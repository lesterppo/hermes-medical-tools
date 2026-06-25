#!/usr/bin/env python3
"""
jamovi integration via R — native, token-efficient statistical analysis.

Instead of requiring the heavy jmv R package (20+ compilation-heavy deps),
this tool uses base R stats functions directly. Same statistical power,
same jamovi-compatible output structure, zero compilation wait.

Available analyses:
  ttestIS, ttestPS, ttestOneS    — t-tests
  anovaOneW                       — one-way ANOVA + Tukey HSD
  linReg                          — linear regression
  logRegBin                       — binary logistic regression
  descriptives                    — N, mean, SD, median, min, max, skew, kurtosis
  corrMatrix                      — Pearson/Spearman correlation matrix
  contTables                      — contingency table + chi-square + Fisher
  wilcoxon                        — Mann-Whitney / Wilcoxon
  kruskal                         — Kruskal-Wallis

Service-gated: only appears when Rscript is on PATH.
"""

from __future__ import annotations

import json
import os
import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import Any, Dict, List, Optional


# ═══════════════════════════════════════════════════════════════════
# Dependency check
# ═══════════════════════════════════════════════════════════════════

_R_LIB = os.path.expanduser("~/R/library")

def _r_available() -> bool:
    return shutil.which("Rscript") is not None


# ═══════════════════════════════════════════════════════════════════
# Helpers
# ═══════════════════════════════════════════════════════════════════

def _ok(result: dict) -> str:
    return json.dumps(result, ensure_ascii=False, separators=(",", ":"))


def _err(msg: str) -> str:
    return _ok({"e": msg})


# ═══════════════════════════════════════════════════════════════════
# R script templates — each analysis gets its own R code block
# ═══════════════════════════════════════════════════════════════════

_R_HEADER = '''
.libPaths(c("{r_lib}", .libPaths()))
data <- read.csv("{data_path}", stringsAsFactors=FALSE, check.names=FALSE)
opts <- jsonlite::fromJSON('{options_json}', simplifyVector=TRUE)
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a
'''

_R_TTESTS = '''
# --- Independent samples t-test ---
grp <- data[[opts$group %||% colnames(data)[1]]]
dep <- data[[opts$dep %||% colnames(data)[2]]]
grp <- as.factor(grp)
lvls <- levels(grp)
if (length(lvls) != 2) stop("Group must have exactly 2 levels")

g1 <- dep[grp == lvls[1]]; g2 <- dep[grp == lvls[2]]

# Welch t-test (default)
tt <- t.test(g1, g2, var.equal=FALSE)
# Effect size (Cohen's d)
n1 <- length(g1); n2 <- length(g2)
sp <- sqrt(((n1-1)*var(g1) + (n2-1)*var(g2)) / (n1+n2-2))
d <- (mean(g1) - mean(g2)) / sp

# Also compute Student's t for comparison
tts <- t.test(g1, g2, var.equal=TRUE)

# Mann-Whitney
mw <- wilcox.test(g1, g2)

res <- list(
  ttest = list(
    welch_t = unname(tt$statistic), welch_df = unname(tt$parameter),
    welch_p = tt$p.value,
    student_t = unname(tts$statistic), student_df = unname(tts$parameter),
    student_p = tts$p.value,
    mean_diff = unname(tt$estimate[1] - tt$estimate[2]),
    ci_lower = tt$conf.int[1], ci_upper = tt$conf.int[2],
    cohens_d = unname(d)
  ),
  descriptives = list(
    group1 = list(n=n1, mean=mean(g1), sd=sd(g1), se=sd(g1)/sqrt(n1)),
    group2 = list(n=n2, mean=mean(g2), sd=sd(g2), se=sd(g2)/sqrt(n2))
  ),
  mann_whitney = list(U=unname(mw$statistic), p=mw$p.value)
)
'''

_R_PAIRED_TTEST = '''
# --- Paired t-test ---
pre <- data[[opts$pre %||% colnames(data)[1]]]
post <- data[[opts$post %||% colnames(data)[2]]]
tt <- t.test(pre, post, paired=TRUE)
d <- mean(pre - post) / sd(pre - post)

res <- list(
  ttest = list(
    t = unname(tt$statistic), df = unname(tt$parameter),
    p = tt$p.value, mean_diff = unname(tt$estimate),
    ci_lower = tt$conf.int[1], ci_upper = tt$conf.int[2],
    cohens_d = unname(d)
  ),
  descriptives = list(
    pre = list(n=length(pre), mean=mean(pre), sd=sd(pre)),
    post = list(n=length(post), mean=mean(post), sd=sd(post))
  )
)
'''

_R_ANOVA = '''
# --- One-way ANOVA ---
dep_name <- opts$dep %||% colnames(data)[2]
grp_name <- opts$group %||% colnames(data)[1]
dep <- data[[dep_name]]
grp <- as.factor(data[[grp_name]])

# ANOVA
fit <- aov(dep ~ grp)
aov_tbl <- summary(fit)[[1]]
# Effect size (eta-squared)
ss_between <- aov_tbl[1, "Sum Sq"]; ss_total <- sum(aov_tbl[, "Sum Sq"])
eta2 <- ss_between / ss_total

# Post-hoc Tukey HSD
tukey <- TukeyHSD(fit)
tukey_df <- as.data.frame(tukey[[1]])

# Homogeneity test — try Levene (car) first, fall back to Bartlett
lv_p <- NA; lv_F <- NA
if (requireNamespace("car", quietly=TRUE)) {
  lv <- car::leveneTest(dep ~ grp)
  lv_F <- lv$`F value`[1]; lv_p <- lv$`Pr(>F)`[1]
} else {
  bt <- bartlett.test(dep ~ grp)
  lv_F <- unname(bt$statistic); lv_p <- bt$p.value
}

# Descriptives per group
desc <- aggregate(dep, list(grp), function(x) c(n=length(x), mean=mean(x), sd=sd(x), se=sd(x)/sqrt(length(x))))

res <- list(
  anova = list(
    F_value = aov_tbl[1, "F value"], df1 = aov_tbl[1, "Df"],
    df2 = aov_tbl[2, "Df"], p = aov_tbl[1, "Pr(>F)"],
    eta_squared = unname(eta2)
  ),
  levene = list(F = lv_F, p = lv_p),
  tukey = list(pairs = rownames(tukey_df), diff = tukey_df[,"diff"],
               lwr = tukey_df[,"lwr"], upr = tukey_df[,"upr"],
               p = tukey_df[,"p adj"]),
  descriptives = desc
)
'''

_R_REGRESSION = '''
# --- Linear regression ---
dep_name <- opts$dep %||% colnames(data)[1]
pred_names <- opts$preds %||% colnames(data)[-1]
fmla <- as.formula(paste(dep_name, "~", paste(pred_names, collapse=" + ")))
fit <- lm(fmla, data=data)
s <- summary(fit)
coefs <- as.data.frame(s$coefficients)
colnames(coefs) <- c("b", "se", "t", "p")

res <- list(
  model = list(
    r_squared = s$r.squared, adj_r_squared = s$adj.r.squared,
    f = unname(s$fstatistic[1]), df1 = unname(s$fstatistic[2]),
    df2 = unname(s$fstatistic[3]), p = unname(pf(s$fstatistic[1], s$fstatistic[2], s$fstatistic[3], lower.tail=FALSE)),
    rmse = sqrt(mean(residuals(fit)^2)), n = nrow(data)
  ),
  coefficients = coefs
)
'''

_R_LOGISTIC = '''
# --- Binary logistic regression ---
dep_name <- opts$dep %||% colnames(data)[1]
pred_names <- opts$preds %||% colnames(data)[-1]
dep <- data[[dep_name]]
dep <- as.factor(dep)

fmla <- as.formula(paste(dep_name, "~", paste(pred_names, collapse=" + ")))
fit <- glm(fmla, data=data, family=binomial)
s <- summary(fit)
coefs <- as.data.frame(s$coefficients)
colnames(coefs) <- c("b", "se", "z", "p")
coefs$OR <- exp(coefs$b)
coefs$ci_lower <- exp(coefs$b - 1.96 * coefs$se)
coefs$ci_upper <- exp(coefs$b + 1.96 * coefs$se)

# Model fit
null_dev <- s$null.deviance; res_dev <- s$deviance
pseudo_r2 <- 1 - res_dev / null_dev
aic <- s$aic

res <- list(
  model = list(
    pseudo_r2 = pseudo_r2, aic = aic,
    null_deviance = null_dev, residual_deviance = res_dev,
    n = nrow(data)
  ),
  coefficients = coefs
)
'''

_R_DESCRIPTIVES = '''
# --- Descriptive statistics ---
vars <- opts$vars %||% colnames(data)
num_data <- data[, sapply(data[vars], is.numeric), drop=FALSE]
# Cap at 15 variables for token efficiency
if (ncol(num_data) > 15) {
  num_data <- num_data[, 1:15, drop=FALSE]
}
desc <- data.frame(
  variable = colnames(num_data),
  n = sapply(num_data, function(x) sum(!is.na(x))),
  missing = sapply(num_data, function(x) sum(is.na(x))),
  mean = sapply(num_data, mean, na.rm=TRUE),
  median = sapply(num_data, median, na.rm=TRUE),
  sd = sapply(num_data, sd, na.rm=TRUE),
  se = sapply(num_data, function(x) sd(x, na.rm=TRUE)/sqrt(sum(!is.na(x)))),
  min = sapply(num_data, min, na.rm=TRUE),
  max = sapply(num_data, max, na.rm=TRUE),
  skew = sapply(num_data, function(x) { n <- sum(!is.na(x)); x <- na.omit(x); m3 <- mean((x-mean(x))^3); m2 <- mean((x-mean(x))^2); m3 / m2^(3/2) * sqrt(n*(n-1))/(n-2) }),
  kurtosis = sapply(num_data, function(x) { x <- na.omit(x); n <- length(x); m4 <- mean((x-mean(x))^4); m2 <- mean((x-mean(x))^2); (n*(n+1)*(n-1)*m4/((n-2)*(n-3)*m2^2)) - 3*(n-1)^2/((n-2)*(n-3)) }),
  shapiro_w = sapply(num_data, function(x) if(length(na.omit(x)) < 3 || length(na.omit(x)) > 5000) NA else shapiro.test(na.omit(x))$statistic),
  shapiro_p = sapply(num_data, function(x) if(length(na.omit(x)) < 3 || length(na.omit(x)) > 5000) NA else shapiro.test(na.omit(x))$p.value)
)
res <- list(descriptives = desc)
'''

_R_CORRELATION = '''
# --- Correlation matrix ---
vars <- opts$vars %||% colnames(data)
num_data <- data[, sapply(data[vars], is.numeric), drop=FALSE]
method <- opts$method %||% "pearson"

cor_mat <- cor(num_data, use="pairwise.complete.obs", method=method)
n <- nrow(num_data)

# P-values
p_mat <- matrix(NA, ncol(cor_mat), ncol(cor_mat))
colnames(p_mat) <- colnames(cor_mat); rownames(p_mat) <- colnames(cor_mat)
for (i in 1:ncol(num_data)) {
  for (j in 1:ncol(num_data)) {
    if (i != j) {
      ct <- cor.test(num_data[[i]], num_data[[j]], method=method)
      p_mat[i,j] <- ct$p.value
    }
  }
}

res <- list(
  correlation = as.data.frame(cor_mat),
  p_values = as.data.frame(p_mat),
  n = n, method = method
)
'''

_R_CHISQ = '''
# --- Contingency table + chi-square ---
row_name <- opts$rows %||% colnames(data)[1]
col_name <- opts$cols %||% colnames(data)[2]
tbl <- table(data[[row_name]], data[[col_name]])
csq <- chisq.test(tbl)
fisher_p <- if (nrow(tbl) == 2 && ncol(tbl) == 2) fisher.test(tbl)$p.value else NA

# Cramer's V
n <- sum(tbl)
chi2 <- unname(csq$statistic)
k <- min(nrow(tbl), ncol(tbl))
cramers_v <- sqrt(chi2 / (n * (k - 1)))

res <- list(
  table = as.data.frame.matrix(tbl),
  chi_square = list(chi2 = chi2, df = unname(csq$parameter), p = csq$p.value),
  fisher_exact = fisher_p,
  cramers_v = unname(cramers_v),
  n = n
)
'''

_R_WILCOXON = '''
# --- Wilcoxon / Mann-Whitney ---
grp <- as.factor(data[[opts$group %||% colnames(data)[1]]])
dep <- data[[opts$dep %||% colnames(data)[2]]]
lvls <- levels(grp)
if (length(lvls) == 2) {
  wt <- wilcox.test(dep ~ grp)
  res <- list(test="Mann-Whitney U", W=unname(wt$statistic), p=wt$p.value,
              n1=sum(grp==lvls[1]), n2=sum(grp==lvls[2]))
} else {
  kw <- kruskal.test(dep ~ grp)
  res <- list(test="Kruskal-Wallis", chi2=unname(kw$statistic), df=unname(kw$parameter),
              p=kw$p.value, n_groups=length(lvls), n=nrow(data))
}
'''

# Map analysis name → R script
_R_SCRIPTS = {
    "ttestIS": _R_TTESTS,
    "ttestPS": _R_PAIRED_TTEST,
    "anovaOneW": _R_ANOVA,
    "linReg": _R_REGRESSION,
    "logRegBin": _R_LOGISTIC,
    "descriptives": _R_DESCRIPTIVES,
    "corrMatrix": _R_CORRELATION,
    "contTables": _R_CHISQ,
    "wilcoxon": _R_WILCOXON,
}


# ═══════════════════════════════════════════════════════════════════
# Tool: jmv — Statistical analysis via R (base R stats, no heavy deps)
# ═══════════════════════════════════════════════════════════════════

def jmv_run(
    analysis: str,
    data: Optional[str] = None,
    data_file: Optional[str] = None,
    options: Optional[str] = None,
) -> str:
    """
    Statistical analysis via R. Base R stats — no heavy package installs.

    analysis: ttestIS, ttestPS, anovaOneW, linReg, logRegBin,
              descriptives, corrMatrix, contTables, wilcoxon

    data:     CSV with header: 'group,score\\n1,10\\n1,12\\n2,20\\n2,22'
    data_file: Path to CSV file (use instead of data for large files)
    options:  JSON: '{"group":"group","dep":"score"}' or omit for auto-detect
    """
    try:
        valid = set(_R_SCRIPTS.keys())
        if analysis not in valid:
            return _err(f"Unknown analysis '{analysis}'. Valid: {', '.join(sorted(valid))}")

        # Handle data
        if data:
            csv_data = data
        elif data_file:
            if not os.path.exists(data_file):
                return _err(f"File not found: {data_file}")
            with open(data_file) as f:
                csv_data = f.read()
        else:
            return _err("Provide data (inline CSV) or data_file (path to CSV)")

        # Parse options
        opts = {}
        if options:
            try:
                opts = json.loads(options)
            except json.JSONDecodeError:
                return _err("options must be valid JSON")

        # Write temp files
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".csv", delete=False, encoding="utf-8"
        ) as f:
            f.write(csv_data)
            data_path = f.name

        try:
            # Build R script
            script = _R_HEADER.replace("{data_path}", data_path).replace(
                "{r_lib}", _R_LIB
            ).replace("{options_json}", json.dumps(opts))
            script += _R_SCRIPTS[analysis]
            script += "\ncat(jsonlite::toJSON(res, dataframe='rows', na='null', digits=6, auto_unbox=FALSE))\n"

            result = subprocess.run(
                ["Rscript", "-"],
                input=script,
                capture_output=True,
                text=True,
                timeout=30,
                env={**os.environ, "R_LIBS_USER": _R_LIB, "R_LIBS_SITE": ""},
            )

            if result.returncode != 0:
                err_lines = [l for l in result.stderr.strip().split("\n")
                           if l.strip() and "Error" in l or "rror" in l]
                err_msg = err_lines[-1][:400] if err_lines else result.stderr.strip()[:400]
                return _err(f"R error: {err_msg}")

            # Check for R errors in output
            stdout = result.stdout.strip()
            if stdout.startswith("Error"):
                return _err(f"R error: {stdout[:400]}")

            try:
                raw = json.loads(stdout)
            except json.JSONDecodeError:
                return _ok({"out": "text", "t": stdout[:3000]})

            return _ok({
                "analysis": analysis,
                "r": raw,
            })

        finally:
            os.unlink(data_path)

    except subprocess.TimeoutExpired:
        return _err("R analysis timed out (30s)")
    except Exception as e:
        return _err(f"jmv error: {e}")


# ═══════════════════════════════════════════════════════════════════
# Schema
# ═══════════════════════════════════════════════════════════════════

JMV_SCHEMA = {
    "name": "jmv",
    "description": (
        "Statistical analysis via R (base stats, no heavy deps). "
        "analysis: ttestIS|ttestPS|anovaOneW|linReg|logRegBin|"
        "descriptives|corrMatrix|contTables|wilcoxon. "
        "data=CSV with header. options=JSON options dict (omit for defaults). "
        "Auto-detects column roles from column order."
    ),
    "parameters": {
        "type": "object",
        "properties": {
            "analysis": {
                "type": "string",
                "description": "ttestIS, ttestPS, anovaOneW, linReg, logRegBin, descriptives, corrMatrix, contTables, wilcoxon",
            },
            "data": {
                "type": "string",
                "description": "CSV data with header: 'group,score\\n1,10\\n1,12\\n2,20\\n2,22'",
            },
            "data_file": {
                "type": "string",
                "description": "Path to CSV file. Use instead of data for large files.",
            },
            "options": {
                "type": "string",
                "description": "JSON: '{\"group\":\"tx\",\"dep\":\"score\"}'. Omit for auto-detection.",
            },
        },
        "required": ["analysis"],
    },
}


# ═══════════════════════════════════════════════════════════════════
# Registry
# ═══════════════════════════════════════════════════════════════════

from tools.registry import registry

registry.register(
    name="jmv",
    toolset="medical",
    schema=JMV_SCHEMA,
    handler=lambda args, **kw: jmv_run(
        analysis=args.get("analysis", ""),
        data=args.get("data"),
        data_file=args.get("data_file"),
        options=args.get("options"),
    ),
    check_fn=_r_available,
    emoji="📈",
)
