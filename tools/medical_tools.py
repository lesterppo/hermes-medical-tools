#!/usr/bin/env python3
"""
Medical research tools — AI-agent-native, token-efficient.

Every tool returns compact JSON with short keys. Schemas use terse
descriptions. All tools are service-gated: they only appear when
dependencies are available (scipy/numpy for stats; pubmed API key
or httpx for PubMed).

Token budget for ALL 5 tool schemas combined: ~1,500 chars.
Output per call: typically 200-600 chars of JSON.
"""

from __future__ import annotations

import json
import os
import sys
import time
from typing import Any, Dict, List, Optional, Tuple

# ═══════════════════════════════════════════════════════════════════
# Dependency checks
# ═══════════════════════════════════════════════════════════════════

def _has_scipy() -> bool:
    try:
        import scipy  # noqa: F401
        return True
    except ImportError:
        return False


def _has_pubmed() -> bool:
    """PubMed needs httpx or urllib (stdlib always available)."""
    # E-utilities only requires urllib (stdlib). httpx is optional for async.
    return True  # Always available — uses stdlib


# ═══════════════════════════════════════════════════════════════════
# Shared helpers
# ═══════════════════════════════════════════════════════════════════

def _sanitize(v):
    """Recursively replace NaN/±inf floats with None for strict JSON."""
    import math
    try:
        import numpy as _np
        if isinstance(v, _np.floating):
            v = float(v)
        elif isinstance(v, _np.integer):
            v = int(v)
        elif isinstance(v, _np.ndarray):
            return [_sanitize(x) for x in v.tolist()]
    except ImportError:
        pass
    if isinstance(v, float):
        return None if (math.isnan(v) or math.isinf(v)) else v
    if isinstance(v, dict):
        return {k: _sanitize(x) for k, x in v.items()}
    if isinstance(v, (list, tuple)):
        return [_sanitize(x) for x in v]
    return v


def _ok(result: dict) -> str:
    """Wrap result dict as compact JSON. No newlines, minimal whitespace."""
    return json.dumps(_sanitize(result), ensure_ascii=False, separators=(",", ":"),
                      allow_nan=False, default=str)


def _err(msg: str) -> str:
    return _ok({"e": msg})


def _parse_numbers(raw: Any) -> Optional[List[float]]:
    """Parse a list of numbers from various input formats."""
    if raw is None:
        return None
    if isinstance(raw, list):
        try:
            return [float(x) for x in raw]
        except (ValueError, TypeError):
            return None
    if isinstance(raw, str):
        parts = raw.replace(",", " ").split()
        try:
            return [float(p) for p in parts]
        except ValueError:
            return None
    return None


# ═══════════════════════════════════════════════════════════════════
# Tool: med_stats — Statistical tests for medical research
# ═══════════════════════════════════════════════════════════════════

def _wilson_ci(k: int, n: int, z: float = 1.96) -> tuple:
    """Wilson score interval — the correct CI for a proportion (never 0-width
    at 0/n or n/n, unlike the Wald interval)."""
    if n <= 0:
        return (0.0, 1.0)
    p = k / n
    d = 1 + z * z / n
    c = p + z * z / (2 * n)
    half = z * ((p * (1 - p) / n + z * z / (4 * n * n)) ** 0.5)
    return (max(0.0, (c - half) / d), min(1.0, (c + half) / d))


def _mean_ci(x, alpha: float = 0.05) -> tuple:
    """95% CI for a mean (t-based)."""
    import numpy as np
    from scipy import stats as sps
    n = len(x)
    if n < 2:
        return (float("nan"), float("nan"))
    se = float(np.std(x, ddof=1)) / (n ** 0.5)
    h = float(sps.t.ppf(1 - alpha / 2, n - 1)) * se
    m = float(np.mean(x))
    return (round(m - h, 3), round(m + h, 3))


def med_stats(
    test: str,
    a: Optional[List[float]] = None,
    b: Optional[List[float]] = None,
    categorical: Optional[List[List[int]]] = None,
    paired: bool = False,
) -> str:
    """
    Run a statistical test, returning compact results.

    test:  "ttest" | "mw" | "wilcoxon" | "chisq" | "fisher" | "anova" | "kw" | "pearson" | "spearman"
    a:     first group data (numbers) or paired pre-values. Omit when using categorical.
    b:     second group data (numbers) or paired post-values; omit for 1-sample
    categorical: 2x2 table [[a,b],[c,d]] for chi-square/fisher; use INSTEAD of a/b
    paired: True for paired t-test or Wilcoxon signed-rank

    Effect sizes and 95% CIs are returned alongside p-values: a p-value alone
    is not reportable in a clinical manuscript.
    """
    try:
        import numpy as np
        from scipy import stats as sps

        # --- Categorical tests (checked FIRST: the documented call passes a=[])
        if categorical is not None:
            if not (isinstance(categorical, list) and len(categorical) == 2
                    and len(categorical[0]) == 2 and len(categorical[1]) == 2):
                return _err("categorical must be [[a,b],[c,d]]")
            tbl = np.array(categorical, dtype=float)
            if tbl.sum() <= 0:
                return _err("categorical table totals zero cells")
            n = int(tbl.sum())
            a_, b_, c_, d_ = (float(tbl[0][0]), float(tbl[0][1]),
                              float(tbl[1][0]), float(tbl[1][1]))
            out: dict = {"n": n}
            # 2x2 effect sizes: OR + risk ratio (with Haldane-Anscombe correction
            # when a cell is zero, so the estimate is not silently inf).
            zero_cell = 0 in (a_, b_, c_, d_)
            aa, bb, cc, dd = (a_ + 0.5, b_ + 0.5, c_ + 0.5, d_ + 0.5) if zero_cell else (a_, b_, c_, d_)
            or_v = (aa * dd) / (bb * cc) if bb and cc else float("inf")
            if (a_ + b_) > 0 and (c_ + d_) > 0:
                out["rr"] = round(((aa / (aa + bb)) / (cc / (cc + dd))), 3)
            out["or"] = round(or_v, 3) if or_v != float("inf") else "inf"
            if zero_cell:
                out["note"] = "zero cell — Haldane-Anscombe (+0.5) correction applied to OR/RR"
            if test == "fisher":
                _, p = sps.fisher_exact(tbl)
                out.update({"t": "fisher", "p": round(float(p), 6)})
                return _ok(out)
            if test == "chisq":
                chi2, p, dof, expected = sps.chi2_contingency(tbl)
                out.update({
                    "t": "chisq", "x2": round(float(chi2), 3),
                    "df": dof, "p": round(float(p), 6),
                    # Cramér's V is the reportable effect size for a table.
                    "v": round(float((chi2 / (n * min(tbl.shape[0] - 1, tbl.shape[1] - 1))) ** 0.5), 3),
                })
                if float(expected.min()) < 5:
                    out["warn"] = f"min expected count {round(float(expected.min()),1)} < 5 — Fisher exact is preferred"
                return _ok(out)
            return _err(f"test='{test}' not for categorical; use chisq or fisher")

        # --- Input validation for continuous tests ---
        if a is None or (isinstance(a, list) and len(a) == 0):
            return _err("a must be a non-empty list of numbers (or pass categorical=[[a,b],[c,d]])")

        # --- Continuous tests ---
        x = np.array(a, dtype=float)
        y = np.array(b, dtype=float) if b is not None else None

        if test == "ttest":
            if paired:
                if y is None:
                    return _err("paired ttest needs both a and b")
                if len(x) != len(y):
                    return _err(f"paired ttest needs equal-length groups (n1={len(x)}, n2={len(y)})")
                if len(x) < 2:
                    return _err("paired ttest needs at least 2 pairs")
                stat, p = sps.ttest_rel(x, y)
            else:
                if y is not None and len(y) < 2:
                    return _err(f"group b has n={len(y)} — need at least 2 per group")
                if y is None and len(x) < 2:
                    return _err(f"group a has n={len(x)} — need at least 2")
                stat, p = sps.ttest_ind(x, y) if y is not None else sps.ttest_1samp(x, 0)
            d_val = _cohens_d(x, y, paired)
            diff = x - y if (paired and y is not None) else (x if y is None else x - y.mean())
            res = {
                "t": "ttest", "s": round(float(stat), 3),
                "p": round(float(p), 6), "n1": len(x),
                "n2": len(y) if y is not None else 0,
                "m1": round(float(x.mean()), 2),
                "m2": round(float(y.mean()), 2) if y is not None else 0,
                "ci1": _mean_ci(x),
                "d": round(d_val, 2) if d_val == d_val else None,
            }
            if y is not None:
                res["m2_ci"] = _mean_ci(y)
                res["md"] = round(float(x.mean() - y.mean()), 3)
            if d_val != d_val:
                res["warn"] = "effect size undefined (zero variance within group)"
            return _ok(res)

        elif test == "mw":
            if y is None:
                return _err("Mann-Whitney needs two groups (a and b)")
            if len(x) < 2 or len(y) < 2:
                return _err(f"Mann-Whitney needs n>=2 per group (n1={len(x)}, n2={len(y)})")
            stat, p = sps.mannwhitneyu(x, y, alternative="two-sided")
            # Rank-biserial correlation (signed: positive = a tends higher).
            u = float(stat)
            rb = 1 - (2 * u) / (len(x) * len(y))
            return _ok({
                "t": "mw", "u": u, "p": round(float(p), 6),
                "n1": len(x), "n2": len(y),
                "md1": round(float(np.median(x)), 2),
                "md2": round(float(np.median(y)), 2),
                "rb": round(rb, 3),
            })

        elif test == "wilcoxon":
            if y is None:
                return _err("Wilcoxon needs a and b")
            if len(x) != len(y):
                return _err(f"Wilcoxon signed-rank needs paired equal-length data (n1={len(x)}, n2={len(y)})")
            if len(x) < 2:
                return _err("Wilcoxon needs at least 2 pairs")
            if np.allclose(x - y, 0):
                return _err("Wilcoxon undefined: all paired differences are zero")
            stat, p = sps.wilcoxon(x, y)
            return _ok({
                "t": "wilcoxon", "s": float(stat), "p": round(float(p), 6),
                "n": len(x),
                "md": round(float(np.median(x - y)), 3),
            })

        elif test == "anova":
            # Multi-group: a = [[g1...], [g2...], ...]
            if isinstance(a, list) and len(a) > 0 and isinstance(a[0], list):
                groups = [np.array(g, dtype=float) for g in a if len(g) > 0]
                if len(groups) < 2:
                    return _err("ANOVA needs at least 2 non-empty groups")
                if any(len(g) < 2 for g in groups):
                    return _err("ANOVA needs n>=2 in every group")
                stat, p = sps.f_oneway(*groups)
                n_total = sum(len(g) for g in groups)
                grand = np.concatenate(groups).mean()
                ss_between = sum(len(g) * (g.mean() - grand) ** 2 for g in groups)
                ss_total = sum(((g - grand) ** 2).sum() for g in groups)
                return _ok({
                    "t": "anova", "f": round(float(stat), 3),
                    "df1": len(groups) - 1, "df2": n_total - len(groups),
                    "p": round(float(p), 6), "k": len(groups),
                    "means": [round(float(g.mean()), 2) for g in groups],
                    "eta2": round(float(ss_between / ss_total), 3) if ss_total else None,
                })
            # Flat a + b for 2-group
            if y is None or len(y) == 0:
                return _err("ANOVA needs at least 2 groups — pass a=[[g1],[g2],...] or a=g1,b=g2")
            if len(x) < 2 or len(y) < 2:
                return _err("ANOVA needs n>=2 per group")
            stat, p = sps.f_oneway(x, y)
            grand = np.concatenate([x, y]).mean()
            ss_between = len(x) * (x.mean() - grand) ** 2 + len(y) * (y.mean() - grand) ** 2
            ss_total = ((x - grand) ** 2).sum() + ((y - grand) ** 2).sum()
            return _ok({
                "t": "anova", "f": round(float(stat), 3),
                "df1": 1, "df2": len(x) + len(y) - 2,
                "p": round(float(p), 6),
                "eta2": round(float(ss_between / ss_total), 3) if ss_total else None,
                "means": [round(float(x.mean()), 2), round(float(y.mean()), 2)],
            })

        elif test == "kw":
            # Multi-group: a = [[g1...], [g2...], ...]
            if isinstance(a, list) and len(a) > 0 and isinstance(a[0], list):
                groups = [np.array(g, dtype=float) for g in a if len(g) > 0]
                if len(groups) < 2:
                    return _err("Kruskal-Wallis needs at least 2 non-empty groups")
                if any(len(g) < 2 for g in groups):
                    return _err("Kruskal-Wallis needs n>=2 in every group")
                stat, p = sps.kruskal(*groups)
                n_total = sum(len(g) for g in groups)
                return _ok({
                    "t": "kw", "h": round(float(stat), 3),
                    "p": round(float(p), 6), "k": len(groups),
                    "n": n_total,
                    # eta-squared H is the standard KW effect size.
                    "eta2h": round(float(stat / (n_total - 1)), 3) if n_total > 1 else None,
                    "medians": [round(float(np.median(g)), 2) for g in groups],
                })
            # Flat a + b for 2-group
            if y is None or len(y) == 0:
                return _err("Kruskal-Wallis needs at least 2 groups — pass a=[[g1],[g2],...] or a=g1,b=g2")
            if len(x) < 2 or len(y) < 2:
                return _err("Kruskal-Wallis needs n>=2 per group")
            stat, p = sps.kruskal(x, y)
            n_total = len(x) + len(y)
            return _ok({
                "t": "kw", "h": round(float(stat), 3),
                "p": round(float(p), 6), "n": n_total,
                "eta2h": round(float(stat / (n_total - 1)), 3) if n_total > 1 else None,
                "medians": [round(float(np.median(x)), 2), round(float(np.median(y)), 2)],
            })

        elif test == "pearson":
            if y is None:
                return _err("Pearson needs x and y")
            if len(x) != len(y):
                return _err(f"Pearson needs paired equal-length data (n1={len(x)}, n2={len(y)})")
            if len(x) < 3:
                return _err(f"Pearson needs n>=3 (got {len(x)})")
            r, p = sps.pearsonr(x, y)
            # Fisher z CI — correlations must be reported with a CI.
            ci = None
            if len(x) > 3 and abs(float(r)) < 1:
                z = np.arctanh(float(r))
                se = 1 / np.sqrt(len(x) - 3)
                ci = (round(float(np.tanh(z - 1.96 * se)), 3), round(float(np.tanh(z + 1.96 * se)), 3))
            return _ok({
                "t": "pearson", "r": round(float(r), 4),
                "p": round(float(p), 6), "n": len(x),
                "ci": ci, "r2": round(float(r) ** 2, 4),
            })

        elif test == "spearman":
            if y is None:
                return _err("Spearman needs x and y")
            if len(x) != len(y):
                return _err(f"Spearman needs paired equal-length data (n1={len(x)}, n2={len(y)})")
            if len(x) < 3:
                return _err(f"Spearman needs n>=3 (got {len(x)})")
            rho, p = sps.spearmanr(x, y)
            return _ok({
                "t": "spearman", "rho": round(float(rho), 4),
                "p": round(float(p), 6), "n": len(x),
            })

        else:
            return _err(f"unknown test='{test}'. Use: ttest, mw, wilcoxon, chisq, fisher, anova, kw, pearson, spearman")

    except Exception as e:
        return _err(f"stats error: {e}")


def _cohens_d(x: "np.ndarray", y: "np.ndarray" = None, paired: bool = False) -> float:
    """Cohen's d effect size."""
    import numpy as np
    if paired or y is None:
        diff = x - y if y is not None else x
        return float(np.mean(diff) / np.std(diff, ddof=1))
    n1, n2 = len(x), len(y)
    sp = np.sqrt(((n1 - 1) * np.var(x, ddof=1) + (n2 - 1) * np.var(y, ddof=1)) / (n1 + n2 - 2))
    return float((np.mean(x) - np.mean(y)) / sp)


# ═══════════════════════════════════════════════════════════════════
# Tool: med_pubmed — PubMed search
# ═══════════════════════════════════════════════════════════════════

_PUBMED_SEARCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
_PUBMED_FETCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"

def med_pubmed(
    q: str = "",
    max_results: int = 5,
    fetch: bool = False,
    pmid: Optional[str] = None,
    retstart: int = 0,
) -> str:
    """
    Search PubMed or fetch article details.

    q:           search query (PubMed syntax: "diabetes AND metformin[TIAB] AND 2024[DP]")
    max_results: articles to return (max 20)
    fetch:       if True, fetches abstracts; if False, returns PMIDs + titles
    pmid:        single PMID to fetch directly; overrides q
    retstart:    pagination offset
    """
    import urllib.request
    import urllib.parse
    import xml.etree.ElementTree as ET

    API_KEY = os.getenv("PUBMED_API_KEY", "")
    TOOL = "hermes-med"
    EMAIL = os.getenv("PUBMED_EMAIL", "user@example.com")

    try:
        if pmid:
            # Direct fetch
            params = {
                "db": "pubmed", "id": str(pmid), "retmode": "xml",
                "tool": TOOL, "email": EMAIL,
            }
            if API_KEY:
                params["api_key"] = API_KEY
            url = _PUBMED_FETCH + "?" + urllib.parse.urlencode(params)
            resp = urllib.request.urlopen(
                urllib.request.Request(url, headers={"User-Agent": f"{TOOL}/1.0"}),
                timeout=15,
            )
            data = resp.read()
            root = ET.fromstring(data)
            article = root.find(".//PubmedArticle")
            if article is None:
                return _ok({"pmid": pmid, "e": "not found"})

            title_el = article.find(".//ArticleTitle")
            title = "".join(title_el.itertext()).strip() if title_el is not None else ""
            # Join ALL AbstractText nodes — structured abstracts
            # (Background/Methods/Results/Conclusions) have several.
            abstract_parts = []
            for ab_el in article.findall(".//AbstractText"):
                label = ab_el.get("Label")
                txt = "".join(ab_el.itertext()).strip()
                if txt:
                    abstract_parts.append(f"{label}: {txt}" if label else txt)
            abstract = " ".join(abstract_parts)
            yr_el = article.find(".//PubDate/Year")
            yr = yr_el.text if yr_el is not None else ""
            jr_el = article.find(".//Journal/Title")
            jr = jr_el.text if jr_el is not None else ""

            # DOI, publication types and retraction status: without these a
            # clinician cannot cite the paper, and can unknowingly cite a
            # retracted one.
            doi = ""
            for el in article.findall(".//ArticleId"):
                if (el.get("IdType") or "").lower() == "doi" and el.text:
                    doi = el.text.strip()
                    break
            if not doi:
                for el in article.findall(".//ELocationID"):
                    if (el.get("EIdType") or "").lower() == "doi" and el.text:
                        doi = el.text.strip()
                        break
            pubs = [p.text.strip() for p in article.findall(".//PublicationType") if p.text]
            warn = ""
            low = {p.lower() for p in pubs}
            if {"retracted publication", "retraction of publication"} & low:
                warn = "RETRACTED"
            elif "expression of concern" in low:
                warn = "EXPRESSION_OF_CONCERN"
            elif low & {"published erratum", "corrected and republished article"}:
                warn = "CORRECTED"

            out = {
                "pmid": pmid, "t": title[:400],
                "ab": abstract[:900] if abstract else "",
                "yr": yr, "jr": jr,
            }
            if doi:
                out["doi"] = doi
            if pubs:
                out["pt"] = pubs[:4]
            if warn:
                out["warn"] = warn
            return _ok(out)

        # Search
        params = {
            "db": "pubmed", "term": q, "retmax": min(max_results, 20),
            "retstart": retstart, "retmode": "xml",
            "sort": "relevance", "tool": TOOL, "email": EMAIL,
        }
        if API_KEY:
            params["api_key"] = API_KEY
        url = _PUBMED_SEARCH + "?" + urllib.parse.urlencode(params)
        resp = urllib.request.urlopen(
            urllib.request.Request(url, headers={"User-Agent": f"{TOOL}/1.0"}),
            timeout=15,
        )
        data = resp.read()
        root = ET.fromstring(data)

        ids = [e.text for e in root.findall(".//Id") if e.text]
        count_el = root.find(".//Count")
        total = int(count_el.text) if count_el is not None else 0

        if not ids:
            return _ok({"q": q, "n": 0, "total": total, "pmids": []})

        if fetch and ids:
            # Batch fetch abstracts
            fparams = {
                "db": "pubmed", "id": ",".join(ids), "retmode": "xml",
                "tool": TOOL, "email": EMAIL,
            }
            if API_KEY:
                fparams["api_key"] = API_KEY
            furl = _PUBMED_FETCH + "?" + urllib.parse.urlencode(fparams)
            fresp = urllib.request.urlopen(
                urllib.request.Request(furl, headers={"User-Agent": f"{TOOL}/1.0"}),
                timeout=20,
            )
            fdata = fresp.read()
            froot = ET.fromstring(fdata)

            articles = []
            for art in froot.findall(".//PubmedArticle"):
                pid_el = art.find(".//PMID")
                pid = pid_el.text if pid_el is not None else ""
                title_el = art.find(".//ArticleTitle")
                title = "".join(title_el.itertext()).strip()[:300] if title_el is not None else ""
                ab_parts = []
                for ab_el in art.findall(".//AbstractText"):
                    label = ab_el.get("Label")
                    txt = "".join(ab_el.itertext()).strip()
                    if txt:
                        ab_parts.append(f"{label}: {txt}" if label else txt)
                abstract = " ".join(ab_parts)[:400]
                yr_el = art.find(".//PubDate/Year")
                yr = yr_el.text if yr_el is not None else ""
                entry = {
                    "pmid": pid, "t": title,
                    "ab": abstract, "yr": yr,
                }
                for el in art.findall(".//ArticleId"):
                    if (el.get("IdType") or "").lower() == "doi" and el.text:
                        entry["doi"] = el.text.strip()
                        break
                low = {(p.text or "").strip().lower() for p in art.findall(".//PublicationType")}
                if {"retracted publication", "retraction of publication"} & low:
                    entry["warn"] = "RETRACTED"
                elif "expression of concern" in low:
                    entry["warn"] = "EXPRESSION_OF_CONCERN"
                articles.append(entry)
            return _ok({"q": q, "n": len(articles), "total": total, "r": articles})

        return _ok({"q": q, "n": len(ids), "total": total, "pmids": ids})

    except Exception as e:
        import traceback
        return _err(f"PubMed error: {e}")


# ═══════════════════════════════════════════════════════════════════
# Tool: med_trial — ClinicalTrials.gov search
# ═══════════════════════════════════════════════════════════════════

_CLINICALTRIALS_API = "https://clinicaltrials.gov/api/v2/studies"

def med_trial(
    q: str,
    max_results: int = 5,
    status: Optional[str] = None,
    fmt: str = "brief",
) -> str:
    """
    Search ClinicalTrials.gov.

    q:           condition or drug (e.g., "pembrolizumab melanoma")
    max_results: max studies to return (max 20)
    status:      filter: "recruiting" | "active" | "completed" | None
    fmt:         output: "brief" (NCT + title + phase + status) | "full" (adds conditions + locations)
    """
    import urllib.request
    import urllib.parse

    try:
        params = {
            "query.term": q,
            "pageSize": min(max_results, 20),
            "format": "json",
            # Without countTotal the v2 API omits totalCount and the tool
            # reported total=0 for every search (wrong screening counts).
            "countTotal": "true",
        }
        if status:
            # 'active' is not a valid v2 overallStatus enum — map the
            # colloquial term to ACTIVE_NOT_RECRUITING (HTTP 400 otherwise).
            _status_map = {"active": "ACTIVE_NOT_RECRUITING"}
            params["filter.overallStatus"] = _status_map.get(status.lower(), status.upper())

        url = _CLINICALTRIALS_API + "?" + urllib.parse.urlencode(params)
        req = urllib.request.Request(url, headers={
            "Accept": "application/json",
            "User-Agent": "hermes-med/1.0",
        })
        resp = urllib.request.urlopen(req, timeout=15)
        data = json.loads(resp.read())

        studies = data.get("studies", [])
        total = data.get("totalCount", 0)
        results = []

        for s in studies:
            prot = s.get("protocolSection", {})
            ident = prot.get("identificationModule", {})
            design = prot.get("designModule", {})
            stat = prot.get("statusModule", {})

            nct = ident.get("nctId", "")
            title = ident.get("briefTitle", "")[:300]
            phase_list = design.get("phases", [])
            phase = phase_list[0] if phase_list else "N/A"
            overall = stat.get("overallStatus", "")

            entry = {
                "nct": nct, "t": title, "ph": phase, "st": overall,
            }
            if fmt == "full":
                cond = prot.get("conditionsModule", {}).get("conditions", [])
                entry["cond"] = cond[:3]
                locs = prot.get("contactsLocationsModule", {}).get("locations", [])
                entry["loc"] = [l.get("facility", "")[:80] for l in locs[:3]]

            results.append(entry)

        return _ok({"q": q, "n": len(results), "total": total, "r": results})

    except Exception as e:
        return _err(f"trial error: {e}")


# ═══════════════════════════════════════════════════════════════════
# Tool: med_power — Sample size / power calculations
# ═══════════════════════════════════════════════════════════════════

def med_power(
    calc: str,
    effect: float = 0.5,
    power: float = 0.80,
    alpha: float = 0.05,
    ratio: float = 1.0,
    p1: Optional[float] = None,
    p2: Optional[float] = None,
    n: Optional[int] = None,
    dropout: float = 0.0,
) -> str:
    """
    Calculate sample size, power, or detectable effect for common designs.

    calc:    "n" (sample size needed) | "power" (actual power given n) | "detect" (detectable effect)
    effect:  Cohen's d (0.2 small, 0.5 medium, 0.8 large). For "detect"/legacy "power" with
             no `n`, this doubles as n per group.
    power:   target power (default 0.80)
    alpha:   significance level (default 0.05, two-sided)
    ratio:   n2/n1 ratio (default 1 for equal groups)
    p1:      proportion in group 1 (proportion tests)
    p2:      proportion in group 2 (proportion tests)
    n:       per-group n for calc="power" (exact noncentral-t power)
    dropout: expected attrition fraction (0.15 = 15%) — inflates the required n
    """
    try:
        from scipy import stats as sps
        import numpy as np

        if dropout and not (0 <= dropout < 1):
            return _err("dropout must be a fraction between 0 and 1 (0.15 = 15%)")

        def _inflate(n1: int) -> tuple:
            """Return (n_per_group_after_dropout, total_after_dropout)."""
            if not dropout:
                return n1, n1 + round(n1 * ratio)
            adj = int(np.ceil(n1 / (1 - dropout)))
            return adj, adj + round(adj * ratio)

        # For proportions
        if p1 is not None and p2 is not None:
            if not (0 < p1 < 1 and 0 < p2 < 1):
                return _err("p1 and p2 must be proportions strictly between 0 and 1")
            # Cohen's h (arcsine transformation)
            h = 2 * (np.arcsin(np.sqrt(p1)) - np.arcsin(np.sqrt(p2)))
            if calc == "n":
                n1 = _n_for_effect(h, power, alpha, ratio)
                adj, tot = _inflate(n1)
                return _ok({
                    "t": "proportions", "calc": "n",
                    "h": round(float(h), 4), "p1": p1, "p2": p2,
                    "n1": adj, "n2": round(adj * ratio), "total": tot,
                    "n1_before_dropout": n1,
                    "dropout": dropout or None,
                    "alpha": alpha, "power": power,
                })
            elif calc == "power":
                n1 = int(n or effect)
                if n1 < 2:
                    return _err("n (per group) must be >= 2 for a power calculation")
                z_alpha = sps.norm.ppf(1 - alpha / 2)
                z_power = abs(h) * np.sqrt(n1 / (1 + 1 / ratio)) - z_alpha
                actual_power = float(sps.norm.cdf(z_power))
                return _ok({"t": "proportions", "calc": "power", "n1": n1, "h": round(float(h), 4),
                            "p": round(actual_power, 4)})
            return _err("calc must be n or power")

        # For continuous (Cohen's d) — exact noncentral-t, not a normal
        # approximation: the normal approximation overstates power noticeably
        # at the small n typical of pilot studies.
        if calc == "n":
            n1 = _n_for_effect(effect, power, alpha, ratio)
            adj, tot = _inflate(n1)
            return _ok({
                "t": "continuous", "calc": "n",
                "d": effect, "n1": adj, "n2": round(adj * ratio), "total": tot,
                "n1_before_dropout": n1,
                "dropout": dropout or None,
                "alpha": alpha, "power": power,
            })
        elif calc == "power":
            n1 = int(n or effect)
            if n1 < 2:
                return _err("n (per group) must be >= 2 for a power calculation")
            n2 = max(2, int(round(n1 * ratio)))
            df = n1 + n2 - 2
            nc = effect * np.sqrt(1.0 / (1.0 / n1 + 1.0 / n2))
            crit = sps.t.ppf(1 - alpha / 2, df)
            actual_power = float(1 - sps.nct.cdf(crit, df, nc) + sps.nct.cdf(-crit, df, nc))
            return _ok({
                "t": "continuous", "calc": "power",
                "n1": n1, "n2": n2, "d": effect,
                "p": round(actual_power, 4),
                "alpha": alpha,
            })
        elif calc == "detect":
            n1 = int(n or effect)
            if n1 < 2:
                return _err("n (per group) must be >= 2 for a detectable-effect calculation")
            # Solve for the smallest d whose noncentral-t power reaches target.
            lo, hi = 1e-4, 5.0
            for _ in range(60):
                mid = (lo + hi) / 2
                df = n1 + max(2, int(round(n1 * ratio))) - 2
                nc_ = mid * np.sqrt(1.0 / (1.0 / n1 + 1.0 / (max(2, int(round(n1 * ratio))))))
                crit = sps.t.ppf(1 - alpha / 2, df)
                pw = float(1 - sps.nct.cdf(crit, df, nc_) + sps.nct.cdf(-crit, df, nc_))
                if pw < power:
                    lo = mid
                else:
                    hi = mid
            return _ok({
                "t": "continuous", "calc": "detect",
                "n1": n1, "d": round(float((lo + hi) / 2), 4),
                "alpha": alpha, "power": power,
            })
        else:
            return _err("calc must be n, power, or detect")

    except Exception as e:
        return _err(f"power error: {e}")


def _n_for_effect(d: float, power: float, alpha: float, ratio: float) -> int:
    """Per-group n for a given effect size.

    Exact solution: bisection on the noncentral-t distribution (the closed-form
    normal approximation under-states the n needed for small effects).
    Uses abs(d): Cohen's h/d is signed but sample size depends on magnitude
    (negative d previously returned n=2 — a dangerous under-powering bug).
    """
    import numpy as np
    from scipy import stats as sps

    d = abs(float(d))
    if d <= 0:
        return 2

    def achieved(n1: int) -> float:
        n2 = max(2, int(round(n1 * ratio)))
        df = n1 + n2 - 2
        nc = d * np.sqrt(1.0 / (1.0 / n1 + 1.0 / n2))
        crit = sps.t.ppf(1 - alpha / 2, df)
        return float(1 - sps.nct.cdf(crit, df, nc) + sps.nct.cdf(-crit, df, nc))

    lo, hi = 2, 4
    while achieved(hi) < power and hi < 200000:
        hi *= 2
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if achieved(mid) >= power:
            hi = mid
        else:
            lo = mid
    return max(int(hi), 2)


# ═══════════════════════════════════════════════════════════════════
# Tool: med_evidence — Evidence-based medicine calculations
# ═══════════════════════════════════════════════════════════════════

def med_evidence(
    tp: Optional[int] = None,
    tn: Optional[int] = None,
    fp: Optional[int] = None,
    fn: Optional[int] = None,
    nnt_calc: Optional[str] = None,
    cer: Optional[float] = None,
    eer: Optional[float] = None,
) -> str:
    """
    Calculate evidence-based medicine metrics from a 2x2 table or event rates.

    Diagnostic metrics (use tp, tn, fp, fn):
      tp=true positives, tn=true negatives, fp=false positives, fn=false negatives

    Treatment metrics (use nnt_calc, cer, eer):
      nnt_calc: "arr" or "nnt"
      cer:      control event rate (proportion, e.g. 0.15 for 15%)
      eer:      experimental event rate (proportion, e.g. 0.10 for 10%)
    """
    try:
        # --- Diagnostic test metrics ---
        if tp is not None and tn is not None and fp is not None and fn is not None:
            if min(tp, tn, fp, fn) < 0:
                return _err("cell counts cannot be negative")
            n_all = tp + tn + fp + fn
            if n_all <= 0:
                return _err("all four cells are zero — nothing to compute")
            sens = tp / (tp + fn) if (tp + fn) > 0 else 0
            spec = tn / (tn + fp) if (tn + fp) > 0 else 0
            ppv = tp / (tp + fp) if (tp + fp) > 0 else 0
            npv = tn / (tn + fn) if (tn + fn) > 0 else 0
            acc = (tp + tn) / n_all
            dor = (tp * tn) / (fp * fn) if fp > 0 and fn > 0 else float("inf")
            prev = (tp + fn) / n_all

            # CIs: Wilson for every proportion, Katz log CI for the likelihood
            # ratios. A point estimate alone is not reportable.
            sens_ci = tuple(round(v, 4) for v in _wilson_ci(tp, tp + fn))
            spec_ci = tuple(round(v, 4) for v in _wilson_ci(tn, tn + fp))
            ppv_ci = tuple(round(v, 4) for v in _wilson_ci(tp, tp + fp)) if (tp + fp) else None
            npv_ci = tuple(round(v, 4) for v in _wilson_ci(tn, tn + fn)) if (tn + fn) else None

            lr_plus = sens / (1 - spec) if spec < 1 else float("inf")
            lr_minus = (1 - sens) / spec if spec > 0 else float("inf")

            lrp_ci = None
            lrn_ci = None
            if fp > 0 and tn > 0 and tp > 0 and fn > 0:
                se_lrp = (1 / tp + 1 / fp - 1 / (tp + fn) - 1 / (tn + fp)) ** 0.5
                se_lrn = (1 / fn + 1 / tn - 1 / (tp + fn) - 1 / (tn + fp)) ** 0.5
                lrp_ci = (round(lr_plus * pow(2.718281828, -1.96 * se_lrp), 3),
                          round(lr_plus * pow(2.718281828, 1.96 * se_lrp), 3))
                lrn_ci = (round(lr_minus * pow(2.718281828, -1.96 * se_lrn), 3),
                          round(lr_minus * pow(2.718281828, 1.96 * se_lrn), 3))

            out: dict = {
                "m": "diagnostic",
                "sens": round(sens, 4), "sens_ci": sens_ci,
                "spec": round(spec, 4), "spec_ci": spec_ci,
                "ppv": round(ppv, 4), "ppv_ci": ppv_ci,
                "npv": round(npv, 4), "npv_ci": npv_ci,
                "acc": round(acc, 4), "prev": round(prev, 4),
                "lr+": round(lr_plus, 2) if lr_plus != float("inf") else "inf",
                "lr+_ci": lrp_ci,
                "lr-": round(lr_minus, 3) if lr_minus != float("inf") else "inf",
                "lr-_ci": lrn_ci,
                "dor": round(dor, 1) if dor != float("inf") else "inf",
                "n": n_all,
            }
            if (tp + fn) == 0 or (tn + fp) == 0:
                out["warn"] = "a disease group has zero members — sensitivity/specificity undefined for that arm"
            return _ok(out)

        # --- Treatment effect metrics ---
        if cer is not None and eer is not None:
            if not (0 <= cer <= 1 and 0 <= eer <= 1):
                return _err("cer and eer must be proportions between 0 and 1")
            arr = cer - eer  # Absolute Risk Reduction
            rrr = arr / cer if cer > 0 else 0  # Relative Risk Reduction
            rr = eer / cer if cer > 0 else 0  # Relative Risk
            or_val = (eer / (1 - eer)) / (cer / (1 - cer)) if cer < 1 and eer < 1 else 0

            # 95% CI for ARR
            out: dict = {
                "m": "treatment",
                "cer": cer, "eer": eer,
                "arr": round(arr, 4),
                "rrr": round(rrr, 4),
                "rr": round(rr, 3),
                "or": round(or_val, 2),
            }

            # ARR CI: Wald on the difference of independent proportions, with the
            # denominator implied by the rates when only rates are supplied.
            n_eff = 100  # documented assumption: percentages behave as n=100
            se = ((cer * (1 - cer) / n_eff) + (eer * (1 - eer) / n_eff)) ** 0.5
            out["arr_ci"] = (round(arr - 1.96 * se, 4), round(arr + 1.96 * se, 4))
            # RR CI (Katz log method) — also an assumed n of 100 per arm.
            if cer > 0 and eer > 0:
                se_lnrr = ((1 / (cer * n_eff) - 1 / n_eff) + (1 / (eer * n_eff) - 1 / n_eff)) ** 0.5
                out["rr_ci"] = (round(rr * pow(2.718281828, -1.96 * se_lnrr), 3),
                                round(rr * pow(2.718281828, 1.96 * se_lnrr), 3))
            out["n_assumed_per_arm"] = n_eff

            if arr > 0:
                nnt = int(1 / arr + 0.9999)  # round up — never under-power the NNT
                out["nnt"] = nnt
                out["nnt_ci"] = (
                    int(1 / out["arr_ci"][1] + 0.9999) if out["arr_ci"][1] > 0 else "inf",
                    int(1 / out["arr_ci"][0] + 0.9999) if out["arr_ci"][0] > 0 else "inf",
                )
            elif arr < 0:
                out["nnh"] = int(-1 / arr + 0.9999)
                out["harm"] = True
            else:
                out["note"] = "no difference in event rates (ARR = 0) — NNT undefined"
            return _ok(out)

        return _err("provide either (tp,tn,fp,fn) for diagnostic or (cer,eer) for treatment")

    except Exception as e:
        return _err(f"evidence error: {e}")


# ═══════════════════════════════════════════════════════════════════
# Schemas (token-efficient — minimal descriptions, short param names)
# ═══════════════════════════════════════════════════════════════════

MED_STATS_SCHEMA = {
    "name": "med_stats",
    "description": (
        "Statistical test for medical data. test: ttest|mw|wilcoxon|chisq|fisher|anova|kw|pearson|spearman. "
        "Continuous: pass a,b arrays. Categorical (chisq/fisher): pass categorical=[[a,b],[c,d]] (no a). "
        "paired=true for paired t-test/Wilcoxon. Returns JSON: {t:test,p:p-value,effect sizes + 95% CIs}"
    ),
    "parameters": {
        "type": "object",
        "properties": {
            "test": {
                "type": "string",
                "description": "Test: ttest, mw (Mann-Whitney), wilcoxon, chisq, fisher, anova, kw (Kruskal-Wallis), pearson, spearman",
            },
            "a": {
                "type": "array",
                "items": {"type": "number"},
                "description": "Group 1 data [n1,n2,...]; or pre-values for paired; or x for correlation. Omit for chisq/fisher.",
            },
            "b": {
                "type": "array",
                "items": {"type": "number"},
                "description": "Group 2 data; or post-values for paired; or y for correlation. Omit for 1-sample ttest.",
            },
            "categorical": {
                "type": "array",
                "items": {"type": "array", "items": {"type": "integer"}},
                "description": "2x2 table [[tp,fp],[fn,tn]] for chi-square/Fisher (use INSTEAD of a,b)",
            },
            "paired": {
                "type": "boolean",
                "description": "True for paired t-test or Wilcoxon signed-rank",
            },
        },
        "required": ["test"],
    },
}

MED_PUBMED_SCHEMA = {
    "name": "med_pubmed",
    "description": (
        "Search PubMed. q=query (PubMed syntax). fetch=true returns abstracts (uses 2 API calls). "
        "pmid=ID fetches single article. max_results: max 20. Returns compact JSON with short keys."
    ),
    "parameters": {
        "type": "object",
        "properties": {
            "q": {
                "type": "string",
                "description": "PubMed search: 'diabetes AND metformin[TIAB] AND 2024[DP]'. Ignored when pmid set.",
            },
            "max_results": {
                "type": "integer",
                "description": "Max articles (default 5, max 20)",
            },
            "fetch": {
                "type": "boolean",
                "description": "Fetch abstracts too (2 API calls); False returns PMIDs only",
            },
            "pmid": {
                "type": "string",
                "description": "Single PMID to fetch directly; overrides query",
            },
            "retstart": {
                "type": "integer",
                "description": "Pagination offset (default 0)",
            },
        },
        "required": [],
    },
}

MED_TRIAL_SCHEMA = {
    "name": "med_trial",
    "description": (
        "Search ClinicalTrials.gov. q=drug/condition. status filter: recruiting|active|completed. "
        "fmt=brief (NCT+title+phase+status) or full (adds conditions+locations). max_results: max 20."
    ),
    "parameters": {
        "type": "object",
        "properties": {
            "q": {
                "type": "string",
                "description": "Search term: 'pembrolizumab melanoma', 'Crohn's disease', etc.",
            },
            "max_results": {
                "type": "integer",
                "description": "Max studies (default 5, max 20)",
            },
            "status": {
                "type": "string",
                "description": "Filter: recruiting, active (=ACTIVE_NOT_RECRUITING), completed. Omit for all.",
            },
            "fmt": {
                "type": "string",
                "description": "brief=NCT+title+phase+status; full=+conditions+locations",
            },
        },
        "required": ["q"],
    },
}

MED_POWER_SCHEMA = {
    "name": "med_power",
    "description": (
        "Sample size/power for medical studies. calc=n (n needed)|power (power given n)|detect (detectable effect). "
        "effect=Cohen's d (0.2/0.5/0.8), or n per group when n is omitted. For proportions set p1+p2 (arcsin). "
        "n=per-group n for calc=power/detect. dropout=expected attrition (0.15) inflates n. ratio=n2/n1."
    ),
    "parameters": {
        "type": "object",
        "properties": {
            "calc": {
                "type": "string",
                "description": "n=sample size needed, power=power given n, detect=detectable effect",
            },
            "effect": {
                "type": "number",
                "description": "Cohen's d (0.2 small, 0.5 medium, 0.8 large). For power/detect, use n instead.",
            },
            "n": {
                "type": "integer",
                "description": "Per-group n for calc=power or detect",
            },
            "dropout": {
                "type": "number",
                "description": "Expected attrition fraction (0.15 = 15%) — inflates the required n",
            },
            "power": {
                "type": "number",
                "description": "Target power (default 0.80)",
            },
            "alpha": {
                "type": "number",
                "description": "Significance level (default 0.05, two-sided)",
            },
            "ratio": {
                "type": "number",
                "description": "n2/n1 ratio (default 1.0 for equal groups)",
            },
            "p1": {
                "type": "number",
                "description": "Proportion in group 1 (for proportion tests)",
            },
            "p2": {
                "type": "number",
                "description": "Proportion in group 2 (for proportion tests)",
            },
        },
        "required": ["calc"],
    },
}

MED_EVIDENCE_SCHEMA = {
    "name": "med_evidence",
    "description": (
        "EBM metrics from 2x2 table or event rates. "
        "Diagnostic: pass tp,tn,fp,fn → sens,specificity,PPV,NPV,LR+,LR-,DOR. "
        "Treatment: pass cer,eer (proportions) → ARR,RRR,RR,OR,NNT."
    ),
    "parameters": {
        "type": "object",
        "properties": {
            "tp": {"type": "integer", "description": "True positives"},
            "tn": {"type": "integer", "description": "True negatives"},
            "fp": {"type": "integer", "description": "False positives"},
            "fn": {"type": "integer", "description": "False negatives"},
            "cer": {
                "type": "number", "description": "Control Event Rate (e.g. 0.15 for 15%)",
            },
            "eer": {
                "type": "number", "description": "Experimental Event Rate (e.g. 0.10 for 10%)",
            },
        },
        "required": [],
    },
}


# ═══════════════════════════════════════════════════════════════════
# Registry
# ═══════════════════════════════════════════════════════════════════

from tools.registry import registry


def _check_stats() -> bool:
    return _has_scipy()


def _check_pubmed() -> bool:
    return _has_pubmed()


def _check_trial() -> bool:
    return True  # Always available (stdlib only)


def _check_power() -> bool:
    return _has_scipy()


def _check_evidence() -> bool:
    return True  # Pure calculation, no deps


registry.register(
    name="med_stats",
    toolset="medical",
    schema=MED_STATS_SCHEMA,
    handler=lambda args, **kw: med_stats(
        test=args.get("test", "ttest"),
        a=args.get("a", []),
        b=args.get("b"),
        categorical=args.get("categorical"),
        paired=args.get("paired", False),
    ),
    check_fn=_check_stats,
    emoji="📊",
)

registry.register(
    name="med_pubmed",
    toolset="medical",
    schema=MED_PUBMED_SCHEMA,
    handler=lambda args, **kw: med_pubmed(
        q=args.get("q", ""),
        max_results=args.get("max_results", 5),
        fetch=args.get("fetch", False),
        pmid=args.get("pmid"),
        retstart=args.get("retstart", 0),
    ),
    check_fn=_check_pubmed,
    emoji="🔬",
)

registry.register(
    name="med_trial",
    toolset="medical",
    schema=MED_TRIAL_SCHEMA,
    handler=lambda args, **kw: med_trial(
        q=args.get("q", ""),
        max_results=args.get("max_results", 5),
        status=args.get("status"),
        fmt=args.get("fmt", "brief"),
    ),
    check_fn=_check_trial,
    emoji="💊",
)

registry.register(
    name="med_power",
    toolset="medical",
    schema=MED_POWER_SCHEMA,
    handler=lambda args, **kw: med_power(
        calc=args.get("calc", "n"),
        effect=args.get("effect", 0.5),
        power=args.get("power", 0.80),
        alpha=args.get("alpha", 0.05),
        ratio=args.get("ratio", 1.0),
        p1=args.get("p1"),
        p2=args.get("p2"),
        n=args.get("n"),
        dropout=args.get("dropout", 0.0),
    ),
    check_fn=_check_power,
    emoji="📐",
)

registry.register(
    name="med_evidence",
    toolset="medical",
    schema=MED_EVIDENCE_SCHEMA,
    handler=lambda args, **kw: med_evidence(
        tp=args.get("tp"),
        tn=args.get("tn"),
        fp=args.get("fp"),
        fn=args.get("fn"),
        cer=args.get("cer"),
        eer=args.get("eer"),
    ),
    check_fn=_check_evidence,
    emoji="⚖️",
)
