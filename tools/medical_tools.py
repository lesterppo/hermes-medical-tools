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

def _ok(result: dict) -> str:
    """Wrap result dict as compact JSON. No newlines, minimal whitespace."""
    return json.dumps(result, ensure_ascii=False, separators=(",", ":"))


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

def med_stats(
    test: str,
    a: List[float],
    b: Optional[List[float]] = None,
    categorical: Optional[List[List[int]]] = None,
    paired: bool = False,
) -> str:
    """
    Run a statistical test, returning compact results.

    test:  "ttest" | "mw" | "wilcoxon" | "chisq" | "fisher" | "anova" | "kw" | "pearson" | "spearman"
    a:     first group data (numbers) or paired pre-values
    b:     second group data (numbers) or paired post-values; omit for 1-sample
    categorical: 2x2 table [[a,b],[c,d]] for chi-square/fisher; use INSTEAD of a/b
    paired: True for paired t-test or Wilcoxon signed-rank
    """
    try:
        import numpy as np
        from scipy import stats as sps

        # --- Input validation ---
        if a is None or (isinstance(a, list) and len(a) == 0):
            return _err("a must be a non-empty list of numbers")

        # --- Categorical tests ---
        if categorical is not None:
            if not (isinstance(categorical, list) and len(categorical) == 2
                    and len(categorical[0]) == 2 and len(categorical[1]) == 2):
                return _err("categorical must be [[a,b],[c,d]]")
            tbl = np.array(categorical)
            if test == "fisher":
                _, p = sps.fisher_exact(tbl)
                return _ok({"t": "fisher", "p": round(float(p), 6)})
            elif test == "chisq":
                chi2, p, dof, expected = sps.chi2_contingency(tbl)
                return _ok({
                    "t": "chisq", "x2": round(float(chi2), 3),
                    "df": dof, "p": round(float(p), 6),
                    "n": int(tbl.sum()),
                })
            return _err(f"test='{test}' not for categorical; use chisq or fisher")

        # --- Continuous tests ---
        x = np.array(a, dtype=float)
        y = np.array(b, dtype=float) if b is not None else None

        if test == "ttest":
            if paired:
                if y is None:
                    return _err("paired ttest needs both a and b")
                stat, p = sps.ttest_rel(x, y)
            else:
                stat, p = sps.ttest_ind(x, y) if y is not None else sps.ttest_1samp(x, 0)
            return _ok({
                "t": "ttest", "s": round(float(stat), 3),
                "p": round(float(p), 6), "n1": len(x),
                "n2": len(y) if y is not None else 0,
                "m1": round(float(x.mean()), 2),
                "m2": round(float(y.mean()), 2) if y is not None else 0,
                "d": round(_cohens_d(x, y, paired), 2),
            })

        elif test == "mw":
            if y is None:
                return _err("Mann-Whitney needs two groups (a and b)")
            stat, p = sps.mannwhitneyu(x, y, alternative="two-sided")
            return _ok({
                "t": "mw", "u": float(stat), "p": round(float(p), 6),
                "n1": len(x), "n2": len(y),
                "md1": round(float(np.median(x)), 2),
                "md2": round(float(np.median(y)), 2),
            })

        elif test == "wilcoxon":
            if y is None:
                return _err("Wilcoxon needs a and b")
            stat, p = sps.wilcoxon(x, y)
            return _ok({
                "t": "wilcoxon", "s": float(stat), "p": round(float(p), 6),
                "n": len(x),
            })

        elif test == "anova":
            # Multi-group: a = [[g1...], [g2...], ...]
            if isinstance(a, list) and len(a) > 0 and isinstance(a[0], list):
                groups = [np.array(g, dtype=float) for g in a if len(g) > 0]
                if len(groups) < 2:
                    return _err("ANOVA needs at least 2 non-empty groups")
                stat, p = sps.f_oneway(*groups)
                n_total = sum(len(g) for g in groups)
                return _ok({
                    "t": "anova", "f": round(float(stat), 3),
                    "df1": len(groups) - 1, "df2": n_total - len(groups),
                    "p": round(float(p), 6), "k": len(groups),
                    "means": [round(float(g.mean()), 2) for g in groups],
                })
            # Flat a + b for 2-group
            if y is None or len(y) == 0:
                return _err("ANOVA needs at least 2 groups — pass a=[[g1],[g2],...] or a=g1,b=g2")
            stat, p = sps.f_oneway(x, y)
            return _ok({
                "t": "anova", "f": round(float(stat), 3),
                "df1": 1, "df2": len(x) + len(y) - 2,
                "p": round(float(p), 6),
            })

        elif test == "kw":
            # Multi-group: a = [[g1...], [g2...], ...]
            if isinstance(a, list) and len(a) > 0 and isinstance(a[0], list):
                groups = [np.array(g, dtype=float) for g in a if len(g) > 0]
                if len(groups) < 2:
                    return _err("Kruskal-Wallis needs at least 2 non-empty groups")
                stat, p = sps.kruskal(*groups)
                return _ok({
                    "t": "kw", "h": round(float(stat), 3),
                    "p": round(float(p), 6), "k": len(groups),
                    "n": sum(len(g) for g in groups),
                })
            # Flat a + b for 2-group
            if y is None or len(y) == 0:
                return _err("Kruskal-Wallis needs at least 2 groups — pass a=[[g1],[g2],...] or a=g1,b=g2")
            stat, p = sps.kruskal(x, y)
            return _ok({
                "t": "kw", "h": round(float(stat), 3),
                "p": round(float(p), 6), "n": len(x) + len(y),
            })

        elif test == "pearson":
            if y is None:
                return _err("Pearson needs x and y")
            r, p = sps.pearsonr(x, y)
            return _ok({
                "t": "pearson", "r": round(float(r), 4),
                "p": round(float(p), 6), "n": len(x),
            })

        elif test == "spearman":
            if y is None:
                return _err("Spearman needs x and y")
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
            abstract_el = article.find(".//AbstractText")
            title = title_el.text if title_el is not None and title_el.text else ""
            abstract = abstract_el.text if abstract_el is not None and abstract_el.text else ""

            return _ok({
                "pmid": pmid, "t": title[:400],
                "ab": abstract[:600] if abstract else "",
            })

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
                title = title_el.text[:300] if title_el is not None and title_el.text else ""
                ab_el = art.find(".//AbstractText")
                abstract = ab_el.text[:400] if ab_el is not None and ab_el.text else ""
                yr_el = art.find(".//PubDate/Year")
                yr = yr_el.text if yr_el is not None else ""
                articles.append({
                    "pmid": pid, "t": title,
                    "ab": abstract, "yr": yr,
                })
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
        }
        if status:
            params["filter.overallStatus"] = status.upper()

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
) -> str:
    """
    Calculate sample size or power for common medical study designs.

    calc:   "n" (sample size needed) | "power" (actual power given n) | "detect" (detectable effect)
    effect: Cohen's d (0.2 small, 0.5 medium, 0.8 large). For calc="detect" with proportions, this is n per group.
    power:  target power (default 0.80)
    alpha:  significance level (default 0.05, two-sided)
    ratio:  n2/n1 ratio (default 1 for equal groups)
    p1:     proportion in group 1 (for proportion tests)
    p2:     proportion in group 2 (for proportion tests). If p1+p2 given, uses arcsin method.
    """
    try:
        from scipy import stats as sps

        # For proportions
        if p1 is not None and p2 is not None:
            import numpy as np
            # Cohen's h (arcsine transformation)
            h = 2 * (np.arcsin(np.sqrt(p1)) - np.arcsin(np.sqrt(p2)))
            if calc == "n":
                n1 = _n_for_effect(h, power, alpha, ratio)
                return _ok({
                    "t": "proportions", "calc": "n",
                    "h": round(float(h), 4), "p1": p1, "p2": p2,
                    "n1": n1, "n2": round(n1 * ratio),
                    "total": n1 + round(n1 * ratio),
                    "alpha": alpha, "power": power,
                })
            elif calc == "power":
                if effect <= 0:
                    return _err("effect (n per group) required for power calc")
                n1 = int(effect)
                z_alpha = sps.norm.ppf(1 - alpha / 2)
                z_power = abs(h) * np.sqrt(n1 / (1 + 1 / ratio)) - z_alpha
                actual_power = float(sps.norm.cdf(z_power))
                return _ok({"t": "proportions", "calc": "power", "p": round(actual_power, 4)})

        # For continuous (Cohen's d)
        if calc == "n":
            n1 = _n_for_effect(effect, power, alpha, ratio)
            return _ok({
                "t": "continuous", "calc": "n",
                "d": effect, "n1": n1, "n2": round(n1 * ratio),
                "total": n1 + round(n1 * ratio),
                "alpha": alpha, "power": power,
            })
        elif calc == "power":
            if effect <= 0:
                return _err("effect must be sample size (n per group) for power calc")
            n1 = int(effect)
            import numpy as np
            z_alpha = sps.norm.ppf(1 - alpha / 2)
            # Non-centrality
            nc = effect * np.sqrt(n1 / (1 + 1/ratio))
            df = n1 + int(n1 * ratio) - 2
            # Use normal approximation for t-test power
            z_power = nc - z_alpha
            actual_power = float(sps.norm.cdf(z_power))
            return _ok({
                "t": "continuous", "calc": "power",
                "n1": n1, "p": round(actual_power, 4),
            })
        elif calc == "detect":
            if effect <= 0:
                return _err("effect must be n per group for detectable effect calc")
            n1 = int(effect)
            import numpy as np
            z_alpha = sps.norm.ppf(1 - alpha / 2)
            z_beta = sps.norm.ppf(power)
            d = (z_alpha + z_beta) * np.sqrt((1 + 1/ratio) / n1)
            return _ok({
                "t": "continuous", "calc": "detect",
                "n1": n1, "d": round(float(d), 4),
                "alpha": alpha, "power": power,
            })
        else:
            return _err("calc must be n, power, or detect")

    except Exception as e:
        return _err(f"power error: {e}")


def _n_for_effect(d: float, power: float, alpha: float, ratio: float) -> int:
    """Sample size per group for given effect size."""
    import numpy as np
    from scipy import stats as sps
    z_alpha = sps.norm.ppf(1 - alpha / 2)
    z_beta = sps.norm.ppf(power)
    n1 = int(np.ceil((z_alpha + z_beta)**2 * (1 + 1/ratio) / d**2))
    return max(n1, 2)


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
            sens = tp / (tp + fn) if (tp + fn) > 0 else 0
            spec = tn / (tn + fp) if (tn + fp) > 0 else 0
            ppv = tp / (tp + fp) if (tp + fp) > 0 else 0
            npv = tn / (tn + fn) if (tn + fn) > 0 else 0
            acc = (tp + tn) / (tp + tn + fp + fn)
            lr_plus = sens / (1 - spec) if spec < 1 else float("inf")
            lr_minus = (1 - sens) / spec if spec > 0 else float("inf")
            dor = (tp * tn) / (fp * fn) if fp > 0 and fn > 0 else float("inf")
            prev = (tp + fn) / (tp + tn + fp + fn)

            return _ok({
                "m": "diagnostic",
                "sens": round(sens, 4), "spec": round(spec, 4),
                "ppv": round(ppv, 4), "npv": round(npv, 4),
                "acc": round(acc, 4), "prev": round(prev, 4),
                "lr+": round(lr_plus, 2) if lr_plus != float("inf") else "inf",
                "lr-": round(lr_minus, 3) if lr_minus != float("inf") else "inf",
                "dor": round(dor, 1) if dor != float("inf") else "inf",
                "n": tp + tn + fp + fn,
            })

        # --- Treatment effect metrics ---
        if cer is not None and eer is not None:
            arr = cer - eer  # Absolute Risk Reduction
            rrr = arr / cer if cer > 0 else 0  # Relative Risk Reduction
            rr = eer / cer if cer > 0 else 0  # Relative Risk
            or_val = (eer / (1 - eer)) / (cer / (1 - cer)) if cer < 1 and eer < 1 else 0
            nnt = int(1 / arr) if arr > 0 else float("inf")

            return _ok({
                "m": "treatment",
                "cer": cer, "eer": eer,
                "arr": round(arr, 4), "rrr": round(rrr, 4),
                "rr": round(rr, 3), "or": round(or_val, 2),
                "nnt": nnt if nnt != float("inf") else "inf",
                "ari": round(-arr, 4) if arr < 0 else 0,  # Absolute Risk Increase (harm)
            })

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
        "Continuous: pass a,b arrays. Categorical (chisq/fisher): pass categorical=[[a,b],[c,d]]. "
        "paired=true for paired t-test/Wilcoxon. Returns JSON: {t:test,p:p-value,...}"
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
                "description": "Group 1 data [n1,n2,...]; or pre-values for paired; or x for correlation",
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
        "required": ["test", "a"],
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
                "description": "Filter: recruiting, active, completed. Omit for all.",
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
        "effect=Cohen's d (0.2/0.5/0.8). For proportions, set p1+p2 (uses arcsin). "
        "For power/detect, effect=n per group. ratio=n2/n1 (default 1)."
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
                "description": "Cohen's d (0.2 small, 0.5 medium, 0.8 large). For power/detect calc, this is n per group.",
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
        "required": ["calc", "effect"],
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
