---
name: hermes-medical-tools
description: Statistics tools for Hermes — PubMed, trials, EBM, power, PSPP, R.
version: 2.0.0
author: lesterppo
license: MIT
metadata:
  hermes:
    category: research
    tags: [medical, statistics, pubmed, clinical-trials, evidence-based-medicine, pspp, jamovi, sample-size]
---

# Hermes Medical Tools Skill

## When to Use

Load when the user asks for:
- Statistical analysis of medical/clinical data
- PubMed literature search (with retraction status)
- Clinical trial lookup
- Sample size / power / detectable-effect calculation
- Evidence-based medicine metrics (NNT, NNH, sensitivity, specificity, LRs)
- Any medical research task requiring statistics

## Install (plugin — survives `hermes update`)

```bash
./install.sh                 # plugin install → ~/.hermes/plugins/hermes_medical_tools
./install.sh --check         # backend availability only
./install.sh --uninstall
```

Do **not** use `./install.sh --legacy`: copying into
`~/.hermes/hermes-agent/tools/` is wiped by the next `hermes update` and leaves
`Unknown toolsets: medical`. The installer also warns when another plugin
already registers the same tool names (`med_pubmed`, `med_stats`, …), because
duplicate registrations overwrite each other.

## Quick Reference

| Task | Tool | Example |
|------|------|---------|
| t-test (+ CI, Cohen's d) | med_stats | `med_stats(test="ttest", a=[120,125,130], b=[140,145,138])` |
| ANOVA (+ eta²) | med_stats, jmv, pspp | `med_stats(test="anova", a=[[10,12],[20,22],[30,32]])` |
| Chi-square / Fisher (+ OR, RR, Cramér's V) | med_stats | `med_stats(test="chisq", categorical=[[45,15],[10,30]])` |
| Mann-Whitney / Kruskal-Wallis (+ effect size) | med_stats | `med_stats(test="mw", a=[1,2,3], b=[5,6,7])` |
| Correlation (+ Fisher-z CI) | med_stats | `med_stats(test="pearson", a=[1,2,3,4,5], b=[2,3,5,4,6])` |
| Regression / logistic | jmv, pspp | `jmv(analysis="linReg", data="x,y\n1,2\n2,3\n3,5")` |
| SPSS syntax | pspp | `pspp(syntax="T-TEST GROUPS=g(1 2) /VARIABLES=s.", data=...)` |
| PubMed search / fetch | med_pubmed | `med_pubmed(q="metformin diabetes[TIAB]", fetch=True)` |
| Retraction check | med_pubmed | `med_pubmed(pmid="9500320")` → `warn: RETRACTED` |
| Clinical trials | med_trial | `med_trial(q="inebilizumab", status="recruiting")` |
| Sample size | med_power | `med_power(calc="n", effect=0.5, dropout=0.15)` |
| Power for a given n | med_power | `med_power(calc="power", effect=0.5, n=30)` |
| Detectable effect | med_power | `med_power(calc="detect", effect=0.5, n=20)` |
| NNT / NNH with CI | med_evidence | `med_evidence(cer=0.15, eer=0.10)` |
| Sens/spec/LR with CI | med_evidence | `med_evidence(tp=85, tn=90, fp=10, fn=15)` |

## Reportable output

Every test returns an **effect size and/or 95% CI next to the p-value**
(Cohen's d, eta², rank-biserial, Cramér's V, OR/RR, Wilson intervals for
proportions). A p-value alone is not reportable in a clinical manuscript —
use these fields directly in the results table.

## Pitfalls

1. **`chisq` / `fisher` take `categorical` only** — do not pass `a=[]`;
   `a` is optional for those tests (this used to fail).
2. **NNT rounds up** (`int(1/ARR + 0.9999)`) so the number never under-powers
   the estimate; `arr_ci`/`nnt_ci` assume n=100 per arm when only rates are
   given (`n_assumed_per_arm` echoes the assumption).
3. **`calc="power"` needs `n`** (per group). Passing n through `effect` still
   works for backward compatibility but is deprecated.
4. **Power is computed from the noncentral t**, not the normal approximation —
   do not compare its numbers against normal-approximation calculators for
   small n; this tool is the accurate one.
5. **`dropout` inflates n before rounding** (`ceil(n/(1-dropout))`); report the
   inflated figure in a protocol.
6. **PSPP needs DATA LIST format**, not CSV: `DATA LIST FREE /a b. BEGIN DATA … END DATA.`
7. **jmv needs CSV with a header row.**
8. **Multi-group ANOVA/KW**: pass `a=[[g1],[g2],[g3]]`.
9. **Service gating**: `pspp` needs the GNU PSPP binary, `jmv` needs `Rscript`,
   `med_stats`/`med_power` need scipy. Unavailable tools simply are not exposed;
   `./install.sh --check` tells you what is missing.
10. **`med_trial` legacy code without `countTotal` reported `total: 0`** — fixed;
    if you port this code elsewhere, keep `countTotal=true`.

## Credits

- [GNU PSPP](https://www.gnu.org/software/pspp/) — free SPSS replacement
- [jamovi](https://github.com/jamovi/jamovi) — free, open-source statistics platform
- [SciPy](https://scipy.org/) — scientific computing
- [PubMed / NCBI](https://pubmed.ncbi.nlm.nih.gov/) — biomedical literature
- [ClinicalTrials.gov](https://clinicaltrials.gov/) — clinical trial registry
- [CEBM Oxford](https://www.cebm.ox.ac.uk/) — evidence-based medicine tools
