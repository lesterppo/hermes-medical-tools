---
name: hermes-medical-tools
description: Native Hermes medical research tools — PSPP, jmv/R, PubMed, trials, stats, EBM.
version: 1.0.0
author: lesterppo
license: MIT
metadata:
  hermes:
    category: research
    tags: [medical, statistics, pubmed, clinical-trials, evidence-based-medicine, pspp, jamovi]
---

# Hermes Medical Tools Skill

## When to Use

Load when the user asks for:
- Statistical analysis of medical/clinical data
- PubMed literature search
- Clinical trial lookup
- Sample size or power calculation
- Evidence-based medicine metrics (NNT, sensitivity, specificity)
- Any medical research task requiring statistics

## How to Run

1. Install prerequisites (see README.md)
2. Copy tools to `~/.hermes/hermes-agent/tools/`
3. Add the `medical` toolset to `~/.hermes/hermes-agent/toolsets.py`
4. Run `hermes tools enable medical`
5. Start a new Hermes session

## Quick Reference

| Task | Tool | Example |
|------|------|---------|
| t-test | pspp, jmv, med_stats | `pspp(syntax="T-TEST GROUPS=g(1 2) /VARIABLES=s.", data=...)` |
| ANOVA | pspp, jmv, med_stats | `jmv(analysis="anovaOneW", data="grp,val\nA,10\nB,20\nC,30")` |
| Regression | pspp, jmv | `jmv(analysis="linReg", data="x,y\n1,2\n2,3\n3,5")` |
| Logistic regression | jmv | `jmv(analysis="logRegBin", data="y,x\n0,1\n0,2\n1,4\n1,5")` |
| Descriptives | pspp, jmv | `jmv(analysis="descriptives", data="a,b,c\n1,2,3\n4,5,6")` |
| Chi-square | pspp, jmv, med_stats | `jmv(analysis="contTables", data="tx,out\n1,1\n1,0\n0,1")` |
| Correlation | jmv | `jmv(analysis="corrMatrix", data="x,y,z\n1,2,3\n2,3,1")` |
| PubMed search | med_pubmed | `med_pubmed(q="metformin AND diabetes[TIAB]")` |
| Clinical trials | med_trial | `med_trial(q="melanoma", status="recruiting")` |
| Sample size | med_power | `med_power(calc="n", effect=0.5)` |
| NNT, sens/spec | med_evidence | `med_evidence(tp=85, tn=90, fp=10, fn=15)` |

## Pitfalls

1. **PSPP needs DATA LIST format** — not CSV. Use `DATA LIST FREE /var1 var2. BEGIN DATA ... END DATA.`
2. **jmv needs CSV with header** — first row must be column names
3. **med_stats ANOVA multi-group** — pass `a=[[g1],[g2],[g3]]` for 3+ groups
4. **Service gating** — tools only appear when backend is installed. Run `hermes doctor` to check.
5. **R packages** — jmv uses base R stats only. No heavy R package installs needed.
6. **PSPP survival analysis** — PSPP supports SURVIVAL command but syntax is complex

## Credits

- [GNU PSPP](https://www.gnu.org/software/pspp/) — free SPSS replacement
- [jamovi](https://github.com/jamovi/jamovi) — free, open-source statistics platform
- [SciPy](https://scipy.org/) — scientific computing
- [PubMed / NCBI](https://pubmed.ncbi.nlm.nih.gov/) — biomedical literature
- [ClinicalTrials.gov](https://clinicaltrials.gov/) — clinical trial registry
- [CEBM Oxford](https://www.cebm.ox.ac.uk/) — evidence-based medicine tools
