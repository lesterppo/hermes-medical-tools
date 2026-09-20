# Hermes Medical Tools

Native, token-efficient medical research tools for [Hermes Agent](https://github.com/NousResearch/hermes-agent).

7 tools covering PubMed search (retraction-aware), clinical trials, statistical
analysis (scipy, GNU PSPP, R/jamovi), exact sample-size/power calculation, and
evidence-based medicine metrics with confidence intervals.

## Install (plugin — survives `hermes update`)

```bash
git clone https://github.com/lesterppo/hermes-medical-tools.git
cd hermes-medical-tools
./install.sh                 # → ~/.hermes/plugins/hermes_medical_tools
./install.sh --check         # report backends only
./install.sh --uninstall
hermes plugins enable hermes_medical_tools   # if not auto-enabled
# restart Hermes / the gateway
```

**Why a plugin and not a copy into `hermes-agent/tools/`:** `hermes update`
resets everything inside the git tree — locally added `tools/*.py`, the
`toolsets.py` entry, `tools_config.py`. The old copy-in install therefore lost
all seven tools (and printed `Unknown toolsets: medical`) after every update.
`./install.sh --legacy` still performs the old copy for anyone who wants it, and
now says so explicitly.

The installer also detects other plugins that register the same tool names
(`med_pubmed`, `med_stats`, …) and tells you which one to disable: duplicate
registrations overwrite each other and the winner is load-order dependent.

Backends:

```bash
sudo apt install pspp          # pspp  — SPSS-syntax analyses
sudo apt install r-base        # jmv   — jamovi-compatible analyses
pip install scipy numpy        # med_stats, med_power
```

`med_pubmed`, `med_trial`, `med_evidence` need no backend (stdlib only).

## What changed in v2 (defects found by running real sessions)

| Symptom observed | Root cause | Fix |
|---|---|---|
| `med_stats(test="chisq", categorical=[[45,15],[10,30]], a=[])` — the call documented in this README — returned `a must be a non-empty list of numbers` | continuous-test validation ran before the categorical branch | categorical branch first; `a` optional; `required` reduced to `["test"]` |
| `med_trial(...)` always reported `total: 0` | ClinicalTrials.gov v2 omits `totalCount` unless `countTotal=true` | parameter added |
| n=1 per group returned `{"t":…,"s":NaN,"p":NaN}` with no warning | no minimum-n guard | explicit "need at least 2 per group" errors across tests |
| 2×2 tests returned a p-value only | no effect sizes | OR + RR (Haldane-Anscombe on zero cells), Cramér's V, expected-count warning |
| t-test returned means but no uncertainty | no CIs | t-based mean CIs, mean difference, Wilson intervals for sens/spec/PPV/NPV, Katz log CIs for LR± |
| `med_evidence(cer=0.15, eer=0.10)` returned a meaningless `ari: 0` | the harm column was emitted even for benefit | `nnh` only when harm; `nnt` rounds **up**; `arr_ci` / `nnt_ci` added |
| power/detect could only be computed by overloading `effect` as n | API gap | `n` parameter; exact **noncentral-t** power instead of the normal approximation |
| protocols under-powered by attrition | no inflation | `dropout` (e.g. `0.15`) inflates n via `ceil(n/(1-dropout))` |
| `med_pubmed(pmid=…)` gave title/abstract only | no citation identifiers, no safety signal | adds `doi`, `pt` (publication types) and `warn: RETRACTED / EXPRESSION_OF_CONCERN / CORRECTED` |

## Tools

### med_stats — statistical tests (scipy)

```python
med_stats(test="ttest", a=[120,125,130,128], b=[140,145,138,142])
# → t, p, n1/n2, means + 95% CIs, mean difference, Cohen's d

med_stats(test="anova", a=[[10,12,11],[20,22,21],[30,32,31]])   # + eta²
med_stats(test="kw", a=[[10,12],[20,22],[30,32]])               # + eta²(H), medians
med_stats(test="chisq", categorical=[[45,15],[10,30]])          # + OR, RR, Cramér's V
med_stats(test="fisher", categorical=[[1,9],[8,2]])
med_stats(test="mw", a=[1,2,3,4], b=[5,6,7,8])                  # + rank-biserial
med_stats(test="pearson", a=[1,2,3,4,5,6], b=[2,3,5,4,6,8])     # + Fisher-z CI, r²
med_stats(test="wilcoxon", a=[100,110,105], b=[95,98,92], paired=True)
```

Tests: `ttest`, `mw`, `wilcoxon`, `chisq`, `fisher`, `anova`, `kw`, `pearson`,
`spearman`. Effect sizes and 95% CIs always accompany the p-value — a p-value
alone is not reportable in a clinical manuscript.

### med_power — sample size, power, detectable effect

```python
med_power(calc="n", effect=0.5)                          # n per group for d=0.5, 80% power
med_power(calc="n", effect=0.5, dropout=0.15)            # inflate for 15% attrition
med_power(calc="n", p1=0.30, p2=0.20)                    # proportions (Cohen's h)
med_power(calc="power", effect=0.5, n=30)                # actual power at n=30/group
med_power(calc="detect", effect=0.5, n=20)               # smallest detectable d
```

Power and n come from the **noncentral t** distribution (normal approximations
overstate power at the small n typical of pilot work). `ratio` sets n2/n1.

Non-inferiority (`design="noninf"`, H0: diff ≤ −margin, one-sided alpha
default 0.025) uses stated normal approximations for both endpoints:

```python
med_power(calc="n", design="noninf", margin=0.5, sd=1.0)        # continuous: 63/group
med_power(calc="n", design="noninf", margin=0.1, p_ctrl=0.30)   # binary risk-diff: 330/group
med_power(calc="detect", design="noninf", sd=1.0, n=63)         # smallest winnable margin
```

Continuous: n1 = (1+1/r)·sd²·(z_a+z_b)²/(true_diff+margin)², power =
Φ((true_diff+margin)/SE − z_a), SE = sd·√(1/n1+1/n2). Binary (risk-difference
scale, pT = p_ctrl+true_diff): n1 =
(z_a+z_b)²·[pC(1−pC)+pT(1−pT)/r]/(true_diff+margin)² with the analogous power.
Wald unpooled variance under the alternative — the constrained
Farrington-Manning variance is NOT used. `true_diff` (T−C) defaults to 0.

### med_evidence — EBM metrics

```python
med_evidence(tp=85, tn=90, fp=10, fn=15)
# sens 0.85 [0.767, 0.907], spec 0.90 [0.826, 0.945], PPV, NPV,
# LR+ 8.5 [4.69, 15.4], LR- 0.167 [0.104, 0.267], DOR, prevalence

med_evidence(cer=0.15, eer=0.10)
# ARR 0.05, RRR 0.333, RR 0.667, OR 0.63, NNT 20, arr_ci, nnt_ci, rr_ci

med_evidence(cer=0.10, eer=0.15)     # → harm: NNH 20, harm: true
```

### med_pubmed — PubMed search / fetch

```python
med_pubmed(q="metformin AND diabetes[TIAB] AND 2024[DP]", max_results=5, fetch=True)
med_pubmed(pmid="9500320")
# → pmid, t, ab (structured labels preserved), yr, jr, doi, pt, warn: RETRACTED
```

Optional: `PUBMED_API_KEY` for higher rate limits, `PUBMED_EMAIL` for
identification.

### med_trial — ClinicalTrials.gov

```python
med_trial(q="inebilizumab IgG4", max_results=5, status="recruiting", fmt="full")
```

### pspp / jmv — SPSS syntax and jamovi-compatible analyses

```python
pspp(syntax="T-TEST GROUPS=group(1 2) /VARIABLES=score.",
     data="DATA LIST FREE /group score.\nBEGIN DATA\n1 85\n1 90\n2 72\n2 75\nEND DATA.")
pspp(syntax="CROSSTABS var1 BY var2 /STATISTICS=CHISQ.", data="...")

jmv(analysis="ttestIS",      data="group,score\n1,10\n1,12\n2,20\n2,22")
jmv(analysis="anovaOneW",    data="grp,val\nA,10\nA,12\nB,20\nB,22\nC,30\nC,32")
jmv(analysis="linReg",       data="x,y\n1,2\n2,3\n3,5\n4,4\n5,7")
jmv(analysis="logRegBin",    data="y,x\n0,1\n0,2\n0,3\n1,4\n1,5\n1,6")
jmv(analysis="descriptives", data="a,b,c\n1,2,3\n4,5,6")
jmv(analysis="corrMatrix",   data="x,y,z\n1,2,3\n2,3,1\n4,5,6")
jmv(analysis="contTables",   data="tx,out\n1,1\n1,0\n0,1\n0,0")
jmv(analysis="wilcoxon",     data="grp,val\n1,10\n1,12\n2,20\n2,22")
jmv(analysis="ttestPS",      data="pre,post\n100,95\n110,98\n105,92")
```

`pspp` needs `DATA LIST FREE` syntax (not CSV). `jmv` needs a CSV header row and
uses base R stats only — no heavyweight R packages.

## Prerequisites by tool

| Tool | Requirement | Install |
|------|-------------|---------|
| `pspp` | `pspp` binary | `sudo apt install pspp` |
| `jmv` | `Rscript` binary | `sudo apt install r-base` |
| `med_stats`, `med_power` | scipy, numpy | `pip install scipy numpy` |
| `med_pubmed`, `med_trial`, `med_evidence` | none (stdlib) | — |

## Credits & reference sources

| Tool | Backend | Reference |
|------|---------|-----------|
| `jmv` | R base stats (jamovi-compatible API) | [jamovi](https://github.com/jamovi/jamovi) |
| `pspp` | GNU PSPP | [GNU PSPP](https://www.gnu.org/software/pspp/) |
| `med_pubmed` | NCBI E-utilities | [PubMed](https://pubmed.ncbi.nlm.nih.gov/) |
| `med_trial` | ClinicalTrials.gov v2 API | [ClinicalTrials.gov](https://clinicaltrials.gov/) |
| `med_stats`, `med_power` | scipy.stats | [SciPy](https://scipy.org/) |
| `med_evidence` | pure Python | [CEBM Oxford](https://www.cebm.ox.ac.uk/) |

## Token efficiency

* short parameter names (`q`, `a`, `b`, `tp`, `tn`, `n`, `calc`)
* compact JSON with 1–3 character keys (`t`, `p`, `d`, `n`, `r`, `ci`, `eta2`)
* service gating via `check_fn`: an unavailable backend costs zero schema tokens
* no MCP layer — the tools call the APIs directly

## Companion project

[`lesterppo/med-search-cli`](https://github.com/lesterppo/med-search-cli) is the
retrieval-side sibling: twin-track PubMed + EuropePMC search, citation chasing,
MeSH lookup, reference export and standing-query surveillance. Use it for the
literature leg of a project and these tools for the analysis leg.

## License

MIT — same as Hermes Agent.
