# Hermes Medical Tools

Native, token-efficient medical research tools for [Hermes Agent](https://github.com/NousResearch/hermes-agent).

7 tools covering PubMed search, clinical trials, statistical analysis (via GNU PSPP, R, and scipy), sample size calculation, and evidence-based medicine metrics. Total schema: ~1,645 tokens.

## Credits & Reference Sources

These tools build on the following excellent projects:

| Tool | Backend | Reference |
|------|---------|-----------|
| `jmv` | R base stats (jamovi-compatible API) | [jamovi](https://github.com/jamovi/jamovi) — free, open-source statistics platform |
| `pspp` | GNU PSPP | [GNU PSPP](https://www.gnu.org/software/pspp/) — free replacement for IBM SPSS |
| `med_pubmed` | NCBI E-utilities | [PubMed](https://pubmed.ncbi.nlm.nih.gov/) — National Library of Medicine |
| `med_trial` | ClinicalTrials.gov v2 API | [ClinicalTrials.gov](https://clinicaltrials.gov/) — US National Library of Medicine |
| `med_stats` | scipy.stats | [SciPy](https://scipy.org/) — scientific computing for Python |
| `med_power` | scipy.stats | [SciPy](https://scipy.org/) |
| `med_evidence` | Pure Python | Evidence-Based Medicine formulas from [CEBM Oxford](https://www.cebm.ox.ac.uk/) |

## Quick Install

```bash
# 1. Clone into your Hermes tools directory
git clone https://github.com/lesterppo/hermes-medical-tools.git /tmp/hermes-medical-tools
cp /tmp/hermes-medical-tools/tools/*.py ~/.hermes/hermes-agent/tools/

# 2. Install backends (choose based on which tools you need)
# PSPP statistics (135 SPSS commands):
sudo apt install pspp          # Debian/Ubuntu
# OR: brew install pspp        # macOS

# R statistics (jamovi-compatible, 9 analyses):
sudo apt install r-base        # Debian/Ubuntu
# OR: brew install r           # macOS

# Python statistics (med_stats, med_power):
pip install scipy numpy        # or: uv pip install scipy numpy

# 3. Add the medical toolset
# Edit ~/.hermes/hermes-agent/toolsets.py
# Find the TOOLSETS dict and add:
#   "medical": {
#       "description": "Medical research — jmv/R, PSPP, PubMed, clinical trials, EBM",
#       "tools": ["jmv", "pspp", "med_pubmed", "med_trial", "med_stats", "med_power", "med_evidence"],
#       "includes": [],
#   },

# Or use the included patch:
python3 -c "
import sys
sys.path.insert(0, '.')
# The install script handles this automatically
"

# 4. Enable and restart
hermes tools enable medical
# Start a new Hermes session
```

## Tools

### pspp — GNU PSPP (SPSS-compatible, 135 commands)
```python
pspp(syntax="DESCRIPTIVES x y /STATISTICS=MEAN STDDEV MIN MAX.",
     data="DATA LIST FREE /x y.\nBEGIN DATA\n1 2\n3 4\nEND DATA.")

pspp(syntax="T-TEST GROUPS=group(1 2) /VARIABLES=score.",
     data="DATA LIST FREE /group score.\nBEGIN DATA\n1 85\n1 90\n2 72\n2 75\nEND DATA.")

pspp(syntax="REGRESSION /DEPENDENT=y /METHOD=ENTER x1 x2.", data="...")
pspp(syntax="CROSSTABS var1 BY var2 /STATISTICS=CHISQ.", data="...")
```
Backend: `pspp` binary. 135 commands: T-TEST, ANOVA, REGRESSION, CROSSTABS, FREQUENCIES, FACTOR, RELIABILITY, ROC, NPAR TESTS, GLM, EXAMINE, CORRELATIONS, LOGISTIC REGRESSION, ONEWAY, MEANS, RANK, GRAPH, SURVIVAL...

### jmv — R statistical analysis (jamovi-compatible)
```python
jmv(analysis="ttestIS", data="group,score\n1,10\n1,12\n2,20\n2,22")
jmv(analysis="anovaOneW", data="grp,val\nA,10\nA,12\nB,20\nB,22\nC,30\nC,32")
jmv(analysis="linReg", data="x,y\n1,2\n2,3\n3,5\n4,4\n5,7")
jmv(analysis="logRegBin", data="y,x\n0,1\n0,2\n0,3\n1,4\n1,5\n1,6")
jmv(analysis="descriptives", data="a,b,c\n1,2,3\n4,5,6")
jmv(analysis="corrMatrix", data="x,y,z\n1,2,3\n2,3,1\n4,5,6")
jmv(analysis="contTables", data="tx,out\n1,1\n1,0\n0,1\n0,0")
jmv(analysis="wilcoxon", data="grp,val\n1,10\n1,12\n2,20\n2,22")
jmv(analysis="ttestPS", data="pre,post\n100,95\n110,98\n105,92")
```
9 analyses: ttestIS, ttestPS, anovaOneW, linReg, logRegBin, descriptives, corrMatrix, contTables, wilcoxon. Backend: Rscript + base R stats. No compilation-heavy R packages needed.

### med_pubmed — PubMed search
```python
med_pubmed(q="metformin AND diabetes[TIAB] AND 2024[DP]", max_results=5, fetch=True)
med_pubmed(pmid="32897388")  # fetch single article by PMID
```
Optional: set `PUBMED_API_KEY` for higher rate limits, `PUBMED_EMAIL` for identification.

### med_trial — ClinicalTrials.gov search
```python
med_trial(q="pembrolizumab melanoma", max_results=5, status="recruiting")
```

### med_stats — Quick statistical tests (scipy)
```python
med_stats(test="ttest", a=[120,125,130], b=[140,145,138])
med_stats(test="anova", a=[[10,12],[20,22],[30,32]])  # multi-group
med_stats(test="chisq", categorical=[[45,15],[10,30]], a=[])
```
Tests: ttest, mw, wilcoxon, chisq, fisher, anova, kw, pearson, spearman.

### med_power — Sample size / power
```python
med_power(calc="n", effect=0.5)           # n for d=0.5 at 80% power
med_power(calc="n", p1=0.3, p2=0.2)       # proportions
```

### med_evidence — EBM metrics
```python
med_evidence(tp=85, tn=90, fp=10, fn=15)  # sensitivity, specificity, PPV, NPV, LR+, LR-
med_evidence(cer=0.15, eer=0.10)          # ARR=0.05, NNT=20, RR=0.667
```

## Manual Toolset Setup

If the automated install doesn't work, add this to `~/.hermes/hermes-agent/toolsets.py`:

```python
# In the TOOLSETS dict, add:
"medical": {
    "description": "Medical research — jmv/R, PSPP, PubMed, clinical trials, EBM",
    "tools": ["jmv", "pspp", "med_pubmed", "med_trial", "med_stats", "med_power", "med_evidence"],
    "includes": [],
},
```

Then enable and restart:
```bash
hermes tools enable medical
# Start a new Hermes session
```

## Token Efficiency

- Total schema: 6,583 chars (~1,645 tokens) for all 7 tools
- Average output: ~390 bytes per call
- Short parameter names throughout (`q`, `a`, `b`, `tp`, `tn`)
- Compact JSON output with 1-2 char keys (`t`, `p`, `d`, `n`, `r`)
- Service-gated: tools only appear when their backend is installed

## Prerequisites by Tool

| Tool | Requirement | Install |
|------|-------------|---------|
| `pspp` | `pspp` binary | `sudo apt install pspp` |
| `jmv` | `Rscript` binary | `sudo apt install r-base` |
| `med_stats` | scipy, numpy | `pip install scipy numpy` |
| `med_power` | scipy, numpy | `pip install scipy numpy` |
| `med_pubmed` | None (stdlib) | — |
| `med_trial` | None (stdlib) | — |
| `med_evidence` | None (stdlib) | — |

## License

MIT — same as Hermes Agent.
