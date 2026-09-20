#!/usr/bin/env python3
"""Offline regression tests for hermes-medical-tools bug fixes.

Run:  python3 tests/test_regression.py   (needs scipy+numpy; no network, no pspp/R)
Exit 0 = all green. Live API tests (PubMed/ClinicalTrials.gov) are NOT here —
they were verified manually against the live endpoints (see task report).
"""
import importlib.util
import json
import sys
import types
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

# Stub tools.registry (only present inside a Hermes install) so the tool
# modules' module-level registry.register(...) calls succeed standalone.
_reg = types.ModuleType("tools.registry")


class _Registry:
    def __init__(self):
        self.e = {}

    def register(self, name, toolset=None, schema=None, handler=None,
                 check_fn=None, emoji=None, **kw):
        from types import SimpleNamespace
        self.e[name] = SimpleNamespace(
            name=name, toolset=toolset, schema=schema, handler=handler,
            check_fn=check_fn, emoji=emoji, is_async=False)

    def get_entry(self, n):
        return self.e.get(n)


_reg.registry = _Registry()
_pkg = types.ModuleType("tools")
_pkg.__path__ = []
sys.modules["tools"] = _pkg
sys.modules["tools.registry"] = _reg


def _load(path, name):
    spec = importlib.util.spec_from_file_location(name, str(path))
    m = importlib.util.module_from_spec(spec)
    sys.modules[name] = m
    spec.loader.exec_module(m)
    return m


mt = _load(ROOT / "tools" / "medical_tools.py", "reg_mt")
pspp = _load(ROOT / "tools" / "pspp_tool.py", "reg_pspp")
jmv = _load(ROOT / "tools" / "jmv_tool.py", "reg_jmv")

G1 = [12.5, 14.1, 13.8, 15.2, 14.9, 13.1, 14.4, 15.0]
G2 = [10.2, 11.0, 9.8, 10.5, 11.3, 10.9, 10.1, 11.5]
PASS = FAIL = 0


def strict(label, s):
    """json.loads(strict=True) + reject NaN/Infinity literals."""
    global PASS, FAIL
    try:
        d = json.loads(s, strict=True)
        json.dumps(d, allow_nan=False)
        assert "NaN" not in s and "Infinity" not in s, "non-strict literal"
        PASS += 1
        return d
    except Exception as e:  # noqa: BLE001
        FAIL += 1
        print(f"FAIL {label}: {e} :: {s[:200]}")
        return None


def check(label, cond, detail=""):
    global PASS, FAIL
    if cond:
        PASS += 1
    else:
        FAIL += 1
        print(f"FAIL {label} {detail}")


# 1. ragged multi-group anova/kw (used to crash in np.array conversion)
d = strict("anova_ragged", mt.med_stats(test="anova", a=[G1, G2, [20.1, 21.3, 19.8, 20.5, 22.0]]))
check("anova_ragged_vals", d and d.get("t") == "anova" and d["k"] == 3, d)
d = strict("kw_ragged", mt.med_stats(test="kw", a=[G1, G2, [20.1, 21.3, 19.8, 20.5, 22.0]]))
check("kw_ragged_vals", d and d.get("t") == "kw" and d["k"] == 3, d)
from scipy import stats as _sps  # noqa: E402
check("anova_matches_scipy",
      abs(strict("x", mt.med_stats(test="anova", a=[G1, G2]))["f"]
          - round(float(_sps.f_oneway(G1, G2)[0]), 3)) < 1e-9)

# 2. rank-biserial sign: positive = group a higher (as documented)
check("rb_pos", strict("x", mt.med_stats(test="mw", a=G1, b=G2))["rb"] == 1.0)
check("rb_neg", strict("x", mt.med_stats(test="mw", a=G2, b=G1))["rb"] == -1.0)

# 3. 1-sample constant data warns on undefined effect size
d = strict("const_1samp", mt.med_stats(test="ttest", a=[5, 5, 5, 5]))
check("const_1samp_warn", d.get("d") is None and "warn" in d, d)

# 4. effect=0 errors; proportions support detect
check("power_d0", "e" in strict("x", mt.med_power(calc="n", effect=0.0)))
check("power_p12", "e" in strict("x", mt.med_power(calc="n", p1=0.2, p2=0.2)))
d = strict("prop_detect", mt.med_power(calc="detect", n=120, p1=0.3, p2=0.15))
check("prop_detect_vals", d.get("calc") == "detect" and abs(d["h"] - 0.3617) < 1e-4, d)
# round-trips n -> power ~= 0.80
for _dd in (0.2, 0.5, 0.8):
    _n = strict("x", mt.med_power(calc="n", effect=_dd))["n1_before_dropout"]
    _p = strict("x", mt.med_power(calc="power", effect=_dd, n=_n))["p"]
    check(f"roundtrip_{_dd}", abs(_p - 0.8) < 0.01, (_n, _p))

# 5/6. OR undefined (not 0, no crash) at 0%/100% event rates
check("or_00", strict("x", mt.med_evidence(cer=0.0, eer=0.0)).get("or") is None)
check("or_11", strict("x", mt.med_evidence(cer=1.0, eer=1.0)).get("or") is None)

# 8. real-n vs assumed-n
da = strict("assumed", mt.med_evidence(cer=0.20, eer=0.12))
dr = strict("real", mt.med_evidence(cer=0.20, eer=0.12, n_c=500, n_e=500))
check("assumed_flag", da.get("n_assumed_per_arm") == 100 and "n_c" not in da, da)
check("real_flag", dr.get("n_c") == 500 and "n_assumed_per_arm" not in dr, dr)
check("real_narrower",
      (dr["arr_ci"][1] - dr["arr_ci"][0]) < (da["arr_ci"][1] - da["arr_ci"][0]))

# 9. unknown trial status -> clean error (no HTTP 400 leak)
d = strict("bad_status", mt.med_trial(q="diabetes", max_results=2, status="nonsense"))
check("bad_status_clean", "e" in d and "HTTP" not in d["e"] and "400" not in d["e"], d)

# 10. pspp/jmv without binaries: install hint + strict-JSON sanitize
check("pspp_hint", "apt install pspp" in strict("x", pspp.pspp_run(syntax="DESCRIPTIVES VARIABLES=ALL."))["e"])
check("jmv_hint", "apt install r-base" in strict("x", jmv.jmv_run(analysis="descriptives", data="a,b\n1,2"))["e"])
check("pspp_san", strict("x", pspp._ok({"x": float("nan")})) == {"x": None})  # noqa: SLF001
check("jmv_san", strict("x", jmv._ok({"x": float("inf")})) == {"x": None})  # noqa: SLF001

# 11. McNemar paired sizing (Miettinen normal approx)
d = strict("mcnemar_n", mt.med_power(calc="n", design="mcnemar", pdisc=0.3, oratio=2.0))
check("mcnemar_n_vals", d and d.get("t") == "mcnemar" and d["n_pairs"] == 234 and d["n_discordant"] == 71, d)
check("mcnemar_auto", strict("x", mt.med_power(calc="n", pdisc=0.3, oratio=2.0))["n_pairs"] == 234)
check("mcnemar_dropout", strict("x", mt.med_power(calc="n", design="mcnemar", pdisc=0.3, oratio=2.0, dropout=0.2))["n_pairs"] == 293)
check("mcnemar_ci", "or_ci" in d and d["or_ci"][0] < 2.0 < d["or_ci"][1], d)
_p = strict("x", mt.med_power(calc="power", design="mcnemar", pdisc=0.3, oratio=2.0, n=234))["p"]
check("mcnemar_roundtrip", abs(_p - 0.8) < 0.02, _p)
_od = strict("x", mt.med_power(calc="detect", design="mcnemar", pdisc=0.3, n=234))["oratio"]
check("mcnemar_detect", abs(_od - 2.0) < 0.05, _od)
check("mcnemar_or1", "e" in strict("x", mt.med_power(calc="n", design="mcnemar", pdisc=0.3, oratio=1.0)))
check("mcnemar_noparam", "e" in strict("x", mt.med_power(calc="n", design="mcnemar", pdisc=0.3)))

# 12. log-rank survival sizing (Schoenfeld/Freedman)
d = strict("logrank_n", mt.med_power(calc="n", design="logrank", hr=0.7, pevent=0.5))
check("logrank_vals", d and d.get("t") == "logrank" and d["events"] == 247 and d["n_total"] == 494, d)
_f = strict("x", mt.med_power(calc="n", design="logrank", hr=0.7, pevent=0.5, method="freedman"))
check("logrank_freedman", _f["events"] == 253 and _f["events"] >= d["events"], (_f["events"], d["events"]))
check("logrank_nopevent", "note" in strict("x", mt.med_power(calc="n", design="logrank", hr=0.7)))
_p = strict("x", mt.med_power(calc="power", design="logrank", hr=0.7, events=247))["p"]
check("logrank_roundtrip", abs(_p - 0.8) < 0.02, _p)
check("logrank_npevent", abs(strict("x", mt.med_power(calc="power", design="logrank", hr=0.7, n=494, pevent=0.5))["p"] - 0.8) < 0.02)
check("logrank_hr1", "e" in strict("x", mt.med_power(calc="n", design="logrank", hr=1.0)))
check("logrank_badcalc", "e" in strict("x", mt.med_power(calc="detect", design="logrank", hr=0.7)))

# 13. diagnostic-accuracy sizing (Wald normal approx)
d = strict("acc_n", mt.med_power(calc="n", design="acc", sens=0.9, width=0.05, prev=0.2))
check("acc_vals", d and d.get("t") == "accuracy" and d["n_diseased"] == 139 and d["n_total"] == 695, d)
_b = strict("x", mt.med_power(calc="n", design="acc", sens=0.9, spec=0.85, width=0.05, prev=0.2))
check("acc_both", _b["n_total"] == 695 and _b["n_total"] == max(_b["n_total_sens"], _b["n_total_spec"]) and _b["driven_by"] == "sens", _b)
check("acc_dropout", strict("x", mt.med_power(calc="n", design="acc", sens=0.9, width=0.05, prev=0.2, dropout=0.1))["n_total"] == 773)
check("acc_calc", "e" in strict("x", mt.med_power(calc="power", design="acc", sens=0.9, width=0.05, prev=0.2)))
check("acc_noparam", "e" in strict("x", mt.med_power(calc="n", design="acc", sens=0.9)))

# 14. registry handler forwards new params (no dropped params)
_h = _reg.registry.get_entry("med_power").handler
_d = json.loads(_h({"calc": "n", "design": "mcnemar", "pdisc": 0.3, "oratio": 2.0}))
check("handler_mcnemar", _d.get("t") == "mcnemar" and _d["n_pairs"] == 234, _d)
_d = json.loads(_h({"calc": "n", "design": "logrank", "hr": 0.7, "pevent": 0.5}))
check("handler_logrank", _d.get("t") == "logrank" and _d["events"] == 247, _d)
_d = json.loads(_h({"calc": "n", "design": "acc", "sens": 0.9, "width": 0.05, "prev": 0.2}))
check("handler_acc", _d.get("t") == "accuracy" and _d["n_total"] == 695, _d)

# 15. non-inferiority sizing (normal approx, one-sided alpha default 0.025)
d = strict("noninf_c_n", mt.med_power(calc="n", design="noninf", margin=0.5, sd=1.0))
check("noninf_c_vals", d and d.get("t") == "noninf" and d.get("endpoint") == "continuous"
      and d["n1"] == 63 and d["alpha"] == 0.025 and d.get("onesided") is True
      and "diff_ci" in d, d)
check("noninf_auto", strict("x", mt.med_power(calc="n", margin=0.5, sd=1.0))["n1"] == 63)
check("noninf_alias", strict("x", mt.med_power(calc="n", design="ni", margin=0.5, sd=1.0))["n1"] == 63)
check("noninf_dropout", strict("x", mt.med_power(calc="n", design="noninf", margin=0.5, sd=1.0, dropout=0.1))["n1"] == 70)
d = strict("noninf_b_n", mt.med_power(calc="n", design="noninf", margin=0.1, p_ctrl=0.3))
check("noninf_b_vals", d and d.get("endpoint") == "binary" and d["n1"] == 330
      and d["p_exp"] == 0.3 and "diff_ci" in d, d)
_p = strict("x", mt.med_power(calc="power", design="noninf", margin=0.5, sd=1.0, n=63))["p"]
check("noninf_c_roundtrip", abs(_p - 0.8) < 0.02, _p)
_p = strict("x", mt.med_power(calc="power", design="noninf", margin=0.1, p_ctrl=0.3, n=330))["p"]
check("noninf_b_roundtrip", abs(_p - 0.8) < 0.02, _p)
_m = strict("x", mt.med_power(calc="detect", design="noninf", sd=1.0, n=63))["margin"]
check("noninf_c_detect", abs(_m - 0.5) < 0.02, _m)
_m = strict("x", mt.med_power(calc="detect", design="noninf", p_ctrl=0.3, n=330))["margin"]
check("noninf_b_detect", abs(_m - 0.1) < 0.01, _m)
check("noninf_truediff", strict("x", mt.med_power(calc="n", design="noninf", margin=0.5, sd=1.0, true_diff=0.1))["n1"] == 44)
check("noninf_explicit_alpha", strict("x", mt.med_power(calc="n", design="noninf", margin=0.5, sd=1.0, alpha=0.05))["n1"] == 50)
check("noninf_margin0", "e" in strict("x", mt.med_power(calc="n", design="noninf", margin=0.0, sd=1.0)))
check("noninf_nomargin", "e" in strict("x", mt.med_power(calc="n", design="noninf", sd=1.0)))
check("noninf_noendpoint", "e" in strict("x", mt.med_power(calc="n", design="noninf", margin=0.5)))
check("noninf_both", "e" in strict("x", mt.med_power(calc="n", design="noninf", margin=0.5, sd=1.0, p_ctrl=0.3)))
check("noninf_badrate", "e" in strict("x", mt.med_power(calc="n", design="noninf", margin=0.1, p_ctrl=1.5)))
check("noninf_beyond", "e" in strict("x", mt.med_power(calc="n", design="noninf", margin=0.5, sd=1.0, true_diff=-0.5)))
check("noninf_badrate2", "e" in strict("x", mt.med_power(calc="n", design="noninf", margin=0.1, p_ctrl=0.3, true_diff=0.8)))

# 16. registry handler forwards noninf params (no dropped params)
_d = json.loads(_h({"calc": "n", "design": "noninf", "margin": 0.5, "sd": 1.0, "true_diff": 0.0}))
check("handler_noninf_c", _d.get("t") == "noninf" and _d["n1"] == 63 and _d["true_diff"] == 0.0, _d)
_d = json.loads(_h({"calc": "n", "design": "noninf", "margin": 0.1, "p_ctrl": 0.3}))
check("handler_noninf_b", _d.get("endpoint") == "binary" and _d["n1"] == 330 and _d["p_ctrl"] == 0.3, _d)
_d = json.loads(_h({"calc": "n", "margin": 0.5, "sd": 1.0}))
check("handler_noninf_auto", _d.get("t") == "noninf" and _d["n1"] == 63, _d)

def test_regression():
    assert FAIL == 0, f"{FAIL} regression checks failed"


if __name__ == "__main__":
    print(f"\n{PASS} passed, {FAIL} failed")
    sys.exit(1 if FAIL else 0)
