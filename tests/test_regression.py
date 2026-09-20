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

def test_regression():
    assert FAIL == 0, f"{FAIL} regression checks failed"


if __name__ == "__main__":
    print(f"\n{PASS} passed, {FAIL} failed")
    sys.exit(1 if FAIL else 0)
