"""End-to-end integration test on synthetic data: builds a synthetic app in a temporary folder (tests/synth_app.py)
and runs the real scripts against it through CLEANX_APP_DIR. Covers every new or modified path: cx_04 --pair,
cx_05 (observed deltas, Holm, sensitivity subsets), cx_06 (d_own / d_other), cx_07 (pure / mixed open set, EC3
fallback), cx_09 --skip-ecbench and --card, cx_10 (certification budget) and cx_08 (figures, including the
fallback to pilot results when the confirmatory folder is incomplete). Never touches the real data/
or results/ folders (checked)."""
import json
import os
import subprocess
import sys
import tempfile
import time
from pathlib import Path

import pandas as pd

APP = Path(__file__).resolve().parents[1]
REAL_RESULTS = APP / "results" / "cleanx"


def _snapshot():
    return {str(p.relative_to(REAL_RESULTS)): p.stat().st_mtime for p in REAL_RESULTS.rglob("*") if p.is_file()} \
        if REAL_RESULTS.is_dir() else {}


def _run(tmp, script, *args, timeout=900):
    env = dict(os.environ)
    env["CLEANX_APP_DIR"] = str(tmp)
    env["PYTHONPATH"] = os.pathsep.join([str(APP / "src"), str(APP), env.get("PYTHONPATH", "")])
    env["PYTHONUNBUFFERED"] = "1"
    t0 = time.time()
    r = subprocess.run([sys.executable, str(APP / script), *args], cwd=str(tmp), env=env, capture_output=True,
                       text=True, timeout=timeout)
    if r.returncode != 0:
        raise AssertionError(f"{script} {' '.join(args)} failed (exit {r.returncode}):\n{r.stdout[-3000:]}\n"
                             f"{r.stderr[-3000:]}")
    print(f"    ran {script} {' '.join(args)} ({time.time() - t0:.0f}s)", flush=True)
    return r.stdout


def test_pipeline_end_to_end_on_synthetic_app():
    import synth_app
    before = _snapshot()
    with tempfile.TemporaryDirectory(ignore_cleanup_errors=True) as td:
        tmp = Path(td) / "app"
        synth_app.build(tmp, seed=0)
        res = tmp / "results" / "cleanx"

        out = _run(tmp, "cx_04_identity.py", "--pair", "sp_temporal_test:sp_temporal_cal", "--topn", "10")
        ident = pd.read_csv(tmp / "data" / "cleanx" / "identity" / "sp_temporal_test_vs_sp_temporal_cal.csv")
        assert len(ident) == 110 and (ident.max_pid_short >= 0.5).sum() >= 10, out[-500:]

        _run(tmp, "cx_05_geometry.py", "--bootstrap", "200")
        d = pd.read_csv(res / "E2_split100" / "E2_deltas.csv")
        assert {"delta_obs", "p_holm", "n"} <= set(d.columns)
        subs = set(d.subset)
        for s in ("closed", "closed|cluster_bootstrap_EC", "closed_cal_disjoint", "closed_PE1"):
            assert s in subs, (s, sorted(subs))
        t = d[d.subset.str.startswith("tier:")]
        assert len(t) and t.p_holm.notna().all() and (t.p_holm >= t.p_two_sided - 1e-12).all()
        closed = d[(d.subset == "closed") & (d.dataset == "sp_temporal_test")]
        m = pd.read_csv(res / "E2_split100" / "E2_metrics_long.csv")
        pt = m[(m.dataset == "sp_temporal_test") & (m.subset_type == "closed") & (m.depth == 4)].set_index("method")
        r = closed[(closed.comparison == "exemplar - centroid") & (closed.metric == "micro_F1")].iloc[0]
        assert abs(r.delta_obs - (pt.loc["exemplar", "micro_F1"] - pt.loc["centroid", "micro_F1"])) < 1e-12
        tuning = json.loads((res / "E2_split100" / "E2_tuning.json").read_text())
        assert tuning["tuning_set"] == "sp_temporal_cal"

        _run(tmp, "cx_06_ec_biology.py", "--smoke")
        st = json.loads((res / "E3_split100_smoke" / "E3_stats.json").read_text())
        assert {"mean_d_own", "mean_d_other"} <= set(st["size_stratified_permutation"]), st.keys()
        assert st["n_nise_like"] >= 2

        _run(tmp, "cx_07_reliability.py")
        os_ = pd.read_csv(res / "E4" / "E4_openset.csv")
        groups = set(os_.group)
        assert {"unseen_EC_enzymes", "unseen_EC_enzymes_pure", "unseen_EC_enzymes_mixed", "non_enzyme"} <= groups
        feats = pd.read_csv(res / "E4" / "E4_features_sp_temporal_test_centroid.csv", keep_default_na=False)
        assert set(feats.unseen_kind) >= {"pure", "mixed", ""}
        assert (res / "E4" / "E4_fallback.csv").is_file()

        _run(tmp, "cx_09_evalcard.py", "--skip-ecbench")
        inv = pd.read_csv(res / "E5_local" / "E5_invariance.csv")
        assert len(inv) and inv.ok.all()
        cards = pd.read_csv(res / "E5_local" / "E5_cards.csv", dtype={"regime": str})
        assert (cards.source == "temporal").any() and (cards.source == "legacy").any()
        assert (res / "E5_local" / "E5_benchmark_card.csv").is_file()

        pred = res / "preds" / "new392_split100_centroid_maxsep.csv"
        _run(tmp, "cx_09_evalcard.py", "--card", str(tmp / "data" / "new.csv"), str(pred))
        c = pd.read_csv(res / "E5_card" / f"E5_card_{pred.stem}.csv")
        assert list(c.depth) == [1, 2, 3, 4] and (c.fp_invisible <= c.fp).all()
        ref = m[(m.dataset == "new392") & (m.method == "centroid") & (m.subset == "all") & (m.depth == 4)]
        assert abs(c.u_micro_F1.iloc[-1] - ref.micro_F1.iloc[0]) < 1e-12
        assert not (res / "E5").exists(), "--card must not write into the confirmatory E5 folder"

        _run(tmp, "cx_10_certbudget.py", "--bootstrap", "20")
        b = pd.read_csv(res / "E6" / "E6_tier_budget.csv")
        assert {"status", "n_needed", "deficit", "years_of_accrual_for_deficit"} <= set(b.columns)
        assert set(b.status.dropna()) <= {"estimated", "risk_above_alpha", "too_few_calibration_proteins"}

        # a results folder without the needed files (e.g. a partial run) must not block the pilot fallback
        (res / "E5").mkdir()
        (res / "E5_pilot").mkdir()
        (res / "E5_pilot" / "E5_cards.csv").write_bytes((res / "E5_local" / "E5_cards.csv").read_bytes())
        (tmp / "pick_check.py").write_text("import cx_08_figures as F\n"
                                           "print('PICK', F._pick(('E5_cards.csv',), 'E5', 'E5_pilot'))\n")
        assert "PICK E5_pilot" in _run(tmp, str(tmp / "pick_check.py"))
        (res / "E5_pilot" / "E5_cards.csv").unlink()
        (res / "E5_pilot").rmdir()

        out = _run(tmp, "cx_08_figures.py")
        assert "skip fig8" in out, out[-800:]
        figs = res / "figures"
        for f in ("fig2_tiers", "fig5_reliability_by_tier", "fig6_risk_coverage", "fig9_certification_budget"):
            for ext in ("png", "pdf", "csv"):
                assert (figs / f"{f}.{ext}").is_file(), f"{f}.{ext}"
    after = _snapshot()
    assert after == before, "the integration test modified the real results folder"
