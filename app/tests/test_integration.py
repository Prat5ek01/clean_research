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
import textwrap
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


def _run_fail(tmp, script, *args, timeout=900):
    env = dict(os.environ)
    env["CLEANX_APP_DIR"] = str(tmp)
    env["PYTHONPATH"] = os.pathsep.join([str(APP / "src"), str(APP), env.get("PYTHONPATH", "")])
    r = subprocess.run([sys.executable, str(APP / script), *args], cwd=str(tmp), env=env, capture_output=True,
                       text=True, timeout=timeout)
    assert r.returncode != 0, f"{script} {' '.join(args)} should have failed:\n{r.stdout[-1500:]}"
    print(f"    ran {script} {' '.join(args)} (failed as expected)", flush=True)
    return r.stdout + r.stderr


def _blast_path(tmp):
    """cx_11 --verify and cx_04 --tool blast on BLAST-format tables with known hits, then a tamper check."""
    import hashlib
    from cleanx import blast as BL
    from cleanx import labels as L
    cx = tmp / "data" / "cleanx"
    tr = pd.read_csv(tmp / "data" / "split100.csv", sep="\t", dtype=str, keep_default_na=False)
    t_ids, t_seq = list(tr.Entry), dict(zip(tr.Entry, tr.Sequence))
    t_ec = dict(zip(tr.Entry, tr["EC number"]))
    paths = {"new392": tmp / "data" / "new.csv", "price149": tmp / "data" / "datasets" / "price.csv",
             "sp_temporal_cal": cx / "sets" / "sp_temporal_cal.csv",
             "sp_temporal_test": cx / "sets" / "sp_temporal_test.csv"}
    sets = {n: pd.read_csv(p, sep="\t", dtype=str, keep_default_na=False) for n, p in paths.items()}
    bdir = cx / "identity" / "blast"
    man = {"blast_version": "synthetic", "created_utc": "test", "sets": {},
           "database": {"name": "split100", "digest": BL.seq_digest(t_ids, t_seq)}}

    def hsp(q, qs, h, hs, pid, bits):
        n = min(len(qs), len(hs))
        nid = max(1, int(round(pid * n)))
        return "\t".join([q, h, f"{100 * nid / n:.3f}", str(n), str(len(qs)), str(len(hs)), f"{bits:.1f}",
                          "1e-30", str(nid)])

    def store(lines, fname):
        BL.write_hits(BL.reduce_hits(BL.parse_lines(lines)), bdir / fname)
        return hashlib.sha256((bdir / fname).read_bytes()).hexdigest()

    expect, pid_exp = {}, {}
    for k, (name, d) in enumerate(sets.items()):
        qs = dict(zip(d.Entry, d.Sequence))
        lines = []
        for j, q in enumerate(d.Entry):
            hits = [t_ids[(7 * j + 3 * k + r) % len(t_ids)] for r in range(3)]
            top = 0.25 + 0.7 * ((37 * j + k) % 100) / 100  # spread over all identity bins
            lines += [hsp(q, qs[q], h, t_seq[h], top - 0.1 * r, 300.0 - 50 * r) for r, h in enumerate(hits)]
            expect[(name, q)] = hits[0]
            n = min(len(qs[q]), len(t_seq[hits[0]]))
            pid_exp[(name, q)] = max(1, int(round(top * n))) / n
        man["sets"][name] = {"digest": BL.seq_digest(list(d.Entry), qs), "file": f"{name}.tsv.gz",
                             "file_sha256": store(lines, f"{name}.tsv.gz")}
    te, ca = sets["sp_temporal_test"], sets["sp_temporal_cal"]
    tq, cq = dict(zip(te.Entry, te.Sequence)), dict(zip(ca.Entry, ca.Sequence))
    c_ids = list(ca.Entry)
    lines = [hsp(q, tq[q], c_ids[j % len(c_ids)], cq[c_ids[j % len(c_ids)]], 0.85 if j < 12 else 0.3, 200.0)
             for j, q in enumerate(te.Entry)]
    man["pair"] = {"query": "sp_temporal_test", "target": "sp_temporal_cal", "file": "pair.tsv.gz",
                   "query_digest": BL.seq_digest(list(te.Entry), tq), "target_digest": BL.seq_digest(c_ids, cq),
                   "file_sha256": store(lines, "pair.tsv.gz")}
    (bdir / "manifest.json").write_text(json.dumps(man), encoding="utf-8")

    out = _run(tmp, "cx_11_blast.py", "--verify")
    assert "FAIL" not in out and "CX11 VERIFY DONE" in out, out[-800:]
    _run(tmp, "cx_04_identity.py", "--tool", "blast", "--sets", *paths)
    for name in paths:
        idt = pd.read_csv(cx / "identity" / f"{name}.csv", dtype={"Entry": str, "top_hit": str})
        assert (idt.tool == "blast").all()
        assert all(expect[(name, q)] == h for q, h in zip(idt.Entry, idt.top_hit)), name
        assert all(abs(pid_exp[(name, q)] - v) < 1e-9 for q, v in zip(idt.Entry, idt.max_pid_short)), name
        assert idt.identity_bin.nunique() >= 4, name
    hp = L.read_pred_csv(tmp / "results" / "cleanx" / "preds" / "new392_homology_tophit.csv")
    for q, items in hp.items():
        assert [e for e, _ in items] == L.split_ecs(t_ec[expect[("new392", q)]]), q
    _run(tmp, "cx_04_identity.py", "--tool", "blast", "--pair", "sp_temporal_test:sp_temporal_cal")
    pr = pd.read_csv(cx / "identity" / "sp_temporal_test_vs_sp_temporal_cal.csv", dtype={"Entry": str})
    assert (pr.tool == "blast").all() and int((pr.max_pid_short >= 0.5).sum()) == 12

    # a changed sequence must be refused by both the verifier and the parser
    p = paths["sp_temporal_test"]
    orig = p.read_bytes()
    d = te.copy()
    d.loc[0, "Sequence"] = d.loc[0, "Sequence"][:-1] + ("A" if d.loc[0, "Sequence"][-1] != "A" else "C")
    d.to_csv(p, sep="\t", index=False)
    assert "FAIL" in _run_fail(tmp, "cx_11_blast.py", "--verify")
    assert "mismatch" in _run_fail(tmp, "cx_04_identity.py", "--tool", "blast", "--sets", "sp_temporal_test")
    p.write_bytes(orig)


def test_drift_end_to_end_synthetic_only():
    from cleanx import labels as L
    from cx_12_drift import synthetic, METHODS
    before = _snapshot()
    with tempfile.TemporaryDirectory() as td:
        tmp = Path(td) / "app"
        ids, truth, raw, _ = synthetic()
        sets = tmp / "data" / "cleanx" / "sets"
        seqs = {q: "ACDEFGHIKLMNPQRSTVWY" for q in ids}
        L.write_dataset(sets / "sp_temporal_test.csv", ids, truth, seqs)
        for bb in ("split100", "split70"):
            L.write_dataset(tmp / "data" / f"{bb}.csv", ids, {q: ["1.1.1.1"] for q in ids}, seqs)
            for m in METHODS:
                L.write_pred_csv(tmp / "results" / "cleanx" / "preds" / f"sp_temporal_test_{bb}_{m}_maxsep.csv",
                                 {q: [(e, None) for e in raw[q]] for q in ids})
        L.write_pred_csv(tmp / "results" / "cleanx" / "preds" / "sp_temporal_test_homology_tophit.csv",
                         {q: [(e, None) for e in raw[q]] for q in ids})
        pd.DataFrame({"accession": ids, "protein_existence": ["Evidence at protein level"] * len(ids)}).to_csv(
            sets / "sp_temporal_meta.csv", index=False)
        nom = tmp / "data" / "cleanx" / "downloads" / "enzyme.dat"
        nom.parent.mkdir(parents=True)
        nom.write_text("ID   1.1.1.1\nDE   Transferred entry: 1.1.1.2 and 1.1.1.3.\n//\n", encoding="utf-8")
        ident = tmp / "data" / "cleanx" / "identity" / "sp_temporal_test.csv"
        ident.parent.mkdir(parents=True)
        pd.DataFrame({"Entry": ids, "max_pid_short": [0.2, 0.6] * 30}).to_csv(ident, index=False)
        out = _run(tmp, "cx_12_drift.py", "--posthoc", "--bootstrap", "20")
        assert "CX12 POSTHOC PASS" in out
        res = tmp / "results" / "cleanx" / "E7_posthoc"
        d = pd.read_csv(res / "E7_deltas.csv")
        assert set(d.backbone) == {"split100", "split70"}
        assert {"closed_PE1", "closed_capped5", "identity:<30%"} <= set(d.subset)
        assert d[(d.subset == "all") & (d.metric == "micro_F1")].delta_obs.min() > 0.4
        man = json.loads((res / "E7_manifest.json").read_text())
        assert man["status"] == "posthoc" and len(man["files"]) == 15
        for ext in ("png", "pdf", "csv"):
            assert (res / f"figE7_drift.{ext}").stat().st_size > 100
        out = _run(tmp, "cx_12_drift.py", "--smoke")
        assert "CX12 SMOKE PASS" in out
        pred = tmp / "results" / "cleanx" / "preds" / "sp_temporal_test_split100_centroid_maxsep.csv"
        pred.write_text("q00,EC:1.1.1.1\n", encoding="utf-8")
        assert "ID mismatch" in _run_fail(tmp, "cx_12_drift.py", "--posthoc", "--bootstrap", "20")
    assert _snapshot() == before, "E7 integration touched real results"


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

        _blast_path(tmp)

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
        assert {"status", "n_needed", "deficit", "years_of_accrual_for_deficit",
                "years_at_deficit_ci_lo"} <= set(b.columns)
        lo = b[b.deficit_ci_lo.notna() & (b.deficit_ci_lo > 0) & (b.deficit_ci_lo < float("inf"))]
        assert (lo.years_at_deficit_ci_lo <= lo.years_of_accrual_for_deficit + 1e-9).all()
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
        assert "does not reproduce" not in out and "WARNING" not in out, out[-800:]
        figs = res / "figures"
        for f in ("fig2_tiers", "fig3_identity", "fig10_primary_endpoint", "fig5_reliability_by_tier",
                  "fig6_risk_coverage", "fig7_tier_risk_control", "fig9_certification_budget"):
            for ext in ("png", "pdf", "csv"):
                assert (figs / f"{f}.{ext}").is_file(), f"{f}.{ext}"
        f2 = pd.read_csv(figs / "fig2_tiers.csv")
        assert f2.plotted.any() and ((0 <= f2.ci_lo) & (f2.ci_lo <= f2.ci_hi) & (f2.ci_hi <= 1)).all()
        f10 = pd.read_csv(figs / "fig10_primary_endpoint.csv")
        prim = f10[(f10.kind == "primary") & (f10.metric == "macro_union_F1")]
        ref = d[(d.dataset == "sp_temporal_test") & (d.subset == "closed") & (d.metric == "macro_union_F1") &
                (d.comparison == "balanced3 - centroid")]
        assert len(prim) == 1 and abs(prim.delta_obs.iloc[0] - ref.delta_obs.iloc[0]) < 1e-12
        f7 = pd.read_csv(figs / "fig7_tier_risk_control.csv")
        assert (f7.errors == (f7.selective_risk.fillna(0) * f7.accepted).round()).all()

        # fig7 with accepted proteins (the synthetic CAL split is too small to certify risk 0.2), from a risk-control
        # table in the format written before error counts were reported: exact Clopper-Pearson intervals expected
        (tmp / "fig7_check.py").write_text(textwrap.dedent("""
            import pandas as pd
            import cx_08_figures as F
            F.style()
            rows = []
            for g in ("centroid", "exemplar"):
                for rule in ("global", "mondrian"):
                    for grp, n, acc, err in (("all", 100, 40, 6), ("pred_tier:common", 60, 30, 3),
                                             ("pred_tier:rare", 40, 10 if rule == "global" else 0,
                                              4 if rule == "global" else 0)):
                        rows.append({"geometry": g, "alpha": 0.2, "delta": 0.1, "rule": rule,
                                     "dataset": "sp_temporal_test", "group": grp, "n": n, "accepted": acc,
                                     "coverage": acc / n, "selective_risk": err / acc if acc else float("nan"),
                                     "theta": 0.5})
            F.fig_tier_risk_control(pd.DataFrame(rows))
            print("FIG7 CHECK DONE")
        """), encoding="utf-8")
        assert "FIG7 CHECK DONE" in _run(tmp, str(tmp / "fig7_check.py"))
        from scipy.stats import beta
        f7 = pd.read_csv(figs / "fig7_tier_risk_control.csv")
        acc = f7[f7.accepted > 0]
        assert len(acc) == 10 and len(f7) == 12
        for k, n in ((6, 40), (3, 30), (4, 10)):
            r = acc[(acc.errors == k) & (acc.accepted == n)]
            assert len(r) and (abs(r.risk_ci_lo - beta.ppf(0.025, k, n - k + 1)) < 1e-12).all()
            assert (abs(r.risk_ci_hi - beta.ppf(0.975, k + 1, n - k)) < 1e-12).all()
    after = _snapshot()
    assert after == before, "the integration test modified the real results folder"
