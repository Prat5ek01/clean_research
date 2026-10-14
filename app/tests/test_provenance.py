"""Known-answer tests for E7: UniProt evidence parsing and label classes (cleanx.provenance), stratum contrasts,
difference-in-differences and identity standardisation (cleanx.contrast), and the pre-registration guard
(cleanx.freeze) on a temporary git repository."""
import shutil
import subprocess
import tempfile
from pathlib import Path

import numpy as np

from cleanx import contrast as C
from cleanx import freeze
from cleanx import metrics as M
from cleanx import provenance as PV

E, S, R, O = (("ECO:0000269", "PubMed", "1"), ("ECO:0000250", "UniProtKB", "P11111-2"),
              ("ECO:0000255", "HAMAP-Rule", "MF_00001"), ("ECO:0000305", "", ""))


def test_evidence_from_names_reactions_and_domains():
    e = PV.synthetic_entry("A1", {"1.1.1.1": [E]}, in_reaction={"1.1.1.1": [S], "2.2.2.2": [R]},
                           contains={"3.3.3.3": [O]})
    ev = PV.ec_evidence(e)
    assert set(ev) == {"1.1.1.1", "2.2.2.2", "3.3.3.3"}
    assert ev["1.1.1.1"]["codes"] == {"ECO:0000269", "ECO:0000250"} and ev["1.1.1.1"]["in_name"] and \
        ev["1.1.1.1"]["in_reaction"]
    assert ev["1.1.1.1"]["sim_sources"] == {"P11111-2"} and PV.source_accession("P11111-2") == "P11111"
    assert ev["2.2.2.2"]["model_sources"] == {"HAMAP-Rule:MF_00001"} and not ev["2.2.2.2"]["in_name"]
    assert ev["3.3.3.3"]["codes"] == {"ECO:0000305"} and ev["3.3.3.3"]["in_name"]


def test_class_precedence():
    assert PV.ec_class(["ECO:0000269", "ECO:0000250", "ECO:0000255"]) == "EXP"
    assert PV.ec_class(["ECO:0000250", "ECO:0000255"]) == "SIM"
    assert PV.ec_class(["ECO:0000256"]) == "RULE" and PV.ec_class(["ECO:0000259"]) == "RULE"
    assert PV.ec_class(["ECO:0000305", "ECO:0000312"]) == "OTHER" and PV.ec_class([]) == "OTHER"
    codes = {"a": ["ECO:0000269"], "b": ["ECO:0000250"], "c": ["ECO:0000255"], "d": ["ECO:0000305"]}
    assert PV.protein_class(["a"], codes) == "EXP"
    assert PV.protein_class(["a", "b"], codes) == "SIM"        # any SIM label
    assert PV.protein_class(["a", "c"], codes) == "RULE"       # not all EXP, a RULE label
    assert PV.protein_class(["a", "d"], codes) == "OTHER"
    assert PV.protein_class(["c", "d"], codes) == "RULE"
    assert PV.protein_class(["a", "z"], codes) == "LABEL_ABSENT"  # a label not on the entry
    assert PV.protein_class(["a"], codes, found=False) == "LABEL_ABSENT"
    assert PV.protein_class([], codes) == "LABEL_ABSENT"
    assert PV.PROPAGATED == {"SIM", "RULE"}


def test_accession_resolution():
    entries = [PV.synthetic_entry("P1", {}, secondary=["S1", "S2"]), PV.synthetic_entry("P2", {}, secondary=["S2"]),
               PV.synthetic_entry("M1", {}, inactive="MERGED", merged_into=["P2"]),
               PV.synthetic_entry("M2", {}, inactive="MERGED", merged_into=["X9"]),
               PV.synthetic_entry("D1", {}, inactive="DELETED"),
               PV.synthetic_entry("P1", {})]  # the same entry read twice (batch and single response)
    idx = PV.index_entries(entries)
    got = {q: (PV.resolve(q, idx)[0] or {}).get("primaryAccession", "-") + "/" + PV.resolve(q, idx)[1]
           for q in ("P1", "S1", "S2", "M1", "M2", "D1", "Q0")}
    assert got == {"P1": "P1/primary", "S1": "-/absent", "S2": "P2/secondary", "M1": "P2/merged",
                   "M2": "-/inactive:MERGED", "D1": "-/inactive:DELETED", "Q0": "-/absent"}, got
    # The later P1 replaces its aliases: S1 is unknown; only P2 claims S2.
    idx2 = PV.index_entries(entries[:2])
    assert PV.resolve("S2", idx2) == (None, "ambiguous_secondary") and PV.resolve("S1", idx2)[1] == "secondary"


def test_tables_rows_and_read_json_round_trip():
    import gzip
    import json
    entries = [PV.synthetic_entry("P1", {"1.1.1.1": [E, S]}), PV.synthetic_entry("D1", {}, inactive="DELETED")]
    with tempfile.TemporaryDirectory() as td:
        p = Path(td) / "r.json.gz"
        p.write_bytes(gzip.compress(json.dumps({"results": entries}).encode("utf-8"), mtime=0))
        back = PV.read_json_results(p)
        q = Path(td) / "single.json"
        q.write_text(json.dumps(entries[0]), encoding="utf-8")
        assert PV.read_json_results(q) == [entries[0]]
    er, cr = PV.tables(["P1", "D1"], back, via={"P1": "batch"})
    assert [(r["query_accession"], r["found"], r["resolution"], r["via"]) for r in er] == \
        [("P1", True, "primary", "batch"), ("D1", False, "inactive:DELETED", "")]
    assert len(cr) == 1 and cr[0]["codes"] == "ECO:0000250;ECO:0000269" and cr[0]["sim_sources"] == "P11111-2"


def _world(n, acc_a, acc_b, seed=0):
    """n single-label proteins; method a correct on the first acc_a, b on the first acc_b (micro-F1 = accuracy)."""
    ids = [f"q{seed}_{k}" for k in range(n)]
    truth = {i: [f"1.1.1.{k}"] for k, i in enumerate(ids)}
    pa = {i: truth[i] if k < acc_a else ["9.9.9.9"] for k, i in enumerate(ids)}
    pb = {i: truth[i] if k < acc_b else ["8.8.8.8"] for k, i in enumerate(ids)}
    return ids, truth, pa, pb


def test_delta_and_did_known_answer():
    i0, t0, a0, b0 = _world(50, 40, 30, 0)   # delta0 = 0.2
    i1, t1, a1, b1 = _world(40, 20, 30, 1)   # delta1 = -0.25
    truth, pa, pb = {**t0, **t1}, {**a0, **a1}, {**b0, **b1}
    o, boot = C.delta(truth, pa, pb, i0, "micro_F1", 300, 5)
    assert abs(o - 0.2) < 1e-12 and boot.shape == (300,)
    r = C.did(truth, pa, pb, i0, i1, "micro_F1", 300, 5)
    assert abs(r["did"]["estimate"] + 0.45) < 1e-12 and abs(r["delta1"]["estimate"] + 0.25) < 1e-12
    assert np.allclose(r["boot1"] - r["boot0"], r["boot1"] - r["boot0"])
    b0_again = M.bootstrap(truth, {"a": pa, "b": pb}, i0, B=300, seed=5, keys=("micro_F1",))
    assert np.allclose(r["boot0"], b0_again["a"]["micro_F1"] - b0_again["b"]["micro_F1"])  # seed for S0
    b1_again = M.bootstrap(truth, {"a": pa, "b": pb}, i1, B=300, seed=6, keys=("micro_F1",))
    assert np.allclose(r["boot1"], b1_again["a"]["micro_F1"] - b1_again["b"]["micro_F1"])  # seed + 1 for S1
    assert r["did"]["ci_lo"] < -0.45 < r["did"]["ci_hi"] and r["did"]["p_two_sided"] < 0.01
    assert abs(C.mde80(0.1) - 0.2801585218) < 1e-9


def test_two_sided_p():
    assert C.two_sided_p(np.array([1.0, 2.0, 3.0])) == 0.0
    assert C.two_sided_p(np.array([-1.0, 1.0])) == 1.0
    assert abs(C.two_sided_p(np.array([-1.0, 1.0, 2.0, 3.0])) - 0.5) < 1e-12
    assert C.two_sided_p(np.zeros(5)) == 1.0


def test_merge_bins_rule():
    order = ["<30%", "30-50%", "50-70%", "70-90%", ">=90%"]
    full = {"a": {b: 20 for b in order}, "b": {b: 20 for b in order}}
    assert C.merge_bins(order, full, 10) == [[b] for b in order]
    top_short = {"a": {**{b: 20 for b in order}, ">=90%": 3}, "b": {b: 20 for b in order}}
    assert C.merge_bins(order, top_short, 10) == [["<30%"], ["30-50%"], ["50-70%"], ["70-90%", ">=90%"]]
    low_short = {"a": {**{b: 20 for b in order}, "<30%": 2}, "b": {b: 20 for b in order}}
    assert C.merge_bins(order, low_short, 10) == [["<30%", "30-50%"], ["50-70%"], ["70-90%"], [">=90%"]]
    mid_short = {"a": {"<30%": 20, "30-50%": 4, "50-70%": 4, "70-90%": 20, ">=90%": 20}, "b": full["b"]}
    # 50-70% joins 30-50% (8 < 10, still short), and that pair joins <30%
    assert C.merge_bins(order, mid_short, 10) == [["<30%", "30-50%", "50-70%"], ["70-90%"], [">=90%"]]
    mid_ok = {"a": {"<30%": 20, "30-50%": 6, "50-70%": 6, "70-90%": 20, ">=90%": 20}, "b": full["b"]}
    assert C.merge_bins(order, mid_ok, 10) == [["<30%"], ["30-50%", "50-70%"], ["70-90%"], [">=90%"]]
    assert C.merge_bins(order, {"a": {"<30%": 9}}, 10) is None


def test_standardized_known_answer():
    # stratum X: bin lo delta 0 (20 proteins), bin hi delta 0.5 (20); stratum Y: lo delta 0.1 (10), hi delta 0.25 (40)
    parts = [("x_lo", 20, 10, 10), ("x_hi", 20, 15, 5), ("y_lo", 10, 6, 5), ("y_hi", 40, 30, 20)]
    truth, pa, pb, bins, strata = {}, {}, {}, {}, {"X": [], "Y": []}
    for k, (name, n, a, b) in enumerate(parts):
        ids, t, x, y = _world(n, a, b, k)
        truth.update(t), pa.update(x), pb.update(y)
        for i in ids:
            bins[i] = "<30%" if name.endswith("lo") else ">=90%"
        strata["X" if name.startswith("x") else "Y"] += ids
    res = C.standardized(truth, pa, pb, strata, bins, ["<30%", ">=90%"], "micro_F1", 200, 3, min_n=10)
    w_lo, w_hi = 30 / 90, 60 / 90
    assert res["groups"] == ["<30%", ">=90%"] and np.allclose(res["weights"], [w_lo, w_hi])
    assert abs(res["strata"]["X"][0] - (w_lo * 0.0 + w_hi * 0.5)) < 1e-12
    assert abs(res["strata"]["Y"][0] - (w_lo * 0.1 + w_hi * 0.25)) < 1e-12
    assert res["strata"]["X"][1].shape == (200,)
    assert C.standardized(truth, pa, pb, strata, bins, ["<30%", ">=90%"], "micro_F1", 50, 3, min_n=100) is None


def _git(cwd, *args):
    subprocess.run(["git", *args], cwd=str(cwd), check=True, capture_output=True, text=True)


def test_freeze_guard_on_a_temporary_repository():
    if shutil.which("git") is None:
        print("    git not found: freeze guard test not run")
        return
    with tempfile.TemporaryDirectory() as td:
        root = Path(td)
        _git(root, "init", "-q")
        _git(root, "config", "user.email", "t@example.org")
        _git(root, "config", "user.name", "t")
        f = root / "analysis.py"
        f.write_text("x = 1\n", encoding="utf-8")
        st = freeze.code_state([f], root)
        assert st["git"] and st["commit"] is None or st["dirty"], st  # no commit yet: not frozen
        try:
            freeze.require_frozen([f], False, "t", root)
            raise AssertionError("uncommitted code must stop a confirmatory run")
        except SystemExit:
            pass
        _git(root, "add", "analysis.py")
        _git(root, "commit", "-q", "-m", "freeze")
        label = freeze.require_frozen([f], False, "t", root)
        assert label.startswith("confirmatory @"), label
        f.write_text("x = 2\n", encoding="utf-8")
        assert freeze.require_frozen([f], True, "t", root).startswith("NOT CONFIRMATORY (uncommitted changes")
        d = freeze.digests([f], root)
        assert list(d) == ["analysis.py"] and len(d["analysis.py"]) == 64
