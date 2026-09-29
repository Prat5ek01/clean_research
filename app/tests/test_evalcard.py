"""Known-answer and property tests for cleanx.evalcard (Props. 1-4) and cleanx.ecbench parsing."""
import tempfile
from pathlib import Path

import numpy as np
from sklearn.metrics import f1_score, precision_score, recall_score
from sklearn.preprocessing import MultiLabelBinarizer

from cleanx import ecbench as EB
from cleanx.evalcard import augment_out_of_space, card, gaming_curve, invariance_check

X, Y, Z, W = "1.1.1.1", "2.2.2.2", "3.3.3.3", "4.4.4.4"


def _random(seed, n=50):
    rng = np.random.default_rng(seed)
    vocab = [f"{a}.{b}.1.{d}" for a in (1, 2, 3) for b in (1, 2) for d in range(1, 8)]
    ids = [f"q{i}" for i in range(n)]
    truth = {i: list(rng.choice(vocab[:20], size=int(rng.integers(1, 3)), replace=False)) for i in ids}
    preds = {i: list(rng.choice(vocab, size=int(rng.integers(0, 4)), replace=False)) for i in ids}
    return truth, preds, ids


def test_hand_computed_example():
    truth = {"a": [X], "b": [Y]}
    preds = {"a": [X, Z], "b": [W]}  # Z, W: wrong and outside the test label set {X, Y}
    c = card(truth, preds, ["a", "b"], 4)
    exp = {"reported_CLEAN_P": 0.5, "reported_CLEAN_R": 0.5, "reported_CLEAN_F1": 0.5,
           "reported_ECBench_F1": 0.5, "u_micro_P": 1 / 3, "u_micro_R": 0.5, "u_micro_F1": 0.4,
           "u_macro_union_F1": 0.25, "u_samples_F1": 1 / 3, "t_micro_P": 1.0, "t_micro_F1": 2 / 3,
           "fp": 2, "fp_invisible": 2, "invisible_share": 1.0, "n_true_classes": 2}
    for k, v in exp.items():
        assert abs(c[k] - v) < 1e-12, (k, c[k], v)


def test_prop1_invariance_and_prop2_sensitivity_random():
    for seed in range(30):
        truth, preds, ids = _random(seed)
        res = invariance_check(truth, preds, ids, 4, n_extra=3, seed=seed)
        assert res["ok"], (seed, res)
        assert res["max_change_invariant_scores"] < 1e-12


def test_lemma3_weighted_recall_equals_micro_recall():
    for seed in range(20):
        truth, preds, ids = _random(seed)
        for depth in (1, 2, 3, 4):
            c = card(truth, preds, ids, depth)
            assert abs(c["reported_ECBench_R"] - c["u_micro_R"]) < 1e-12
            assert abs(c["reported_CLEAN_R"] - c["t_micro_R"]) < 1e-12


def test_card_matches_both_published_protocols_in_sklearn():
    for seed in range(15):
        truth, preds, ids = _random(seed)
        c = card(truth, preds, ids, 4)
        true_space = sorted({e for i in ids for e in truth[i]})
        mlb = MultiLabelBinarizer(classes=true_space).fit([true_space])  # CLEAN: fitted on test labels only
        T = mlb.transform([truth[i] for i in ids])
        P = mlb.transform([[e for e in preds[i] if e in true_space] for i in ids])
        for k, f in (("P", precision_score), ("R", recall_score), ("F1", f1_score)):
            assert abs(c[f"reported_CLEAN_{k}"] - f(T, P, average="weighted", zero_division=0)) < 1e-12
        union = sorted({e for i in ids for e in truth[i] + preds[i]})  # EC-Bench: union label space
        mlb = MultiLabelBinarizer(classes=union).fit([union])
        Tu, Pu = mlb.transform([truth[i] for i in ids]), mlb.transform([preds[i] for i in ids])
        assert abs(c["reported_ECBench_F1"] - f1_score(Tu, Pu, average="weighted", zero_division=0)) < 1e-12
        assert abs(c["u_micro_F1"] - f1_score(Tu, Pu, average="micro", zero_division=0)) < 1e-12


def test_augmentation_only_adds_phantom_labels():
    truth, preds, ids = _random(0)
    aug, added = augment_out_of_space(preds, ids, n_extra=3, seed=0)
    for i in ids:
        assert aug[i][:len(preds[i])] == preds[i]
        assert len(aug[i]) - len(preds[i]) == added[i]
        assert all(e.startswith("99.99.99.") for e in aug[i][len(preds[i]):])


def test_gaming_curve_recall_non_decreasing():
    truth, preds, ids = _random(5)
    ranked = {i: list(dict.fromkeys(truth[i][:1] + preds[i] + [f"9.9.9.{k}" for k in range(10)])) for i in ids}
    rows = gaming_curve(truth, ranked, ids, ks=(1, 2, 3, 5, 10))
    rec = [r["u_micro_R"] for r in rows]
    assert all(b >= a - 1e-15 for a, b in zip(rec, rec[1:]))
    assert rows[-1]["n_pred_labels"] > rows[0]["n_pred_labels"]


def test_ecbench_parsing_and_synthetic_file():
    assert EB.parse_cell(" 1.1.1.1, - ,2.2.2.2,1.1.1.1,") == ["1.1.1.1", "2.2.2.2"]
    assert EB.parse_cell("-") == [] and EB.parse_cell("nan") == []
    with tempfile.TemporaryDirectory() as td:
        p = EB.write_synthetic(Path(td) / "syn.csv", n=30, seed=2)
        ids, truth, preds = EB.read_results(p)
        assert len(ids) == 30 and set(preds) == set(EB.METHODS)
        assert all(truth[i] for i in ids)


def test_ecbench_pinned_hashes_fail_loudly():
    assert set(EB.SHA256) == set(EB.FILES.values())
    with tempfile.TemporaryDirectory() as td:
        p = Path(td) / "results_test_100_go_3.csv"
        p.write_text("id,true\n")
        try:
            EB.verify({("test", "100"): p})
        except ValueError as e:
            assert "pinned" in str(e)
        else:
            raise AssertionError("a file with the wrong SHA-256 was accepted")
    cached = {k: EB.DEST / v for k, v in EB.FILES.items()}
    if all(q.is_file() for q in cached.values()):  # the real downloads, when present, match the pinned bytes
        assert EB.verify(cached) == cached
