"""Known-answer tests for the additions to cleanx.metrics (Holm, BH, cluster bootstrap, observed delta)."""
import numpy as np

from cleanx.metrics import bh, bootstrap, evaluate, holm, observed_delta


def _toy(n=40, seed=1):
    rng = np.random.default_rng(seed)
    vocab = [f"1.1.1.{k}" for k in range(12)]
    ids = [f"q{i}" for i in range(n)]
    truth = {i: [vocab[int(rng.integers(0, 8))]] for i in ids}
    pa = {i: ([truth[i][0]] if rng.random() < 0.6 else [vocab[int(rng.integers(0, 12))]]) for i in ids}
    pb = {i: ([truth[i][0]] if rng.random() < 0.4 else [vocab[int(rng.integers(0, 12))]]) for i in ids}
    return truth, pa, pb, ids


def test_holm_known_answer():
    assert np.allclose(holm([0.01, 0.04, 0.03, 0.005]), [0.03, 0.06, 0.06, 0.02])
    assert np.allclose(holm([0.001, 0.01, 0.02, 0.04, 0.08]), [0.005, 0.04, 0.06, 0.08, 0.08])
    assert np.allclose(holm([0.9, 0.8]), [1.0, 1.0])


def test_bh_known_answer():
    assert np.allclose(bh([0.01, 0.04, 0.03, 0.005]), [0.02, 0.04, 0.04, 0.02])
    assert np.allclose(bh([0.01, 0.02, 0.03, 0.04, 0.05]), [0.05] * 5)
    assert np.allclose(bh([0.5]), [0.5])


def test_cluster_bootstrap_with_singleton_clusters_equals_protein_bootstrap():
    truth, pa, pb, ids = _toy()
    b1 = bootstrap(truth, {"a": pa, "b": pb}, ids, B=300, seed=7)
    b2 = bootstrap(truth, {"a": pa, "b": pb}, ids, B=300, seed=7, clusters={i: i for i in ids})
    for m in b1:
        for k in b1[m]:
            assert np.array_equal(b1[m][k], b2[m][k]), (m, k)


def test_cluster_bootstrap_is_wider_for_duplicated_proteins():
    truth, pa, _, base = _toy(n=20, seed=3)
    ids = [f"{i}_{r}" for i in base for r in range(5)]
    t = {f"{i}_{r}": truth[i] for i in base for r in range(5)}
    p = {f"{i}_{r}": pa[i] for i in base for r in range(5)}
    sd_prot = bootstrap(t, {"m": p}, ids, B=800, seed=1, keys=("micro_F1",))["m"]["micro_F1"].std()
    sd_clus = bootstrap(t, {"m": p}, ids, B=800, seed=1, keys=("micro_F1",),
                        clusters={f"{i}_{r}": i for i in base for r in range(5)})["m"]["micro_F1"].std()
    assert sd_clus > 1.6 * sd_prot, (sd_prot, sd_clus)  # ~sqrt(5) expected


def test_observed_delta_is_the_point_difference():
    truth, pa, pb, ids = _toy()
    d = observed_delta(truth, pa, pb, ids, "macro_union_F1")
    ref = evaluate(truth, pa, ids)["macro_union_F1"] - evaluate(truth, pb, ids)["macro_union_F1"]
    assert abs(d - ref) < 1e-15


def test_cluster_bootstrap_samples_f1_is_a_weighted_mean():
    # perfect predictions: per-protein F1 is 1 for every protein, so every resample must give exactly 1, whatever
    # the cluster sizes (the resample weights sum to N only for the protein bootstrap)
    truth, _, _, ids = _toy(n=30, seed=5)
    clusters = {i: ("big" if k < 20 else f"s{k}") for k, i in enumerate(ids)}
    bs = bootstrap(truth, {"m": truth}, ids, B=200, seed=2, keys=("samples_F1",), clusters=clusters)
    assert np.allclose(bs["m"]["samples_F1"], 1.0), bs["m"]["samples_F1"][:5]
