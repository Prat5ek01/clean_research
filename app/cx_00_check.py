"""E0: environment + exact-replica checks. Every later experiment depends on these passing."""
import hashlib
import sys
import time
import warnings

warnings.filterwarnings("ignore")

from cleanx import config

config.enter_app_dir()
results = []


def check(name, ok, detail=""):
    results.append(ok)
    print(f"[{'PASS' if ok else 'FAIL'}] {name}" + (f" -- {detail}" if detail else ""), flush=True)
    return ok


def md5(path, n=1 << 22):
    h = hashlib.md5()
    with open(path, "rb") as fh:
        while True:
            b = fh.read(n)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


def main():
    t0 = time.time()
    import numpy, pandas, scipy, sklearn, torch
    print(f"python {sys.version.split()[0]} | numpy {numpy.__version__} | pandas {pandas.__version__} | "
          f"scipy {scipy.__version__} | sklearn {sklearn.__version__} | torch {torch.__version__} | "
          f"threads {torch.get_num_threads()}")
    try:
        import Bio
        print(f"biopython {Bio.__version__} (identity fallback available)")
    except ImportError:
        print("biopython not installed (only needed for cx_04_identity.py --tool python)")

    needed = [config.BACKBONES["split100"][k] for k in ("csv", "emb", "pth")] + [
        config.GMM_PKL, config.LEGACY_SETS["new392"], config.LEGACY_SETS["price149"],
        config.APP_DIR / "results" / "new_maxsep.csv", config.APP_DIR / "results" / "inputs" / "price_maxsep.csv"]
    missing = [str(p) for p in needed if not p.is_file()]
    if not check("required files present", not missing, "; ".join(missing)):
        sys.exit(1)
    for p in needed[1:4]:
        print(f"       md5 {p.name}: {md5(p)}")

    from cleanx import labels as L
    from cleanx.index import TrainIndex
    from cleanx.metrics import evaluate
    from cleanx.scoring import centroid_distances, maxsep_predict

    idx = TrainIndex("split100")
    check("split100 index", len(idx.ecs) == 5242 and idx.E.shape == (241025, 128),
          f"{len(idx.ecs)} ECs, embeddings {idx.E.shape}")

    ids, truth, _ = L.read_dataset(config.LEGACY_SETS["new392"])
    Q = idx.embed(ids)
    preds = maxsep_predict(centroid_distances(Q, idx), idx.ecs, ids)
    ref = L.read_pred_csv(config.APP_DIR / "results" / "new_maxsep.csv")
    ec_mis = sum([e for e, _ in preds[q]] != [e for e, _ in ref[q]] for q in ids)
    dmax = max(abs(round(d1, 4) - d2) for q in ids for (_, d1), (_, d2) in zip(preds[q], ref[q]))
    check("New-392 max-separation replica vs results/new_maxsep.csv", ec_mis == 0 and dmax <= 1.5e-4,
          f"EC-list mismatches {ec_mis}/392, max |d diff| {dmax:.5f}")

    r = evaluate(truth, L.labels_only(preds), ids)
    ok = (abs(r["weighted_P"] - 0.5955) < 1e-4 and abs(r["weighted_R"] - 0.4791) < 1e-4
          and abs(r["weighted_F1"] - 0.4967) < 1e-4)
    check("CLEAN published metric (support-weighted) reproduced", ok,
          f"P {r['weighted_P']:.4f} R {r['weighted_R']:.4f} F1 {r['weighted_F1']:.4f} (paper 0.596/0.479/0.497)")
    from CLEAN.evaluate import get_eval_metrics, get_pred_probs, get_true_labels, get_pred_labels
    pl = get_pred_labels("results/new", pred_type="_maxsep")
    pp = get_pred_probs("results/new", pred_type="_maxsep")
    tl, al = get_true_labels("data/new")
    pre, rec, f1, auc, _ = get_eval_metrics(pl, pp, tl, al)
    check("cleanx metric == CLEAN.evaluate.get_eval_metrics",
          abs(pre - r["weighted_P"]) < 1e-9 and abs(rec - r["weighted_R"]) < 1e-9 and abs(f1 - r["weighted_F1"]) < 1e-9,
          f"CLEAN P {pre:.4f} R {rec:.4f} F1 {f1:.4f} AUC {auc:.4f}")
    check("blind spot is visible in micro metrics", r["micro_P"] < r["weighted_P"],
          f"micro P {r['micro_P']:.4f}; {r['fp_outside_true_labelset']}/{r['fp']} false positives are ECs "
          f"outside the test label set")

    from cleanx.calibration import gmm_confidence, load_gmms
    pids, _, _ = L.read_dataset(config.LEGACY_SETS["price149"])
    Dp = centroid_distances(idx.embed(pids), idx)
    ppred = maxsep_predict(Dp, idx.ecs, pids)
    gmms = load_gmms()
    pref = L.read_pred_csv(config.APP_DIR / "results" / "inputs" / "price_maxsep.csv")
    mis, cmax = 0, 0.0
    for q in pids:
        mine = [e for e, _ in ppred[q]]
        if mine != [e for e, _ in pref[q]]:
            mis += 1
            continue
        conf = gmm_confidence([d for _, d in ppred[q]], gmms)
        cmax = max(cmax, max(abs(round(c, 4) - s) for c, (_, s) in zip(conf, pref[q])))
    check("Price-149 GMM-confidence replica vs results/inputs/price_maxsep.csv", mis == 0 and cmax <= 1.5e-4,
          f"EC-list mismatches {mis}/149, max |conf diff| {cmax:.5f}")

    freq = config.DATA_DIR / "ec_train_frequency.csv"
    if freq.is_file():
        import pandas as pd
        f = pd.read_csv(freq)
        mine = dict(zip(idx.ecs, idx.tiers))
        bad = sum(mine.get(e) != t for e, t in zip(f.ec, f.tier))
        check("tier thresholds match data/ec_train_frequency.csv", bad == 0 and len(f) == len(idx.ecs),
              f"{bad} mismatches over {len(f)} ECs")

    sp70 = config.BACKBONES["split70"]
    if all(sp70[k].is_file() for k in ("csv", "emb", "pth")):
        try:
            i70 = TrainIndex("split70")
            check("split70 backbone loads (second backbone for robustness)", True,
                  f"{len(i70.ecs)} ECs, embeddings {i70.E.shape}")
        except Exception as e:
            check("split70 backbone loads", False, str(e))
    n_fail = results.count(False)
    print(f"\nCX00 SUMMARY: {results.count(True)}/{len(results)} checks passed in {time.time() - t0:.0f}s")
    sys.exit(1 if n_fail else 0)


if __name__ == "__main__":
    main()
