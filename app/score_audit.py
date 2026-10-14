"""Score CLEAN predictions by EC depth and training-frequency tier.

Defaults to results/new_maxsep.csv; use --preds and --tag for other predictions.

See docs/code_reference.md#score-audit."""

import argparse
import os
import sys
import pandas as pd


# ----------------------------------------------------------------------
# loading
# ----------------------------------------------------------------------

def load_truth(path):
    """Entry -> set of true EC numbers."""
    df = pd.read_csv(path, sep="\t")
    if "Entry" not in df.columns or "EC number" not in df.columns:
        sys.exit(f"[!] {path} does not have Entry / 'EC number' columns. Found: {list(df.columns)}")
    truth = {}
    for _, row in df.iterrows():
        ecs = {
            e.strip()
            for e in str(row["EC number"]).split(";")
            if e.strip() and e.strip().lower() != "nan"
        }
        truth[str(row["Entry"]).strip()] = ecs
    return truth


def load_preds(path):
    """
    Entry -> (set of predicted ECs, {ec: score}).

    CLEAN output lines look like:
        Q9RYA6,EC:5.1.1.20/7.4553
        Q5TZ07,EC:3.6.1.43/8.0610,EC:3.1.3.4/8.0627
    The trailing number is a DISTANCE when infer_maxsep is called without gmm=,
    and a GMM CONFIDENCE (0-1) when called with it. Semantics differ - see notes.
    """
    preds, scores = {}, {}
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            parts = line.split(",")
            entry = parts[0].strip()
            ec_set, sc = set(), {}
            for chunk in parts[1:]:
                chunk = chunk.strip()
                if not chunk:
                    continue
                body = chunk[3:] if chunk.startswith("EC:") else chunk
                if "/" in body:
                    ec, val = body.rsplit("/", 1)
                    ec = ec.strip()
                    try:
                        sc[ec] = float(val)
                    except ValueError:
                        pass
                else:
                    ec = body.strip()
                if ec:
                    ec_set.add(ec)
            preds[entry] = ec_set
            scores[entry] = sc
    return preds, scores


def load_freq(path):
    df = pd.read_csv(path)
    return dict(zip(df["ec"], df["train_count"])), dict(zip(df["ec"], df["tier"]))


# ----------------------------------------------------------------------
# metrics
# ----------------------------------------------------------------------

def truncate(ec, depth):
    """3.1.1.3 at depth 2 -> '3.1'"""
    return ".".join(ec.split(".")[:depth])


def micro_prf(pairs, depth=4):
    """
    Micro-averaged multi-label precision / recall / F1.
    pairs: iterable of (true_set, pred_set)
    """
    tp = fp = fn = 0
    for true_set, pred_set in pairs:
        t = {truncate(e, depth) for e in true_set}
        p = {truncate(e, depth) for e in pred_set}
        tp += len(t & p)
        fp += len(p - t)
        fn += len(t - p)
    prec = tp / (tp + fp) if (tp + fp) else 0.0
    rec = tp / (tp + fn) if (tp + fn) else 0.0
    f1 = 2 * prec * rec / (prec + rec) if (prec + rec) else 0.0
    return prec, rec, f1, tp, fp, fn


def sample_tier(true_ecs, counts, tiers):
    """
    A sample's difficulty is driven by its RAREST true EC.
    An EC absent from the training index is 'unseen'.
    """
    best_tier, best_count = "unseen", None
    order = {"unseen": 0, "singleton": 1, "very_rare": 2, "rare": 3, "common": 4}
    for ec in true_ecs:
        c = counts.get(ec)
        t = tiers.get(ec, "unseen")
        if best_count is None or order[t] < order[best_tier]:
            best_tier, best_count = t, c
    return best_tier


# ----------------------------------------------------------------------
# reporting
# ----------------------------------------------------------------------

def line(width=74):
    print("-" * width)


def report_depth(pairs, tag):
    print()
    line()
    print(f"  DEPTH BREAKDOWN  ({tag})   -- does it get the family but miss the reaction?")
    line()
    print(f"{'depth':<24}{'precision':>11}{'recall':>10}{'F1':>10}{'TP':>8}{'FP':>7}{'FN':>7}")
    labels = {1: "1  (broad class)", 2: "2  (sub-class)",
              3: "3  (sub-sub-class)", 4: "4  (exact reaction)"}
    rows = []
    for d in (1, 2, 3, 4):
        p, r, f, tp, fp, fn = micro_prf(pairs, depth=d)
        print(f"{labels[d]:<24}{p:>11.3f}{r:>10.3f}{f:>10.3f}{tp:>8}{fp:>7}{fn:>7}")
        rows.append({"tag": tag, "depth": d, "precision": p, "recall": r,
                     "f1": f, "tp": tp, "fp": fp, "fn": fn})
    return rows


def report_tier(by_tier, tag, depth=4):
    print()
    line()
    print(f"  RARITY BREAKDOWN  ({tag}, depth={depth})  -- how far does it fall on rare classes?")
    line()
    print(f"{'tier':<14}{'n':>6}{'precision':>11}{'recall':>10}{'F1':>10}{'TP':>8}{'FP':>7}{'FN':>7}")
    rows = []
    for t in ("common", "rare", "very_rare", "singleton", "unseen"):
        pairs = by_tier.get(t, [])
        if not pairs:
            continue
        p, r, f, tp, fp, fn = micro_prf(pairs, depth=depth)
        print(f"{t:<14}{len(pairs):>6}{p:>11.3f}{r:>10.3f}{f:>10.3f}{tp:>8}{fp:>7}{fn:>7}")
        rows.append({"tag": tag, "tier": t, "n": len(pairs), "depth": depth,
                     "precision": p, "recall": r, "f1": f})
    return rows


def report_cross(by_tier, tag):
    print()
    line()
    print(f"  CORE AUDIT TABLE  ({tag})  --  F1 by rarity tier x digit depth")
    line()
    print(f"{'tier':<14}{'n':>6}{'d1':>9}{'d2':>9}{'d3':>9}{'d4':>9}")
    rows = []
    for t in ("common", "rare", "very_rare", "singleton", "unseen"):
        pairs = by_tier.get(t, [])
        if not pairs:
            continue
        f1s = [micro_prf(pairs, depth=d)[2] for d in (1, 2, 3, 4)]
        print(f"{t:<14}{len(pairs):>6}" + "".join(f"{v:>9.3f}" for v in f1s))
        rows.append({"tag": tag, "tier": t, "n": len(pairs),
                     "f1_d1": f1s[0], "f1_d2": f1s[1],
                     "f1_d3": f1s[2], "f1_d4": f1s[3]})
    return rows


def report_confidence(truth, preds, scores, tag):
    """
    Is the model confident exactly when it is wrong?
    NOTE: the score is a DISTANCE (maxsep without gmm) or a CONFIDENCE 0-1
    (CLEAN_infer_fasta.py, which passes gmm=). Read the direction accordingly.
    """
    correct, wrong = [], []
    for entry, pred_set in preds.items():
        true_set = truth.get(entry)
        if not true_set:
            continue
        for ec in pred_set:
            s = scores.get(entry, {}).get(ec)
            if s is None:
                continue
            (correct if ec in true_set else wrong).append(s)
    if not correct and not wrong:
        return []
    print()
    line()
    print(f"  SCORE SEPARATION  ({tag})  -- confident when right, unsure when wrong?")
    line()
    def stat(v):
        return (f"n={len(v):<6} mean={sum(v)/len(v):<9.4f}"
                f"min={min(v):<9.4f}max={max(v):<9.4f}") if v else "n=0"
    print(f"  correct predictions : {stat(correct)}")
    print(f"  wrong   predictions : {stat(wrong)}")
    print("  (maxsep w/o gmm -> distance; CLEAN_infer_fasta.py -> GMM confidence 0-1)")
    return [{"tag": tag, "group": "correct", "n": len(correct),
             "mean": sum(correct)/len(correct) if correct else None},
            {"tag": tag, "group": "wrong", "n": len(wrong),
             "mean": sum(wrong)/len(wrong) if wrong else None}]


# ----------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--truth", default="data/new.csv")
    ap.add_argument("--preds", default="results/new_maxsep.csv")
    ap.add_argument("--freq",  default="data/ec_train_frequency.csv")
    ap.add_argument("--tag",   default="maxsep")
    ap.add_argument("--outdir", default="audit_out")
    args = ap.parse_args()

    for path in (args.truth, args.preds, args.freq):
        if not os.path.exists(path):
            sys.exit(f"[!] missing file: {path}\n"
                     f"    (run with --preds/--truth/--freq to point at the right paths)")

    truth = load_truth(args.truth)
    preds, scores = load_preds(args.preds)
    counts, tiers = load_freq(args.freq)

    scored = [e for e in truth if e in preds]
    missing = [e for e in truth if e not in preds]

    print()
    line()
    print(f"  CLEAN AUDIT  --  {args.tag}")
    line()
    print(f"  truth file        : {args.truth}   ({len(truth)} sequences)")
    print(f"  predictions file  : {args.preds}   ({len(preds)} sequences)")
    print(f"  scored            : {len(scored)}")
    if missing:
        print(f"  !! no prediction for {len(missing)} sequences (counted as all-FN below)")

    # include missing-prediction sequences as empty predictions - honest recall
    pairs_all = [(truth[e], preds.get(e, set())) for e in truth]

    by_tier = {}
    for e in truth:
        t = sample_tier(truth[e], counts, tiers)
        by_tier.setdefault(t, []).append((truth[e], preds.get(e, set())))

    d_rows = report_depth(pairs_all, args.tag)
    t_rows = report_tier(by_tier, args.tag, depth=4)
    c_rows = report_cross(by_tier, args.tag)
    s_rows = report_confidence(truth, preds, scores, args.tag)

    os.makedirs(args.outdir, exist_ok=True)
    pd.DataFrame(d_rows).to_csv(f"{args.outdir}/{args.tag}_by_depth.csv", index=False)
    pd.DataFrame(t_rows).to_csv(f"{args.outdir}/{args.tag}_by_tier.csv", index=False)
    pd.DataFrame(c_rows).to_csv(f"{args.outdir}/{args.tag}_tier_x_depth.csv", index=False)
    if s_rows:
        pd.DataFrame(s_rows).to_csv(f"{args.outdir}/{args.tag}_score_separation.csv", index=False)

    print()
    line()
    print(f"  wrote CSVs to {os.path.abspath(args.outdir)}")
    line()
    print()


if __name__ == "__main__":
    main()