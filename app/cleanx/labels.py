import csv
import json
import re
from pathlib import Path

import pandas as pd

from . import config

EC_COMPLETE = re.compile(r"^\d+\.\d+\.\d+\.n?\d+$")


def split_ecs(field):
    """Split an EC field ('1.1.1.1;2.3.4.5' or UniProt '1.1.1.1; 2.3.4.-') into unique ECs, order kept."""
    out = []
    for tok in str(field).split(";"):
        t = tok.strip()
        if t and t.lower() != "nan" and t not in out:
            out.append(t)
    return out


def complete_ecs(ecs):
    return [e for e in ecs if EC_COMPLETE.match(e)]


def truncate(ec, depth):
    return ".".join(ec.split(".")[:depth])


def read_dataset(path):
    """Read a CLEAN-format tab-separated file (Entry, EC number[, Sequence]).

    Returns (ids, labels, seqs) where labels maps id -> list of ECs."""
    df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    missing = {"Entry", "EC number"} - set(df.columns)
    if missing:
        raise ValueError(f"{path}: missing columns {sorted(missing)}; found {list(df.columns)}")
    ids = [s.strip() for s in df["Entry"]]
    if len(set(ids)) != len(ids):
        raise ValueError(f"{path}: duplicate Entry ids")
    labels = {i: split_ecs(e) for i, e in zip(ids, df["EC number"])}
    seqs = dict(zip(ids, df["Sequence"])) if "Sequence" in df.columns else {}
    return ids, labels, seqs


def load_named_dataset(name):
    path = config.dataset_csv(name)
    if not Path(path).is_file():
        raise FileNotFoundError(f"dataset '{name}' not found at {path} (build it with cx_02_build_temporal.py)")
    return read_dataset(path)


def write_dataset(path, ids, labels, seqs):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["Entry", "EC number", "Sequence"])
        for i in ids:
            w.writerow([i, ";".join(labels.get(i, [])), seqs[i]])


def write_fasta(path, ids, seqs):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="\n") as fh:
        for i in ids:
            fh.write(f">{i}\n{seqs[i]}\n")


def read_fasta(path):
    recs, name, buf = [], None, []
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    recs.append((name, "".join(buf)))
                name, buf = line[1:].split()[0], []
            else:
                buf.append(line)
    if name is not None:
        recs.append((name, "".join(buf)))
    return recs


def read_pred_csv(path):
    """CLEAN prediction file: 'id,EC:1.2.3.4/score,...'. Returns id -> [(ec, score or None)]."""
    preds = {}
    with open(path, encoding="utf-8") as fh:
        for row in csv.reader(fh):
            if not row or not row[0].strip():
                continue
            items = []
            for cell in row[1:]:
                cell = cell.strip()
                if not cell:
                    continue
                body = cell[3:] if cell.startswith("EC:") else cell
                ec, sep, sc = body.rpartition("/")
                if not sep:
                    ec, sc = body, None
                items.append((ec.strip(), float(sc) if sc not in (None, "") else None))
            preds[row[0].strip()] = items
    return preds


def write_pred_csv(path, preds):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, lineterminator="\n")
        for q, items in preds.items():
            w.writerow([q] + [f"EC:{ec}/{sc:.4f}" if sc is not None else f"EC:{ec}" for ec, sc in items])


def labels_only(preds):
    return {q: [it[0] if isinstance(it, (tuple, list)) else it for it in items] for q, items in preds.items()}


def query_tier(true_ecs, ec_counts):
    tiers = [config.tier_of_count(ec_counts.get(e, 0)) for e in true_ecs] or ["unseen"]
    return min(tiers, key=lambda t: config.TIER_RANK[t])


def harmonize(pred_ecs, truth_set, fmap):
    """Map training-era (2022) EC predictions onto current EC numbers.

    A prediction becomes the truth label it now maps to when any successor matches;
    otherwise its first successor (or itself if deleted), so one prediction stays one label."""
    out = []
    for ec in pred_ecs:
        succ = fmap.get(ec, [ec]) if fmap else [ec]
        hit = [s for s in succ if s in truth_set]
        lab = hit[0] if hit else (succ[0] if succ else ec)
        if lab not in out:
            out.append(lab)
    return out


def harmonize_all(preds_labels, truth, fmap):
    return {q: harmonize(p, set(truth.get(q, [])), fmap) for q, p in preds_labels.items()}


def save_json(obj, path):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(obj, fh, indent=2, default=str)


def load_json(path):
    with open(path, encoding="utf-8") as fh:
        return json.load(fh)
