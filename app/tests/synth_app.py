"""Builds a small synthetic app/ folder (data + weights) in a temporary directory so the real cx_* scripts can run
end to end without CLEAN's weights, ESM-1b or any download. Nothing here reads or writes the real repository data.

The world is coherent: training ESM vectors are clustered by EC, 100.pt is those vectors passed through a
(randomly initialised, seeded) CLEAN LayerNormNet saved as split100.pth, and query vectors lie near their EC's
centre, so max-separation, calibration and risk control have something real to find."""
import csv
import gzip
import json
import pickle
from pathlib import Path

import numpy as np
import torch

AA = "ACDEFGHIKLMNPQRSTVWY"
SIZES = [1] * 10 + [2, 2, 3, 3, 2, 3, 2, 3] + [4, 5, 6, 7, 8, 9, 10, 5, 6, 7] + \
        [12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 15, 25]
TRANSFERRED, SUCCESSOR, DELETED = "2.1.1.40", "2.1.1.99", "2.2.2.30"


def _seq(rng, n=None):
    return "".join(rng.choice(list(AA), size=int(n or rng.integers(60, 160))))


def _mutate(rng, s, rate=0.2):
    s = list(s)
    for k in range(len(s)):
        if rng.random() < rate:
            s[k] = AA[int(rng.integers(0, 20))]
    return "".join(s)


def _write_tsv(path, header, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(header)
        w.writerows(rows)


def build(root, seed=0):
    import sys
    from sklearn.mixture import GaussianMixture
    app_src = Path(__file__).resolve().parents[1] / "src"
    if str(app_src) not in sys.path:
        sys.path.insert(0, str(app_src))
    from CLEAN.model import LayerNormNet
    from CLEAN.utils import get_ec_id_dict

    rng = np.random.default_rng(seed)
    root = Path(root)
    data = root / "data"
    cx = data / "cleanx"
    for d in (data / "pretrained", data / "esm_data", data / "datasets", cx / "downloads", cx / "sets", cx / "fasta",
              cx / "identity", cx / "cache", root / "results" / "cleanx" / "preds"):
        d.mkdir(parents=True, exist_ok=True)

    # vocabulary: 40 training ECs (long tail); TRANSFERRED is a common one, DELETED is rare
    ecs = [f"{1 + k % 3}.{1 + (k // 3) % 2}.{1 + (k // 6) % 2}.{10 + k}" for k in range(len(SIZES))]
    ecs[-3], ecs[20] = TRANSFERRED, DELETED
    centers = {e: rng.normal(0, 1.0, 1280) for e in ecs}
    nise = set(ecs[-6:-3])  # common ECs whose members form two separate sub-centres
    offset = {e: rng.normal(0, 0.8, 1280) for e in nise}

    entries, esm_train = [], {}
    for e, n in zip(ecs, SIZES):
        for j in range(n):
            acc = f"T{len(entries):05d}"
            v = centers[e] + rng.normal(0, 0.25, 1280)
            if e in nise:
                v = v + (offset[e] if j % 2 else -offset[e])
            entries.append([acc, e, _seq(rng)])
            esm_train[acc] = v
    for k in range(6):  # multi-EC training entries
        entries[40 + 7 * k][1] += ";" + ecs[-1 - k]
    _write_tsv(data / "split100.csv", ["Entry", "EC number", "Sequence"], entries)

    torch.manual_seed(seed)
    model = LayerNormNet(512, 128, "cpu", torch.float32)
    model.eval()
    torch.save(model.state_dict(), data / "pretrained" / "split100.pth")
    _, ec_id = get_ec_id_dict(str(data / "split100.csv"))
    rows = []
    with torch.no_grad():
        for e in ec_id:
            for acc in sorted(ec_id[e]):
                rows.append(model(torch.tensor(esm_train[acc], dtype=torch.float32)[None])[0])
    torch.save(torch.stack(rows), data / "pretrained" / "100.pt")

    gmms = []
    x = np.r_[rng.normal(3.0, 1.0, 400), rng.normal(8.0, 1.5, 400)].reshape(-1, 1)
    for k in range(3):
        gmms.append(GaussianMixture(n_components=2, random_state=k).fit(x))
    with open(data / "pretrained" / "gmm_ensumble.pkl", "wb") as fh:
        pickle.dump(gmms, fh)

    counts = {}
    for _, ecf, _ in entries:
        for e in ecf.split(";"):
            counts[e] = counts.get(e, 0) + 1
    tier = lambda n: "singleton" if n == 1 else ("very_rare" if n <= 3 else ("rare" if n <= 10 else "common"))
    with open(data / "ec_train_frequency.csv", "w", encoding="utf-8", newline="") as fh:
        fh.write("ec,train_count,tier\n" + "".join(f"{e},{n},{tier(n)}\n" for e, n in counts.items()))

    with open(cx / "downloads" / "enzyme.dat", "w", encoding="utf-8") as fh:
        fh.write("CC   Release of 02-Sep-2026\n//\n")
        fh.write(f"ID   {TRANSFERRED}\nDE   Transferred entry: {SUCCESSOR}.\n//\n")
        fh.write(f"ID   {DELETED}\nDE   Deleted entry.\n//\n")
        for e in ecs:
            if e not in (TRANSFERRED, DELETED):
                fh.write(f"ID   {e}\nDE   synthetic enzyme {e}.\n//\n")

    weights = np.array([n for n in SIZES], float)
    weights /= weights.sum()
    esm_query, seqs = {}, {}

    def query(prefix, n, unseen_frac=0.0, mixed_frac=0.0, successor_frac=0.0):
        ids, lab = [], {}
        for k in range(n):
            q = f"{prefix}{k:04d}"
            u = rng.random()
            e = ecs[int(rng.choice(len(ecs), p=weights))]
            if u < unseen_frac:  # pure unseen: an EC number CLEAN never saw, new region of space
                lab[q] = [f"3.3.3.{90 + k % 5}"]
                v = rng.normal(0, 1.0, 1280)
            elif u < unseen_frac + mixed_frac:  # mixed: one seen EC + one unseen EC
                lab[q] = [e, f"3.3.3.{90 + k % 5}"]
                v = centers[e] + rng.normal(0, 0.5, 1280)
            elif u < unseen_frac + mixed_frac + successor_frac:  # labelled with the successor of a transferred EC
                lab[q] = [SUCCESSOR]
                v = centers[TRANSFERRED] + rng.normal(0, 0.5, 1280)
            else:
                lab[q] = [e]
                v = centers[e] + rng.normal(0, 0.5, 1280) + (offset[e] if e in nise and k % 2 else 0)
            esm_query[q] = v
            seqs[q] = _seq(rng)
            ids.append(q)
        return ids, lab

    sets = {"new392": query("N", 50), "price149": query("P", 30),
            "sp_temporal_cal": query("C", 80, 0.15, 0.1, 0.05), "sp_temporal_test": query("E", 110, 0.15, 0.1, 0.05)}
    cal_ids = sets["sp_temporal_cal"][0]
    for k, q in enumerate(sets["sp_temporal_test"][0][:12]):  # near-copies of CAL sequences
        seqs[q] = _mutate(rng, seqs[cal_ids[k]], 0.15)
    non_ids = [f"X{k:04d}" for k in range(40)]
    for q in non_ids:
        esm_query[q] = rng.normal(0, 1.0, 1280) * 1.5
        seqs[q] = _seq(rng)

    paths = {"new392": data / "new.csv", "price149": data / "datasets" / "price.csv"}
    for name, (ids, lab) in sets.items():
        p = paths.get(name, cx / "sets" / f"{name}.csv")
        _write_tsv(p, ["Entry", "EC number", "Sequence"], [[i, ";".join(lab[i]), seqs[i]] for i in ids])
    _write_tsv(cx / "sets" / "sp_nonenzyme.csv", ["Entry", "EC number", "Sequence"], [[i, "", seqs[i]] for i in non_ids])

    for q, v in esm_query.items():
        torch.save({"label": q, "mean_representations": {33: torch.tensor(v, dtype=torch.float32)}},
                   data / "esm_data" / f"{q}.pt")

    from cleanx import labels as L
    from cleanx.nomenclature import forward_counts, forward_map, parse_enzyme_dat
    fmap = forward_map(parse_enzyme_dat(cx / "downloads" / "enzyme.dat"))
    cur = forward_counts(counts, fmap)
    meta = []
    for name in ("sp_temporal_cal", "sp_temporal_test"):
        ids, lab = sets[name]
        for q in ids:
            meta.append([q, name, L.query_tier(lab[q], cur), any(cur.get(e, 0) == 0 for e in lab[q]),
                         f"202{2 + int(rng.integers(0, 4))}-0{1 + int(rng.integers(0, 8))}-15", "Synthetica",
                         "Evidence at protein level" if rng.random() < 0.6 else "Evidence at transcript level",
                         len(seqs[q]), "; ".join(lab[q]), ""])
    for q in non_ids:
        meta.append([q, "sp_nonenzyme", "non_enzyme", False, "2024-03-27", "Synthetica",
                     "Evidence at protein level", len(seqs[q]), "", ""])
    with open(cx / "sets" / "sp_temporal_meta.csv", "w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, lineterminator="\n")
        w.writerow(["accession", "set", "tier", "novel_EC", "date_created", "organism_name", "protein_existence",
                    "len", "ec", "xref_pfam"])
        w.writerows(meta)
    with open(cx / "sets" / "sp_temporal_manifest.json", "w", encoding="utf-8") as fh:
        json.dump({"uniprot": {"uniprot_release": "2026_03", "uniprot_release_date": "02-September-2026"}}, fh)

    for name, (ids, _) in sets.items():  # identity to training (for the identity-stratified subsets)
        with open(cx / "identity" / f"{name}.csv", "w", encoding="utf-8", newline="") as fh:
            fh.write("Entry,max_pid_short,max_pid_aln,top_hit,top_score,top_hit_pid_short,identity_bin,tool\n")
            for q in ids:
                pid = float(rng.uniform(0.1, 0.99))
                fh.write(f"{q},{pid:.3f},{pid:.3f},,0,0,,synthetic\n")

    fam = {e: f"PF{10000 + k:05d}" for k, e in enumerate(ecs)}
    clan = {f: (f"CL{k % 7:04d}" if k % 3 else "") for k, f in enumerate(sorted(fam.values()))}
    nise_fam = {e: f"PF{20000 + k:05d}" for k, e in enumerate(sorted(nise))}
    with gzip.open(cx / "downloads" / "uniprot_sp_ec_all.tsv.gz", "wt", encoding="utf-8") as fh:
        fh.write("Entry\tEC number\tDate of creation\tPfam\tLength\n")
        per_ec = {}
        for acc, ecf, s in entries:
            e = ecf.split(";")[0]
            j = per_ec.get(e, 0)
            per_ec[e] = j + 1
            pf = nise_fam[e] if e in nise and j % 2 else fam[e]
            fh.write(f"{acc}\t{ecf.replace(';', '; ')}\t2020-01-01\t{pf};\t{len(s)}\n")
    with gzip.open(cx / "downloads" / "Pfam-A.clans.tsv.gz", "wt", encoding="utf-8") as fh:
        for f, c in clan.items():
            fh.write(f"{f}\t{c}\t{c or ''}\t{f}\tsynthetic family\n")
        for f in nise_fam.values():
            fh.write(f"{f}\tCL9999\tCL9999\t{f}\tsynthetic family (other clan)\n")
    return root
