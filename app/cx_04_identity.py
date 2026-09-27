"""Max sequence identity of each evaluation protein to CLEAN's split100 training set, plus the
homology-transfer baseline (EC of the best alignment hit).

--tool python : 3-mer prefilter + Smith-Waterman (Biopython), runs anywhere, slower
--tool m8     : parse MMseqs2 results (recommended for the paper; see RUNBOOK for the command)
--write-fasta : write split100 + query FASTAs for running MMseqs2 elsewhere"""
import argparse
import time
import warnings

warnings.filterwarnings("ignore")

from cleanx import config

config.enter_app_dir()

import pandas as pd

from cleanx import identity as I
from cleanx import labels as L
from cleanx.index import TrainIndex


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sets", nargs="+", default=["new392", "price149", "sp_temporal_cal", "sp_temporal_test",
                                                  "sp_nonenzyme"])
    ap.add_argument("--tool", choices=["python", "m8"], default="python")
    ap.add_argument("--m8-dir", default=str(config.CX_IDENTITY / "mmseqs"),
                    help="folder holding <set>.m8 files from MMseqs2")
    ap.add_argument("--topn", type=int, default=30)
    ap.add_argument("--limit", type=int, default=0, help="first N proteins per set (smoke)")
    ap.add_argument("--write-fasta", action="store_true")
    args = ap.parse_args()
    config.ensure_dirs()
    t0 = time.time()

    train = pd.read_csv(config.BACKBONES["split100"]["csv"], sep="\t", dtype=str, keep_default_na=False)
    t_ids, t_seqs = list(train.Entry), list(train.Sequence)
    id_ec = dict(zip(train.Entry, train["EC number"]))
    if args.write_fasta:
        L.write_fasta(config.CX_FASTA / "split100.fasta", t_ids, dict(zip(t_ids, t_seqs)))
        for name in args.sets:
            ids, _, seqs = L.load_named_dataset(name)
            L.write_fasta(config.CX_FASTA / f"{name}.fasta", ids, seqs)
        print(f"FASTA files written to {config.CX_FASTA}")
        return

    for name in args.sets:
        try:
            ids, truth, seqs = L.load_named_dataset(name)
        except FileNotFoundError as e:
            print(f"  skip {name}: {e}")
            continue
        if args.limit:
            ids = ids[:args.limit]
        print(f"== {name}: {len(ids)} proteins")
        if args.tool == "python":
            rows = I.python_search(ids, [seqs[i] for i in ids], t_ids, t_seqs, topn=args.topn)
        else:
            m8 = f"{args.m8_dir}/{name}.m8"
            rows = I.parse_m8(m8, ids)
        for r in rows:
            r["top_hit_EC"] = id_ec.get(r["top_hit"], "")
            r["identity_bin"] = I.identity_bin(r["max_pid_short"])
            r["tool"] = args.tool
        suffix = "_smoke" if args.limit else ""
        df = pd.DataFrame(rows)
        df.to_csv(config.CX_IDENTITY / f"{name}{suffix}.csv", index=False)
        preds = {r["Entry"]: ([(e, r["top_score"]) for e in L.split_ecs(r["top_hit_EC"])] if r["top_hit"] else [])
                 for r in rows}
        L.write_pred_csv(config.CX_PREDS / f"{name}_homology_tophit{suffix}.csv", preds)
        print(df.identity_bin.value_counts().reindex(I.IDENTITY_BINS, fill_value=0).to_string())
        if truth and any(truth.values()):
            ok = sum(bool(set(L.split_ecs(r["top_hit_EC"])) & set(truth[r["Entry"]])) for r in rows)
            print(f"  homology top-hit exact-EC agreement: {ok}/{len(rows)}")
    print(f"\nCX04 DONE in {time.time() - t0:.0f}s -> {config.CX_IDENTITY}")


if __name__ == "__main__":
    main()
