"""Max sequence identity of each evaluation protein to CLEAN's split100 training set, plus the
homology-transfer baseline (EC of the best alignment hit).

--tool python : 3-mer prefilter + Smith-Waterman (Biopython), runs anywhere, slower
--tool m8     : parse MMseqs2 results (recommended for the paper; see RUNBOOK for the command)
--write-fasta : write split100 + query FASTAs for running MMseqs2 elsewhere
--pair Q:T    : max identity of every protein of set Q to set T (e.g. sp_temporal_test:sp_temporal_cal), written
                to data/cleanx/identity/Q_vs_T.csv; used by cx_05 for the CAL-disjoint TEST sensitivity analysis
                (sequences only: no prediction or label of either set is used)"""
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
    ap.add_argument("--pair", default=None, help="QUERY_SET:TARGET_SET, e.g. sp_temporal_test:sp_temporal_cal")
    args = ap.parse_args()
    config.ensure_dirs()
    t0 = time.time()

    if args.pair:
        qn, tn = args.pair.split(":")
        q_ids, _, q_seqs = L.load_named_dataset(qn)
        t_ids, _, t_seqs = L.load_named_dataset(tn)
        if args.limit:
            q_ids = q_ids[:args.limit]
        overlap = set(q_ids) & set(t_ids)
        if overlap:
            raise SystemExit(f"{qn} and {tn} share {len(overlap)} accessions; they must be disjoint")
        print(f"== {qn} vs {tn}: {len(q_ids)} x {len(t_ids)} proteins")
        rows = I.python_search(q_ids, [q_seqs[i] for i in q_ids], t_ids, [t_seqs[i] for i in t_ids], topn=args.topn)
        df = pd.DataFrame(rows).rename(columns={"top_hit": "top_hit_in_target"})
        df["identity_bin"] = df.max_pid_short.map(I.identity_bin)
        out = config.CX_IDENTITY / f"{qn}_vs_{tn}{'_smoke' if args.limit else ''}.csv"
        df.to_csv(out, index=False)
        print(df.identity_bin.value_counts().reindex(I.IDENTITY_BINS, fill_value=0).to_string())
        print(f"  {int((df.max_pid_short >= 0.5).sum())} of {len(df)} {qn} proteins are >= 50% identical to a {tn} protein")
        print(f"\nCX04 PAIR DONE in {time.time() - t0:.0f}s -> {out}")
        return

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
