"""E-data: build the leakage-controlled temporal benchmark (Swiss-Prot entries created after the
CLEAN training snapshot) and the non-enzyme open-set set; download ENZYME nomenclature and,
with --with-pfam, the full Swiss-Prot EC/Pfam table used by E1 (label drift) and E3."""
import argparse
import time
import warnings

warnings.filterwarnings("ignore")

from cleanx import config

config.enter_app_dir()

import numpy as np
import pandas as pd

from cleanx import downloads as dl
from cleanx import labels as L
from cleanx.nomenclature import forward_counts, forward_map, parse_enzyme_dat

ENZYME_CLASS_KW = {"KW-0560", "KW-0808", "KW-0378", "KW-0456", "KW-0413", "KW-0436", "KW-1278"}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--snapshot", default=config.SNAPSHOT_DATE)
    ap.add_argument("--min-len", type=int, default=config.MIN_LEN)
    ap.add_argument("--max-len", type=int, default=config.MAX_LEN)
    ap.add_argument("--cal-frac", type=float, default=0.4)
    ap.add_argument("--nonenz-max", type=int, default=2000)
    ap.add_argument("--seed", type=int, default=config.SEED)
    ap.add_argument("--with-pfam", action="store_true", help="also download all Swiss-Prot EC entries + Pfam (~8 min)")
    ap.add_argument("--smoke", action="store_true", help="first 400 UniProt records only; writes smk_* sets")
    ap.add_argument("--force", action="store_true", help="re-download even if cached")
    args = ap.parse_args()
    config.ensure_dirs()
    t0 = time.time()
    prefix = "smk" if args.smoke else "sp"

    info_enz = dl.plain_file(dl.ENZYME_DAT_URL, config.CX_DOWNLOADS / "enzyme.dat", force=args.force)
    nom = parse_enzyme_dat(config.CX_DOWNLOADS / "enzyme.dat")
    fmap = forward_map(nom)
    print(f"ENZYME {nom['release']}: {len(nom['active'])} active, {len(nom['transfers'])} transferred, "
          f"{len(nom['deleted'])} deleted ECs")

    query = f"(reviewed:true) AND (date_created:[{args.snapshot} TO *])"
    new_path = config.CX_DOWNLOADS / (f"uniprot_new_since_{args.snapshot}" + ("_smoke" if args.smoke else "") + ".tsv.gz")
    info_up = dl.uniprot_table(query, dl.NEW_ENTRY_FIELDS, new_path, max_records=400 if args.smoke else None,
                               force=args.force)
    if args.with_pfam:
        dl.uniprot_table("(reviewed:true) AND (ec:*)", dl.ALL_EC_FIELDS,
                         config.CX_DOWNLOADS / "uniprot_sp_ec_all.tsv.gz", force=args.force)
        dl.plain_file(dl.PFAM_CLANS_URL, config.CX_DOWNLOADS / "Pfam-A.clans.tsv.gz", force=args.force)
        dl.plain_file(dl.PFAM_RELNOTES_URL, config.CX_DOWNLOADS / "Pfam_relnotes.txt", force=args.force)

    df = dl.read_uniprot_tsv(new_path, dl.NEW_ENTRY_FIELDS)
    n_raw = len(df)
    train_acc, train_seq = set(), set()
    for split in ("split10", "split30", "split50", "split70", "split100"):
        p = config.DATA_DIR / f"{split}.csv"
        if p.is_file():
            t = pd.read_csv(p, sep="\t", dtype=str, keep_default_na=False)
            train_acc |= set(t.iloc[:, 0].str.strip())
            train_seq |= set(t.iloc[:, 2].str.strip().str.upper())
    legacy = set()
    for p in list(config.LEGACY_SETS.values()) + [config.DATA_DIR / "datasets" / "halogenase.csv"]:
        if p.is_file():
            legacy |= set(pd.read_csv(p, sep="\t", dtype=str).Entry)

    reasons = {}

    def drop(mask, why):
        nonlocal df
        reasons[why] = int(mask.sum())
        df = df[~mask]

    df["seq"] = df.sequence.str.upper()
    df["len"] = df.seq.str.len()
    drop(df.fragment.str.strip() != "", "fragment")
    drop((df.len < args.min_len) | (df.len > args.max_len), f"length_outside_{args.min_len}-{args.max_len}")
    drop(df.accession.isin(train_acc), "accession_in_CLEAN_training_csv")
    drop(df.accession.isin(legacy), "accession_in_legacy_benchmarks")
    drop(df.seq.isin(train_seq), "sequence_identical_to_training")
    df = df.sort_values("accession")
    drop(df.duplicated("seq", keep="first"), "duplicate_sequence_within_new_entries")

    df["all_ecs"] = df.ec.map(L.split_ecs)
    df["complete"] = df.all_ecs.map(L.complete_ecs)
    kw = df.keywordid.map(lambda s: {k.strip() for k in s.split(";") if k.strip()})
    is_enz = df.complete.map(len) > 0
    partial_only = (df.all_ecs.map(len) > 0) & ~is_enz
    no_ec = df.all_ecs.map(len) == 0
    enz_kw = kw.map(lambda s: bool(s & ENZYME_CLASS_KW))
    reasons["partial_EC_only(excluded)"] = int(partial_only.sum())
    reasons["no_EC_but_enzyme_class_keyword(excluded)"] = int((no_ec & enz_kw).sum())
    enz = df[is_enz].copy()
    non = df[no_ec & ~enz_kw].copy()

    train100 = pd.read_csv(config.BACKBONES["split100"]["csv"], sep="\t", dtype=str, keep_default_na=False)
    counts = {}
    for ecs in train100["EC number"]:
        for e in L.split_ecs(ecs):
            counts[e] = counts.get(e, 0) + 1
    cur_counts = forward_counts(counts, fmap)
    enz["tier"] = enz.complete.map(lambda e: L.query_tier(e, cur_counts))
    enz["novel_EC"] = enz.complete.map(lambda e: any(cur_counts.get(x, 0) == 0 for x in e))

    rng = np.random.default_rng(args.seed)
    enz["split"] = "test"
    for t, g in enz.groupby("tier"):
        perm = rng.permutation(g.index.to_numpy())
        enz.loc[perm[:int(round(args.cal_frac * len(perm)))], "split"] = "cal"
    if len(non) > args.nonenz_max:
        non = non.loc[rng.choice(non.index.to_numpy(), size=args.nonenz_max, replace=False)]
    non = non.sort_values("accession")

    names = {"cal": f"{prefix}_temporal_cal", "test": f"{prefix}_temporal_test"}
    for sp, name in names.items():
        g = enz[enz.split == sp]
        ids = list(g.accession)
        L.write_dataset(config.CX_SETS / f"{name}.csv", ids, dict(zip(ids, g.complete)), dict(zip(ids, g.seq)))
        L.write_fasta(config.CX_FASTA / f"{name}.fasta", ids, dict(zip(ids, g.seq)))
    nid = list(non.accession)
    L.write_dataset(config.CX_SETS / f"{prefix}_nonenzyme.csv", nid, {}, dict(zip(nid, non.seq)))
    L.write_fasta(config.CX_FASTA / f"{prefix}_nonenzyme.fasta", nid, dict(zip(nid, non.seq)))
    meta = pd.concat([
        enz.assign(set=enz.split.map(names))[["accession", "set", "tier", "novel_EC", "date_created",
                                             "organism_name", "protein_existence", "len", "ec", "xref_pfam"]],
        non.assign(set=f"{prefix}_nonenzyme", tier="non_enzyme", novel_EC=False)[
            ["accession", "set", "tier", "novel_EC", "date_created", "organism_name", "protein_existence",
             "len", "ec", "xref_pfam"]]])
    meta.to_csv(config.CX_SETS / f"{prefix}_temporal_meta.csv", index=False)

    by_tier = enz.groupby(["split", "tier"]).size().unstack(fill_value=0)
    manifest = {"snapshot_created_on_or_after": args.snapshot, "uniprot": info_up, "enzyme_dat": info_enz,
                "enzyme_release": nom["release"], "raw_records": n_raw, "excluded": reasons,
                "n_enzymes": len(enz), "n_cal": int((enz.split == "cal").sum()),
                "n_test": int((enz.split == "test").sum()), "n_nonenzyme": len(non),
                "n_novel_EC_proteins": int(enz.novel_EC.sum()),
                "tiers": {sp: {k: int(v) for k, v in by_tier.loc[sp].items()} for sp in by_tier.index},
                "params": vars(args)}
    L.save_json(manifest, config.CX_SETS / f"{prefix}_temporal_manifest.json")
    print(f"raw {n_raw} records | excluded {reasons}")
    print(f"enzymes {len(enz)} (cal {manifest['n_cal']}, test {manifest['n_test']}; "
          f"{manifest['n_novel_EC_proteins']} with an EC unseen by CLEAN) | non-enzymes {len(non)}")
    print(by_tier.to_string())
    print(f"\nCX02 DONE in {time.time() - t0:.0f}s -> {config.CX_SETS}")


if __name__ == "__main__":
    main()
