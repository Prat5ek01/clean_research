"""E3: why the centroid assumption fails - EC classes whose training members split into
non-homologous groups (no shared Pfam clan; candidate non-homologous isofunctional enzymes, NISE).

Per-EC geometry from CLEAN's training embeddings (split100 or split70):
  loo_misassign = share of members closer to another EC's centroid than to their own
                  leave-one-out centroid; rms_radius = within-EC spread.
Per-EC biology from current Swiss-Prot Pfam annotations of the same training accessions:
  clan groups = connected components of members linked by a shared Pfam clan (a family without
                a clan counts as its own clan); nise_like = >=2 groups each covering >=20% of the
                Pfam-annotated members.
E3b (--write-nise-fasta, cx_03_embed.py --sets nise_members, then --subgroups): are those
subgroups separate clusters in CLEAN space, i.e. is the single EC centroid between them?"""
import argparse
import json
import time
import warnings
from collections import defaultdict

warnings.filterwarnings("ignore")

from cleanx import config

config.enter_app_dir()

import numpy as np
import pandas as pd
import torch
from scipy import stats
from sklearn.metrics import roc_auc_score, silhouette_score

from cleanx import labels as L
from cleanx.index import TrainIndex, load_esm_means
from cleanx.pipeline import fmt


def coannotated(idx):
    """EC column -> columns of ECs that share at least one training sequence with it."""
    co = defaultdict(set)
    for ecs in idx.id_ec.values():
        cols = [idx.col[e] for e in ecs if e in idx.col]
        for c in cols:
            co[c].update(x for x in cols if x != c)
    return {c: np.array(sorted(v)) for c, v in co.items() if v}


def loo_geometry(idx, chunk=20000):
    """A member's own EC and every EC co-annotated with it are excluded from 'other' centroids,
    because rows cannot be mapped back to accessions (a multi-EC sequence sits in several blocks)."""
    n_row = idx.sizes[idx.row_ec].astype(float)
    co = coannotated(idx)
    mis = np.zeros(len(idx.E), dtype=bool)
    margin = np.full(len(idx.E), np.nan)
    for s in range(0, len(idx.E), chunk):
        X = idx.E[s:s + chunk]
        ri = idx.row_ec[s:s + chunk]
        nn = n_row[s:s + chunk]
        d2 = np.einsum("ij,ij->i", X, X)[:, None] + idx.C_sq[None, :] - 2.0 * X @ idx.centroids.T
        d = np.sqrt(np.maximum(d2, 0))
        loo = (nn[:, None] * idx.centroids[ri] - X) / np.maximum(nn[:, None] - 1, 1)
        d_own = np.sqrt(((X - loo) ** 2).sum(1))
        d[np.arange(len(X)), ri] = np.inf
        for c in np.unique(ri):
            if c in co:
                rows_c = np.nonzero(ri == c)[0]
                d[np.ix_(rows_c, co[c])] = np.inf
        d_other = d.min(1)
        ok = nn > 1
        mis[s:s + chunk] = ok & (d_other < d_own)
        margin[s:s + chunk] = np.where(ok, d_other - d_own, np.nan)
    rate = np.bincount(idx.row_ec, weights=mis, minlength=len(idx.ecs)) / idx.sizes
    mean_margin = np.array([np.nanmean(margin[a:a + n]) if n > 1 else np.nan
                            for a, n in zip(idx.starts, idx.sizes)])
    rms = np.sqrt(np.bincount(idx.row_ec, weights=((idx.E - idx.centroids[idx.row_ec]) ** 2).sum(1),
                              minlength=len(idx.ecs)) / idx.sizes)
    rate[idx.sizes < 2] = np.nan
    rms[idx.sizes < 2] = np.nan
    return rate, mean_margin, rms


def clan_groups(accs, pf, clan_of):
    """acc -> group rank (0 = largest group) for members with Pfam annotation; plus group sizes."""
    clans = {}
    for a in accs:
        fam = [x for x in pf.get(a, "").split(";") if x]
        if fam:
            clans[a] = sorted({clan_of.get(x, x) for x in fam})
    parent = {}

    def find(x):
        while parent.setdefault(x, x) != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for cl in clans.values():
        for c in cl[1:]:
            ra, rb = find(cl[0]), find(c)
            if ra != rb:
                parent[max(ra, rb)] = min(ra, rb)
    root = {a: find(cl[0]) for a, cl in clans.items()}
    size = defaultdict(int)
    for r in root.values():
        size[r] += 1
    order = sorted(size, key=lambda r: (-size[r], r))
    rank = {r: k for k, r in enumerate(order)}
    return {a: rank[r] for a, r in root.items()}, [size[r] for r in order]


def subgroup_analysis(idx, ec_members, pf, clan_of, geo, out, max_n, perm, seed):
    rng = np.random.default_rng(seed)
    rows, skipped = [], 0
    for ec in geo[(geo.nise_like) & (geo.n <= max_n)].ec:
        grp, sizes = clan_groups(ec_members[ec], pf, clan_of)
        n_pf = len(grp)
        keep = {g for g, s in enumerate(sizes) if s >= max(2, 0.2 * n_pf)}
        accs = sorted(a for a, g in grp.items() if g in keep)
        if len(keep) < 2 or any(not (config.ESM_DIR / f"{a}.pt").is_file() for a in accs):
            skipped += 1
            continue
        with torch.no_grad():
            Z = idx.model()(torch.from_numpy(load_esm_means(accs))).numpy().astype(float)
        lab = np.array([grp[a] for a in accs])
        sil = silhouette_score(Z, lab)
        null = np.array([silhouette_score(Z, rng.permutation(lab)) for _ in range(perm)])
        d_ec = np.linalg.norm(Z - Z.mean(0), axis=1).mean()
        d_grp = np.mean(np.concatenate([np.linalg.norm(Z[lab == g] - Z[lab == g].mean(0), axis=1)
                                        for g in np.unique(lab)]))
        rows.append({"ec": ec, "n_members": len(accs), "n_groups": len(keep), "silhouette": sil,
                     "perm_p": float((1 + np.sum(null >= sil)) / (1 + perm)),
                     "mean_dist_to_EC_centroid": d_ec, "mean_dist_to_subgroup_centroid": d_grp})
    df = pd.DataFrame(rows)
    df.to_csv(out / "E3b_nise_subgroups.csv", index=False)
    if len(df):
        ratio = (df.mean_dist_to_EC_centroid / df.mean_dist_to_subgroup_centroid).median()
        msg = (f"E3b: {len(df)} NISE-like ECs analysed ({skipped} skipped: missing embeddings or <2 groups); "
               f"median silhouette {df.silhouette.median():.3f}; {(df.perm_p < 0.05).mean():.0%} have separated "
               f"subgroups (permutation p<0.05); members are {ratio:.2f}x farther from the EC centroid than "
               f"from their own subgroup centroid (median)")
    else:
        msg = "E3b: no NISE-like EC has all member embeddings yet (run --write-nise-fasta, then cx_03_embed.py)"
    print(msg)
    return msg


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--backbone", default="split100", choices=list(config.BACKBONES))
    ap.add_argument("--min-n", type=int, default=5)
    ap.add_argument("--min-pfam-frac", type=float, default=0.8)
    ap.add_argument("--perm", type=int, default=10000)
    ap.add_argument("--seed", type=int, default=config.SEED)
    ap.add_argument("--smoke", action="store_true", help="1000 permutations")
    ap.add_argument("--write-nise-fasta", action="store_true",
                    help="write FASTA of members of NISE-like ECs (n <= --nise-max-n) for cx_03_embed.py")
    ap.add_argument("--subgroups", action="store_true", help="E3b subgroup analysis (needs those embeddings)")
    ap.add_argument("--nise-max-n", type=int, default=150)
    ap.add_argument("--subgroup-perm", type=int, default=200)
    args = ap.parse_args()
    if args.smoke:
        args.perm = 1000
    out = config.CX_RESULTS / (f"E3_{args.backbone}" + ("_smoke" if args.smoke else ""))
    out.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    up_path = config.CX_DOWNLOADS / "uniprot_sp_ec_all.tsv.gz"
    clan_path = config.CX_DOWNLOADS / "Pfam-A.clans.tsv.gz"
    for p in (up_path, clan_path):
        if not p.is_file():
            raise SystemExit(f"{p} missing: run  cx.bat cx_02_build_temporal.py --with-pfam")

    idx = TrainIndex(args.backbone)
    rate, margin, rms = loo_geometry(idx)
    print(f"LOO centroid misassignment (co-annotated ECs masked): {np.nanmean(rate):.3f} mean per EC; "
          f"{time.time() - t0:.0f}s")

    up = pd.read_csv(up_path, sep="\t", dtype=str, keep_default_na=False, compression="gzip")
    up.columns = ["accession", "ec", "date_created", "xref_pfam", "length"]
    pf = dict(zip(up.accession, up.xref_pfam))
    cl = pd.read_csv(clan_path, sep="\t", header=None, dtype=str, keep_default_na=False, compression="gzip")
    clan_of = {a: (c if c else a) for a, c in zip(cl[0], cl[1])}

    ec_members = defaultdict(list)
    for acc, ecs in idx.id_ec.items():
        for e in ecs:
            ec_members[e].append(acc)
    rows = []
    for j, e in enumerate(idx.ecs):
        accs = ec_members[e]
        grp, sizes = clan_groups(accs, pf, clan_of)
        n_pf = len(grp)
        archs = pd.Series([tuple(sorted(x for x in pf.get(a, "").split(";") if x)) for a in grp])
        ac = archs.value_counts() if n_pf else pd.Series(dtype=int)
        big = [s for s in sizes if s >= max(2, 0.2 * n_pf)]
        rows.append({"ec": e, "ec_class": e.split(".")[0], "n": int(idx.sizes[j]), "tier": idx.tiers[j],
                     "n_with_pfam": n_pf, "loo_misassign": rate[j], "mean_margin": margin[j],
                     "rms_radius": rms[j], "n_arch": int(len(ac)),
                     "top_arch_share": float(ac.iloc[0] / n_pf) if n_pf else np.nan,
                     "n_clan_groups": len(sizes), "n_big_clan_groups": len(big), "nise_like": len(big) >= 2})
    df = pd.DataFrame(rows)
    mp_file = config.CX_RESULTS / f"E2_{args.backbone}" / "multiproto_ecs.csv"
    if mp_file.is_file():
        mp = pd.read_csv(mp_file)
        df = df.merge(mp[["ec", "k", "silhouette"]].rename(columns={"k": "multiproto_k",
                                                                     "silhouette": "multiproto_silhouette"}),
                      on="ec", how="left")
    df.to_csv(out / "E3_ec_geometry.csv", index=False)

    sel = df[(df.n >= args.min_n) & (df.n_with_pfam >= args.min_pfam_frac * df.n)].copy()
    a, b = sel[sel.nise_like], sel[~sel.nise_like]
    res = {"backbone": args.backbone, "n_EC_analysed": int(len(sel)), "n_nise_like": int(len(a)),
           "share_nise_like": float(len(a) / max(len(sel), 1))}
    cols = ["loo_misassign", "rms_radius", "mean_margin"] + (["multiproto_silhouette"] if "multiproto_silhouette" in sel else [])
    for col in cols:
        m = sel[col].notna()
        if m.sum() < 10 or sel.nise_like[m].nunique() < 2:
            continue
        u = stats.mannwhitneyu(a[col].dropna(), b[col].dropna(), alternative="two-sided")
        res[col] = {"median_nise": float(a[col].median()), "median_other": float(b[col].median()),
                    "mean_nise": float(a[col].mean()), "mean_other": float(b[col].mean()),
                    "mannwhitney_p": float(u.pvalue), "auroc": float(roc_auc_score(sel.nise_like[m], sel[col][m]))}
    rng = np.random.default_rng(args.seed)
    sel["size_bin"] = pd.qcut(np.log(sel.n), q=10, duplicates="drop", labels=False)
    lab = sel.nise_like.to_numpy()
    groups = [np.nonzero(sel.size_bin.to_numpy() == g)[0] for g in np.unique(sel.size_bin)]
    perms = []
    for _ in range(args.perm):
        p = lab.copy()
        for g in groups:
            p[g] = rng.permutation(p[g])
        perms.append(p)
    res["size_stratified_permutation"] = {}
    for col in ("loo_misassign", "rms_radius", "mean_margin"):
        y = sel[col].to_numpy()
        ok = ~np.isnan(y)
        obs = y[lab & ok].mean() - y[~lab & ok].mean()
        null = np.array([y[p & ok].mean() - y[~p & ok].mean() for p in perms])
        res["size_stratified_permutation"][col] = {
            "observed_mean_diff_nise_minus_other": float(obs), "n_perm": args.perm,
            "p_two_sided": float((1 + np.sum(np.abs(null) >= abs(obs))) / (1 + args.perm))}
    rho = stats.spearmanr(sel.top_arch_share, sel.loo_misassign, nan_policy="omit")
    res["spearman_top_arch_share_vs_loo"] = {"rho": float(rho.correlation), "p": float(rho.pvalue)}
    res["nise_like_by_EC_class"] = {k: int(v) for k, v in a.ec_class.value_counts().sort_index().items()}

    e2 = config.CX_RESULTS / f"E2_{args.backbone}" / "E2_per_query.csv"
    if e2.is_file():
        pq = pd.read_csv(e2, dtype=str, keep_default_na=False)
        geo = df.set_index("ec")
        link = []
        if {"f1_exemplar", "f1_centroid"} <= set(pq.columns):
            for _, r in pq[pq.tier != "unseen"].iterrows():
                tes = [t for t in r.true_EC.split(";") if t in geo.index]
                if tes:
                    g = geo.loc[tes]
                    link.append({"dataset": r.dataset, "nise_like": bool(g.nise_like.any()),
                                 "gain": float(r.f1_exemplar) - float(r.f1_centroid)})
        if link:
            lk = pd.DataFrame(link)
            out_link = {}
            for ds, g in lk.groupby("dataset"):
                x, z = g[g.nise_like].gain, g[~g.nise_like].gain
                out_link[ds] = {"n_nise": int(len(x)), "gain_nise": float(x.mean()) if len(x) else None,
                                "n_other": int(len(z)), "gain_other": float(z.mean()) if len(z) else None,
                                "mannwhitney_p": float(stats.mannwhitneyu(x, z).pvalue) if len(x) > 4 and len(z) > 4 else None}
            res["exemplar_gain_by_true_EC_nise_status"] = out_link

    if args.write_nise_fasta:
        ecs_sel = df[(df.nise_like) & (df.n <= args.nise_max_n)].ec
        accs = sorted({x for e in ecs_sel for x in ec_members[e]})
        seqs = idx.train_sequences()
        L.write_fasta(config.CX_FASTA / "nise_members.fasta", accs, seqs)
        print(f"wrote {len(accs)} member sequences of {len(ecs_sel)} NISE-like ECs -> "
              f"{config.CX_FASTA / 'nise_members.fasta'}\n  next: cx.bat cx_03_embed.py --sets nise_members")
    if args.subgroups:
        res["E3b"] = subgroup_analysis(idx, ec_members, pf, clan_of, df, out, args.nise_max_n,
                                       args.subgroup_perm, args.seed)

    with open(out / "E3_stats.json", "w", encoding="utf-8") as fh:
        json.dump(res, fh, indent=2, default=str)
    text = (f"# E3 summary ({args.backbone})\n\n"
            f"- ECs analysed (n>={args.min_n}, >={args.min_pfam_frac:.0%} of members with Pfam): {res['n_EC_analysed']}; "
            f"NISE-like (>=2 non-homologous Pfam-clan groups): {res['n_nise_like']} ({res['share_nise_like']:.1%}); "
            f"by EC class {res['nise_like_by_EC_class']}\n")
    for col, label in (("rms_radius", "within-EC spread (RMS radius)"),
                       ("mean_margin", "margin to nearest other centroid"),
                       ("loo_misassign", "LOO centroid misassignment")):
        s, pm = res.get(col), res["size_stratified_permutation"][col]
        if s:
            text += (f"- {label}, NISE-like vs other: median {fmt(s['median_nise'])} vs {fmt(s['median_other'])}, "
                     f"Mann-Whitney p={s['mannwhitney_p']:.3g}, AUROC {fmt(s['auroc'])}; size-stratified permutation "
                     f"(mean diff {pm['observed_mean_diff_nise_minus_other']:+.3f}) p={pm['p_two_sided']:.4g}\n")
    text += f"- Spearman(top Pfam-architecture share, misassignment) rho={res['spearman_top_arch_share_vs_loo']['rho']:.3f}\n"
    for ds, v in res.get("exemplar_gain_by_true_EC_nise_status", {}).items():
        text += (f"- {ds}: exemplar-minus-centroid per-protein F1 gain, true EC NISE-like "
                 f"{v['gain_nise']} (n={v['n_nise']}) vs other {v['gain_other']} (n={v['n_other']}), "
                 f"Mann-Whitney p={v['mannwhitney_p']}\n")
    if "E3b" in res:
        text += f"- {res['E3b']}\n"
    with open(out / "E3_summary.md", "w", encoding="utf-8") as fh:
        fh.write(text)
    print(text)
    print(f"E3 DONE in {time.time() - t0:.0f}s -> {out}")


if __name__ == "__main__":
    main()
