"""E7: label provenance. Is BLASTp's advantage over CLEAN on SP-temporal a property of the methods, or partly of how
the benchmark's EC labels were assigned?

PRE-REGISTRATION (drafted 2026-09-29, finalised 2026-10-07, before any UniProt evidence had been downloaded for
any benchmark protein).
The analysis code (this file, cleanx/provenance.py, cleanx/contrast.py, cleanx/metrics.py) is frozen by committing
it; the confirmatory run refuses to start from uncommitted code (cleanx/freeze.py). Known when this was written:
the overall and per-identity-bin TEST results of C3 (BLASTp top hit 0.735 vs CLEAN 0.641 micro-F1 on the 876 closed
TEST proteins; confirmatory run of 2026-09-28). Not known: the provenance of any label, the class sizes, and every
class-specific result.

Background. Swiss-Prot curators assign many EC numbers by transferring them from a characterised homolog ('by
similarity', ECO:0000250, which cites the source protein) or from a family rule such as HAMAP (ECO:0000255). A
homology-transfer predictor re-enacts that process, so on such labels part of its agreement with the benchmark is
built in (circularity). Only labels with experimental evidence (ECO:0000269) test prediction against a reference that
was not itself derived from sequence similarity.

Data. Current UniProtKB entries of the benchmark proteins (cx_12_uniprot_evidence.py; one release, recorded in
data/cleanx/provenance/manifest.json). Predictions and truth from the confirmatory E2 run
(results/cleanx/E2_split100/E2_per_query.csv: pred_homology_tophit = BLASTp top hit against split100, pred_centroid
= CLEAN; ENZYME transfers applied as in E2). BLASTp tables in data/cleanx/identity/. Price-149 is not included: its
proteins are not UniProtKB entries and its labels come from mutant fitness experiments.
Classes (cleanx/provenance.py). EC level: EXP (ECO:0000269) > SIM (ECO:0000250) > RULE (ECO:0000255, ECO:0000256,
ECO:0000259) > OTHER. Protein level: LABEL_ABSENT if the entry is not found or any benchmark label is not an EC on
the current entry; else SIM if any label is SIM; else EXP if every label is EXP; else RULE if any label is RULE;
else OTHER. PROPAGATED = SIM or RULE (labels transferred by homology).
Population. SP-temporal TEST proteins whose true ECs are all in CLEAN's vocabulary (closed; n = 876).

Primary endpoints (one family of two, Holm). EC depth 4, micro-F1, paired bootstrap over proteins, B = 2000.
  P7a  delta_EXP = micro-F1(BLASTp) - micro-F1(CLEAN) on the closed EXP proteins; two-sided; seed 20260927.
  P7b  DiD = delta_PROPAGATED - delta_EXP; EXP and PROPAGATED proteins resampled independently (seeds 20260927 and
       20260928); two-sided.
Gate G7. The primaries are confirmatory only if n_EXP >= 100 and n_PROPAGATED >= 50; otherwise every E7 result is
descriptive. Effects are observed full-sample differences with 95% percentile intervals; MDE80 = 2.80 x bootstrap SE.

Secondary (CIs reported; no claim beyond the primaries).
  S7a  DiD for SIM and for RULE separately (Holm over the two).
  S7b  identity-standardised P7a and P7b: per-bin deltas over the five BLASTp identity bins (maximum identity to
       split100), weighted by the pooled EXP + PROPAGATED bin distribution; bins with fewer than 10 proteins in either
       class are merged by the fixed rule in cleanx/contrast.merge_bins; every (class, bin) cell resampled
       independently.
  S7c  P7a and P7b with mean per-protein F1 instead of micro-F1.
  S7d  EC-Bench test set (100% regime; proteins whose EC-Bench labels are all complete ECs; classes from the same
       UniProt download against EC-Bench's labels): DiD for blastp - clean-supcon, micro-F1 at depth 4 (seeds as in
       P7b); descriptive unless n_EXP >= 30 and n_PROPAGATED >= 30.
Mechanism (closed TEST SIM proteins). Share whose BLASTp top hit is a protein cited as the source of one of their SIM
labels (M7a), whose cited source is among the kept BLASTp hits (M7b), or is a split100 training protein (M7c), with
exact 95% CIs; BLASTp and CLEAN micro-F1 when the top hit is a cited source and when it is not (M7d).
Sensitivity. P7a and P7b with (a) proteins resampled in clusters of their rarest true EC (as in E2), (b) at most 5
proteins per rarest true EC (seeded sample, as in E2).
Descriptive. Class composition of New-392, SP-temporal CAL and TEST and EC-Bench test; ECO code frequencies; deltas
by class x identity bin and class x superkingdom. LABEL_ABSENT and OTHER are reported, never pooled with a class.
Interpretation (fixed in advance; q = Holm-adjusted p, alpha = 0.05):
  G7 fails                               -> descriptive only.
  P7a q<.05, delta>0 and P7b q>=.05      -> BLASTp's advantage holds on experimentally supported labels; no
                                            detectable provenance effect (report the DiD CI and MDE80).
  P7a q<.05, delta>0 and P7b q<.05, DiD>0 -> the advantage holds on experimental labels and is larger on propagated
                                            labels; the provenance-matched estimate is delta_EXP.
  P7a q>=.05 and P7b q<.05, DiD>0         -> the advantage is concentrated in propagated labels and is not
                                            significant on experimental labels; C3's BLASTp advantage is qualified.
  P7a q>=.05 and P7b q>=.05              -> no significant difference on experimental labels and no detectable
                                            provenance effect.
  a significant negative effect          -> reported as such (CLEAN better on experimental labels, or a smaller
                                            BLASTp advantage on propagated labels).
Pilot. New-392 (development set) through the same code (--pilot). It may change parsing and class rules, as a dated
amendment recorded here before the confirmatory run; it may not change the endpoints.

Modes: default (confirmatory; committed code required), --pilot (New-392; writes E7_pilot), --smoke (synthetic data
with known answers; writes E7_smoke), --allow-dirty (development and tests only: outputs are labelled NOT
CONFIRMATORY)."""
import argparse
import sys
import tempfile
import time
import warnings
from pathlib import Path

warnings.filterwarnings("ignore")

from cleanx import config

config.enter_app_dir()

import numpy as np
import pandas as pd

from cleanx import blast as BL
from cleanx import contrast as C
from cleanx import ecbench as ECB
from cleanx import freeze
from cleanx import labels as L
from cleanx import provenance as PV
from cleanx.calibration import risk_ci
from cleanx.identity import IDENTITY_BINS
from cleanx.metrics import evaluate, holm
from cleanx.nomenclature import forward_counts
from cleanx.pipeline import fmt, load_fmap, write_rows

PROV_DIR = config.CX_DATA / "provenance"
E2_PER_QUERY = config.CX_RESULTS / "E2_split100" / "E2_per_query.csv"
BLAST_M, CLEAN_M = "homology_tophit", "centroid"
EB_BLAST, EB_CLEAN = "blastp", "clean-supcon"
GATE = {"EXP": 100, "PROPAGATED": 50}
GATE_ECB = 30
MIN_BIN = 10
CAP = 5
FROZEN = ["cx_13_provenance.py", "cleanx/provenance.py", "cleanx/contrast.py", "cleanx/metrics.py",
          "cleanx/freeze.py"]
fails = []


def check(name, ok, detail=""):
    print(f"[{'PASS' if ok else 'FAIL'}] {name}" + (f" -- {detail}" if detail else ""), flush=True)
    if not ok:
        fails.append(name)
    return ok


# ----------------------------------------------------------------------------------------------- inputs
def read_provenance(tag):
    """(entries DataFrame, {query: {ec: [codes]}}, {query: {ec: [source accessions]}}) from cx_12's tables."""
    e_path, c_path = PROV_DIR / f"{tag}_entries.csv", PROV_DIR / f"{tag}_ec_evidence.csv"
    if not e_path.is_file() or not c_path.is_file():
        raise SystemExit(f"{e_path} missing: run cx_12_uniprot_evidence.py on a machine with internet access")
    ent = pd.read_csv(e_path, dtype=str, keep_default_na=False)
    ecs = pd.read_csv(c_path, dtype=str, keep_default_na=False)
    return tables_to_maps(ent, ecs)


def tables_to_maps(ent, ecs):
    ent = ent.copy()
    ent["found"] = ent["found"].astype(str).str.lower() == "true"
    codes, sources = {}, {}
    for r in ecs.itertuples(index=False):
        codes.setdefault(r.query_accession, {})[r.ec] = [c for c in str(r.codes).split(";") if c]
        sources.setdefault(r.query_accession, {})[r.ec] = [PV.source_accession(s) for s in
                                                            str(r.sim_sources).split(";") if s]
    return ent, codes, sources


def classify(ids, truth, ent, codes, sources):
    """Per protein: class, label classes, cited sources of its SIM labels, superkingdom, resolution."""
    meta = ent.set_index("query_accession")
    out = {}
    for i in ids:
        found = bool(meta.found.get(i, False)) if i in meta.index else False
        c = codes.get(i, {})
        cls = PV.protein_class(truth[i], c, found)
        lab = [PV.ec_class(c[e]) for e in truth[i]] if cls != "LABEL_ABSENT" else []
        src = sorted({s for e, k in zip(truth[i], lab) if k == "SIM" for s in sources.get(i, {}).get(e, [])})
        out[i] = {"class": cls, "label_classes": ";".join(lab), "sources": ";".join(src),
                  "superkingdom": str(meta.superkingdom.get(i, "")) if i in meta.index and found else "",
                  "resolution": str(meta.resolution.get(i, "")) if i in meta.index else "not_queried"}
    return out


def load_e2(dataset):
    if not E2_PER_QUERY.is_file():
        raise SystemExit(f"{E2_PER_QUERY} missing (confirmatory E2 run)")
    df = pd.read_csv(E2_PER_QUERY, dtype=str, keep_default_na=False)
    df = df[df.dataset == dataset]
    need = {f"pred_{BLAST_M}", f"pred_{CLEAN_M}", "true_EC", "tier", "identity_bin"}
    if df.empty or need - set(df.columns):
        raise SystemExit(f"{E2_PER_QUERY}: no {dataset} rows with {sorted(need)}")
    ids = list(df.Entry)
    truth = {i: L.split_ecs(t) for i, t in zip(ids, df.true_EC)}
    preds = {m: {i: L.split_ecs(p) for i, p in zip(ids, df[f"pred_{m}"])} for m in (BLAST_M, CLEAN_M)}
    return ids, truth, preds, dict(zip(ids, df.tier)), dict(zip(ids, df.identity_bin))


def rarest_ec(ids, truth, temporal):
    """id -> rarest true EC in split100 (training counts after ENZYME transfers for temporal sets), as in E2."""
    from CLEAN.utils import get_ec_id_dict
    _, ec_id = get_ec_id_dict(str(config.BACKBONES["split100"]["csv"]))
    counts = {e: len(v) for e, v in ec_id.items()}
    if temporal:
        counts = forward_counts(counts, load_fmap(required=True))
    return {i: min(truth[i], key=lambda e: (counts.get(e, 0), e)) if truth[i] else "none" for i in ids}


def blast_context(dataset):
    """(top hit per query, set of kept hits per query, split100 accessions)."""
    ident = pd.read_csv(config.CX_IDENTITY / f"{dataset}.csv", dtype=str, keep_default_na=False)
    if "tool" not in ident.columns or (ident.tool != "blast").any():
        raise SystemExit(f"{dataset}: identity table is not from BLASTp (run cx_04_identity.py --tool blast)")
    top = dict(zip(ident.Entry, ident.top_hit))
    hits = BL.read_hits(config.CX_IDENTITY / "blast" / f"{dataset}.tsv.gz")
    kept = hits.groupby("qseqid").sseqid.apply(set).to_dict()
    train = set(pd.read_csv(config.BACKBONES["split100"]["csv"], sep="\t", dtype=str, usecols=["Entry"]).Entry)
    return top, kept, train


# ----------------------------------------------------------------------------------------------- analysis
def capped(ids, cluster, k, seed):
    """At most k proteins per cluster, chosen with a fixed seed; order of ids kept (same rule as E2)."""
    rng = np.random.default_rng(seed)
    groups = {}
    for i in ids:
        groups.setdefault(cluster[i], []).append(i)
    keep = set()
    for members in groups.values():
        keep.update(members if len(members) <= k else list(rng.choice(members, size=k, replace=False)))
    return [i for i in ids if i in keep]


def did_rows(names, family, population, truth, pa, pb, exp_ids, other_ids, key, B, seed, clusters=None, gate=True):
    """Two rows: delta_EXP (named names[0]) and DiD = delta(other) - delta_EXP (named names[1])."""
    r = C.did(truth, pa, pb, exp_ids, other_ids, key, B, seed, clusters=clusters)
    base = {"family": family, "population": population, "n0": len(exp_ids), "n1": len(other_ids), "metric": key,
            "gate_passed": gate, "B": B, "seed": seed}
    return [{"endpoint": names[0], "contrast": "delta_EXP", **base, **r["delta0"]},
            {"endpoint": names[1], "contrast": "DiD", **base, **r["did"]}]


def verdict(p7a, p7b, gate):
    """The pre-registered interpretation for P7a (delta_EXP) and P7b (DiD), from Holm-adjusted p-values."""
    if not gate:
        return "G7 failed: descriptive only"
    sa, sb = p7a["p_holm"] < 0.05, p7b["p_holm"] < 0.05
    if (sa and p7a["estimate"] < 0) or (sb and p7b["estimate"] < 0):
        parts = []
        if sa and p7a["estimate"] < 0:
            parts.append("CLEAN is better than BLASTp on experimentally supported labels")
        if sb and p7b["estimate"] < 0:
            parts.append("BLASTp's advantage is smaller on propagated labels")
        return "significant negative effect: " + "; ".join(parts)
    if sa and not sb:
        return ("BLASTp's advantage holds on experimentally supported labels; no detectable provenance effect "
                f"(DiD {p7b['estimate']:+.3f} [{p7b['ci_lo']:+.3f}, {p7b['ci_hi']:+.3f}], MDE80 {p7b['mde80']:.3f})")
    if sa and sb:
        return ("BLASTp's advantage holds on experimental labels and is larger on propagated labels; the "
                f"provenance-matched estimate is delta_EXP = {p7a['estimate']:+.3f}")
    if sb:
        return ("BLASTp's advantage is concentrated in propagated labels and is not significant on experimental "
                "labels; C3's BLASTp advantage is qualified")
    return "no significant difference on experimental labels and no detectable provenance effect"


def mechanism(sim_ids, cls, top, kept, train, truth, pa, pb, B, seed):
    rows, withsrc = [], [i for i in sim_ids if cls[i]["sources"]]
    src = {i: set(cls[i]["sources"].split(";")) for i in withsrc}
    for mid, what, hit in [("M7a", "top hit is a cited source", lambda i: top.get(i, "") in src[i]),
                           ("M7b", "a cited source is among the kept hits", lambda i: bool(src[i] & kept.get(i, set()))),
                           ("M7c", "a cited source is a split100 training protein", lambda i: bool(src[i] & train))]:
        k = sum(1 for i in withsrc if hit(i))
        lo, hi = risk_ci(k, len(withsrc)) if withsrc else (np.nan, np.nan)
        rows.append({"measure": mid, "what": what, "k": k, "n": len(withsrc), "n_sim": len(sim_ids),
                     "share": k / len(withsrc) if withsrc else np.nan, "ci_lo": lo, "ci_hi": hi})
    is_src = [i for i in withsrc if top.get(i, "") in src[i]]
    not_src = [i for i in sim_ids if i not in set(is_src)]
    for label, sub in (("top hit is a cited source", is_src), ("top hit is not a cited source", not_src)):
        row = {"measure": "M7d", "what": label, "n": len(sub)}
        if sub:
            ea, eb = evaluate(truth, pa, sub, 4)["micro_F1"], evaluate(truth, pb, sub, 4)["micro_F1"]
            row.update({"blastp_micro_F1": ea, "clean_micro_F1": eb, "delta": ea - eb})
            if len(sub) >= 10:
                s = C.summarize(*C.delta(truth, pa, pb, sub, "micro_F1", B, seed))
                row.update({"ci_lo": s["ci_lo"], "ci_hi": s["ci_hi"]})
        rows.append(row)
    return rows


def strata_rows(dataset, ids, cls, groups, truth, pa, pb, B, seed, gname):
    rows = []
    for c in PV.PROTEIN_CLASSES:
        for g in sorted({groups[i] for i in ids}):
            sub = [i for i in ids if cls[i]["class"] == c and groups[i] == g]
            if not sub:
                continue
            row = {"dataset": dataset, "class": c, "grouping": gname, "group": g, "n": len(sub),
                   "blastp_micro_F1": evaluate(truth, pa, sub, 4)["micro_F1"],
                   "clean_micro_F1": evaluate(truth, pb, sub, 4)["micro_F1"]}
            row["delta"] = row["blastp_micro_F1"] - row["clean_micro_F1"]
            if len(sub) >= 10:
                o, b = C.delta(truth, pa, pb, sub, "micro_F1", B, seed)
                s = C.summarize(o, b)
                row.update({"ci_lo": s["ci_lo"], "ci_hi": s["ci_hi"]})
            rows.append(row)
    return rows


def composition(dataset, subset, ids, cls, truth, codes):
    rows = []
    n = len(ids)
    for c in PV.PROTEIN_CLASSES:
        k = sum(1 for i in ids if cls[i]["class"] == c)
        rows.append({"dataset": dataset, "subset": subset, "level": "protein", "class": c, "n": k,
                     "share": k / n if n else np.nan, "total": n})
    labs = [(i, e) for i in ids for e in truth[i] if e in codes.get(i, {})]
    missing = sum(1 for i in ids for e in truth[i] if e not in codes.get(i, {}))
    total = len(labs) + missing
    for c in PV.EC_CLASSES:
        k = sum(1 for i, e in labs if PV.ec_class(codes[i][e]) == c)
        rows.append({"dataset": dataset, "subset": subset, "level": "label", "class": c, "n": k,
                     "share": k / total if total else np.nan, "total": total})
    rows.append({"dataset": dataset, "subset": subset, "level": "label", "class": "LABEL_ABSENT", "n": missing,
                 "share": missing / total if total else np.nan, "total": total})
    return rows


def code_rows(dataset, ids, truth, codes):
    cnt = {}
    for i in ids:
        for e in truth[i]:
            for c in codes.get(i, {}).get(e, []) or ["(none)"]:
                cnt[c] = cnt.get(c, 0) + 1
    return [{"dataset": dataset, "eco_code": c, "n_labels": n} for c, n in sorted(cnt.items(), key=lambda x: -x[1])]


def analyse(dataset, data, B, seed, cap, confirmatory):
    """All E7 tables for one benchmark. data: dict with ids, truth, preds, tier, ibin, cls, top, kept, train,
    rarest, codes (see load_dataset / synthetic)."""
    ids, truth, cls = data["ids"], data["truth"], data["cls"]
    pa, pb = data["preds"][BLAST_M], data["preds"][CLEAN_M]
    closed = [i for i in ids if data["tier"][i] != "unseen"]
    by = {c: [i for i in closed if cls[i]["class"] == c] for c in PV.PROTEIN_CLASSES}
    exp_ids = by["EXP"]
    prop_ids = [i for i in closed if cls[i]["class"] in PV.PROPAGATED]
    gate = len(exp_ids) >= GATE["EXP"] and len(prop_ids) >= GATE["PROPAGATED"]
    pop = f"{dataset} closed"
    ep, sens = [], []
    if exp_ids and prop_ids:
        rows = did_rows(("P7a", "P7b"), "primary", pop, truth, pa, pb, exp_ids, prop_ids, "micro_F1", B, seed,
                        gate=gate)
        ph = holm([rows[0]["p_two_sided"], rows[1]["p_two_sided"]])
        rows[0]["p_holm"], rows[1]["p_holm"] = float(ph[0]), float(ph[1])
        ep += rows
        sub = []
        for k in ("SIM", "RULE"):
            if by[k]:
                sub.append(did_rows(("-", f"S7a_{k}"), "S7a", pop, truth, pa, pb, exp_ids, by[k], "micro_F1", B,
                                    seed, gate=gate)[1])
        for r_, q in zip(sub, holm([x["p_two_sided"] for x in sub])):
            r_["p_holm"] = float(q)
        ep += sub
        std = C.standardized(truth, pa, pb, {"EXP": exp_ids, "PROPAGATED": prop_ids}, data["ibin"], IDENTITY_BINS,
                             "micro_F1", B, seed, min_n=MIN_BIN)
        if std is None:
            ep.append({"endpoint": "S7b", "family": "S7b", "population": pop, "contrast": "DiD_standardised",
                       "note": f"not estimable: fewer than {MIN_BIN} proteins per class after merging all bins"})
        else:
            (oe, be), (op, bp) = std["strata"]["EXP"], std["strata"]["PROPAGATED"]
            note = "bins " + ", ".join(f"{g} (w={w:.3f})" for g, w in zip(std["groups"], std["weights"]))
            ep.append({"endpoint": "S7b_delta_EXP", "family": "S7b", "population": pop,
                       "contrast": "delta_EXP_standardised", "n0": len(exp_ids), "n1": len(prop_ids),
                       "metric": "micro_F1", **C.summarize(oe, be), "note": note, "B": B, "seed": seed})
            ep.append({"endpoint": "S7b", "family": "S7b", "population": pop, "contrast": "DiD_standardised",
                       "n0": len(exp_ids), "n1": len(prop_ids), "metric": "micro_F1",
                       **C.summarize(op - oe, bp - be), "note": note, "B": B, "seed": seed})
        ep += did_rows(("S7c_delta_EXP", "S7c"), "S7c", pop, truth, pa, pb, exp_ids, prop_ids, "samples_F1", B,
                       seed, gate=gate)
        for label, clusters, a_ids, p_ids in [
                ("cluster_bootstrap_rarest_EC", data["rarest"], exp_ids, prop_ids),
                (f"capped{cap}_per_rarest_EC", None, capped(exp_ids, data["rarest"], cap, seed),
                 capped(prop_ids, data["rarest"], cap, seed))]:
            sens += did_rows(("P7a", "P7b"), "sensitivity", f"{pop} | {label}", truth, pa, pb, a_ids, p_ids,
                             "micro_F1", B, seed, clusters=clusters, gate=gate)
        p7a, p7b = ep[0], ep[1]
        v = verdict(p7a, p7b, gate) if confirmatory else "not confirmatory: " + verdict(p7a, p7b, gate)
    else:
        v = f"not estimable: n_EXP = {len(exp_ids)}, n_PROPAGATED = {len(prop_ids)}"
    mech = mechanism(by["SIM"], cls, data["top"], data["kept"], data["train"], truth, pa, pb, B, seed)
    strata = (strata_rows(dataset, closed, cls, data["ibin"], truth, pa, pb, B, seed, "identity_bin") +
              strata_rows(dataset, closed, cls, {i: cls[i]["superkingdom"] or "(unknown)" for i in closed}, truth,
                          pa, pb, B, seed, "superkingdom"))
    comp = composition(dataset, "all", ids, cls, truth, data["codes"])
    if len(closed) < len(ids):
        comp += composition(dataset, "closed", closed, cls, truth, data["codes"])
    counts = {c: len(v_) for c, v_ in by.items()}
    return {"endpoints": ep, "sensitivity": sens, "mechanism": mech, "strata": strata, "composition": comp,
            "codes": code_rows(dataset, ids, truth, data["codes"]), "verdict": v, "gate": gate, "counts": counts,
            "n_closed": len(closed), "n_prop": len(prop_ids)}


def ecbench_secondary(B, seed):
    """S7d: EC-Bench test (100% regime), complete labels, blastp - clean-supcon, DiD PROPAGATED - EXP."""
    path = ECB.DEST / ECB.FILES[("test", "100")]
    if not path.is_file():
        ECB.fetch()
    ECB.verify({("test", "100"): path})
    ids, truth, preds = ECB.read_results(path)
    ids = [i for i in ids if truth[i] and all(L.EC_COMPLETE.match(e) for e in truth[i])]
    ent, codes, sources = read_provenance("ecbench_test")
    cls = classify(ids, truth, ent, codes, sources)
    exp_ids = [i for i in ids if cls[i]["class"] == "EXP"]
    prop_ids = [i for i in ids if cls[i]["class"] in PV.PROPAGATED]
    gate = len(exp_ids) >= GATE_ECB and len(prop_ids) >= GATE_ECB
    comp = composition("ecbench_test", "complete_labels", ids, cls, truth, codes)
    rows = []
    if exp_ids and prop_ids:
        rows = did_rows(("S7d_delta_EXP", "S7d"), "S7d", "EC-Bench test (100%), complete labels", truth,
                        preds[EB_BLAST], preds[EB_CLEAN], exp_ids, prop_ids, "micro_F1", B, seed, gate=gate)
    return rows, comp, code_rows("ecbench_test", ids, truth, codes), gate


# ----------------------------------------------------------------------------------------------- real data
def load_dataset(dataset):
    ids, truth, preds, tier, ibin = load_e2(dataset)
    ent, codes, sources = read_provenance(dataset)
    missing = [i for i in ids if i not in set(ent.query_accession)]
    if missing:
        raise SystemExit(f"{dataset}: {len(missing)} proteins have no provenance row (e.g. {missing[:3]}); "
                         f"re-run cx_12_uniprot_evidence.py")
    cls = classify(ids, truth, ent, codes, sources)
    top, kept, train = blast_context(dataset)
    return {"ids": ids, "truth": truth, "preds": preds, "tier": tier, "ibin": ibin, "cls": cls, "top": top,
            "kept": kept, "train": train, "rarest": rarest_ec(ids, truth, config.is_temporal(dataset)),
            "codes": codes}


def composition_only(tag):
    """Composition and ECO codes of a set whose model results E7 does not use (SP-temporal CAL)."""
    ids, truth, _ = L.load_named_dataset(tag)
    ent, codes, sources = read_provenance(tag)
    cls = classify(ids, truth, ent, codes, sources)
    return composition(tag, "all", ids, cls, truth, codes), code_rows(tag, ids, truth, codes)


def fmt_p(p, B):
    """A bootstrap p-value of 0 means below 1/B."""
    return f"<{1 / B:.2g}" if p == 0 else f"{p:.3g}"


def write_all(out, res, extra_comp, extra_codes, header):
    write_rows(res["endpoints"], out / "E7_endpoints.csv")
    write_rows(res["sensitivity"], out / "E7_sensitivity.csv")
    write_rows(res["mechanism"], out / "E7_mechanism.csv")
    write_rows(res["strata"], out / "E7_strata.csv")
    write_rows(res["composition"] + extra_comp, out / "E7_composition.csv")
    write_rows(res["codes"] + extra_codes, out / "E7_eco_codes.csv")
    lines = [header, "", f"- verdict: **{res['verdict']}**",
             f"- gate G7 (n_EXP >= {GATE['EXP']}, n_PROPAGATED >= {GATE['PROPAGATED']}): "
             f"{'passed' if res['gate'] else 'FAILED'}; closed proteins by class: {res['counts']}"]
    for r in res["endpoints"]:
        if "estimate" in r:
            lines.append(f"- {r['endpoint']} ({r['contrast']}, {r['metric']}, n0={r.get('n0')}, n1={r.get('n1')}): "
                         f"{r['estimate']:+.3f} [{r['ci_lo']:+.3f}, {r['ci_hi']:+.3f}] "
                         f"p={fmt_p(r['p_two_sided'], r['B'])}"
                         + (f" q_holm={fmt_p(r['p_holm'], r['B'])}" if "p_holm" in r else "")
                         + f" MDE80={r['mde80']:.3f}" + (f" ({r['note']})" if r.get("note") else ""))
        else:
            lines.append(f"- {r['endpoint']}: {r.get('note', '')}")
    for r in res["sensitivity"]:
        lines.append(f"- sensitivity {r['population']}: {r['endpoint']} {r['estimate']:+.3f} "
                     f"[{r['ci_lo']:+.3f}, {r['ci_hi']:+.3f}] p={fmt_p(r['p_two_sided'], r['B'])}")
    for r in res["mechanism"]:
        if "share" in r:
            lines.append(f"- {r['measure']} {r['what']}: {r['k']}/{r['n']} = {fmt(r['share'])} "
                         f"[{fmt(r['ci_lo'])}, {fmt(r['ci_hi'])}]")
        else:
            lines.append(f"- {r['measure']} {r['what']}: n={r['n']}, BLASTp {fmt(r.get('blastp_micro_F1', np.nan))},"
                         f" CLEAN {fmt(r.get('clean_micro_F1', np.nan))}")
    (out / "E7_summary.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


# ----------------------------------------------------------------------------------------------- smoke
CELLS = {  # (class, identity bin): (n proteins, BLASTp correct, CLEAN correct); single label, single prediction
    ("EXP", "<30%"): (50, 35, 35), ("EXP", ">=90%"): (50, 45, 35),
    ("SIM", "<30%"): (10, 5, 5), ("SIM", ">=90%"): (50, 49, 31),
    ("RULE", "<30%"): (10, 5, 5), ("RULE", ">=90%"): (30, 25, 25),
    ("OTHER", "30-50%"): (20, 10, 18), ("LABEL_ABSENT", "30-50%"): (10, 10, 0)}
CODE = {"EXP": ("ECO:0000269", "PubMed", "1"), "SIM": ("ECO:0000250", "UniProtKB", "SRC"),
        "RULE": ("ECO:0000255", "HAMAP-Rule", "MF_00001"), "OTHER": ("ECO:0000305", "", "")}


def synthetic():
    """Proteins with known per-cell accuracies, UniProt-shaped entries run through cx_12's parser, BLAST context."""
    ids, truth, pa, pb, tier, ibin, entries, top, kept = [], {}, {}, {}, {}, {}, [], {}, {}
    n, s = 0, 0
    for (c, b), (k, ka, kb) in CELLS.items():
        for j in range(k):
            n += 1
            i = f"Q{n:05d}"
            ec = f"1.1.{n % 7 + 1}.{n}"
            ids.append(i)
            truth[i] = [ec]
            pa[i] = [ec] if j < ka else ["9.9.9.9"]
            pb[i] = [ec] if j < kb else ["8.8.8.8"]
            tier[i], ibin[i] = "common", b
            if c == "SIM":  # the s-th SIM protein: s < 40 source is the top hit; 40 <= s < 50 a lower kept hit
                src = f"S{n:05d}"
                code = (CODE["SIM"][0], CODE["SIM"][1], src)
                top[i] = src if s < 40 else f"T{n:05d}"
                kept[i] = {top[i], src} if s < 50 else {top[i]}
                s += 1
            else:
                code = CODE.get(c, CODE["EXP"])
                top[i] = f"T{n:05d}"
                kept[i] = {top[i]}
            if c == "LABEL_ABSENT":
                entries.append(PV.synthetic_entry(i, {"7.7.7.7": [code]}))
            else:
                entries.append(PV.synthetic_entry(i, {ec: [code]}, lineage=("Eukaryota" if j % 2 else "Bacteria",)))
    for k in range(5):  # unseen-tier proteins: excluded from every closed analysis
        n += 1
        i = f"U{n:05d}"
        ids.append(i)
        truth[i], pa[i], pb[i], tier[i], ibin[i] = ["6.6.6.6"], [], [], "unseen", "<30%"
        entries.append(PV.synthetic_entry(i, {"6.6.6.6": [CODE["EXP"]]}))
        top[i], kept[i] = "", set()
    erows, crows = PV.tables(ids, entries)
    ent, codes, sources = tables_to_maps(pd.DataFrame(erows).astype(str), pd.DataFrame(crows).astype(str))
    cls = classify(ids, truth, ent, codes, sources)
    train = {f"S{m:05d}" for m in range(0, n + 1) if m % 10 != 0}  # 1 in 10 cited sources not in training
    return {"ids": ids, "truth": truth, "preds": {BLAST_M: pa, CLEAN_M: pb}, "tier": tier, "ibin": ibin,
            "cls": cls, "top": top, "kept": kept, "train": train, "rarest": {i: truth[i][0] for i in ids},
            "codes": codes}


def smoke(out, B):
    data = synthetic()
    counts = {c: sum(1 for i in data["ids"] if data["cls"][i]["class"] == c and data["tier"][i] != "unseen")
              for c in PV.PROTEIN_CLASSES}
    check("synthetic proteins are classified as designed",
          counts == {"EXP": 100, "SIM": 60, "RULE": 40, "OTHER": 20, "LABEL_ABSENT": 10}, str(counts))
    res = analyse("synthetic", data, B, config.SEED, CAP, confirmatory=True)
    ep = {(r["endpoint"], r.get("contrast")): r for r in res["endpoints"]}
    p7a, p7b = ep[("P7a", "delta_EXP")], ep[("P7b", "DiD")]
    # hand values: single label + single prediction, so micro-F1 = accuracy
    check("P7a observed delta_EXP = 80/100 - 70/100", abs(p7a["estimate"] - 0.10) < 1e-12, f"{p7a['estimate']:.6f}")
    check("P7b observed DiD = (84-66)/100 - 0.10", abs(p7b["estimate"] - 0.08) < 1e-12, f"{p7b['estimate']:.6f}")
    check("bootstrap intervals contain the observed effects",
          p7a["ci_lo"] <= 0.10 <= p7a["ci_hi"] and p7b["ci_lo"] <= 0.08 <= p7b["ci_hi"],
          f"P7a [{p7a['ci_lo']:.3f}, {p7a['ci_hi']:.3f}], P7b [{p7b['ci_lo']:.3f}, {p7b['ci_hi']:.3f}]")
    q = holm([p7a["p_two_sided"], p7b["p_two_sided"]])
    check("Holm over the two primaries", abs(p7a["p_holm"] - q[0]) < 1e-12 and abs(p7b["p_holm"] - q[1]) < 1e-12 and
          p7a["p_holm"] >= p7a["p_two_sided"])
    check("gate G7 passes at n_EXP = 100, n_PROPAGATED = 100", res["gate"] is True)
    sim, rule = ep[("S7a_SIM", "DiD")], ep[("S7a_RULE", "DiD")]
    check("S7a: SIM DiD = 0.30 - 0.10, RULE DiD = 0.00 - 0.10",
          abs(sim["estimate"] - 0.20) < 1e-12 and abs(rule["estimate"] + 0.10) < 1e-12,
          f"{sim['estimate']:.4f}, {rule['estimate']:.4f}")
    sb, sbe = ep[("S7b", "DiD_standardised")], ep[("S7b_delta_EXP", "delta_EXP_standardised")]
    # merged bins: <30%..70-90% (EXP 50, PROP 20) and >=90% (EXP 50, PROP 80); weights 70/200 and 130/200
    exp_std, prop_std = 0.35 * 0.0 + 0.65 * (45 - 35) / 50, 0.35 * 0.0 + 0.65 * (74 - 56) / 80
    check("S7b: identity-standardised deltas by direct standardisation",
          abs(sbe["estimate"] - exp_std) < 1e-12 and abs(sb["estimate"] - (prop_std - exp_std)) < 1e-12 and
          sb["note"].startswith("bins <30%+30-50%+50-70%+70-90% (w=0.350), >=90% (w=0.650)"),
          f"{sbe['estimate']:.5f}, {sb['estimate']:.5f}; {sb['note']}")
    check("S7c: per-protein F1 equals micro-F1 for single-label single-prediction proteins",
          abs(ep[("S7c", "DiD")]["estimate"] - 0.08) < 1e-12)
    mech = {r["measure"]: r for r in res["mechanism"] if "share" in r}
    check("mechanism shares: top hit = source 40/60, source among kept hits 50/60, source in training 54/60",
          (mech["M7a"]["k"], mech["M7b"]["k"], mech["M7c"]["k"], mech["M7a"]["n"]) == (40, 50, 54, 60),
          str({k: (v["k"], v["n"]) for k, v in mech.items()}))
    lo, hi = mech["M7a"]["ci_lo"], mech["M7a"]["ci_hi"]
    check("exact CI of a share brackets it", lo < 40 / 60 < hi and 0 < lo and hi < 1, f"[{lo:.3f}, {hi:.3f}]")
    sens = {r["population"].split("| ")[1]: r for r in res["sensitivity"] if r["contrast"] == "DiD"}
    check("sensitivity analyses run (cluster bootstrap with singleton clusters, capping with no effect)",
          all(abs(r["estimate"] - 0.08) < 1e-12 for r in sens.values()) and len(sens) == 2)
    comp = {(r["subset"], r["level"], r["class"]): r["n"] for r in res["composition"]}
    check("composition counts closed proteins and labels by class",
          comp[("closed", "protein", "EXP")] == 100 and comp[("closed", "label", "LABEL_ABSENT")] == 10 and
          comp[("all", "protein", "EXP")] == 105)
    check("verdict follows the pre-registered table", res["verdict"] == verdict(p7a, p7b, True))
    fake = {"p_holm": 0.01, "estimate": 0.1, "ci_lo": 0.05, "ci_hi": 0.15, "mde80": 0.05}
    nul = {"p_holm": 0.5, "estimate": 0.01, "ci_lo": -0.05, "ci_hi": 0.07, "mde80": 0.08}
    neg = {"p_holm": 0.01, "estimate": -0.1, "ci_lo": -0.15, "ci_hi": -0.05, "mde80": 0.05}
    check("interpretation table: all six outcomes",
          verdict(fake, nul, True).startswith("BLASTp's advantage holds on experimentally supported labels") and
          verdict(fake, fake, True).startswith("BLASTp's advantage holds on experimental labels and is larger") and
          verdict(nul, fake, True).startswith("BLASTp's advantage is concentrated") and
          verdict(nul, nul, True).startswith("no significant difference") and
          verdict(neg, nul, True).startswith("significant negative effect: CLEAN") and
          verdict(fake, fake, False) == "G7 failed: descriptive only")
    small = dict(data)
    small["cls"] = {i: (dict(v, **{"class": "OTHER"}) if v["class"] == "EXP" and int(i[1:]) > 60 else v)
                    for i, v in data["cls"].items()}
    check("gate G7 fails below 100 EXP proteins", analyse("synthetic", small, 50, config.SEED, CAP, True)["gate"]
          is False)
    check("merge_bins: unmergeable strata give None", C.merge_bins(IDENTITY_BINS, {"a": {"<30%": 3}}, 10) is None)
    with tempfile.TemporaryDirectory() as td:  # the freeze guard labels uncommitted code
        p = Path(td) / "x.py"
        p.write_text("x = 1\n", encoding="utf-8")
        lab = freeze.require_frozen([p], True, "smoke")
        check("freeze guard: a file outside the repository is not confirmatory", lab.startswith("NOT CONFIRMATORY"),
              lab)
    write_all(out, res, [], [], "# E7 summary (SMOKE: synthetic data)")


# ----------------------------------------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bootstrap", type=int, default=2000)
    ap.add_argument("--seed", type=int, default=config.SEED)
    ap.add_argument("--cap-per-ec", type=int, default=CAP)
    ap.add_argument("--pilot", action="store_true", help="New-392 (development set) only; writes E7_pilot")
    ap.add_argument("--smoke", action="store_true", help="synthetic data with known answers; writes E7_smoke")
    ap.add_argument("--allow-dirty", action="store_true",
                    help="development and tests only: run from uncommitted code, outputs labelled NOT CONFIRMATORY")
    ap.add_argument("--skip-ecbench", action="store_true",
                    help="tests only: leave out the EC-Bench secondary S7d (its files are pinned downloads)")
    args = ap.parse_args()
    config.ensure_dirs()
    t0 = time.time()
    tag = "E7" + ("_smoke" if args.smoke else "_pilot" if args.pilot else "")
    out = config.CX_RESULTS / tag
    out.mkdir(parents=True, exist_ok=True)
    if args.smoke:
        smoke(out, min(args.bootstrap, 400))
        print(f"\nCX13 SMOKE {'FAILED' if fails else 'DONE'} in {time.time() - t0:.0f}s -> {out}")
        sys.exit(1 if fails else 0)

    if args.pilot:
        label = "PILOT on New-392 (development set; descriptive)"
        state = freeze.code_state(freeze.app_paths(*FROZEN))
        dataset = "new392"
    else:
        label = freeze.require_frozen(freeze.app_paths(*FROZEN), args.allow_dirty, "E7")
        state = freeze.code_state(freeze.app_paths(*FROZEN))
        dataset = "sp_temporal_test"
    confirmatory = label.startswith("confirmatory")
    data = load_dataset(dataset)
    res = analyse(dataset, data, args.bootstrap, args.seed, args.cap_per_ec, confirmatory)
    extra_comp, extra_codes = [], []
    if not args.pilot:
        for other in ("new392", "sp_temporal_cal"):
            c, k = composition_only(other)
            extra_comp += c
            extra_codes += k
        if args.skip_ecbench:
            res["endpoints"].append({"endpoint": "S7d", "family": "S7d", "note": "skipped (--skip-ecbench)"})
        else:
            rows, comp, codes, gate_ecb = ecbench_secondary(args.bootstrap, args.seed)
            if rows and not gate_ecb:
                for r in rows:
                    r["note"] = f"descriptive: fewer than {GATE_ECB} EXP or PROPAGATED proteins"
            res["endpoints"] += rows
            extra_comp += comp
            extra_codes += codes
    manifest = L.load_json(PROV_DIR / "manifest.json") if (PROV_DIR / "manifest.json").is_file() else {}
    rel = sorted({r for s in manifest.get("sets", {}).values() for r in s.get("uniprot_releases", [])})
    digests = freeze.digests(freeze.app_paths(*FROZEN))
    header = (f"# E7 summary ({dataset}; {label})\n\n- code commit {state['commit']}; analysis code SHA-256 "
              + ", ".join(f"{k} {v[:12]}" for k, v in digests.items())
              + f"\n- UniProt release(s) of the evidence: {rel or 'unknown'}; B = {args.bootstrap}, seed {args.seed}")
    write_all(out, res, extra_comp, extra_codes, header)
    L.save_json({"label": label, "commit": state["commit"], "dirty": state["dirty"], "digests": digests,
                 "uniprot_releases": rel, "B": args.bootstrap, "seed": args.seed, "verdict": res["verdict"],
                 "gate": res["gate"], "counts": res["counts"]}, out / "E7_run.json")
    print(open(out / "E7_summary.md", encoding="utf-8").read())
    print(f"CX13 DONE in {time.time() - t0:.0f}s -> {out}")


if __name__ == "__main__":
    main()
