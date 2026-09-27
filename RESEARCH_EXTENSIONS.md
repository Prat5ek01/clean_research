# Research Investigation: Novel, Publishable Extensions to CLEAN

> Companion to [`PROJECT_OVERVIEW.md`](PROJECT_OVERVIEW.md) (architecture/repo map) and
> [`app/AUDIT_NOTES.md`](app/AUDIT_NOTES.md) + [`app/audit_out/`](app/audit_out/) (this repo's own
> evaluation findings). Written 2026-09-27 from (1) a full read of this repo's source/data/audit
> outputs and (2) a live web literature search conducted the same day. **The field around CLEAN
> moves fast — several of the closest competitor papers below were posted within the last few
> months of this writing. Re-verify status/citations directly before committing real work to any
> direction.**

---

## 0. Headline finding, stated up front

**The "obvious" extensions to CLEAN are almost all already published**, several within the last
year, some explicitly building on or benchmarking against CLEAN by name. Before this
investigation, the natural reading of this repo's own audit (`AUDIT_NOTES.md`,
`audit_out/*.csv`) suggests fruitful directions: "add EC-hierarchy structure to the loss,"
"fix singleton/rare classes," "calibrate the confidence score," "add structural/chemical
information." **All four are already done, in print, by named methods** (§2). This doesn't mean
there's nothing left — it means the *entry points* need to shift from "do the obvious thing" to
"find the specific slice of the obvious thing that the existing work has explicitly left open, or
attack it from an angle (interpretability, evaluation methodology, theoretical assumption) that
the architecture-race hasn't touched." §5 and §6 identify three such slices.

---

## 1. What this repo's audit actually established (and what it didn't)

From `AUDIT_NOTES.md` and `audit_out/`, this local audit found, on New-392 with max-separation:

| Finding | This repo's number | Is it a genuine finding? |
|---|---|---|
| Depth-1 (family) F1 ≫ depth-4 (exact reaction) F1 | 0.939 vs 0.498 | **Real pattern, but not novel to report** — independently confirmed by multiple 2024–2026 papers (§2), though their exact numbers differ from this repo's because of a methodological choice explained in §1.1 below. |
| Singleton-tier depth-4 F1 collapses | 0.152 (vs 0.497 aggregate) | **Real and important**, and conceptually corroborated by EnzPlacer's whole premise and GloEC's stated limitation of "insufficient rare EC coverage" — but no one I found reproduces this *exact* tier-crossed-with-depth table. It's solid supporting evidence, not yet a contribution on its own. |
| Score separation is weak (correct mean dist 4.33 vs wrong mean 5.99, overlapping) | — | **Real**, and this is the one finding that points most directly at a still-open, specific gap (§5.3). |

### 1.1 A concrete, actionable correction to make before doing anything else

`app/score_audit.py`'s `micro_prf()` computes **micro-averaged** precision/recall/F1 only. The
closest published comparator, **GloEC** (Huang et al. 2024, *Briefings in Bioinformatics*,
[10.1093/bib/bbae365](https://doi.org/10.1093/bib/bbae365)), reports **macro-averaged** F1 by EC
level and finds CLEAN's depth-4 macro-F1 is **23.9%** on their New-432 benchmark — dramatically
worse than this repo's own micro-F1 of 49.8% on New-392. Micro-averaging is dominated by the
common-EC samples (which do fine); macro-averaging weights every EC number equally, so it's
punished hard by exactly the long tail this repo's own `ec_train_frequency.csv` shows dominates
the label space (67% of training EC numbers are seen ≤10 times). **Before any new method work,
add macro-F1 (per depth, per tier) to `score_audit.py`.** This is a ~30-minute change, makes the
audit directly comparable to GloEC/EnzPlacer/CARE's own reported numbers, and on its own will
likely make the singleton-collapse story *more* dramatic, not less — worth knowing before it's
built into a paper's headline claim.

---

## 2. The landscape: who has already done what (as of this search)

| Method / paper | What it does | Relation to CLEAN | Status |
|---|---|---|---|
| **CLEAN-Contact** — Wang et al., *Communications Biology* 2024 (published from [bioRxiv 2024.05.14.594148](https://www.biorxiv.org/content/10.1101/2024.05.14.594148); DOI `s42003-024-07359-z`); code [PNNL-CompBio/CLEAN-Contact](https://github.com/PNNL-CompBio/CLEAN-Contact) | Adds a structural signal (2D residue **contact maps**, not full 3D structure) alongside the ESM sequence embedding into the contrastive objective. | Direct fork/extension of CLEAN, same authorship lineage area (PNNL). | **Closes** the naive "add structure" idea (Dimension 6). Its own stated limitation: contact maps are a lossy 2D proxy for 3D structure — a gap, but a narrow one. |
| **GloEC** — Huang et al., *Briefings in Bioinformatics* 2024, [10.1093/bib/bbae365](https://doi.org/10.1093/bib/bbae365); code [hyr0771/GloEC](https://github.com/hyr0771/GloEC) | Models the EC hierarchy explicitly as a **directed graph** with a 3-layer **Hierarchy-GCN** encoder (bottom-up + top-down label-dependency propagation), not contrastive learning. | Direct empirical comparator to CLEAN; beats it on macro-F1 at every EC level on New-432 and COFACTOR-237. | **Closes** most of "use the EC hierarchy structure" (Dimension 2) via a non-contrastive route. Self-stated limits: single-label only (no multifunctional enzymes), "insufficient coverage of rare EC numbers." |
| **EnzPlacer** ("How not to be seen: predicting unseen enzyme functions using contrastive learning") — Ma, Joshi, Friedberg, Li, *Bioinformatics* 42(Suppl. 1), 2026 (ISMB proceedings), [10.1093/bioinformatics/btag215](https://doi.org/10.1093/bioinformatics/btag215); code [drxiangma/EnzPlacer](https://github.com/drxiangma/EnzPlacer) | A **hierarchical exemplar contrastive loss (HiNCE)**: standard SupCon instance term *plus* a term pulling embeddings toward leave-one-out centroids computed at every EC depth (1–4). Singleton classes get synthetic positives via **residue masking + re-embedding** (different from CLEAN's sequence-mutation approach). Periodic hard-negative refresh from the evolving embedding space. | **This is the closest possible paper to "extend CLEAN's contrastive loss to fix hierarchy + rare classes."** Explicitly built as a CLEAN successor; beats CLEAN on unseen-EC4 EC2/EC3 accuracy (0.435 vs 0.385 at EC2) and on low-identity (<10%) subsets. | **Effectively closes Dimensions 2 and 3 combined** as a "modify the contrastive loss for hierarchy + rare classes" idea. Its own explicitly stated open problems (straight from the paper): (a) **no calibrated uncertainty/abstention mechanism for open-set use**, (b) multi-label/promiscuous enzymes not modeled, (c) no structure/substrate/cofactor signal, (d) EC3 macro-F1 "low across the board" even for EnzPlacer itself. |
| **EnzymeHunter** — Cao, Ouyang, Xiong, Wu et al., 2024–2026 ([bioRxiv 2024.07.07.602424](https://www.biorxiv.org/content/10.1101/2024.07.07.602424); code [cgxbio/EnzymeHunter](https://github.com/cgxbio/EnzymeHunter)) | Hierarchically-aware contrastive loss **plus** structural information together, explicitly targeting "low-homology" and "rare enzyme classes" performance; validated with a real wet-lab-confirmed novel annotation in *T. thermophilus*. | Same lineage of "CLEAN + hierarchy + structure" idea. | **Further closes** the combined hierarchy+structure+rare-class space — this is essentially the maximal version of the "obvious" extension, already built and validated with a real biological finding. |
| **CPEC** — Ding, Luo, Luo, *PLOS Computational Biology* 2024, [10.1371/journal.pcbi.1012135](https://doi.org/10.1371/journal.pcbi.1012135) | Wraps a contrastive embedding-distance model (**PenLight2**, not CLEAN) in a **conformal-prediction / "Learn-then-Test" FDR-control** layer: picks a distance threshold guaranteeing predicted-set false-discovery rate ≤ α with high probability. Uses an **adaptive cascade** (try EC4, fall back to EC3 under strict α). | Different backbone, but the method is backbone-agnostic and directly transferable to CLEAN's own distance outputs. | **Substantially narrows** "add calibrated confidence / abstention to CLEAN" (Dimension 5). What it does **not** do: stratify calibration by training-set rarity tier — global FDR control is proven only under an i.i.d./exchangeability assumption the authors explicitly flag as violated under covariate shift (e.g., rare/unseen classes, exactly this repo's weak point). |
| **HIT-EC** — *Nature Communications* (from [bioRxiv 2025.02.01.635810](https://www.biorxiv.org/content/10.1101/2025.02.01.635810)) and a separate **evidential-deep-learning EC model**, Han et al. (ECPICK), *Briefings in Bioinformatics* 2024, [10.1093/bib/bbad401](https://doi.org/10.1093/bib/bbad401) | Hierarchical transformer / CNN classifiers (not contrastive-embedding models) aimed at "trustworthy"/interpretable EC prediction, with residue-level importance scores. | Different model family (classifier, not nearest-cluster-center). I was only able to confirm ECPICK's method in detail (CNN + hierarchical global/local flow + residue importance scores as interpretability evidence, **not** formal Bayesian/Dirichlet evidential uncertainty despite the name overlap) — **verify HIT-EC's exact mechanism directly before relying on it**, I only have search-snippet-level confidence there. | Relevant prior art for "trustworthy EC prediction" generally, but a different architectural family from CLEAN's embedding-distance approach — leaves room for a **distance/embedding-native** calibration method specifically. |
| **CARE benchmark** — Yang, Mora, Liu, Wittmann, Anandkumar, Arnold, Yue (Caltech/Arnold lab), [arXiv:2406.15669](https://arxiv.org/abs/2406.15669) 2024; code [jsunn-y/CARE](https://github.com/jsunn-y/CARE) | A standardized benchmark: EC classification at 4 identity-stratified splits (<30%, 30–50%, previously-misclassified, promiscuous) **plus** an entirely separate enzyme-*from*-reaction retrieval task (easy/medium/hard by reaction similarity). Reports CLEAN's numbers directly: 55.1% (Level-4 acc, <30% identity), 80.2% (30–50%), 31.8% (misclassified), 69.4% (promiscuous). | Rigorous, large-scale, well-resourced independent benchmark of CLEAN and peers. | **Substantially narrows** "do a more rigorous stratified evaluation of CLEAN" (Dimension 7) as a standalone contribution — this has already been done, more rigorously and at larger scale, by a strong group. **But CARE explicitly flags an open gap that this repo is well-positioned to fill**: *"reaction equivalence wasn't thoroughly analyzed"* — i.e., no existing benchmark weights a wrong EC prediction by how chemically close it was to the true reaction (see §5.1). |
| **SOLVE** — bioRxiv 2025.01.28.635225 / *Chemical Science* 2025, [10.1039/D5SC04513D](https://doi.org/10.1039/D5SC04513D) | Ensemble (RF + LightGBM + decision tree, focal loss) that explicitly does **enzyme-vs-non-enzyme discrimination** before EC assignment, plus SHAP-based interpretability. Beats CLEAN/ECPred on an independent 2023+ UniProt set (F1 0.82 vs CLEAN's lower score). | Solves a real, separately-documented CLEAN weakness ("CLEAN indiscriminately assigns an EC number" to non-enzyme input) — but via a completely different (non-embedding, tree-ensemble) backbone. | **Narrows but doesn't close** an embedding-native version of open-set/non-enzyme rejection (i.e., doing it *inside* CLEAN's own distance space rather than a separate classifier). |
| **Multi-Alignment Contrastive Learning for Enzyme–Reaction Retrieval** — Zhou et al., [arXiv:2512.08508](https://arxiv.org/abs/2512.08508), and the related **TIGER** ([arXiv:2605.24489](https://arxiv.org/pdf/2605.24489)) | Contrastive alignment **between enzyme embeddings and reaction embeddings** (not just enzyme-to-enzyme), with a Gromov–Wasserstein-style geometric-consistency regularizer, for enzyme-for-reaction retrieval. | Different task framing (retrieval, not EC classification) but directly relevant to any "bring reaction chemistry into the loss" idea (Dimension 6). | Shows the "align with reaction-side chemistry" space is **already being worked by a well-resourced group**, very recently. Any "reaction-informed hard negatives" direction (§5.3-alt) needs to be positioned carefully against this. |
| Sparse-autoencoder interpretability of ESMC features — [arXiv:2606.12209](https://arxiv.org/pdf/2606.12209) | Applies mechanistic-interpretability tooling (SAEs) to protein-LM features "across the microbial protein universe," including enzyme function. | Orthogonal to all of the above — an *interpretability*, not *accuracy*, contribution. | **Thin/fresh territory.** I found no paper applying SAE-style mechanistic interpretability specifically to *why* CLEAN-style embeddings fail at exact-reaction-level or singleton-class prediction. Flagged as "potentially underexplored" per the hedging this investigation was asked to use — not verified as absent, just not found. |
| Convergent evolution of enzymes — Riziotis et al., *FEBS Journal* 2025 (from [bioRxiv 2024.04.08.588552](https://www.biorxiv.org/content/10.1101/2024.04.08.588552)); background: M-CSA database | Documents (with real curated cases) that enzymes sharing an EC number can have **no detectable sequence, fold, or active-site similarity** ("non-homologous isofunctional enzymes") — a real biological phenomenon, not a data artifact. | Not a CLEAN paper at all — a biology paper that undermines a **modeling assumption** every EC-cluster-center method (CLEAN, GloEC, EnzPlacer, EnzymeHunter) makes: one centroid per EC number. | I found **no** EC-prediction paper that explicitly tests or corrects for this (§5.2). This is a theoretical gap, not an engineering one — worth checking directly before committing, since absence-from-search is weak evidence. |

**Read as a whole:** the "architecture/loss" race (hierarchy-aware contrastive losses, structure
fusion, rare-class augmentation) is dense and being actively pushed by multiple groups with more
compute and more recent embeddings than this repo has access to. Competing there head-on (Level 3
in the requested framework) is high-risk for a resource-constrained solo effort. The gaps that
survive this search are narrower and more specific: (a) an **evaluation-methodology** gap CARE
itself names as unsolved, (b) a **theoretical/modeling-assumption** gap nobody in this search
addresses, (c) a **calibration** gap that exists on a different backbone (PenLight2/CPEC) than
CLEAN's own, explicitly under-addressed for rarity-stratified exchangeability.

---

## 3. Per-dimension status (mapping to the 7 dimensions requested)

| Dimension | Status | Key closing work | What's left |
|---|---|---|---|
| 1. Representation learning (alt. PLMs, fine-tuning, multi-level) | Partially explored | BMC Bioinformatics comparative-PLM study; GloEC/EnzymeHunter use ESM-1b too | Swapping in ESM2/ESM-C alone is a Level-1 hyperparameter-style change, not a contribution by itself |
| 2. EC hierarchy / reaction-level joint prediction | **Saturated** | GloEC (GCN), EnzPlacer (HiNCE), MAPred (autoregressive digit prediction), HIT-EC (hierarchical transformer) | Little room left for "just add hierarchy to the loss/architecture" |
| 3. Rare/singleton/unseen EC classes | **Saturated as a target, but not as a diagnosis** | EnzPlacer is built entirely around this | The *mechanism* of the singleton failure (is it embedding geometry? insufficient augmentation? label noise?) is under-diagnosed everywhere I looked |
| 4. Hard-negative mining | Partially explored, chemistry-blind everywhere found | CLEAN's own `mine_hard_negative`, EnzPlacer's periodic refresh — both purely embedding-space | **No paper found** that mines negatives using independent reaction/substrate chemical similarity rather than the model's own current embedding distances (see §5.3-alt; needs verification against the very recent Multi-Alignment CL paper before committing) |
| 5. Inference / calibration / uncertainty / open-set | Partially explored, backbone-fragmented | CPEC (conformal, PenLight2 backbone), SOLVE (enzyme/non-enzyme gate, tree ensemble), HIT-EC/ECPICK (classifier-based "trustworthy" framing) | **No paper found doing rarity/depth-stratified, exchangeability-aware conformal calibration natively on CLEAN's own max-separation distances** — this is the most concrete surviving gap (§5.3) |
| 6. Biological/chemical information (structure, active sites, reactions) | **Saturated** for structure; active for reactions | CLEAN-Contact, EnzymeHunter (structure); Multi-Alignment CL, TIGER (reaction retrieval); ECPICK/PARSE/ProtDETR/EC-LMGraph (catalytic-residue interpretability) | Reaction-chemistry-*informed evaluation* (not retrieval) is the one sub-slice CARE names as unsolved |
| 7. Rigorous/leakage-free evaluation | Done at scale by CARE | CARE (Caltech/Arnold group) | CARE's own named gap: reaction-similarity-weighted error analysis ("reaction equivalence wasn't thoroughly analyzed") — see §5.1 |

---

## 4. Level classification of concrete candidate ideas

**Level 1 (incremental — useful, not a paper alone):**
- Swapping ESM-1b → ESM2/ESM-C embeddings and re-running the pipeline.
- Adding macro-F1 alongside the existing micro-F1 in `score_audit.py` (do this regardless — see §1.1).
- Running `score_audit.py --preds results/new_pvalue.csv --tag pvalue` to close this repo's own
  noted gap between the two inference algorithms.
- Tuning `p_value`/`nk_random` or the GMM ensemble size.

**Level 2 (methodological extension — a real, scoped contribution):**
- §5.1 Chemistry-weighted evaluation metric for EC mispredictions (reaction/substrate similarity
  rather than exact/nothing).
- §5.3 Rarity/depth-stratified conformal risk control on CLEAN's own distances, positioned
  explicitly against CPEC.
- A CLEAN-native (embedding-space) enzyme/non-enzyme open-set gate, positioned against SOLVE.

**Level 3 (stronger, riskier, higher-compute contribution):**
- §5.2 Multi-prototype-per-EC modeling to handle convergent-evolution/non-homologous isofunctional
  enzymes.
- Chemistry-informed (not purely embedding-space) hard-negative mining requiring retraining —
  needs careful positioning against the Multi-Alignment CL / TIGER line of work before committing
  real compute.

---

## 5. Shortlist for deeper investigation

### 5.1 A chemistry-weighted "how wrong was the wrong answer" evaluation framework

1. **Research question.** When CLEAN (or any EC classifier) predicts the wrong exact EC number,
   is the error a chemically-near-miss (adjacent reaction, same superfamily) or an arbitrary
   mistake — and does today's exact-match scoring (this repo's own `score_audit.py` included)
   systematically overstate how "wrong" depth-4 errors really are?
2. **Existing limitation it addresses.** CARE ([arXiv:2406.15669](https://arxiv.org/abs/2406.15669))
   explicitly states: *"Similarity ambiguities: Similar reactions may carry different ECs;
   reaction equivalence wasn't thoroughly analyzed."* This repo's own audit_out tables also treat
   every off-target EC at depth 4 as an equally bad false positive.
3. **Related published work.** CARE (evaluation benchmark, names the gap); DRFP — Probst & Reymond,
   *Digital Discovery* 2022, [10.1039/D1DD00006C](https://doi.org/10.1039/D1DD00006C) (a
   data-independent, cheap reaction fingerprint usable without training a model); RHEA/BRENDA as
   sources of EC-to-reaction-SMILES mappings; Multi-Alignment CL / TIGER (reaction embeddings, but
   for retrieval, not error-weighting).
4. **Proposed method (conceptual).** For each EC number in the training/test data, pull its
   catalyzed reaction(s) (RHEA has curated reaction-EC mappings), compute a reaction fingerprint
   (DRFP is a strong, training-free default), and define a continuous **chemical distance between
   the true EC and the predicted EC**. Recompute the existing depth/tier tables not just as
   exact-match F1 but as a **distance-weighted score** (e.g., mean reaction-fingerprint distance of
   wrong predictions, or a soft-credit F1 that gives partial credit within some similarity radius).
   Apply this to CLEAN's existing `results/new_maxsep.csv`/`new_pvalue.csv` predictions already in
   this repo — no retraining required for the diagnostic.
5. **Code changes needed.** A new module (e.g. `app/reaction_distance.py`) mapping EC → reaction
   SMILES → DRFP vector; extend `score_audit.py`'s reporting functions with a
   `report_chemical_distance()` alongside the existing `report_depth`/`report_tier`. No changes to
   `src/CLEAN/` needed at all for the diagnostic version.
6. **Data/compute needed.** RHEA (or BRENDA) EC→reaction mapping (free download), `rdkit` +
   `drfp` (pip-installable). No GPU required. This is the cheapest direction on the list.
7. **Experimental validation plan.** (a) Compute the metric on this repo's existing predictions;
   (b) sanity-check against known convergent-evolution cases (§5.2) — those *should* show up as
   chemically-close "near-miss" errors if the metric works; (c) compare CLEAN vs. a naive baseline
   (e.g. BLAST top-hit) under the new metric to see whether the *ranking* of methods changes once
   near-misses are credited — that reordering (or lack of it) is itself the finding.
8. **Main risks / failure modes.** Reaction-SMILES coverage is incomplete for rare/singleton ECs
   (exactly the tier this audit cares about) — the metric could end up undefined precisely where
   it's most needed; EC-to-reaction mapping is sometimes many-to-many (promiscuous enzymes,
   multi-step reactions), so the "distance" isn't always well-defined either.
9. **How to distinguish from existing work.** CARE stratifies by *sequence* identity; this
   stratifies wrongness by *chemical* similarity of the predicted vs. true reaction — a
   genuinely different axis that CARE itself flags as unaddressed. Frame explicitly as "the
   evaluation gap CARE names, not a new model."
10. **Evidence needed before committing.** Check RHEA/BRENDA EC-reaction coverage against this
    repo's `split100.csv`/`new.csv` EC list first — if coverage is under, say, 60–70%, the metric
    won't be trustworthy enough to build a paper on, and this becomes a smaller methods-section
    contribution inside a larger paper rather than a standalone one.

### 5.2 Multi-prototype-per-EC modeling for convergent-evolution / non-homologous isofunctional enzymes

1. **Research question.** Does representing each EC number by a *single* cluster center (as
   `distance_map.py`'s `get_cluster_center` does today) systematically hurt exactly the EC numbers
   known to contain non-homologous isofunctional enzymes — and does detecting and
   splitting those into multiple prototypes recover accuracy without hurting the (majority of) EC
   numbers that really are one coherent cluster?
2. **Existing limitation it addresses.** A modeling assumption, not a benchmark result: CLEAN (and
   every hierarchy-aware successor found in §2) still uses one embedding centroid per EC number.
   Riziotis et al. (*FEBS Journal* 2025, from [bioRxiv 2024.04.08.588552](https://www.biorxiv.org/content/10.1101/2024.04.08.588552))
   and the M-CSA database document real, curated cases of enzymes sharing an EC number with **no
   detectable sequence/fold/active-site similarity**.
3. **Related published work.** General (non-enzyme-specific) multi-prototype / long-tailed
   prototype learning exists in computer vision (e.g., "Learning Prototype Classifiers for
   Long-Tailed Recognition," IJCAI 2023, [arXiv:2302.00491](https://arxiv.org/pdf/2302.00491)) but
   I found **no enzyme-function paper applying this to the convergent-evolution phenomenon
   specifically** — flagged per the required hedging as "potentially underexplored, requires
   further novelty verification," not confirmed absent.
4. **Proposed method (conceptual).** For every EC number with ≥ some minimum training count, run
   a within-class clustering pass (e.g., Gaussian mixture or k-means with model selection via
   BIC/silhouette) over the *existing* training embeddings already computed in this repo's
   `data/esm_data/` + `data/pretrained/100.pt`. EC numbers whose training sequences form one tight
   cluster keep a single centroid (status quo); EC numbers whose embeddings are genuinely
   multimodal get **k > 1 centroids**, and `infer_maxsep`'s nearest-center step matches against
   the nearest of the k centroids for that EC, not a single average.
5. **Code changes needed.** Modify `get_cluster_center` in `src/CLEAN/distance_map.py` (or add a
   parallel `get_cluster_centers_multi`), and adjust `infer_maxsep`'s distance-map construction in
   `src/CLEAN/infer.py` to handle a variable number of centroids per EC. This is a moderate,
   surgical change — no retraining of the projection head required, since it operates on already-
   trained embeddings.
6. **Data/compute needed.** None beyond what's already in this repo (`data/pretrained/100.pt`,
   `data/split100.csv`). CPU-only, cheap to run. Optionally cross-reference against M-CSA's curated
   convergent-evolution EC list to have ground-truth "known-multimodal" ECs for validation.
7. **Experimental validation plan.** (a) Flag which training EC numbers this repo's embeddings
   show as multimodal (silhouette/BIC test); (b) check overlap with M-CSA's documented convergent
   cases — if the *unsupervised* multimodality detector recovers *known* biological convergent-
   evolution cases, that's strong validation the signal is real, not noise; (c) re-run
   `score_audit.py`-style depth/tier evaluation with multi-prototype inference vs. the existing
   single-prototype baseline, specifically watching whether flagged EC numbers' accuracy improves
   without regressing the rest.
8. **Main risks / failure modes.** Embedding multimodality could reflect **within-species /
   technical clustering artifacts** (e.g., taxonomic sampling bias in training data) rather than
   genuine convergent evolution — the M-CSA cross-check in step 7 is the guard against this, and a
   negative result there (no overlap) would be an important, reportable negative finding, not a
   dead end.
9. **How to distinguish from existing work.** Every hierarchy/structure-aware CLEAN successor
   found (GloEC, EnzPlacer, EnzymeHunter, CLEAN-Contact) still collapses each EC number to one
   representation (a GCN node, a single contrastive class, etc.) — none were found to test
   *within-EC* multimodality as a first-class phenomenon.
10. **Evidence needed before committing.** A cheap 1–2 day spike: run the clustering pass on the
    existing `100.pt` embeddings for, say, the 50 EC numbers with the most training sequences, and
    check by eye/silhouette whether *any* show genuine bimodality before investing in the full
    M-CSA cross-reference and inference-pipeline changes.

### 5.3 Rarity/depth-stratified, exchangeability-aware conformal calibration for CLEAN's own distances

1. **Research question.** Can CLEAN's max-separation distances be turned into a **calibrated,
   abstention-capable** prediction with a formal false-discovery guarantee that holds *per rarity
   tier and per EC depth* — specifically targeting the singleton-tier failure this repo's own
   audit found (F1 0.152, and weak overall score separation: correct mean 4.33 vs. wrong mean
   5.99) — rather than one global threshold?
2. **Existing limitation it addresses.** This repo's own `score_audit.py`'s
   `report_confidence()` finding (weak, overlapping score separation) plus EnzPlacer's explicitly
   stated open problem ("need for calibrated uncertainty/abstention mechanisms for open-set
   applications").
3. **Related published work.** CPEC — Ding, Luo, Luo, *PLOS Comp Bio* 2024,
   [10.1371/journal.pcbi.1012135](https://doi.org/10.1371/journal.pcbi.1012135) — already does
   conformal FDR control on embedding distances with an EC-depth fallback cascade, but on a
   **different backbone** (PenLight2) and under a **single global** exchangeability assumption
   that the authors themselves flag as questionable under covariate shift (their own stated
   limitation). HIT-EC and ECPICK approach "trustworthiness" from a classifier-probability angle,
   not a distance/conformal one.
4. **Proposed method (conceptual).** Apply a CPEC-style Learn-then-Test procedure **natively to
   CLEAN's `infer_maxsep`/GMM distances**, but calibrate **separately per rarity tier**
   (common/rare/very_rare/singleton, using this repo's own `ec_train_frequency.csv` scheme) rather
   than globally — directly testing CPEC's own stated open concern (non-exchangeability across
   strata) as the central hypothesis, with the specific goal of surfacing when the model *should
   abstain* on likely-wrong singleton-tier calls rather than confidently emitting them.
5. **Code changes needed.** No changes to `src/CLEAN/` model code. A new calibration module
   consuming `results/*_maxsep.csv` distances + a held-out calibration split, implementing
   per-tier threshold search (Learn-then-Test / Hoeffding-based, following CPEC's published
   recipe) and an abstention option in the reporting layer of `score_audit.py`.
6. **Data/compute needed.** Needs a calibration split distinct from `new.csv`'s test role
   (currently used purely for evaluation) — likely a held-out slice of `split100.csv` or a fresh
   Swiss-Prot pull, since conformal calibration requires its own exchangeable sample. CPU-only.
7. **Experimental validation plan.** Report realized FDR vs. target α, **broken out by tier**,
   directly testing whether a single global threshold under- or over-covers specific tiers (the
   CPEC paper's own flagged risk) compared to the tier-stratified version; report
   coverage/abstention-rate trade-off curves per tier.
8. **Main risks / failure modes.** Splitting calibration by tier shrinks the calibration sample
   size per stratum severely for singleton/very-rare tiers (by definition, few examples) —
   conformal guarantees weaken or become vacuous exactly where they're needed most; this is a
   real, citable risk (small-sample conformal calibration is a known hard problem) and should be
   treated as a central experimental question, not glossed over.
9. **How to distinguish from existing work.** CPEC already did backbone-agnostic conformal
   distance calibration; the differentiator is (a) a different, more widely-used backbone (CLEAN
   itself, not PenLight2) and (b) treating tier-stratified exchangeability as the explicit research
   question CPEC raised but left open, rather than adding stratification as an afterthought.
10. **Evidence needed before committing.** Read the CPEC paper's methods section in full (this
    investigation only has an AI-summarized read) to confirm exactly how their Learn-then-Test
    procedure would or wouldn't already handle stratification via covariates — it's possible a
    straightforward extension already covers this and the gap is narrower (or wider) than it
    currently appears.

---

## 6. Comparison table

| Research Direction | Scientific Problem | Existing Literature | Proposed Contribution | Implementation Difficulty | Data Requirements | Experimental Test | Novelty Confidence | Publication Potential |
|---|---|---|---|---|---|---|---|---|
| §5.1 Chemistry-weighted error/evaluation metric | Exact-match scoring can't tell a near-miss from an arbitrary error | CARE (names this exact gap, doesn't fill it); DRFP | Reaction-fingerprint-distance-weighted re-scoring of existing EC predictions | Low | RHEA/BRENDA EC→reaction mapping, no GPU | Re-score existing `results/*.csv`; check against known convergent cases | Strong existing evidence of a research gap (CARE names it directly) | Promising as a methods contribution within a larger paper; possibly standalone as a short evaluation paper |
| §5.2 Multi-prototype-per-EC (convergent evolution) | Single centroid per EC is a modeling assumption known to be biologically false for some ECs | M-CSA / Riziotis et al. document the phenomenon; general (non-enzyme) multi-prototype ML literature exists | Detect within-EC embedding multimodality, use k-prototype nearest-center inference | Moderate | Existing repo embeddings only; optional M-CSA cross-reference | Unsupervised multimodality detection validated against curated convergent-evolution cases; before/after accuracy on flagged ECs | Promising but requires further investigation — absence of prior art not fully confirmed | Feasibility uncertain until the 1–2 day spike (§5.2 point 10) runs; could be a strong Level-2/3 paper if the spike shows signal |
| §5.3 Tier-stratified conformal calibration | Score separation is weak, especially for singleton/rare tiers; global calibration may not hold across strata | CPEC (same idea, different backbone, global-only, explicitly flags the stratification question as open) | Apply/extend CPEC's method to CLEAN's own distances, calibrated per rarity tier | Moderate | A held-out calibration split beyond `new.csv`; CPU-only | Realized-FDR-vs-target-α broken out by tier; coverage/abstention trade-off | Promising but requires further investigation (need to read CPEC methods in full first) | Feasible; risk is that small per-tier calibration samples weaken the guarantee — itself a reportable finding either way |
| (For contrast) Hierarchy-aware contrastive loss + rare-class augmentation | Depth-4 and singleton performance is weak | **Already solved**: EnzPlacer, EnzymeHunter, GloEC | — | — | — | — | Significant overlap with existing work | Low — would need a very sharp differentiator to be publishable now |
| (For contrast) Add structural information | Sequence-only embeddings miss 3D context | **Already solved**: CLEAN-Contact, EnzymeHunter | — | — | — | — | Significant overlap with existing work | Low, unless combined with one of the above in a genuinely new way |

---

## 7. Recommended staged workflow

1. **Repository understanding** — done (this document + `PROJECT_OVERVIEW.md`).
2. **Baseline validation** — before anything else: (a) add macro-F1 to `score_audit.py` (§1.1);
   (b) run the already-noted missing `--tag pvalue` audit pass; (c) don't just trust this repo's
   New-392 numbers as "the" baseline — CARE's numbers on the same model under identity-stratified
   splits are the more defensible baseline to cite going forward.
3. **Literature confirmation** — before committing to any of §5's three directions, do a fresh,
   targeted search (this investigation's searches are a snapshot, not exhaustive) specifically for:
   the exact phrase "conformal prediction CLEAN enzyme" and "multi-prototype EC number" to catch
   anything published in the weeks around this writing.
4. **Research hypothesis** — pick one of §5.1/5.2/5.3 (or a combination) and write it as a single
   falsifiable sentence before writing code (each direction's point 1 above is a draft of this).
5. **Small proof-of-concept** — for §5.1: score the *existing* `results/new_maxsep.csv` under the
   new metric (no new code beyond the metric itself, no GPU). For §5.2: the 1–2 day clustering
   spike on `100.pt`. For §5.3: a single global (non-stratified) conformal pass first, to confirm
   the CPEC recipe reproduces on CLEAN's distances before adding stratification.
6. **Controlled experiments** — compare against CLEAN's own unmodified max-separation/p-value
   outputs (already in `results/`) as the baseline in every case; for §5.1/5.3 also compare against
   the relevant published comparator's own reported numbers where reproducible.
7. **Ablation/robustness** — §5.1: does the finding hold if BRENDA is used instead of RHEA for the
   reaction mapping? §5.2: does the multimodality signal survive using ESM2 embeddings instead of
   ESM-1b? §5.3: does the tier-stratified result hold with a different calibration split size?
8. **Novelty/publication assessment** — re-run the exact literature searches in §2 one more time
   immediately before writing up, given how fast this specific sub-field has moved even during the
   course of this one investigation (three of the closest papers found were posted within months
   of this writing).

---

## 8. Guardrails carried over from the investigation brief (self-check)

- No paper, DOI, or result above was fabricated — every claim is attributed to a specific fetched
  or searched source; where I could only get search-snippet-level confidence (HIT-EC's exact
  mechanism), that's flagged explicitly rather than stated as fact.
- Higher-F1-alone is **not** treated as evidence of novelty anywhere above — the shortlist was
  selected because each item addresses a *specific stated gap* in a *specific* paper, not because
  it "sounds advanced."
- The biological meaning of EC numbers and the depth-1-vs-depth-4 distinction is treated as central
  throughout, not incidental.
- Singleton/rare-class performance is treated as a first-class concern (§5.2, §5.3), not an
  afterthought.
- Nothing here is guaranteed publishable; §6 uses descriptive confidence language, not fabricated
  scores, per the requested format.
