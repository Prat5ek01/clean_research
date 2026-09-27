# Publication plan: CLEAN audit and extension for *Computational Biology and Chemistry*

> Status: 2026-09-27. This supersedes the direction ranking in `RESEARCH_EXTENSIONS.md` (§15 lists
> what changed and why). How to run everything: [`app/cleanx/RUNBOOK.md`](app/cleanx/RUNBOOK.md).
>
> ✅ = run and verified in this repository today. ⏳ = code written and smoke-tested; waits only
> for the ESM-1b embedding step (RUNBOOK Step B, ~2–4 h CPU). Numbers marked *pilot* come from
> New-392, which was used for development and is disclosed as such in the paper.
> Nothing here guarantees acceptance. It is the most defensible paper the evidence supports, with
> decision gates so it survives negative results.

---

## 0. The paper in one paragraph

**Working title:** *"Weighted metrics, single centroids and silent overconfidence: auditing and
extending contrastive enzyme function prediction (CLEAN)"* (alternatives in §9.1).

We show that **(C1)** the metric CLEAN is reported with (support-weighted precision/F1 over the
test label set) does not see most false positives: 67% on New-392 and 96% on Price-149. It can
also be gamed: predicting CLEAN's top-10 ECs for every Price-149 protein gives a *reported* F1 of
0.618, above CLEAN's own 0.495. **(C2)** A new leakage-controlled temporal benchmark (2,099
Swiss-Prot enzymes created after CLEAN's training snapshot plus 2,000 non-enzymes) shows that
**30% of new enzymes carry an EC number CLEAN has never seen**, a failure mode New-392 and
Price-149 cannot expose (every EC in them is in CLEAN's vocabulary). **(C3)** Inference geometry
matters: nearest-member (exemplar) retrieval in CLEAN's own frozen embedding space helps common EC
classes but hurts rare ones, an order-statistic density bias. We test a class-balanced correction
(equal-prior kNN density) tuned only on a calibration split. **(C4)** CLEAN's released GMM
confidence is most overconfident exactly on rare EC numbers (pilot: 0.84 confidence at 34%
accuracy for very-rare-tier predictions). Tier-conditional calibration and certified selective
risk control (SGR) give per-tier guarantees and principled abstention, including on unseen-EC
enzymes and non-enzymes.

C1 and C2 are already verified and stand whatever C3/C4 show. C3 and C4 are fully implemented,
each with a pre-registered primary endpoint and a stated fallback framing.

---

## 1. Target journal: facts and requirements

| Item | Value | Status |
|---|---|---|
| Journal | *Computational Biology and Chemistry* (Elsevier; *Computers & Chemistry* until 2003), ISSN 1476-9271 | verified |
| Editors-in-Chief | Wentian Li; Donald Hamelberg | Wikipedia snapshot, re-check |
| Metrics | Impact factor 3.1 (2024); CiteScore 4.3 (2024) | verified via ScienceDirect/Wikipedia |
| Scope fit | "all areas of computational life sciences", explicitly including protein sequence research and ML in bioinformatics | verified |
| Access | Hybrid. Subscription route free; open access ≈ US$3,320 list price | third-party index, re-check |
| Abstract | **≤ 250 words** | verified (Guide for Authors, via search) |
| Highlights | Elsevier standard: 3–5 bullets, ≤ 85 characters each incl. spaces | **verify on the Guide** |
| Keywords, graphical abstract, reference style, length | not retrievable automatically | **verify on the Guide** |
| Declarations | competing interests; CRediT author roles; data availability; generative-AI use in writing (Elsevier-wide policy) | **verify exact wording on the Guide** |

The official Guide for Authors (`sciencedirect.com/journal/computational-biology-and-chemistry/publish/guide-for-authors`)
blocks automated access, so open it in your browser and tick §13's checklist against it.

**Why this journal fits:** CBAC publishes applied computational-biology method and evaluation work,
often ML on protein sequences. An evaluation + reliability + method paper with a new benchmark and
open code is a natural fit. A pure "new architecture beats SOTA" paper would compete with
better-resourced groups (§3).

---

## 2. Evidence already in hand

| # | Finding | Numbers | Source (in repo) |
|---|---|---|---|
| 1 | Exact reproduction of CLEAN | New-392 predictions 0/392 mismatches; P/R/F1 0.5955/0.4791/0.4967 (paper 0.596/0.479/0.497), bit-identical to CLEAN's `get_eval_metrics`; Price-149 GMM confidences identical to 4 d.p. | `cx_00_check.py` ✅ 9/9 PASS |
| 2 | **Blind spot of the reported precision** | False positives invisible to support-weighted precision: New-392 150/224 (67%); Price-149 106/110 (96%); p-value method 247/398 (62%) | `results/cleanx/E1/E1_blind_spot.csv` ✅ |
| 3 | Reported vs micro scores | New-392: P 0.596 vs 0.518, F1 0.497 vs 0.498. Price-149: P 0.584 vs 0.392, F1 0.495 vs 0.426. Union-macro F1: 0.371 / 0.238 | E1 ✅ |
| 4 | Metric ranking distortion | p-value vs max-separation precision gap: 0.04 reported, 0.14 micro | E1 ✅ |
| 5 | **Over-prediction games the metric** | Price-149, top-10 ECs per protein: reported F1 **0.618** (> CLEAN 0.495) while micro F1 0.127. New-392 top-3: reported F1 0.565 vs micro 0.370. True-label macro-F1 is gameable too (0.48 → 0.58) | `E1_overprediction.csv` ✅ |
| 6 | Leakage / label noise in New-392 | 5 proteins identical to a training sequence (2 correct; B0XPE4 wrong despite an identical training copy with the right EC; 2 Swiss-Prot annotation inconsistencies, P9WNP2/Q7U2S5 vs P9WNP3); 43 (11%) ≥ 90% identity to training | `E1_exact_duplicates.csv`, `data/cleanx/identity/new392.csv` ✅ |
| 7 | Label drift (complete ECs) | New-392: 385 unchanged, 3 changed, 2 gained an EC, 1 lost one, 1 no longer an enzyme entry. Price-149 IDs are RefSeq, so not applicable. ENZYME now lists 1,314 transferred and 178 deleted ECs | `E1_label_drift.csv` ✅ |
| 8 | **Temporal benchmark** | 8,920 Swiss-Prot entries created ≥ 2022-03-01 → 2,099 enzymes (839 calibration / 1,260 test) + 2,000 non-enzymes (sampled from 4,058 eligible); **640 enzymes (30%) carry ≥ 1 EC absent from CLEAN's training vocabulary** (353 distinct ECs, all valid ENZYME entries) | `data/cleanx/sets/sp_temporal_manifest.json` ✅ |
| 9 | Identity profile (fallback aligner) | Temporal test: 36% < 30% identity, 7.5% ≥ 90%. Price-149: 38% < 30%. New-392: 28% < 30%, 11% ≥ 90% | `data/cleanx/identity/*.csv` ✅ |
| 10 | **Exemplar trade-off** (dev + confirmatory) | Exemplar − centroid micro-F1: New-392 +0.035 [+0.002, +0.068]; Price-149 −0.012 [−0.048, +0.024]. By tier on Price-149: common **+0.062** (p=0.016), rare −0.085 (p=0.034), very rare **−0.132** (p<0.001), singleton −0.043 (p=0.024). New-392 common +0.059 (p=0.004) | `results/cleanx/E2_split100/E2_deltas.csv` ✅ |
| 11 | Local-mean (3-NN) inference | Union-macro F1: New-392 +0.032 [+0.007, +0.056]; Price-149 +0.006 (ns). K-means sub-prototypes slightly *worse* (negative result) | E2 ✅ |
| 12 | Homology baseline | Best local-alignment hit: ties CLEAN on New-392 (micro +0.007, ns), clearly worse on low-identity Price-149 (−0.074, p=0.038; union-macro −0.080, p=0.001) | E2 ✅ |
| 13 | **NISE-like ECs are dispersed in CLEAN space** | 183 of 2,493 ECs (7.3%) have members in ≥ 2 non-homologous Pfam-clan groups, consistent with the literature estimate of up to ~10% of reactions (Omelchenko et al. 2010) and with its class profile (mostly oxidoreductases 55 and hydrolases 54). **Controlling for EC size** (10,000 size-stratified permutations): within-EC RMS radius +0.84 (median 4.14 vs 3.44; AUROC 0.80) and margin to the nearest other centroid −1.00, both p < 10⁻⁴. Training-member misassignment is not significant after size control (p = 0.14). Test-set link on New-392 underpowered (25 NISE proteins; gain 0.035 vs 0.041, p = 0.58) | `results/cleanx/E3_split100/E3_stats.json` ✅ |
| 14 | **Tier-dependent overconfidence** (*pilot*, New-392 halves) | CLEAN GMM confidence, by predicted-EC tier: common ECE 0.18 (conf 0.64, acc 0.64); rare 0.36 (0.74 vs 0.48); very rare 0.50 (0.84 vs 0.34); singleton 0.49 (0.85 vs 0.36). Tier-conditional isotonic calibration lowers ECE in every tier | `results/cleanx/E4_smoke/` ✅ (pilot) |
| 15 | Certification limits (*pilot*) | Even CLEAN's 20 most confident predictions contain 2 errors. Smallest certifiable risk (centroid): common 0.35, rare 0.60, very rare 0.80, singleton 0.95 (tiny n; the temporal CAL split will tighten these) | E4 smoke ✅ (pilot) |

**What we could *not* verify and will not claim:** exact reported numbers of other methods
re-computed with our metrics (we only re-score CLEAN and our variants); the internals of EC-Bench's
metric code; and whether other groups' benchmark numbers suffer the same blind spot (we can only say
that any support-weighted average over the test label space does, by construction).

---

## 3. Novelty and gap analysis (literature checked 2026-09-27)

| Gap | Closest prior work | What it does / does not do | Our contribution | Novelty confidence |
|---|---|---|---|---|
| **G1 evaluation blind spot and gameability** | CLEAN (Yu et al. 2023) reports support-weighted P/R/F1; CLEAN-Contact's code copies CLEAN's `get_eval_metrics` verbatim (checked on GitHub); EC-Bench (Davoudi et al. 2026) reports weighted F1 at 4 EC levels; CARE (Yang et al. 2024) uses k=1 accuracy | None quantifies how many false positives support-weighted scores ignore, or shows they can be raised by over-predicting | Formal statement + quantification on two benchmarks + over-prediction stress test + reporting recommendations | **Strong**: no report found; verified numerically |
| **G2 temporal hold-out incl. unseen ECs & non-enzymes, for the released model** | New-392/Price-149 (CLEAN paper); EC-Bench temporal split (Swiss-Prot 2018 → 2023, models retrained); CARE identity splits; SOLVE's post-2023 set (505 enzymes + 694 non-enzymes) | None evaluates the *released* CLEAN on post-2022 entries with IUBMB-transfer mapping, identity strata and an explicit unseen-EC share | 2,099 + 2,000 benchmark, fully scripted and versioned (UniProt 2026_03, ENZYME 02-Sep-2026, Pfam 38.2) | Moderate (resource contribution; related sets exist) |
| **G3 single-centroid assumption** | CLEAN (centroids); EnzPlacer (kNN transfer + hierarchical loss, retrained; Ma et al. 2026); FEDKEA (distance-weighted kNN on ESM-2); LMKNN classifiers in general ML (Mitani & Hamamoto 2006) | No controlled comparison of inference geometries on the same frozen contrastive space; no analysis of the common-vs-rare trade-off or of NISE classes | Controlled geometry comparison; density-bias mechanism; class-balanced kNN-density scoring; NISE link | Moderate. kNN itself is not new; the controlled analysis, trade-off and correction are |
| **G4 rarity-dependent miscalibration; per-tier guarantees** | CPEC (Ding et al. 2024; conformal FDR on PenLight2, global threshold); HIT-EC, ECPICK (classifier "trust"); EnzPlacer names calibrated abstention as open | Nobody tests CLEAN's own released confidence by tier or gives per-tier guarantees | Tier-conditional calibration + SGR risk control + open-set abstention on unseen ECs and non-enzymes | Moderate–strong |
| G5 chemistry-weighted errors (earlier idea) | **EC-Bench already scores partial/additional EC recommendations with RDKit reaction similarity** | Substantial overlap | **Dropped** as a main claim; optional supplementary analysis only | Weak |
| G6 new training losses (hierarchy / rare classes) | EnzPlacer, EnzymeHunter, GloEC | Already done, and needs ESM-1b for all 227k training sequences | **Not pursued** (future work) | Weak |

Re-run these searches in the week before submission; this sub-field moves monthly:
`"CLEAN" enzyme weighted F1 precision evaluation`, `enzyme commission temporal benchmark 2025 2026`,
`conformal OR selective enzyme function prediction rare classes`,
`nearest neighbour vs centroid enzyme contrastive embedding`.

---

## 4. Research questions, hypotheses, endpoints

| RQ | Hypothesis | Endpoint (pre-registered) | Status |
|---|---|---|---|
| RQ1 How much of CLEAN's error does the reported metric see, and can it be gamed? | **H1** Support-weighted P/F1 ignore most FPs and rise under over-prediction; micro, union-macro and per-protein F1 do not | Share of FPs outside the test label space; score vs number of ECs predicted (1–10) | ✅ verified on legacy sets; ⏳ repeat on temporal test |
| RQ2 How does CLEAN generalise past 2022? | **H2** Exact-EC performance falls steeply with rarity tier and identity; a large share of new enzymes carry unseen ECs | Micro/union-macro F1 by depth × tier × identity bin; unseen-EC share | unseen share ✅ (30%); ⏳ performance |
| RQ3 Is one centroid per EC the right way to use CLEAN's embedding? | **H3a (secondary)** Exemplar 1-NN helps common and hurts rare tiers (replicates Price-149). **H3b (primary)** Class-balanced density scoring `balanced3` (λ tuned on CAL only) beats centroids on union-macro F1. **H3c (mechanism)** NISE-like ECs are more dispersed in CLEAN space (✅ size-controlled, p < 10⁻⁴) and exemplar gains concentrate in them (⏳ needs temporal power) | **Primary:** union-macro F1 at EC depth 4, `balanced3` − `centroid`, SP-temporal TEST, seen-EC proteins, paired bootstrap (B=2000), two-sided α=0.05. **Key secondary:** micro-F1, same comparison | ⏳ (H3c dispersion ✅) |
| RQ4 Is CLEAN's confidence reliable, and can abstention be certified per tier? | **H4a** GMM-confidence ECE rises from common to singleton predicted tiers. **H4b** Per-tier (Mondrian) SGR thresholds keep realized risk ≤ α in every tier on TEST; one global threshold does not. **H4c** The smallest certifiable risk differs by tier | ECE by predicted tier; realized selective risk and coverage per tier at α ∈ {0.1, 0.2, 0.3}; min certifiable risk; abstention on unseen-EC enzymes and non-enzymes | pilot ✅; ⏳ |

### Decision gates (what each outcome means for the paper)

- **G-A (after E2 on TEST):** if the primary CI excludes 0 in favour of `balanced3`, C3 is a positive
  method contribution. If not, C3 becomes a controlled negative result plus mechanism: "the
  exemplar/centroid trade-off is real and not fixed by a density correction". Put the emphasis on
  C1/C2/C4; the paper still stands.
- **G-B (after E4):** if Mondrian thresholds hold per-tier risk where the global threshold fails,
  C4 is a positive guarantee result. If rare tiers simply cannot be certified at useful α, report
  the *smallest certifiable risk per tier* as the finding: CLEAN should abstain on rare-EC calls
  unless more calibration data exist. That is still a clear, practical contribution.
- **G-C (after E3):** a significant NISE gain difference makes it a mechanism figure; otherwise it
  goes to the Supplement as a well-powered null.

---

## 5. Datasets

| Set | Size | Role | Notes |
|---|---|---|---|
| split100 (CLEAN training) | 227,362 sequences, 5,242 ECs, 241,025 (EC, seq) pairs | index for all inference geometries (released `100.pt` embeddings) | UniProt 2022_01 snapshot (newest regular entries created 2022-02-23) |
| split70 (CLEAN training, 70% clusters) | 79,172 pairs, 5,093 ECs | second backbone (robustness) | released `70.pt` + `split70.pth` |
| New-392 | 392 | **development** (used in piloting; disclosed) | 5 exact training duplicates; all ECs seen in training |
| Price-149 | 149 | confirmatory for the frozen centroid/exemplar/kmean3/multiproto; **post hoc** for `balanced*` | RefSeq IDs; low identity; all ECs seen |
| SP-temporal CAL | 839 enzymes | λ tuning (E2) and calibration (E4) **only** | stratified random 40% by tier |
| SP-temporal TEST | 1,260 enzymes (876 seen-EC, 384 unseen-EC) | **confirmatory** | untouched by any model selection |
| SP non-enzymes | 2,000 (of 4,058 eligible) | open-set rejection | no EC and none of the 7 enzyme-class keywords |
| NISE members (optional) | 3,817 training sequences of 168 NISE-like ECs (n ≤ 150) | E3b mechanism | FASTA already written |

**Temporal construction rules** (`cx_02_build_temporal.py`, recorded in the manifest): Swiss-Prot,
created ≥ 2022-03-01; not a fragment; length 30–1022; accession not in any CLEAN split, New-392,
Price-149 or halogenase; sequence not identical to any CLEAN training sequence; within-set duplicate
sequences removed; ≥ 1 complete EC (partial-only entries excluded). Labels are current complete
ECs. CLEAN's 2022 predictions are mapped through the ENZYME transfer table before scoring
(`labels.harmonize`). Rarity tier = count of the protein's rarest true EC in split100 (after
transfer mapping); `unseen` = count 0.

**Known limitation to state:** ENZYME records no transfer for some reclassifications (e.g. DNA
helicases: 3.6.4.12 stays active while UniProt moved many proteins to 5.6.2.3/5.6.2.4). Such
proteins land in the unseen tier and are analysed as an open-set problem, never as closed-set
errors.

---

## 6. Experiments

| ID | Question | Script | Status | Runtime (this PC) | Key outputs |
|---|---|---|---|---|---|
| E0 | Exact reproduction | `cx_00_check.py` | ✅ | 15 s | 9/9 PASS |
| E1 | Evaluation validity (RQ1) | `cx_01_audit.py` | ✅ legacy; ⏳ add `sp_temporal_test` | 10 s | blind spot, over-prediction, tier×depth, duplicates, drift |
| D | Temporal benchmark | `cx_02_build_temporal.py --with-pfam` | ✅ | 20 s (+8 min first Pfam table) | sets, FASTA, manifest |
| B | ESM-1b embeddings | `cx_03_embed.py` | ⏳ **your step** | smoke 5 min; full ≈ 4,099 × (printed s/protein) | `data/esm_data/*.pt` |
| I | Identity + homology baseline | `cx_04_identity.py` | ✅ (fallback); MMseqs2 optional | 13 min | identity bins, top-hit predictions |
| E2 | Inference geometry (RQ2/RQ3) | `cx_05_geometry.py` (+ `--backbone split70`) | ✅ legacy; ⏳ temporal | ~5 min (+3 min first multi-prototype build) | metrics, deltas, λ tuning, per-protein table |
| E3 | NISE mechanism (H3c) | `cx_06_ec_biology.py` (+ E3b `--subgroups`) | ✅ legacy link; ⏳ re-run after temporal E2 | 1–2 min | stats JSON, per-EC geometry |
| E4 | Reliability (RQ4) | `cx_07_reliability.py` | pilot ✅; ⏳ full | ~5 min | calibration, risk control, min certifiable risk, open set |
| F | Figures | `cx_08_figures.py` | partial ✅ | 30 s | PNG + PDF + CSV twins |

**How to run (details in the RUNBOOK):** from `app\`, run `.\cx_smoke.bat` (all 7 smoke tests,
≈2 min), then `.\cx.bat cx_03_embed.py --verify-only` (ESM-1b equivalence, needs ~8–10 GB free
RAM), then `.\cx_full.bat` (≈2–4 h, resumable).

**Verification done on 2026-09-27:** `cx_smoke.bat` → ALL 7 SMOKE TESTS PASSED; an integration
test with synthetic embeddings (temp folders only) exercised every temporal-only path (unseen ECs,
EC-transfer mapping, λ tuning on CAL, non-enzyme open set, reliability figures) → 5/5 PASS; unit
tests of the embedding orchestration (PASS and FAIL branches), the MMseqs2 parser and E3b → PASS.
Bugs found and fixed by these tests: an E3 summary key error and figure crashes on empty panels.
**Not executed here:** loading ESM-1b itself (only ~4 GB RAM free); its equivalence test is
Command 2 in the RUNBOOK and gates the full run.

**Statistics:** paired non-parametric bootstrap over proteins (B = 2,000; seed 20260927) for all
method differences; 95% percentile CIs; two-sided p. Per-tier secondary tests are Holm-corrected
across the 4 tiers. E3 uses a size-stratified permutation test (10,000 permutations). SGR uses
δ = 0.1, split by Bonferroni across the 4 predicted tiers for Mondrian thresholds.

**Metrics reported everywhere:** micro P/R/F1; union-macro F1 (every false-positive class counts);
per-protein (samples) F1; subset accuracy; and, for comparability with the literature, the
support-weighted and true-label macro scores, **explicitly marked as blind to out-of-label-space
false positives**.

---

## 7. Pre-registration and anti-leakage protocol

1. **Freeze before the first TEST run.** Commit the `app/cleanx/` package and `app/cx_*.py`, then
   write the commit hash into the Methods ("analysis code frozen at commit …"). The primary and
   secondary endpoints are written in `cx_05_geometry.py`'s docstring (dated 2026-09-27).
2. **CAL only for fitting**: λ (E2) and every calibrator/threshold (E4) are fitted on
   SP-temporal CAL. TEST is only scored.
3. **Disclosure:** New-392 was used for piloting (1-NN vs centroid, dispersion scaling (rejected),
   BIC multi-prototypes (rejected), GMM calibration). Price-149 was first scored on 2026-09-27 with
   frozen methods, and its per-tier result motivated `balanced*`, so `balanced*` on Price-149 is
   post hoc.
4. **Leakage controls:** accession + exact-sequence exclusion (done); identity strata (done;
   MMseqs2 confirmation recommended); PE-level sensitivity (Swiss-Prot "evidence at protein level":
   1,663 of 2,099 enzymes). Per-EC cap sensitivity if one EC dominates.
5. **No tuning on anything that is later reported as a test result.**
6. **Recorded endpoint change (2026-09-27, before any SP-temporal TEST result existed):** the
   primary metric and the λ-tuning criterion were switched from true-label macro-F1 to
   union-label macro-F1 after E1 showed true-label macro-F1 also rises under over-prediction
   (0.48 → 0.58 on New-392). State this in the Methods; it strengthens the protocol.

---

## 8. Figures and tables

| Item | Content | Produced by |
|---|---|---|
| Fig. 1 | Schematic / graphical abstract: frozen CLEAN encoder → four inference geometries → calibration/abstention; benchmark timeline | draw manually (PowerPoint/Inkscape) |
| Fig. 2 | Blind spot: reported vs micro precision; visible vs invisible FPs; + over-prediction curve (reported vs micro F1 against ECs predicted) | `fig1_blind_spot` (+ add over-prediction panel from `E1_overprediction.csv`) |
| Fig. 3 | Exact-EC F1 by true tier for centroid / exemplar / homology (3 datasets) | `fig2_tiers` |
| Fig. 4 | F1 vs identity bin | `fig3_identity` |
| Fig. 5 | NISE gain (exemplar − centroid) | `fig4_nise_gain` (+ E3b silhouettes) |
| Fig. 6 | Reliability diagrams by predicted tier: CLEAN GMM vs tier-conditional | `fig5_reliability_by_tier` |
| Fig. 7 | Risk–coverage curves; realized per-tier risk, global vs per-tier thresholds | `fig6_risk_coverage`, `fig7_tier_risk_control` |
| Table 1 | Datasets: sizes, tiers, identity bins, unseen share | manifest + identity CSVs |
| Table 2 | Main results with CIs (methods × datasets × metric families) | `E2_bootstrap_ci.csv` |
| Table 3 | Calibration (ECE) and min certifiable risk by tier | `E4_calibration.csv`, `E4_min_certifiable_risk.csv` |
| Supplement | full metric tables; duplicates; drift; over-prediction; split70 replication; MMseqs2 vs fallback identity; λ tuning curve; depth ablation (`figS1`); E3b; PE-level sensitivity | E-scripts |

Figure style is already fixed in `cx_08_figures.py`: validated colour-blind-safe palette (worst
CVD ΔE 9.2), fixed method → colour mapping, one y-axis per panel, hairline grids, and a CSV twin
per figure for the data statement.

---

## 9. Manuscript blueprint

### 9.1 Title options
1. *Weighted metrics, single centroids and silent overconfidence: auditing and extending contrastive enzyme function prediction*
2. *What reported scores hide: evaluation blind spots, inference geometry and tier-aware abstention in CLEAN*
3. *Beyond the centroid: class-balanced retrieval and certified abstention for contrastive EC-number prediction* (use only if G-A and G-B are both positive)

### 9.2 Structure and word budget (~6,500 words main text)
- **Introduction (~800):** why EC prediction matters; CLEAN and successors (CLEAN-Contact,
  EnzPlacer, GloEC, EnzymeHunter); benchmark efforts (CARE, EC-Bench); calibration (CPEC). Gap
  statement (§3). Numbered contributions C1–C4.
- **Methods (~2,000):** reproduction; formal metric definitions, including the one-line proof that
  support-weighted precision gives weight 0 to classes without true instances; temporal benchmark
  rules; identity; the four geometries and the balanced-density derivation; NISE definition;
  calibration, Mondrian thresholds, SGR; statistics; software versions and data releases.
- **Results (~2,500):** 3.1 reproduction + blind spot + over-prediction; 3.2 temporal
  generalisation (tier × depth × identity; 30% unseen); 3.3 inference geometry (trade-off, primary
  endpoint, split70 replication); 3.4 NISE mechanism; 3.5 reliability (ECE by tier, certified
  per-tier risk, min certifiable risk, open-set abstention).
- **Discussion (~900):** a boxed **reporting checklist** for EC-prediction papers (micro +
  union-macro + FP counts + per-tier + identity + unseen-EC share); implications for CLEAN users;
  limitations (frozen model, one model family, Swiss-Prot label noise, EC reclassification, small
  rare tiers); future work (retraining with balanced objectives; applying the audit to
  CLEAN-Contact/EnzPlacer embeddings).
- **Conclusion (~150).**

### 9.3 Abstract draft (fill the ⟦ ⟧ after E2/E4; keep ≤ 250 words)
> Contrastive learning models such as CLEAN assign Enzyme Commission (EC) numbers by comparing a
> protein's embedding with one centroid per EC class, and are evaluated with support-weighted
> precision and F1. We show that this metric ignores 67–96% of CLEAN's false positives on standard
> benchmarks and can be inflated by over-prediction: predicting ten EC numbers per protein raises
> the reported F1 on Price-149 from 0.495 to 0.618. We built a leakage-controlled temporal benchmark
> of 2,099 Swiss-Prot enzymes created after CLEAN's training snapshot and 2,000 non-enzymes; 30% of
> the new enzymes carry an EC number absent from CLEAN's vocabulary. On proteins with known EC
> numbers, exact-EC performance ⟦falls from … for common to … for singleton classes⟧. Replacing
> centroids by nearest-member retrieval in CLEAN's frozen embedding space improves common classes
> but harms rare ones, a class-size bias that ⟦a class-balanced density score corrects / does not
> correct⟧ (union-macro F1 ⟦…⟧). CLEAN's released confidence is most overconfident on rare EC
> numbers (expected calibration error ⟦…⟧ vs ⟦…⟧). Tier-conditional calibration with certified
> selective risk control ⟦keeps error ≤ α in every rarity tier / shows that rare-tier calls cannot
> be certified below …⟧ and abstains on ⟦…⟧% of unseen-EC enzymes and ⟦…⟧% of non-enzymes. We
> release the benchmark and analysis code and propose a reporting checklist for EC prediction.

### 9.4 Highlights (character counts checked; all ≤ 85)
1. `Reported precision hides 67-96% of CLEAN's false positives on standard benchmarks` (81)
2. `A 2022-2026 Swiss-Prot hold-out: 30% of new enzymes carry ECs CLEAN never saw` (77)
3. `CLEAN's released confidence is most overconfident on rare EC numbers` (68), after E4 confirms
4. `Tier-conditional risk control certifies abstention separately per rarity tier` (77), if G-B positive
5. `Exemplar retrieval in CLEAN's own embedding space improves exact-EC F1` (70). **Use only if
   G-A is positive**; otherwise replace with a trade-off highlight.

### 9.5 Keywords (≤ 6, confirm the limit)
Enzyme Commission number; enzyme function prediction; contrastive learning; evaluation metrics;
calibration; selective prediction.

---

## 10. Timeline (working days, with gates)

| Day | Work | Gate |
|---|---|---|
| 1 ✅ | E0, E1 legacy, benchmark build, identity, E2 legacy, E3 legacy, pilots, code | done |
| 2 | RUNBOOK Step B: ESM-1b verify (must PASS) → full embeddings (2–4 h) → E2 on TEST + split70 → E3 re-run | **G-A** |
| 3 | E4 full; E1 on TEST; figures; optional E3b embeddings (1.5–2 h) and MMseqs2 identity | **G-B, G-C** |
| 4–5 | Freeze numbers; Tables 1–3; Supplement; sensitivity analyses (PE level, per-EC cap) | numbers frozen |
| 6–12 | Write Methods + Results first, then Introduction/Discussion, abstract, highlights, graphical abstract | internal read-through |
| 13–15 | Licence steps (§13), Zenodo archive, final literature re-check, cover letter, submission | submit |

---

## 11. Risks and mitigations

| Risk | Likelihood | Mitigation |
|---|---|---|
| ESM-1b load fails (RAM) or is too slow | medium | close apps (another Claude session also runs jobs here); resumable; or run `cx_03` on a GPU server (copy `app/`, same commands) |
| Primary endpoint (H3b) negative | medium | pre-registered fallback framing (G-A); the trade-off itself is a replicated finding |
| Rare tiers cannot be certified | high | report the smallest certifiable risk per tier as the result (G-B) |
| Reviewer: "the metric point is trivial" | medium | quantify invisibility, show gameability (Price-149 0.618 > 0.495), show the ranking distortion, and give a checklist; the definition is standard, its consequence here was unreported |
| New overlapping paper appears | medium | re-run §3 searches before submission; adjust claims |
| Identity via fallback aligner questioned | medium | MMseqs2 re-run (RUNBOOK §3); report agreement |
| Temporal labels partly inferred | low–medium | PE-level sensitivity (1,663 protein-level entries) |
| **Licence** (CLEAN code redistribution) | real | §13: private repo; add-on release only; UIUC confirmation |
| Git mishaps (the aborted pull on 2026-09-27) | real | commit often; never `git pull --rebase` between unrelated histories; keep `data/cleanx/` and `results/` out of git |

---

## 12. Reviewer questions and prepared answers

- **"Why not retrain CLEAN?"** Scope: we keep the released model fixed to isolate evaluation,
  inference geometry and calibration. Retraining needs ESM-1b for all 227k sequences, and
  retraining-based fixes already exist (EnzPlacer, EnzymeHunter). Our analyses apply to any
  embedding model.
- **"Isn't 1-NN just BLAST?"** We compare against a homology top-hit baseline directly. It ties
  CLEAN at high identity and loses clearly at low identity (Price-149: −0.074 micro, −0.080
  union-macro).
- **"Weighted F1 is standard."** We don't call it wrong; we show what it cannot see and that it
  rewards over-prediction, and we recommend reporting it *alongside* micro/union-macro and FP
  counts.
- **"Why SGR and not conformal FDR like CPEC?"** Annotation pipelines ask "accept this top-1 call
  or abstain?", which is selective risk. CPEC's set-valued FDR is complementary; we discuss it and
  can add it as a supplementary analysis.
- **"Was anything tuned on the test set?"** No: §7, with the frozen commit hash and dated
  pre-registered endpoints.
- **"Is the temporal set leaking?"** Accession and exact-sequence exclusion, identity strata, 7.5%
  of TEST ≥ 90% identity is reported and analysed separately.

---

## 13. Submission, licence and ethics checklist

**CLEAN licence (read `NON-EXCLUSIVE RESEARCH USE LICENSE FOR CLEAN SOFTWARE.pdf`; not legal advice):**
- [ ] §1.1 counts *"distribution or otherwise making the Software available to a third party"* as
      Commercial Use, which is not permitted without a UIUC licence. **Keep `Prat5ek01/clean_project`
      private**, and release only the `cleanx` add-on (no CLEAN source, no weights, no `split*.csv`).
      Ask UIUC (otm@illinois.edu) before any public release; mention that users obtain CLEAN from
      the official repository.
- [ ] §2.5–2.6: derivative work renamed (`cleanx`) and carries the required notice (already in
      `app/cleanx/__init__.py`).
- [ ] §2.7.1: put this sentence in Acknowledgements verbatim: *"The CLEAN Software, used by [your
      name], was developed by the Department of Chemical and Biomolecular Engineering at the
      University of Illinois Urbana-Champaign."*
- [ ] §2.7.2: email Prof. Huimin Zhao (zhao5@illinois.edu) about the publication (on
      acceptance at the latest; a courtesy note at submission is reasonable).
- [ ] `app/build.py` says "Apache-2.0"; the PDF licence governs. Don't cite Apache.

**Journal (verify each against the Guide):** abstract ≤ 250 words · highlights · keywords ·
graphical abstract (recommended) · CRediT statement · competing-interest declaration · funding ·
data availability · **declaration of generative-AI use in the writing process** (Elsevier policy;
disclose AI assistance truthfully, and describe AI-assisted code in Methods if applicable) ·
cover letter · suggested reviewers (avoid anyone with a conflict).

**Reproducibility:** environment versions (RUNBOOK §0) · data releases (UniProt 2026_03, ENZYME
02-Sep-2026, Pfam 38.2, all recorded in `data/cleanx/downloads/*.json`) · weight md5s printed by
`cx_00_check.py` · frozen commit hash · Zenodo DOI for the `cleanx` add-on, benchmark accession
lists and result CSVs (UniProt data are CC BY 4.0: cite UniProt).

**Ethics:** no human or animal data; public databases only.

---

## 14. Literature to cite

Identifiers marked ✔ were seen on the source during this session; others are standard references to
double-check before submission.

| Reference | Identifier |
|---|---|
| Yu et al., Enzyme function prediction using contrastive learning, *Science* 2023 (CLEAN) | 10.1126/science.adf2465 ✔ |
| CLEAN-Contact, *Communications Biology* 2024 | 10.1038/s42003-024-07359-z ✔ (bioRxiv 10.1101/2024.05.14.594148) |
| Ma, Joshi, Friedberg, Li, How not to be seen (EnzPlacer), *Bioinformatics* 2026 | 10.1093/bioinformatics/btag215 ✔ |
| Huang et al., GloEC, *Brief. Bioinform.* 2024 | 10.1093/bib/bbae365 ✔ |
| Davoudi et al., EC-Bench, *Bioinformatics Advances* 2026 | 10.1093/bioadv/vbag004 ✔ |
| Yang et al., CARE benchmark, 2024 | arXiv:2406.15669 ✔ |
| Ding, Luo, Luo, CPEC conformal enzyme annotation, *PLOS Comput. Biol.* 2024 | 10.1371/journal.pcbi.1012135 ✔ |
| SOLVE, *Chemical Science* 2025 | 10.1039/D5SC04513D ✔ |
| Han et al., ECPICK, *Brief. Bioinform.* 2024 | 10.1093/bib/bbad401 ✔ |
| HIT-EC, *Nature Communications* 2026 | s41467-026-68727-3 ✔ |
| EnzymeHunter, 2026 | PMID 42453691 ✔ |
| Omelchenko et al., Non-homologous isofunctional enzymes, *Biology Direct* 2010 | PMID 20433725 ✔ |
| Riziotis et al., Paradigms of convergent evolution in enzymes, *FEBS J.* 2025 | 10.1111/febs.17332 ✔ |
| Geifman & El-Yaniv, Selective classification for deep neural networks, NeurIPS 2017 (SGR) | arXiv:1705.08500 (check) |
| Angelopoulos et al., Learn then Test, 2021 | arXiv:2110.01052 (check) |
| Vovk, Gammerman, Shafer, *Algorithmic Learning in a Random World* (Mondrian CP) | book (check edition) |
| Ding et al., Class-conditional conformal prediction with many classes, 2023 | arXiv:2306.09335 ✔ |
| Guo et al., On calibration of modern neural networks, ICML 2017 (ECE) | arXiv:1706.04599 (check) |
| Mitani & Hamamoto, A local mean-based nonparametric classifier, *Pattern Recogn. Lett.* 2006 (LMKNN) | (check DOI) |
| Rives et al., ESM-1b, *PNAS* 2021 | 10.1073/pnas.2016239118 (check) |
| Steinegger & Söding, MMseqs2, *Nat. Biotechnol.* 2017 | 10.1038/nbt.3988 (check) |
| UniProt Consortium (current NAR database-issue paper); Bairoch, ENZYME database; Pfam/InterPro (Pfam 38.2) | check latest |
| FEDKEA (distance-weighted kNN on ESM-2), bioRxiv 2024 | 10.1101/2024.08.12.604109 ✔ |

---

## 15. What changed vs `RESEARCH_EXTENSIONS.md` (and why)

- **Chemistry-weighted evaluation dropped as a main claim**: EC-Bench (2026) already scores EC
  recommendations by RDKit reaction similarity.
- **The GloEC macro-F1 comparison there needs care**: GloEC's CLEAN numbers are macro-averaged on a
  different set (New-432); our own micro/union-macro/true-macro numbers are now computed directly.
- **Training-space "misassignment" pilot (p = 3×10⁻⁷) was an artifact** of multi-EC sequences.
  With co-annotated ECs masked, the effect is weak. The NISE analysis therefore relies on test-set
  behaviour.
- **Dispersion-scaled centroids and BIC multi-prototypes** failed in pilots and are not in the
  method set (BIC split 1,614 of 1,811 ECs in 128-D; silhouette-gated k-means is kept as a
  documented negative ablation).
- **New since then:** the metric's gameability, the temporal benchmark (30% unseen ECs), the
  exemplar common-vs-rare trade-off, and the tier-dependent overconfidence of CLEAN's GMM
  confidence.
