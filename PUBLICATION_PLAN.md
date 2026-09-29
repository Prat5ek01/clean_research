# Publication plan: CLEAN audit and extension for *Computational Biology and Chemistry*

> Status: 2026-09-28, updated after the confirmatory run of the pre-registration commit c8e42f7 (results in
> §2b; post-run changes in §7 and §15). Supersedes the 2026-09-27 version and the direction ranking in
> `RESEARCH_EXTENSIONS.md`. How to run everything: [`app/cleanx/RUNBOOK.md`](app/cleanx/RUNBOOK.md).
>
> ✅ = run and verified. ⏳ = code written and tested (unit + synthetic end-to-end), waiting for the full run.
> *pilot* = development data, disclosed as such. Evidence levels in §3/§14: **[code]** = checked in the
> public source code at the stated commit; **[search]** = from search-engine results/abstracts only (full texts
> of PMC, OUP, Nature, bioRxiv and arXiv were not reachable from the analysis environment; re-check them).
> Nothing here guarantees acceptance. It is the most defensible paper the evidence supports, with decision
> gates so it survives negative results.

---

## 0. The paper in one paragraph

**Working title:** *"What reported scores hide: metric blind spots, inference geometry and the calibration
budget for trusting rare enzyme-function calls"* (alternatives in §9.1).

**(C1) Metric blind spot, formalised and shown field-wide.** We prove that EC scores computed after projecting
predictions onto the test label set, and support-weighted or true-label-macro scores on any label space, are
invariant to predictions of EC numbers absent from the test labels, and that support-weighted recall equals
micro recall (all distortion is in precision). The same scoring is used in the public code of CLEAN,
CLEAN-Contact, ProtDETR, EnzymeHunter and the EC-Bench benchmark [code]. It hides 67–96% of CLEAN's false
positives, is gameable on small benchmarks (top-10 predictions give Price-149 a reported F1 of 0.618 vs CLEAN's
0.495), and across EC-Bench's 16 published predictors the inflation tracks how many ECs a method predicts
(confirmatory ρ = 0.75) and reverses one published method comparison (pre-registered, q = 0.026).
**(C2) A leakage-controlled temporal benchmark** (2,099 Swiss-Prot enzymes created after CLEAN's snapshot +
2,000 non-enzymes): 30% of new enzymes carry an EC number CLEAN has never seen, and a model's training
vintage must be matched per model (a 2025 successor's training set already contains half of it).
**(C3) Inference geometry and the homology baseline:** on post-2022 enzymes, member-level inference in CLEAN's
frozen space beats its one-centroid-per-EC rule (pre-registered primary: `balanced3` +0.043 union-macro F1
[+0.030, +0.060]; exemplar +0.079, descriptive), yet a standard BLASTp top hit against the same training set beats
every CLEAN-space rule except on remote homologs (seen-EC micro-F1 0.735 vs 0.641 for CLEAN; crossover by identity
⟦final⟧). **(C4) Reliability
and its price:** CLEAN's released confidence is most overconfident on rare ECs (expected calibration error 0.32 in
the very-rare and singleton predicted tiers vs 0.19 in the common tier; tier-conditional isotonic recalibration on
the calibration split brings all four to 0.05–0.10); exact binomial certification budgets show that per-tier risk
control at 20% needs at least 550–630 calibration proteins in each rarer tier (lower 95% bounds; 90–122 exist),
at least 6–10 years of Swiss-Prot curation (pre-registered H6 and G-E met).

C1 and C2 stand whatever C3/C4 show. Every experiment has a pre-registered primary endpoint and a stated
fallback framing (§4, §7).

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

The official Guide for Authors blocks automated access: open it in your browser and tick §13 against it.
**Why this journal fits:** applied computational-biology method and evaluation work, often ML on protein
sequences. An evaluation + reliability + method paper with a new benchmark, a formal result and open code is a
natural fit; a "new architecture beats SOTA" paper would compete with better-resourced groups (§3).

---

## 2. Evidence already in hand

| # | Finding | Numbers | Source |
|---|---|---|---|
| 1 | Exact reproduction of CLEAN | New-392 0/392 prediction mismatches; P/R/F1 0.5955/0.4791/0.4967 (paper 0.596/0.479/0.497), bit-identical to CLEAN's `get_eval_metrics`; Price-149 GMM confidences identical to 4 d.p. | `cx_00_check.py` ✅ 9/9 (your PC, 2026-09-27) |
| 2 | **Blind spot of the reported precision** | False positives invisible to support-weighted precision: New-392 150/224 (67%); Price-149 106/110 (96%); p-value method 247/398 (62%) | `E1/E1_blind_spot.csv` ✅; re-derived by the E5 card ✅ |
| 3 | Reported vs micro scores | New-392: P 0.596 vs 0.518, F1 0.497 vs 0.498. Price-149: P 0.584 vs 0.392, F1 0.495 vs 0.426. Union-macro F1: 0.371 / 0.238 | E1 ✅ |
| 4 | Metric ranking distortion (CLEAN's two algorithms) | p-value vs max-separation precision gap: 0.04 reported, 0.14 micro | E1 ✅ |
| 5 | **Over-prediction games the metric** | Price-149, top-10 ECs per protein: reported F1 **0.618** (> CLEAN 0.495) while micro F1 0.127. New-392 top-3: reported F1 0.565 vs micro 0.370. True-label macro-F1 is gameable too (0.48 → 0.58) | `E1_overprediction.csv` ✅ |
| 6 | **Formal result** | Props. 1–2, Lemma 3, Prop. 4 (§9.2; `cleanx/evalcard.py`); verified numerically on 45 real prediction tables and in 21 unit tests | `cx_09_evalcard.py` ✅ |
| 7 | **The blind spot is field-wide** [code] | Scoring invariant to out-of-label predictions in CLEAN, CLEAN-Contact, ProtDETR, EnzymeHunter (verbatim copies of `get_eval_metrics`) and EC-Bench (support-weighted on the union label space); all false positives counted by CPEC (samples), CARE (top-1 accuracy), GloEC and SOLVE (macro on single-label outputs) | §3 table, file:line at commit |
| 8 | **EC-Bench re-scoring validated** | Our scorer reproduces EC-Bench's printed weighted F1 for all 16 methods at EC levels 1 and 4: max \|diff\| 0.0 | E5 gate ✅ |
| 9 | **Inflation tracks over-prediction** (*pilot*, EC-Bench Price-149) | Spearman ρ(reported − micro F1, ECs predicted per protein) = 0.94 (16 methods, 100% regime; 0.93 for 30%); median invisible share of FPs 0.96 (range 0.83–0.98; 0.93 for 30%); Kendall τ_b(reported, micro) 0.88 [0.80, 0.95]; 7 directional, 0 significant reversals among 120 pairs | `results/cleanx/E5_pilot/` ✅ (pilot) |
| 10 | Leakage / label noise in New-392 | 5 proteins identical to a training sequence (2 correct; B0XPE4 wrong despite an identical training copy with the right EC; 2 Swiss-Prot annotation inconsistencies, P9WNP2/Q7U2S5 vs P9WNP3); 43 (11%) ≥ 90% identity to training | `E1_exact_duplicates.csv` ✅ |
| 11 | Label drift (complete ECs) | New-392: 385 unchanged, 3 changed, 2 gained an EC, 1 lost one, 1 no longer an enzyme entry. ENZYME lists 1,314 transferred and 178 deleted ECs | `E1_label_drift.csv` ✅ |
| 12 | **Temporal benchmark** | 8,920 Swiss-Prot entries created ≥ 2022-03-01 → 2,099 enzymes (839 CAL / 1,260 TEST) + 2,000 non-enzymes; **640 enzymes (30.5%) carry ≥ 1 EC absent from CLEAN's vocabulary** | `sp_temporal_manifest.json` ✅ |
| 13 | **What "unseen" means** | TEST: 331 pure (no seen EC) + 53 mixed (also a seen EC, e.g. LPMOs annotated 1.14.99.54; 1.14.99.56); 209 carry an EC numbered above anything CLEAN knows in its sub-subclass (likely created after 2022), 173 an existing EC without a training member; largest single EC 1.14.99.56 (3.6% of TEST) | `E5_*/E5_benchmark_card.csv` ✅ |
| 14 | **CAL–TEST family overlap** | BLASTp: 345/1,260 TEST proteins (27%) ≥ 50% identical to a CAL protein, 59 ≥ 90% (fallback aligner: 355, 59) | `identity/blast/`, `identity/sp_temporal_test_vs_sp_temporal_cal.csv` ✅ |
| 15 | **Model vintage matters** | EnzPlacer's released training set contains 631/1,260 TEST (50.1%) and 430/839 CAL accessions, and 316/392 New-392 sequences; newest overlapping entry created 2025-04-09; 363 TEST proteins (261 seen-EC) were created ≥ 2025-06-01 | computed from the EnzPlacer repo ✅ |
| 16 | Identity profile (BLASTp; identity over the shorter sequence) | TEST: 22.6% < 30% identity to training, 40.9% ≥ 50%, 7.4% ≥ 90% (1 protein without a hit). CAL: 25.4% / 40.9% / 9.1%. Price-149: 17.4% / 32.2% / 1.3%. New-392: 15.1% / 54.8% / 11.0%. The fallback aligner missed hits and put 36–38% of TEST and Price-149 below 30% | `identity/blast/` ✅ |
| 17 | **Exemplar trade-off** (dev + confirmatory) | Exemplar − centroid micro-F1 (observed): New-392 +0.035 [+0.002, +0.068]; Price-149 −0.013 [−0.048, +0.024]. Price-149 by tier: common +0.062 (raw p 0.016), rare −0.085 (0.034), very rare −0.132 (< 0.001), singleton −0.043 (0.024); **Holm-adjusted: 0.048, 0.048, 0.002, 0.048** (borderline at B = 2000); the very-rare tier is 18 proteins but only 7 ECs | `E2_deltas.csv` ✅ (Holm added 2026-09-28) |
| 18 | Local-mean (3-NN) inference | Union-macro F1 New-392 **+0.036 observed** (bootstrap mean +0.032) [+0.007, +0.056]; Price-149 +0.006 (ns). K-means sub-prototypes slightly *worse* (negative result) | E2 ✅ |
| 19 | **Homology baseline (BLASTp, standard settings)** | EC of the best BLASTp hit against CLEAN's own training set **beats CLEAN** on New-392 (micro-F1 0.584 vs 0.498; CLEAN's own reported metric 0.574 vs 0.497) and on SP-temporal TEST (seen-EC micro-F1 0.735 vs 0.641, above every CLEAN-space rule, best 0.703; all of TEST, reported metric 0.495 vs 0.443); CLEAN wins on Price-149 (micro 0.426 vs 0.375; reported 0.495 vs 0.392). BLASTp on TEST seen-EC by identity: 0.46 (< 30%, n = 166), 0.69 (30–50%), 0.82 (50–70%), 0.87 (70–90%), 0.94 (≥ 90%). The old fallback aligner understated homology (New-392 0.505, TEST 0.600). Consistent with EC-Bench, where blastp (weighted F1 0.502 at EC4) beats its CLEAN retraining (0.251) [code]. Paired CIs and CLEAN per identity bin: ⟦final run⟧ | scratch analysis on the committed tables ✅; E2 ⟦final⟧ |
| 20 | **NISE-like ECs are dispersed in CLEAN space** | 183 of 2,493 ECs (7.3%) have members in ≥ 2 non-homologous Pfam-clan groups; size-stratified permutation: RMS radius +0.84, p < 10⁻⁴. The margin (−1.00) is *not* independent evidence (it is d_other − d_own); d_own and d_other are now tested separately (⏳ re-run E3) | `E3_stats.json` ✅ |
| 21 | **Tier-dependent overconfidence** (*pilot*, New-392 halves) | CLEAN GMM confidence by predicted tier: common ECE 0.18 (conf 0.64, acc 0.64); rare 0.36 (0.74 vs 0.48); very rare 0.50 (0.84 vs 0.34); singleton 0.49 (0.85 vs 0.36) | `E4_smoke/` ✅ (pilot) |
| 22 | **Certification budget** (*pilot*, New-392 half as calibration) | To certify selective risk 0.2 (power 0.8): common tier needs ≥ 630 calibration proteins (106 available; deficit 524 [34, 6,744] ≈ 2 years of Swiss-Prot accrual); rare, very-rare and singleton tiers have no operating point with estimated risk below α even at α = 0.3 (bootstrap resamples certifiable: 0–9% at α = 0.2, 0–29% at α = 0.3). SGR's own search needs ≈ 2× the oracle budget (simulation) | `E6_pilot/` ✅ (pilot) |
| 23 | Swiss-Prot accrual | 466 new enzymes/year pass the benchmark filters (2022-03 → 2026-09) | E6 ✅ |

**What we could *not* verify and will not claim:** other methods' published numbers re-computed from their own
releases (only EC-Bench ships per-protein predictions; the EC-Bench models are EC-Bench's retrainings, not the
original authors' checkpoints); full texts of most 2025–2026 papers (only abstracts/search results and code).

## 2b. Confirmatory results (run of the pre-registration commit c8e42f7 on the analysis PC, 2026-09-28)

Gates passed: smoke 10/10 (CX00 9/9); ESM-1b equivalence max |diff| 1.76e-07, min cosine 1.0000000; EC-Bench
reproduction gate 32/32 exact; Props. 1–2 / Lemma 3 hold on 82/82 prediction tables. The EC-Bench part of E5,
re-run in a Linux copy, gives identical numbers (cross-platform reproducibility). Numbers marked ⟦final⟧ are
re-read from the CSVs of the final run (BLASTp baseline, §7); no pre-registered endpoint depends on them.

| Endpoint (pre-registered) | Rule | Result | Verdict |
|---|---|---|---|
| **E2 primary (H3b)** | `balanced3` − `centroid`, union-macro F1 (EC4), SP-temporal TEST seen-EC (n = 876), paired bootstrap B = 2000, two-sided | **+0.043 [+0.030, +0.060]**, p < 0.001 (0.432 → 0.475) | **met** → G-A: C3 is a method contribution |
| E2 key secondary | micro-F1, same comparison | +0.047 [+0.032, +0.063] (0.641 → 0.689) | met |
| E2 robustness (split70 backbone) | same | +0.036 [+0.022, +0.052] | consistent |
| **E5 primary (H5)** | ≥ 1 of 120 EC-Bench pairs (test, 100%, EC4) ranked one way by the reported F1 and significantly the other way by micro-F1 (BH q < 0.05, sign agreement) | **1** significant (6 directional): DeepECtransformer-learnt ranks above PRIAM by the reported F1 (+0.016) but below by micro-F1 (−0.041 [−0.075, −0.006], q = 0.026). Mechanism: DeepECtransformer predicts ≈1.1 ECs/protein and 393 of its 452 FPs (87%) are invisible to the reported metric; PRIAM abstains on 402/468 proteins with micro precision 0.64. Kendall τ_b 0.90 [0.83, 0.93] | **met, by one pair** → G-D: "the blind spot changes a published method comparison"; rankings otherwise largely preserved |
| E5 key secondary (S4) | ρ(inflation, ECs/protein) > 0, one-sided | ρ = 0.75, p = 0.0004 | met |
| E5 S6 (complete-EC-only variant) | same count | 1 significant (DeepECtransformer-regular vs PRIAM, q = 0.015) | consistent |
| **E6 primary (H6)** | α = 0.2: budget > CAL proteins available in the very-rare and singleton predicted tiers | both tiers: estimated selective risk ≥ 0.2 at every usable coverage → not certifiable at any budget (46% and 21% of bootstrap resamples could certify); deficit 95% CI [528, ∞) and [515, ∞) proteins (102 and 90 available), i.e. ≥ 9.3 and ≥ 10.3 years of Swiss-Prot accrual at the lower bound. Rare tier: deficit 49,218 (unstable point estimate), CI [428, ∞), ≥ 6.3 years. Common tier: feasible now (deficit 110, CI [−250, 2,055]) | **met**; **G-E met** (deficit CI above 0 in all three rarer tiers) |
| H4a | GMM ECE rises common → singleton | 0.19 / 0.22 / 0.32 / 0.32 (confidence − accuracy −0.02 / +0.15 / +0.30 / +0.28) | met |
| H4b | per-tier (Mondrian) risk ≤ α on TEST; a global threshold not | global thresholds exceed α = 0.2 in the singleton tier: 4/14 = 0.29 [0.08, 0.58], 5/15 = 0.33 [0.12, 0.62], 11/24 = 0.46 [0.26, 0.67] errors (exact 95% CI) for centroid, exemplar, balanced3; only balanced3's CI excludes 0.2. Overall, each global threshold holds (0.12–0.14). Mondrian never exceeds α but accepts nothing in the rare tiers, and nothing at all for centroid at α = 0.2 | validity met, usefulness not → G-B: report the smallest certifiable risk per tier (centroid 0.25 common vs 0.60–0.65 rarer tiers) |
| H3c / G-C | NISE exemplar gain > other | 0.073 (n = 87) vs 0.062 (n = 787), p = 0.36 | not met → Supplement; NISE ECs are dispersed (d_own AUROC 0.82, p < 10⁻⁴), not crowded (d_other p = 0.36) |
| S7 | invisible share of CLEAN's FPs on TEST < New-392's 0.67 | 0.74 (616/830) | **not met** (the prediction was wrong: label-set size alone does not set the invisible share) |
| H1 (second half) | reported precision rises under over-prediction | only on Price-149 (0.52 → 0.64, top-1 → top-10); it falls on New-392 (0.56 → 0.48) and TEST (0.45 → 0.28) while micro precision collapses (→ 0.09, 0.07, 0.07) | partly met; exactly Prop. 4's condition (a small label set) |
| H3a | exemplar helps common, hurts rare tiers | exemplar − centroid, micro-F1 by true tier on TEST (Holm across 4 tiers): common +0.082 [+0.060, +0.106], p_Holm < 0.005; rare +0.025 [−0.011, +0.063], very rare +0.048 [−0.020, +0.109], singleton −0.020 [−0.068, +0.026] (p_Holm 0.48 each). `balanced3`: common +0.070 [+0.051, +0.089], p_Holm < 0.005; rarer tiers +0.002 / +0.016 / −0.022, ns. Price-149 shows the full pattern (common +0.063, rare −0.083, very rare −0.130, singleton −0.042; p_Holm ≤ 0.048), but that pattern was seen during development (§7.3), so it is not confirmatory | **half met**: helps common; no significant harm in rarer tiers on TEST. The member-level gain comes from common ECs |

**Homology baseline (BLASTp, after the run; §7):** the best BLASTp hit against CLEAN's training set beats CLEAN on
New-392 and on TEST (seen-EC micro-F1 0.735 vs 0.641; above every CLEAN-space rule), and loses on Price-149 (0.375
vs 0.426); details in §2 row 19. This does not touch a pre-registered endpoint, but it changes what C3 can claim.

**Descriptive, not pre-registered (reported as such):** on TEST every member-level geometry beats the centroid
(union-macro: exemplar +0.079, kmean3 +0.063, balanced1 +0.055, balanced3 +0.043; multiproto +0.007, ns), and
the class-size correction *lowers* union-macro F1 relative to uncorrected kmean3 (0.475 vs 0.494). The gain
therefore comes from local, member-level inference, not from correcting a class-size bias; C3 is framed that way.
On Price-149 no geometry differs from the centroid (balanced3 +0.008 [−0.016, +0.027]).

---

## 3. Novelty and gap analysis (literature and code checked 2026-09-28)

Conclusion per contribution: **no prior report found / partially reported / already reported**, with evidence.

| Contribution | Closest prior work (evidence level) | Verdict |
|---|---|---|
| **C1a** invariance of CLEAN-style / support-weighted EC scores to out-of-label predictions, with proof | No paper found stating or quantifying it for EC prediction. *The Accuracy Paradox* (arXiv:2609.07897, Sept 2026) [search] critiques default 0.5 thresholds for 6 top-level EC classes on compounds: a different mechanism. sklearn documents zero-support handling, not its consequence | **no prior report found** (the mathematics is elementary; the contribution is the consequence and its prevalence) |
| **C1b** field-wide prevalence | [code] CLEAN `app/src/CLEAN/evaluate.py` (`get_eval_metrics`, binarizer fitted on test labels, `average='weighted'`); CLEAN-Contact `src/CLEAN/evaluate.py:320-342` @0e423db; ProtDETR `util/clean.py:311-331`, `engine.py:217-218` @aa49809; EnzymeHunter `src/evaluate.py:241-262` @fd2aa18; EC-Bench `code/results.ipynb` cell 60 (weighted over union) and `CLEAN/app/src/CLEAN/evaluate.py:301-322` (micro over test labels) @5e98a2c | **no prior report found** |
| **C1c** gameability / inflation vs over-prediction | E1 (ours) + E5 pilot; EC-Bench notes "additional EC numbers predicted" as a diagnostic [code, notebook cell 99] but does not relate it to the metric | **no prior report found** |
| **C2** temporal hold-out for the *released* CLEAN with unseen-EC share, ENZYME transfer mapping and non-enzymes | EC-Bench (temporal split with retraining) [search]; ECPred temporal hold-out (2018) [search]; SOLVE post-2023 set [search]; EnzPlacer unseen-EC4 experimental set [code: data/] | **partially reported** (resource contribution; the vintage-overlap finding is new) |
| **C3** controlled inference-geometry comparison on one frozen contrastive space | kNN head-class bias and prototype classifiers for long tails are known in ML (Learning Prototype Classifiers for Long-Tailed Recognition, IJCAI 2023, doi:10.24963/ijcai.2023/151; logit adjustment) [search]; FEDKEA and EnzPlacer use kNN transfer [search/code] | **partially reported** (mechanism known in general ML; the controlled EC analysis and tier trade-off are new) |
| **C4a** tier-dependent miscalibration of CLEAN's released confidence | Boger et al. 2025 (Nat. Commun.) apply conformal risk control to CLEAN's embedding distances (max-hierarchical loss, marginal guarantee) [code: docs/VERIFICATION_NOTES.md, "CLEAN Enzyme Classification (Tables 1-2)"]; CPEC (FDR, PenLight2) [code]; ECPICK/HIT-EC "trustworthy" classifiers [search] | **partially reported** — conformal guarantees on CLEAN exist; per-tier miscalibration of the released GMM confidence is not reported |
| **C4b** per-tier certification budget (exact binomial) and curation-time conversion | Small-sample conformal (SSBC, arXiv:2509.15349) and class-conditional conformal with many classes (arXiv:2306.09335) [search] give general theory; nothing EC-specific found | **no prior report found** for EC prediction; general theory exists |
| Chemistry-weighted errors (earlier idea) | EC-Bench scores partial/additional ECs with RDKit reaction similarity [search] | **already reported** — dropped |
| New training losses | EnzPlacer, EnzymeHunter, GloEC | **already reported** — not pursued |
| PLM > BLAST at low identity; BLAST strong at high identity | *PLMs outperform BLAST for evolutionarily distant enzymes* (bioRxiv 10.64898/2026.03.31.715487) [search]; EC-Bench (blastp above its CLEAN retraining) [code]; CLEAN's paper reports CLEAN ahead of BLASTp on New-392 (⟦read its BLASTp configuration in the Methods before writing⟧) | **partially reported** — new here: with standard settings BLASTp beats the *released* CLEAN on New-392 and on a temporal hold-out, and the identity-stratified crossover |
| Non-functional decoys | EnzymARC (bioRxiv 10.64898/2026.08.21.746242): CLEAN > 90% false-positive rate on active-site-disrupted decoys [search] | **already reported** — cite for the open-set motivation |

Re-run before submission (the field moves monthly): `CLEAN enzyme weighted F1 evaluation`,
`enzyme commission benchmark 2026`, `conformal selective enzyme function rare`, `EnzPlacer`, `EC-Bench`.

---

## 4. Research questions, hypotheses, endpoints

| RQ | Hypothesis | Endpoint (pre-registered) | Status |
|---|---|---|---|
| RQ1 How much of CLEAN's error does the reported metric see, and can it be gamed? | **H1** Support-weighted P/F1 ignore most FPs and rise under over-prediction; micro, union-macro and per-protein F1 do not | Share of FPs outside the test label space; score vs number of ECs predicted (1–10) | first half ✅ everywhere; the rise only on Price-149; S7 not met (§2b) |
| RQ2 How does CLEAN generalise past 2022? | **H2** Exact-EC performance falls steeply with rarity tier and identity; a large share of new enzymes carry unseen ECs | Micro/union-macro F1 by depth × tier × identity bin; unseen-EC share | unseen share ✅; strata ⟦final run, BLASTp identity⟧ |
| RQ3 Is one centroid per EC the right way to use CLEAN's embedding? | **H3a** exemplar helps common, hurts rare tiers. **H3b (primary)** `balanced3` beats centroids on union-macro F1. **H3c** NISE-like ECs are more dispersed (✅) and exemplar gains concentrate in them | **Primary:** union-macro F1 (EC4), `balanced3` − `centroid`, SP-temporal TEST, seen-EC proteins, paired bootstrap (B = 2000), two-sided α = 0.05. **Key secondary:** micro-F1. Sensitivity (2026-09-28): EC-cluster bootstrap, ≤ 5 proteins per EC, CAL-disjoint TEST, protein-level evidence | **H3b met** (+0.043 [+0.030, +0.060]); H3a half met (helps common; no significant harm in rarer tiers); H3c not met |
| RQ4 Is CLEAN's confidence reliable, and can abstention be certified per tier? | **H4a** GMM ECE rises from common to singleton. **H4b** Mondrian SGR keeps per-tier risk ≤ α on TEST; a global threshold does not. **H4c** smallest certifiable risk differs by tier | ECE by predicted tier; realized selective risk/coverage per tier at α ∈ {0.1, 0.2, 0.3}; min certifiable risk; abstention on pure vs mixed unseen-EC enzymes and non-enzymes | H4a met; H4b valid but vacuous in rare tiers; H4c met |
| **RQ5** Does the reported metric change conclusions across methods? (E5) | **H5 (primary)** ≥ 1 of 120 method pairs on EC-Bench's test set (100% regime, EC4) is ordered one way by the reported weighted F1 and *significantly* the other way by micro-F1 | Number of such pairs; paired bootstrap B = 2000, BH q < 0.05, sign agreement. **Key secondary (amended after pilot):** Spearman ρ(inflation, ECs/protein) > 0, one-sided α = 0.05 | **H5 met (1 pair)**; S4 met (ρ 0.75) |
| **RQ6** How much calibration data would per-tier certification need? (E6) | **H6 (primary)** For α = 0.2 the budget exceeds the CAL proteins available in the very-rare and singleton predicted tiers | Calibration deficit per tier with bootstrap 95% CI (1000 resamples); secondaries: α 0.1/0.3, other geometries, years of accrual, SGR simulation, realized TEST certification | **H6 met**; G-E met (deficit CI lower bounds 428–528 proteins, ≥ 6.3–10.3 years) |

### Decision gates

- **G-A (E2):** primary CI excludes 0 in favour of `balanced3` → C3 is a method contribution. Otherwise C3 is a
  controlled negative result plus mechanism; emphasis on C1/C2/C4.
- **G-B (E4):** Mondrian thresholds hold per-tier risk where a global threshold fails → positive guarantee
  result. Otherwise report the smallest certifiable risk per tier.
- **G-C (E3):** significant NISE gain difference → mechanism figure; otherwise Supplement. The d_own / d_other
  split decides whether NISE ECs are "dispersed", "crowded", or both.
- **G-D (E5):** ≥ 1 significant reversal → "the blind spot changes method comparisons in a published
  benchmark". Otherwise, if S4 holds → "inflation tracks over-prediction; rankings largely preserved". If both
  are null → C1 rests on CLEAN's benchmarks, the proof and the code audit (still publishable).
- **G-E (E6):** deficit CI above 0 in the rare tiers → C4 reports "rare-EC calls cannot be certified with the
  calibration data that exist: N proteins ≈ Y years of curation; use tier pooling or hierarchical abstention".
  Otherwise → certification is feasible now (G-B expected positive).

---

## 5. Datasets

| Set | Size | Role | Notes |
|---|---|---|---|
| split100 (CLEAN training) | 227,362 sequences, 5,242 ECs, 241,025 (EC, seq) pairs | index for all inference geometries (released `100.pt`) | UniProt 2022_01 snapshot |
| split70 | 79,172 pairs, 5,093 ECs | second backbone | released `70.pt` + `split70.pth` |
| New-392 | 392 | **development** (disclosed) | 5 exact training duplicates; all ECs seen |
| Price-149 | 149 | confirmatory for the frozen centroid/exemplar/kmean3/multiproto; **post hoc** for `balanced*` | RefSeq IDs; low identity |
| SP-temporal CAL | 839 enzymes | λ tuning (E2), calibration (E4), budgets (E6) **only** | stratified random 40% by tier |
| SP-temporal TEST | 1,260 enzymes (876 seen-EC, 384 unseen-EC: 331 pure, 53 mixed) | **confirmatory** | untouched by any model selection; 28% ≥ 50% identical to a CAL protein |
| SP non-enzymes | 2,000 (of 4,058 eligible) | open-set rejection (descriptive; no guarantee) | not split into CAL/TEST |
| EC-Bench predictions | 16 methods × {test 468, Price-149 149} × {100%, 30%} | E5 (test = confirmatory; Price-149 = pilot) | MIT licence, commit 5e98a2c, SHA-256 in manifests |
| NISE members (optional) | 3,817 training sequences of 168 NISE-like ECs | E3b | FASTA written |

**Temporal construction rules** (`cx_02_build_temporal.py`): Swiss-Prot, created ≥ 2022-03-01; not a fragment;
length 30–1022; accession not in any CLEAN split, New-392, Price-149 or halogenase; sequence not identical to any
CLEAN training sequence; within-set duplicates removed; ≥ 1 complete EC. Labels are current complete ECs;
CLEAN's predictions are mapped through the ENZYME transfer table (`labels.harmonize`). Tier = count of the
protein's rarest true EC in split100 after transfer mapping; `unseen` = 0.

**Known limitations to state:** ENZYME records no transfer for some reclassifications (e.g. 3.6.4.12 →
5.6.2.3/5.6.2.4); such proteins are analysed as open-set. One family dominates the unseen tier (fungal LPMOs,
78 proteins with 1.14.99.56, 73 created in one UniProt release) → per-EC cap and EC-cluster bootstrap
sensitivity analyses. CAL and TEST are a random split of one pool → CAL-disjoint TEST sensitivity analysis.

---

## 6. Experiments

| ID | Question | Script | Status | Runtime (your PC) | Key outputs |
|---|---|---|---|---|---|
| E0 | Exact reproduction | `cx_00_check.py` | ✅ (your PC) | 15 s | 9/9 PASS |
| E1 | Evaluation validity (RQ1) | `cx_01_audit.py` | ✅ legacy + TEST (confirmatory) | 10 s | blind spot, over-prediction, tier×depth, duplicates, drift |
| D | Temporal benchmark | `cx_02_build_temporal.py --with-pfam` | ✅ | 20 s | sets, FASTA, manifest |
| B | ESM-1b embeddings | `cx_03_embed.py` | ✅ (your PC; equivalence PASS, 0 missing) | 1.5 min now | `data/esm_data/*.pt` |
| BL | BLASTp search | `cx_11_blast.py` | ✅ (Linux copy, BLAST+ 2.12.0) | 21 min | hit tables + SHA-256 manifest |
| I | Identity + homology baseline | `cx_04_identity.py --tool blast` | ✅ tested; ⏳ final run | seconds | identity bins, top-hit predictions |
| I2 | TEST vs CAL identity | `cx_04_identity.py --tool blast --pair sp_temporal_test:sp_temporal_cal` | ✅ tested; ⏳ final run | seconds | CAL-disjoint subset |
| E2 | Inference geometry (RQ3) | `cx_05_geometry.py` (+ `--backbone split70`) | ✅ confirmatory; ⏳ final run (only homology, identity and CAL-disjoint rows change) | ~5 min | metrics, deltas (observed + bootstrap), Holm, sensitivity subsets |
| E3 | NISE mechanism | `cx_06_ec_biology.py` | ✅ confirmatory | 30 s | stats JSON |
| E4 | Reliability (RQ4) | `cx_07_reliability.py` | ✅ confirmatory | 50 s | calibration, risk control (with exact CIs), open set (pure/mixed), fallback |
| **E5** | **Metric validity across methods (RQ5)** | `cx_09_evalcard.py` | ✅ confirmatory (re-run identical in a Linux copy) | 45 s | cards, pairs, rank agreement, invariance, benchmark card |
| **E6** | **Certification budget (RQ6)** | `cx_10_certbudget.py` | ✅ confirmatory (re-run identical in a Linux copy) | 3 min | n* grid, tier budgets with CIs, accrual, SGR simulation |
| F | Figures | `cx_08_figures.py` | ✅ revised and rendered on the confirmatory results; ⏳ final run | 30 s | PNG + PDF + CSV twins, CIs |
| T | Unit + integration tests | `tests/run_tests.py` | ✅ 29/29 | ~1 min | PASS/FAIL per test |

**Verification (2026-09-28, Linux replica with your exact pins:** Python 3.10.20, numpy 1.22.3, pandas 1.4.2,
scipy 1.7.3, scikit-learn 1.2.0, torch 1.11.0, fair-esm 2.0.0, biopython 1.83, matplotlib 3.7.0): 29/29 tests (28
known-answer unit tests + an end-to-end run of cx_04/05/06/07/08/09/10/11 on a synthetic app in a temp folder);
`cx_10 --smoke` and `cx_11 --smoke` PASS; the EC-Bench part of E5 and all of E6 re-run from your confirmatory
features give your numbers exactly (apart from the NaN → ∞ reporting fix in E6); the revised figures were rendered on your confirmatory results (every F1
reproduces `E2_metrics_long.csv`). **Not run there:** anything that loads CLEAN's weights (cx_00, E1–E4 on real
data), which stay on your PC (licence).

**Statistics:** paired non-parametric bootstrap over proteins (B = 2,000; seed 20260927); effects reported as the
**observed** full-sample difference with the 95% percentile CI; two-sided bootstrap p; per-tier tests Holm-
corrected across tiers; E5 pairs BH-corrected across 120 pairs; E3 size-stratified permutation (10,000); SGR
with δ = 0.1, Bonferroni-split over 4 tiers; E6 exact binomial (Clopper–Pearson duality).

**Metrics reported everywhere:** micro P/R/F1; union-macro F1; per-protein F1; subset accuracy; and, for
comparability, the support-weighted and true-label macro scores **explicitly marked as blind to out-of-label-set
false positives** (Prop. 1). Every table carries the invisible-FP count (the "EC evaluation card").

---

## 7. Pre-registration and anti-leakage protocol

1. **Freeze before the first confirmatory run.** Done: the Methods say "analysis code frozen at commit c8e42f7
   (c8e42f704a406104a6fb61868392e73991f25698), 2026-09-28, before the first confirmatory run; its results are
   commit 49a6747; the changes after that run (§7.7) are in later commits". Confirmatory = anything on
   SP-temporal TEST (E1–E6) and E5 on the EC-Bench test set.
2. **CAL only for fitting:** λ (E2), every calibrator/threshold (E4) and every selective-risk estimate for the
   budgets (E6) come from SP-temporal CAL. TEST is only scored.
3. **Disclosure:** New-392 was used for piloting (1-NN vs centroid, dispersion scaling (rejected), BIC
   multi-prototypes (rejected), GMM calibration, E6 pilot). Price-149's per-tier result motivated `balanced*`
   (post hoc on Price-149). EC-Bench's Price-149 predictions were the E5 pilot.
4. **Leakage controls:** accession + exact-sequence exclusion; identity strata (BLASTp); CAL-disjoint TEST;
   protein-level-evidence subset; per-EC cap; EC-cluster bootstrap.
5. **No tuning on anything later reported as a test result.**
6. **Dated record of changes, all before any SP-temporal TEST result existed:**
   - *2026-09-27:* E2 primary metric and λ criterion switched from true-label to union-label macro-F1 after E1
     showed the former rises under over-prediction.
   - *2026-09-28:* (a) E2 reports observed deltas and Holm-adjusted per-tier p (the bootstrap mean of a
     macro-averaged metric is biased: New-392 kmean3 union-macro +0.036 observed vs +0.032 bootstrap mean);
     sensitivity analyses added (EC-cluster bootstrap, ≤ 5 proteins per rarest true EC, CAL-disjoint TEST at
     < 50% identity, protein-level evidence). (b) E3 tests d_own and d_other separately (the margin is their
     difference). (c) E4 splits unseen-EC enzymes into pure and mixed, and certifies the EC3 fallback on the
     calibration proteins the EC4 rule rejects. (d) **E5 pre-registered** (docstring of `cx_09_evalcard.py`)
     before any EC-Bench metric other than the reported one was computed; **amendment** after the Price-149
     pilot, before confirmatory scoring: S4 promoted to key secondary with a one-sided direction; primary
     unchanged. (e) **E6 pre-registered** (docstring of `cx_10_certbudget.py`) before any CAL feature existed;
     the Jeffreys shrinkage and the ≥ 10-accepted-protein rule were set after the disclosed pilot.
7. **Changes after the confirmatory run (2026-09-28; no pre-registered endpoint or decision rule changed; both runs
   are in the repository history):**
   - *Homology baseline and identity:* BLASTp (BLAST+ 2.12.0, default scoring, E ≤ 10, best hit by bit score)
     replaces the fallback aligner (3-mer prefilter + Smith–Waterman on 30 candidates), for the baseline, the
     identity strata and the CAL-disjoint subset. The switch to a standard aligner was planned before the run
     (RUNBOOK §4 at c8e42f7: "recommended for final paper numbers"); BLASTp rather than MMseqs2 because it is the
     baseline of CLEAN's paper and of EC-Bench. Hit tables are stored with SHA-256 digests of the sequences searched.
   - *Bug fix:* the EC-cluster bootstrap divided per-protein F1 by N instead of each resample's total weight;
     only the `samples_F1` rows of that sensitivity analysis were affected (micro and union-macro F1 are ratios).
   - *Reporting:* realized selective risks carry error counts and exact 95% CIs (E4); E6 bootstrap bounds between
     two non-certifiable resamples are +∞ (numpy's interpolation returned NaN), and E6 adds the years of accrual
     at the deficit's lower bound (`years_at_deficit_ci_lo`), because the point estimate is unstable near the
     certifiability boundary (rare tier: 49,218 proteins, "727 years", from a risk estimate of 0.18 on 2 errors in 13 accepted proteins).
   - *Figures (reporting only):* 95% CIs on every F1 and risk; reliability diagrams compare the two confidences
     on the same CLEAN predictions (before: CLEAN GMM on centroid vs recalibrated exemplar, two changes at once);
     risk–coverage adds CLEAN with tier-conditional confidence; the per-tier risk figure shows all three
     geometries; the EC-Bench figure marks the significant reversal; the budget figure shows the deficit CI and
     the years at its lower bound instead of the unstable point estimate; new forest plot of the primary endpoint.
   - *Interpretation:* claims that the confirmatory data did not support are withdrawn or narrowed (§2b: S7,
     over-prediction reward, the class-size mechanism of C3); the descriptive geometry comparisons are labelled as
     such.

---

## 8. Figures and tables

| Item | Content | Produced by |
|---|---|---|
| Fig. 1 | Schematic / graphical abstract: frozen CLEAN encoder → geometries → calibration/abstention → budget | draw manually |
| Fig. 2 | Blind spot: reported vs micro precision; visible vs invisible FPs (+ over-prediction curve) | `fig1_blind_spot` |
| Fig. 3 | **EC-Bench re-scoring: reported vs micro F1 for 16 methods; inflation vs ECs/protein** | `fig8_ecbench_rescoring` (pilot: `figS2_…_pilot`) |
| Fig. 4 | **E2 primary endpoint (forest plot):** `balanced3` − centroid with the sensitivity analyses, the descriptive geometries, BLASTp and the other test sets; union-macro and micro F1 with 95% CIs | `fig10_primary_endpoint` |
| Fig. 5 | Exact-EC micro-F1 by true tier with 95% CIs: centroid / `balanced3` / exemplar / BLASTp top hit; * = Holm-significant pre-registered per-tier test | `fig2_tiers` |
| Fig. 6 | Micro-F1 vs BLASTp identity bin with 95% CIs (bins with n < 10 marked) | `fig3_identity` |
| Fig. 7 | Reliability diagrams by predicted tier on the same CLEAN predictions: released GMM confidence vs tier-conditional isotonic (fitted on CAL); exact CIs per bin; ECE | `fig5_reliability_by_tier` |
| Fig. 8 | (a) Risk–coverage: centroid with GMM vs tier-conditional confidence, exemplar and `balanced3` tier-conditional (separates the calibration and geometry effects); (b) realized per-tier selective risk with exact CIs, global vs per-tier thresholds, all three geometries | `fig6_risk_coverage`, `fig7_tier_risk_control` |
| Fig. 9 | **Certification budget: n* vs true selective risk; needed (point and 95% CI) vs available per tier, years of accrual at the lower bound** | `fig9_certification_budget` (pilot: `figS3_…_pilot`) |
| Supplement | NISE gain (G-C not met), with Mann–Whitney p | `fig4_nise_gain` |
| Table 1 | Datasets (sizes, tiers, identity bins, unseen pure/mixed, EC concentration) | manifest, `E5_benchmark_card.csv` |
| Table 2 | Code audit: 9 code bases, averaging, label space, blind or not (file:line@commit) | §3 |
| Table 3 | Main results with CIs (methods × datasets × metric families, evaluation cards) | `E2_bootstrap_ci.csv`, `E5_cards.csv` |
| Table 4 | Calibration (ECE) and certification budget by tier | `E4_calibration.csv`, `E6_tier_budget.csv` |
| Supplement | full tables; duplicates; drift; split70; MMseqs2 vs fallback identity; λ curve; depth ablation (`figS1`); E3b; sensitivity subsets; E5 levels 1–3, 30% regime, complete-EC variant; E5/E6 pilots | E-scripts |

Figure style fixed in `cx_08_figures.py`: validated colour-blind-safe palette (re-validated 2026-09-28 with the
dataviz validator: CVD ΔE 9.2, normal-vision ΔE 27.6; the green needs visible labels, which every figure has; the
violet added for `balanced3` passes the all-pairs check next to blue, orange and green, worst CVD ΔE 9.2),
fixed method → colour mapping (grey for anything that is not a method), one y-axis per panel, direct labels, CSV
twin per figure. Every interval is a 95% interval (paired bootstrap B = 2000 as in E2, or exact Clopper–Pearson);
F1 point estimates are recomputed from `E2_per_query.csv` and a figure is skipped if they differ from
`E2_metrics_long.csv`. Captions under the plots are wrapped so that every figure is at most 7.5 in wide (double
column).

---

## 9. Manuscript blueprint

### 9.1 Title options
1. *What reported scores hide: metric blind spots, inference geometry and the calibration budget for trusting rare enzyme-function calls*
2. *Weighted metrics, single centroids and silent overconfidence: auditing contrastive enzyme function prediction*
3. *Beyond the centroid: class-balanced retrieval and certified abstention for contrastive EC-number prediction* (only if G-A and G-B are both positive)

### 9.2 Structure and word budget (~6,500 words main text)
- **Introduction (~800):** EC prediction; CLEAN and successors (CLEAN-Contact, EnzPlacer, EnzymeHunter, ProtDETR,
  GloEC, HIT-EC, MAPred); benchmarks (CARE, EC-Bench); calibration (CPEC, Boger et al. 2025). Gap (§3).
  Contributions C1–C4.
- **Methods (~2,100):** reproduction; formal definitions and **Props. 1–4 with proofs** (≈ half a page; text in
  `cleanx/evalcard.py`); code audit protocol; EC-Bench re-scoring (gate, pairs, BH); temporal benchmark;
  identity; the geometries and the balanced-density derivation; NISE; calibration, Mondrian SGR; **certification
  budget** (duality, power, SGR δ', fixed point, Jeffreys estimate, lower-bound argument); statistics; versions.
- **Results (~2,500):** 3.1 reproduction + blind spot + gaming; 3.2 the blind spot across the field (code audit,
  EC-Bench); 3.3 temporal generalisation (unseen decomposition, vintage); 3.4 inference geometry; 3.5 NISE;
  3.6 reliability and the certification budget.
- **Discussion (~900):** boxed **EC evaluation card** checklist (micro + union-macro + per-protein F1, FP count and
  invisible share, label-space size, tier and identity strata, unseen-EC share, model vintage); implications for
  CLEAN users and annotation pipelines; limitations (frozen model, one model family, EC-Bench retrainings,
  Swiss-Prot label noise, EC reclassification, small rare tiers, oracle budget); future work (EnzPlacer
  vintage-matched audit, balanced training objectives).
- **Conclusion (~150).**

### 9.3 Abstract draft (244 words by whitespace count with the two ⟦ ⟧ as written, counted by script on 2026-09-28; filled from the confirmatory run except the two ⟦ ⟧, which depend on the BLASTp run: read them from `E2_metrics_long.csv` and `fig3_identity.csv`)
> Contrastive models such as CLEAN assign Enzyme Commission (EC) numbers by the nearest class centroid and are
> compared by support-weighted precision and F1. We prove that these scores, as computed in the public code of
> CLEAN, four successors and the EC-Bench benchmark, cannot see predicted EC numbers absent from the test labels,
> and that support-weighted recall equals micro recall. On standard benchmarks the reported precision misses
> 67–96% of CLEAN's false positives, and predicting ten ECs per protein raises CLEAN's reported F1 on Price-149
> from 0.495 to 0.618. Across 16 methods re-scored from EC-Bench's published predictions, the inflation tracks how
> many ECs a method predicts (Spearman ρ = 0.75) and one published comparison reverses. On a leakage-controlled
> temporal benchmark of 2,099 Swiss-Prot enzymes created after CLEAN's training snapshot and 2,000 non-enzymes,
> 30% of new enzymes carry an EC number absent from CLEAN's vocabulary. In CLEAN's frozen embedding space, a
> class-balanced nearest-member rule raises the pre-registered union-label macro-F1 by 0.043 (95% CI
> 0.030–0.060), mainly on common ECs, yet a standard BLASTp top hit on the same training set scores higher
> (micro-F1 ⟦0.735 vs 0.641⟧) ⟦except below 30% identity⟧. CLEAN's released confidence is most overconfident on
> rare EC numbers (expected calibration error 0.32 vs 0.19). Exact binomial certification budgets show that
> certifying rare-EC calls at 20% risk needs at least 550–630 calibration proteins per rarity tier (90–122 exist),
> at least 6–10 years of Swiss-Prot curation. We release an EC evaluation card, the benchmark and all analysis
> code.

### 9.4 Highlights (character counts include spaces, counted by script on 2026-09-28; all ≤ 85)
1. `Five public EC-prediction code bases use scores blind to out-of-label predictions` (81)
2. `Reported precision hides 67-96% of CLEAN's false positives on standard benchmarks` (81)
3. `A 2022-2026 Swiss-Prot hold-out: 30% of new enzymes carry ECs CLEAN never saw` (77)
4. `CLEAN's released confidence is most overconfident on rare EC numbers` (68) — confirmed (H4a)
5. `Certifying rare-EC calls needs more calibration data than Swiss-Prot adds yearly` (80) — confirmed (H6, G-E)

### 9.5 Keywords (≤ 6, confirm the limit)
Enzyme Commission number; enzyme function prediction; evaluation metrics; contrastive learning; calibration;
selective prediction.

---

## 10. Timeline (working days, with gates)

| Day | Work | Gate |
|---|---|---|
| 1 ✅ | E0, E1 legacy, benchmark, identity, E2/E3 legacy, pilots | done |
| 2 ✅ | Literature + code audit; E5 and E6 implemented, pre-registered, piloted; rigor fixes; 22 tests | done (commit pending) |
| 3 ✅ | Confirmatory run of c8e42f7 on your PC (smoke 10/10, full pipeline); audit; BLASTp baseline (§7) | **G-A met, G-B → min certifiable risk, G-C → Supplement, G-D met, G-E met** |
| 4 | **You:** commit the confirmatory results ✅ (49a6747); pull the audit; `.\cx_smoke.bat` (11/11); `.\cx_full.bat` (≈20 min; embeddings skipped); commit the final results | final numbers |
| 5 | Fill every ⟦final⟧ marker; Tables 1–4; Supplement; optional EnzPlacer vintage-matched audit (§11) | numbers frozen |
| 6–12 | Write Methods + Results, then Introduction/Discussion, abstract, highlights, graphical abstract | internal read-through |
| 13–15 | Licence steps (§13), Zenodo archive, final literature re-check, cover letter, submission | submit |

---

## 11. Risks and mitigations

| Risk | Likelihood | Mitigation |
|---|---|---|
| ESM-1b load fails (RAM) or is too slow | medium | close apps; resumable; or a GPU server (embedding is CPU-only in `esm_embed.py`: add a device option first) |
| E2 primary (H3b) negative | medium | pre-registered fallback (G-A) |
| E5 primary met by a single pair (q = 0.026) | certain | report it as exactly that, with its mechanism (abstaining vs over-predicting method) and the S6 replication; the key secondary (ρ = 0.75) carries the general claim |
| Rare tiers cannot be certified | high | this *is* E6's result: report the budget (G-E) |
| Reviewer: "the metric point is trivial" | medium | the proof is short by design; the contribution is prevalence (5 code bases), size (67–96%), gameability, and the evaluation card |
| EC-Bench models are retrainings, not original checkpoints | certain | state it; the claim is about the metric's response across 16 heterogeneous predictors, validated by exact reproduction |
| Boger et al. 2025 already put conformal guarantees on CLEAN | certain | cite prominently; ours is per tier, on the released confidence, with open-set rejection and the data budget |
| Model vintage leakage for any 2023+ model | real | never evaluate EnzPlacer on the full TEST; only the 363 proteins created ≥ 2025-06-01 (+ Price-149) |
| CAL–TEST family overlap (28% ≥ 50% identity; LPMO cluster) | real | CAL-disjoint TEST, EC-cluster bootstrap, per-EC cap (pre-registered sensitivity) |
| E6 budget is an oracle lower bound | certain | stated; makes "needed > available" conservative; SGR simulation reports the gap (~2×) |
| New overlapping paper appears | medium | re-run §3 searches before submission |
| **BLASTp beats CLEAN on New-392 and TEST** (reviewer: "why use CLEAN?") | certain | report it prominently as a finding; show the identity-stratified crossover (where CLEAN still helps: low identity, Price-149); document the BLAST configuration exactly (manifest) and contrast it with the BLASTp setup in CLEAN's paper; C3 becomes "how much of the gap member-level inference closes" |
| Identity/baseline via the fallback aligner | resolved | replaced by BLASTp (§7); the fallback understated homology |
| **Licence** | real | §13 |

**Optional next extension (design only, pre-register before running): EnzPlacer vintage-matched audit.** Run the
released `EnzPlacer_seen.pth` (GPL-3.0; reference embeddings on Zenodo 10.5281/zenodo.18110452) on Price-149 (no
overlap) and on the 363 TEST proteins created ≥ 2025-06-01 (261 seen-EC), reusing our ESM-1b embeddings (its
extraction matches CLEAN's: layer 33, mean over residues 1..L), and score it with the evaluation card
(`cx_09_evalcard.py --card`). Needs ~0.8 GB download and CPU only.

---

## 12. Reviewer questions and prepared answers

- **"Isn't Proposition 1 trivial?"** As mathematics, yes, and deliberately stated as such. The contribution is that
  five public code bases (including the benchmark the field uses) report scores it applies to, that the hidden
  share is 67–96% of false positives, that the inflation tracks over-prediction across 16 methods (ρ = 0.75,
  confirmatory) and reverses a published comparison, and a drop-in evaluation card that fixes the reporting.
- **"EC-Bench's models were retrained — are these the real tools?"** The question is how the metric responds across
  heterogeneous predictors; we reproduce EC-Bench's own published numbers exactly before re-scoring.
- **"Boger et al. already certified CLEAN."** They control a marginal, set-valued hierarchical loss. We show why
  marginal guarantees hide rare-tier failure (tier-dependent overconfidence) and quantify the per-tier
  calibration data that tier-level guarantees need; the two are complementary.
- **"Why not retrain CLEAN / evaluate a second model?"** Retraining needs ESM-1b for 227k sequences; our analyses
  apply to any embedding model via the evaluation card. EnzPlacer's training data contain half of our TEST set,
  so we specify a vintage-matched subset (§11).
- **"Is the certification budget real or an SGR artifact?"** It is an exact binomial lower bound for certifying a
  fixed coverage with a one-sided Clopper–Pearson bound; SGR's own search needs about twice it (simulation).
- **"Are CAL and TEST independent?"** A random split of one pool: exchangeable, but 27% of TEST has a ≥ 50% identical
  CAL relative (BLASTp); pre-registered CAL-disjoint and EC-cluster sensitivity analyses.
- **"Is your BLASTp baseline fair to CLEAN?"** Same training sequences as CLEAN's index, BLAST+ 2.12.0 default scoring,
  E ≤ 10, best hit by bit score, same label mapping and the same scorer; exact commands and SHA-256 fingerprints are in
  `data/cleanx/identity/blast/manifest.json`. It is the configuration a user gets by default, and it also beats CLEAN
  under CLEAN's own reported metric on New-392 (0.574 vs 0.497).
- **"Weighted F1 is standard."** We don't call it wrong; we show what it cannot see and that it rewards
  over-prediction, and we recommend reporting it alongside micro/union-macro and the invisible-FP count.
- **"Isn't 1-NN just BLAST?"** Compared directly: ties CLEAN at high identity, loses at low identity.
- **"Why SGR and not conformal FDR like CPEC?"** Annotation pipelines ask "accept this top-1 call or abstain?", which
  is selective risk; set-valued FDR/hierarchical-loss control (CPEC, Boger et al.) is complementary.
- **"Was anything tuned on the test set?"** No: §7, frozen commit hash, dated pre-registrations and amendments.

---

## 13. Submission, licence and ethics checklist

**CLEAN licence (read the PDF; not legal advice):**
- [ ] §1.1 counts *"distribution or otherwise making the Software available to a third party"* as Commercial Use.
      **Keep `Prat5ek01/clean_research` private** (verified private on 2026-09-28), and release only the `cleanx`
      add-on (no CLEAN source, weights, `split*.csv` or CLEAN-derived embeddings). Ask UIUC (otm@illinois.edu)
      before any public release.
- [ ] §2.5–2.6: derivative work renamed (`cleanx`), notice in `app/cleanx/__init__.py`.
- [ ] §2.7.1: Acknowledgements verbatim: *"The CLEAN Software, used by [your name], was developed by the
      Department of Chemical and Biomolecular Engineering at the University of Illinois Urbana-Champaign."*
- [ ] §2.7.2: email Prof. Huimin Zhao (zhao5@illinois.edu) about the publication.
- [ ] `app/build.py` says "Apache-2.0"; the PDF licence governs.
- [ ] **EC-Bench** (MIT): cite Davoudi et al. 2026 and the commit; its files are downloaded at run time (not
      redistributed). **EnzPlacer** (GPL-3.0) only if the optional audit is run: call it, do not copy its code.

**Journal (verify against the Guide):** abstract ≤ 250 words · highlights · keywords · graphical abstract ·
CRediT · competing interests · funding · data availability · **generative-AI declaration** (disclose AI
assistance in writing and code truthfully) · cover letter · suggested reviewers.

**Reproducibility:** environment versions (RUNBOOK §0) · data releases (UniProt 2026_03, ENZYME 02-Sep-2026,
Pfam 38.2, EC-Bench 5e98a2c, all in `data/cleanx/downloads/**/*.json`) · weight md5s (`cx_00_check.py`) · frozen
commit hash · Zenodo DOI for the add-on, accession lists and result CSVs (UniProt CC BY 4.0).

**Ethics:** no human or animal data; public databases only.

---

## 14. Literature (identifiers; ✔ = identifier seen on a search result or the code repository on the stated date)

| Reference | Identifier(s) | Code / weights | Evidence |
|---|---|---|---|
| Yu et al., Enzyme function prediction using contrastive learning, *Science* 379:1358 (2023) — CLEAN | 10.1126/science.adf2465 ✔ | github.com/tttianhao/CLEAN (research licence), weights on Google Drive | [code] |
| CLEAN-Contact, *Commun. Biol.* (2024) | 10.1038/s42003-024-07359-z ✔; bioRxiv 10.1101/2024.05.14.594148 ✔ | github.com/pnnl-predictive-phenomics/clean-contact @0e423db | [code] metric; [search] text |
| Ma, Joshi, Friedberg, Li, How not to be seen (EnzPlacer), *Bioinformatics* 42 (Suppl. 1) btag215 (2026) | 10.1093/bioinformatics/btag215 ✔; PMID 42412796 ✔; bioRxiv 10.64898/2026.02.23.707489 ✔ | github.com/drxiangma/EnzPlacer @97a8f0a (GPL-3.0; weights in repo; reference embeddings Zenodo 10.5281/zenodo.18110452) | [code] data/extraction; no eval code |
| Huang, Lin, Lan, Huang, Zhong, GloEC, *Brief. Bioinform.* 25(5) bbae365 (2024) | 10.1093/bib/bbae365 ✔; PMID 39073830 ✔ | github.com/hyr0771/GloEC @b58647c | [code] macro |
| Cao et al., EnzymeHunter (2026; Cell Press journal, pii S2666389926000760 — journal inferred from ISSN, confirm) | PMID 42453691 ✔; Zenodo 18896678 ✔; SSRN 5928085 ✔ | github.com/cgxbio/EnzymeHunter @fd2aa18 | [code] metric |
| Davoudi, Henry, Miller, Banaei-Kashani, EC-Bench, *Bioinform. Adv.* 6(1) vbag004 (2026) | 10.1093/bioadv/vbag004 ✔; PMID 41675694 ✔ | github.com/dsaeedeh/EC-Bench @5e98a2c (MIT; predictions published) | [code] metric + data |
| Yang et al., CARE benchmark (2024) | arXiv:2406.15669 ✔ | github.com/jsunn-y/CARE @e570c36 | [code] top-1 accuracy |
| Ding, Luo, Luo, CPEC, *PLoS Comput. Biol.* 20:e1012135 (2024) | 10.1371/journal.pcbi.1012135 ✔; PMID 38809942 ✔ | github.com/luo-group/CPEC @3c9d7a5 | [code] samples average |
| Boger, Chithrananda, Angelopoulos, Yoon, Jordan, Doudna, Functional protein mining with conformal guarantees, *Nat. Commun.* (2025) | 10.1038/s41467-024-55676-y ✔ | github.com/ronboger/conformal-protein-retrieval @509888f | [code] CLEAN experiments; [search] text |
| Dhibar, Basak, Jana, SOLVE, *Chem. Sci.* 16:18438 (2025) | 10.1039/D5SC04513D ✔; PMID 40951780 ✔ | github.com/saikat-ai/Enzyme_prediction @618681a | [code] macro |
| HIT-EC, *Nat. Commun.* 17:1146 (2026) | 10.1038/s41467-026-68727-3 ✔; bioRxiv 10.1101/2025.02.01.635810 ✔ | github.com/datax-lab/HIT-EC @7b45c7a (no eval code) | [search] |
| Han et al., Evidential deep learning for trustworthy EC prediction (ECPICK), *Brief. Bioinform.* 25(1) bbad401 (2024) | 10.1093/bib/bbad401 ✔ | github.com/datax-lab/ECPICK | [search] |
| FEDKEA — bioRxiv (2024); published as *Large language models enhance annotation of enzymes in metagenomes*, *Sci. Adv.* 12(34) (2026) | bioRxiv 10.1101/2024.08.12.604109 ✔; 10.1126/sciadv.aee4389 ✔; PMID 42616879 ✔ | github.com/Stevenleizheng/FEDKEA @602a11d (no eval code) | [search] |
| Rong, Zhong, Zheng, Hong, Liu, MAPred, *Brief. Bioinform.* 26(5) bbaf476 (2025) | 10.1093/bib/bbaf476 ✔; PMID 40966652 ✔; arXiv:2408.06391 ✔ | none found | [search] |
| RXNRECer (2026) | arXiv:2603.12694 ✔ | github.com/kingstdio/RXNRECer | [search] |
| ProtDETR, Interpretable enzyme function prediction via residue-level detection (ICLR 2025 workshop) | arXiv:2501.05644 ✔ | github.com/yangzhao1230/ProtDETR @aa49809 | [code] metric |
| Zhou et al., Multi-alignment contrastive learning for enzyme–reaction retrieval (FGW-CLIP) | arXiv:2512.08508 ✔ | — | [search] |
| Zhang, Ding et al., TIGER, ACL 2026 | arXiv:2605.24489 ✔ | — | [search] |
| The Accuracy Paradox … multi-label EC prediction (Sept 2026) | arXiv:2609.07897 ✔ | "with code" | [search] |
| Limitations of current ML models in predicting enzymatic functions for uncharacterized proteins, *G3* 15(10) jkaf169 (2025) | 10.1093/g3journal/jkaf169 ✔; bioRxiv 10.1101/2024.07.01.601547 ✔ | — | [search] |
| EnzymARC decoy benchmark (Aug 2026) | bioRxiv 10.64898/2026.08.21.746242 ✔ | — | [search] |
| PLMs outperform BLAST for evolutionarily distant enzymes (Apr 2026) | bioRxiv 10.64898/2026.03.31.715487 ✔ | — | [search] |
| ENZYMAP: modelling EC number changes in Swiss-Prot, *PLoS ONE* (2014) | 10.1371/journal.pone.0089162 ✔ | — | [search] |
| Omelchenko et al., Non-homologous isofunctional enzymes, *Biol. Direct* (2010) | PMID 20433725 ✔ | — | [search] |
| Riziotis et al., Convergent evolution in enzymes, *FEBS J.* (2025) | 10.1111/febs.17332 ✔ | — | [search] |
| Learning prototype classifiers for long-tailed recognition, IJCAI 2023 | 10.24963/ijcai.2023/151 ✔; arXiv:2302.00491 ✔ | — | [search] |
| Probabilistic conformal coverage guarantees in small-data settings (SSBC) | arXiv:2509.15349 ✔ | — | [search] |
| Ding et al., Class-conditional conformal prediction with many classes (2023) | arXiv:2306.09335 ✔ | — | [search] |
| Geifman & El-Yaniv, Selective classification for deep neural networks, NeurIPS 2017 (SGR) | arXiv:1705.08500 (check) | — | — |
| Angelopoulos et al., Learn then Test (2021) | arXiv:2110.01052 (check) | — | — |
| Guo et al., On calibration of modern neural networks, ICML 2017 | arXiv:1706.04599 (check) | — | — |
| Rives et al., ESM-1b, *PNAS* (2021); Steinegger & Söding, MMseqs2 (2017); UniProt, ENZYME, Pfam 38.2 | check latest | — | — |

---

## 15. What changed on 2026-09-28 (and why)

- **New E5** (evaluation card + EC-Bench re-scoring) turns C1 from a CLEAN-specific observation into a proof plus
  a field-wide, code-verified finding; pre-registered, piloted on EC-Bench Price-149, amended once (dated).
- **New E6** (certification budget) turns C4's likely negative gate (rare tiers not certifiable) into a
  quantitative, actionable result; pre-registered, piloted.
- **Boger et al. 2025 added** as the closest prior work for C4 (conformal guarantees on CLEAN distances).
- **Rigor fixes before the confirmatory run:** Holm (claimed before, now implemented), observed deltas,
  sensitivity analyses for clustered data (LPMO family, 7-EC tiers on Price-149, CAL–TEST overlap), E3 margin
  decomposition, E4 pure/mixed open set and a correctly certified EC3 fallback.
- **Model vintage:** EnzPlacer's training set contains half of the temporal TEST set; any second-model audit must
  be vintage-matched.
- **Safeguards added the same day:** the four EC-Bench files are pinned by SHA-256 in `cleanx/ecbench.py` (the
  confirmatory run must read the bytes the gate and pilot used); `cx_09 --card` writes to `E5_card/`, never the
  confirmatory `E5/` folder; `cx_08` draws a figure from confirmatory results only when their files exist, else
  from the pilot; test suite 22/22.
- **After the confirmatory run (same day):** results and verdicts recorded (§2b); audit of all analysis code found
  one bug (EC-cluster bootstrap, per-protein F1 only) and two reporting gaps (realized-risk CIs; infinite E6
  bounds), all fixed with tests (29/29); the fallback aligner replaced by BLASTp (§7), which revealed that a
  standard BLASTp top hit beats CLEAN on New-392 and the temporal TEST set; claims narrowed accordingly (S7,
  over-prediction reward, class-size mechanism, per-tier certification usefulness). Figures revised so that none
  can mislead (CIs everywhere, like-for-like comparisons, the primary endpoint drawn; §7.7); H3a, H4b and the E6
  confidence intervals filled in from the confirmatory CSVs (they do not depend on the BLASTp change).
- **Dropped earlier:** chemistry-weighted evaluation (EC-Bench), training-space misassignment (multi-EC
  artifact), dispersion-scaled centroids and BIC multi-prototypes (failed pilots).
