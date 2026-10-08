# Prompt 1 of 2: Introduction and Related Work (literature review)

Copy everything below the line into the writing assistant, and attach `paper/references.bib`.

---

## Your role

You are a co-author writing for *Computational Biology and Chemistry* (Elsevier; author–year citations with
natbib, `\citet{}` / `\citep{}`). Write **Section 1 (Introduction)** and **Section 2 (Related work)** of a research
paper, in LaTeX, in plain and precise scientific English. You know the code repository, but rely only on the
facts in this prompt. Do not add any number, claim or reference that is not given here.

Length: Introduction 900–1,100 words; Related work 1,600–1,900 words in five subsections. No bullet lists in
the running text, except for the contributions list at the end of the Introduction.

---

## 1. The project, in detail

**Working title:** *What reported scores hide: metric blind spots, test-set vintage, label provenance and the
calibration budget for trusting rare enzyme-function calls.*

**Author:** Prativa Khanal, Butwal Multiple Campus, Butwal, Nepal.

**What the paper is.** It is an audit and extension of **CLEAN**, a contrastive-learning EC-number predictor
published in *Science* in 2023. CLEAN embeds each protein sequence with the ESM-1b protein language model, maps
it through a contrastively trained network, and assigns the EC numbers whose training-set centroid (the mean
embedding of all training proteins with that EC number) is closest, using a "maximum-separation" rule. CLEAN's
paper reported better accuracy than BLASTp. Several later models reuse its design and its evaluation code. We
did not retrain CLEAN. We reproduced the released model exactly: 0 of 392 prediction mismatches on its New-392
test set, and the published precision, recall and F1 (0.596 / 0.479 / 0.497). We then asked what CLEAN's
reported numbers do and do not measure. The paper is not a new architecture. It is an evaluation, reliability
and method paper, with a new benchmark, a formal result and open (private-until-publication) code.

**The central story.** A reported score is a product of four things, not one: **the model, the metric, the test
set, and the labels.** For deployment there is a fifth: **whether the model's confidence can be trusted where
it matters most, on rare EC numbers.** For EC predictors each of the four non-model factors had been taken for
granted. We test each one, and each changes the picture.

1. **The metric (C1).** CLEAN's public evaluation code, three later methods that copy it (CLEAN-Contact,
   ProtDETR, EnzymeHunter) and the EC-Bench benchmark use support-weighted (sklearn `average="weighted"`)
   precision and F1. We prove that these scores give zero weight to any predicted EC number that does not occur
   among the test labels, so such false positives are never counted. Support-weighted recall equals micro
   recall, so all the distortion is in precision. On CLEAN's own test sets the reported precision hides 67%
   (New-392) and 96% (Price-149) of CLEAN's false positives (74% on our new benchmark). The score can be gamed:
   predicting the top 10 EC numbers per protein raises the reported F1 on Price-149 to 0.618, above CLEAN's
   0.495, while micro-F1 falls to 0.127. When we re-scored the published per-protein predictions of the 16
   methods in EC-Bench, the inflation tracked how many EC numbers a method predicts per protein (Spearman
   ρ = 0.75, pre-registered). One published pairwise ranking reversed significantly: DeepECtransformer ranks
   above PRIAM by the reported F1, but below it by micro-F1 (−0.041 [−0.075, −0.006], BH q = 0.026). Rankings
   are otherwise largely preserved (Kendall τ_b 0.90).
2. **The test set (C2).** We built a leakage-controlled temporal benchmark: 2,099 Swiss-Prot enzymes created
   after CLEAN's training snapshot (2022-03-01), with no shared accession or identical sequence, split into a
   calibration set (CAL, 839) and a test set (TEST, 1,260), plus 2,000 non-enzymes. **30.5% of the new enzymes
   carry an EC number that CLEAN has never seen**, which caps any closed-set model's attainable recall. A
   model's training "vintage" has to be matched model by model: the training set of a 2026 successor (EnzPlacer)
   already contains 631 of the 1,260 TEST proteins.
3. **Inference geometry and the homology baseline (C3).** Using CLEAN's frozen embedding space, member-level
   (nearest-neighbour) inference beats CLEAN's one-centroid-per-EC rule. The pre-registered primary endpoint, a
   class-balanced 3-nearest-member rule, gives +0.043 union-label macro-F1 [+0.030, +0.060]. **However, a
   standard BLASTp top hit against CLEAN's own training sequences beats every embedding-based rule**: micro-F1
   0.735 vs 0.641 for CLEAN on the 876 TEST proteins with seen EC numbers (+0.093 [+0.070, +0.116]). Its
   advantage lies at 30–90% sequence identity. Below 30% identity the two do not differ significantly (+0.048
   [−0.016, +0.113], n = 166). This contradicts CLEAN's own BLASTp comparison, which used a different set-up.
4. **The labels (C3, extension E7, pre-registered, result pending).** Swiss-Prot curators assign many EC numbers
   by transferring them from a characterised homologue ("by similarity", evidence code ECO:0000250) or from a
   family rule such as HAMAP (ECO:0000255). Only labels with experimental evidence (ECO:0000269) test a
   predictor against an independent reference. A homology-transfer baseline re-enacts the curators' process, so
   on propagated labels part of its agreement is built in (circularity). Our pre-registered test asks whether
   BLASTp's advantage holds on experimentally supported labels, and whether it is larger on propagated labels.
   **Do not state a result. Write: "we test whether…; results in Section 4".**
5. **Reliability and its price (C4).** CLEAN's released confidence (a Gaussian-mixture score on embedding
   distance) is most overconfident for rare EC numbers. The expected calibration error rises from 0.19 for
   common EC numbers to 0.32 for very-rare and singleton ones. Tier-conditional isotonic recalibration on the
   calibration split brings all tiers to 0.05–0.10. Exact binomial calculations show that certifying rare-EC
   predictions at a 20% error rate is not possible with the calibration proteins that exist. In the very-rare
   and singleton tiers the deficit's 95% lower bound is 528 and 515 proteins, while only 102 and 90 exist. At the
   current rate of curation (≈ 466 new enzymes per year that pass our filters) Swiss-Prot would need at least
   6–10 years to supply them.

**Planned extensions in the same paper** (mention them only as scope in the contributions paragraph, and only
if the author confirms they were run; otherwise leave them out): E8, a benchmark card reporting each
benchmark's minimum detectable effect and its exact overlap with each model's training data; E9, the deepest EC
level at which predictions can be certified, per rarity tier; E10, a vintage-matched comparison of a second
model (EnzPlacer) on TEST proteins created after its training data.

**Rigour that the Introduction should convey in one or two sentences.**
- Every confirmatory analysis was pre-registered in the code and frozen by a git commit before it was run.
- New-392 was used for development and is disclosed as such.
- Calibration and tuning used only the CAL split. TEST was never used for any choice.
- Paired bootstrap (B = 2,000), with Holm or Benjamini–Hochberg correction.

---

## 2. The narrative arc, paragraph by paragraph

### Introduction (6 paragraphs)

- **P1 — Why EC prediction matters.** The EC system classifies enzyme-catalysed reactions in a four-level
  hierarchy. Most newly sequenced proteins have no experimentally characterised function, so computational EC
  assignment is the default. Cite: `McDonald2023`.
- **P2 — The rise of learned predictors, and CLEAN.** CLEAN uses ESM-1b embeddings (`Rives2021`) and contrastive
  learning. It reported higher accuracy than BLASTp and became a reference model (`Yu2023`). Say what CLEAN does
  in two sentences.
- **P3 — The tension.** Independent, expert-curated evaluations of EC predictors have been less favourable than
  the original reports. For example, models largely failed to predict functions that were unseen in training
  (`deCrecyLagard2025`). Machine-learning-based science in general is vulnerable to evaluation leakage
  (`Kapoor2023`). Pose the paper's question: how much of a reported EC score reflects the model, and how much
  reflects the metric, the test set and the labels?
- **P4 — Deployment.** In practice users need to know when to trust a call, and modern neural networks are often
  miscalibrated (`Guo2017`). The calls most worth making, on rare EC numbers, are the ones with the least data
  to calibrate them.
- **P5 — What we do.** One paragraph. We reproduce CLEAN exactly and audit each factor in turn: metric, test
  set, inference rule and labels, then confidence. No citations.
- **P6 — Contributions.** A numbered list of 4–5 items (C1–C4, with E7 inside C3), each with **one** headline
  number from Section 1 above. Close with one sentence on pre-registration and code availability.

### Related work (5 subsections)

- **2.1 Homology transfer and the reliability of EC labels.** EC function is less conserved between homologues
  than assumed (`Rost2002`). Conservation depends on identity (`Tian2003`: roughly 40% identity for the first
  three EC digits and about 60% for all four; say "roughly" and attribute the figures). Misannotation is rare in
  Swiss-Prot but common in automatically annotated databases (`Schnoes2009`). Labels can be wrong at scale even
  in a single EC class (`Rembeza2021`: at least 78% of 122 sequences annotated as EC 1.1.3.15 were inferred to
  be misannotated). **Bridge to E7:** benchmark labels are themselves partly products of homology transfer. That
  matters for comparing a homology baseline with a learned model. (UniProt evidence codes need a UniProt or ECO
  citation that is not in the verified list: insert `[CITATION NEEDED: UniProt/ECO — verify]`.)
- **2.2 Learned EC predictors.** Give the arc from per-EC classifiers (`Dalkiran2018`, ECPred) to deep sequence
  models (`Ryu2019`, DeepEC; `Sanderson2023`, ProteInfer; `Kim2023`, DeepECtransformer, which assigned EC numbers
  to 464 unannotated *E. coli* genes and validated three). Then the CLEAN family: contact maps (`YangY2024`,
  CLEAN-Contact) and hierarchy-aware contrastive learning (`Cao2026`, EnzymeHunter). Note neutrally that the
  public code of CLEAN-Contact and EnzymeHunter reuses CLEAN's support-weighted evaluation; this is our code
  audit, so say "their public code", not "their paper". Close with models aimed at under-represented EC numbers
  (`Dumontet2026`, HIT-EC).
- **2.3 Independent tests and the homology baseline.** `Capela2025`: BLASTp was marginally better overall than
  protein-language-model classifiers, the two were complementary, and the language models helped below about
  25% identity. `Sathyamoorthy2026` (a preprint; say so): with a train-set-matched BLASTp, performance was at
  parity overall, but language models were well ahead for distant eukaryotic enzymes. **This partly
  contradicts our tie below 30% identity.** State the contradiction explicitly and say the outcome appears to
  depend on the identity profile and the baseline configuration. `Sartori2026` (a preprint): CLEAN and DIAMOND
  assigned the original EC number to more than 90% of decoys with weakly perturbed active sites, so both
  homology transfer and CLEAN over-assign.
- **2.4 How function predictors are evaluated: metrics, benchmarks and time.**
  - CAFA established time-delayed, prospective evaluation (`Radivojac2013`); LAFA runs it continuously
    (`Phan2026`).
  - Some GO evaluation metrics give erroneously high scores to sets of false positives (`Plyusnin2019`). This
    is the closest precedent for C1; say so plainly and state that we found no equivalent analysis for EC
    predictors ("to our knowledge").
  - Incomplete ground truth means that a scored false positive may be an unannotated true function
    (`Jiang2014`). This is a caveat to C1, so state it.
  - Two macro-F1 formulas in common use can rank classifiers differently (`Opitz2019`, an arXiv note).
  - EC-Bench unified the comparison of EC predictors and adopted support-weighted F1 because of class imbalance
    (`Davoudi2026`). CARE controls leakage with sequence-identity splits (`YangJ2024`).
  - Contrast these with our temporal hold-out of the *released* model, its unseen-EC share and the vintage
    problem.
- **2.5 Rare and unseen EC numbers, and trustworthy predictions.** EnzPlacer targets EC numbers absent from
  training (`Ma2026`). Conformal methods give guarantees for EC annotation: false-discovery-rate control in
  CPEC (`DingK2024`), and user-chosen risk guarantees in protein search including CLEAN's embeddings
  (`Boger2025`). These guarantees hold on average over test proteins, not per rarity tier, and none asks how
  much calibration data a per-tier guarantee would need. That is our C4 gap (hedge with "to our knowledge").
  End the section with a 3–4 sentence paragraph listing the open questions this paper answers, in the order
  of Section 1.

---

## 3. The 30 references: each cited exactly once in Sections 1–2

Use only these keys. **Each key appears exactly once** across the Introduction and Related work. Claim only what
the "verified claim" column says (claims were checked against abstracts, publisher summaries and code, not full
texts).

| # | Key | Place | Verified claim you may make |
|---|---|---|---|
| 1 | McDonald2023 | 1-P1 | EC numbers classify enzyme-catalysed reactions in a four-level hierarchy |
| 2 | Rives2021 | 1-P2 | ESM-1b protein language model (CLEAN's input embeddings) |
| 3 | Yu2023 | 1-P2 | CLEAN: contrastive learning on ESM-1b embeddings, nearest EC centroids; reported higher accuracy, reliability and sensitivity than BLASTp |
| 4 | deCrecyLagard2025 | 1-P3 | Expert dissection of DeepECtransformer predictions for >450 *E. coli* proteins of unknown function: ML largely fails to predict functions unseen in training; basic logic errors |
| 5 | Kapoor2023 | 1-P3 | Leakage undermines reproducibility in ML-based science |
| 6 | Guo2017 | 1-P4 | Modern neural networks are often poorly calibrated |
| 7 | Rost2002 | 2.1 | Enzyme function is less conserved between homologues than anticipated |
| 8 | Tian2003 | 2.1 | Roughly 40% identity is needed for >90% conservation of the first three EC digits; about 60% for all four |
| 9 | Schnoes2009 | 2.1 | Across 37 enzyme families, misannotation is low in Swiss-Prot and high in TrEMBL, GenBank NR and KEGG |
| 10 | Rembeza2021 | 2.1 | Screen of 122 sequences: ≥78% of EC 1.1.3.15 inferred misannotated; four alternative activities confirmed |
| 11 | Dalkiran2018 | 2.2 | ECPred: per-EC classifiers; its time-separated hold-out was very difficult (F1 0.14 at the substrate level) |
| 12 | Ryu2019 | 2.2 | DeepEC: deep (CNN) EC prediction, high quality and throughput |
| 13 | Sanderson2023 | 2.2 | ProteInfer: deep networks for protein function inference |
| 14 | Kim2023 | 2.2 | DeepECtransformer uses transformer layers; EC numbers for 464 unannotated *E. coli* genes, 3 validated |
| 15 | YangY2024 | 2.2 | CLEAN-Contact: CLEAN plus predicted contact maps; reported gains |
| 16 | Cao2026 | 2.2 | EnzymeHunter: hierarchy-aware contrastive learning for fine-grained EC prediction |
| 17 | Dumontet2026 | 2.2 | HIT-EC: hierarchical interpretable transformer, motivated by under-represented ECs and incomplete annotations |
| 18 | Capela2025 | 2.3 | BLASTp marginally better overall than PLM classifiers; complementary; PLMs useful below ~25% identity |
| 19 | Sathyamoorthy2026 | 2.3 | **Preprint.** Train-set-matched BLASTp: parity overall; PLMs well ahead for distant eukaryotes |
| 20 | Sartori2026 | 2.3 | **Preprint.** CLEAN and DIAMOND assign the original EC to >90% of weakly perturbed active-site decoys |
| 21 | Radivojac2013 | 2.4 | CAFA: time-delayed community evaluation of function prediction |
| 22 | Phan2026 | 2.4 | LAFA: continuous longitudinal evaluation under evolving ground truth |
| 23 | Plyusnin2019 | 2.4 | Of 37 GO metrics, some give erroneously high scores to false-positive sets |
| 24 | Jiang2014 | 2.4 | Incomplete annotation biases evaluation; a scored FP may be a true, unannotated function |
| 25 | Opitz2019 | 2.4 | **arXiv note.** Two macro-F1 formulas can differ and rank classifiers differently |
| 26 | Davoudi2026 | 2.4 | EC-Bench: unified protocol for EC predictors; adopts support-weighted F1 for class imbalance |
| 27 | YangJ2024 | 2.4 | CARE: EC classification with identity-based out-of-distribution splits (<30%, 30–50%) |
| 28 | Ma2026 | 2.5 | EnzPlacer: contrastive learning to predict EC numbers absent from training |
| 29 | DingK2024 | 2.5 | CPEC: conformal EC annotation with limited false positives (FDR control) |
| 30 | Boger2025 | 2.5 | Conformal guarantees for protein mining with user-chosen risk, including CLEAN embeddings |

The remaining keys in `references.bib` belong to the Methods (Prompt 2). Do not use them here.

---

## 4. Rules (non-negotiable)

1. **No invented references, numbers, quotations or claims.** If a sentence needs support that is not in the
   table, write `[CITATION NEEDED: what]`. Never fill the gap yourself.
2. **No repeated citations.** Each of the 30 keys is cited exactly once in Sections 1–2, in the place assigned.
   Refer back with words ("the benchmark above", "CLEAN"), not a second citation.
3. **Label preprints** (Sartori2026, Sathyamoorthy2026) as preprints, and label Opitz2019 as an arXiv note.
4. **Hedge absence claims** with "to our knowledge". Three are allowed: that no EC study has analysed whether
   the averaging counts every false positive; that no study has compared centroid and neighbour rules in one
   frozen EC embedding space; and that no study has estimated the calibration data a per-tier guarantee needs.
5. **No overselling.** No "first", "novel paradigm", "groundbreaking" or "state of the art". The mathematics of
   C1 is elementary: the contribution is its consequence and its prevalence. C2 is partly anticipated by
   EC-Bench and CARE, so say what is different. The BLASTp result has a peer-reviewed precedent (`Capela2025`),
   so present ours as confirming it on the released CLEAN and a temporal hold-out.
6. **Discuss the contradiction** with Sathyamoorthy2026 honestly. Do not hide it.
7. **E7–E10:** describe as research questions, never as findings. Use exactly
   `[RESULT PENDING: E7]` wherever a result would go.
8. Use the same names throughout: "support-weighted F1" (never just "weighted F1"), "union-label macro-F1",
   "micro-F1", "per-protein F1", "seen-EC" and "unseen-EC" proteins, "rarity tier"
   (common / rare / very rare / singleton / unseen).
9. Never say CLEAN is "wrong" or "flawed". Say what its reported numbers do and do not measure.
10. Tone: neutral and specific; short sentences; one idea per paragraph.

## 5. Output

1. The LaTeX for `\section{Introduction}` and `\section{Related work}` (subsections `\subsection{}`), with
   `\citet{}`/`\citep{}` and the keys above.
2. A table that maps each of the 30 keys to the sentence that cites it. Use it to verify that every key is used
   exactly once.
3. A list of every `[CITATION NEEDED]` and `[RESULT PENDING]` placeholder.
4. A three-line self-check: word counts per section; whether any key is missing or repeated; whether any number
   is not taken from Section 1 of this prompt.
