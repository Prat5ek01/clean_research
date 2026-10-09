# Prompt 2 of 2: Methodology (Materials and methods)

Copy everything below the line into the writing assistant, and attach `paper/references.bib`. This prompt is
independent of Prompt 1. Its references do not overlap with those of Prompt 1.

---

## Your role

You are a co-author writing for *Computational Biology and Chemistry* (Elsevier; author–year citations with
natbib, `\citet{}` / `\citep{}`). Write **Section 3, "Materials and methods"**, in LaTeX. Make it detailed enough
that a competent reader could re-implement every analysis without our code. You know the repository, but use
only the facts in this prompt. Do not add any number, setting or reference that is not given here.

Length: 2,800–3,600 words in the main text. Move long definitions and tables to an "Appendix / Supplementary
methods" block at the end of your output. Use past tense for what was done and present tense for definitions.

---

## 1. What the study does (context for the writer)

The paper audits CLEAN, a contrastive-learning EC-number predictor that embeds sequences with ESM-1b and assigns
the EC numbers of the nearest training-set centroids. It asks what CLEAN's reported numbers do and do not
measure, along four factors:

- the **metric** (C1);
- the **test set**, including leakage and how old the model's training data is (C2);
- the **inference rule and the homology baseline**, plus the **provenance of the labels** (C3, with E7);
- **reliability**: calibration, certified abstention and the calibration data that certification would need
  (C4).

CLEAN was **not retrained**. All analyses use the released weights, the released training embeddings and
predictions computed from them. Sections 2 and 4 of the paper (Related work and Results) are written
separately. Refer to models and benchmarks introduced in Section 2 by name ("CLEAN", "EC-Bench", "EnzPlacer",
"CARE") **without citing them again**.

---

## 2. Required structure and content

### 3.1 Study design, pre-registration and data roles

- **Pre-registration.** Every confirmatory analysis was pre-registered in the analysis code (the docstring of
  each script states hypotheses, endpoints and decision rules). The code was frozen by a git commit before the
  first confirmatory run: commit `c8e42f7`, 28 September 2026, for E1–E6. Changes made after that run are
  recorded with dates, and none altered a pre-registered endpoint.
  - E7 has its own dated pre-registration, drafted 29 September 2026 and finalised 7 October 2026, before any
    label evidence was downloaded. Its script refuses to run in confirmatory mode unless the analysis code is
    committed, and it writes the commit hash and the SHA-256 digests of the analysis files into every output.
- **Data roles (state them exactly; put them in a table).**
  - New-392 is the **development set** (disclosed).
  - Price-149 is confirmatory for methods frozen before 27 September 2026, and post hoc for the class-balanced
    rules designed later.
  - The SP-temporal **calibration split (CAL)** is the only data used for tuning or calibration.
  - The SP-temporal **test split (TEST)** was never used for any choice.
  - EC-Bench's Price-149 predictions served as the pilot for E5; its test-set predictions are confirmatory.
- **Metrics allowed for endpoints:** micro-F1, union-label macro-F1 and per-protein F1. Support-weighted and
  true-label macro scores are reported only for comparability, labelled as blind to out-of-label-set false
  positives and gameable.

### 3.2 The audited model and exact reproduction

- CLEAN, as released, uses the split100 model.
- **Inputs:** ESM-1b mean embeddings (layer 33, 1,280 dimensions; sequences truncated to 1,022 residues). They
  were recomputed locally; an equivalence test against the released embeddings gave a maximum absolute
  difference of 1.76 × 10⁻⁷.
- A learned projection to a 128-dimensional space, followed by the **maximum-separation** rule, which assigns
  the EC numbers whose centroid distance is separated from the rest.
- **Confidence:** CLEAN's released two-component Gaussian-mixture score on embedding distance.
- **Reproduction:** 0 of 392 prediction mismatches on New-392, and precision, recall and F1 of
  0.5955 / 0.4791 / 0.4967, bit-identical to CLEAN's own evaluation function. The Price-149 confidences were
  identical to four decimal places.
- Our code imports CLEAN's functions rather than copying them. CLEAN is licensed for research use only.

### 3.3 Datasets (with a Table 1)

- **Training reference: split100.** 227,362 sequences, 5,242 EC numbers and 241,025 (EC, sequence) pairs,
  from UniProt release 2022_01.
  - The released embedding file stores rows in per-EC blocks, so the training accession behind a row inside a
    block is not recoverable.
  - A protein with several EC numbers appears in several blocks, so co-annotated EC numbers were masked in
  every "other EC" computation.
  - Robustness backbone: split70 (79,172 pairs, 5,093 EC numbers).
- **New-392** and **Price-149** are CLEAN's original test sets. Price-149 uses RefSeq identifiers and has low
  identity to the training data. New-392 contains 5 sequences identical to training sequences.
- **SP-temporal benchmark.** Construction rules, in order, applied to Swiss-Prot entries created on or after
  1 March 2022 (UniProt release 2026_03; ENZYME release of 2 September 2026):
  - not a fragment; length 30–1,022;
  - accession not in any CLEAN split, New-392, Price-149 or CLEAN's other released sets;
  - sequence not identical to any CLEAN training sequence; within-set duplicates removed;
  - at least one complete EC number. Entries with only partial EC numbers, or with an enzyme keyword but no EC
    number, were excluded.

  The result: 8,920 raw entries gave 2,099 enzymes. Split those 2,099 into CAL (839) and TEST (1,260) by a
  stratified random 40/60 split by tier, with seed 20260927. TEST has 876 seen-EC proteins and 384 unseen-EC
  proteins (331 pure, 53 mixed). Add 2,000 non-enzymes, sampled from 4,058 eligible entries.
  - Labels are current complete EC numbers.
  - CLEAN's 2022-era predictions were mapped to current nomenclature through the ENZYME transfer table before
    scoring. Some reclassifications have no recorded transfer; such proteins are analysed as open-set.
- **Rarity tier** of a protein is the training count of its rarest true EC number in split100, after transfer
  mapping:
  - unseen: 0
  - singleton: 1
  - very rare: 2–3
  - rare: 4–10
  - common: more than 10

  A protein is "seen-EC" (closed) if all its true EC numbers are in CLEAN's vocabulary. "Pure" unseen proteins
  have no seen EC number; "mixed" unseen proteins have at least one.
- **EC-Bench predictions.** Per-protein predictions of 16 methods on EC-Bench's test set (468 proteins) and on
  Price-149, under its 100% and 30% identity regimes.
  - Fetched from a pinned commit (`5e98a2c`) and verified by SHA-256 digests.
  - Labels parsed exactly as in EC-Bench's own scoring code.
  - Reproduction gate: our scorer reproduced EC-Bench's printed support-weighted F1 for all 16 methods at EC
    levels 1 and 4, with maximum absolute difference 0.
- **Known limitations to state in this subsection:**
  - one family dominates the unseen tier: fungal lytic polysaccharide monooxygenases, 78 proteins with EC
    1.14.99.56;
  - CAL and TEST are a random split of one pool, so 27% of TEST proteins are ≥ 50% identical to a CAL protein.
    A CAL-disjoint sensitivity subset addresses this.

### 3.4 Metrics, the formal blind spot and the EC evaluation card (C1)

- Define, for EC depth d ∈ {1, 2, 3, 4} (prefix truncation):
  - micro-precision, micro-recall and micro-F1 (pooled TP/FP/FN);
  - **union-label macro-F1**: the mean over every class that is true or predicted;
  - **per-protein F1**: the mean over proteins of 2|T∩P| / (|T| + |P|);
  - support-weighted F1: the class F1 weighted by true support, as in sklearn's `average="weighted"`;
  - true-label macro-F1.

  Cite `Zhang2014` for the multi-label averaging schemes and `Sokolova2009` for invariance properties of
  classification measures.
- **Proposition 1 (state and prove in the Appendix).** Support-weighted precision and F1, and true-label macro
  scores, on any label space give zero weight to a class with no true instance. They are therefore invariant to
  every prediction of an EC number absent from the test labels. Scores computed after projecting predictions
  onto the test label set are invariant to such predictions by construction. Support-weighted recall equals
  micro-recall, so all of the distortion is in precision.
  - Verified numerically on 82 real prediction tables and in unit tests.
- **Invisible false positives:** the count of false positives on EC numbers outside the test label set.
- **Over-prediction test (E1):** emit the top-k EC numbers per protein, k = 1–10, and track each score.
- **Code audit.** We read the public evaluation code of CLEAN, CLEAN-Contact, ProtDETR, EnzymeHunter and
  EC-Bench at pinned commits. Report file and line numbers in a supplementary table. No citations are needed
  here; the papers are cited in Section 2.
- **EC evaluation card.** For any prediction file, at each depth, report:
  - micro-F1, union-label macro-F1 and per-protein F1;
  - the support-weighted score with its invisible false-positive count;
  - abstentions;
  - the label-space sizes.
- **E5 (EC-Bench re-scoring), pre-registered.**
  - **H5:** at least one of the 120 method pairs (test set, 100% regime, EC level 4) is ordered one way by
    support-weighted F1 and significantly the other way by micro-F1. "Significantly" means a paired bootstrap
    with B = 2,000, Benjamini–Hochberg adjustment over 120 pairs, q < 0.05, and the bootstrap median agreeing in
    sign with the observed difference.
  - **Key secondary (amended after the disclosed pilot, before confirmatory scoring):** Spearman ρ between
    inflation (support-weighted F1 minus micro-F1) and mean EC numbers predicted per protein; one-sided.
  - **Other secondaries:** Kendall τ_b with a bootstrap CI; invisible-FP shares; the 30% regime; EC levels 1–3;
    a complete-EC-only variant.
  - PRIAM is one of the 16 methods; cite it here as `ClaudelRenard2003`.

### 3.5 Inference geometry in CLEAN's frozen space (C3, E2)

- Using the same embeddings and the same maximum-separation rule throughout, compare these scoring rules:
  - **centroid** (CLEAN's own rule): distance to the mean of the EC's training members. This is the prototype
    classifier of `Snell2017`. CLEAN's embedding is trained with a supervised contrastive objective
    (`Khosla2020`).
  - **exemplar**: the distance to the nearest member (1-NN).
  - **kmean3**: the distance to the mean of the 3 nearest members.
  - **multiproto**: silhouette-gated k-means sub-prototypes, used only for EC numbers with at least 20 members:
    k ≤ 4, silhouette ≥ 0.25, each cluster at least 10% of members and at least 3 members.
  - **balanced1 / balanced3**: exemplar or kmean3 distances multiplied by n_c^λ, where n_c is the EC's training
    size. This is an equal-prior density correction motivated by prototype methods for long-tailed recognition
    (`Sharma2023`). λ was chosen on **CAL only**, maximising union-label macro-F1 at EC level 4 on seen-EC
    proteins, over the grid {0, 0.01, 0.02, 0.03, 0.05, 0.075, 0.1, 0.15, 0.2}.
- **Primary endpoint (H3b, pre-registered 27 September 2026, metric switched to union-label macro-F1 before any
  TEST result existed):**
  - balanced3 − centroid in union-label macro-F1 at EC level 4;
  - on TEST seen-EC proteins (n = 876);
  - paired bootstrap, B = 2,000, two-sided α = 0.05.

  Key secondary: micro-F1. Per-tier tests are Holm-corrected across tiers. Robustness: the split70 backbone.
- **Sensitivity analyses:**
  - cluster bootstrap by each protein's rarest true EC number, because proteins that share an EC share a
    centroid;
  - at most 5 proteins per rarest EC (seeded);
  - the CAL-disjoint TEST subset (< 50% identity to every CAL protein);
  - proteins with evidence at protein level.
- **E3 (Supplement only):** dispersion of EC numbers whose training members fall into non-homologous Pfam-clan
  groups, i.e. candidate non-homologous isofunctional enzymes (`Omelchenko2010`). Use a size-stratified
  permutation test (10,000 permutations), and test d_own and d_other separately.

### 3.6 Homology baseline and identity strata

- **BLASTp:** the BLAST+ 2.12.0 implementation (`Camacho2009`) of gapped BLAST (`Altschul1997`).
  - Default scoring: BLOSUM62, gap open 11, gap extend 1, composition-based statistics.
  - E-value ≤ 10, up to 500 subjects per query.
  - Database: CLEAN's split100 training sequences, i.e. the same training set as CLEAN.
- **Prediction:** the EC numbers of the best hit by bit score. Ties are broken by lower E-value, then by subject
  accession. A query with no hit abstains.
- **Identity:** identical residues divided by the length of the shorter sequence, maximised over HSPs.
  - Five bins: < 30%, 30–50%, 50–70%, 70–90% and ≥ 90%.
  - The same search, run against CAL, defines the CAL-disjoint subset.
  - Hit tables are stored with SHA-256 digests of the sequences searched. The analysis refuses tables whose
    sequences differ from the local data.
- State plainly that this replaced an approximate aligner after the confirmatory run. The replacement had been
  planned before the run, and it changes a comparator, not a pre-registered endpoint. The BLASTp-versus-CLEAN
  comparisons are paired bootstraps and are descriptive (not pre-registered).

### 3.7 Label provenance (E7, pre-registered; results pending)

- **Data.** The current UniProtKB entry of every benchmark protein: New-392, CAL, TEST and the EC-Bench test set.
  - Retrieved through the UniProt REST API as JSON, in batches of 100, from a single release recorded with
    SHA-256 digests.
  - Accession resolution order: primary accession; else the unique active entry listing the query as a
    secondary accession; else an entry the query was merged into. Everything else is "not found".
  - Price-149 is excluded: its proteins are not UniProtKB entries, and its labels come from mutant-fitness
    experiments.
- **Evidence of a label.** The Evidence and Conclusion Ontology (ECO) codes attached to that EC number in the
  protein names (including domains) and in catalytic-activity annotations.
- **EC-level classes**, in order of precedence:
  - **EXP:** ECO:0000269, experimental;
  - **SIM:** ECO:0000250, "by similarity", which cites the source protein;
  - **RULE:** ECO:0000255, ECO:0000256 or ECO:0000259, a match to a sequence model such as a HAMAP rule;
  - **OTHER:** everything else.
- **Protein-level class:**
  - LABEL_ABSENT if the entry is not found, or any benchmark label is missing from the current entry;
  - else SIM if any label is SIM;
  - else EXP if every label is EXP;
  - else RULE if any label is RULE;
  - else OTHER.

  **PROPAGATED** = SIM ∪ RULE. If you want a citation for using experimental evidence codes to define reference
  labels, the CAFA reports `Jiang2016` and `Zhou2019` may be cited for time-delayed evaluation on newly
  accumulated annotations. Do not claim more than that. UniProt and ECO need citations that are not in the
  verified list: insert `[CITATION NEEDED: UniProt Consortium; ECO — verify]`.
- **Population:** TEST seen-EC proteins (n = 876). The comparison is BLASTp top hit versus CLEAN, at EC level 4.
- **Primary endpoints** (a family of two, Holm-corrected; paired bootstrap B = 2,000):
  - **P7a:** Δ_EXP = micro-F1(BLASTp) − micro-F1(CLEAN) on EXP proteins; two-sided; seed 20260927.
  - **P7b:** the difference-in-differences DiD = Δ_PROPAGATED − Δ_EXP. EXP and PROPAGATED proteins are
    resampled independently, with seeds 20260927 and 20260928.
- **Gate G7.** The primaries are confirmatory only if n_EXP ≥ 100 and n_PROPAGATED ≥ 50; otherwise every E7
  result is descriptive. Effects are observed full-sample differences with 95% percentile intervals, plus the
  minimum detectable effect at 80% power, MDE80 = 2.80 × bootstrap SE.
- **Secondaries:**
  - **S7a:** SIM and RULE DiDs separately, Holm over the two.
  - **S7b:** identity-standardised P7a and P7b. Per-bin deltas are weighted by the pooled EXP + PROPAGATED bin
    distribution (direct standardisation). Bins with fewer than 10 proteins in either class are merged by a fixed
    rule: from the highest bin downward into the next lower bin, then the lowest upward. Each (class, bin) cell
    is resampled independently.
  - **S7c:** per-protein F1 instead of micro-F1.
  - **S7d:** the EC-Bench test set (100% regime; proteins whose EC-Bench labels are all complete): DiD for
    BLASTp minus EC-Bench's CLEAN retraining. Descriptive unless both classes have ≥ 30 proteins.
- **Mechanism (SIM proteins):**
  - the share whose BLASTp top hit is a source protein cited by their SIM label;
  - the share with the cited source among the kept hits;
  - the share with the cited source in CLEAN's training set;
  - BLASTp and CLEAN micro-F1 when the top hit is or is not a cited source.

  Use exact Clopper–Pearson 95% CIs (`Clopper1934`).
- **Sensitivity:** EC-cluster bootstrap; at most 5 proteins per rarest EC number.
- **Interpretation table** (fixed before the run; reproduce it as a table):
  - G7 fails: descriptive only.
  - P7a significant positive, P7b not: BLASTp's advantage holds on experimental labels; no detectable
    provenance effect.
  - Both significant positive: the advantage holds on experimental labels and is larger on propagated labels;
    the provenance-matched estimate is Δ_EXP.
  - P7a not significant, P7b significant positive: the advantage is concentrated in propagated labels.
  - Neither significant: no significant difference on experimental labels and no detectable provenance effect.
  - Any significant negative effect is reported as such.
- **Pilot:** New-392 (the development set), through the same code. The pilot may change parsing and class rules
  (as a dated amendment), but not the endpoints.

### 3.8 Calibration, certified abstention and the certification budget (C4: E4, E6)

- Use top-1 predictions under the centroid, exemplar and balanced3 rules, and group proteins by the tier of the
  *predicted* EC number.
- **Calibration:**
  - expected calibration error (ECE) with 10 equal-width bins, and the Brier score;
  - CLEAN's released confidence versus isotonic regression (`Zadrozny2002`), fitted on CAL either globally or
    per predicted tier (tier-conditional);
  - reliability diagrams with exact per-bin CIs.
- **Certified abstention:**
  - Selection with Guaranteed Risk (SGR; `Geifman2017`) on the raw distance, with δ = 0.1 and α ∈ {0.1, 0.2,
    0.3};
  - one global threshold versus one threshold per predicted tier (Mondrian; δ Bonferroni-split over 4 tiers),
    following the class-conditional guarantees literature (`DingT2023`);
  - the realized selective risk on TEST, with error counts and exact CIs;
  - open-set behaviour: abstention on pure and mixed unseen-EC proteins and on non-enzymes; AUROC of each score
    for seen-EC enzymes versus non-enzymes;
  - an EC-level-3 fallback, certified on the CAL proteins that the EC-level-4 rule rejects.
- **E6 certification budget (pre-registered H6).** For tier t and risk α, the budget n_needed is the smallest
  number of calibration proteins at which the exact binomial power to certify selective risk α reaches 0.8.
  - It uses the most favourable coverage on a 0.1–0.9 grid and SGR's per-test δ.
  - Selective risks are estimated on CAL only, as the Jeffreys posterior mean (errors + 0.5)/(accepted + 1),
    using coverages with at least 10 accepted proteins. This rule was set after the disclosed pilot.
  - Because the coverage is chosen with hindsight, n_needed is a lower bound on what SGR itself needs, which
    makes "needed > available" conservative.
  - **H6:** at α = 0.2 the budget exceeds the CAL proteins available in the very-rare and singleton predicted
    tiers.
  - **Endpoint:** the deficit (n_needed − n_available) with a bootstrap 95% CI (1,000 resamples; resamples that
    cannot certify count as infinite).
  - Convert it to years of Swiss-Prot curation: about 466 new enzymes per year pass the benchmark filters
    (March 2022 to September 2026), multiplied by the tier's share of CAL predictions.
  - Add an SGR simulation check.

### 3.9 Planned extensions (write these subsections ONLY if the author confirms that they were run)

Mark each subsection with `% [PLANNED — remove if not run]`.

- **E8, the benchmark card:**
  - for every benchmark and method pair, the minimum detectable effect, MDE80 ≈ 2.80 × paired-bootstrap SE, and
    its scaling with sample size, MDE(n) = MDE(N)·√(N/n);
  - exact accession and exact-sequence overlap between each benchmark (New-392, Price-149, the EC-Bench test,
    CAL, TEST, CARE's test sets, EnzPlacer's experimental unseen test) and each model's training set (CLEAN
    split100, EnzPlacer's seen/unseen training sets, CARE's training set);
  - label drift of New-392's false positives (descriptive).
- **E9, certifiable depth:** for each tier, the deepest EC level d ∈ {4, 3, 2, 1} whose certification deficit
  has an upper 95% bound ≤ 0 at α = 0.2. Pre-registered before any depth-level CAL statistic is computed.
  Secondary: the realized TEST risk at that depth.
- **E10, a vintage-matched second model:**
  - TEST proteins created on or after the first day of the second month after the newest TEST/CAL entry found
    in EnzPlacer's training data (expected cut-off 1 June 2025; about 363 proteins, 261 seen-EC), and not
    themselves in that training data;
  - EnzPlacer's released seen-EC model with its k-nearest-neighbour rule (L2 distance on unit-normalised
    projected ESM-1b embeddings, k = 7, nearest-neighbour vote), re-implemented and verified against its own
    command-line output on its example data;
  - **primary:** micro-F1 of BLASTp against EnzPlacer's training sequences minus micro-F1 of EnzPlacer; paired
    bootstrap; two-sided.

### 3.10 Statistical analysis

- Paired nonparametric bootstrap over proteins, with the same resamples for every method compared: B = 2,000,
  seed 20260927.
- Effects are the **observed** full-sample differences, with 95% percentile intervals and two-sided bootstrap
  p-values. A p-value of 0 is reported as p < 1/B. The bootstrap mean of a macro-averaged metric is biased, so
  it is never reported as the effect.
- **Multiplicity:**
  - Holm (step-down, family-wise) across tiers and within each pre-registered family;
  - Benjamini–Hochberg across the 120 EC-Bench pairs;
  - for both, insert `[CITATION NEEDED: Holm 1979; Benjamini & Hochberg 1995 — verify]`.
- **Exact binomial (Clopper–Pearson) intervals** for risks and shares.
- Cluster bootstrap for proteins that share an EC number. This was corrected after the confirmatory run so that
  per-protein F1 is divided by each resample's total weight; say so.
- Report everything as recommended for supervised ML validation in biology (`Walsh2021`): data, optimisation,
  model and evaluation.

### 3.11 Software, reproducibility and licensing

- **Versions:** Python 3.10, numpy 1.22.3, pandas 1.4.2, scipy 1.7.3, scikit-learn 1.2.0, torch 1.11.0 (CPU),
  fair-esm 2.0.0, biopython 1.83, matplotlib 3.7.0, BLAST+ 2.12.0.
- Every script has a smoke mode with known-answer checks. Unit and integration tests run the real scripts on
  synthetic data in temporary folders.
- External files (EC-Bench, UniProt responses) are pinned and verified by SHA-256 digests.
- The analyses ran on a Windows workstation with 15 GB of RAM and no GPU. The EC-Bench part of E5 and all of E6
  were re-run on Linux with identical results.
- CLEAN's licence permits research use only, so no CLEAN weights, embeddings or split files are redistributed.
  The analysis code will be released at publication (state the repository policy as `[TO CONFIRM]`).

---

## 3. References allowed in this section (16; each cited once)

These do not overlap with the 30 references of the Introduction and Related work. Use each **exactly once**, at
the place assigned, and claim only what is stated.

| # | Key | Place | What you may cite it for |
|---|---|---|---|
| 1 | Zhang2014 | 3.4 | Multi-label metric definitions (micro, macro, example-based averaging) |
| 2 | Sokolova2009 | 3.4 | Classification measures differ in which confusion-matrix changes they are invariant to |
| 3 | ClaudelRenard2003 | 3.4 | PRIAM: enzyme-specific profiles for genome annotation (one of the EC-Bench methods) |
| 4 | Snell2017 | 3.5 | The class-mean (prototype) classifier, which CLEAN's centroid rule instantiates |
| 5 | Khosla2020 | 3.5 | Supervised contrastive learning (CLEAN's training objective family) |
| 6 | Sharma2023 | 3.5 | Prototype classifiers for long-tailed recognition (motivation for class-balanced rules) |
| 7 | Omelchenko2010 | 3.5 (E3) | Non-homologous isofunctional enzymes: unrelated proteins catalysing the same reaction |
| 8 | Altschul1997 | 3.6 | Gapped BLAST algorithm |
| 9 | Camacho2009 | 3.6 | The BLAST+ software suite (version used: 2.12.0) |
| 10 | Jiang2016 | 3.7 | CAFA2: time-delayed evaluation of function prediction on newly accumulated annotations |
| 11 | Zhou2019 | 3.7 | CAFA3: the continuation of that evaluation, with experimental screens |
| 12 | Clopper1934 | 3.7 | Exact binomial confidence limits |
| 13 | Zadrozny2002 | 3.8 | Isotonic regression for calibrating classifier scores |
| 14 | Geifman2017 | 3.8 | Selective classification with a guaranteed risk (SGR), holding with high probability |
| 15 | DingT2023 | 3.8 | Class-conditional conformal prediction is hard with many classes and few examples per class |
| 16 | Walsh2021 | 3.10 | DOME recommendations for supervised ML validation in biology |

Not verified, and therefore not in `references.bib`: UniProt, ENZYME, the Evidence and Conclusion Ontology,
Holm (1979), Benjamini–Hochberg (1995), Efron's bootstrap and scikit-learn. Write `[CITATION NEEDED: …]`; never
supply bibliographic details yourself.

---

## 4. Rules (non-negotiable)

1. Use only the facts and numbers above. Do not invent a parameter, a software version, a dataset size or a
   reference. If something is missing, write `[TO CONFIRM: …]`.
2. **No citation is repeated**, and none of Prompt 1's 30 references is cited here. Models and benchmarks from
   Section 2 are named, not re-cited.
3. Present **no results** except the reproduction and gate checks listed in 3.2 and 3.3, which are part of the
   method. Everything else belongs to Section 4.
4. Keep "pre-registered", "post hoc", "descriptive" and "pilot/development" exactly as labelled here. Never
   upgrade a descriptive analysis to a confirmatory one.
5. Use the same terms as the rest of the paper: support-weighted F1, union-label macro-F1, micro-F1,
   per-protein F1, seen-EC / unseen-EC (pure, mixed), rarity tier, CAL, TEST, PROPAGATED.
6. Equations: number the metric definitions, Proposition 1, the DiD, the direct standardisation, the
   certification budget and MDE80. Put the proof in the Appendix.
7. Include **Table 1** (datasets and roles) and **Table 2** (pre-registered endpoints: ID, hypothesis, metric,
   population, test, correction, decision rule, date frozen).

## 5. Output

1. The LaTeX for `\section{Materials and methods}` with subsections 3.1–3.11, plus an Appendix with the proof
   of Proposition 1 and the full E7 interpretation table.
2. Tables 1 and 2 in LaTeX.
3. A map from each of the 16 keys to the sentence that cites it, confirming each is used exactly once.
4. A list of every `[CITATION NEEDED]`, `[TO CONFIRM]` and `[PLANNED]` marker.
5. Word count of the main text.
