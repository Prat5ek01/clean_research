# Literature package: Related Work for the CLEAN evaluation paper

**Target journal:** *Computational Biology and Chemistry* (Elsevier). **Drafted:** 2026-09-29.

**Files:**
- `paper/literature_review.tex`: Section 2, *Related work*, 1,640 words, ready to `\input`.
- `paper/references.bib`: 50 entries, every one cited, every one checked.
- This file: the package A–O your prompt asked for, plus integration notes (section P).

Nothing here is committed. It will be committed only when you say "commit and push".

---

## 0. Read this first

### 0.1 Scope

Your prompt template is written for a wet-lab chemistry BSc project, with materials, pH, temperature and local Nepali context. The manuscript it is meant for is the computational CLEAN evaluation paper for *Computational Biology and Chemistry*: Prativa Khanal (corresponding author) and the 239-word abstract in `paper/manuscript.tex`. So I adapted the template instead of forcing sections that don't apply:
- "Materials" became data sets.
- "Experimental variables" became scoring rules, inference rules, identity strata and frequency tiers.

If you meant a different project, send its title and description; nothing here would transfer.

### 0.2 How the sources were checked, and what could not be checked

The analysis environment cannot reach scholarly APIs or publisher pages. Crossref, doi.org, PubMed E-utilities, Semantic Scholar, OpenAlex, arXiv, bioRxiv and journal sites are all blocked. Every reference was therefore checked through web-search results: the publisher or preprint page, plus a second independent index where one was available (PubMed/PMC, OSTI, dblp, Semantic Scholar, institutional repositories).

- **V2:** a publisher or preprint page and an independent index agree on title, authors, year, venue and DOI.
- **V1:** only one source.

What each paper is cited *for* was checked against its abstract, the publisher's summary or a search excerpt. I did not read full texts. Your prompt asks that the most important papers be read in full. That was impossible here, so L3 lists nine papers to open before submission.

No LaTeX installation is available here, so the section was not compiled. Instead, a script (pybtex in strict mode plus a key checker) confirmed four things:
- the .bib parses;
- every `\cite` key exists;
- every entry is cited;
- braces balance.

The last automated run came before three small wording edits, which changed no citation: "strongly imbalanced" became "imbalanced", "Later methods add" became "Later methods introduce", and "Two lessons from this work" became "Two lessons from the evaluation of Gene Ontology predictors". The shell was unavailable after that, so section K was brought up to date by hand from the checked rendering.

### 0.3 Weaknesses to know before you use the review

Your prompt's rules 18 and 19 ask for these up front.

1. **The BLASTp result is not new as a general claim, and one preprint partly contradicts it.**
   - A peer-reviewed comparison already found BLASTp marginally better overall than protein-language-model classifiers (Capela et al., 2025, *BMC Bioinformatics*).
   - A 2026 preprint found language-model classifiers *ahead* of BLASTp for evolutionarily distant enzymes (Sathyamoorthy and Puri, 2026). We see only a tie below 30% identity.
   - CLEAN's own paper reports the opposite of our result: CLEAN better than BLASTp.
   - Present our BLASTp finding as confirmation in a leakage-controlled temporal setting for the *released* CLEAN (Type 3, section G), not as a discovery.
   - The CLEAN paper's BLASTp configuration still needs reading in its Methods (full text not reachable here). A reviewer will ask why the two results differ.
2. **C1's mathematics is elementary.** A support-weighted average gives zero weight to labels with zero support; scikit-learn's documentation describes the zero-support handling.
   - The closest published precedent is in Gene Ontology evaluation: Plyusnin et al. (2019) showed that some metrics give erroneously high scores to deliberately false annotation sets. Opitz and Burst (2019) showed that two macro-F1 formulas can reorder classifiers.
   - C1 is therefore **Type 2**. Its value lies in three things: the EC-specific formal statement, the prevalence in five public code bases, and the consequences for published rankings. Do not present it as a new theory of metrics.
3. **"Invisible false positives" are defined against Swiss-Prot labels, which are incomplete and sometimes wrong.**
   - Incomplete: Jiang et al. (2014). Sometimes wrong: Schnoes et al. (2009); Rembeza and Engqvist (2021).
   - Some predictions counted as false positives may be true functions nobody has annotated yet.
   - Where precision matters, say "predictions not supported by the reference annotation", and keep this in the Limitations.
   - The E1 label-drift file (`E1_label_drift.csv`) is the natural place to show how many such predictions later became annotations.
4. **The abstract says "the published EC-Bench predictions of 16 methods".**
   - EC-Bench describes *ten* existing methods. Our 16 are *predictors*: the ten methods, "learnt" and "regular" threshold variants of three of them, and two ensembles (majority voting with two variants, and stacking).
   - Suggested wording: "the published EC-Bench predictions of 16 predictors (ten methods and their variants)".
   - I have **not** changed the abstract; you asked me to change it only when needed, so the call is yours. It is a one-word fix that a reviewer comparing with EC-Bench would catch.
5. **The field is moving fast.**
   - 12 of the 50 references are from 2025–2026, and two of them are unrefereed preprints: EnzymARC and Sathyamoorthy and Puri.
   - EnzymeHunter, LAFA and HIT-EC were published in 2026.
   - Re-run the searches in section M in the week before submission.
6. **One model family.** The study evaluates the released, frozen CLEAN.
   - Reviewers may ask whether the conclusions carry over to CLEAN-Contact or EnzymeHunter.
   - The review makes claims about CLEAN only, and says "commonly" and "five public code bases" rather than "the field".
7. **Some bibliographic fields are incomplete.** L2 lists them: full author lists for eight entries, volume numbers for two 2026 articles, and the published HIT-EC author list. They are marked in `references.bib` and must be completed from the publisher pages; guessing them is not safe.
8. **Journal register.**
   - *Computational Biology and Chemistry* readers are computational biologists and chemists.
   - 10 of the 50 references are machine-learning venues (NeurIPS, ICML, IJCAI, KDD, arXiv). Each one carries a specific statistical or ML concept used in the paper (prototype classifiers, calibration, selective risk, class-conditional conformal prediction, macro-F1 variants).
   - That is defensible, but do not add more of them.
9. **Nepal/South Asia.**
   - A targeted search found no Nepal-specific literature on EC-number prediction or its evaluation. The topic has no local material, environmental or regional dimension, so no local literature is included for its own sake.
   - The one South Asian primary study in this area, SOLVE (Dhibar, Basak and Jana, 2025, IACS Kolkata, *Chemical Science*), is cited on its scientific merit.

---

## A. Research understanding

| Item | This study |
|---|---|
| Research problem | EC-number predictors are ranked with scores that may not count every wrong prediction. The released CLEAN model has not been evaluated on a leakage-controlled temporal set, with its inference rule, its relation to homology search and the calibration of its confidence for rare EC numbers examined together. |
| Research questions | (1) Do the standard EC scores count all false positives, and does it matter for published rankings? (2) How does the released CLEAN generalise to enzymes added to Swiss-Prot after its training snapshot? How much of that set is out of reach because of unseen EC numbers? How do inference rules and BLASTp compare? (3) Is CLEAN's confidence calibrated across EC-frequency tiers, and how much calibration data would certification need? |
| Hypotheses | Pre-registered, with primary endpoints (`PUBLICATION_PLAN.md` §4 and §7; the pre-registration commit c8e42f7). |
| "Materials" (data) | New-392 (development set, disclosed as such). Price-149 (confirmatory for methods frozen before 2026-09-27). SP-temporal: 2,099 Swiss-Prot enzymes created after CLEAN's snapshot plus 2,000 non-enzymes, split into CAL (839 enzymes, the only set used for tuning and calibration) and TEST (1,260 enzymes, never used for any choice). The published EC-Bench predictions of 16 predictors. |
| Independent variables | Scoring rule: support-weighted vs micro, union-label macro and per-protein. Inference rule in CLEAN's frozen space: centroid, exemplar (1-NN), kmean3, multiproto, balanced3. Method: CLEAN vs BLASTp. Sequence-identity stratum. EC-frequency tier. |
| Dependent variables | Micro-F1, union-label macro-F1, per-protein F1; the share of false positives the weighted score cannot see; rank changes (Kendall τ_b, reversed pairs); expected calibration error; selective risk; the certification budget (calibration proteins needed, years of curation). |
| Controls | CLEAN's own centroid rule. A BLASTp top hit against the same training sequences. CAL/TEST separation. Pre-registration. Paired bootstrap (B = 2000), with Holm or Benjamini–Hochberg correction. |
| Mechanism | A support-weighted average gives each label a weight equal to its number of true instances. A predicted EC number that never occurs in the test labels has zero weight, so its false positives never enter the reported precision. Support-weighted recall equals micro recall, so all of the distortion is in precision. |
| Main results (final run 532bc1f) | 67–96% of CLEAN's false positives go uncounted on New-392/Price-149. Inflation tracks predicted-EC count across EC-Bench predictors (ρ = 0.75). One significant ranking reversal (q = 0.026). 30.5% of temporal enzymes carry an EC number CLEAN never saw. balanced3 − centroid union-macro F1 +0.043 [0.030, 0.060]. BLASTp micro-F1 0.735 vs 0.641 for CLEAN, a tie below 30% identity. ECE 0.32 in the rare/singleton tiers vs 0.19 in the common tier. Certification at 20% risk would need at least four times the available calibration proteins, at least six years of curation. |
| Application | Benchmark design and reporting (an EC evaluation card), annotation pipelines, curation priorities. |
| Type of science | Computational biology and bioinformatics: evaluation methodology for enzyme-function prediction, with a statistical-learning component. The chemistry link is that EC numbers classify reactions. |
| Geographic relevance | None scientifically (section 0.3, item 9). |
| Claimed novelty | C1–C4 (`PUBLICATION_PLAN.md` §0); assessed in section G. |
| Missing information | The Introduction is not written yet. Two things need the CLEAN paper's full text: its BLASTp configuration and the provenance of Price-149. The journal's length limits could not be retrieved (the Guide for Authors blocks automated access). |

---

## B. Scientific field map

| # | Area | What it contributes | References used | Link to our contributions |
|---|---|---|---|---|
| 1 | EC system and annotation quality | EC numbers classify reactions. Unrelated proteins can share an EC number. Labels contain errors. | McDonald2023, Omelchenko2010, Schnoes2009, Rembeza2021 | Defines the labels. NISE analysis. Limits what "false positive" means. |
| 2 | Homology transfer | The default baseline and its identity dependence | Altschul1997, Camacho2009, ClaudelRenard2003, Rost2002, Tian2003 | C3: the BLASTp baseline and identity strata |
| 3 | Learned EC predictors | Model families from per-EC classifiers to contrastive PLM embeddings | Dalkiran2018, Ryu2019, Sanderson2023, Gligorijevic2021, Kim2023, Dhibar2025, Rives2021, Yu2023, YangY2024, YangZ2025, Cao2026, Ma2026, Dumontet2026 | The object of study. C1b: code audit of CLEAN, CLEAN-Contact, ProtDETR, EnzymeHunter. |
| 4 | Independent evaluations | Where published claims fail | deCrecyLagard2025, Sartori2026, Capela2025, Sathyamoorthy2026 | Motivation, and the context for C3 |
| 5 | Evaluation of function prediction | Prospective (time-delayed) evaluation; metric pathologies; incomplete ground truth | Radivojac2013, Jiang2016, Zhou2019, Phan2026, Plyusnin2019, Jiang2014 | Framing for C1 and C2 |
| 6 | Multi-label metrics | Definitions, invariance, formula variants | Zhang2014, Sokolova2009, Opitz2019 | The C1 proof and its framing |
| 7 | EC benchmarks and leakage | Tasks, identity splits, temporal hold-outs, reporting guidance | YangJ2024, Davoudi2026, Dalkiran2018, Walsh2021, Kapoor2023 | C1c (EC-Bench re-scoring), C2 |
| 8 | Long tail and prototypes | Inference rules in embedding spaces | Snell2017, Khosla2020, Sharma2023 (+ Dumontet2026, Ma2026) | C3a |
| 9 | Calibration and guarantees | Calibration, selective prediction, conformal guarantees, exact binomial limits | Guo2017, Zadrozny2002, Geifman2017, Yu2023, Han2024, DingK2024, Boger2025, DingT2023, Clopper1934 | C4 |

---

## C. Verified papers

The review cites 50 works. The 30 below carry the argument, in the evidence-table format of your prompt (section 6). The other 20 (table C2) are foundations and methods references a reviewer expects to see. Their roles are narrower, and several can go if space is tight (section M, redundancy check).

Abbreviations in the "Where" column:
- 2.1-P1 means subsection 2.1, paragraph 1.
- Types: **PR** primary research, **Bench** benchmark or resource, **Pre** preprint (not peer reviewed), **Conf** peer-reviewed conference paper, **Rev** review, **Note** arXiv note.

### C1. Core papers (30)

| # | Citation | Year | DOI / ID | Type | What the study did | Key finding (as verified) | Why it matters to this study | Where | Ev. |
|---|---|---|---|---|---|---|---|---|---|
| 1 | Yu et al., *Science* 379:1358 (CLEAN) | 2023 | 10.1126/science.adf2465 | PR | Contrastive learning on ESM-1b embeddings for EC assignment | Reported better accuracy, reliability and sensitivity than BLASTp; annotates understudied and promiscuous enzymes | The model we audit. Its public evaluation code uses support-weighted scores. Its BLASTp claim is not reproduced on our temporal set. | 2.1-P2, 2.3, 2.4 | V2 |
| 2 | Tian & Skolnick, *J Mol Biol* 333:863 | 2003 | 10.1016/j.jmb.2003.08.057 | PR | EC conservation vs pairwise identity among enzyme relatives | ~40% identity needed for >90% conservation of the first three EC digits; ~60% for all four | Sets expectations for homology transfer; motivates identity strata | 2.1-P1 | V2 |
| 3 | Rost, *J Mol Biol* 318:595 | 2002 | 10.1016/S0022-2836(02)00016-5 | PR | Conservation of EC numbers among sequence-similar enzymes | Function less conserved than anticipated | Early evidence that homology transfer errs even for similar pairs | 2.1-P1 | V2 |
| 4 | Schnoes et al., *PLoS Comput Biol* 5:e1000605 | 2009 | 10.1371/journal.pcbi.1000605 | PR | Misannotation in 37 enzyme families across public databases | Low in Swiss-Prot; high in TrEMBL, GenBank NR and KEGG | Justifies Swiss-Prot as reference; caveat for false-positive counts | 2.1-P1 | V2 |
| 5 | Rembeza & Engqvist, *PLoS Comput Biol* 17:e1009446 | 2021 | 10.1371/journal.pcbi.1009446 | PR | High-throughput assay of 122 sequences annotated EC 1.1.3.15 | ≥78% of the class inferred misannotated; four alternative activities confirmed | Experimental proof that EC labels can be wrong at scale | 2.1-P1 | V2 |
| 6 | Omelchenko et al., *Biol Direct* 5:31 | 2010 | 10.1186/1745-6150-5-31 | PR | Systematic analysis of non-homologous isofunctional enzymes (NISE) | Unrelated proteins catalyse the same reactions; complicates homology-based pathway reconstruction | Why sequence similarity cannot capture every EC relation; basis of our NISE analysis | 2.1-P1 | V2 |
| 7 | Camacho et al., *BMC Bioinformatics* 10:421 (BLAST+) | 2009 | 10.1186/1471-2105-10-421 | PR (software) | Rewritten BLAST suite | New command-line applications, speed gains | The BLASTp implementation we use (BLAST+ 2.12.0) | 2.1-P1; Methods | V2 |
| 8 | Dalkiran et al., *BMC Bioinformatics* 19:334 (ECPred) | 2018 | 10.1186/s12859-018-2368-y | PR | Per-EC ensemble classifiers; time-separated hold-out of 30 enzymes + 30 non-enzymes | F1 0.14 at the substrate level on the temporal set; such cases are very difficult | Early temporal hold-out showing much lower performance (C2) | 2.1-P2, 2.2-P2 | V2 |
| 9 | Ryu, Kim & Lee, *PNAS* 116:13996 (DeepEC) | 2019 | 10.1073/pnas.1821905116 | PR | CNN-based EC prediction | High-quality, high-throughput EC prediction (title-level) | Representative deep predictor; one of EC-Bench's methods | 2.1-P2 | V2 |
| 10 | Kim et al., *Nat Commun* 14:7370 (DeepECtransformer) | 2023 | 10.1038/s41467-023-43216-z | PR | Transformer EC predictor applied to *E. coli* K-12 | EC numbers for 464 unannotated genes; 3 validated | Discovery use case; dissected by #16; its EC-Bench variant is in our significant ranking reversal | 2.1-P2, P3 | V2 |
| 11 | Yang Y. et al., *Commun Biol* 7:1690 (CLEAN-Contact) | 2024 | 10.1038/s42003-024-07359-z | PR | CLEAN plus predicted contact maps | Reported gains over CLEAN | Successor whose public code reuses the support-weighted evaluation (C1b) | 2.1-P2 | V2 |
| 12 | Cao et al., *Patterns* 101567 (EnzymeHunter) | 2026 | 10.1016/j.patter.2026.101567 | PR | Hierarchy-aware contrastive learning with sequence + structure | Fine-grained fourth-level prediction; robust at low homology (their report) | Successor; public code uses the support-weighted metric (C1b) | 2.1-P2 | V2 |
| 13 | Ma et al., *Bioinformatics* 42(S1):btag215 (EnzPlacer) | 2026 | 10.1093/bioinformatics/btag215 | PR | Contrastive learning to predict unseen EC functions | Addresses EC numbers absent from training (title-level) | Directly addresses the unseen-EC share we quantify (C2) | 2.1-P2, 2.3 | V2 |
| 14 | Dumontet et al., *Nat Commun* 17:1146 (HIT-EC) | 2026 | 10.1038/s41467-026-68727-3 | PR | Four-level hierarchical interpretable transformer | Targets under-represented EC numbers and incomplete annotations | Shows the rare-EC problem is recognised; context for C3/C4 | 2.1-P2, 2.3 | V2 (authors V1) |
| 15 | Capela et al., *BMC Bioinformatics* 26:68 | 2025 | 10.1186/s12859-025-06081-9 | PR | ESM2, ESM1b and ProtBERT classifiers vs BLASTp for EC prediction | BLASTp marginally better overall; complementary; LLMs useful below 25% identity; ESM2 best LLM | Peer-reviewed precedent for our BLASTp result; lowers C3b to Type 3 | 2.1-P3 | V2 |
| 16 | de Crécy-Lagard et al., *G3* 15(10):jkaf169 | 2025 | 10.1093/g3journal/jkaf169 | PR | Expert dissection of DeepECtransformer predictions for >450 *E. coli* unknowns | ML largely fails to predict functions unseen in training; basic logic errors | Independent evidence that benchmark scores overstate novel-function performance | 2.1-P3 | V2 |
| 17 | Sartori et al., bioRxiv (EnzymARC) | 2026 | 10.64898/2026.08.21.746242 | Pre | Decoys with active sites disrupted at 5–15 Å radii; DIAMOND, CLEAN, DeepEC tested | >90% false-positive rate for DIAMOND and CLEAN on low-perturbation decoys | Both homology and CLEAN over-assign EC numbers; supports the open-set concern | 2.1-P2, P3 | V1 |
| 18 | Sathyamoorthy & Puri, bioRxiv | 2026 | 10.64898/2026.03.31.715487 | Pre | 3 PLMs × 9 architectures vs a train-set-matched BLASTp, by identity threshold | Parity in distribution; PLMs +31.8 pp for distant eukaryotes | **Partly contradicts our tie below 30% identity; must be discussed** | 2.1-P3 | V1 |
| 19 | Radivojac et al., *Nat Methods* 10:221 (CAFA1) | 2013 | 10.1038/nmeth.2340 | PR (community assessment) | First CAFA: 54 methods, 866 targets, time-delayed | Established prospective evaluation of function prediction | Precedent for temporal evaluation (C2) | 2.2-P1 | V2 |
| 20 | Phan et al., *Bioinform Adv* vbag221 (LAFA) | 2026 | 10.1093/bioadv/vbag221 | Bench | Continuous longitudinal evaluation server | Up-to-date comparison under evolving ground truth | Partial overlap with C2 (temporal evaluation) | 2.2-P1 | V2 |
| 21 | Plyusnin, Holm & Törönen, *PLoS Comput Biol* 15:e1007419 | 2019 | 10.1371/journal.pcbi.1007419 | PR (methods) | Compared 37 GO metrics using artificial dilution series and false-positive sets | Some metrics give erroneously high scores to false-positive sets | **Closest precedent for C1**; makes C1 Type 2 | 2.2-P1 | V2 |
| 22 | Jiang et al., *Bioinformatics* 30:i609 | 2014 | 10.1093/bioinformatics/btu472 | PR (methods) | Impact of incomplete annotations on evaluation | Evaluation is biased when ground truth is incomplete (title-level) | Caveat: uncounted "false positives" may be true functions | 2.2-P1 | V2 |
| 23 | Sokolova & Lapalme, *Inf Process Manag* 45:427 | 2009 | 10.1016/j.ipm.2009.03.002 | PR (methods) | Systematic analysis of classification measures and their invariance properties | Measures differ in which confusion-matrix changes they ignore | Formal frame for the C1 invariance result | 2.2-P2 | V2 |
| 24 | Opitz & Burst, arXiv:1911.03347 | 2019 | arXiv:1911.03347 | Note | Compared the two macro-F1 formulas in use | Can differ by up to 0.5 and reorder classifiers | Averaging choices change rankings; supports C1c | 2.2-P2 | V2 |
| 25 | Davoudi et al., *Bioinform Adv* 6:vbag004 (EC-Bench) | 2026 | 10.1093/bioadv/vbag004 | Bench | 10 methods; exact prediction, completion and recommendation; accuracy and efficiency metrics | Unified protocol; support-weighted F1 adopted for class imbalance | The benchmark we re-score (C1c); partial overlap with C2 | 2.1-P2, 2.2-P2, 2.3 | V2 |
| 26 | Yang J. et al., NeurIPS 2024 D&B (CARE) | 2024 | arXiv:2406.15669 | Bench / Conf | EC classification with identity-based OOD splits (<30%, 30–50%); reaction→EC retrieval | Standardised OOD evaluation for EC prediction | Leakage control by identity; contrast to our temporal split | 2.2-P2 | V2 |
| 27 | Ding K., Luo & Luo, *PLoS Comput Biol* 20:e1012135 (CPEC) | 2024 | 10.1371/journal.pcbi.1012135 | PR | Conformal prediction for EC annotation with limited false positives | FDR-controlled EC annotation | Closest EC-specific guarantee work; C4 differs (per-tier budget) | 2.4 | V2 |
| 28 | Boger et al., *Nat Commun* 16:85 | 2025 | 10.1038/s41467-024-55676-y | PR | Conformal guarantees for protein search, including CLEAN embeddings | User-specified risk guarantees; strong enzyme classification without new training | Guarantees on CLEAN exist but are marginal, not per tier (C4 partial overlap) | 2.4 | V2 |
| 29 | Ding T. et al., NeurIPS 2023 | 2023 | arXiv:2306.09335 | Conf | Class-conditional conformal prediction with many classes (clustered conformal) | Standard methods fail with little labelled data per class | General theory behind the rare-tier budget (C4b) | 2.4 | V2 |
| 30 | Guo et al., ICML 2017 | 2017 | PMLR 70:1321–1330 | Conf | Calibration of modern neural networks | Modern networks are miscalibrated; temperature scaling | Motivates checking CLEAN's calibration (C4a) | 2.4 | V2 |

### C2. Supporting references also cited (20)

| Key | Reference | DOI / ID | Type | Role |
|---|---|---|---|---|
| McDonald2023 | McDonald & Tipton, *FEBS J* 290:2214 | 10.1111/febs.16274 | Rev | Defines the EC system |
| Altschul1997 | Altschul et al., *Nucleic Acids Res* 25:3389 | 10.1093/nar/25.17.3389 | PR | BLAST algorithm (gapped BLAST) |
| ClaudelRenard2003 | Claudel-Renard et al., *Nucleic Acids Res* 31:6633 (PRIAM) | 10.1093/nar/gkg847 | PR | Profile-based EC transfer; an EC-Bench method |
| Gligorijevic2021 | Gligorijević et al., *Nat Commun* 12:3168 (DeepFRI) | 10.1038/s41467-021-23303-9 | PR | Structure-aware learned predictor (a reviewer will expect it) |
| Dhibar2025 | Dhibar, Basak & Jana, *Chem Sci* 16:18438 (SOLVE) | 10.1039/D5SC04513D | PR | Interpretable non-deep predictor; published in a chemistry journal |
| Rives2021 | Rives et al., *PNAS* 118:e2016239118 (ESM-1b) | 10.1073/pnas.2016239118 | PR | CLEAN's input embeddings |
| YangZ2025 | Yang Z. et al., ICLR 2025 workshop (ProtDETR) | arXiv:2501.05644 | Conf (workshop) | Successor in the code audit (C1b) |
| Jiang2016 | Jiang et al., *Genome Biol* 17:184 (CAFA2) | 10.1186/s13059-016-1037-6 | PR | CAFA series |
| Zhou2019 | Zhou et al., *Genome Biol* 20:244 (CAFA3) | 10.1186/s13059-019-1835-8 | PR | CAFA series |
| Zhang2014 | Zhang & Zhou, *IEEE TKDE* 26:1819 | 10.1109/TKDE.2013.39 | Rev | Multi-label metric definitions |
| Walsh2021 | Walsh et al., *Nat Methods* 18:1122 (DOME) | 10.1038/s41592-021-01205-4 | Guidelines | ML validation in biology |
| Kapoor2023 | Kapoor & Narayanan, *Patterns* 4:100804 | 10.1016/j.patter.2023.100804 | PR (meta-research) | Leakage in ML-based science |
| Snell2017 | Snell et al., NIPS 2017 | NIPS 30 | Conf | Prototype (class-mean) classifier = CLEAN's centroid rule |
| Khosla2020 | Khosla et al., NeurIPS 2020 | NeurIPS 33 | Conf | Supervised contrastive learning (CLEAN's objective family) |
| Sharma2023 | Sharma et al., IJCAI-23 | 10.24963/ijcai.2023/151 | Conf | Prototype classifiers for long-tailed recognition (C3a) |
| Zadrozny2002 | Zadrozny & Elkan, KDD '02 | 10.1145/775047.775151 | Conf | Isotonic calibration (our recalibration method) |
| Geifman2017 | Geifman & El-Yaniv, NIPS 2017 | NIPS 30, pp. 4878–4887 | Conf | Selective classification with guaranteed risk (SGR) |
| Han2024 | Han et al., *Brief Bioinform* 25:bbad401 (ECPICK) | 10.1093/bib/bbad401 | PR | Evidential uncertainty for EC prediction |
| Clopper1934 | Clopper & Pearson, *Biometrika* 26:404 | 10.1093/biomet/26.4.404 | PR | Exact binomial limits (certification budget) |
| Sanderson2023 | Sanderson et al., *eLife* 12:e80942 (ProteInfer) | 10.7554/eLife.80942 | PR | Representative deep predictor |

The C1 table has 30 rows and C2 has 20, for 50 in total.

---

## D. Evidence matrix

The columns are your prompt's section 7 matrix, adapted to a computational study.

| Study | Data / setting | What was varied | Method | Main finding | Limitation *for our question* | What remains unanswered |
|---|---|---|---|---|---|---|
| Tian & Skolnick 2003 | Enzyme relatives with EC annotation | Pairwise identity | Pairwise comparison of EC digits | 40% (3 digits) / 60% (4 digits) for >90% conservation | Pre-dates learned predictors; no test-time scores | How homology compares with learned embeddings on recent proteins |
| Yu et al. 2023 (CLEAN) | New-392, Price-149 and others | Method | Contrastive learning vs baselines incl. BLASTp | CLEAN better than BLASTp (reported) | Evaluation code uses support-weighted averaging (our audit); BLASTp configuration unverified | Whether the advantage holds on a leakage-controlled temporal set with complete false-positive counting |
| Capela et al. 2025 | EC data set (split design not checked; full text) | Language model; identity to reference | PLM embeddings + classifiers vs BLASTp | BLASTp marginally better; complementary below 25% identity | Not the released CLEAN; not a temporal hold-out as far as the abstract shows | Whether the pattern holds for CLEAN on post-snapshot enzymes |
| Sathyamoorthy & Puri 2026 (preprint) | 1,296 trained models; identity thresholds | PLM, architecture, identity, EC level | MLPs etc. vs train-set-matched BLASTp | Parity in distribution; PLMs far ahead for distant eukaryotes | Trained classifiers, not a frozen contrastive model; preprint | Why a frozen CLEAN shows only a tie at low identity (model, data or taxonomy?) |
| de Crécy-Lagard et al. 2025 | >450 *E. coli* unknowns, DeepECtransformer | — | Expert curation of predictions | Mostly no correct novel predictions; logic errors | One model, one organism | How often scores on benchmarks hide such errors |
| Sartori et al. 2026 (preprint) | Active-site-disrupted decoys | Perturbation radius | DIAMOND, CLEAN, DeepEC | >90% FP rate on low-perturbation decoys | Decoys are synthetic; preprint | How many real false positives current scores fail to count |
| Plyusnin et al. 2019 | GO annotations | Metric (37 of them) | Artificial dilution series + false-positive sets | Some metrics reward false predictions | GO, not EC; does not examine support-weighted flat averaging in EC code | Whether EC scoring has an analogous blind spot, and how large it is |
| Opitz & Burst 2019 | Theory + examples | Macro-F1 formula | Analytical | Up to 0.5 difference; rankings change | Macro-F1 only; not EC | Effect of the *label set* over which averaging runs |
| Davoudi et al. 2026 (EC-Bench) | 10 methods (16 predictors in its release) | Method, task | Uniform retraining and evaluation | Unified comparison; support-weighted F1 adopted | Support-weighted scores cannot count out-of-label predictions (our C1) | How rankings change under micro / union-macro / per-protein F1 (our E5) |
| Yang J. et al. 2024 (CARE) | EC classification + retrieval | Identity split, promiscuity, reaction split | Benchmark suite | Standardised OOD tasks | Identity splits, not temporal; top-1 accuracy (plan §14 [code]) | Temporal generalisation of released models; unseen-EC share |
| Dalkiran et al. 2018 (ECPred) | 30 + 30 time-separated proteins | Method | Temporal hold-out | F1 0.14 at substrate level | Very small set | A large temporal set with unseen-EC decomposition |
| Radivojac 2013 / Phan 2026 (CAFA / LAFA) | GO/HPO targets | Time | Prospective / continuous evaluation | Standard for time-delayed evaluation | GO-centric; EC is not the focus | EC-specific temporal benchmark for released models; model vintage |
| Ding K. et al. 2024 (CPEC) | EC annotation | Risk level | Conformal FDR control | Limited false positives with guarantees | Marginal (average) guarantees | Per-tier guarantees for rare EC numbers and their data cost |
| Boger et al. 2025 | Protein search incl. CLEAN | Loss / risk level | Conformal risk control | Guarantees for any search model | Marginal guarantees | Calibration of CLEAN's own confidence by tier |
| Ding T. et al. 2023 | General ML benchmarks | Classes, data per class | Clustered conformal | Class-conditional coverage fails with few examples per class | Not EC-specific | How far EC tiers are from the needed sample size |
| Dumontet et al. 2026 / Ma et al. 2026 | EC prediction | Rare / unseen EC numbers | New architectures | Recognise rare and unseen ECs as the hard cases | Method papers; do not quantify how much a released model cannot reach | The unseen share of a real temporal set and its effect on headline scores |

**Patterns.**
1. Homology transfer is well characterised, but its comparison with learned embeddings gives conflicting results, which depend on the identity profile of the test set and on how the baseline is configured.
2. Metric pathologies are documented for GO evaluation but not for flat EC scoring.
3. Temporal evaluation is standard in CAFA and appeared early in EC prediction (ECPred), but on small sets and without an unseen-EC decomposition.
4. Guarantees for EC predictions exist, but they are marginal; the per-class data cost has not been measured.

---

## E. Research gaps identified

| # | Candidate gap | Verdict | Evidence |
|---|---|---|---|
| 1 | Support-weighted EC scores cannot count predictions of EC numbers absent from the test labels | **Real; addressed (C1a)** | No EC-specific report found; GO analogue in Plyusnin 2019 |
| 2 | How common this scoring is in public EC code | **Real; addressed (C1b)** | Our audit of five code bases; no prior report found |
| 3 | Whether published EC rankings change under counting scores | **Real; addressed (C1c)** | Opitz & Burst 2019 show the general phenomenon; EC-Bench not re-scored before |
| 4 | Temporal evaluation of the *released* CLEAN with an unseen-EC decomposition | **Partially addressed by others; our addition is the decomposition and vintage matching (C2)** | CAFA, ECPred, EC-Bench (retrained), LAFA; CLEAN's own New-392 is temporal |
| 5 | Training-vintage overlap between successor models and test sets | **Real; addressed (C2 finding)** | Not found elsewhere |
| 6 | Centroid vs neighbour inference in one frozen EC space | **Real but modest; addressed (C3a)** | Known in general ML (Snell 2017; Sharma 2023); not tested for EC |
| 7 | Whether PLM-based predictors beat homology search | **Not a gap: contested and actively studied** (Capela 2025; Sathyamoorthy 2026). We add evidence (C3b). | — |
| 8 | Calibration of CLEAN's released confidence by EC-frequency tier | **Real; addressed (C4a)** | CPEC, Boger 2025, ECPICK, HIT-EC address guarantees or trustworthiness, not tier-resolved calibration of CLEAN's confidence |
| 9 | Calibration data needed to certify rare-EC predictions | **Real; addressed (C4b)** | General theory (Clopper 1934; Ding T. 2023); not quantified for EC |
| 10 | Rejecting non-functional look-alikes (open-set failure) | Real, **but addressed by EnzymARC (Sartori 2026), not us** | Cite; do not claim |
| 11 | Experimental validation of predictions | Real; **not addressed by us** (limitation) | de Crécy-Lagard 2025; Kim 2023 |
| 12 | Separating true errors from unannotated true functions | Real; **only partly addressed** (label drift in E1) | Jiang 2014 |

**Gaps that must *not* be claimed:**
- "No independent evaluation of CLEAN exists." False: EC-Bench, EnzymARC, Boger 2025.
- "No temporal benchmarks exist." False: CAFA, ECPred, LAFA, CLEAN's New-392.
- "Calibration has not been studied for EC prediction." False: CPEC, ECPICK, Boger 2025.
- "BLASTp beating PLMs is new." False: Capela 2025.

---

## F. The actual research gap

The most defensible gap has two parts:
- **The scores.** EC-number predictors are ranked with support-weighted scores that, by construction, never count a prediction of an EC number absent from the test labels. No study has quantified what this does to reported precision or to published rankings, even though the same scoring appears in the public code of CLEAN, three successors and the EC-Bench benchmark.
- **The released model.** For the released CLEAN, four things have not been examined together on a leakage-controlled temporal set with pre-registered endpoints: the share of new enzymes that are out of reach (unseen EC numbers), whether its centroid inference rule suits rare classes, how it compares with a train-set-matched BLASTp search, and whether its confidence could be certified for rare EC numbers.

**Suggested wording for the Introduction:**

> "Previous work has established that metric choice can change the ranking of function predictors (Plyusnin et al., 2019; Opitz and Burst, 2019) and that time-separated test sets are much harder than random splits (Dalkiran et al., 2018). However, the scoring used for EC predictors has not been examined for errors it cannot count, and the released CLEAN has not been evaluated on a leakage-controlled temporal benchmark that separates unseen EC numbers, compares inference rules and a homology baseline, and tests the calibration of its confidence for rare EC numbers."

---

## G. Novelty assessment

These are the types defined in your prompt:
- **Type 1:** completely unexplored.
- **Type 2:** incremental but meaningful.
- **Type 3:** contextual or application novelty.

| Claim | Closest prior work | Type | What can be claimed | What should NOT be claimed |
|---|---|---|---|---|
| C1a: invariance of support-weighted and true-label-macro EC scores to out-of-label predictions (proof) | Plyusnin 2019 (GO); Opitz & Burst 2019; Sokolova & Lapalme 2009 | **2** | "We formalise, for EC prediction, a blind spot of support-weighted averaging" | "A new metric theory"; "first to show that metrics can mislead" |
| C1b: prevalence in five public code bases | None found | **2** (a new, concrete variable: evaluation code) | "To our knowledge, not previously reported" | "Field-wide" or "all EC studies" (GloEC and SOLVE use macro averaging; CPEC uses sample averaging; plan §14 [code]) |
| C1c: EC-Bench re-scoring (inflation tracks predicted-EC count, ρ = 0.75; one significant reversal) | Opitz & Burst 2019 (rankings change with averaging) | **2** | "Re-scoring the published predictions changes one ranking significantly" | "Published rankings are wrong" (only 1 of 120 pairs is significant) |
| C2: leakage-controlled temporal benchmark with unseen-EC share and model-vintage matching | CAFA; ECPred 2018; LAFA 2026; EC-Bench; CLEAN's New-392 | **2**, with a Type 3 element (known design applied to the released CLEAN) | "A temporal benchmark that separates unseen EC numbers and matches training vintage" | "The first temporal EC benchmark" |
| C3a: balanced3 vs centroid in one frozen space (+0.043 union-macro F1) | Snell 2017; Sharma 2023 (prototype vs neighbour rules in general ML) | **2** | "A modest, pre-registered improvement from class-balanced neighbour inference" | "A new method that outperforms CLEAN" in general; BLASTp beats it |
| C3b: BLASTp beats the released CLEAN on a temporal set; tie below 30% identity | **Capela 2025 (same direction, peer reviewed)**; contradicts CLEAN's own report; partly contradicts Sathyamoorthy & Puri 2026 (preprint) | **3** | "Consistent with Capela et al. (2025), a train-set-matched BLASTp search remains a strong baseline on post-2022 enzymes" | Anything implying the BLASTp result is a discovery |
| C4a: tier-dependent miscalibration of CLEAN's released confidence | CPEC 2024; Boger 2025; ECPICK 2024; HIT-EC 2026 | **2** | "Calibration of CLEAN's released confidence has not been reported by EC-frequency tier" | "Nobody has studied uncertainty in EC prediction" |
| C4b: exact-binomial certification budget and curation-time conversion | Clopper & Pearson 1934; Ding T. 2023 (general) | **2** | "To our knowledge, the first EC-specific estimate of the calibration data needed to certify rare-tier predictions" | That certification is impossible (it is a lower bound under stated assumptions) |
| EC evaluation card | DOME (Walsh 2021); model info sheets (Kapoor 2023) | **3** | A practical reporting checklist for EC prediction | A new framework |

**Overall:** the study is **Type 2**: incremental but meaningful. It is an evaluation and reliability paper with one formal result, one resource and one confirmatory homology finding. No claim can be defended as Type 1.

The strongest honest framing: *"Existing EC scores, benchmarks and confidence estimates each leave a specific blind spot; we measure each one for the most widely used contrastive model and provide the reporting needed to close them."*

---

## H. Recommended structure

The Related Work section is about 1,650 words, has 50 references and sits after a short Introduction.

**2.1 From homology transfer to learned predictors** (3 paragraphs)
- **Purpose:** establish what EC numbers are, why homology transfer is the baseline and where it fails, what learned models added, and what independent tests found.
- **Papers:** McDonald2023, Omelchenko2010, Altschul1997, Camacho2009, ClaudelRenard2003, Rost2002, Tian2003, Schnoes2009, Rembeza2021, Dalkiran2018, Ryu2019, Sanderson2023, Gligorijevic2021, Kim2023, Dhibar2025, Yu2023, Rives2021, YangY2024, YangZ2025, Cao2026, Ma2026, Dumontet2026, deCrecyLagard2025, Sartori2026, Capela2025, Sathyamoorthy2026, Davoudi2026.
- **Argument:** both homology and learned methods have documented limits. Whether embeddings beat homology depends on test-set identity and baseline configuration.
- **Connects to 2.2:** if conclusions depend on the test design, the way methods are scored matters.

**2.2 Metrics and benchmarks** (2 paragraphs)
- **Purpose:** show what the function-prediction community learned about evaluation, then how EC predictors are actually scored.
- **Papers:** Radivojac2013, Jiang2016, Zhou2019, Phan2026, Plyusnin2019, Jiang2014, Zhang2014, Sokolova2009, Opitz2019, Davoudi2026, YangJ2024, Dalkiran2018, Walsh2021, Kapoor2023.
- **Argument:** metric choice can decide rankings. EC studies commonly use support-weighted averaging, and whether it counts every false positive has not been asked (gap 1).
- **Connects to 2.3:** temporal sets bring unseen EC numbers, which scores and models handle badly.

**2.3 Rare and unseen EC numbers** (1 paragraph)
- **Purpose:** explain why unseen and rare classes cap and bias performance, and why the inference rule matters.
- **Papers:** Davoudi2026, Dumontet2026, Ma2026, Snell2017, Khosla2020, Sharma2023.
- **Argument:** a closed-set model cannot recall unseen labels. Centroid vs neighbour rules for rare classes are untested in EC space (gaps 4 and 6).
- **Connects to 2.4:** rare classes are also where confidence is least reliable.

**2.4 Confidence, calibration and guarantees** (1 paragraph)
- **Purpose:** review calibration, selective prediction and conformal guarantees, including those for EC and CLEAN.
- **Papers:** Guo2017, Zadrozny2002, Geifman2017, Yu2023, Han2024, DingK2024, Boger2025, DingT2023, Clopper1934.
- **Argument:** guarantees exist but are marginal. Per-tier certification needs data whose amount has not been measured (gaps 8 and 9).
- **Connects to 2.5:** leads to the synthesis.

**2.5 Open questions** (1 paragraph)
- **Purpose:** state the three open questions that the Results answer, in the order of the Results.

**Why this length.** About 1,650 words is roughly three pages in the single-column cas-sc layout; that page count is my estimate, since the section was not compiled. The Guide for Authors' length limit could not be retrieved, so check it.

The plan's main-text budget is about 6,500 words, with an 800-word Introduction. With this section, keep the Introduction to **about 500–600 words**, covering motivation, CLEAN, the three questions in one sentence each, and the contributions. The total then stays near 7,000 words.

If the editor asks for cuts, work in this order:
1. Merge 2.3 into 2.4.
2. Drop the optional references (section N).
3. Shorten 2.1-P2 to one sentence per model family.

---

## I. Paper-by-paper placement

The table follows order of appearance. "Why here" explains why the paper comes at this point and not earlier or later.

| Order | Key | Location | Supports | Why here |
|---|---|---|---|---|
| 1 | McDonald2023 | 2.1-P1 s1 | EC = four-level reaction classification | Must define the object before anything else |
| 2 | Omelchenko2010 | 2.1-P1 s2 | Unrelated proteins share EC numbers | Qualifies the definition before homology transfer is introduced |
| 3–5 | Altschul1997, Camacho2009, ClaudelRenard2003 | 2.1-P1 s3 | Homology/profile transfer | The baseline comes before its limits |
| 6–7 | Rost2002, Tian2003 | 2.1-P1 s4–5 | Transfer reliability vs identity | Limits of the baseline; sets up identity strata |
| 8–9 | Schnoes2009, Rembeza2021 | 2.1-P1 s6–7 | Label errors | Last point of P1: the labels themselves are imperfect (needed again in 2.2) |
| 10–13 | Dalkiran2018, Ryu2019, Sanderson2023, Gligorijevic2021 | 2.1-P2 s1 | Learned model families | Learned models are introduced as the answer to homology's limits |
| 14–15 | Kim2023, Dhibar2025 | 2.1-P2 s2 | Transformer discovery use; simple ensembles compete | Before CLEAN, to show both the ambition and the strength of simple methods |
| 16–17 | Yu2023, Rives2021 | 2.1-P2 s3–4 | CLEAN and its input | The model under study, placed at the peak of the development |
| 18 | Davoudi2026, Sartori2026 | 2.1-P2 s4 | CLEAN as a reference method | Establishes why auditing CLEAN matters |
| 19–23 | YangY2024, YangZ2025, Cao2026, Ma2026, Dumontet2026 | 2.1-P2 s5 | Successors | After CLEAN, because they build on or compete with it (and three share its metric code) |
| 24–25 | deCrecyLagard2025, Sartori2026 | 2.1-P3 s1–2 | Independent failures | After the claims, to contrast them |
| 26–27 | Capela2025, Sathyamoorthy2026 | 2.1-P3 s3–5 | PLM vs BLAST disagreement | Ends 2.1 with the unresolved question that our BLASTp analysis addresses |
| 28–31 | Radivojac2013, Jiang2016, Zhou2019, Phan2026 | 2.2-P1 s1–2 | Prospective evaluation | Opens the evaluation subsection with the established gold standard |
| 32–33 | Plyusnin2019, Jiang2014 | 2.2-P1 s3–4 | Metric pathologies; incomplete truth | Lessons from GO that carry over to EC, placed before EC practice |
| 34–36 | Zhang2014, Sokolova2009, Opitz2019 | 2.2-P2 s1–3 | Averaging definitions and invariance | Technical base immediately before stating EC practice |
| 37 | Davoudi2026 | 2.2-P2 s4 | EC-Bench uses support-weighted F1 | The concrete EC practice that gap 1 targets |
| 38 | YangJ2024 (+ Davoudi2026) | 2.2-P2 s5–6 | Benchmarks and leakage control | Credit the benchmarks before stating what they leave open |
| 39 | Dalkiran2018 | 2.2-P2 s7 | Early temporal hold-out | Precedent for C2, placed with the benchmarks |
| 40–41 | Walsh2021, Kapoor2023 | 2.2-P2 s8 | General validation and leakage guidance | Wider frame just before the gap statement |
| 42–44 | Davoudi2026, Dumontet2026, Ma2026 | 2.3 s1–2 | Imbalance; rare and unseen ECs | Opens the rare-EC argument with the methods that target it |
| 45–47 | Snell2017, Khosla2020, Sharma2023 | 2.3 s6–8 | CLEAN's centroid = prototype classifier; long-tail bias | Explains the inference-rule question after the recall ceiling |
| 48–50 | Guo2017, Zadrozny2002, Geifman2017 | 2.4 s2–3 | Calibration and selective prediction | General tools first |
| 51–52 | Yu2023, Han2024 | 2.4 s4 | EC confidence methods | EC-specific confidence next |
| 53–54 | DingK2024, Boger2025 | 2.4 s5–6 | Conformal guarantees for EC and CLEAN | Strongest existing guarantees, just before their limitation |
| 55–56 | DingT2023, Clopper1934 | 2.4 s8–9 | Per-class guarantees need data; exact limits | Leads directly to the certification-budget gap |

"s" means sentence.

---

## J. Paragraph-by-paragraph blueprint

**P0: Opening (48 words)**
- **Purpose:** announce the four strands.
- **Main argument:** how predictors should be compared depends on homology, metrics, rare classes and confidence.
- **Papers:** none.
- **Transition:** into the first strand.

**P1: Homology transfer and its limits**
- **Purpose:** define EC numbers and establish the baseline and its failure modes.
- **Main argument:** homology transfer is long-standing but identity-dependent, and it propagates errors. The labels themselves are imperfect.
- **Papers:** McDonald2023, Omelchenko2010, Altschul1997, Camacho2009, ClaudelRenard2003, Rost2002, Tian2003, Schnoes2009, Rembeza2021.
- **Evidence:**
  - the 40%/60% identity thresholds;
  - Swiss-Prot has low error while automatic databases have high error;
  - at least 78% misannotation in EC 1.1.3.15.
- **Transition:** learned models were built to go beyond close homologues.

**P2: Learned predictors and CLEAN**
- **Purpose:** trace the model families up to CLEAN and its successors.
- **Main argument:** the field moved to contrastive PLM embeddings. CLEAN became a reference method and spawned successors.
- **Papers:** Dalkiran2018, Ryu2019, Sanderson2023, Gligorijevic2021, Kim2023, Dhibar2025, Yu2023, Rives2021, Davoudi2026, Sartori2026, YangY2024, YangZ2025, Cao2026, Ma2026, Dumontet2026.
- **Evidence:**
  - 464 genes proposed and 3 validated (DeepECtransformer);
  - CLEAN reported to be better than BLASTp;
  - CLEAN used in EC-Bench and EnzymARC.
- **Transition:** original reports vs independent tests.

**P3: Independent evaluations and the homology debate**
- **Purpose:** show the counter-evidence and the unresolved embedding-vs-homology question.
- **Main argument:** independent tests find failures on novel functions and on decoys. The PLM-vs-BLAST result depends on the test set's identity and on the baseline configuration.
- **Papers:** deCrecyLagard2025, Sartori2026, Capela2025, Sathyamoorthy2026.
- **Evidence:**
  - more than 450 unknowns, with mostly no novel correct predictions;
  - a false-positive rate above 90% on decoys;
  - BLASTp marginally better overall, with PLMs useful below 25% identity;
  - parity in distribution, with a PLM advantage for distant eukaryotes.
- **Transition:** our train-set-matched BLASTp by identity stratum. This leads into 2.2, because conclusions depend on evaluation design.

**P4: Lessons from prospective function-prediction assessment**
- **Purpose:** import CAFA's lessons.
- **Main argument:** time-delayed evaluation is the gold standard. Metrics can decide rankings, and incomplete ground truth blurs false positives.
- **Papers:** Radivojac2013, Jiang2016, Zhou2019, Phan2026, Plyusnin2019, Jiang2014.
- **Evidence:** 37 metrics; erroneously high scores for false-positive sets.
- **Transition:** EC predictors are scored differently.

**P5: How EC predictors are scored and benchmarked**
- **Purpose:** state EC practice and credit the benchmarks.
- **Main argument:** averaging choices are not interchangeable. EC studies commonly use support weighting; benchmarks now control leakage and tasks; temporal hold-outs are hard. Nobody has asked whether support weighting counts every false positive.
- **Papers:** Zhang2014, Sokolova2009, Opitz2019, Davoudi2026, YangJ2024, Dalkiran2018, Walsh2021, Kapoor2023.
- **Evidence:**
  - macro-F1 formulas differ by up to 0.5 (the paper's own claim; not repeated in the text);
  - EC-Bench adopts weighted F1 for imbalance;
  - CARE uses <30% and 30–50% identity splits;
  - ECPred reaches F1 0.14 on its temporal set.
- **Transition:** the gap statement leads to C1. Temporal sets bring unseen classes (2.3).

**P6: Rare and unseen EC numbers**
- **Purpose:** explain the recall ceiling and the inference-rule question.
- **Main argument:** unseen labels cap the recall of closed-set models. CLEAN's centroid rule is a prototype classifier, and its suitability for rare classes is untested in EC space.
- **Papers:** Davoudi2026, Dumontet2026, Ma2026, Snell2017, Khosla2020, Sharma2023.
- **Evidence:** HIT-EC's stated motivation; EnzPlacer's unseen-EC focus; prototype classifiers in long-tailed recognition.
- **Transition:** rare classes are also where confidence matters most.

**P7: Confidence, calibration and guarantees**
- **Purpose:** review the calibration and guarantee tools, including those for CLEAN.
- **Main argument:** calibration and conformal guarantees exist for EC prediction. They are marginal, and per-tier guarantees need calibration data whose amount is unknown.
- **Papers:** Guo2017, Zadrozny2002, Geifman2017, Yu2023, Han2024, DingK2024, Boger2025, DingT2023, Clopper1934.
- **Evidence:** the GMM confidence in CLEAN; FDR control (CPEC); risk control on CLEAN embeddings (Boger); class-conditional failure with few examples (Ding T.).
- **Transition:** into the synthesis.

**P8: Open questions (79 words)**
- **Purpose:** set up the Results.
- **Main argument:** three questions (scores, temporal reach with inference rule and homology, certification data).
- **Papers:** none.
- **Transition:** into Methods.

---

## K. Full literature-review draft

The LaTeX source is `paper/literature_review.tex`. The plain-text copy below is generated from it by `check_review.py` in the scratchpad, so the two cannot drift apart. The citations are rendered as natbib author-year with cas-model2-names; "2024a/b" is the style's automatic disambiguation of two "Yang et al., 2024" entries.

<!-- K:BEGIN -->

*1,640 words of prose in 9 paragraphs; 58 citations of 50 references.*

#### 2. Related work

Four lines of earlier work shape how enzyme function predictors should be compared: the limits of homology transfer and the learned methods that followed it, the metrics and benchmarks used to rank these methods, the handling of rare and unseen EC numbers, and the reliability of prediction confidence.

##### 2.1 From homology transfer to learned predictors

EC numbers classify enzyme-catalysed reactions in a four-level hierarchy (McDonald and Tipton, 2023). Because the scheme classifies reactions rather than proteins, one protein can carry several EC numbers, and unrelated proteins can catalyse the same reaction (Omelchenko et al., 2010). A long-standing way to annotate a new sequence is to transfer EC numbers from its closest annotated homologue, found with BLAST (Altschul et al., 1997; Camacho et al., 2009) or with EC-specific sequence profiles such as PRIAM (Claudel-Renard et al., 2003). The reliability of this transfer falls with sequence identity. Rost (2002) found that enzyme function is less conserved between homologues than had been assumed, and Tian and Skolnick (2003) estimated that about 40% identity is needed before more than 90% of homologous pairs share the first three EC digits, and about 60% before they share all four. Transfer also copies existing errors. Across 37 enzyme families, misannotation was rare in the manually reviewed Swiss-Prot but frequent in automatically annotated databases (Schnoes et al., 2009), and an experimental screen of 122 sequences from a single class, EC 1.1.3.15, implied that at least 78% of the sequences annotated to that class were misannotated (Rembeza and Engqvist, 2021).

Machine-learning predictors were developed to reach beyond close homologues. ECPred trains an independent classifier for each EC number (Dalkiran et al., 2018), DeepEC and ProteInfer apply convolutional neural networks to the sequence (Ryu et al., 2019; Sanderson et al., 2023), and DeepFRI combines language-model features with protein structure in a graph convolutional network (Gligorijević et al., 2021). DeepECtransformer, which uses transformer layers, proposed EC numbers for 464 unannotated *Escherichia coli* genes and confirmed three of them experimentally (Kim et al., 2023), while SOLVE showed that interpretable ensembles built on short sequence fragments can compete with deep models (Dhibar et al., 2025). CLEAN (Yu et al., 2023) takes embeddings from the protein language model ESM-1b (Rives et al., 2021), trains a contrastive embedding on them and assigns the EC numbers whose cluster centres lie nearest to the query. It was reported to be more accurate and sensitive than BLASTp and has since been used as a reference method in independent benchmarks (Davoudi et al., 2026; Sartori et al., 2026). Later methods introduce predicted contact maps (CLEAN-Contact; Yang et al., 2024b), residue-level detection (ProtDETR; Yang et al., 2025), hierarchy-aware contrastive losses with structural information (EnzymeHunter; Cao et al., 2026), the prediction of EC numbers absent from the training data (EnzPlacer; Ma et al., 2026) and a hierarchical transformer (HIT-EC; Dumontet et al., 2026).

Independent evaluations have been less favourable than the original reports. An expert analysis of DeepECtransformer predictions for more than 450 *E. coli* enzymes of unknown function concluded that current models mostly fail to make novel predictions and make basic logic errors that human annotators avoid (de Crécy-Lagard et al., 2025). In a 2026 preprint, CLEAN and the homology tool DIAMOND still assigned the original EC number to more than 90% of weakly perturbed decoys whose active sites had been deliberately disrupted (Sartori et al., 2026). Whether learned embeddings outperform homology search at all is also contested. A comparison of ESM2, ESM-1b and ProtBERT classifiers with BLASTp found BLASTp marginally better overall, although the two approaches were complementary: the language models, ESM2 in particular, predicted some EC numbers better and were more useful for enzymes with less than 25% identity to the reference database (Capela et al., 2025). Another 2026 preprint, which matched the BLASTp database to each model's training data, found language-model classifiers level with BLASTp on in-distribution proteins and well ahead of it on evolutionarily distant eukaryotic enzymes (Sathyamoorthy and Puri, 2026). The outcome therefore appears to depend on the identity profile of the test set and on how the homology baseline is configured. For this reason we compare CLEAN with a standard BLASTp search against CLEAN's own training sequences and report the results by identity stratum.

##### 2.2 Metrics and benchmarks

Community assessments of protein function prediction evaluate methods prospectively. In the CAFA challenges, methods predict functions for proteins that lack experimental annotations and are scored on the annotations that accumulate afterwards (Radivojac et al., 2013; Jiang et al., 2016; Zhou et al., 2019); LAFA now runs this time-delayed evaluation continuously (Phan et al., 2026). Two lessons from the evaluation of Gene Ontology predictors apply directly to EC prediction. First, the metric can decide the ranking: Plyusnin et al. (2019) compared 37 evaluation metrics for Gene Ontology classifiers and found that some of them gave erroneously high scores to deliberately constructed false-positive annotation sets. Second, the ground truth is incomplete, so a prediction scored as a false positive may be a true function that has not yet been annotated (Jiang et al., 2014).

EC predictors, by contrast, are usually compared with standard multi-label precision, recall and F1. These can be averaged over proteins, averaged over labels, or pooled over all decisions by micro-averaging (Zhang and Zhou, 2014). These choices are not interchangeable. Measures differ in which changes to the confusion matrix they are invariant to (Sokolova and Lapalme, 2009), and even the two formulas in use for macro-F1 can rank the same classifiers differently (Opitz and Burst, 2019). EC studies commonly average over EC classes with weights proportional to class support; EC-Bench, for example, adopts the support-weighted F1 to account for class imbalance (Davoudi et al., 2026). Recent benchmarks have also standardised tasks and controlled leakage. CARE tests EC classification on sequences with less than 30% or 30–50% identity to the training data and adds the retrieval of EC numbers from reactions (Yang et al., 2024a), and EC-Bench compares ten existing methods on exact EC prediction, EC completion and EC recommendation under uniform assumptions (Davoudi et al., 2026). Temporal hold-outs have a longer history in EC prediction. ECPred was tested on a small set of proteins that received an EC number only in a later Swiss-Prot release and reached an F1 of 0.14 at the substrate level; its authors noted that scores on such proteins can differ greatly from published values because these cases are very difficult (Dalkiran et al., 2018). These efforts are in line with general guidance on validation and data leakage in biological machine learning (Walsh et al., 2021; Kapoor and Narayanan, 2023). To our knowledge, however, no study has asked whether the averaging used to score EC predictors counts every false positive. We show that the support-weighted form does not.

##### 2.3 Rare and unseen EC numbers

EC labels are imbalanced (Davoudi et al., 2026), and several recent methods target the rare classes. HIT-EC was motivated by the low performance of existing models on under-represented EC numbers and by incomplete annotations (Dumontet et al., 2026), and EnzPlacer predicts EC numbers that are absent from the training set altogether (Ma et al., 2026). A closed-set predictor such as CLEAN can only output EC numbers seen in training. The test labels it has never seen therefore cap its attainable recall, and a single headline score on a time-split test set mixes proteins with unseen EC numbers and proteins whose EC numbers were seen. For the seen classes, the inference rule matters. CLEAN scores each EC number by the distance between the query and the mean embedding of that EC number's training sequences. This is the prototype classifier of few-shot learning (Snell et al., 2017), applied to an embedding trained with a supervised contrastive objective (Khosla et al., 2020). In long-tailed recognition, prototype classifiers such as the nearest class mean have been proposed to counter the bias of standard classifiers towards frequent classes (Sharma et al., 2023). To our knowledge, centroid and nearest-neighbour rules have not been compared within one frozen EC embedding space, so it is not known which of them serves rare EC numbers better.

##### 2.4 Confidence, calibration and guarantees

For a curator, a predicted EC number is most useful when its confidence can be trusted. Modern neural networks are often poorly calibrated (Guo et al., 2017), and their scores can be mapped to calibrated probabilities after training, for example by isotonic regression (Zadrozny and Elkan, 2002). Selective classification goes one step further: the model abstains when its confidence is low, and the error rate of the predictions it keeps is guaranteed with high probability (Geifman and El-Yaniv, 2017). For EC prediction, CLEAN expresses confidence by fitting a two-component Gaussian mixture to embedding distances (Yu et al., 2023), and ECPICK uses evidential deep learning (Han et al., 2024). Conformal methods add distribution-free guarantees. CPEC limits the false discovery rate of EC annotations (Ding et al., 2024), and Boger et al. (2025) applied conformal methods to protein search with user-chosen risk guarantees, including search with CLEAN embeddings. These guarantees hold on average over test proteins. A guarantee for each EC class, or for each tier of rare classes, needs calibration data from that class, and class-conditional conformal prediction becomes difficult when there are many classes with few examples each (Ding et al., 2023). Exact binomial confidence limits (Clopper and Pearson, 1934) determine the minimum number of calibration proteins that such a guarantee needs. To our knowledge, nobody has measured how far the rare EC tiers are from that minimum, or how long curation would take to close the gap.

##### 2.5 Open questions

The field therefore has strong predictors and increasingly careful benchmarks, but three questions remain open for EC prediction. Do the standard scores count all false positives? How much of a temporal test set lies beyond the reach of a closed-set model, and how does a frozen contrastive embedding compare with homology search on such a set? And how much calibration data do rare EC numbers need before their predictions can be certified? This study addresses these questions in turn.

##### References (approximate rendering of the Elsevier author-year list)

- Altschul, S. F., Madden, T. L., Schäffer, A. A., Zhang, J., Zhang, Z., Miller, W., Lipman, D. J., 1997. Gapped BLAST and PSI-BLAST: a new generation of protein database search programs. Nucleic Acids Research 25(17), 3389–3402. doi:10.1093/nar/25.17.3389.
- Boger, R. S., Chithrananda, S., Angelopoulos, A. N., Yoon, P. H., Jordan, M. I., Doudna, J. A., 2025. Functional protein mining with conformal guarantees. Nature Communications 16, 85. doi:10.1038/s41467-024-55676-y.
- Camacho, C., Coulouris, G., Avagyan, V., Ma, N., Papadopoulos, J., Bealer, K., Madden, T. L., 2009. BLAST+: architecture and applications. BMC Bioinformatics 10, 421. doi:10.1186/1471-2105-10-421.
- Cao, G., Ouyang, J., Xiong, X., Liu, C., Zhang, Y., Yang, S., Shi, T., Wu, J., 2026. EnzymeHunter: achieving fine-grained enzyme function prediction with a hierarchically aware contrastive learning framework. Patterns, 101567. doi:10.1016/j.patter.2026.101567.
- Capela, J., Zimmermann-Kogadeeva, M., van Dijk, A. D. J., et al., 2025. Comparative assessment of protein large language models for enzyme commission number prediction. BMC Bioinformatics 26, 68. doi:10.1186/s12859-025-06081-9.
- Claudel-Renard, C., Chevalet, C., Faraut, T., Kahn, D., 2003. Enzyme-specific profiles for genome annotation: PRIAM. Nucleic Acids Research 31(22), 6633–6639. doi:10.1093/nar/gkg847.
- Clopper, C. J., Pearson, E. S., 1934. The use of confidence or fiducial limits illustrated in the case of the binomial. Biometrika 26(4), 404–413. doi:10.1093/biomet/26.4.404.
- Dalkiran, A., Rifaioglu, A. S., Martin, M. J., et al., 2018. ECPred: a tool for the prediction of the enzymatic functions of protein sequences based on the EC nomenclature. BMC Bioinformatics 19, 334. doi:10.1186/s12859-018-2368-y.
- Davoudi, S., Henry, C. S., Miller, C. S., Banaei-Kashani, F., 2026. EC-Bench: a benchmark for enzyme commission number prediction. Bioinformatics Advances 6(1), vbag004. doi:10.1093/bioadv/vbag004.
- de Crécy-Lagard, V., Dias, R., Sexson, N., Friedberg, I., Yuan, Y., Swairjo, M. A., 2025. Limitations of current machine learning models in predicting enzymatic functions for uncharacterized proteins. G3 Genes|Genomes|Genetics 15(10), jkaf169. doi:10.1093/g3journal/jkaf169.
- Dhibar, S., Basak, S., Jana, B., 2025. Prediction of enzyme function using an interpretable optimized ensemble learning framework. Chemical Science 16(39), 18438–18449. doi:10.1039/D5SC04513D.
- Ding, K., Luo, J., Luo, Y., 2024. Leveraging conformal prediction to annotate enzyme function space with limited false positives. PLoS Computational Biology 20(5), e1012135. doi:10.1371/journal.pcbi.1012135.
- Ding, T., Angelopoulos, A. N., Bates, S., Jordan, M. I., Tibshirani, R. J., 2023. Class-conditional conformal prediction with many classes. In: Advances in Neural Information Processing Systems 36 (NeurIPS 2023). arXiv:2306.09335.
- Dumontet, L., Han, S.-R., Oh, T.-J., Kang, M., 2026. Trustworthy prediction of enzyme commission numbers using a hierarchical interpretable transformer. Nature Communications 17, 1146. doi:10.1038/s41467-026-68727-3.
- Geifman, Y., El-Yaniv, R., 2017. Selective classification for deep neural networks. In: Advances in Neural Information Processing Systems 30 (NIPS 2017), pp. 4878–4887.
- Gligorijević, V., Renfrew, P. D., Kosciolek, T., et al., 2021. Structure-based protein function prediction using graph convolutional networks. Nature Communications 12, 3168. doi:10.1038/s41467-021-23303-9.
- Guo, C., Pleiss, G., Sun, Y., Weinberger, K. Q., 2017. On calibration of modern neural networks. In: Proceedings of the 34th International Conference on Machine Learning, vol. 70, pp. 1321–1330.
- Han, S.-R., Park, M., Kosaraju, S., Lee, J., Lee, H., Lee, J. H., Oh, T.-J., Kang, M., 2024. Evidential deep learning for trustworthy prediction of enzyme commission number. Briefings in Bioinformatics 25(1), bbad401. doi:10.1093/bib/bbad401.
- Jiang, Y., Clark, W. T., Friedberg, I., Radivojac, P., 2014. The impact of incomplete knowledge on the evaluation of protein function prediction: a structured-output learning perspective. Bioinformatics 30(17), i609–i616. doi:10.1093/bioinformatics/btu472.
- Jiang, Y., Oron, T. R., Clark, W. T., et al., 2016. An expanded evaluation of protein function prediction methods shows an improvement in accuracy. Genome Biology 17, 184. doi:10.1186/s13059-016-1037-6.
- Kapoor, S., Narayanan, A., 2023. Leakage and the reproducibility crisis in machine-learning-based science. Patterns 4(9), 100804. doi:10.1016/j.patter.2023.100804.
- Khosla, P., Teterwak, P., Wang, C., Sarna, A., Tian, Y., Isola, P., Maschinot, A., Liu, C., Krishnan, D., 2020. Supervised contrastive learning. In: Advances in Neural Information Processing Systems 33 (NeurIPS 2020).
- Kim, G. B., Kim, J. Y., Lee, J. A., Norsigian, C. J., Palsson, B. O., Lee, S. Y., 2023. Functional annotation of enzyme-encoding genes using deep learning with transformer layers. Nature Communications 14, 7370. doi:10.1038/s41467-023-43216-z.
- Ma, X., Joshi, P., Friedberg, I., Li, Q., 2026. How not to be seen: predicting unseen enzyme functions using contrastive learning. Bioinformatics 42(Suppl. 1), btag215. doi:10.1093/bioinformatics/btag215.
- McDonald, A. G., Tipton, K. F., 2023. Enzyme nomenclature and classification: the state of the art. The FEBS Journal 290(9), 2214–2231. doi:10.1111/febs.16274.
- Omelchenko, M. V., Galperin, M. Y., Wolf, Y. I., Koonin, E. V., 2010. Non-homologous isofunctional enzymes: a systematic analysis of alternative solutions in enzyme evolution. Biology Direct 5, 31. doi:10.1186/1745-6150-5-31.
- Opitz, J., Burst, S., 2019. Macro F1 and macro F1. arXiv preprint arXiv:1911.03347. Preprint, not peer reviewed.
- Phan, A., Wang, Y., Boadu, F., Cheng, J., Radivojac, P., Friedberg, I., 2026. LAFA: a framework for reproducible longitudinal assessment of protein function annotation models. Bioinformatics Advances, vbag221. doi:10.1093/bioadv/vbag221.
- Plyusnin, I., Holm, L., Törönen, P., 2019. Novel comparison of evaluation metrics for gene ontology classifiers reveals drastic performance differences. PLoS Computational Biology 15(11), e1007419. doi:10.1371/journal.pcbi.1007419. Correction: doi:10.1371/journal.pcbi.1010249.
- Radivojac, P., et al., 2013. A large-scale evaluation of computational protein function prediction. Nature Methods 10, 221–227. doi:10.1038/nmeth.2340.
- Rembeza, E., Engqvist, M. K. M., 2021. Experimental and computational investigation of enzyme functional annotations uncovers misannotation in the EC 1.1.3.15 enzyme class. PLoS Computational Biology 17, e1009446. doi:10.1371/journal.pcbi.1009446.
- Rives, A., Meier, J., Sercu, T., Goyal, S., Lin, Z., Liu, J., Guo, D., Ott, M., Zitnick, C. L., Ma, J., Fergus, R., 2021. Biological structure and function emerge from scaling unsupervised learning to 250 million protein sequences. Proceedings of the National Academy of Sciences of the USA 118(15), e2016239118. doi:10.1073/pnas.2016239118.
- Rost, B., 2002. Enzyme function less conserved than anticipated. Journal of Molecular Biology 318, 595–608. doi:10.1016/S0022-2836(02)00016-5.
- Ryu, J. Y., Kim, H. U., Lee, S. Y., 2019. Deep learning enables high-quality and high-throughput prediction of enzyme commission numbers. Proceedings of the National Academy of Sciences of the USA 116(28), 13996–14001. doi:10.1073/pnas.1821905116.
- Sanderson, T., Bileschi, M. L., Belanger, D., Colwell, L. J., 2023. ProteInfer, deep neural networks for protein functional inference. eLife 12, e80942. doi:10.7554/eLife.80942.
- Sartori, J., et al., 2026. A novel benchmark dataset for enzyme function prediction reveals the limitations of state-of-the-art models. bioRxiv preprint. doi:10.64898/2026.08.21.746242. Preprint, not peer reviewed.
- Sathyamoorthy, R., Puri, M., 2026. Protein language models outperform BLAST for evolutionarily distant enzymes: a systematic benchmark of EC number prediction. bioRxiv preprint. doi:10.64898/2026.03.31.715487. Preprint, not peer reviewed.
- Schnoes, A. M., Brown, S. D., Dodevski, I., Babbitt, P. C., 2009. Annotation error in public databases: misannotation of molecular function in enzyme superfamilies. PLoS Computational Biology 5(12), e1000605. doi:10.1371/journal.pcbi.1000605.
- Sharma, S., Xian, Y., Yu, N., Singh, A., 2023. Learning prototype classifiers for long-tailed recognition. In: Proceedings of the Thirty-Second International Joint Conference on Artificial Intelligence (IJCAI-23). doi:10.24963/ijcai.2023/151.
- Snell, J., Swersky, K., Zemel, R., 2017. Prototypical networks for few-shot learning. In: Advances in Neural Information Processing Systems 30 (NIPS 2017).
- Sokolova, M., Lapalme, G., 2009. A systematic analysis of performance measures for classification tasks. Information Processing & Management 45(4), 427–437. doi:10.1016/j.ipm.2009.03.002.
- Tian, W., Skolnick, J., 2003. How well is enzyme function conserved as a function of pairwise sequence identity? Journal of Molecular Biology 333, 863–882. doi:10.1016/j.jmb.2003.08.057.
- Walsh, I., Fishman, D., Garcia-Gasulla, D., et al., 2021. DOME: recommendations for supervised machine learning validation in biology. Nature Methods 18, 1122–1127. doi:10.1038/s41592-021-01205-4.
- Yang, J., Mora, A., Liu, S., Wittmann, B. J., Anandkumar, A., Arnold, F. H., Yue, Y., 2024a. CARE: a benchmark suite for the classification and retrieval of enzymes. In: Advances in Neural Information Processing Systems 37 (NeurIPS 2024), Datasets and Benchmarks Track. arXiv:2406.15669.
- Yang, Y., Jerger, A., Feng, S., Wang, Z., Brasfield, C., Cheung, M. S., Zucker, J., Guan, Q., 2024b. Improved enzyme functional annotation prediction using contrastive learning with structural inference. Communications Biology 7, 1690. doi:10.1038/s42003-024-07359-z.
- Yang, Z., Su, B., Chen, J., Wen, J.-R., 2025. Interpretable enzyme function prediction via residue-level detection. In: ICLR 2025 Workshop on Learning Meaningful Representations of Life (LMRL). arXiv:2501.05644.
- Yu, T., Cui, H., Li, J. C., Luo, Y., Jiang, G., Zhao, H., 2023. Enzyme function prediction using contrastive learning. Science 379(6639), 1358–1363. doi:10.1126/science.adf2465.
- Zadrozny, B., Elkan, C., 2002. Transforming classifier scores into accurate multiclass probability estimates. In: Proceedings of the Eighth ACM SIGKDD International Conference on Knowledge Discovery and Data Mining, pp. 694–699. doi:10.1145/775047.775151.
- Zhang, M.-L., Zhou, Z.-H., 2014. A review on multi-label learning algorithms. IEEE Transactions on Knowledge and Data Engineering 26(8), 1819–1837. doi:10.1109/TKDE.2013.39.
- Zhou, N., et al., 2019. The CAFA challenge reports improved protein function prediction and new functional annotations for hundreds of genes through experimental screens. Genome Biology 20, 244. doi:10.1186/s13059-019-1835-8.

<!-- K:END -->

---

## L. Citation audit

### L1. Claim by claim

"Abstract" in the table means checked against the abstract, the publisher summary or a search excerpt; no full texts were read.

| # | Claim in draft | Citation | Supported? | Problem found | Correction made |
|---|---|---|---|---|---|
| 1 | EC numbers classify enzyme-catalysed reactions in a four-level hierarchy | McDonald2023 | Yes (review of EC nomenclature; EC-Bench's abstract states the same) | — | — |
| 2 | One protein can carry several EC numbers; unrelated proteins catalyse the same reaction | Omelchenko2010 | Yes for the second clause. The first clause is common knowledge; it is also stated in CLEAN's abstract ("promiscuous enzymes with two or more EC numbers"). | First clause not cited | None needed. Add Yu2023 there if a reviewer asks. |
| 3 | Transfer from homologues via BLAST or profiles (PRIAM) | Altschul1997, Camacho2009, ClaudelRenard2003 | Yes | The original draft said "the standard way … is still", an unsupported generalisation | Now "A long-standing way …" |
| 4 | Function less conserved between homologues than assumed | Rost2002 | Yes (title-level claim of the paper) | — | — |
| 5 | ~40% identity for >90% conservation of the first three digits, ~60% for all four | Tian2003 | Yes (abstract as quoted in a citing source) | Exact figures should be read in the paper | Listed in L3 |
| 6 | Across 37 families, misannotation rare in Swiss-Prot, frequent in automatic databases | Schnoes2009 | Yes | — | — |
| 7 | Screen of 122 sequences; ≥78% of EC 1.1.3.15 misannotated | Rembeza2021 | Yes | — | — |
| 8 | ECPred trains an independent classifier per EC number | Dalkiran2018 | Yes | — | — |
| 9 | DeepEC and ProteInfer use CNNs on sequence | Ryu2019, Sanderson2023 | Yes (well-known architecture; not re-read here) | — | Low risk; check in L3 if desired |
| 10 | DeepFRI combines language-model features with structure in a GCN | Gligorijevic2021 | Yes | The EC-number task was not confirmed from the abstract | Sentence claims only "function prediction" |
| 11 | DeepECtransformer: 464 genes, 3 confirmed | Kim2023 | Yes | The original draft said "adds transformer layers", which implied building on DeepEC | Now "uses transformer layers" |
| 12 | SOLVE: interpretable ensembles on short fragments can compete with deep models | Dhibar2025 | Yes ("surpassing most existing algorithms" per the publisher summary) | "Can compete" is deliberately weaker than the paper's own claim | — |
| 13 | CLEAN uses ESM-1b embeddings, a contrastive embedding and nearest cluster centres | Yu2023, Rives2021 | Yes (paper plus our code audit) | — | — |
| 14 | CLEAN reported more accurate and sensitive than BLASTp | Yu2023 | Yes (abstract) | Attributed as "reported" | — |
| 15 | CLEAN used as reference in independent benchmarks | Davoudi2026, Sartori2026 | Yes (EC-Bench predictions include CLEAN; EnzymARC tests CLEAN) | One source is a preprint | Stated as preprint in 2.1-P3 |
| 16 | Successors: contact maps, residue-level detection, hierarchy-aware losses, unseen ECs, hierarchical transformer | YangY2024, YangZ2025, Cao2026, Ma2026, Dumontet2026 | Yes | "Later methods add" implied they all extend CLEAN | Now "introduce" |
| 17 | Independent evaluations less favourable than original reports | (synthesis of rows 18–21) | Yes as a synthesis of the cited studies | "Independent tests are" was too categorical | Now "have been" |
| 18 | >450 *E. coli* unknowns; models mostly fail at novel predictions and make logic errors | deCrecyLagard2025 | Yes (abstract wording) | — | — |
| 19 | CLEAN and DIAMOND assign original EC to >90% of weakly perturbed decoys | Sartori2026 | Yes (abstract); **preprint** | — | Labelled "a 2026 preprint" |
| 20 | BLASTp marginally better overall; complementary; LLMs useful below 25% identity | Capela2025 | Yes (abstract) | — | — |
| 21 | Train-set-matched BLASTp; parity in distribution; PLMs well ahead for distant eukaryotes | Sathyamoorthy2026 | Yes (abstract); **preprint** | — | Labelled "Another 2026 preprint" |
| 22 | Outcome depends on identity profile and baseline configuration | synthesis of rows 20–21 | Reasonable inference, stated as "appears to" | — | — |
| 23 | We compare with BLASTp against CLEAN's own training sequences, by identity | our Methods | Yes (plan §0 C3) | — | — |
| 24 | CAFA: predict for unannotated proteins, scored on later annotations | Radivojac2013, Jiang2016, Zhou2019 | Yes | — | — |
| 25 | LAFA runs time-delayed evaluation continuously | Phan2026 | Yes | — | — |
| 26 | 37 GO metrics; some give erroneously high scores to false-positive sets | Plyusnin2019 | Yes (abstract wording) | — | Correction notice cited in bib note |
| 27 | Incomplete ground truth: a scored false positive may be an unannotated true function | Jiang2014 | Yes (premise of the paper) | — | — |
| 28 | Averages over proteins, over labels, or micro | Zhang2014 | Yes | — | — |
| 29 | Measures differ in invariance to confusion-matrix changes | Sokolova2009 | Yes | — | — |
| 30 | Two macro-F1 formulas can rank classifiers differently | Opitz2019 | Yes; arXiv note | Not peer reviewed | Kept because the finding is analytical and widely cited; labelled in the bib |
| 31 | EC studies commonly use support-weighted averaging; EC-Bench adopts it for imbalance | Davoudi2026 (+ our audit) | Yes (EC-Bench text; five code bases in our audit) | "Commonly" rests on our audit | Worded "commonly", not "always" |
| 32 | CARE: <30% and 30–50% identity splits; reaction→EC retrieval | YangJ2024 | Yes | — | — |
| 33 | EC-Bench: ten methods; exact prediction, completion, recommendation; uniform assumptions | Davoudi2026 | Yes | Conflicts with our abstract's "16 methods" | See 0.3, item 4 |
| 34 | ECPred temporal set: F1 0.14 at substrate level; such cases are very difficult | Dalkiran2018 | Yes | The authors' remark was about scores differing from published values | Paraphrase now matches |
| 35 | General guidance on validation and leakage | Walsh2021, Kapoor2023 | Yes | "All of this follows" implied a chronology | Now "in line with" |
| 36 | No study has asked whether EC averaging counts every false positive | — (absence claim) | Hedged "to our knowledge"; searches in M found none | Absence claims can be overturned | Re-run the M searches before submission |
| 37 | EC labels are imbalanced | Davoudi2026 | Yes | The original draft said "strongly" | Adjective removed |
| 38 | HIT-EC motivated by under-represented ECs and incomplete annotations | Dumontet2026 | Yes | — | — |
| 39 | EnzPlacer predicts ECs absent from training | Ma2026 | Yes (title) | — | — |
| 40 | A closed-set model cannot recall unseen labels | logical | Yes | The original draft said "share of proteins … sets an upper limit", imprecise for multi-EC proteins | Now "the test labels it has never seen cap its attainable recall" |
| 41 | CLEAN's rule = distance to mean embedding of training sequences | Yu2023 (+ code) | Yes | — | — |
| 42 | This is the prototype classifier; the embedding is trained with a supervised contrastive objective | Snell2017, Khosla2020 | Yes (CLEAN's losses are supervised contrastive: triplet and SupCon-Hard) | — | — |
| 43 | Prototype classifiers (incl. nearest class mean) counter bias towards frequent classes | Sharma2023 | Yes (they address the biased softmax in long-tailed recognition) | — | — |
| 44 | Centroid vs neighbour rules not compared in one frozen EC space | — (absence claim) | Hedged "to our knowledge" | EnzPlacer and FEDKEA use kNN transfer, but not a controlled comparison | Be ready to cite them if a reviewer asks (section O) |
| 45 | Modern neural networks often poorly calibrated | Guo2017 | Yes | — | — |
| 46 | Isotonic regression maps scores to calibrated probabilities | Zadrozny2002 | Yes | — | — |
| 47 | Selective classification guarantees kept-prediction risk with high probability | Geifman2017 | Yes | "Guaranteed" alone was too strong | "with high probability" added |
| 48 | CLEAN confidence = two-component GMM on embedding distances | Yu2023 | Yes (description of CLEAN seen in search results; also the released `gmm` file) | — | — |
| 49 | ECPICK uses evidential deep learning | Han2024 | Yes | — | — |
| 50 | CPEC limits the false discovery rate of EC annotations | DingK2024 | Yes (title plus our code notes) | — | — |
| 51 | Boger et al. applied conformal methods with user-chosen risk, including CLEAN embeddings | Boger2025 | Yes (abstract plus the repository's CLEAN experiments) | The original draft said "conformal risk control", which may be too specific | Now "conformal methods … with user-chosen risk guarantees" |
| 52 | These guarantees hold on average over test proteins | synthesis | Yes (standard property of conformal and risk-control guarantees) | — | — |
| 53 | Class-conditional conformal is hard with many classes and few examples | DingT2023 | Yes | — | — |
| 54 | Exact binomial limits determine the minimum calibration sample | Clopper1934 | Yes (the method; the application is ours) | "fix" was too strong | Now "determine" |
| 55 | Nobody has measured the distance of rare tiers from that minimum | — (absence claim) | Hedged "to our knowledge" | — | — |

**Result:** no citation was found that fails to support its sentence. Eleven wordings were tightened (rows 3, 11, 16, 17, 34, 35, 37, 40, 47, 51, 54). Three absence claims (rows 36, 44 and 55) are all hedged with "to our knowledge" and rest on the searches in section M.

### L2. Bibliographic fields to complete from the publisher pages

These fields were left out or shortened rather than guessed.

| Key | What to complete |
|---|---|
| Dalkiran2018, Capela2025, Radivojac2013, Jiang2016, Zhou2019, Walsh2021, Gligorijevic2021 | Full author lists ("and others" is used now) |
| Sartori2026 | Full author list. The first author is J. Sartori; co-authors seen include A. C. R. Guimarães and L. de A. Machado, but the order is unconfirmed. |
| Dumontet2026 | The author list is from the bioRxiv version (Dumontet, Han, Oh, Kang); confirm it on the *Nature Communications* page |
| Cao2026, Phan2026 | Volume (and issue) numbers for 2026 |
| Tian2003, Rost2002, Rembeza2021, Radivojac2013 | Issue numbers (optional in the Elsevier style) |
| Snell2017, Khosla2020, Sharma2023, DingT2023, YangJ2024 | Page numbers. Sources disagreed for Snell 2017 (4077–4087 vs 4080–4090), so pages were omitted. |
| YangZ2025 | This is an ICLR 2025 *workshop* paper (LMRL). Keep it as a code-audit reference; do not present it as a main-conference paper. |

### L3. Open these in full before submission

Your prompt asks for this for key papers; it was impossible here.

1. **Yu et al. 2023 (CLEAN):** the BLASTp configuration (database, E-value, top-hit rule), how the test sets were built (New-392, Price-149 provenance) and the exact metric wording.
2. **Davoudi et al. 2026 (EC-Bench):** metric definitions, the list of methods and predictors, the temporal split.
3. **Capela et al. 2025:** split design, identity bins, how BLASTp was run.
4. **Sathyamoorthy & Puri 2026:** what "train-set-matched" means, taxonomic composition, EC levels.
5. **Plyusnin et al. 2019** and its correction: the exact false-positive-set construction.
6. **Tian & Skolnick 2003:** the 40%/60% figures and conditions.
7. **Boger et al. 2025:** which CLEAN quantities are used and the form of the guarantee.
8. **Ding K. et al. 2024 (CPEC):** the guarantee (marginal FDR) and the base model.
9. **Sartori et al. 2026 (EnzymARC):** decoy construction and CLEAN settings.

---

## M. Missing-literature audit: second search, 2025–2026 overlap check

### M1. What was searched

The searches were done on 2026-09-29 and in the earlier sessions recorded in `PUBLICATION_PLAN.md` §3 and §14.

- **Exact terms and synonyms:**
  - "enzyme commission number prediction"
  - "EC number prediction benchmark"
  - "weighted F1 enzyme"
  - "support-weighted"
  - "macro F1"
  - "multi-label evaluation metrics"
  - "gene ontology evaluation metrics"
  - "incomplete annotation evaluation"
- **Methods:**
  - "contrastive learning enzyme function"
  - "protein language model EC prediction BLAST"
  - "prototype classifier long-tailed"
  - "conformal prediction enzyme"
  - "class-conditional conformal"
  - "selective classification"
  - "calibration neural networks"
  - "Clopper-Pearson"
- **Named works:** CLEAN, CLEAN-Contact, EnzymeHunter, ProtDETR, EnzPlacer, HIT-EC, ECPICK, SOLVE, DeepEC, DeepECtransformer, ProteInfer, DeepFRI, ECPred, PRIAM, EC-Bench, CARE, CAFA1–3, LAFA, CPEC, EnzymARC, "The Accuracy Paradox", "Rethinking Benchmarks and Models for Enzyme Specificity Prediction".
- **Backward chaining:** the methods CLEAN compares against (DeepEC, ProteInfer, ECPred, BLASTp) and the evaluation literature behind CAFA (Clark & Radivojac 2013, Jiang 2014).
- **Forward chaining** (work that uses or evaluates CLEAN): EC-Bench, EnzymARC, Boger 2025, CLEAN-Contact, EnzymeHunter, EnzPlacer, HIT-EC, Capela 2025.
- **Regional:** "EC number prediction Nepal"; "enzyme function prediction India".

**Limits.**
- No Google Scholar or Scopus access, so there was no systematic citation count.
- Forward chaining relied on search engines, not "cited by" lists.
- Theses were not searched.

### M2. Recent work that overlaps (2025–2026), flagged prominently

| Work | Status | Overlap with | Verdict |
|---|---|---|---|
| **Capela et al. 2025**, *BMC Bioinformatics* | Peer reviewed | C3b (BLASTp vs PLMs) | **Reduces novelty.** Same direction as our result, so C3b becomes Type 3. Cited. |
| **Sathyamoorthy & Puri 2026**, bioRxiv | Preprint | C3b | **Partly contradicts** our tie below 30% identity. Cited; discuss possible reasons (trained classifiers vs frozen CLEAN; eukaryotic distant set; baseline configuration). |
| **Sartori et al. 2026 (EnzymARC)**, bioRxiv | Preprint | Open-set false positives | **Strengthens the rationale** for counting all false positives. Cited. |
| **Phan et al. 2026 (LAFA)**, *Bioinform Adv* | Peer reviewed | C2 (temporal evaluation) | **Partial overlap.** GO-centred continuous evaluation. Cited as precedent; our novelty is EC-specific decomposition and vintage. |
| **Davoudi et al. 2026 (EC-Bench)**, *Bioinform Adv* | Peer reviewed | C1c, C2 | **Partial overlap.** It is the benchmark we re-score; it motivates weighted F1. Cited. |
| **Dumontet et al. 2026 (HIT-EC)**, *Nat Commun* | Peer reviewed | Rare ECs, "trustworthy" | **Strengthens the rationale** (rare ECs are a recognised problem). Cited. |
| **Ma et al. 2026 (EnzPlacer)**, *Bioinformatics* | Peer reviewed | Unseen ECs (C2) | **Partial overlap.** It solves unseen ECs; we measure how many a released model cannot reach. Cited. |
| **Cao et al. 2026 (EnzymeHunter)**, *Patterns* | Peer reviewed | C1b | **Strengthens C1b.** Its public code uses the same weighted evaluation. Cited. |
| *The Accuracy Paradox* (Ahmad, arXiv:2609.07897, Sept 2026) | Preprint, single author | C1 | **No overlap in mechanism.** It is about default 0.5 thresholds for top-level EC classes predicted from molecular descriptors. Not cited; mention only if a reviewer raises it. |
| *Rethinking Benchmarks and Models for Enzyme Specificity Prediction* (Mahood et al., arXiv:2607.05084) | Preprint | Homology baselines | **Tangential.** It covers enzyme–substrate specificity, not EC numbers; most models do not beat BLAST baselines. Optional support for C3b; not cited. |
| PAD dataset (bioRxiv 10.64898/2025.12.17.694800) | Preprint | Benchmarks | Resource only; no overlap. Not cited. |
| EC-Reason-Bench (arXiv:2607.26397) | Preprint | LLM enzyme classification | Different question. Not cited. |

**None of these eliminates the C1 gap.** Two of them (Capela, Sathyamoorthy) change how C3b must be framed. LAFA, EC-Bench and EnzPlacer mean C2 must be framed as a partial extension.

### M3. Papers a reviewer might expect but not cited in Section 2

Some are for other sections; some are optional.

| Work | Why a reviewer might expect it | Recommendation |
|---|---|---|
| UniProt Consortium, *NAR* 53:D609 (2025), doi:10.1093/nar/gkae1010 (V2) | Data source (Swiss-Prot) | **Cite in Methods** |
| Price et al., *Nature* 557:503 (2018), doi:10.1038/s41586-018-0124-0 (V2) | Likely origin of Price-149 | Cite in Methods **after confirming the provenance in CLEAN's paper** |
| scikit-learn (Pedregosa et al. 2011) | The weighted average is scikit-learn's implementation | Cite in Methods (not verified this session) |
| ENZYMAP, *PLoS ONE* 2014, doi:10.1371/journal.pone.0089162 | EC reclassification over time (our transfer mapping) | Methods; verify first |
| FEDKEA (*Sci Adv* 2026, doi:10.1126/sciadv.aee4389) and EnzPlacer (kNN transfer) | Neighbour transfer in embedding space (C3a) | Optional in 2.3 if a reviewer questions the absence claim |
| TopEC (*Nat Commun* 2025), KAN-based EC prediction (npj AI 2025), EnzymeNet (2023), ECRECer (2023), GloEC (2024), MAPred (2025) | Other EC predictors | Optional; they add nothing to the argument |
| Conformal foundations (Vovk et al. 2005; Angelopoulos & Bates), Learn-then-Test | Background for conformal statements | Methods, if conformal tools are used there |
| ESM-2 (Lin et al. 2023) | Mentioned via Capela / Sathyamoorthy | Optional |
| Clark & Radivojac 2013 (*Bioinformatics* 29:i53, doi:10.1093/bioinformatics/btt228, V2) | Information-theoretic CAFA metrics | Optional; dropped to keep 2.2 short |

### M4. Redundancy check

- **CAFA trio** (Radivojac2013, Jiang2016, Zhou2019): Radivojac alone carries the point. Cut Jiang2016 and Zhou2019 first if space is tight.
- **Walsh2021 and Kapoor2023** overlap (validation and leakage guidance). Keep Walsh (biology-specific) if one must go.
- **Five learned-predictor examples** (ECPred, DeepEC, ProteInfer, DeepFRI, SOLVE): ECPred stays (it is also the temporal precedent). DeepFRI and SOLVE are the first cuts; SOLVE's value is that it is a chemistry-journal reference.
- **Altschul1997 and Camacho2009**: different roles (algorithm vs the software version used). Keep both, or move Camacho to Methods.
- **Capela2025 and Sathyamoorthy2026**: *not* redundant. They disagree, and the disagreement is the point.
- **Han2024 and Dumontet2026**: same group but different roles (uncertainty vs rare-EC motivation). Keep both.

### M5. Nepal / South Asia

- No literature from Nepal on EC-number prediction or its evaluation was found.
- The topic has no local material, environmental or regional dimension, so the absence is not a gap to claim.
- From South Asia, SOLVE (IACS Kolkata, India; *Chemical Science* 2025) is a relevant primary study and is cited on its merit.

---

## N. Must-cite papers

**Essential, the 15 core citations for the central argument:**
1. Yu et al. 2023 (CLEAN): the model under study.
2. Davoudi et al. 2026 (EC-Bench): the benchmark re-scored; uses support-weighted F1.
3. Plyusnin et al. 2019: closest precedent for C1.
4. Sokolova & Lapalme 2009: the invariance framework for the proof.
5. Jiang et al. 2014: incomplete ground truth (the main caveat for C1).
6. Radivojac et al. 2013 (CAFA): prospective evaluation standard.
7. Dalkiran et al. 2018 (ECPred): the early EC temporal hold-out.
8. Yang J. et al. 2024 (CARE): identity-controlled EC benchmark.
9. Tian & Skolnick 2003: identity thresholds for EC transfer.
10. Capela et al. 2025: BLASTp vs PLMs (same direction as C3b).
11. de Crécy-Lagard et al. 2025: independent failure on novel functions.
12. Ma et al. 2026 (EnzPlacer): unseen EC numbers.
13. Ding K. et al. 2024 (CPEC): EC-specific conformal guarantees.
14. Boger et al. 2025: guarantees on CLEAN embeddings.
15. Ding T. et al. 2023: class-conditional guarantees need per-class data.

**Strongly useful:**
- Schnoes2009, Rembeza2021, Omelchenko2010, Rost2002
- Camacho2009, Altschul1997, Rives2021
- Ryu2019 (DeepEC, also an EC-Bench method), Kim2023, YangY2024, Cao2026, YangZ2025 (the last three are needed for the code audit), Dumontet2026
- Sartori2026, Sathyamoorthy2026 (preprints, but the contradiction must be discussed)
- Phan2026, Opitz2019, Zhang2014
- Guo2017, Geifman2017, Zadrozny2002, Clopper1934
- Walsh2021, Snell2017, Sharma2023, Han2024

**Optional:** McDonald2023, ClaudelRenard2003, Gligorijevic2021, Dhibar2025, Sanderson2023, Jiang2016, Zhou2019, Khosla2020, Kapoor2023.

---

## O. Reviewer perspective

### 1. What would a reviewer question?

- **"Isn't the metric issue trivial or known?"** Answer: the zero-weight property is elementary and documented for scikit-learn. What was not known is its size for EC benchmarks (67–96% of false positives uncounted), its presence in five public code bases, and its effect on EC-Bench rankings. Plyusnin et al. (2019) show the same *class* of problem in GO. Cite them rather than wait to be told.
- **"Are the uncounted predictions really errors?"** Incomplete labels (Jiang 2014) mean some are not. Use the E1 label-drift analysis, report how many of those predictions later appeared in Swiss-Prot, and phrase precision claims relative to the reference annotation.
- **"Why does your BLASTp beat CLEAN when CLEAN's paper reports the reverse?"** Read CLEAN's BLASTp configuration (L3.1) and state the difference. Our result agrees with Capela et al. (2025).
- **"Only CLEAN?"** Say so in the title/abstract scope and in the Limitations; the evaluation card is model-agnostic.
- **"Is +0.043 union-macro F1 meaningful?"** It is modest and pre-registered, and it is smaller than the BLASTp gap. Say so; the plan already does.
- **"What do the certification-budget numbers assume?"** Exchangeability between calibration and future data, a fixed target risk, and a stable curation rate. They are lower bounds, not predictions.

### 2. Where is the review weak?

- It relies on abstracts, not full texts (L3).
- Two preprints carry specific points: EnzymARC and the PLM-vs-BLAST contradiction.
- The chemistry framing is thin for a *Computational Biology and Chemistry* audience. One sentence on why EC errors matter chemically (wrong reactions in metabolic reconstructions or biocatalyst searches) would help, but it needs a citation; Omelchenko 2010 covers pathway reconstruction.
- Hierarchical (partial-credit) EC evaluation is mentioned only through EC-Bench's tasks.

### 3. Which novelty claims could be challenged?

- C3b (Capela 2025).
- C1a ("elementary").
- C2 (temporal hold-outs exist: CAFA, ECPred, LAFA, New-392).
- C4a (conformal guarantees for CLEAN already exist: Boger 2025).

Section G gives defensible wordings for each.

### 4. What literature might still be missing?

- Methods references (UniProt, scikit-learn, Price 2018, ENZYMAP).
- Neighbour-transfer tools (FEDKEA, EnzPlacer's kNN) for the C3a absence claim.
- Conformal foundations, if conformal tools appear in Methods.
- Any late-2026 successor that reports temporal or calibrated EC results. Re-run M1 before submission.

### 5. What would strengthen the paper?

1. The label-drift check on uncounted false positives (item 1 above).
2. Report BLASTp settings sensitivity (E-value threshold, top-k) and optionally DIAMOND, which EnzymARC and EC-Bench use.
3. If a successor's released weights run on CPU, apply the evaluation card to it (for example, EnzymeHunter or CLEAN-Contact). Even a partial check answers "only CLEAN?".
4. Add a Results table that summarises, for each audited benchmark or code base, the averaging used and whether out-of-label predictions are counted (from the plan §3/§14 audit).

### 6. What statements should be softened?

- Say "is not counted by" rather than "hides" where precision matters. The abstract's "hides 67–96% of CLEAN's false positives" is defensible but vivid; keep it if space is short.
- Say "five public code bases" rather than "field-wide".
- Say "16 predictors" rather than "16 methods" (0.3, item 4).
- Any "first" claim should become "to our knowledge", with the specific scope.

### 7. What would make the Introduction and review more convincing?

- Open the Introduction with the practical stake: EC annotations drive metabolic reconstruction and enzyme discovery, and wrong EC numbers propagate (Schnoes 2009).
- End it with the three open questions from 2.5 mapped one-to-one to C1–C4 and to the Results subsections.
- Keep the Introduction short (about 500–600 words) so that it does not repeat Section 2.

---

## P. How to put the section into the manuscript

These edits are to `paper/manuscript.tex`; I did not make them, since you are writing the other sections.

1. After the Introduction, add `\input{literature_review}`.
2. Replace `\bibliography{cas-refs}` with `\bibliography{references}`.
3. **Recommended:** change `\usepackage[authoryear,longnamesfirst]{natbib}` to `\usepackage[authoryear]{natbib}`. With `longnamesfirst`, the first citation of every work prints its full author list in the text; Rives et al. would print 11 names.
4. Compile on your PC with the Elsevier cas-sc bundle (`cas-sc.cls`, `cas-model2-names.bst`). There is no TeX installation here, so the section has been checked by script but not compiled.
5. The two "Yang et al., 2024" entries will print as 2024a (CARE) and 2024b (CLEAN-Contact); that is the style's normal behaviour.
