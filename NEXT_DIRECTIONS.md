# Next directions: novel follow-up ideas and work plan

**Written:** 2026-09-29.

**Scope:** follow-up research after the CLEAN evaluation paper for *Computational Biology and Chemistry*. That paper is described in `PUBLICATION_PLAN.md`; its contributions C1–C4 are referred to below.

**Method:**
- About 30 web searches, on 2026-09-29 and in earlier sessions.
- Feasibility checks against the data in this repository.
- One external repository inspected: EnzPlacer at commit 97a8f0a.

**Status:** plan only. No code has been written, no test-set result was computed for any new idea, and nothing is committed.

**Novelty types** follow your literature prompt:
- **Type 1:** completely unexplored.
- **Type 2:** incremental but meaningful.
- **Type 3:** contextual or application novelty.

No idea here is Type 1. I say so rather than oversell.

---

## 0. Bottom line

| Rank | Idea | Type | Why it is strong | Effort | Where it goes |
|---|---|---|---|---|---|
| **1** | **Label provenance.** Do Swiss-Prot benchmarks favour homology search because many test labels were themselves assigned "by similarity"? | 2 | See below the table. | ~1 week for the core check; ~4 weeks for the full study | A small pre-registered check (E7-lite) **before** the CBC submission, then paper 2 |
| **2** | **Benchmark hygiene.** Training-vintage overlap and statistical power across the EC-prediction literature. | 2 | Already motivated by a concrete finding (section 3.1). A minimum-detectable-effect table costs one day and strengthens the current paper. | 1 day for power; 6–8 weeks for the audit | Power table in the current paper; the audit joins idea 1 in paper 2 |
| **3** | **Certified open-set annotation.** Conformal detection of enzymes with new EC numbers, a certifiable-depth atlas and prospective validation. | 2 | Extends C4 into a positive method with formal guarantees, and adds a truly prospective test. | 6–10 weeks, plus waiting for new Swiss-Prot releases | Paper 3 |

**Why idea 1 is strong:**
- It tests the most contested claim in our paper (BLASTp beats CLEAN) at its root.
- Any outcome is publishable.
- It needs one download plus predictions already in the repo.
- It generalises to EC-Bench's 16 predictors.

**Recommendation:**
1. Do the E7-lite check (section 2.5) before submitting the CBC paper. Its result could change how C3 is worded.
2. Then build paper 2 from ideas 1 and 2 together: *label provenance, training-vintage overlap and statistical power in enzyme-function benchmarks*.
3. Start idea 3 in parallel, because its prospective test set needs months to accumulate.

---

## 1. How the ideas were chosen

**Excluded, because they are done in our paper, already published, or need a GPU:**
- Chemistry-weighted error metrics: EC-Bench already scores partial and additional ECs.
- Hierarchy-aware or rare-class training losses: GloEC, EnzPlacer, EnzymeHunter, HIT-EC.
- Structure fusion: CLEAN-Contact, EnzymeHunter.
- Anything that needs retraining CLEAN (no GPU confirmed).
- Rarity-stratified calibration and certification budgets: our E4 and E6.

**Kept if they met four tests:**
1. They answer a question our results leave open.
2. They can run on CPU with the data in this repo or public downloads.
3. They respect the data protocol (New-392 for development, CAL for tuning, TEST never used for choices, new dated pre-registrations).
4. Their outcome is informative in either direction.

**Checked and rejected** (details in section 6):
- EnzPlacer's "unseen" test set as an independent open-set test. It is mostly inside CLEAN's training data.
- Open-world precision correction. It is already studied for Gene Ontology.
- Taxonomy-only stratification. The EC-Bench test set is 93% eukaryotic.
- Sparse-autoencoder interpretability. Its claims are hard to validate.
- Curation-allocation optimisation. It is speculative.

---

## 2. Idea 1 (flagship): label provenance

### 2.1 The question

Swiss-Prot records *how* each annotation was made, using Evidence and Conclusion Ontology (ECO) codes:
- **ECO:0000269** means experimental evidence used in a manual assertion.
- **ECO:0000250** means sequence-similarity evidence used in a manual assertion: the annotation was propagated from a related, experimentally characterised protein.

If a test protein's EC number was assigned by similarity to a protein that is in the training set, then a BLASTp top hit essentially repeats the curator's own inference.

**Question:** is BLASTp's advantage over learned predictors concentrated in similarity-labelled proteins, and does it survive on experimentally labelled ones?

### 2.2 Why it matters now

**Our own results make the question sharp.** On the temporal TEST set, BLASTp beats CLEAN: seen-EC micro-F1 0.735 vs 0.641, difference +0.093 [0.070, 0.116]. On New-392 BLASTp is also ahead (+0.086). On Price-149, whose labels come from experimental data (to confirm in the CLEAN and ProteInfer papers), the direction flips: −0.051 [−0.116, +0.015], not significant.

**It speaks to the contradiction in the literature** (literature package, 0.3):
- CLEAN's paper reports CLEAN ahead of BLASTp.
- Capela et al. (2025) find BLASTp marginally ahead.
- Sathyamoorthy & Puri (2026, preprint) find protein-language-model classifiers far ahead for distant enzymes.

Differences in how the test labels were produced could explain part of this.

**The field knows the problem but has not measured it:**
- CAFA builds its benchmarks only from GO annotations with experimental or curator-statement evidence codes, excluding computationally inferred ones. This avoids scoring methods against labels that computational inference produced.
- For EC benchmarks, a recent paper remarks in passing that "almost all existing database annotations were assigned via homology". This was seen in search excerpts; confirm the source (EC-Bench or arXiv:2606.12209) in the full text.
- EC-Bench responds with a 30%-identity training regime, not with evidence stratification.
- **No study was found that stratifies EC-benchmark performance by label evidence.** The searches covered evidence-code-stratified EC evaluation, annotation circularity, homology-inflated benchmarks, experimental-only EC benchmarks, and Swiss-Prot evidence proportions.

**Novelty: Type 2.** The concept is known from Gene Ontology and CAFA practice. The quantification for EC benchmarks, the mechanism test (next item) and the reconciliation of conflicting BLAST-vs-PLM results are new, to our knowledge.

**A direct mechanism test is available.** A UniProt ECO:0000250 evidence usually names its source protein (e.g. `ECO:0000250|UniProtKB:P12345`). Our BLAST results already store each query's top-hit accession (`data/cleanx/identity/*.csv`, column `top_hit`) and the 10 best hits (`identity/blast/*.tsv.gz`). So we can measure how often BLASTp's top hit is literally the protein the curator copied from. I found no report of this measurement.

### 2.3 Feasibility, as checked on 2026-09-29

**Model outputs already exist:**
- CLEAN (all inference rules) and BLASTp top-hit predictions per protein, for TEST, New-392 and Price-149: `results/cleanx/E2_split100/E2_per_query.csv`.
- For CAL: CLEAN in `E4_features_cal_centroid.csv`, BLASTp in `identity/sp_temporal_cal.csv`.
- The EC-Bench test set: 468 UniProt accessions with the predictions of 16 predictors (`data/cleanx/downloads/ec-bench/results_test_*.csv`).

**Sample sizes:**
- TEST has 1,260 enzymes, 876 of them with seen ECs.
- A proxy for evidence mix: 79% of TEST (992/1,260) and 80% of CAL (671/839) are "Evidence at protein level".
- That is *protein-existence* evidence, not EC evidence, so the real EC-evidence split is unknown until downloaded. If most labels are experimental, the similarity stratum may be small (decision gate G7, section 2.5).

**What is missing:** the EC-level evidence codes and sources.
- They come from the UniProt REST API as JSON: EC numbers in the protein description and in catalytic-activity comments, each with its evidences.
- The API is blocked from this analysis machine (HTTP 403), so the download runs on your PC.
- About 3,000 accessions (New-392 + CAL + TEST + EC-Bench test) take a few minutes.
- Record the UniProt release and a SHA-256 digest in the manifest, as the existing downloads do.

**Compute:** seconds. No CLEAN weights are needed for the analysis, since it uses stored predictions.

### 2.4 Design (draft pre-registration, E7)

Commit this in the script docstring before any TEST number exists, as was done for E5 and E6.

**Label evidence per EC:**
- Collect the ECO codes attached to that EC in the protein name and in catalytic-activity reactions with that EC.
- **EXP** if any code is ECO:0000269. Otherwise **SIM** if any code is ECO:0000250. Otherwise **OTHER** (e.g. ECO:0000305 curator inference, rule-based codes, no evidence).

**Label evidence per protein:**
- **EXP** if all its true ECs are EXP.
- **SIM** if at least one is SIM.
- **OTHER** otherwise.

The pilot may refine these mapping rules. The endpoints below may not change.

**Population:** SP-temporal TEST proteins whose true ECs are all seen in split100 after transfer mapping (n = 876). BLASTp top-hit transfer is closed-set too, so unseen-EC proteins are uninformative for this contrast.

**Estimands:** Δ(class) = micro-F1(BLASTp top hit) − micro-F1(CLEAN centroid, max-separation), EC level 4, within an evidence class.

**Primary endpoints** (Holm-adjusted across the two, α = 0.05):
- **P7a:** Δ(EXP). Does BLASTp still beat CLEAN on experimentally supported labels?
- **P7b:** Δ(SIM) − Δ(EXP). Does provenance modulate the gap?
- **Method:** stratified paired bootstrap over proteins, resampled within class (B = 2,000, seed 20260927), 95% percentile CI.

**Key secondary endpoints:**
1. The same contrasts after identity standardisation. Δ is computed within BLAST identity bins (<30, 30–50, 50–70, 70–90, ≥90%) and re-weighted to the SIM class's bin distribution. SIM labels should sit at higher identity, so this separates "provenance" from "identity".
2. EC-Bench test (100% regime, EC4): Δ(SIM) − Δ(EXP) for blastp vs enzbert-learnt and for blastp vs clean-supcon, Holm-adjusted over these two. EC-Bench's labels are from Swiss-Prot 2023_01; ECs no longer present in the current entry form a separate "label changed" class and are excluded.
3. Per-protein F1 versions of P7a and P7b.

**Mechanism (descriptive, with exact CIs):** among SIM ECs with a cited source accession, report:
- the share where BLASTp's top hit is the source;
- the share where the source is among the 10 kept hits;
- the share where the source is in split100;
- BLASTp's micro-F1 when the top hit is the source vs when it is not.

**Sensitivity analyses:**
- EC-cluster bootstrap.
- At most 5 proteins per EC.
- Excluding OTHER.
- Price-149 as a fully experimental set. This is **post hoc** under the data protocol, since it was designed after 2026-09-27.

**Decision gate G7:** if n_EXP < 100 or n_SIM < 50 among seen-EC TEST proteins, the affected endpoint is reported descriptively (estimate and CI, no test). The EC-Bench and pooled analyses then carry the argument.

**How each outcome would be read** (fixed in advance):

| P7a Δ(EXP) CI | P7b DiD CI | Claim |
|---|---|---|
| > 0 | includes 0 | "BLASTp's advantage holds on experimentally supported labels; not a provenance artefact". This strengthens C3 as written. |
| > 0 | > 0 | "BLASTp's advantage is real but inflated where labels were transferred by similarity". Report both numbers in C3. |
| includes 0 or < 0 | > 0 | "BLASTp's advantage on Swiss-Prot benchmarks is largely a label-provenance effect". **C3 and the abstract must be reworded before submission.** |

**Pilot (allowed by the protocol):** New-392 (development set) and CAL (the tuning set), to check parsing, class sizes and the mapping rules.

### 2.5 E7-lite: the part to do before the CBC submission

E7-lite is P7a and P7b plus the mechanism counts on TEST, and nothing else. It answers whether C3 is safe as worded.

| Step | Where | Time |
|---|---|---|
| 1. Write `cx_12_provenance.py`. It downloads UniProt JSON for listed accessions in batches, parses EC-level evidence and sources, and writes `data/cleanx/provenance/<set>.csv` plus a manifest (release, date, SHA-256). It needs a `--smoke` mode (5 accessions) and an `--offline` mode that parses cached JSON. | here (code and tests) | 1 day |
| 2. Write `cx_13_provenance_eval.py`, with the pre-registration in its docstring. It builds the stratified bootstrap, the identity standardisation and the mechanism counts, and writes `results/cleanx/E7/` (CSV plus `E7_summary.md`). | here | 1 day |
| 3. Add tests. The unit test is a known-answer difference-in-differences on a hand-built table. The integration test uses synthetic UniProt JSON fixtures in a temp folder (never writing to `data/esm_data` or the real `results/`). Add entries to `cx_smoke.bat` and `cx_full.bat`, and keep `cx_00_check.py` at 9/9. | here | ½ day |
| 4. Commit the frozen code, which registers E7 (needs your "commit and push"). | — | — |
| 5. Run `.\cx.bat cx_12_provenance.py`, then the pilot, then the freeze, then `.\cx.bat cx_13_provenance_eval.py`. | your PC | 30 min |
| 6. Audit the results. If needed, reword C3 and the abstract, and add one paragraph plus a supplementary table. | here | ½–1 day |

### 2.6 The full study (paper 2, section 3.3)

The full study adds three things:
1. EC-Bench across all 16 predictors, which generalises beyond CLEAN.
2. The 30% training regime of EC-Bench.
3. A time analysis: do annotations created by similarity become experimentally supported later?

The last one uses the `date_created` field and the evidence dates in UniProt, if available. Check before planning on it.

### 2.7 Risks

- **Small SIM class.** Mitigated by gate G7 and the EC-Bench test set.
- **Confounding by identity, organism or family.** Mitigated by identity standardisation, the EC-cluster bootstrap and the per-EC cap. Organism lineage from the same download gives a taxonomy secondary.
- **Evidence codes change over time.** We use the codes in the current release, matching the current labels used for TEST, and state this.
- **ECO codes are not recorded consistently for EC numbers.** Parse both the protein-name and the catalytic-activity evidences, and report the share with no evidence at all.

---

## 3. Idea 2: benchmark hygiene (training vintage and statistical power)

### 3.1 Motivating findings, checked on 2026-09-29

- **EnzPlacer's "unseen" test set is not unseen for the released CLEAN.** Its `test_unseen_experimental.csv` holds 9,901 proteins and 3,682 ECs (repository commit 97a8f0a). **8,348 of those proteins (84%) are in CLEAN's split100 training set**, and 3,363 of the ECs (91%) are in CLEAN's vocabulary. The set is "unseen" only relative to EnzPlacer's own seen/unseen split. This is not a flaw in EnzPlacer, whose split is internally consistent. But a comparison that used *released* CLEAN weights on this set would be contaminated. Whether any published comparison did so has to be checked in the papers; I have not checked it and make no claim.
- **Practice varies between papers.** Both points below come from search excerpts; verify them in the full texts.
  - A recent benchmark paper notes that DeepECtransformer, which it could not retrain, may have had data leakage. Confirm whether this is EC-Bench.
  - SOLVE removed all New-392 sequences from its training set (UniProt up to 2023).
  - Other successors that report New-392 and Price-149 results, such as ProtDETR, still need their training vintage checked.
- **Standard benchmarks are small.** For the BLASTp–CLEAN pair, our paired bootstrap CI on Price-149 (n = 149) has a half-width of about 0.066. That implies a minimum detectable micro-F1 difference of about **0.09** at 80% power and two-sided α = 0.05, using MDE ≈ 2.8 × SE. On SP-temporal TEST (n = 1,260) the same pair gives about **0.03**. Many published "improvements" on Price-149 and New-392 are a few points; how many are detectable is an open, checkable question.
- **Precedent:** power analysis is established in NLP (Card et al., EMNLP 2020, "With Little Power Comes Great Responsibility"). No EC-specific power analysis was found.

### 3.2 Quick win for the current paper (1 day)

Add a **minimum-detectable-effect (MDE) field** to the EC evaluation card: an MDE row for New-392, Price-149, EC-Bench test and SP-temporal TEST.

- Compute it from the paired per-protein outcomes that already exist (`E2_per_query.csv`, EC-Bench predictions) with the existing bootstrap.
- This is descriptive, not a new hypothesis test, so it needs no pre-registration. It still goes into the Methods as a dated addition.

### 3.3 Full audit (paper 2, joined with idea 1)

1. **Vintage-overlap matrix.** Rows are models with public training data: CLEAN splits, CLEAN-Contact and ProtDETR (CLEAN splits), EnzPlacer (`train_seen.csv`), EnzymeHunter, GloEC, SOLVE, and others where released. Columns are benchmarks: New-392, Price-149, EC-Bench test, CARE splits, SP-temporal TEST, EnzPlacer unseen. Each cell reports exact accession/sequence overlap and ≥90%/≥50% identity (BLAST on your PC). Evidence is recorded as file@commit, as in C1b.
2. **Claims table.** For each paper that reports an improvement over CLEAN, record:
   - the dataset, metric and claimed difference;
   - whether the metric is blind to out-of-label predictions (C1);
   - whether the difference exceeds the benchmark's MDE;
   - whether the training data overlap the benchmark;
   - whether baselines were retrained or released weights were used.

   This needs the papers read in full, which is your part, since full texts are blocked here.
3. **Deliverable:** paper 2, *"Label provenance, training vintage and statistical power: how much of the progress in enzyme-function prediction is measurable?"*

**Novelty: Type 2.** Leakage and power are known issues in general machine learning (Kapoor & Narayanan 2023; Card et al. 2020). A systematic, file-level audit for EC prediction was not found.

### 3.4 Risks

- **Manual effort.** Mitigated by limiting the audit to the 10–12 most-cited successors.
- **The audit may find little overlap.** It is still publishable as "the field is clean on these axes", with the provenance result as the lead.
- **Tone.** Report what the files and papers state (file@commit), never intent, and send each group its row before submission.

---

## 4. Idea 3: certified open-set annotation

### 4.1 The question

Can we flag, with a controlled false-discovery rate, proteins whose EC number is **not among those a model knows**? About 30% of new enzymes in our temporal set carry such an EC number. Can we then say how deep an annotation can be certified?

### 4.2 Literature status

**General methods exist:**
- Conformal outlier p-values with Benjamini–Hochberg FDR control (Bates et al., *Ann. Stat.* 51:149, 2023).
- Open-set conformal prediction with unseen labels and imbalance (Xie, Zhou, Liang, Favaro & Sesia, arXiv:2510.13037, Oct 2025).
- Hierarchical conformal classification (arXiv:2508.13288; Mortier et al., arXiv:2501.19038).

**For EC prediction:**
- CPEC falls back from EC4 to EC3.
- Boger et al. (2025) use a hierarchical loss.
- EnzPlacer places unseen EC4 proteins at higher levels. Our earlier reading of its paper noted that it has no calibrated open-set abstention; re-verify that.

**Not found:** an application of conformal novelty detection to EC prediction with temporal ground truth. **Type 2.**

### 4.3 Test-data caveat

This decides the design.

**The TEST set is compromised for this idea.** Our E4 already reported open-set abstention and AUROC on TEST descriptively. A new method designed now and evaluated on TEST would be partly informed by those numbers. Report TEST results with that disclosure.

**External sets:**
- **EnzPlacer's unseen set cannot serve** (section 3.1).
- **A prospective set is the clean test.** It consists of Swiss-Prot entries created after 2026-09-02 (release 2026_03), built with the `cx_02` rules. Release 2026_03 added only **247 sequences** since 2026_02, per the Swiss-Prot release statistics. At that pace, about 150–250 new enzymes accumulate in 6–9 months.
- **Pre-register a fixed evaluation date.** For example, the first release on or after 2027-03-31, or n ≥ 150 enzymes, whichever comes later. Freeze all thresholds now, with the commit hash in the pre-registration.

### 4.4 Plan

1. **Novelty detection.** Build conformal p-values from CAL seen-EC enzymes (inliers). Pre-fix the score: CLEAN's distance to the nearest centroid, with BLAST's maximum identity as a secondary. Apply Benjamini–Hochberg at q ∈ {0.1, 0.2}.
   - **Endpoints:** realised FDR and power, with exact CIs.
   - **Evaluation sets:** TEST (disclosed), the 2,000 non-enzymes, and the prospective set.
   - **Endpoint rule:** these are reliability endpoints, as in E4 and E6. Your rule restricting endpoints to micro, union-macro and per-protein F1 governs classification endpoints. Tell me if you want it applied here too.
2. **Open-set prediction sets** following Xie et al. (2025), after you read the paper. The algorithm details are in the full text, which is not reachable here.
3. **Certifiable-depth atlas.** Extend E6's budget to EC3, EC2 and EC1 predictions, and map the deepest level that can be certified today in each part of EC space, with projected years for deeper levels.
4. **Prospective run** at the pre-registered date. It needs ESM-1b embeddings of about 200 new sequences on your PC (minutes).

### 4.5 Risks

- **CLEAN's distances may separate unseen-EC proteins poorly.** That would be a valid negative result, with power reported.
- **Family concentration.** The fungal LPMO family (1.14.99.56) dominates the unseen tier, so use a per-EC cap and the EC-cluster bootstrap.
- **The prospective sample will be small.** Pre-register the minimum n and report exact CIs.

---

## 5. Quick wins before the CBC submission

1. **E7-lite (section 2.5).** This is the most important; it could change C3's wording.
2. **MDE row in the evaluation card (section 3.2).**
3. **Open-world check of the "uncounted false positives".** Use `E1_label_drift.csv`, which compares 2022 benchmark ECs with current complete ECs: how many of CLEAN's New-392 false positives are now annotated ECs? New-392 is the development set, so this is descriptive. It answers the reviewer point about incomplete labels (Jiang et al. 2014).
4. **Abstract.** Say "16 predictors" instead of "16 methods" (literature package, 0.3 item 4).

---

## 6. Ideas considered and rejected

| Idea | Why not now |
|---|---|
| EnzPlacer's unseen-EC4 set as an external open-set test for CLEAN | 84% of its proteins are in CLEAN's training set (checked, section 3.1) |
| Open-world ("annotation delay") precision correction as a main idea | Studied for GO: Jiang et al. 2014; "The impact of incomplete knowledge on evaluation: an experimental benchmark for protein function prediction", *Bioinformatics* 25(18):2404 (2009); an "Open World Assumption" analysis repository. Kept only as quick win 3. |
| Taxonomy-stratified BLAST vs PLM as its own study | The EC-Bench test set is 434 eukaryotic vs 34 bacterial proteins, too unbalanced. SP-temporal needs lineage, which comes free with the E7 download; it is kept as an E7 secondary. |
| EC-Bench's "deep learning beats BLASTp at 30%" as a metric artefact | Checked in our E5 cards: at EC4 in the 30% regime, EnzBERT-learnt vs BLASTp is 0.463 vs 0.430 on weighted F1, 0.528 vs 0.509 on micro-F1, 0.165 vs 0.169 on union-macro F1. Mostly robust, so no story there. |
| Sparse-autoencoder interpretability of CLEAN failures | A 2026 preprint already applies SAEs to ESMC features for EC prediction (arXiv:2606.12209). Interpretability claims are hard to validate. |
| Chemistry-informed hard negatives, new losses, ESM-2 swaps | Need retraining (GPU) or are already published |
| Curation-allocation optimisation ("which proteins to curate first") | Speculative, with no validation data. Revisit after idea 3. |
| A living EC benchmark server | LAFA already runs continuous evaluation (GO-centred). Heavy maintenance for a solo project. |

---

## 7. Timeline and gates

| When | Work | Gate |
|---|---|---|
| Week 1 | E7-lite code and tests (here) → commit/register (your go-ahead) → download and run (your PC) → audit | **G7 and the outcome table in 2.4 decide C3's wording** |
| Week 1 | MDE row; open-world check; abstract fix | — |
| Week 2 | Submit the CBC paper | Literature package L3 items read by you |
| Weeks 2–4 | Idea 3: pre-register novelty detection and freeze thresholds; build the prospective pipeline | Freeze commit before 2026_04 is released |
| Weeks 3–8 | Paper 2: full E7 (EC-Bench, time analysis) + vintage matrix + claims table | Few overlaps found → provenance leads the paper |
| ~2027-03/04 | Idea 3 prospective evaluation at the pre-registered date | n ≥ 150 prospective enzymes, else report with the smaller n as pre-registered |

---

## 8. Guardrails carried over

**Data protocol:**
- New-392 is the development set; CAL is the only set for tuning and pilots.
- TEST is never used for choices; every new endpoint gets a dated pre-registration committed before TEST is scored.
- Price-149 is post hoc for anything designed after 2026-09-27.

**Metrics:**
- Classification endpoints use only micro-F1, union-macro F1 or per-protein F1.
- Support-weighted and true-label-macro scores are reported for comparability only, labelled as gameable.
- Paired bootstrap B = 2,000; Holm correction for per-tier or per-stratum tests.

**Code:**
- Follow cleanx conventions: paths from `cleanx/config.py`, argparse, `--smoke` writing `*_smoke` outputs, fixed seeds, `results/cleanx/<EXP>/` with `<EXP>_summary.md`.
- Every new script needs smoke and integration tests and entries in `cx_smoke.bat` and `cx_full.bat`.
- `cx_00_check.py` must stay at 9/9.
- Never copy CLEAN code.

**Licence and privacy:** the repository stays private. No CLEAN weights, embeddings or split CSVs go into public artefacts. UniProt data are CC BY 4.0; cite the release.

**Git:** nothing is committed or pushed without your explicit "commit and push".

---

## 9. Sources checked for this plan

These are search results, not full texts.

**Our own paper and data:**
- Our results: `results/cleanx/E2_split100/`, `E4/`, `E5/`, `data/cleanx/identity/`, `data/cleanx/sets/sp_temporal_meta.csv`.

**Evaluation practice and power:**
- EC-Bench: <https://academic.oup.com/bioinformaticsadvances/article/6/1/vbag004/8417619>
- Card et al. 2020: <https://aclanthology.org/2020.emnlp-main.745/>
- LiveProteinBench (contamination-free protein benchmark for LLMs): <https://arxiv.org/abs/2512.22257>

**Evidence codes and incomplete knowledge:**
- ECO:0000250 <https://evidenceontology.org/term/ECO:0000250/> and ECO:0000269 <https://evidenceontology.org/term/ECO:0000269/>
- Open World Assumption analysis: <https://github.com/arkatebi/OpenWorld-problem>
- Incomplete-knowledge benchmark, *Bioinformatics* 2009: <https://academic.oup.com/bioinformatics/article/25/18/2404/196123>

**Conformal and open-set methods:**
- Bates et al. 2023: <https://projecteuclid.org/journals/annals-of-statistics/volume-51/issue-1/Testing-for-outliers-with-conformal-p-values/10.1214/22-AOS2244.full>
- Xie et al. 2025: <https://arxiv.org/abs/2510.13037>
- Hierarchical conformal: <https://arxiv.org/abs/2508.13288>, <https://arxiv.org/abs/2501.19038>

**EC predictors and data sources:**
- EnzPlacer data (commit 97a8f0a): <https://github.com/drxiangma/EnzPlacer>
- SOLVE (New-392 removed from training): <https://pmc.ncbi.nlm.nih.gov/articles/PMC12424445/>
- Swiss-Prot release statistics (2026_03: 247 sequences added since 2026_02): <https://web.expasy.org/docs/relnotes/relstat.html>
- ESMC sparse-autoencoder EC preprint: <https://arxiv.org/abs/2606.12209>
