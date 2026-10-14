# Analysis protocols and historical notes

These module docstrings were relocated unchanged from commit `874abff997bc766e160f0b45869e21a4dbc3e187` during a
documentation-only cleanup. The original wording, dates and qualifications are preserved
below; this relocation is not a new registration or a revision of the recorded protocols.
Historical interpretations may differ from the final manuscript. Existing result digests
refer to the earlier code versions; no analyses were rerun for this cleanup.

<a id="cx-05-geometry"></a>

## `app/cx_05_geometry.py`

```text
E2: does the single-centroid-per-EC assumption limit CLEAN? Same frozen embeddings, different
inference geometry, all with CLEAN's max-separation rule:
  centroid   CLEAN (one mean per EC)
  exemplar   per-EC 1-NN
  kmean3     local mean of the 3 nearest members (LMKNN)
  multiproto silhouette-gated k-means sub-prototypes
  balanced1/balanced3  exemplar / kmean3 x n_c**lambda (equal-prior kNN density), lambda tuned on
             the calibration split only (--tune-on), never on a test set

Pre-registered endpoints (fixed 2026-09-27, before any sp_temporal_test result existed). Only
metrics that count every false positive are used (E1 shows support-weighted and true-label macro
scores rise when a method over-predicts):
  primary   macro-F1 over the union label space (EC depth 4), balanced3 vs centroid, sp_temporal_test,
            proteins whose ECs are all in CLEAN's vocabulary ('closed'), paired bootstrap, two-sided
  secondary micro-F1 for the same comparison; per-true-EC-tier micro-F1 differences

Reporting and sensitivity analyses added 2026-09-28, before any sp_temporal_test result existed (the primary and
secondary endpoints above are unchanged):
  - every delta row also carries delta_obs, the observed full-sample difference (reported as the effect; the
    bootstrap mean 'delta' is biased for macro-averaged metrics); per-tier rows carry p_holm (Holm across tiers)
  - closed|cluster_bootstrap_EC: the same comparisons with proteins resampled in clusters of their rarest true
    EC (proteins sharing an EC share a centroid, so the protein bootstrap is anti-conservative)
  - closed_capped<k>: at most k proteins per rarest true EC (seeded sample; --cap-per-ec, default 5)
  - closed_cal_disjoint: TEST proteins below 50% identity to every CAL protein (needs
    cx_04_identity.py --pair sp_temporal_test:sp_temporal_cal)
  - closed_PE1: temporal proteins with UniProt evidence at protein level
```

<a id="cx-06-ec-biology"></a>

## `app/cx_06_ec_biology.py`

```text
E3: why the centroid assumption fails - EC classes whose training members split into
non-homologous groups (no shared Pfam clan; candidate non-homologous isofunctional enzymes, NISE).

Per-EC geometry from CLEAN's training embeddings (split100 or split70):
  loo_misassign = share of members closer to another EC's centroid than to their own
                  leave-one-out centroid; rms_radius = within-EC spread.
Per-EC biology from current Swiss-Prot Pfam annotations of the same training accessions:
  clan groups = connected components of members linked by a shared Pfam clan (a family without
                a clan counts as its own clan); nise_like = >=2 groups each covering >=20% of the
                Pfam-annotated members.
E3b (--write-nise-fasta, cx_03_embed.py --sets nise_members, then --subgroups): are those
subgroups separate clusters in CLEAN space, i.e. is the single EC centroid between them?

Added 2026-09-28: mean_margin = mean(d_other - d_own) is not independent of rms_radius (a larger spread raises
d_own and so lowers the margin mechanically). The two components are therefore tested separately:
mean_d_own (distance of members to their own leave-one-out centroid: "dispersed") and mean_d_other (distance to
the nearest other, non-co-annotated centroid: "crowded").
```

<a id="cx-07-reliability"></a>

## `app/cx_07_reliability.py`

```text
E4: reliability. Is CLEAN's released confidence calibrated across rarity tiers, and can
tier-conditional (Mondrian) risk control give certified abstention per tier, including on
unseen-EC enzymes and non-enzymes?

Top-1 prediction per protein under each inference geometry (centroid = CLEAN; exemplar; balanced3
with the lambda E2 tuned on the calibration split).
- calibration: CLEAN's GMM confidence as released; isotonic maps fitted on the calibration split
  (global, and per tier of the predicted EC); ECE/Brier on the test sets
- risk control: SGR thresholds (Geifman & El-Yaniv 2017) on the raw distance, global vs one
  threshold per predicted-EC tier (delta split by Bonferroni); realized risk/coverage on test
- open set: share of unseen-EC enzymes and non-enzymes abstained at the certified thresholds;
  AUROC of each score for seen-EC enzymes vs non-enzymes
- hierarchical fallback: when EC4 abstains, emit the EC3 prefix if its own certified threshold passes

Changes recorded 2026-09-28, before any SP-temporal result existed:
- unseen-EC enzymes are split into 'pure' (no true EC in CLEAN's vocabulary: abstaining is the right answer) and
  'mixed' (at least one true EC is in the vocabulary, so a correct call exists; e.g. proteins annotated with both
  1.14.99.54 and 1.14.99.56); abstention is reported for both, plus accuracy among accepted mixed calls
- the EC3 fallback threshold is certified on the calibration proteins that the EC4 rule rejects, the population it
  is applied to (it was certified on all calibration proteins before, which does not cover the rejected subset)
Reporting added 2026-09-28 after the confirmatory run (no decision changes): every realized selective risk carries its
error count and an exact 95% Clopper-Pearson interval.
```

<a id="cx-09-evalcard"></a>

## `app/cx_09_evalcard.py`

```text
E5: does the reported EC metric change conclusions across the field? (evaluation card + EC-Bench re-scoring)

PRE-REGISTRATION (written 2026-09-28, before any EC-Bench metric other than the reported one was computed; the
only EC-Bench numbers computed before this text were the reproduction gate below).

Data. EC-Bench (Davoudi et al. 2026, doi:10.1093/bioadv/vbag004, MIT licence) publishes the predictions of 16 EC
predictors, commit 5e98a2c (cleanx/ecbench.py): its test set (468 proteins) and Price-149, for models trained on
its 100%- and 30%-identity clusterings. Labels are parsed exactly as in EC-Bench's own scoring code.
Gate. Our re-implementation must reproduce EC-Bench's printed support-weighted F1 (EC levels 1 and 4, test set,
100% regime, all 16 methods) to 1e-9; otherwise the run stops with FAIL.
Theory. Props. 1-4 in cleanx/evalcard.py: scores computed after projecting predictions onto the test label set,
and support-weighted / true-label-macro scores on any label space, cannot see false positives on EC numbers that
are absent from the evaluated labels; support-weighted recall equals micro recall.

H5 (primary). On the EC-Bench test set (100% regime, EC level 4) at least one pair of methods is ordered one way
by the reported metric (support-weighted F1, EC-Bench's scoring) and significantly the other way by micro-F1.
Primary endpoint: the number of such pairs among the 120, where "significantly" = two-sided paired-bootstrap p
(B = 2000 resamples of proteins, seed 20260927) with Benjamini-Hochberg adjustment across the 120 pairs,
q < 0.05, and the bootstrap median agrees in sign with the observed micro-F1 difference.
Secondary. (S1) the same count with union-label macro-F1 and with per-protein F1; (S2) Kendall tau-b between
the method rankings under the reported metric and each metric that counts every false positive, with a
bootstrap 95% CI; (S3) the share of each method's false positives that the reported metric cannot see;
(S4) Spearman rho between inflation (reported F1 - micro-F1) and mean predicted ECs per protein; (S5) S1-S4 for
the 30% regime, EC levels 1-3, and Price-149 (disclosed: Price-149 was used as the pilot, --pilot); (S6)
complete-EC-only sensitivity (tokens containing '-' dropped); (S7) after the full run, the invisible share of
CLEAN max-separation false positives on SP-temporal TEST, expected below New-392's 0.67 because its label set is
larger (Prop. 4).
Decision rules. >= 1 significant reversal: C1 becomes "the blind spot changes method comparisons in a published
benchmark". No reversal but tau-b < 1: "rankings largely preserved; reported precision inflated by X". tau-b = 1
and no reversal: report the benchmark as robust and keep C1 to CLEAN's own benchmarks and gameability.

AMENDMENT (2026-09-28, after the disclosed Price-149 pilot, before any scoring of the EC-Bench test set beyond
the gate). Pilot: rho(inflation, ECs/protein) = 0.94 (100% regime) and 0.93 (30%); 7 directional but no
significant reversals among 120 pairs (n = 149). Therefore S4 becomes the key secondary endpoint with a
directional hypothesis (rho > 0, one-sided alpha = 0.05, test set, 100% regime, EC level 4). The primary
endpoint is unchanged; the pilot suggests it may be underpowered, and a null result will be reported as such.

Modes: default (confirmatory, downloads the EC-Bench files at the pinned commit), --pilot (EC-Bench Price-149
only; writes E5_pilot), --smoke (synthetic EC-Bench-format data + New-392 CLEAN predictions; writes E5_smoke),
--card TRUTH PREDS (evaluation card for any prediction file; writes E5_card/E5_card_<name>.csv, outside the
confirmatory E5 folder).
```

<a id="cx-10-certbudget"></a>

## `app/cx_10_certbudget.py`

```text
E6: certification budget. How many calibration proteins does each rarity tier need before CLEAN's rare-EC
calls can be accepted with a certified selective risk, and how long does Swiss-Prot take to supply them?

PRE-REGISTRATION (written 2026-09-28, before any SP-temporal CAL or TEST feature existed).
Setting. Top-1 predictions of CLEAN (centroid geometry), score = -distance, groups = tier of the predicted EC,
Mondrian SGR with delta = 0.1 split over the 4 groups (as in cx_07_reliability.py). The budget n_needed(tier,
alpha) is the smallest number of calibration proteins in the tier for which the exact binomial power of certifying
selective risk alpha, at the most favourable coverage on a 0.1-0.9 grid and with SGR's per-test delta, reaches
0.8 (cleanx/certify.py). Selective risks at each coverage are estimated on SP-temporal CAL only, as the Jeffreys
posterior mean (errors + 0.5) / (accepted + 1), using only coverages with >= 10 accepted calibration proteins
(set after the disclosed pilot showed raw 0/2-0/5 error counts driving the budget). Because the coverage is
chosen with hindsight, n_needed is a lower bound on what SGR's own threshold search needs, which makes
"n_needed > n_available" a conservative claim. A tier with no usable coverage is reported as not estimable,
distinct from not certifiable (every usable coverage has risk >= alpha).
H6 (primary). For alpha = 0.2, n_needed exceeds the number of CAL proteins available in the very-rare and in the
singleton predicted tiers. Primary endpoint: the calibration deficit (n_needed - n_available) per predicted tier,
with a bootstrap 95% CI (1000 resamples of the tier's CAL proteins; resamples that cannot certify count as
infinite).
Secondary. (a) alpha = 0.1 and 0.3; (b) exemplar and balanced3 geometries; (c) years of Swiss-Prot accrual the
deficit represents: yearly number of new Swiss-Prot enzymes passing the benchmark filters, 2022-03-01 to the
UniProt release date, times the tier's share of CAL predictions; (d) the universal curve n*(r, alpha, delta');
(e) a simulation check that SGR certifies at the budgeted size with about the targeted probability; (f) after E4,
the realized per-tier certification on TEST compared with the budget's prediction.
Decision rules. Deficit CI above 0 in the rare tiers: C4 reports "rare-EC calls cannot be certified with the
calibration data that exist; they need N proteins, about Y years of curation, or tier pooling / hierarchical
abstention". Deficit <= 0: certification is feasible with current data (G-B expected positive).

Reporting changes 2026-09-28, after the confirmatory run (no estimate or decision changes): a bootstrap quantile that
falls between two non-certifiable (infinite) resamples is reported as inf (numpy's interpolation returned NaN); the
years of accrual at the deficit's lower 95% bound are added (years_at_deficit_ci_lo), since the point estimate of the
deficit is unstable when a tier is close to not certifiable.

Modes: default (confirmatory; needs results/cleanx/E4/E4_features_cal_<geometry>.csv from the full run),
--pilot (New-392 half used as calibration in the E4 smoke run; disclosed pilot; writes E6_pilot),
--smoke (synthetic features and known-answer checks; writes E6_smoke).
```

<a id="cx-13-provenance"></a>

## `app/cx_13_provenance.py`

```text
E7: label provenance. Is BLASTp's advantage over CLEAN on SP-temporal a property of the methods, or partly of how
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
CONFIRMATORY).
```
