# E7 summary (new392; PILOT on New-392 (development set; descriptive))

- code commit 8b050b63367aadab977f43a379873f8161756425; analysis code SHA-256 cx_13_provenance.py b20c13233437, cleanx/provenance.py 2ba6a87b1473, cleanx/contrast.py 6e9cda7466e3, cleanx/metrics.py 130a00d469bc, cleanx/freeze.py bead0e9751a6
- UniProt release(s) of the evidence: ['2026_03']; B = 2000, seed 20260927

- verdict: **not confirmatory: BLASTp's advantage holds on experimentally supported labels; no detectable provenance effect (DiD +0.060 [-0.016, +0.131], MDE80 0.107)**
- gate G7 (n_EXP >= 100, n_PROPAGATED >= 50): passed; closed proteins by class: {'EXP': 198, 'SIM': 133, 'RULE': 11, 'OTHER': 45, 'LABEL_ABSENT': 5}
- P7a (delta_EXP, micro_F1, n0=198, n1=144): +0.075 [+0.024, +0.127] p=0.003 q_holm=0.006 MDE80=0.074
- P7b (DiD, micro_F1, n0=198, n1=144): +0.060 [-0.016, +0.131] p=0.13 q_holm=0.13 MDE80=0.107
- S7a_SIM (DiD, micro_F1, n0=198, n1=133): +0.056 [-0.017, +0.127] p=0.153 q_holm=0.306 MDE80=0.103
- S7a_RULE (DiD, micro_F1, n0=198, n1=11): +0.099 [-0.095, +0.358] p=0.443 q_holm=0.443 MDE80=0.332
- S7b_delta_EXP (delta_EXP_standardised, micro_F1, n0=198, n1=144): +0.084 [+0.031, +0.141] p=0.002 MDE80=0.077 (bins <30% (w=0.146), 30-50% (w=0.292), 50-70% (w=0.260), 70-90% (w=0.199), >=90% (w=0.102))
- S7b (DiD_standardised, micro_F1, n0=198, n1=144): +0.033 [-0.038, +0.109] p=0.376 MDE80=0.105 (bins <30% (w=0.146), 30-50% (w=0.292), 50-70% (w=0.260), 70-90% (w=0.199), >=90% (w=0.102))
- S7c_delta_EXP (delta_EXP, samples_F1, n0=198, n1=144): +0.111 [+0.063, +0.165] p=<0.0005 MDE80=0.073
- S7c (DiD, samples_F1, n0=198, n1=144): -0.012 [-0.082, +0.054] p=0.693 MDE80=0.100
- sensitivity new392 closed | cluster_bootstrap_rarest_EC: P7a +0.075 [+0.012, +0.144] p=0.015
- sensitivity new392 closed | cluster_bootstrap_rarest_EC: P7b +0.060 [-0.039, +0.165] p=0.254
- sensitivity new392 closed | capped5_per_rarest_EC: P7a +0.077 [+0.029, +0.126] p=0.002
- sensitivity new392 closed | capped5_per_rarest_EC: P7b +0.068 [-0.012, +0.146] p=0.09
- M7a top hit is a cited source: 22/133 = 0.165 [0.107, 0.240]
- M7b a cited source is among the kept hits: 48/133 = 0.361 [0.279, 0.449]
- M7c a cited source is a split100 training protein: 62/133 = 0.466 [0.379, 0.555]
- M7d top hit is a cited source: n=22, BLASTp 0.958, CLEAN 0.809
- M7d top hit is not a cited source: n=111, BLASTp 0.605, CLEAN 0.481
