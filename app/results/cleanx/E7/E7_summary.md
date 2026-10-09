# E7 summary (sp_temporal_test; confirmatory @8b050b63367a)

- code commit 8b050b63367aadab977f43a379873f8161756425; analysis code SHA-256 cx_13_provenance.py b20c13233437, cleanx/provenance.py 2ba6a87b1473, cleanx/contrast.py 6e9cda7466e3, cleanx/metrics.py 130a00d469bc, cleanx/freeze.py bead0e9751a6
- UniProt release(s) of the evidence: ['2026_03']; B = 2000, seed 20260927

- verdict: **BLASTp's advantage holds on experimentally supported labels; no detectable provenance effect (DiD +0.012 [-0.037, +0.059], MDE80 0.068)**
- gate G7 (n_EXP >= 100, n_PROPAGATED >= 50): passed; closed proteins by class: {'EXP': 477, 'SIM': 218, 'RULE': 61, 'OTHER': 120, 'LABEL_ABSENT': 0}
- P7a (delta_EXP, micro_F1, n0=477, n1=279): +0.082 [+0.052, +0.114] p=<0.0005 q_holm=<0.0005 MDE80=0.044
- P7b (DiD, micro_F1, n0=477, n1=279): +0.012 [-0.037, +0.059] p=0.607 q_holm=0.607 MDE80=0.068
- S7a_SIM (DiD, micro_F1, n0=477, n1=218): +0.029 [-0.023, +0.078] p=0.289 q_holm=0.289 MDE80=0.073
- S7a_RULE (DiD, micro_F1, n0=477, n1=61): -0.058 [-0.117, +0.006] p=0.067 q_holm=0.134 MDE80=0.089
- S7b_delta_EXP (delta_EXP_standardised, micro_F1, n0=477, n1=279): +0.078 [+0.050, +0.108] p=<0.0005 MDE80=0.041 (bins <30% (w=0.188), 30-50% (w=0.329), 50-70% (w=0.208), 70-90% (w=0.169), >=90% (w=0.106))
- S7b (DiD_standardised, micro_F1, n0=477, n1=279): +0.017 [-0.032, +0.066] p=0.489 MDE80=0.071 (bins <30% (w=0.188), 30-50% (w=0.329), 50-70% (w=0.208), 70-90% (w=0.169), >=90% (w=0.106))
- S7c_delta_EXP (delta_EXP, samples_F1, n0=477, n1=279): +0.077 [+0.046, +0.108] p=<0.0005 MDE80=0.045
- S7c (DiD, samples_F1, n0=477, n1=279): +0.004 [-0.043, +0.053] p=0.856 MDE80=0.067
- S7d_delta_EXP (delta_EXP, micro_F1, n0=98, n1=46): +0.065 [+0.001, +0.136] p=0.049 MDE80=0.096
- S7d (DiD, micro_F1, n0=98, n1=46): -0.086 [-0.205, +0.018] p=0.108 MDE80=0.157
- sensitivity sp_temporal_test closed | cluster_bootstrap_rarest_EC: P7a +0.082 [+0.050, +0.113] p=<0.0005
- sensitivity sp_temporal_test closed | cluster_bootstrap_rarest_EC: P7b +0.012 [-0.036, +0.061] p=0.654
- sensitivity sp_temporal_test closed | capped5_per_rarest_EC: P7a +0.081 [+0.049, +0.112] p=<0.0005
- sensitivity sp_temporal_test closed | capped5_per_rarest_EC: P7b +0.019 [-0.029, +0.072] p=0.436
- M7a top hit is a cited source: 45/218 = 0.206 [0.155, 0.266]
- M7b a cited source is among the kept hits: 128/218 = 0.587 [0.519, 0.653]
- M7c a cited source is a split100 training protein: 164/218 = 0.752 [0.689, 0.808]
- M7d top hit is a cited source: n=45, BLASTp 0.942, CLEAN 0.796
- M7d top hit is not a cited source: n=173, BLASTp 0.853, CLEAN 0.751
