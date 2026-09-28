# E3 summary (split100)

- ECs analysed (n>=5, >=80% of members with Pfam): 2493; NISE-like (>=2 non-homologous Pfam-clan groups): 183 (7.3%); by EC class {'1': 55, '2': 38, '3': 54, '4': 15, '5': 9, '6': 6, '7': 6}
- within-EC spread (RMS radius), NISE-like vs other: median 4.138 vs 3.435, Mann-Whitney p=1.06e-40, AUROC 0.796; size-stratified permutation (mean diff +0.839) p=9.999e-05
- margin to nearest other centroid, NISE-like vs other: median 4.212 vs 5.258, Mann-Whitney p=5.52e-26, AUROC 0.266; size-stratified permutation (mean diff -1.001) p=9.999e-05
- LOO centroid misassignment, NISE-like vs other: median 0.000 vs 0.000, Mann-Whitney p=0.0356, AUROC 0.536; size-stratified permutation (mean diff +0.005) p=0.1371
- distance to own LOO centroid (dispersed), NISE-like vs other: median 4.368 vs 3.538, Mann-Whitney p=1.51e-46, AUROC 0.818; size-stratified permutation (mean diff +0.920) p=9.999e-05
- distance to nearest other centroid (crowded), NISE-like vs other: median 8.639 vs 8.754, Mann-Whitney p=0.133, AUROC 0.467; size-stratified permutation (mean diff -0.081) p=0.3644
- Spearman(top Pfam-architecture share, misassignment) rho=-0.227
- new392: exemplar-minus-centroid per-protein F1 gain, true EC NISE-like 0.034666666666666665 (n=25) vs other 0.04064054322909909 (n=367), Mann-Whitney p=0.5762531764717762
- price149: exemplar-minus-centroid per-protein F1 gain, true EC NISE-like None (n=0) vs other -0.016554809843400447 (n=149), Mann-Whitney p=None
- sp_temporal_test: exemplar-minus-centroid per-protein F1 gain, true EC NISE-like 0.07252326217843459 (n=87) vs other 0.061608277364312944 (n=787), Mann-Whitney p=0.35920334015892064
