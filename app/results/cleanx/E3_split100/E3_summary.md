# E3 summary (split100)

- ECs analysed (n>=5, >=80% of members with Pfam): 2493; NISE-like (>=2 non-homologous Pfam-clan groups): 183 (7.3%)
- training-space LOO centroid misassignment, NISE-like vs other: median 0.000 vs 0.000, mean 0.020 vs 0.015, Mann-Whitney p=0.0356, AUROC 0.536; size-stratified permutation p=0.07799
- Spearman(top Pfam-architecture share, misassignment) rho=-0.227
- new392: exemplar-minus-centroid per-protein F1 gain, true EC NISE-like 0.034666666666666665 (n=25) vs other 0.04064054322909909 (n=367), Mann-Whitney p=0.5762531764717762
- price149: exemplar-minus-centroid per-protein F1 gain, true EC NISE-like None (n=0) vs other -0.016554809843400447 (n=149), Mann-Whitney p=None
