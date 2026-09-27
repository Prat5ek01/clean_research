# E3 summary (split100)

- ECs analysed (n>=5, >=80% of members with Pfam): 2493; NISE-like (>=2 non-homologous Pfam-clan groups): 183 (7.3%)
- training-space LOO centroid misassignment, NISE-like vs other: median 0.000 vs 0.000, mean 0.020 vs 0.015, Mann-Whitney p=0.0356, AUROC 0.536; size-stratified permutation p=0.08092
- Spearman(top Pfam-architecture share, misassignment) rho=-0.227
- E3b: no NISE-like EC has all member embeddings yet (run --write-nise-fasta, then cx_03_embed.py)
