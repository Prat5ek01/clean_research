# E4 summary (smoke: New-392 halves)

geometries: ['centroid', 'exemplar', 'balanced3']; lambda: {'balanced1': 0.02, 'balanced3': 0.02}

- new392_half ECE by predicted-EC tier, centroid/clean_gmm(raw): very_rare 0.495 (acc 0.34, conf 0.84, n=32), common 0.179 (acc 0.64, conf 0.64, n=107), rare 0.359 (acc 0.48, conf 0.74, n=46), singleton 0.490 (acc 0.36, conf 0.85, n=11)
- new392_half ECE by predicted-EC tier, centroid/mondrian: very_rare 0.296 (acc 0.34, conf 0.64, n=32), common 0.114 (acc 0.64, conf 0.64, n=107), rare 0.129 (acc 0.48, conf 0.55, n=46), singleton 0.349 (acc 0.36, conf 0.68, n=11)
- new392_half ECE by predicted-EC tier, exemplar/mondrian: common 0.095 (acc 0.67, conf 0.63, n=132), rare 0.305 (acc 0.49, conf 0.73, n=41), very_rare 0.298 (acc 0.38, conf 0.67, n=16)
- new392_half ECE by predicted-EC tier, balanced3/mondrian: common 0.096 (acc 0.68, conf 0.64, n=124), rare 0.267 (acc 0.58, conf 0.58, n=36), very_rare 0.281 (acc 0.33, conf 0.61, n=21), singleton 0.428 (acc 0.27, conf 0.63, n=15)
- risk control alpha=0.2 global/centroid: all risk nan cov 0.00, pred_tier:common risk nan cov 0.00, pred_tier:rare risk nan cov 0.00, pred_tier:very_rare risk nan cov 0.00, pred_tier:singleton risk nan cov 0.00
- risk control alpha=0.2 global/exemplar: all risk nan cov 0.00, pred_tier:common risk nan cov 0.00, pred_tier:rare risk nan cov 0.00, pred_tier:very_rare risk nan cov 0.00, pred_tier:singleton risk nan cov 0.00
- risk control alpha=0.2 global/balanced3: all risk nan cov 0.00, pred_tier:common risk nan cov 0.00, pred_tier:rare risk nan cov 0.00, pred_tier:very_rare risk nan cov 0.00, pred_tier:singleton risk nan cov 0.00
- risk control alpha=0.2 mondrian/centroid: all risk nan cov 0.00, pred_tier:common risk nan cov 0.00, pred_tier:rare risk nan cov 0.00, pred_tier:very_rare risk nan cov 0.00, pred_tier:singleton risk nan cov 0.00
- risk control alpha=0.2 mondrian/exemplar: all risk nan cov 0.00, pred_tier:common risk nan cov 0.00, pred_tier:rare risk nan cov 0.00, pred_tier:very_rare risk nan cov 0.00, pred_tier:singleton risk nan cov 0.00
- risk control alpha=0.2 mondrian/balanced3: all risk nan cov 0.00, pred_tier:common risk nan cov 0.00, pred_tier:rare risk nan cov 0.00, pred_tier:very_rare risk nan cov 0.00, pred_tier:singleton risk nan cov 0.00
- smallest certifiable selective risk (calibration split), centroid: all (global threshold) 0.4 (cal coverage 0.24, n=196), pred_tier:common 0.35 (cal coverage 0.76, n=106), pred_tier:rare 0.6 (cal coverage 0.16, n=45), pred_tier:very_rare 0.8 (cal coverage 0.97, n=34), pred_tier:singleton 0.95 (cal coverage 0.73, n=11)
- smallest certifiable selective risk (calibration split), exemplar: all (global threshold) 0.3 (cal coverage 0.62, n=196), pred_tier:common 0.3 (cal coverage 0.63, n=130), pred_tier:rare 0.75 (cal coverage 0.97, n=40), pred_tier:very_rare 0.65 (cal coverage 0.50, n=20), pred_tier:singleton 0.8 (cal coverage 0.50, n=6)
- smallest certifiable selective risk (calibration split), balanced3: all (global threshold) 0.3 (cal coverage 0.23, n=196), pred_tier:common 0.25 (cal coverage 0.25, n=116), pred_tier:rare 0.7 (cal coverage 0.97, n=38), pred_tier:very_rare 0.8 (cal coverage 0.62, n=29), pred_tier:singleton 0.95 (cal coverage 0.62, n=13)
