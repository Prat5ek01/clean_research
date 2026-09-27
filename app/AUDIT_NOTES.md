# CLEAN Audit - Reproduction & Training-Set Evidence
_Generated: 2026-07-22 23:43_

## 1. Reproduction of published CLEAN results (New-392, split100 pretrained)

| Algorithm | Metric | README | Mine |
|---|---|---|---|
| max-separation | precision | 0.596 | 0.596 |
| max-separation | recall | 0.479 | 0.479 |
| max-separation | F1 | 0.497 | 0.497 |
| max-separation | AUC | 0.739 | 0.739 |
| p-value (1e-5, 20k) | precision | 0.558 | 0.558 |
| p-value (1e-5, 20k) | recall | 0.477 | 0.477 |
| p-value (1e-5, 20k) | F1 | 0.482 | 0.482 |
| p-value (1e-5, 20k) | AUC | 0.737 | 0.737 |

Both algorithms reproduced exactly. 392 samples | 177 EC | 5242 cluster centers.

## 2. Training-set (split100.csv) frequency analysis

- Unique EC numbers: **5242** (matches 5242 inference cluster centers)
- Total (EC, sequence) pairs: **241025** (matches 241025 loaded embeddings)
- Seen exactly once: **1360** (25.9%)
- Seen <=3 times: **2415** (46.1%)
- Seen <=10 times: **3520** (67.1%)
- Median times seen: **4**
- Max times seen: **2158**

Caveat: counts are for the released split100.csv, which the README notes is not
guaranteed identical to the data behind the pretrained weights.