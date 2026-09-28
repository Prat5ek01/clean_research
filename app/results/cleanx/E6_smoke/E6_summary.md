# E6_smoke summary (delta=0.1, target power 0.8, source synthetic)

Synthetic features: code-path test only, no scientific content.

- centroid common alpha=0.1: available 481, top-1 acc 0.77; needed 1930 (coverage 0.1, selective risk 0.03); deficit 1449 [CI 380, 12131]; ~5.2 years of Swiss-Prot accrual
- centroid common alpha=0.2: available 481, top-1 acc 0.77; needed 354 (coverage 0.3, selective risk 0.07); deficit -127 [CI -301, 221]; ~0.0 years of Swiss-Prot accrual
- centroid common alpha=0.3: available 481, top-1 acc 0.77; needed 148 (coverage 0.4, selective risk 0.10); deficit -333 [CI -401, -269]; ~0.0 years of Swiss-Prot accrual
- centroid rare alpha=0.1: available 159, top-1 acc 0.60; risk above alpha; bootstrap share of certifiable resamples 0.22
- centroid rare alpha=0.2: available 159, top-1 acc 0.60; needed 9865 (coverage 0.2, selective risk 0.17); deficit 9706 [CI 93, inf]; ~104.8 years of Swiss-Prot accrual
- centroid rare alpha=0.3: available 159, top-1 acc 0.60; needed 573 (coverage 0.4, selective risk 0.19); deficit 414 [CI 21, 7504]; ~4.5 years of Swiss-Prot accrual
- centroid very_rare alpha=0.1: available 71, top-1 acc 0.37; risk above alpha; bootstrap share of certifiable resamples 0.02
- centroid very_rare alpha=0.2: available 71, top-1 acc 0.37; risk above alpha; bootstrap share of certifiable resamples 0.08
- centroid very_rare alpha=0.3: available 71, top-1 acc 0.37; risk above alpha; bootstrap share of certifiable resamples 0.42
- centroid singleton alpha=0.1: available 89, top-1 acc 0.44; risk above alpha; bootstrap share of certifiable resamples 0.04
- centroid singleton alpha=0.2: available 89, top-1 acc 0.44; risk above alpha; bootstrap share of certifiable resamples 0.20
- centroid singleton alpha=0.3: available 89, top-1 acc 0.44; risk above alpha; bootstrap share of certifiable resamples 0.60
- accrual: 466 new Swiss-Prot enzymes/year pass the benchmark filters
- SGR simulation common alpha=0.1: certification rate 0.50 at the oracle budget (1930), 0.96 at twice it (target 0.8)
- SGR simulation common alpha=0.2: certification rate 0.69 at the oracle budget (354), 1.00 at twice it (target 0.8)
- SGR simulation common alpha=0.3: certification rate 0.54 at the oracle budget (148), 0.95 at twice it (target 0.8)
- SGR simulation rare alpha=0.3: certification rate 0.48 at the oracle budget (573), 0.94 at twice it (target 0.8)
