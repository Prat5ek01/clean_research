# E6_pilot summary (delta=0.1, target power 0.8, source pilot_new392_half)

- centroid common alpha=0.1: available 106, top-1 acc 0.64; needed 3050 (coverage 0.1, selective risk 0.04); deficit 2944 [CI 628, inf]; ~11.7 years of Swiss-Prot accrual (>= 2.5 at the lower bound)
- centroid common alpha=0.2: available 106, top-1 acc 0.64; needed 630 (coverage 0.1, selective risk 0.04); deficit 524 [CI 34, 6744]; ~2.1 years of Swiss-Prot accrual (>= 0.1 at the lower bound)
- centroid common alpha=0.3: available 106, top-1 acc 0.64; needed 200 (coverage 0.7, selective risk 0.16); deficit 94 [CI -41, 474]; ~0.4 years of Swiss-Prot accrual (>= 0.0 at the lower bound)
- centroid rare alpha=0.1: available 45, top-1 acc 0.47; risk above alpha (not certifiable at any budget); deficit CI [inf, inf] (>= nan years of accrual at the lower bound); bootstrap share of certifiable resamples 0.01
- centroid rare alpha=0.2: available 45, top-1 acc 0.47; risk above alpha (not certifiable at any budget); deficit CI [595, inf] (>= 5.6 years of accrual at the lower bound); bootstrap share of certifiable resamples 0.09
- centroid rare alpha=0.3: available 45, top-1 acc 0.47; risk above alpha (not certifiable at any budget); deficit CI [314, inf] (>= 2.9 years of accrual at the lower bound); bootstrap share of certifiable resamples 0.23
- centroid very_rare alpha=0.1: available 34, top-1 acc 0.44; risk above alpha (not certifiable at any budget); deficit CI [inf, inf] (>= nan years of accrual at the lower bound); bootstrap share of certifiable resamples 0.00
- centroid very_rare alpha=0.2: available 34, top-1 acc 0.44; risk above alpha (not certifiable at any budget); deficit CI [1153, inf] (>= 14.3 years of accrual at the lower bound); bootstrap share of certifiable resamples 0.05
- centroid very_rare alpha=0.3: available 34, top-1 acc 0.44; risk above alpha (not certifiable at any budget); deficit CI [240, inf] (>= 3.0 years of accrual at the lower bound); bootstrap share of certifiable resamples 0.29
- centroid singleton alpha=0.1: available 11, top-1 acc 0.27; risk above alpha (not certifiable at any budget); deficit CI [inf, inf] (>= nan years of accrual at the lower bound); bootstrap share of certifiable resamples 0.00
- centroid singleton alpha=0.2: available 11, top-1 acc 0.27; risk above alpha (not certifiable at any budget); deficit CI [inf, inf] (>= nan years of accrual at the lower bound); bootstrap share of certifiable resamples 0.00
- centroid singleton alpha=0.3: available 11, top-1 acc 0.27; risk above alpha (not certifiable at any budget); deficit CI [inf, inf] (>= nan years of accrual at the lower bound); bootstrap share of certifiable resamples 0.00
- accrual: 466 new Swiss-Prot enzymes/year pass the benchmark filters
- SGR simulation common alpha=0.1: certification rate 0.42 at the oracle budget (3050), 0.92 at twice it (target 0.8)
- SGR simulation common alpha=0.2: certification rate 0.52 at the oracle budget (630), 0.97 at twice it (target 0.8)
- SGR simulation common alpha=0.3: certification rate 0.62 at the oracle budget (200), 0.97 at twice it (target 0.8)
