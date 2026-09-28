# E6 summary (delta=0.1, target power 0.8, source sp_temporal_cal)

- centroid common alpha=0.1: available 525, top-1 acc 0.61; needed 15400 (coverage 0.2, selective risk 0.08); deficit 14875 [CI 455, inf]; ~51.0 years of Swiss-Prot accrual (>= 1.6 at the lower bound)
- centroid common alpha=0.2: available 525, top-1 acc 0.61; needed 635 (coverage 0.2, selective risk 0.08); deficit 110 [CI -250, 2055]; ~0.4 years of Swiss-Prot accrual (>= 0.0 at the lower bound)
- centroid common alpha=0.3: available 525, top-1 acc 0.61; needed 225 (coverage 0.2, selective risk 0.08); deficit -300 [CI -380, -50]; ~0.0 years of Swiss-Prot accrual (>= 0.0 at the lower bound)
- centroid rare alpha=0.1: available 122, top-1 acc 0.36; risk above alpha (not certifiable at any budget); deficit CI [2248, inf] (>= 33.2 years of accrual at the lower bound); bootstrap share of certifiable resamples 0.13
- centroid rare alpha=0.2: available 122, top-1 acc 0.36; needed 49340 (coverage 0.1, selective risk 0.18); deficit 49218 [CI 428, inf]; ~726.7 years of Swiss-Prot accrual (>= 6.3 at the lower bound)
- centroid rare alpha=0.3: available 122, top-1 acc 0.36; needed 1820 (coverage 0.1, selective risk 0.18); deficit 1698 [CI 178, inf]; ~25.1 years of Swiss-Prot accrual (>= 2.6 at the lower bound)
- centroid very_rare alpha=0.1: available 102, top-1 acc 0.25; risk above alpha (not certifiable at any budget); deficit CI [2948, inf] (>= 52.1 years of accrual at the lower bound); bootstrap share of certifiable resamples 0.15
- centroid very_rare alpha=0.2: available 102, top-1 acc 0.25; risk above alpha (not certifiable at any budget); deficit CI [528, inf] (>= 9.3 years of accrual at the lower bound); bootstrap share of certifiable resamples 0.46
- centroid very_rare alpha=0.3: available 102, top-1 acc 0.25; needed 3340 (coverage 0.1, selective risk 0.21); deficit 3238 [CI 198, inf]; ~57.2 years of Swiss-Prot accrual (>= 3.5 at the lower bound)
- centroid singleton alpha=0.1: available 90, top-1 acc 0.19; risk above alpha (not certifiable at any budget); deficit CI [13515, inf] (>= 270.5 years of accrual at the lower bound); bootstrap share of certifiable resamples 0.04
- centroid singleton alpha=0.2: available 90, top-1 acc 0.19; risk above alpha (not certifiable at any budget); deficit CI [515, inf] (>= 10.3 years of accrual at the lower bound); bootstrap share of certifiable resamples 0.21
- centroid singleton alpha=0.3: available 90, top-1 acc 0.19; risk above alpha (not certifiable at any budget); deficit CI [135, inf] (>= 2.7 years of accrual at the lower bound); bootstrap share of certifiable resamples 0.42
- accrual: 466 new Swiss-Prot enzymes/year pass the benchmark filters
- SGR simulation common alpha=0.2: certification rate 0.56 at the oracle budget (635), 0.92 at twice it (target 0.8)
- SGR simulation common alpha=0.3: certification rate 0.58 at the oracle budget (225), 0.91 at twice it (target 0.8)
- SGR simulation rare alpha=0.3: certification rate 0.55 at the oracle budget (1820), 0.93 at twice it (target 0.8)
- SGR simulation very_rare alpha=0.3: certification rate 0.53 at the oracle budget (3340), 0.91 at twice it (target 0.8)
