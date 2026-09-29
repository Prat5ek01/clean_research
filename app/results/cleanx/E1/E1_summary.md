# E1 summary (split100, bootstrap B=2000)

- **new392** (n=392): reported (support-weighted) P/R/F1 = 0.596/0.479/0.497; micro P/R/F1 = 0.518/0.479/0.498; macro(true) F1 = 0.533; 150/224 FPs invisible to the reported precision; exact training duplicates: 5
- label drift price149: not applicable (ids are not Swiss-Prot accessions)
- **price149** (n=149): reported (support-weighted) P/R/F1 = 0.584/0.467/0.495; micro P/R/F1 = 0.392/0.467/0.426; macro(true) F1 = 0.447; 106/110 FPs invisible to the reported precision; exact training duplicates: 0
- **sp_temporal_test** (n=1260): reported (support-weighted) P/R/F1 = 0.459/0.451/0.443; micro P/R/F1 = 0.441/0.451/0.446; macro(true) F1 = 0.404; 616/830 FPs invisible to the reported precision; exact training duplicates: 0
- over-prediction stress test new392: top-1 -> top-10 ECs per protein changes reported P 0.560 -> 0.482 but micro P 0.541 -> 0.093
- over-prediction stress test price149: top-1 -> top-10 ECs per protein changes reported P 0.519 -> 0.636 but micro P 0.403 -> 0.070
- over-prediction stress test sp_temporal_test: top-1 -> top-10 ECs per protein changes reported P 0.453 -> 0.279 but micro P 0.479 -> 0.067
- label drift new392: same=385, EC_changed=3, EC_added=2, no_longer_has_EC_or_absent=1, EC_removed=1
- consistency with the earlier audit (audit_out/maxsep_tier_x_depth.csv): max |F1 diff| = 1.11e-16 over 16 cells
