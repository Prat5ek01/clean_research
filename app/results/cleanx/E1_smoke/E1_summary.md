# E1 summary (split100, bootstrap B=200)

- **new392** (n=392): reported (support-weighted) P/R/F1 = 0.596/0.479/0.497; micro P/R/F1 = 0.518/0.479/0.498; macro(true) F1 = 0.533; 150/224 FPs invisible to the reported precision; exact training duplicates: 5
- over-prediction stress test new392: top-1 -> top-10 ECs per protein changes reported P 0.560 -> 0.482 but micro P 0.541 -> 0.093
- label drift new392: same=385, EC_changed=3, EC_added=2, no_longer_has_EC_or_absent=1, EC_removed=1
- consistency with the earlier audit (audit_out/maxsep_tier_x_depth.csv): max |F1 diff| = 1.11e-16 over 16 cells
