# Selection and rescore results (metrics only, no sequences or design ids)

Run 2026-11-01. All numbers are predictions (Boltz-2 refold + ipSAE), not measurements. Costs are controller or script estimates.

## Pool
- `kras_set5` top-up (cap $60): 14 shards, 74 designs refolded, est $53.57, 0 designs passed the tier A/B gates.
- Pooled with `kras_set2` (32 designs): 106 scored designs, 65 unique sequences with G12D ipSAE >= 0.5.

## Method
1. Rank by `delta` = G12D ipSAE (mean of 3 samples) minus the higher of WT and peptide-free (max of 3 samples). The top 20 were the same with no identity filter, an 80% filter and a 95% filter.
2. Rescore the top 20 alone on all 3 targets with 8 samples each (60 folds, 8.8 min, est $2.72): blocks 10-11.
3. Rescore the other 45 the same way (135 folds, expected ~$6): blocks 12-13.
4. Merge and rank all 65 by score = (delta_mean + delta_worst) / 2. delta_mean uses the mean avoid score, delta_worst the worst avoid sample. pep_delta = G12D peptide-chain ipSAE minus WT's.

## Final top 20 (8-sample means)
| rank | len | G12D ipSAE | WT ipSAE | PF ipSAE | delta_mean | delta_worst | pep_delta | Asp6 min (A) |
|---|---|---|---|---|---|---|---|---|
| 1 | 78 | 0.736 | 0.279 | 0.170 | 0.457 | 0.383 | 0.361 | 6.47 |
| 2 | 66 | 0.890 | 0.465 | 0.388 | 0.425 | 0.291 | 0.485 | 6.81 |
| 3 | 98 | 0.846 | 0.558 | 0.239 | 0.287 | 0.180 | 0.292 | 13.15 |
| 4 | 66 | 0.912 | 0.527 | 0.591 | 0.321 | 0.035 | 0.448 | 5.90 |
| 5 | 66 | 0.916 | 0.662 | 0.346 | 0.254 | 0.096 | 0.288 | 6.70 |
| 6 | 66 | 0.926 | 0.573 | 0.733 | 0.193 | 0.033 | 0.410 | 6.21 |
| 7 | 62 | 0.761 | 0.558 | 0.314 | 0.203 | 0.020 | 0.407 | 5.16 |
| 8 | 64 | 0.721 | 0.444 | 0.294 | 0.277 | -0.091 | 0.323 | 5.09 |
| 9 | 62 | 0.647 | 0.522 | 0.324 | 0.124 | 0.052 | 0.100 | 9.46 |
| 10 | 66 | 0.912 | 0.766 | 0.337 | 0.146 | 0.021 | 0.247 | 6.48 |
| 11 | 68 | 0.944 | 0.849 | 0.860 | 0.084 | 0.032 | 0.154 | 5.99 |
| 12 | 99 | 0.657 | 0.513 | 0.359 | 0.144 | -0.044 | 0.089 | 7.06 |
| 13 | 64 | 0.604 | 0.380 | 0.412 | 0.192 | -0.101 | 0.170 | 5.11 |
| 14 | 64 | 0.967 | 0.914 | 0.628 | 0.052 | 0.036 | 0.141 | 6.55 |
| 15 | 99 | 0.975 | 0.563 | 0.935 | 0.040 | 0.028 | 0.540 | 5.95 |
| 16 | 61 | 0.951 | 0.424 | 0.906 | 0.044 | 0.017 | 0.422 | 8.62 |
| 17 | 64 | 0.962 | 0.924 | 0.374 | 0.038 | -0.005 | 0.127 | 6.54 |
| 18 | 64 | 0.969 | 0.951 | 0.776 | 0.018 | 0.007 | 0.086 | 9.09 |
| 19 | 73 | 0.950 | 0.940 | 0.663 | 0.010 | 0.001 | 0.054 | 9.08 |
| 20 | 79 | 0.952 | 0.942 | 0.571 | 0.010 | -0.005 | 0.074 | 6.25 |

Reading: only ranks 1 and 2 show a gap that holds on the worst sample and a positive peptide gap. Rank 3 does not touch the peptide (Asp6 13 A). Ranks 15 and 16 have a large peptide gap but bind peptide-free HLA at ~0.9, so they are not peptide-dependent. Ranks 14-20 are close to pan-HLA binders. No design passed the strict gates.
