# Run log (honest record)

Costs are the controller's own estimates unless stated. Counts are what was reported at the time.

| Tag | What it was | Result |
|---|---|---|
| first smoke | Original smoke entrypoint | Killed: it called `read_progress` after `spawn`, and `--detach` keeps only the last call alive. Fixed by removing that line |
| `kras_smoke2110_set1` | Smoke, forced targeting on | Refused: avoid targets had no hotspots ("Epitope-focused design needs target hotspots"). Fixed by adding hotspots |
| `kras_smoke4259_set1` | Smoke, forced targeting on, 100 min | est $6.58, 11 trajectories, 0 accepted (G12D ipTM collapse at harden) |
| `kras_smoke3220_set1` | Smoke, forced targeting off, 50 min | est $3.29, 5 trajectories, 0 accepted |
| `kras_set1` | 4 shards, only final-stage floors loosened (0.5/0.55) | est $13.17, 18 trajectories, 0 accepted. Funnel: 17 passed screen, 16 refine, 4 anneal, 1 harden, 0 accepted |
| `kras_set2` | G12D-only config, hotspot set 2, 4 shards | est $16.28, 32 designs refolded, 0 passing. G12D ipSAE about 0.94-0.98 but WT 0.94-0.98 and peptide-free 0.7-0.95: pan-HLA binders |
| `kras_set4` | Peptide-free-only avoid config | Stopped after about 10 minutes, about $1-2 |
| `kras_set5` | **Bind G12D, avoid both, all floors 0.3, peptide-only hotspots, 2 shards** | est $7.27, 6 designs refolded, 0 passing the gates. Best by G12D ipSAE: 99 aa, G12D 0.72 but WT 0.43 and peptide-free 0.95 (not specific). Best specificity: 78 aa, G12D ipSAE 0.512, WT 0.43, peptide-free 0.0, delta 0.08, Asp6 distance 6.92 A |
| `kras_set5` top-up | Cap $60 total, up to 14 shards | done: 14 shards, 74 designs refolded (includes the first 6), est $53.57, 0 passing the tier A/B gates. Best by gap: 66 aa, G12D 0.87, WT 0.44, peptide-free 0.25 |
| `kras_rescore1` | Top 20 pooled designs refolded alone on all 3 targets, 8 samples each | 60 folds, 8.8 min, est compute $2.72. Top gaps held up (see RESCORE_RESULTS.md) |
| `kras_rescore2` | The other 45 designs with G12D ipSAE >= 0.5, same method | 135 folds; finished, changed the ranking little: the final top 20 is mostly the first 20 |

Bug found on the way: whole-complex ipTM was always about 0.8 because of the HLA-b2m pair; fixed with `binder_iptm` (pair_chains_iptm).

Cumulative BindCraft2 spend before the top-up was about $48 by the controller estimates; the top-up added est $53.57 for set5 as a whole. Rescoring added an estimated ~$9. None of this is confirmed against the Modal dashboard.

Selection: pooled set2 + set5 (106 designs, 65 unique with G12D ipSAE >= 0.5), ranked by G12D ipSAE minus the higher avoid score, rescored with 8 samples; details in RESCORE_RESULTS.md.

Conclusion so far: bind-G12D-only design yields pan-HLA binders; the three-target design at default floors yields almost nothing; at floors 0.3 it yields designs, but none has yet shown G12D specificity above the gates.


## Selection and rescoring (2026-11-01)

- Pooled set2 (32) and set5 (74) refolds: 106 designs, 65 unique sequences with G12D ipSAE >= 0.5. Ranked by delta = G12D ipSAE (mean of 3 samples) minus the higher of WT and peptide-free (max of 3). Identity filters of 80% and 95% were tried; the unfiltered top 20 had the same sequences as the 95% version.
- `kras_rescore1`: top 20 refolded with 8 samples per target. 60 folds, 8.8 min, compute est $2.72.
- `kras_rescore2`: the other 45 designs, 8 samples per target, 135 folds. Cost not recorded here; check the Modal dashboard.
- Final 20 chosen by score = (delta_mean + delta_worst) / 2 over all 65. The second rescore did not add a better design than the top 5 from the first.

Final top 20 metrics (8 samples; sequences not published here). ipSAE values are means over samples; delta_worst uses the worst avoid sample.

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

Reading it: only ranks 1 and 2 show a clear gap that holds on the worst sample. Rank 3 is 13 A from Asp6. Ranks 15 and 16 discriminate WT by peptide-chain ipSAE but bind peptide-free HLA at 0.9+, so they are not peptide-dependent. Ranks 14-20 are close to pan-HLA binders. All numbers are predictions; none has been tested in the lab.
