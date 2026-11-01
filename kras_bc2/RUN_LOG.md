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
