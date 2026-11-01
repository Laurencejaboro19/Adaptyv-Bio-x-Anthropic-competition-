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
| `kras_set5` top-up | Cap $60 total, up to 14 shards | launched 2026-11-01 00:04 UK, results not included here |

Bug found on the way: whole-complex ipTM was always about 0.8 because of the HLA-b2m pair; fixed with `binder_iptm` (pair_chains_iptm).

Cumulative BindCraft2 spend before the top-up was about $48 by the controller estimates (not confirmed against the Modal dashboard).

Conclusion so far: bind-G12D-only design yields pan-HLA binders; the three-target design at default floors yields almost nothing; at floors 0.3 it yields designs, but none has yet shown G12D specificity above the gates.
