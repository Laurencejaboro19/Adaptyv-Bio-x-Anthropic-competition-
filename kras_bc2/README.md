# Challenge 5: KRAS G12D / HLA-A*11:01 (VVVGADGVGK), BindCraft2 -> Boltz-2 -> ipSAE on Modal

Goal: a binder that binds the G12D peptide-HLA complex and avoids the wild-type complex (VVVGAGGVGK) and the peptide-free HLA + b2m.
Run from Windows PowerShell. Everything runs on [Modal](https://modal.com). **All evidence is computational. No design has passed the strict gates; see RUN_LOG.md and RESCORE_RESULTS.md.**
No design sequences, API keys or run outputs are in this folder.

## What is here

| File | What it does | Status |
|---|---|---|
| `kras_app.py` | Modal app: BindCraft2 design (H100), Boltz-2 refold of each design on G12D, WT and peptide-free (L40S, 3 samples), ipSAE per target chain, contact check, spend-capped controller | final version incl. the rescore functions (fold_n, rescore_ctl, rescore_start, rescore_get); see note below |
| `00_modal_login.ps1` | Install Modal client, log in | one-time |
| `01_setup_prep_warm_cache.ps1` | Upload the 3 PDBs, cache Boltz-2 weights and HLA/b2m MSAs | ran, worked |
| `02_smoke_test.ps1` | 1 shard x 30 trajectories, refold on all 3 | ran several times, see RUN_LOG |
| `03_refold_test.ps1` | Refold a dummy 78-aa sequence on all 3 targets to test the scoring path | ran |
| `04_status_logs.ps1` | Progress, shard log tails, file listing | ran |
| `05_launch_set5_about_7usd.ps1` | **The run that worked for about $7**: bind G12D, avoid both, all floors 0.3 | ran, worked as a pipeline |
| `06_topup_cap60.ps1` | Resume the same tag with a $60 total cap | ran: 14 shards, 74 designs, est $53.57, 0 passed tier A/B |
| `07_fetch_results.ps1` | Download `designs.csv`, show top rows | ran |
| `08_stop_app.ps1` | Stop a Modal app | ran once (kras_set4) |
| `09_select_top20.ps1` | Pool set2 + set5, rank by G12D-minus-worst-avoid ipSAE gap, save top 20 | ran |
| `10_rescore_top20_detached.ps1` | Detached rescore of the top 20 on all 3 targets, 8 samples, adds peptide-chain ipSAE and Asp6 distance | ran, 8.8 min, est $2.72 |
| `11_rescore_get.ps1` | Check and download a rescore run | ran |
| `12_rescore_remaining45_start.ps1` | Same rescore for the other 45 designs with G12D ipSAE >= 0.5 | ran (launched 10:50 UK) |
| `13_get_and_merge_final_top20.ps1` | Merge both rescores, rank all 65, save the final top 20 | ran |
| `RESCORE_RESULTS.md` | Selection method and final top-20 metrics (no sequences) | |
| `configs/campaign_set5_shard0.json` | The exact BindCraft2 campaign JSON for set5 (campaign() merged with LOOSE), shard 0 | generated from `kras_app.py` |
| `inputs/` | The 3 target PDBs | see below |

Note on `kras_app.py`: the file you run on the PC was written by a PowerShell here-string and then patched in place several times. This copy was rebuilt from the original source plus those same patches in order, then compiled and its set5 config output checked. It should match what ran, but it is a reconstruction, not a byte copy of the PC file (comments differ).

## Inputs

`inputs/` holds the three target files used, named as the app expects: `target_G12D_ABC.pdb` (G12D complex, chains A HLA 1-276, B b2m 1-99, C peptide 1-10), `avoid_wildtype_ABC.pdb` (WT complex, A 2-276), `avoid_peptide_free_AB.pdb` (HLA + b2m only). They are the 7OW6 / 7OW5 complexes with the TCR removed; `prep` checks chains and peptide sequences against the expected values. Peptide residue 6 (chain C) is the Asp6 / Gly6 difference.

## Run order

1. `00` once, then `01`.
2. Optional cheap checks: `03`, then `02`.
3. `05` for the $7 run, `04` to watch it, `07` to read results.
4. `06` for the top-up, then `07` again.
5. `09` to pick the top 20, `10` then `11` to rescore them, `12` then `13` to rescore the rest and choose the final 20.

## The settings that mattered (all in `kras_app.py`)

- Three targets in one BindCraft2 campaign: G12D `target` (weight 1.0, hotspots, coldspots A183-276 and B1-99), WT `detarget` (-0.3), peptide-free `detarget` (-1.0).
- Set5 hotspots: peptide residues C4,C5,C6,C7,C8 (Asp6 is C6). The WT target gets the same hotspots; the peptide-free target gets only the chain-A hotspots, which for set5 is an empty string.
- `forced_targeting` False (`FORCED`). With it on, BindCraft2 restored the wild-type surface at the harden stage and G12D ipTM collapsed.
- `LOOSE` dict: every stage floor (plddt, iptm, ptm, monomer plddt) 0.3, `max_ipae_final` 0.7, `min_hotspot_contact_final` 0.3, `max_detarget_interface_residues_final` 30, detarget ipTM ceilings 0.75. Lowering only the final-stage floors did nothing; the per-stage floors (anneal, harden, mutate) were what rejected designs.
- Boltz-2 refold: `boltz==2.2.1`, recycling 3, sampling 100, 3 diffusion samples, binder and peptide `msa: empty`, HLA and b2m MSAs cached. Binder-vs-target ipTM is taken from `pair_chains_iptm` because whole-complex ipTM is dominated by the HLA-b2m pair.
- ipSAE (`ipsae==1.0.1`, cutoff 15, dist 10) per target chain; G12D uses the mean of its samples, avoid targets the max.
- Gates: G12D ipTM >= 0.65, ipSAE >= 0.55, contact checks incl. Asp6 within 4.5 A; tier A needs WT <= 0.30 and peptide-free <= 0.20, tier B needs delta >= 0.15.
- Spend: the controller's cap is its own estimate (H100 $3.95/h, L40S $1.95/h). The Modal dashboard is the real number.

Throughput seen: about 10 minutes per trajectory with three targets, about $0.66 per trajectory.

## Caveats

- BindCraft2 commit pinned to 301efdd1937fb963cc40a5b0ecc1bc2f9b2b2d15.
- A random non-binder once scored ipSAE 0.375 on WT in 1 of 3 samples, so ipSAE has noise at this level.
- Rotate any API key that was ever pasted in a chat.
