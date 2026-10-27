# TfR1 (Challenge 4): BindCraft2 -> Boltz-2 -> ipSAE -> PROPKA on Modal

Seven PowerShell blocks, in the order they were run on a Windows PC, each one that worked. Run them in Windows PowerShell opened from the app; paths are absolute under `Downloads`.

**Not here yet:** `tfr1_app.py` (the Modal app) and the `tfr1_inputs` folder, which Block 1 unzips from `TfR1_pdfa.zip`. The app file is being added. **Not published while the challenge is open:** design sequences, hotspot and coldspot residue numbers, per-design results. Block 4 has a placeholder for the spots file.

## The blocks

| # | File | What it does | Outcome |
|---|---|---|---|
| 1 | `01_unzip_install_login.ps1` | Unzips the package, sets binder length 60-120 and budget $70, installs Modal, logs in | Worked |
| 2 | `02_gemmi_patch_and_prep.ps1` | Adds `gemmi` to the cloud image, uploads and validates the human and mouse targets | Worked (`PREP OK`). The first `prep` without this patch failed |
| 3 | `03_check.ps1` | Plumbing test: refolds and scores a throwaway control on both species | Worked (`CHECK DONE`) |
| 4 | `04_patch_spots_prep_launch.ps1` | Sets binder length 60-100, 24 shards, pLDDT floors 0.4/0.5; writes the spots file; re-uploads targets; starts the cloud run (tag `tfr1_pB`, $70 cap, goal 24) | Worked |
| 5 | `05_status.ps1` | Shows run progress | Worked |
| 6 | `06_fetch_results.ps1` | Downloads the result CSVs and FASTA | Worked |
| 7 | `07_top20_one_csv.ps1` | Ranks the passing designs, saves the top 20 as one CSV | Worked |

Left out: two smoke tests (20 trajectories each) that accepted 0 of 20 designs, status and log checks on them, and a first-row preview of the results file whose output was never seen.

## Run numbers (the app's own estimates, not Modal billing)

6 shards in 3 rounds, 151 sequences refolded and scored, 37 passed every gate, 3.67 h, about $28.77.

## Pipeline and settings

BindCraft2 dual-target (human and mouse apical domain, 60-100 aa) on an H100; each sequence refolded with Boltz-2 on both species; ipSAE and contact scoring; PROPKA with a Wyman-linkage pH 6.0 to 7.4 estimate. Gates, both species: ipTM 0.75, ipSAE 0.60, contact PAE 7, binder pLDDT 0.80, hotspot fraction 0.40, coldspot fraction 0.30, His fold 0.30 log10.

## Lessons

- BindCraft2's default pLDDT floors (0.6) rejected every trajectory; Block 4 lowers them to 0.4 (screen) and 0.5 (refine).
- A first hotspot face accepted 0 of 40 trajectories; a second face worked.
- `fetch` downloads only the final CSVs and FASTA, not the BindCraft2 trajectory logs.
- Rank score in Block 7 is the lower ipSAE across species plus 0.2 times the lower His fold (capped at 1.5).

## spots.json format

`{"human": {"hot": "A<start>-<end>,...", "cold": "...", "hot_list": [ints], "cold_list": [ints]}, "mouse": {...same keys...}}` with UniProt residue numbers on chain A.
