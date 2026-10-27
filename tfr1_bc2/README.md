# TfR1 (Challenge 4): BindCraft2 -> Boltz-2 -> ipSAE -> PROPKA on Modal (work in progress)

PowerShell blocks that were run on a Windows PC and worked, in order. `tfr1_app.py` (the Modal app) is being added; until it is here these blocks cannot be run from this folder alone.

- **Not included while the challenge is open:** design sequences, hotspot and coldspot residue numbers, per-design results. Block 2 has a placeholder for the spots file.
- **Run totals (the app's own estimates, not Modal billing):** 6 shards, 3 rounds, 151 sequences refolded and scored, 37 passed every gate, 3.67 h, about $28.77.
- **Pipeline:** BindCraft2 dual-target (human and mouse) design on an H100, Boltz-2 refold of each sequence on both species, ipSAE and contact scoring, then PROPKA with a Wyman-linkage pH 6.0 to 7.4 estimate.
- **Lessons:** with BindCraft2's default pLDDT floors (0.6) two 20-trajectory smoke tests accepted 0 of 20, so the floors were lowered to 0.4 (screen) and 0.5 (refine) in Block 2. A first hotspot face accepted 0 of 40 trajectories; a second face worked.
- **Note:** `fetch` (Block 6) downloads only the final CSVs and FASTA, not the BindCraft2 trajectory logs.

Run in Windows PowerShell opened from the app, in the order of the file numbers.
