# Challenge 3 (GLP1R): BindCraft2 -> Boltz-2 refold -> ipSAE on Modal (bridge, 30-50 aa)

This section is the **BindCraft2 (BC2) campaign** for Challenge 3, run on Modal so it keeps going with the PC switched off. It is separate from the BoltzProt-1 side test in [`../glp1r_boltzprot1/`](../glp1r_boltzprot1/).

- **What was run:** the **bridge** campaign (30-50 aa). BC2 designed binders on an H100; each sequence was refolded with open-source Boltz-2 (2.2.1) as a receptor + binder complex, then scored with ipSAE and a contact check. All designs are de novo (no GLP-1 sequence, no seeding, no motif scaffolding, no peptides or VHHs).
- **What is in this folder:** `glp1r_app.py` (the Modal app) and the PowerShell blocks that were run on a Windows PC and worked, in the order they were run.
- **What is not in this folder:** design sequences and per-design results. Challenge 3 closes on 25 Oct 2026 and this repo is public.
- **Not validated:** the **deep** campaign (35-60 aa) is in `glp1r_app.py` but was never run, so there is no deep block here.

Everything runs in **Windows PowerShell opened from the app** (not from a project folder). All paths are absolute and sit under `Downloads`, so you can paste each block from anywhere. Block 1 writes the app file itself, so you create nothing by hand. Modal must be installed and you must be able to log in (Block 0).

## Settings used (documentation, not code)

- **Target:** 6X18 chain R (the active GLP-1 bound GLP1R-Gs structure), cropped to 203 residues and renamed chain A. Residue numbers are UniProt P43220 numbers. The crop is built inside the Modal image from the RCSB file, and the build stops if the crop is not exactly 203 residues.
- **Bridge hotspots (18):** A32, A39, A69, A88, A91, A121 (extracellular domain groove) and A197, A201, A205, A210, A211, A214, A215, A221, A298, A299, A300, A306 (extracellular rim: top of TM2, ECL1, ECL2, top of TM5).
- **Coldspots (22):** A63, A82, A115, A146, A207, A208, A279, A282, A283, A286, A290, A303, A304, A307, A311, A371, A373, A374, A376, A378, A381, A382.
- **BindCraft2:** commit `301efdd1937fb963cc40a5b0ecc1bc2f9b2b2d15`, H100, 3 workers per GPU, binder 30-50 aa, `aa_bias` Cys 0, `min_hotspot_contact_final` 0.55, `max_coldspot_contact_final` 0.1, `min_iptm_final` 0.55, `min_monomer_plddt_final` 0.65, Interface_Residues >= 10, and a final `C_Terminus_Away_Cosine` filter at 0.0. Each shard runs 25 trajectories and is stopped after 100 minutes.
- **Refold:** `boltz==2.2.1` on an L40S, `--model boltz2 --no_kernels --write_full_pae --diffusion_samples 1 --recycling_steps 3 --sampling_steps 100`. The receptor MSA comes from the ColabFold server and is computed once and cached; the binder MSA is empty.
- **Scoring:** `ipsae==1.0.1`, PAE cutoff 15, distance cutoff 10, reported as the larger of the two directions. The contact check counts a pair as in contact at 4.5 Angstrom.
- **"Good" means all of:** ipTM >= 0.60, ipSAE >= 0.50, hotspot fraction >= 0.30, at most 2 coldspots hit, at least 10 binder residues in contact, and the last three binder residues not in the interface.
- **Cost control:** the controller runs inside Modal with its own spend cap. The dollar figures it prints are **its own estimates** from GPU time, not Modal billing, so check the Modal billing page.

## Run order

Run the blocks in order, in one PowerShell window.

### Block 0: log in to Modal (`00_modal_login.ps1`)

**What it does:** opens a browser page to approve the login and saves a token on the PC.
**How long it takes:** about 30 seconds, once per PC.
**What you should see:** `Token verified successfully!` and `Token written to ... .modal.toml`.
**Why it works:** the token is stored in your user profile, so later blocks need no login.

### Block 1: write the app and cache the receptor MSA (`01_write_app_and_warm_cache.ps1`)

**What it does:** creates `Downloads\glp1r_modal`, writes `glp1r_app.py` into it (the same file as in this folder), and runs `warm_cache`. That builds the Modal images and saves the receptor MSA on a Modal volume.
**How long it takes:** the first run builds the cloud images and fetches the MSA (not timed). Reruns finish in seconds because the images and MSA are cached.
**What you should see:** the list of created Modal functions, then `warm OK`. On a rerun it says `warm OK (already done)`. In the image build log the crop step prints `crop residues 203`.
**Why it works:** nothing is downloaded to your PC; Modal builds the images in the cloud and reuses cached layers. The weights and compile cache live on the Modal volume `pdfa-bc2-weights`, the MSA on `pdfa-boltz-cache`, and the results on `pdfa-bc2-results`.

### Block 2: smoke test (`02_smoke_test.ps1`)

**What it does:** runs one H100 shard of 20 BC2 trajectories, harvests the sequences, refolds them with Boltz-2 and scores them. It checks the whole chain before the real spend.
**How long it takes:** about 20-25 minutes (the run recorded 0.36 hours). **Keep the PC on until the result prints**; this block was run with the PC on and the call waits for the result.
**What you should see:** a final line like `STOPPING: shard ceiling reached | shards 1 | refolded N | good M | est $X`, then a progress block.
**What good looks like:** `shards_launched` 1, `sequences_seen` above 0 and `refolded` above 0. `good` can be 0 on 20 trajectories.
**Result from this run:** 20 trajectories, 6 sequences refolded, 0 good, about $1.59 by the app's estimate. This smoke run used the same settings as the campaign except the final `C_Terminus_Away_Cosine` filter, which was added afterwards and was tested only in the campaign itself.
**Note:** running it again with the same tag resumes the saved state and does nothing. Change the tag text in the file to run it fresh.
**Why it works:** the controller in `glp1r_app.py` runs the design shard, then `harvest`, then `refold`, then `score_one`, and saves progress after each stage.

### Block 3: launch the campaign in the cloud (`03_launch_bridge_campaign.ps1`)

**What it does:** starts the controller inside Modal. It runs rounds of BC2 design shards, harvests the sequences, refolds them with Boltz-2, and scores contacts and ipSAE. It stops at 14 good designs, 12 shards, the cap (`--budget 45`) or 20 hours.
**How long it takes:** the launch takes a minute. The campaign itself ran for **1.67 hours** here.
**What you should see:** `CONTROLLER STARTED IN THE CLOUD  tag: glp1r_bridge  goal: 14  budget cap (est): $ 45.0`, then `=== round 1: shards [0, 1, 2, 3] ...`. A `Timed out waiting for final app logs` line afterwards is harmless.
**When to switch the PC off:** as soon as the `CONTROLLER STARTED` line prints.
**Why it works:** the controller is launched with `spawn`, so it keeps running in Modal after the PC disconnects, and it starts the design, refold and scoring jobs from there.
**Result from this run:** 12 shards in 4 rounds, 132 sequences refolded, **4 passed every gate**, about $23.94 by the app's estimate. It stopped at the 12-shard limit before reaching 14 good designs.

### Block 3b (optional): resume for more designs (`03b_resume_for_more_designs.ps1`)

**What it does:** restarts the same campaign with a higher shard limit. The saved state is reused: earlier sequences are not refolded again and new shards get new seeds. `--budget` is a **total** cap that includes what was already spent.
**What you should see:** `=== round 1: shards [12, 13, ...] | good 4/14 | est $23.94`, which shows the earlier results were picked up.
**Status:** the launch and resume were confirmed. Its results are not reported here.

### Block 4: progress (`04_status.ps1`)

**What it does:** prints the saved progress for the smoke test, the bridge campaign and the deep campaign (deep says `no progress file yet`).
**How long it takes:** about 30 seconds. It works while the campaign runs and the PC is off.
**What good looks like:** `state` moving through `designing`, `harvesting`, `refolding`, `scoring`, `done`; `sequences_seen` and `refolded` rising; the estimated spend growing slowly.

### Block 5: BindCraft2's log and files (`05_logs_and_files.ps1`)

**What it does:** prints the last lines of BC2's own log for each shard and lists the files BC2 wrote for a run.
**How long it takes:** under a minute.
**What to look for:** trajectory headers of the form `accepted N/999`, and `rejected at ... design stage` lines with i_pTM and pLDDT. This is the quickest way to see whether BC2 is accepting designs.

### Block 6: download the CSV and pull out the passing designs (`06_fetch_csv_and_good_designs.ps1`)

**What it does:** downloads `glp1r_bc2_refolded_ipsae.csv` (one row per refolded sequence, good designs first), prints the top 15 rows, and saves the passing designs to `glp1r_bridge_good_designs.csv` in `Downloads\glp1r_bc2_results`.
**How long it takes:** about a minute.
**What you should see:** `TOTAL ROWS: 132`, `GOOD: 4`, a table, and `saved 4 passing designs to ...`.
**What good looks like:** ipTM above 0.6, ipSAE above 0.5, `contact_PASS` True and `cterm_in_interface` False on the good rows.

**Result from this run (4 passing designs, sequences not published here):**

| Length | ipTM | ipSAE | Hotspots hit | Coldspots hit |
|---|---|---|---|---|
| 47 | 0.94 | 0.95 | 13 of 18 | 2 |
| 38 | 0.93 | 0.93 | 15 of 18 | 2 |
| 48 | 0.87 | 0.83 | 7 of 18 | 2 |
| 50 | 0.77 | 0.64 | 13 of 18 | 2 |

## Caveats

- **All evidence is computational.** A high refold ipTM or ipSAE says the binder is predicted to bind the receptor. It does not say it activates the receptor (cAMP signalling) in the assay.
- **BC2's own check accepted none of the first 100 trajectories** (0 of 25 in each of the first four shards). The harvest step takes accepted sequences when BC2 writes them and otherwise falls back to every sequence of 25-80 aa in BC2's output files, so these designs are judged mainly by the Boltz-2 refold, not by BC2's own validation.
- **Single-sample Boltz-2 scores are noisy.** Differences of about 0.05-0.1 should not be over-read, and all four passing designs sit at the coldspot limit of 2.
- **The refold setup differs from the BoltzProt-1 section.** This section uses open-source Boltz-2 2.2.1 with a cached receptor MSA; the BoltzProt-1 section uses the hosted Boltz-2.1 API with an empty MSA, so numbers from the two are not strictly comparable.
- **"Good" counts designs, not independent backbones.** Up to two sequences per backbone can both pass.
- **Cost figures are the app's estimates, not Modal billing.**
