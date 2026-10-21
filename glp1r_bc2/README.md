# Challenge 3 (GLP1R): BindCraft2 -> Boltz-2 refold -> ipSAE on Modal (bridge 30-50 aa and deep 35-60 aa)

This section is the **BindCraft2 (BC2) campaigns** for Challenge 3, run on Modal so it keeps going with the PC switched off. It is separate from the BoltzProt-1 side test in [`../glp1r_boltzprot1/`](../glp1r_boltzprot1/).

- **What was run:** the **bridge** campaign (30-50 aa, extracellular groove and rim) and the **deep** campaign (35-60 aa, upper transmembrane pocket). BC2 designed binders on an H100; each sequence was refolded with open-source Boltz-2 (2.2.1) as a receptor + binder complex, then scored with ipSAE and a contact check. All designs are de novo (no GLP-1 sequence, no seeding, no motif scaffolding, no peptides or VHHs).
- **What is in this folder:** `glp1r_app.py` (the Modal app) and the PowerShell blocks that were run on a Windows PC and worked, in the order they were run.
- **What is not in this folder:** design sequences and per-design results. Challenge 3 closes on 25 Oct 2026 and this repo is public.
- **Final numbers (app's own estimates):** bridge 27 shards, 523 sequences refolded, **14 passed every gate**, about $62.23, 2.41 h. Deep 6 shards, 195 refolded, **4 passed**, about $20.80, 2.14 h. The 4 deep passes are not 4 independent backbones: 3 of them share one backbone.

Everything runs in **Windows PowerShell opened from the app** (not from a project folder). All paths are absolute and sit under `Downloads`, so you can paste each block from anywhere. Block 1 writes the app file itself, so you create nothing by hand. Modal must be installed and you must be able to log in (Block 0).

## Settings used (documentation, not code)

- **Target:** 6X18 chain R (the active GLP-1 bound GLP1R-Gs structure), cropped to 203 residues and renamed chain A. Residue numbers are UniProt P43220 numbers. The crop is built inside the Modal image from the RCSB file, and the build stops if the crop is not exactly 203 residues.
- **Deep hotspots (9):** A298, A299, A300, A306, A230, A233, A234, A309, A310. Deep binders are 35-60 aa with `min_hotspot_contact_final` 0.70; everything else matches the bridge settings.
- **Bridge hotspots (18):** A32, A39, A69, A88, A91, A121 (extracellular domain groove) and A197, A201, A205, A210, A211, A214, A215, A221, A298, A299, A300, A306 (extracellular rim: top of TM2, ECL1, ECL2, top of TM5).
- **Coldspots (22):** A63, A82, A115, A146, A207, A208, A279, A282, A283, A286, A290, A303, A304, A307, A311, A371, A373, A374, A376, A378, A381, A382.
- **BindCraft2:** commit `301efdd1937fb963cc40a5b0ecc1bc2f9b2b2d15`, H100, 3 workers per GPU, binder 30-50 aa (bridge), `aa_bias` Cys 0, `min_hotspot_contact_final` 0.55 (bridge), `max_coldspot_contact_final` 0.1, `min_iptm_final` 0.55, `min_monomer_plddt_final` 0.65, Interface_Residues >= 10, and a final `C_Terminus_Away_Cosine` filter at 0.0. `sequence_candidates` 30, forced targeting off, termini-accessible off. A C-terminus design loss was tried and removed because BC2 then accepted nothing. Each shard runs 25 trajectories and is stopped after 100 minutes.
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

**What it does:** runs one H100 shard of 30 BC2 trajectories in the cloud, harvests the sequences, refolds them with Boltz-2 and scores them. It checks the whole chain before the real spend. Set `$camp` to `bridge` or `deep`.
**How long it takes:** about 20-25 minutes. It is launched with `spawn`, so the PC can be switched off once the function call prints.
**What you should see:** a function call id, then the progress block. Check it later with Block 4.
**What good looks like:** `shards_launched` 1, `sequences_seen` above 0 and `refolded` above 0. `good` can be 0 on 30 trajectories.
**Results:** an early bridge smoke (20 trajectories, before the C-terminus filter) refolded 6 sequences, 0 good, about $1.59. The deep smoke refolded 22 sequences and 1 passed. That passing sequence is the same one as a design in the deep campaign, and its scores differ between the two runs (ipSAE 0.886 in the smoke test, 0.851 in the campaign), which shows how noisy a single Boltz-2 sample is.
**Note:** running it again with the same tag resumes the saved state and does nothing. Change the tag text (`glp1r_smoke4_`) in the file to run it fresh.
**Why it works:** the controller in `glp1r_app.py` runs the design shard, then `harvest`, then `refold`, then `score_one`, and saves progress after each stage.

### Block 3: launch the campaign in the cloud (`03_launch_bridge_campaign.ps1`)

**What it does:** starts the controller inside Modal. It runs rounds of BC2 design shards, harvests the sequences, refolds them with Boltz-2, and scores contacts and ipSAE. It stops at 14 good designs, 12 shards, the cap (`--budget 45`) or 20 hours (the shard limit applies unless changed).
**How long it takes:** the launch takes a minute. The first launch ran for **1.67 hours** here.
**What you should see:** `CONTROLLER STARTED IN THE CLOUD  tag: glp1r_bridge  goal: 14  budget cap (est): $ 45.0`, then `=== round 1: shards [0, 1, 2, 3] ...`. A `Timed out waiting for final app logs` line afterwards is harmless.
**When to switch the PC off:** as soon as the `CONTROLLER STARTED` line prints.
**Why it works:** the controller is launched with `spawn`, so it keeps running in Modal after the PC disconnects, and it starts the design, refold and scoring jobs from there.
**Result from the first launch:** 12 shards in 4 rounds, 132 sequences refolded, **4 passed every gate**, about $23.94 by the app's estimate. It stopped at the 12-shard limit before reaching 14 good designs.

### Block 3b (optional): resume for more designs (`03b_resume_for_more_designs.ps1`)

**What it does:** restarts the same campaign with a higher shard limit. The saved state is reused: earlier sequences are not refolded again and new shards get new seeds. `--budget` is a **total** cap that includes what was already spent.
**What you should see:** `=== round 1: shards [12, 13, ...] | good 4/14 | est $23.94`, which shows the earlier results were picked up.
**Result:** the resume reached the goal: 27 shards in total, 523 sequences refolded, **14 passed**, about $62.23 by the app's estimate, 2.41 h of campaign time.

### Block 3c: launch the deep campaign (`03c_launch_deep_campaign.ps1`)

**What it does:** the same controller for the deep site (35-60 aa, 9 hotspots). It stops at 10 good designs, 12 shards, the `$25` cap or 20 hours.
**When to switch the PC off:** once `CONTROLLER STARTED IN THE CLOUD` prints.
**Result:** it stopped at the cap after 6 shards: 195 sequences refolded, **4 passed**, about $20.80, 2.14 h. BC2 accepted few deep designs, so the goal of 10 was not reached.

### Block 4: progress (`04_status.ps1`)

**What it does:** prints the saved progress for the smoke test, the bridge campaign and the deep campaign. `04b_progress_bridge_and_deep.ps1` prints a shorter summary for the two campaigns.
**How long it takes:** about 30 seconds. It works while the campaign runs and the PC is off.
**What good looks like:** `state` moving through `designing`, `harvesting`, `refolding`, `scoring`, `done`; `sequences_seen` and `refolded` rising; the estimated spend growing slowly.

### Block 5: BindCraft2's log and files (`05_logs_and_files.ps1`)

**What it does:** prints the last lines of BC2's own log for each shard and lists the files BC2 wrote for a run.
**How long it takes:** under a minute.
**What to look for:** trajectory headers of the form `accepted N/999`, and `rejected at ... design stage` lines with i_pTM and pLDDT. This is the quickest way to see whether BC2 is accepting designs.

### Block 6: download the CSV and pull out the passing designs (`06_fetch_csv_and_good_designs.ps1`)

**What it does:** downloads `glp1r_bc2_refolded_ipsae.csv` (one row per refolded sequence) for both campaigns, prints the passing designs, and saves them to `glp1r_good_designs_bridge_and_deep.csv` in `Downloads\glp1r_final_results`.
**How long it takes:** about a minute.
**What you should see:** `TOTAL ROWS: 718` (523 bridge, 195 deep), `GOOD (bridge + deep): 18`, a table, and the save path.
**What good looks like:** ipTM above 0.6, ipSAE above 0.5, `contact_PASS` True and `cterm_in_interface` False on the good rows.

### Block 7 (optional): one CSV of the passing designs (`07_good_designs_one_csv.ps1`)

**What it does:** puts the passing bridge and deep designs, plus the passing deep smoke-test design, in one CSV in `Downloads`, best ipSAE first within each campaign. **How long:** seconds. The smoke-test design duplicates a campaign deep design, so count it once.

**Results (passing designs only, sequences not published here):** bridge 14, deep 4 (distinct). The four best bridge designs from the first launch were lengths 47, 38, 48 and 50 aa with ipTM 0.77-0.94 and ipSAE 0.64-0.95; the later resume added ten more. All passing designs sit at or near the coldspot limit of 2.

## Caveats

- **All evidence is computational.** A high refold ipTM or ipSAE says the binder is predicted to bind the receptor. It does not say it activates the receptor (cAMP signalling) in the assay.
- **BC2's own check accepted none of the first 100 trajectories** (0 of 25 in each of the first four shards), and 0 of the first 30 deep trajectories. The harvest step takes accepted sequences when BC2 writes them and otherwise falls back to every sequence of 25-80 aa in BC2's output files, so these designs are judged mainly by the Boltz-2 refold, not by BC2's own validation.
- **Single-sample Boltz-2 scores are noisy.** Differences of about 0.05-0.1 should not be over-read, and the passing designs sit at or near the coldspot limit of 2.
- **The refold setup differs from the BoltzProt-1 section.** This section uses open-source Boltz-2 2.2.1 with a cached receptor MSA; the BoltzProt-1 section uses the hosted Boltz-2.1 API with an empty MSA, so numbers from the two are not strictly comparable.
- **"Good" counts designs, not independent backbones.** Up to two sequences per backbone can both pass.
- **Cost figures are the app's estimates, not Modal billing.**
- **Deep diversity is limited:** 3 of the 4 deep passes share one backbone.
- **Not done here:** no GLP-1 positive-control refold and no multi-sample re-refold of the final list, and no activation (cAMP) evidence.
