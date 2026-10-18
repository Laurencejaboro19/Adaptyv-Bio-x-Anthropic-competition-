# Challenge 3 (GLP1R): BoltzProt-1 designs, Boltz-2.1 refold, contact check, ipSAE

This section is the **BoltzProt-1 side test** for Challenge 3. It is separate from the BindCraft2 campaigns in the rest of this repo.

- **What was run:** 50 "bridge" designs (30-50 aa) and 50 "deep" designs (35-60 aa) from BoltzProt-1, all de novo (no GLP-1 sequence, no seeding, no motif scaffolding, no peptides or VHHs). Every design was refolded with Boltz-2.1 as a receptor + binder complex.
- **What is in this folder:** only the PowerShell that was run on a Windows PC and worked, in the order it was run. Nothing here is untested.
- **What is not in this folder:** design sequences and refold metrics. Challenge 3 closes on 25 Oct 2026 and this repo is public.

Everything runs in **Windows PowerShell opened from the app** (not from a project folder). All paths are absolute and sit under `Downloads`, so you can paste each block from anywhere. Each script is one paste: it creates its own folders and files, so you create nothing by hand.

## What is NOT included, and why

| Step | Why it is not here |
|---|---|
| Submitting the 100 BoltzProt-1 design jobs | Run through the Boltz API from a Claude session, not from a PowerShell script you ran. No tested PowerShell exists for it. |
| Submitting the 100 Boltz-2.1 refold jobs | Same: run through the Boltz API from a Claude session. |
| Merging refold metrics into the 100-design CSVs | Done in a Claude sandbox from the API results, not on your PC. |
| ipSAE on all 100 designs | Never run on the PC. Only the 4 contact-passing designs were scored (Block 3). |
| First-version download scripts and the first key prompt | The first key prompt captured one character and all 100 downloads failed. Block 2 uses the corrected prompt that worked. |

The settings used for the two submitted steps are recorded below so the run can be repeated.

### Settings used (documentation, not code)

- **Design:** BoltzProt-1 protein design against a C-alpha-only crop of 6X18 chain R (203 residues, the active GLP-1 bound GLP1R-Gs structure). 50 designs at 30-50 aa ("bridge") and 50 at 35-60 aa ("deep"). The exact epitope residues sent to the design jobs could not be recovered afterwards. The contact check in Block 2 therefore uses the BindCraft2 hotspot and coldspot lists for the same two campaigns as its working definition.
- **Refold:** Boltz-2.1 structure-and-binding job (`protein_protein_binding`). Chain A is the receptor, 409 aa, UniProt P43220 residues 24-432. Chain B is the binder. Empty MSA, 1 sample.
- **Numbering:** UniProt residue = position in the Boltz receptor sequence + 23. Block 2 checks this on every file by requiring residue 29 to be THR.

## Run order

Run the blocks in order, in one PowerShell window.

### Block 1: get the designed sequences (`01_get_designed_sequences.ps1`)

**What it does:** reads `glp1r_v2_links_scores.csv` from your Downloads folder, downloads each of the 100 designed structures, reads the binder sequence out of each one, and writes three CSVs.

**Needs first:** `glp1r_v2_links_scores.csv` in `Downloads`. It lists the 100 design results with these columns: `design, campaign, length_expected, iptm, min_interaction_pae, binding_confidence, structure_confidence, helix_fraction, sheet_fraction, loop_fraction, result_id, url`. The `url` column holds signed download links that **expire about 30 minutes after they are created**, so run the script straight away.

**How long it takes:** not timed. It is 100 small downloads one after another, so expect a few minutes.

**What you should see:** the last two lines print the counts (total, bridge, deep, length mismatches, download failures) and `Wrote glp1r_100_designs_v2.csv, bridge_30-50_designs_v2.csv, deep_35-60_designs_v2.csv to Downloads`. Good looks like 100 total, 50 bridge, 50 deep, 0 length mismatches and 0 download failures. It also creates `Downloads\bp2_cifs\` with the 100 CIFs.

**Why it works:** the binder is the only chain of 60 residues or fewer in each file (the receptor is far longer). The script reads only the C-alpha line of each residue (`label_atom_id = CA`), turns the three-letter code into one letter, and builds one string per chain. `length_expected` comes from the design request, so a mismatch with the length it reads back would flag a bad parse.

### Block 2: download the refolds and run the contact check (`02_download_refolds_and_contact_check.ps1`)

**What it does:** downloads the Boltz-2.1 refold of all 100 designs (the structure CIF and the full prediction archive), lists what is inside one archive, then checks which binder residues touch the receptor hotspots and coldspots.

**Needs first:** your Boltz API key. The script asks for it in a hidden prompt (right-click to paste, then Enter). **The key is never written into the script or a file**, and it is cleared from the session at the end. If you pasted a key anywhere shared, revoke it and make a new one.

**How long it takes:** not timed. It is 200 small downloads one after another (a CIF and an archive for each design), so expect a few minutes. The contact check itself takes seconds.

**What you should see, in order:**

1. `Key characters received: <length of your key>`. If it is under 10 the script stops and calls nothing.
2. `Test call status: succeeded`, which proves the key and the endpoint work.
3. `ok bridge_01` ... `ok deep_50`, one line per design.
4. `cifs downloaded: 100 / 100; failed: 0`.
5. `Files inside one archive:` followed by `\prediction\metrics.json`, `\prediction\sample_0_pae.npz` and `\prediction\sample_0_predicted_structure.cif`.
6. `wrote C:\Users\<you>\Downloads\glp1r_refold\contact_check.csv 100 designs; 4 pass`, then `bridge 50 pass 3` and `deep 50 pass 1`.

**Why it works:**

- The signed download links in the API response last about 30 minutes, so the script asks the API for fresh links for each design immediately before downloading.
- The archive holds the PAE matrix (`sample_0_pae.npz`), which Block 3 needs.
- The contact check counts a receptor and a binder residue as in contact when any heavy atoms are within 4.5 Angstrom. A design passes when it:
  - hits at least 30% of the hotspots for its campaign,
  - touches at most 2 coldspots,
  - has at least 10 binder residues in contact, and
  - keeps its last three residues (the C-terminus) out of the interface.
- The hotspot and coldspot lists are in the script (`BRIDGE_HOT`, `DEEP_HOT`, `COLD`). The list used is chosen from the file name (`bridge_*` or `deep_*`).

**Result from this run:** 4 of 100 designs passed (bridge 3, deep 1). The thresholds were chosen by us and are fairly strict, so a pass is a useful filter, not a verdict.

### Block 3: ipSAE on the designs that passed (`03_ipsae_on_contact_passes.ps1`)

**What it does:** reads the PASS rows from `contact_check.csv`, unpacks only those archives, computes ipSAE for each one, and writes `ipsae_pass4.csv`.

**How long it takes:** the run for 4 designs finished in about a couple of minutes, including `pip install`. It is CPU only, with no GPU and no cloud.

**What you should see:** `PASS designs: ...` with the names, then one line per design, then `wrote C:\Users\<you>\Downloads\glp1r_refold\ipsae_pass4.csv`. Each line has `ipsae_A_to_B`, `ipsae_B_to_A`, `ipsae_max` and the contact columns, sorted by `ipsae_max`, highest first.

**Why it works:** it uses the `ipsae` Python package with a PAE cutoff of 15 and a distance cutoff of 10 (the setting for a large target), and reports the larger of the two directions (receptor to binder, binder to receptor). Two details matter:

- The package's own file-saving step crashes on these files, so the script switches it off. The scores are unaffected.
- The package expects the PAE array under the key `pae`; the script copies it to a temporary file under that name.

**Result from this run (4 contact-passing designs):**

| Design | ipSAE | Hotspots hit | Coldspots hit |
|---|---|---|---|
| bridge_45 | 0.88 | 13 of 18 | 1 |
| deep_15 | 0.74 | 3 of 9 | 1 |
| bridge_25 | 0.67 | 12 of 18 | 1 |
| bridge_07 | 0.53 | 12 of 18 | 0 |

## Caveats

- **All evidence is computational.** A high refold ipTM or ipSAE says the binder is predicted to bind the receptor. It does not say it activates the receptor (cAMP signalling) in the assay.
- **The bridge designs reach the lower rim of the receptor only.** The deep designs reach the pocket without a groove anchor. ECD-only GLP1R binders reported by Baker's lab were antagonists.
- **Boltz-2.1 is the hosted API model.** The BindCraft2 campaigns elsewhere in this repo use the open-source Boltz-2 (2.2.1), so numbers from the two will not match exactly.
- **Prediction archives expire about 25 Oct 2026.** Download them before then.
