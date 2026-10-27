# Adaptyv Bio x Anthropic protein design competition: PDFA code and notes

Code, configs and notes from Protein Design for Africa (PDFA) for the Anthropic x Adaptyv competition on Proteinbase.

| Challenge | Target | Folder |
|---|---|---|
| 1 | EGFR, conditional (pH-switch) binder, human and mouse | [`egfr/`](egfr/) |
| 2 | TNF-alpha, conditional (pH-switch) binder, human and mouse | [`tnf/`](tnf/) |
| 3 | GLP1R agonist, BoltzProt-1 side test (PowerShell only) | [`glp1r_boltzprot1/`](glp1r_boltzprot1/) |
| 3 | GLP1R agonist, BindCraft2 -> Boltz-2 -> ipSAE on Modal (bridge 30-50 aa, deep 35-60 aa) | [`glp1r_bc2/`](glp1r_bc2/) |
| 4 | TfR1, RFD4-Proteina generator on Modal (install, smoke test, design, read outputs; no scoring yet) | [`tfr1_rfd4/`](tfr1_rfd4/) |
| 4 | TfR1, BindCraft2 -> Boltz-2 -> ipSAE -> PROPKA on Modal (PowerShell blocks; sequences and hotspots withheld until the challenge closes) | [`tfr1_bc2/`](tfr1_bc2/) |

The design approach for both is in [`DESIGN_STRATEGY.md`](DESIGN_STRATEGY.md).

**All evidence is computational.** No wet-lab validation has been done. The pH preference is a PROPKA3 / Wyman-linkage estimate on a predicted structure, used to rank designs. It is not a measured affinity ratio.

## Pipeline

BindCraft2 (design) -> Boltz-2 (refold, source of truth) -> PROPKA3 + Wyman linkage (pH switch) -> ranking. Everything ran on [Modal](https://modal.com) and was launched from Windows PowerShell.

## Layout

```
egfr/
  modal/      egfr_dual10.py (design + Boltz-2 refold + PROPKA controller), modal_bc2_dual.py (design shards only),
              and the two EGFR target PDBs they load (keep them in the same folder as the scripts)
  scoring/    score_ph_switch.py (PROPKA3 wrapper, Wyman linkage, His geometry), runs offline on any complex
  configs/    BindCraft2 campaign JSONs, controller config, measured funnel and cost records (reference)
  notes/      EGFR_methodology.md
tnf/
  modal/      tnf_app.py (Modal app), tnf_client.py (launcher), tnf_smoke.py, and the four TNF-alpha target PDBs
  notes/      methodology, configs and design/validation notes
glp1r_boltzprot1/
  01_get_designed_sequences.ps1, 02_download_refolds_and_contact_check.ps1, 03_ipsae_on_contact_passes.ps1,
  README.md with run order, timings, expected output and why each block works (Challenge 3, BoltzProt-1 only)
glp1r_bc2/
  glp1r_app.py (Modal app: BindCraft2 design, Boltz-2 refold, ipSAE + contact check, spend-capped controller),
  00_modal_login.ps1 ... 07_good_designs_one_csv.ps1, README.md with run order, timings, expected output (bridge and deep, no sequences)
tfr1_rfd4/
  rfd4_app.py (Modal app: RFD4-Proteina install, weights, smoke test, design, fetch), 01_write_app_and_setup.ps1 ... 06_sequences_and_chains.ps1, README.md (generator only; no sequences or scores)
tfr1_bc2/
  01_unzip_install_login.ps1 ... 07_top20_one_csv.ps1, README.md (what each block does and whether it worked; tfr1_app.py to follow)
```

## How to run (PowerShell, in this order)

Each run spends real GPU money on Modal. Do the cheap check first, then a small run, then the full run.

### Step 0: one-time setup

What it does: installs the Modal client and logs your PowerShell into your Modal account.
What you should see next: a browser window opens to approve the login, then PowerShell says the token was verified.

```powershell
pip install modal
python -m modal token new
```

### EGFR (dual-species BindCraft2 campaign)

Step 1: go to the folder. Everything below runs from here, because the script loads `egfr_d3ext_human.pdb` and `egfr_d3ext_mouse.pdb` from beside it.

```powershell
cd egfr\modal
```

Step 2: self-check (about $0.30). What it does: builds the cloud image and checks the software and GPU are visible. What you should see next: labelled sections (`===== name =====`) and a `jax devices:` line listing the GPU. If anything errors, stop here and do not launch.

```powershell
python -m modal run modal_bc2_dual.py::selfcheck
```

Step 3: launch the campaign (10 shards, up to 30 trajectories each, on H100). What it does: runs BindCraft2 against human and mouse EGFR domain III at once and writes results to the Modal Volume `pdfa-bc2-results`. `--detach` keeps it running if you close PowerShell. What you should see next: a Modal run link, then (if you stay connected) a summary with `tag`, `shards clean: n/10`, minutes per shard, GPU-hours and the command to pull results.

```powershell
python -m modal run --detach modal_bc2_dual.py::launch
```

Step 4: check progress any time. What you should see next: one line per shard with its status JSON, or `nothing under /results/egfr_dual yet` if no shard has finished.

```powershell
python -m modal run modal_bc2_dual.py::status
```

Step 5: download results when shards are done. What you should see next: an `egfr_dual` folder under `.\results` with one `shardNN_seedXXXX` folder per shard, each holding its `campaign.json`.

```powershell
python -m modal volume get pdfa-bc2-results egfr_dual .\results
```

Step 6 (alternative to Steps 3 to 5): the full controller. What it does: runs design shards, refolds every candidate with Boltz-2 on human and mouse, scores the pH switch with PROPKA3 and Wyman linkage, and writes a ranked table. It stops at the budget you give it, and the default budget inside the script is $150, so always pass `--budget`. What you should see next: `CONTROLLER STARTED IN THE CLOUD`, a call id and the tag. You can then close PowerShell.

```powershell
python -m modal run --detach egfr_dual10.py::start --tag dual10 --goal 10 --budget 40
python -m modal run egfr_dual10.py::status --tag dual10
python -m modal run egfr_dual10.py::fetch --tag dual10
```

`status` prints the progress JSON, or `no progress file yet` early on. `fetch` downloads the ranked designs table.

Step 7 (optional, offline, no GPU): score the pH switch on any predicted complex. What you should see next: one JSON line with the interface histidines, the PROPKA linkage values and the anchor-residue pKa in the complex versus free. Needs `pip install biopython numpy propka`. Run from the repo root.

```powershell
python egfr\scoring\score_ph_switch.py complex.cif --target A --binder B --target-offset 304 --anchor 409
```

### TNF-alpha

Step 1: go to the folder. The four TNF PDBs sit beside the scripts here.

```powershell
cd tnf\modal
```

Step 2: check the input files. What it does: counts residues per chain in each PDB. What you should see next: `OK` lines for each file, then `ALL FILES GOOD`. If you see `BAD` or `PROBLEM`, fix the files before deploying.

```powershell
python tnf_client.py check
```

Step 3: deploy the app to Modal. What it does: uploads the app `pdfa-tnf-dual`. What you should see next: Modal prints the deployed app and its functions, with no error.

```powershell
python -m modal deploy tnf_app.py
```

Step 4: smoke test (about 10 minutes, under $1). What it does: sets up Boltz-2 weights and target MSAs, then refolds and scores one known control binder on human and mouse TNF-alpha. The control is not part of any submission. What you should see next: `1/2 Boltz setup`, `2/2 refold + score`, one line of metrics per species (`iptm`, `ipsae`, `pass` and so on), then `SMOKE DONE`.

```powershell
python tnf_smoke.py
```

Step 5: small test run. What it does: starts the cloud controller with 1 shard. What you should see next: `CONTROLLER STARTED IN THE CLOUD` with a call id. You can then close PowerShell.

```powershell
python tnf_client.py start --tag tnf --max-shards 1 --per-round 1 --goal 1 --budget 12
```

Step 6: watch it. `status` shows progress and `logs` shows the log tail.

```powershell
python tnf_client.py status --tag tnf
python tnf_client.py logs --tag tnf
```

Step 7: full run (only after Step 5 looks healthy).

```powershell
python tnf_client.py start --tag tnf --max-shards 14 --per-round 3 --goal 24 --budget 85
```

Step 8: download results. What you should see next: `wrote tnf_designs.csv`, `tnf_all_scored.csv`, `tnf_top20.fasta` and `tnf_progress.json`, or `not yet: <file>` for any that isn't ready.

```powershell
python tnf_client.py fetch --tag tnf
```

## Notes and caveats

- Third-party tools are not included and must be obtained under their own licences: BindCraft2, Boltz-2, ProteinMPNN and PROPKA. The Modal apps install them in the cloud image.
- Only scripts that load on their own are included. I imported every Modal script with the Modal client and checked that each app defines its functions, and I ran `score_ph_switch.py` on a real predicted complex. The Modal pipelines have not been run end to end from here, so do the cheap self-check and smoke test first.
- Earlier EGFR scripts (the first design shards, the first refold script, a campaign controller and a results-ingest script) depended on files in a private repository or on one person's local paths, so they are not included. `modal_bc2_dual.py` supersedes the first EGFR design script. The Stage 1 numbers in the notes came from those earlier scripts.
- The notes in `egfr/notes/` and `tnf/notes/` are the methodology texts as submitted.
- Scripts mentioned in the TNF notes that are not in this repository (for example `bc2_campaigns.py`, `refold_multichain.py` and `bc2_submission.py`) live in a separate private repo.
