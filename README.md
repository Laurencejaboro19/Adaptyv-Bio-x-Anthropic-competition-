# Adaptyv Bio x Anthropic protein design competition: PDFA code and notes

Code, configs and notes from Protein Design for Africa (PDFA) for the Anthropic x Adaptyv competition on Proteinbase.

| Challenge | Target | Folder |
|---|---|---|
| 1 | EGFR, conditional (pH-switch) binder, human and mouse | [`egfr/`](egfr/) |
| 2 | TNF-alpha, conditional (pH-switch) binder, human and mouse | [`tnf/`](tnf/) |

The design approach for both is in [`DESIGN_STRATEGY.md`](DESIGN_STRATEGY.md).

**All evidence is computational.** No wet-lab validation has been done. The pH preference is a PROPKA3 / Wyman-linkage estimate on a predicted structure, used to rank designs. It is not a measured affinity ratio.

## Pipeline

BindCraft2 (design) -> Boltz-2 (refold, source of truth) -> PROPKA3 + Wyman linkage (pH switch) -> ranking. Everything ran on [Modal](https://modal.com) and was launched from Windows PowerShell.

## Layout

```
egfr/
  modal/      Modal apps and the two EGFR target PDBs they load (keep them in the same folder as the scripts)
  scoring/    score_ph_switch.py (PROPKA3 wrapper, Wyman linkage, His geometry) and epitope analysis
  configs/    BindCraft2 campaign JSONs, controller config, measured funnel and cost records
  notes/      EGFR_methodology.md
tnf/
  modal/      tnf_app.py (Modal app), tnf_client.py (launcher), tnf_smoke.py, and the four TNF-alpha target PDBs
  notes/      methodology, configs and design/validation notes
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

After this, refold the candidates with Boltz-2 (`modal_bc2_egfr_refold.py`) and score the pH switch with `egfr/scoring/score_ph_switch.py`. The stage-2 controller `egfr_dual10_stage2_controller.py` automates design, refold and scoring with a budget cap (`BUDGET_USD` inside it).

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

- Third-party tools are not included and must be obtained under their own licences: BindCraft2, Boltz-2, ProteinMPNN and PROPKA. Some Modal scripts read an adapter path from the `RB_ADAPTER_DIR` environment variable.
- The notes in `egfr/notes/` and `tnf/notes/` are the methodology texts as submitted.
- The EGFR scripts in `egfr/modal/` are the versions kept in the submission bundle (see its manifest) plus the dual-species stage-2 pipeline. `modal_bc2_dual.py` supersedes `modal_bc2_egfr.py` for the dual campaign.
- Scripts mentioned in the notes that are not in this repository (for example `bc2_campaigns.py`, `refold_multichain.py` and `bc2_submission.py`) live in a separate private repo.
