# RFD4-Proteina on Modal (generator for the next challenge)

This section installs and runs **RFD4-Proteina** (`nvidia/RFD4-Proteina`, the RosettaCommons/NVIDIA latent flow-matching generator, BSD-3 licence for code and weights) on Modal from Windows PowerShell. It is a **generator only**: it makes backbones and sequences. It does not score them.

- **What ran:** the install and weights download, the authors' own PD-L1 example (2 designs), and a 4-design test against a 193-residue cropped target with 5 hotspot residues. All three finished with exit code 0 and wrote their `.cif` files, then the files were downloaded and read in PowerShell.
- **What is in this folder:** `rfd4_app.py` (the Modal app) and the PowerShell blocks that were run on a Windows PC, in order. No design sequences, target residue numbers or results are included, because the TfR1 challenge was still open when this was written.
- **What was NOT done:** no refolding, no ipSAE, no pH filter on RFD4 output. The metrics CSV from a `~metrics` run holds run settings only (see Block 5). Quality scoring is the Boltz-2 -> ipSAE -> PROPKA steps used elsewhere in this repo. Design quality from RFD4 is therefore **unknown**.

Run everything in **Windows PowerShell opened from the app**. Paths are absolute under `Downloads`, so each block can be pasted from anywhere. Modal must be installed and logged in (see the top-level README, Step 0).

## Settings used

- **Code and weights pinned:** repo commit `aaccdd2ad0720c180468df5f9cf037a713d403e5`, Hugging Face revision `a06e2509c7033b010fc43614e9fd840303588502`.
- **Install:** mirrors the repo's own Dockerfile: Python 3.12, `uv`, the `all-cuda12` extra with Torch `cu126` (the safer profile; `RFD4_ACCEL=cu13` is untested here), `tmol` from wheels only. `gemmi` is added to convert a PDB target to mmCIF.
- **GPU:** L40S by default (`RFD4_GPU` overrides; A100 and H100 are listed as tested by the authors but were not run here).
- **Run:** `rfd4-generate` with `compile_lengths=null` (eager, no compile wait), `~metrics`, `batch_size=1`.
- **Target conditioning (what the app writes for you):** `C_CRDGRP` target backbone only, `C_SEQ` full target sequence, `C_HOT` on the hotspot residues you give, and a contig of the form `<binder range>, /0, A<target range>`.
- **Cost control:** each function has a time limit, not a dollar limit: setup 1 h, smoke 3 h, design 4 h. At roughly $2/h for an L40S (check Modal pricing) 4 h is about $8. Observed: 2 designs took about 24 s of generation; the 4-design TfR1 test took a few minutes.

## Two bugs hit and fixed (both are already in `rfd4_app.py`)

1. `CIF polymer label_asym_id 'Axp' has missing label_seq_id values`: gemmi's PDB -> mmCIF conversion left `label_seq_id` empty. Fixed with `assign_label_seq_id`.
2. `Preserved segment 'A191-383' matched no atoms in the input structure`: gemmi had renamed the chain to `Axp` and renumbered residues. RFD4 matches the contig against the chain and residue numbers, so the app now copies `auth_asym_id` and `auth_seq_id` into the label columns. Checked by loading the converted file with the atomworks parser (chain `A`, original residue numbers).

Input constraints worth knowing: the spec reader only accepts `.cif`, `.cif.gz` or `.bcif`, so a PDB must be converted (the app does it), and the contig `A<start>-<end>` assumes the target residues in that range are contiguous.

## Run order

### Block 1: write the app, build the image, download weights (`01_write_app_and_setup.ps1`)
**What it does:** creates `Downloads\tfr1_rfd4`, writes `rfd4_app.py` there, builds the cloud image, downloads the weights to a Modal volume (`rfd4-cache`) and runs `rfd4-generate --help`.
**What you should see:** a first run builds the image (not timed here), then `torch ... cuda ok: True NVIDIA L40S`, the weights listed (about 1.7 GB and 1.0 GB), and `{'help_rc': 0}`. Reruns reuse the cached image and weights.

### Block 2: smoke test (`02_smoke_test.ps1`)
**What it does:** runs the authors' PD-L1 example for 2 designs.
**What you should see:** `Generation complete`, `rc: 0`, and `OUTPUT FILES: 3` (2 structures plus a metrics CSV). Harmless warnings about `CCD_MIRROR_PATH`, `PDB_MIRROR_PATH` and "sentinel insertion" appear.

### Block 3: design against your own target (`03_launch_design.ps1`)
**What it does:** uploads your target PDB and launches a detached run (the PC can be switched off once the run link prints). Fill in `$tag`, `$pdb`, `$hot`, `$range` first; the block refuses to run until they are set.
**What you should see:** the generated YAML, then (in the Modal log) a progress bar and `rc: 0` with the output file count. The 4-design test finished in a few minutes.
**If it fails:** get the error with `python -m modal app logs <app id>` and look for the last `Error` lines; both failures above were found this way.

### Block 4: download and unpack (`04_fetch_and_unpack.ps1`)
**What you should see:** `saved ...tar.gz`, then the paths of the `*_generated.cif` files and `per_example_metrics.csv`.

### Block 5: read the metrics CSV (`05_metrics_csv.ps1`)
**What you should see:** one row per design and columns that are run settings (`example_id`, `ckpt_path`, `seed`, `generation.args.nsteps`, ...). **There are no ipTM, pLDDT or ipSAE values**, because generation ran with `~metrics`.

### Block 6: chains, lengths and sequences (`06_sequences_and_chains.ps1`)
**What it does:** reads each generated `.cif` and prints, per chain, the residue count, residue-number range and sequence.
**What you should see:** the target chain `A` renumbered `1..N` (a 193-residue crop gives 193), and the designed binder as another chain with the length you asked for.

## Notes from the 4-design test

RFD4 does not know about pH, so histidines for a pH switch are not designed in. Generated sequences can also look low-complexity or antibody-like, so check each sequence for novelty and fold plausibility (refold it) before using it.
