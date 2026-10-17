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
  modal/      Modal apps: BindCraft2 design shards, Boltz-2 co-folding, campaign controllers, cost accounting
  scoring/    score_ph_switch.py (PROPKA3 wrapper, Wyman linkage, His geometry) and epitope analysis
  configs/    BindCraft2 campaign JSONs, controller config, measured funnel and cost records
  notes/      EGFR_methodology.md
tnf/
  modal/      tnf_app.py (Modal app), tnf_client.py (launcher), tnf_smoke.py
  targets/    human and mouse TNF-alpha structures used as design targets
  notes/      methodology, configs and design/validation notes
```

## Running on Modal from PowerShell

One-time setup:

```powershell
pip install modal
python -m modal token new
```

EGFR (dual-species BindCraft2 campaign; `--detach` keeps it running after PowerShell closes, and results go to a Modal Volume):

```powershell
cd egfr\modal
python -m modal run modal_bc2_dual.py::selfcheck
python -m modal run --detach modal_bc2_dual.py::launch
python -m modal run modal_bc2_dual.py::status
python -m modal volume get pdfa-bc2-results egfr_dual .\results
```

TNF-alpha (deploy the app, then drive it with the client; run from `tnf\modal` with the PDBs from `tnf\targets` alongside):

```powershell
cd tnf\modal
python -m modal deploy tnf_app.py
python tnf_client.py check
python tnf_client.py start --tag tnf --max-shards 1 --per-round 1 --goal 1 --budget 12
python tnf_client.py status --tag tnf
python tnf_client.py fetch --tag tnf
```

Each run costs real GPU money on Modal. Do a smoke run first and set a budget.

## Notes and caveats

- Third-party tools are not included and must be obtained under their own licences: BindCraft2, Boltz-2, ProteinMPNN and PROPKA. Some Modal scripts read an adapter path from the `RB_ADAPTER_DIR` environment variable.
- The notes in `egfr/notes/` and `tnf/notes/` are the methodology texts as submitted.
- The EGFR scripts in `egfr/modal/` are the versions kept in the submission bundle (see its manifest) plus the dual-species stage-2 pipeline. `modal_bc2_dual.py` supersedes `modal_bc2_egfr.py` for the dual campaign.
- Scripts mentioned in the notes that are not in this repository (for example `bc2_campaigns.py`, `refold_multichain.py` and `bc2_submission.py`) live in a separate private repo.
