# Block 2: smoke test. 1 H100 shard x 20 BindCraft2 trajectories, then Boltz-2 refold + ipSAE on what it produced.
# Keep the PC on until the final progress block prints.
cd "$env:USERPROFILE\Downloads\glp1r_modal"
python -m modal run --detach glp1r_app.py::smoke --camp bridge
