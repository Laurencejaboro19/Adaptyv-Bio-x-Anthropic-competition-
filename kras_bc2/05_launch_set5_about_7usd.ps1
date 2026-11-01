# Block 5: THE RUN THAT WORKED FOR ABOUT $7 (estimate): bind G12D, avoid wild-type and peptide-free, every BC2 stage floor 0.3,
# peptide-only hotspots (C4-C8), forced_targeting off. 2 shards, 30 trajectories each. Cap $50 (the controller stopped on its own at est $7.27).
# Only launch if status shows "no progress file" for kras_set5, otherwise this resumes the saved state.
Set-Location "$env:USERPROFILE\Downloads\kras_modal"
python -m modal run --detach kras_app.py::start --camp set5 --budget 50 --goal 3 --max-shards 2 --per-round 2 --max-traj 30
