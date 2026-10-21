# Block 2: smoke test. 1 H100 shard x 30 BindCraft2 trajectories, then Boltz-2 refold + ipSAE on what it produced.
# It runs in the cloud, so the PC can be switched off once it prints a function call. Set $camp to "deep" for the deep site.
$camp = "bridge"
cd "$env:USERPROFILE\Downloads\glp1r_modal"
python -m modal run --detach glp1r_app.py::smoke --camp $camp
