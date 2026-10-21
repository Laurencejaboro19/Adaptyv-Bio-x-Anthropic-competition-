# Block 3b (optional): resume the same campaign for more shards. The saved state is reused, so nothing is refolded twice.
# --budget is a TOTAL cap that includes what the campaign already spent.
cd "$env:USERPROFILE\Downloads\glp1r_modal"
python -m modal run --detach glp1r_app.py::start --camp bridge --goal 14 --max-shards 40 --per-round 6 --budget 70
