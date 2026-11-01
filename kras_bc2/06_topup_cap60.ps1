# Block 6: top-up on the same tag (resumes saved state; only works once the first run's state is "done").
# Raises the total cap to $60 (estimate), asks for up to 14 shards in rounds of 6, goal 999 so it runs until the cap.
Set-Location "$env:USERPROFILE\Downloads\kras_modal"
python -m modal run --detach kras_app.py::start --camp set5 --budget 60 --goal 999 --max-shards 14 --per-round 6 --max-traj 30
