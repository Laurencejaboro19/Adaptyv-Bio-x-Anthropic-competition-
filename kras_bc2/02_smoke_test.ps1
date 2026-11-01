# Block 2: smoke test, 1 shard x 30 trajectories then refold on all 3 targets. Runs in the cloud (--detach), PC can be off.
# Change the tag in the smoke() entrypoint (kras_smoke9999_) to rerun: a used tag resumes saved state and does nothing.
Set-Location "$env:USERPROFILE\Downloads\kras_modal"
python -m modal run --detach kras_app.py::smoke --camp set1
