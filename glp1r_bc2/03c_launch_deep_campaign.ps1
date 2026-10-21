# Block 3c: launch the deep campaign (35-60 aa) in the cloud. The PC can be switched off once it prints
# "CONTROLLER STARTED IN THE CLOUD". Stops at 10 good designs, 12 shards, the $25 estimated cap, or 20 hours.
cd "$env:USERPROFILE\Downloads\glp1r_modal"
python -m modal run --detach glp1r_app.py::start --camp deep --goal 10 --max-shards 12 --per-round 4 --budget 25
