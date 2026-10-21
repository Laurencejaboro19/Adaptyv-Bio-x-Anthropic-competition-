# Block 3: launch the bridge campaign (30-50 aa) in the cloud. The PC can be switched off once it prints
# "CONTROLLER STARTED IN THE CLOUD". Stops at 14 good designs, 12 shards, the $45 estimated cap, or 20 hours.
cd "$env:USERPROFILE\Downloads\glp1r_modal"
python -m modal run --detach glp1r_app.py::start --camp bridge --budget 45
