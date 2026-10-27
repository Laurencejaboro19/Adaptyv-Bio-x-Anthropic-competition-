# Block 4: start the controller in the cloud (detached). Wait for CONTROLLER STARTED IN THE CLOUD before switching the PC off.
Set-Location "$env:USERPROFILE\Downloads\TfR1_pdfa"
python -m modal run --detach tfr1_app.py::start --tag tfr1_pB --budget 70 --goal 24
