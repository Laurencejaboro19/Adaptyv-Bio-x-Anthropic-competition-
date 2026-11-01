# Block 4: check a run. status prints progress.json for each tag listed in status().
Set-Location "$env:USERPROFILE\Downloads\kras_modal"
python -m modal run kras_app.py::status
# last 25 lines of each shard's BindCraft2 log for a tag:
python -m modal run kras_app.py::logs --tag kras_set5
# list the files a tag wrote:
python -m modal run kras_app.py::inspect_out --tag kras_set5
