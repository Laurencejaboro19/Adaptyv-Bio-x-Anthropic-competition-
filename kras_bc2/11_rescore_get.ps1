# Block 11: check a rescore and download it. Prints "not finished yet" until the run is done.
Set-Location "$env:USERPROFILE\Downloads\kras_modal"
python -m modal run kras_app.py::rescore_get --outdir "$env:USERPROFILE\Downloads\kras_final" --tag kras_rescore1
