# Block 3: cheap check of the Boltz-2 + ipSAE path. Refolds a made-up 78-aa sequence on G12D, WT and peptide-free and prints the combined row.
# A harmless StopAsyncIteration traceback from the Modal client may print at the end.
Set-Location "$env:USERPROFILE\Downloads\kras_modal"
python -m modal run kras_app.py::refold_test
