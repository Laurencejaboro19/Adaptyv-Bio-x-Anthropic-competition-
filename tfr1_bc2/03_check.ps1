# Block 3: plumbing test (about 10 minutes, under $1). Refolds a throwaway control binder on human and mouse and scores it.
# Keep PowerShell open and the PC on until it prints CHECK DONE. A failing score for the control is expected; an error is not.
Set-Location "$env:USERPROFILE\Downloads\TfR1_pdfa"
python -m modal run tfr1_app.py::check
