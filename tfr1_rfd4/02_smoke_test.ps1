# Block 2: smoke test with the authors' own PD-L1 example (2 designs, eager mode, no metrics). About 24 s of generation after the container starts.
Set-Location "$env:USERPROFILE\Downloads\tfr1_rfd4"
$env:PYTHONUTF8="1"; $env:PYTHONIOENCODING="utf-8"
python -m modal run rfd4_app.py::run_smoke --n 2
