# Block 5: progress. Reads from Modal, starts no GPU work.
$env:PYTHONUTF8 = "1"
$env:PYTHONIOENCODING = "utf-8"
Set-Location "$env:USERPROFILE\Downloads\TfR1_pdfa"
python -m modal run tfr1_app.py::status --tags tfr1_pB
