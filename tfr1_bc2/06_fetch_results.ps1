# Block 6: download the result files into Downloads\TfR1_pdfa\tfr1_pB_results.
Set-Location "$env:USERPROFILE\Downloads\TfR1_pdfa"
$env:PYTHONUTF8="1"; $env:PYTHONIOENCODING="utf-8"
python -m modal run tfr1_app.py::fetch --outdir "$env:USERPROFILE\Downloads\TfR1_pdfa\tfr1_pB_results" --tag tfr1_pB
Get-ChildItem "$env:USERPROFILE\Downloads\TfR1_pdfa\tfr1_pB_results" -Recurse -File | Select-Object -First 20 -ExpandProperty FullName
