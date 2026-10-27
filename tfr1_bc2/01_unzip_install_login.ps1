# Block 1: unzip the TfR1 package, set binder length and budget, install the Modal client, log in to Modal.
# Needs TfR1_pdfa.zip (tfr1_app.py + tfr1_inputs) in Downloads. The browser opens once to approve the Modal login.
$d = "$env:USERPROFILE\Downloads"
Expand-Archive -Force "$d\TfR1_pdfa.zip" $d
$f = "$d\TfR1_pdfa\tfr1_app.py"
(Get-Content $f -Raw) -replace 'BINDER_LENGTHS = \[50, 100\]','BINDER_LENGTHS = [60, 120]' -replace 'BUDGET_USD = 60\.0','BUDGET_USD = 70.0' | Set-Content $f -Encoding utf8
Select-String -Path $f -Pattern '^BINDER_LENGTHS|^BUDGET_USD'
python -m pip install --upgrade modal biopython
python -m modal setup
