# Block 2: apply the run settings, write the hotspot/coldspot file, and upload + validate the targets in Modal.
# Needs, in $env:USERPROFILE\Downloads\TfR1_pdfa : tfr1_app.py, and the folder tfr1_inputs with the two cropped apical PDBs and two cropped MSAs.
# Replace the $json text with your own spots.json (format in README.md). It is not published here while the challenge is open.
$d = "$env:USERPROFILE\Downloads\TfR1_pdfa"
$f = "$d\tfr1_app.py"
$t = Get-Content $f -Raw
$t = $t -replace 'BINDER_LENGTHS = \[\d+, \d+\]','BINDER_LENGTHS = [60, 100]'
$t = $t.Replace('MAX_SHARDS = 14','MAX_SHARDS = 24')
if ($t -notmatch 'min_plddt_screen') { $t = $t.Replace('"min_monomer_plddt_final": 0.70,','"min_monomer_plddt_final": 0.70, "min_plddt_screen": 0.4, "min_plddt_refine": 0.5,') }
Set-Content $f $t -Encoding utf8
Select-String -Path $f -Pattern '^BINDER_LENGTHS|^MAX_SHARDS|min_plddt_screen'
$json = '<PASTE YOUR spots.json TEXT HERE, ON ONE LINE>'
[System.IO.File]::WriteAllText("$d\tfr1_inputs\spots.json", $json, (New-Object System.Text.UTF8Encoding($false)))
Set-Location $d
python -m modal run tfr1_app.py::prep
