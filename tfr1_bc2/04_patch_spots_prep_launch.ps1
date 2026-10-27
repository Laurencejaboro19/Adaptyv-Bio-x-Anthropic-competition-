# Block 4: final settings, hotspot/coldspot file, target upload, and the real run (detached; wait for CONTROLLER STARTED IN THE CLOUD before switching the PC off).
# Replace the $json text with your own spots.json (format in README.md). The residue numbers are not published here while the challenge is open.
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
python -m modal run --detach tfr1_app.py::start --tag tfr1_pB --budget 70 --goal 24
