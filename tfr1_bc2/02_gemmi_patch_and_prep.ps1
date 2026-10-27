# Block 2: add gemmi to the cloud image (the first prep failed without it), then upload and validate the targets.
$f = "$env:USERPROFILE\Downloads\TfR1_pdfa\tfr1_app.py"
$t = Get-Content $f -Raw
$t = $t.Replace('"propka==3.5.1", "biopython", "numpy", "pandas")','"propka==3.5.1", "biopython", "numpy", "pandas", "gemmi")')
Set-Content $f $t -Encoding utf8
Select-String -Path $f -Pattern 'pip_install\("propka'
Set-Location "$env:USERPROFILE\Downloads\TfR1_pdfa"
python -m modal run tfr1_app.py::prep
