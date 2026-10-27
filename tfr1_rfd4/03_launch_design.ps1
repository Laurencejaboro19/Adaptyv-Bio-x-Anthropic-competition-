# Block 3: design binders against your own target. Detached, so the PC can be switched off once it prints the run link.
# Fill in the four variables. Target PDB: one chain, real residue numbers. $hot = comma-separated residue numbers on the target (no chain letter).
$dir    = "$env:USERPROFILE\Downloads\tfr1_rfd4"
$tag    = "my_first_run"                          # letters, digits, underscore; used for the output folder name
$pdb    = "$dir\target.pdb"                       # put your target PDB here
$hot    = ""                                      # e.g. "12,15,40,41"
$binder = "60-100"                                # binder length range
$range  = ""                                      # target residue range exactly as in the PDB, e.g. "191-383" (must have no missing residues)
$n      = 4                                       # designs to make
if (-not ($hot -and $range -and (Test-Path $pdb))) { Write-Host "Set `$hot, `$range and put target.pdb in $dir first" -ForegroundColor Red; return }
Set-Location $dir
$env:PYTHONUTF8="1"; $env:PYTHONIOENCODING="utf-8"
python -m modal run --detach rfd4_app.py::run_design --tag $tag --pdb $pdb --hot $hot --contig-binder $binder --target-range $range --n $n
