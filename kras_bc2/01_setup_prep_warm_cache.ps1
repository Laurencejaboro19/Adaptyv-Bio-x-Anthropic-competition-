# Block 1: put kras_app.py and the 3 target PDBs in $env:USERPROFILE\Downloads (kras_app.py goes in a kras_modal subfolder),
# then upload the PDBs to the Modal volume and cache the Boltz-2 weights + the HLA and b2m MSAs (done once).
# prep prints "PREP OK" and the chain/peptide check for the 3 files. warm_cache prints "warm OK | hla: MSA saved ... | b2m: MSA saved ...".
Set-Location "$env:USERPROFILE\Downloads\kras_modal"
python -m modal run kras_app.py::prep --indir "$env:USERPROFILE\Downloads"
python -m modal run kras_app.py::warm_cache
