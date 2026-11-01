# Block 10: rescore the top 20 from block 9 with 8 Boltz-2 samples per target, in the cloud (PC can be off).
# kras_app.py in this folder already contains the rescore code (fold_n, rescore_ctl, rescore_start, rescore_get).
# Per design and target it gives mean/min/max ipSAE, mean ipTM, peptide-chain ipSAE, closest Asp6 distance, and
# delta_mean (G12D mean minus higher avoid mean), delta_worst (G12D mean minus the worst avoid sample), pep_delta (G12D minus WT, peptide chain).
# Ran 2026-11-01: 20 designs, 60 folds, finished in 8.8 min, compute est $2.72.
Set-Location "$env:USERPROFILE\Downloads\kras_modal"
python -m modal run --detach kras_app.py::rescore_start --csv-path "$env:USERPROFILE\Downloads\kras_final\kras_submit_top20_v3_no_diversity.csv" --tag kras_rescore1 --samples 8 --top 20 --max-usd 15
