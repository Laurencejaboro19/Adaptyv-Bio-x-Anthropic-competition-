# Block 11: check a rescore run and download its table. Prints "not finished yet" until done. Columns: g_/wt_/pf_ = mean over 8 samples,
# d_mean = G12D ipSAE minus the higher avoid mean, d_worst = G12D mean minus the worst avoid sample, pep_d = G12D peptide-chain ipSAE minus WT's.
Set-Location "C:\Users\HP\Downloads\kras_modal"
python -m modal run kras_app.py::rescore_get --outdir "C:\Users\HP\Downloads\kras_final" --tag kras_rescore1
