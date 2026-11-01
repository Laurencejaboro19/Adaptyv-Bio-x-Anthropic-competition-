# Block 7: download designs.csv (one row per refolded design, all 3 targets scored) and read the top rows.
Set-Location "$env:USERPROFILE\Downloads\kras_modal"
python -m modal run kras_app.py::fetch --outdir "$env:USERPROFILE\Downloads\kras_results5" --tags "kras_set5"
Import-Csv "$env:USERPROFILE\Downloads\kras_results5\kras_set5_designs.csv" |
  Sort-Object { [double]$_.g_ipsae } -Descending |
  Select-Object -First 20 design_id, binder_len, g_iptm, g_ipsae, wt_ipsae_max, pf_ipsae_max, delta, asp6_min_dist, tier |
  Format-Table -AutoSize
