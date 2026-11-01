# Block 12a: rescore the 45 pooled designs (G12D ipSAE >= 0.5) that were not in the first rescore. Detached, PC can be off.
# Ran 2026-11-01: 45 designs, 135 folds, 8 samples each. Cost not recorded here (estimate was about $6 from block 10's rate; check the Modal dashboard).
$d = "$env:USERPROFILE\Downloads\kras_final"
$all = Import-Csv "$d\kras_all_designs.csv" | Where-Object { $_.g_ipsae -ne '' -and $_.sequence -ne '' }
$cand = $all | Sort-Object sequence -Unique | Where-Object { [double]$_.g_ipsae -ge 0.5 }
$done = (Import-Csv "$d\kras_rescore1.csv").sequence
$rest = $cand | Where-Object { $done -notcontains $_.sequence }
$rest | Select-Object design_id, sequence, campaign | Export-Csv "$d\kras_rest45.csv" -NoTypeInformation -Encoding utf8
Write-Host "remaining designs to rescore:" @($rest).Count
Set-Location "$env:USERPROFILE\Downloads\kras_modal"
python -m modal run --detach kras_app.py::rescore_start --csv-path "$d\kras_rest45.csv" --tag kras_rescore2 --samples 8 --top 60 --max-usd 40

# Block 12b (run later): download the second rescore, merge it with the first (65 designs), rank by
# score = (delta_mean + delta_worst) / 2, keep the top 20. Prints "not finished yet" until the cloud run is done.
# python -m modal run kras_app.py::rescore_get --outdir "$env:USERPROFILE\Downloads\kras_final" --tag kras_rescore2
# $m = @(Import-Csv "$d\kras_rescore1.csv") + @(Import-Csv "$d\kras_rescore2.csv")
# $m = $m | Where-Object { $_.delta_mean -ne '' -and [double]$_.g_ipsae -ge 0.5 }
# foreach ($x in $m) { $x | Add-Member -NotePropertyName score -NotePropertyValue (([double]$x.delta_mean + [double]$x.delta_worst) / 2) -Force }
# $top = $m | Sort-Object { [double]$_.score } -Descending | Select-Object -First 20
# $r = 0
# $final = $top | ForEach-Object { $r++; $_ | Select-Object @{n='final_rank';e={$r}}, campaign, design_id, binder_len, sequence, g_iptm, g_ipsae, wt_ipsae, pf_ipsae, delta_mean, delta_worst, score, pep_delta, g_asp6_min }
# $final | Export-Csv "$d\kras_FINAL_top20.csv" -NoTypeInformation -Encoding utf8
# $final | Select-Object final_rank, binder_len, g_ipsae, wt_ipsae, pf_ipsae, delta_mean, delta_worst, pep_delta, g_asp6_min | Format-Table -AutoSize
