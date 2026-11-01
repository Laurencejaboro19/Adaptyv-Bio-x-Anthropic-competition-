# Block 12: rescore the 45 designs (G12D ipSAE >= 0.5) that were not in the first 20, same 8-sample method, detached.
# Ran 2026-11-01: 45 designs, 135 folds. Expected cost about $6 (estimated from the first rescore: 60 folds ~ $2.72). Needs block 10 run first.
$d = "C:\Users\HP\Downloads\kras_final"
$all = Import-Csv "$d\kras_all_designs.csv" | Where-Object { $_.g_ipsae -ne '' -and $_.sequence -ne '' }
$cand = $all | Sort-Object sequence -Unique | Where-Object { [double]$_.g_ipsae -ge 0.5 }
$done = (Import-Csv "$d\kras_rescore1.csv").sequence
$rest = $cand | Where-Object { $done -notcontains $_.sequence }
$rest | Select-Object design_id, sequence, campaign | Export-Csv "$d\kras_rest45.csv" -NoTypeInformation -Encoding utf8
Write-Host "remaining designs to rescore:" @($rest).Count
Set-Location "C:\Users\HP\Downloads\kras_modal"
python -m modal run --detach kras_app.py::rescore_start --csv-path "$d\kras_rest45.csv" --tag kras_rescore2 --samples 8 --top 60 --max-usd 40
