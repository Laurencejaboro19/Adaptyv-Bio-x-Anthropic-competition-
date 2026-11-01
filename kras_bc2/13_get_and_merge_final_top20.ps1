# Block 13: download the second rescore, merge it with the first 20, keep designs with mean G12D ipSAE >= 0.5 and rank all 65 by
# score = (delta_mean + delta_worst) / 2, so a design needs a good gap on both the mean and the worst sample. Saves the final top 20.
Set-Location "C:\Users\HP\Downloads\kras_modal"
python -m modal run kras_app.py::rescore_get --outdir "C:\Users\HP\Downloads\kras_final" --tag kras_rescore2 | Select-String -Pattern 'state|not finished|saved'
$d = "C:\Users\HP\Downloads\kras_final"
if (Test-Path "$d\kras_rescore2.csv") {
  $m = @(Import-Csv "$d\kras_rescore1.csv") + @(Import-Csv "$d\kras_rescore2.csv")
  $m = $m | Where-Object { $_.delta_mean -ne '' -and [double]$_.g_ipsae -ge 0.5 }
  foreach ($x in $m) { $x | Add-Member -NotePropertyName score -NotePropertyValue (([double]$x.delta_mean + [double]$x.delta_worst) / 2) -Force }
  $top = $m | Sort-Object { [double]$_.score } -Descending | Select-Object -First 20
  $r = 0
  $final = $top | ForEach-Object { $r++; $_ | Select-Object @{n='final_rank';e={$r}}, campaign, design_id, binder_len, sequence, g_iptm, g_ipsae, wt_ipsae, pf_ipsae, delta_mean, delta_worst, score, pep_delta, g_asp6_min }
  $final | Export-Csv "$d\kras_FINAL_top20.csv" -NoTypeInformation -Encoding utf8
  $final | Select-Object final_rank, binder_len, g_ipsae, wt_ipsae, pf_ipsae, delta_mean, delta_worst, pep_delta, g_asp6_min | Format-Table -AutoSize
  Write-Host "saved:" "$d\kras_FINAL_top20.csv"
}
