# Block 7: rank the passing designs and save the top 20 with every metric column as one CSV in Downloads.
$in  = "$env:USERPROFILE\Downloads\TfR1_pdfa\tfr1_pB_results\tfr1_pB_designs.csv"
$out = "$env:USERPROFILE\Downloads\trf1 20 bc2 run.csv"
$rows = Import-Csv $in
$cols = $rows[0].PSObject.Properties.Name
foreach ($r in $rows) {
  if ($cols -contains "score") { $s = [double]$r.score }
  elseif (($cols -contains "min_ipsae") -and ($cols -contains "min_his_log10_fold")) { $s = [double]$r.min_ipsae + 0.2 * [math]::Min(1.5, [double]$r.min_his_log10_fold) }
  else { Write-Host "Expected columns missing. Your columns are:" -ForegroundColor Red; $cols -join ", "; return }
  $r | Add-Member -NotePropertyName rank_score -NotePropertyValue $s -Force
}
$top = @($rows | Sort-Object { [double]$_.rank_score } -Descending | Select-Object -First 20)
$i = 0
foreach ($r in $top) { $i++; $r | Add-Member -NotePropertyName rank -NotePropertyValue $i -Force }
$first = @("rank","design_id","rank_score") + @("min_ipsae","min_his_log10_fold","length","sequence" | Where-Object { $cols -contains $_ })
$rest  = $cols | Where-Object { $_ -notin $first }
$csv = $top | Select-Object ($first + $rest) | ConvertTo-Csv -NoTypeInformation
[System.IO.File]::WriteAllLines($out, $csv, (New-Object System.Text.UTF8Encoding($true)))
"Saved: $out"
"Rows: $($top.Count)   Columns: $(($first + $rest).Count)"
