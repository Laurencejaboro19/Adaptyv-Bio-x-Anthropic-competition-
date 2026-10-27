# Block 5: read per_example_metrics.csv. With '~metrics' (as run here) this file holds run settings only (model, seed, steps), NOT quality scores.
$tag = "my_first_run"
$csv = "$env:USERPROFILE\Downloads\tfr1_rfd4\out\$tag\inference_metrics\per_example_metrics.csv"
$rows = Import-Csv $csv
"Rows: $($rows.Count)"
"Columns:"; $rows[0].PSObject.Properties.Name -join ", "
$rows | Format-List
