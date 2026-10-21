# Block 6: download the results for bridge and deep, print the passing designs, and save them to their own CSV.
$out = "$env:USERPROFILE\Downloads\glp1r_final_results"
cd "$env:USERPROFILE\Downloads\glp1r_modal"
python -m modal run glp1r_app.py::fetch --outdir $out --tags "glp1r_bridge,glp1r_deep"
$csv = Import-Csv "$out\glp1r_bc2_refolded_ipsae.csv"
"TOTAL ROWS: " + $csv.Count
$good = $csv | Where-Object { $_.good -eq 'True' }
"GOOD (bridge + deep): " + @($good).Count
$good | Export-Csv "$out\glp1r_good_designs_bridge_and_deep.csv" -NoTypeInformation -Encoding UTF8
$good | Format-Table campaign, rank, binder_len, iptm, ipsae, hotspots_hit, hotspots_total, coldspots_hit, cterm_in_interface, contact_PASS -AutoSize | Out-String -Width 250
"saved to $out\glp1r_good_designs_bridge_and_deep.csv"
