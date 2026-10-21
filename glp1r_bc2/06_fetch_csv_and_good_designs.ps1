# Block 6: download the results CSV, print the top rows, and save the passing designs to their own CSV.
$out = "$env:USERPROFILE\Downloads\glp1r_bc2_results"
cd "$env:USERPROFILE\Downloads\glp1r_modal"
python -m modal run glp1r_app.py::fetch --outdir $out --tags glp1r_bridge
$csv = Import-Csv "$out\glp1r_bc2_refolded_ipsae.csv"
"TOTAL ROWS: " + $csv.Count
"GOOD: " + @($csv | Where-Object { $_.good -eq 'True' }).Count
$csv | Select-Object -First 15 | Format-Table rank, good, binder_len, iptm, ipsae, complex_plddt, hotspots_hit, hotspot_frac, coldspots_hit, n_binder_contacts, cterm_in_interface, contact_PASS -AutoSize | Out-String -Width 250
$good = $csv | Where-Object { $_.good -eq 'True' }
$good | Export-Csv "$out\glp1r_bridge_good_designs.csv" -NoTypeInformation -Encoding UTF8
"saved " + @($good).Count + " passing designs to $out\glp1r_bridge_good_designs.csv"
