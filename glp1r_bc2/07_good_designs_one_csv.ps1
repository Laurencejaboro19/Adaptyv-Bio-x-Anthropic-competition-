# Block 7: one CSV in Downloads with the passing bridge and deep designs, including the passing design from the deep smoke test.
# Needs the results from Block 6 and the smoke-test results fetched with: ::fetch --outdir "$env:USERPROFILE\Downloads\glp1r_smoke4_deep_results" --tags glp1r_smoke4_deep
$main = Import-Csv "$env:USERPROFILE\Downloads\glp1r_final_results\glp1r_bc2_refolded_ipsae.csv" | Where-Object { $_.good -eq 'True' }
$main | ForEach-Object { Add-Member -InputObject $_ -NotePropertyName run -NotePropertyValue 'campaign' -Force }
$smoke = Import-Csv "$env:USERPROFILE\Downloads\glp1r_smoke4_deep_results\glp1r_bc2_refolded_ipsae.csv" | Where-Object { $_.good -eq 'True' }
$smoke | ForEach-Object { Add-Member -InputObject $_ -NotePropertyName campaign -NotePropertyValue 'deep' -Force; Add-Member -InputObject $_ -NotePropertyName run -NotePropertyValue 'smoke test' -Force }
$rows = @($main) + @($smoke)
$rows | Sort-Object campaign, @{Expression = { [double]$_.ipsae }; Descending = $true} | Select-Object campaign, run, rank, binder_len, sequence, iptm, ptm, ipsae, complex_plddt, hotspots_hit, hotspots_total, hotspot_frac, coldspots_hit, n_binder_contacts, cterm_in_interface, contact_PASS, design_id | Export-Csv "$env:USERPROFILE\Downloads\glp1r_good_designs.csv" -NoTypeInformation -Encoding UTF8
"bridge: " + @($rows | Where-Object { $_.campaign -eq 'bridge' }).Count + "   deep: " + @($rows | Where-Object { $_.campaign -eq 'deep' }).Count
"saved to $env:USERPROFILE\Downloads\glp1r_good_designs.csv"
