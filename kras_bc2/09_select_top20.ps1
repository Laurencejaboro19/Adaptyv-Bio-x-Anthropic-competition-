# Block 9: pool every refolded design (set2 + set5), keep unique sequences with G12D ipSAE >= 0.5, rank by delta
# (G12D ipSAE minus the higher of WT and peptide-free), tiebreak closer Asp6. Writes the top 20 with sequences (kept local, not in this repo).
# History: version 1 skipped designs >=80% identical to one already picked, version 2 used 95%, version 3 (this file) has no identity filter.
# All three gave the same 20 sequences. Pooled 106 scored designs, 65 unique with G12D ipSAE >= 0.5.
Set-Location "C:\Users\HP\Downloads\kras_modal"
$out = "C:\Users\HP\Downloads\kras_final"
python -m modal run kras_app.py::fetch --outdir $out --tags "kras_set2,kras_set5" 2>&1 | Select-String -Pattern 'wrote|not yet'
$all = Import-Csv "$out\kras_all_designs.csv" | Where-Object { $_.g_ipsae -ne '' -and $_.sequence -ne '' }
Write-Host "all scored designs:" @($all).Count
$cand = $all | Sort-Object sequence -Unique | Where-Object { [double]$_.g_ipsae -ge 0.5 }
Write-Host "unique with G12D ipSAE >= 0.5:" @($cand).Count
$top = $cand | Sort-Object @{e={[double]$_.delta};Descending=$true}, @{e={ if ($_.asp6_min_dist -eq '') {99} else {[double]$_.asp6_min_dist} };Descending=$false} | Select-Object -First 20
$rank = 0
$final = $top | ForEach-Object { $rank++; $_ | Select-Object @{n='rank';e={$rank}}, campaign, design_id, binder_len, sequence, g_iptm, g_ipsae, wt_ipsae_max, pf_ipsae_max, delta, asp6_min_dist, hotspots_hit, wt_site_overlap }
$final | Export-Csv "$out\kras_submit_top20_v3_no_diversity.csv" -NoTypeInformation -Encoding utf8
$final | Select-Object rank, campaign, binder_len, g_ipsae, wt_ipsae_max, pf_ipsae_max, delta, asp6_min_dist, hotspots_hit | Format-Table -AutoSize
