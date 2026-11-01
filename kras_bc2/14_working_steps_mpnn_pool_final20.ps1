# KRAS G12D - ProteinMPNN redesign results, pooling and final-20 selection. Steps that ran and worked.
# Replace C:\Users\HP with your own user folder. No sequences or keys are stored in this file.
# The redesign + Boltz-2 refold + ipSAE code lives in kras_app.py (ProteinMPNN with hotspot + Asp6 residues frozen).

# ---- STEP 1: launch the redesign/refold (cloud, detached). Not re-run here; set the ranked CSV you used.
# Set-Location "C:\Users\HP\Downloads\kras_modal"
# python -m modal run --detach kras_app.py::mpnn_hot_start --csv-path "<ranked csv>" --rows "1,2" --tag kras_mpnn3

# ---- STEP 2: download the results
Set-Location "C:\Users\HP\Downloads\kras_modal"
python -m modal run kras_app.py::mpnn_get --outdir "C:\Users\HP\Downloads\kras_final" --tag kras_mpnn3

# ---- STEP 3: pool every run, dedupe by sequence (most-sampled files first), rank by strict d_worst
$dir = "C:\Users\HP\Downloads\kras_final"
$inv = [Globalization.CultureInfo]::InvariantCulture
function Pick($r, $names) { foreach ($n in $names) { $p = $r.PSObject.Properties[$n]; if ($p -and "$($p.Value)" -ne "") { return $p.Value } }; return $null }
function Num($v) { $o = 0.0; if ($null -ne $v -and [double]::TryParse("$v", [Globalization.NumberStyles]::Float, $inv, [ref]$o)) { return $o } else { return $null } }
$prio = @('kras_mpnn3_stage2_all','kras_mpnn1_stage2_all','kras_rescore1','kras_rescore2','kras_FINAL_top20','kras_mpnn3_stage1_all','kras_mpnn1_stage1_all','kras_all_designs','kras_set5_designs','kras_set2_designs','kras_rest45','kras_submit_top20_v3_no_diversity','kras_submit_top20_v2','kras_submit_top20')
$files = Get-ChildItem $dir -Filter *.csv | Where-Object { $_.Name -notmatch 'ALL_ranked|variants|top20_new|submit_top20_v[45]|FINAL20' }
$files = $files | Sort-Object { $b = [IO.Path]::GetFileNameWithoutExtension($_.Name); $ix = [array]::IndexOf($prio, $b); if ($ix -lt 0) { 99 } else { $ix } }
$seen = @{}; $keep = @()
foreach ($f in $files) {
  foreach ($r in (Import-Csv $f.FullName)) {
    $g  = Num (Pick $r @('g_ipsae','g_ips'))
    $wt = Num (Pick $r @('wt_ipsae_max','wt_ips'))
    $pf = Num (Pick $r @('pf_ipsae_max','pf_ips'))
    $seq = Pick $r @('sequence','seq','binder_seq')
    if ($null -eq $g -or $null -eq $wt -or $null -eq $pf -or -not $seq) { continue }
    if ($seen.ContainsKey($seq)) { continue }
    $seen[$seq] = 1
    $dw = Num (Pick $r @('d_worst'))
    if ($null -eq $dw) { $dw = [math]::Round($g - [math]::Max($wt,$pf), 4) }
    $keep += [pscustomobject]@{
      source=$f.Name; id=(Pick $r @('design_id','id','name','variant')); parent=(Pick $r @('parent','seed','from'))
      len=$seq.Length; g_iptm=(Num (Pick $r @('g_iptm','g_ipt'))); g_ipsae=$g; wt=$wt; pf=$pf
      d_worst=$dw; pep_delta=(Num (Pick $r @('pep_d','pep_delta'))); asp6=(Num (Pick $r @('asp6_min_dist','asp6'))); sequence=$seq }
  }
}
$ranked = $keep | Where-Object { $_.g_ipsae -ge 0.5 } | Sort-Object d_worst -Descending
$i = 0; $ranked = $ranked | ForEach-Object { $i++; $_ | Add-Member rank $i -PassThru }
$ranked | Export-Csv "$dir\kras_ALL_ranked_v2.csv" -NoTypeInformation

# ---- STEP 4: final 20 = tier A (strict) then tier B (relaxed fillers)
$r = Import-Csv "$dir\kras_ALL_ranked_v2.csv"
$t1 = $r | Where-Object { [double]$_.g_ipsae -ge 0.7 -and [double]$_.pf -le 0.6 -and [double]$_.wt -le 0.8 } | Sort-Object { [double]$_.d_worst } -Descending
$ids = @{}; foreach ($x in $t1) { $ids[$x.sequence] = 1 }
$t2 = $r | Where-Object { -not $ids.ContainsKey($_.sequence) -and [double]$_.g_ipsae -ge 0.6 -and [double]$_.pf -le 0.75 -and [double]$_.wt -le 0.85 } | Sort-Object { [double]$_.d_worst } -Descending
$t1 | ForEach-Object { $_ | Add-Member tier "A_strict" -Force -PassThru } | Out-Null
$t2 | ForEach-Object { $_ | Add-Member tier "B_relaxed" -Force -PassThru } | Out-Null
$final = @($t1) + @($t2) | Select-Object -First 20
$i = 0; $final = $final | ForEach-Object { $i++; $_.rank = $i; $_ }
$out = "C:\Users\HP\Downloads\kras_FINAL20_with_sequences.csv"
$final | Select-Object rank, tier, id, parent, len, sequence, g_iptm, g_ipsae,
  @{n='wt_ipsae_max';e={$_.wt}}, @{n='pf_ipsae_max';e={$_.pf}}, d_worst, pep_delta, asp6, source |
  Export-Csv $out -NoTypeInformation -Encoding UTF8
"Final 20 saved: $out"
