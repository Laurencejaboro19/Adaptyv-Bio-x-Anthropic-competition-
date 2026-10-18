$ErrorActionPreference = 'Stop'
$csv = Get-ChildItem "$HOME\Downloads\glp1r_v2_links_scores*.csv" | Sort-Object LastWriteTime -Descending | Select-Object -First 1
if (-not $csv) { Write-Host "glp1r_v2_links_scores.csv not found in Downloads"; return }
$rows = Import-Csv $csv.FullName
$aa = @{ALA='A';ARG='R';ASN='N';ASP='D';CYS='C';GLN='Q';GLU='E';GLY='G';HIS='H';ILE='I';LEU='L';LYS='K';MET='M';PHE='F';PRO='P';SER='S';THR='T';TRP='W';TYR='Y';VAL='V'}
$dir = "$HOME\Downloads\bp2_cifs"
New-Item -ItemType Directory -Force -Path $dir | Out-Null
$out = New-Object System.Collections.ArrayList
$bad = 0
foreach ($r in $rows) {
  $p = Join-Path $dir ($r.design + '.cif')
  try { Invoke-WebRequest -Uri $r.url -OutFile $p -UseBasicParsing } catch { Write-Host ("download failed: " + $r.design); $bad++; continue }
  $cols = @(); $inLoop = $false; $chains = @{}; $order = @()
  foreach ($line in Get-Content $p) {
    if ($line -like '_atom_site.*') { $cols += $line.Trim(); $inLoop = $true; continue }
    if ($inLoop -and $line -match '^(ATOM|HETATM)') {
      $f = $line -split '\s+'
      $iAtom = $cols.IndexOf('_atom_site.label_atom_id'); $iRes = $cols.IndexOf('_atom_site.label_comp_id')
      $iCh = $cols.IndexOf('_atom_site.auth_asym_id'); if ($iCh -lt 0) { $iCh = $cols.IndexOf('_atom_site.label_asym_id') }
      $iSeq = $cols.IndexOf('_atom_site.label_seq_id'); $iMod = $cols.IndexOf('_atom_site.pdbx_PDB_model_num')
      if ($iMod -ge 0 -and $f[$iMod] -ne '1') { continue }
      if ($f[$iAtom] -ne 'CA') { continue }
      $c = $f[$iCh]
      if (-not $chains.ContainsKey($c)) { $chains[$c] = New-Object System.Text.StringBuilder; $order += $c }
      $one = $aa[$f[$iRes]]; if (-not $one) { $one = 'X' }
      [void]$chains[$c].Append($one)
    }
  }
  $seq = ''
  foreach ($c in $order) { if ($chains[$c].Length -le 60) { $seq = $chains[$c].ToString() } }
  $o = [pscustomobject]@{ design=$r.design; campaign=$r.campaign; sequence=$seq; length=$seq.Length; length_expected=$r.length_expected
    iptm=$r.iptm; min_interaction_pae=$r.min_interaction_pae; binding_confidence=$r.binding_confidence; structure_confidence=$r.structure_confidence
    helix_fraction=$r.helix_fraction; sheet_fraction=$r.sheet_fraction; loop_fraction=$r.loop_fraction; result_id=$r.result_id }
  [void]$out.Add($o)
}
$out | Where-Object { $_.campaign -eq 'bridge' } | Export-Csv "$HOME\Downloads\bridge_30-50_designs_v2.csv" -NoTypeInformation
$out | Where-Object { $_.campaign -eq 'deep' } | Export-Csv "$HOME\Downloads\deep_35-60_designs_v2.csv" -NoTypeInformation
$out | Export-Csv "$HOME\Downloads\glp1r_100_designs_v2.csv" -NoTypeInformation
$mm = ($out | Where-Object { [int]$_.length -ne [int]$_.length_expected }).Count
Write-Host ("total: " + $out.Count + "  bridge: " + ($out | Where-Object { $_.campaign -eq 'bridge' }).Count + "  deep: " + ($out | Where-Object { $_.campaign -eq 'deep' }).Count + "  length mismatches: " + $mm + "  download failures: " + $bad)
Write-Host "Wrote glp1r_100_designs_v2.csv, bridge_30-50_designs_v2.csv, deep_35-60_designs_v2.csv to Downloads"
