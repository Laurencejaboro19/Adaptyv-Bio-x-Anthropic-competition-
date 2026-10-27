# Block 6: read each generated .cif and print every chain's length, residue-number range and sequence (CA atoms, plain PowerShell).
# The target comes back renumbered 1..N; the designed binder is the other chain.
$tag  = "my_first_run"
$base = "$env:USERPROFILE\Downloads\tfr1_rfd4\out\$tag\inference_samples"
$aa = @{ALA='A';ARG='R';ASN='N';ASP='D';CYS='C';GLN='Q';GLU='E';GLY='G';HIS='H';ILE='I';LEU='L';LYS='K';MET='M';PHE='F';PRO='P';SER='S';THR='T';TRP='W';TYR='Y';VAL='V'}
Get-ChildItem $base -Recurse -Filter "*_generated.cif" | Sort-Object FullName | ForEach-Object {
  $lines = Get-Content $_.FullName
  $cols = @($lines | Where-Object { $_.StartsWith('_atom_site.') } | ForEach-Object { $_.Trim() })
  $iAtom = [array]::IndexOf($cols,'_atom_site.label_atom_id')
  $iComp = [array]::IndexOf($cols,'_atom_site.label_comp_id')
  $iCh   = [array]::IndexOf($cols,'_atom_site.auth_asym_id')
  $iNum  = [array]::IndexOf($cols,'_atom_site.auth_seq_id')
  $chains = [ordered]@{}
  foreach ($l in $lines) {
    if ($l -match '^ATOM\s') {
      $p = $l.Trim() -split '\s+'
      if ($p[$iAtom] -eq 'CA') {
        $c = $p[$iCh]
        if (-not $chains.Contains($c)) { $chains[$c] = @{seq=''; nums=@()} }
        $one = $aa[$p[$iComp]]; if (-not $one) { $one = 'X' }
        $chains[$c].seq += $one; $chains[$c].nums += [int]$p[$iNum]
      }
    }
  }
  "=== $($_.Name)"
  foreach ($c in $chains.Keys) {
    $n = $chains[$c].nums
    "  chain {0}: {1} residues (numbers {2}-{3})" -f $c, $n.Count, ($n | Measure-Object -Minimum).Minimum, ($n | Measure-Object -Maximum).Maximum
    "    " + $chains[$c].seq
  }
}
