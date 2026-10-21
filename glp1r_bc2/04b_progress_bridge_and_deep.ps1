# Block 4b: short progress summary for the bridge and deep campaigns (read-only; works with the PC off).
cd "$env:USERPROFILE\Downloads\glp1r_modal"
foreach ($tag in "glp1r_bridge", "glp1r_deep") {
  $f = "$env:USERPROFILE\Downloads\progress_$tag.json"
  if (Test-Path $f) { Remove-Item $f -Force }
  python -m modal volume get pdfa-bc2-results "/$tag/progress.json" $f --force 2>$null | Out-Null
  "===== $tag"
  if ((Test-Path $f) -and (Get-Item $f).Length -gt 0) { Get-Content $f -Raw | ConvertFrom-Json | Format-List state, round, shards_launched, sequences_seen, refolded, good, goal, est_spend_usd, budget_cap_usd, stop_reason, elapsed_h, updated } else { "no progress file yet" }
}
