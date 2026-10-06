# Restore specific cells from squeezing-gastruloids.ipynb.bak and drop stray
# cells, preserving everything else in the current notebook.
#
#   powershell -File scripts/repair_notebook.ps1
#
# Written after a Jupyter autosave clobbered two cells with a stale in-memory
# state (one was split in two and lost ~1260 characters). Reads and writes UTF-8
# explicitly and BOM-less; PowerShell 5.1's defaults mangle non-ASCII.

param(
  [string]$Notebook = "squeezing-gastruloids_claude.ipynb",
  [string[]]$RestoreIds = @('9f98a77f', '0a275ce7'),
  [string[]]$DropIds    = @('a3630dfa')
)
$ErrorActionPreference = "Stop"
$U8 = [System.Text.Encoding]::UTF8

$cur = [System.IO.File]::ReadAllText((Resolve-Path $Notebook), $U8) | ConvertFrom-Json
$bak = [System.IO.File]::ReadAllText((Resolve-Path "$Notebook.bak"), $U8) | ConvertFrom-Json

foreach ($id in $RestoreIds) {
  $b = $bak.cells | Where-Object { $_.id -eq $id }
  $c = $cur.cells | Where-Object { $_.id -eq $id }
  if (-not $b) { throw "cell $id not in backup" }
  if (-not $c) { throw "cell $id not in current" }
  $c.source = $b.source
  Write-Output ("restored {0}: {1} -> {2} chars" -f $id, ($c.source -join '').Length, ($b.source -join '').Length)
}

if ($DropIds.Count -gt 0) {
  $before = $cur.cells.Count
  $cur.cells = @($cur.cells | Where-Object { $_.id -notin $DropIds })
  Write-Output ("dropped {0} stray cell(s): {1} -> {2}" -f ($before - $cur.cells.Count), $before, $cur.cells.Count)
}

$json = $cur | ConvertTo-Json -Depth 100
$noBom = New-Object System.Text.UTF8Encoding($false)
[System.IO.File]::WriteAllText((Resolve-Path $Notebook), $json, $noBom)

$chk = [System.IO.File]::ReadAllText((Resolve-Path $Notebook), $U8) | ConvertFrom-Json
Write-Output ("cells now: {0}" -f $chk.cells.Count)
foreach ($id in $RestoreIds) {
  $c = $chk.cells | Where-Object { $_.id -eq $id }
  Write-Output ("  {0}: {1} chars" -f $id, ($c.source -join '').Length)
}
