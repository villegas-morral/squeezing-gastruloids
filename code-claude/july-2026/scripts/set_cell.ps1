# Replace the source of one notebook cell from a .jl file.
#
#   powershell -File scripts/set_cell.ps1 -CellId 4cd7bbef -SourceFile scripts/newcell_4cd7bbef.jl
#
# Writes a .bak first and verifies that no OTHER cell's source changed.

param(
  [Parameter(Mandatory=$true)][string]$CellId,
  [Parameter(Mandatory=$true)][string]$SourceFile,
  [string]$Notebook = "squeezing-gastruloids_claude.ipynb"
)
$ErrorActionPreference = "Stop"

# PowerShell 5.1 reads a BOM-less file with the ANSI codepage unless told
# otherwise, which mangles every non-ASCII character in the notebook (an em dash
# becomes "a€"", a subscript zero becomes "a,€"). Read and write UTF-8
# explicitly, and write WITHOUT a BOM -- Jupyter and git both expect none.
$raw = [System.IO.File]::ReadAllText((Resolve-Path $Notebook), [System.Text.Encoding]::UTF8)
$nb  = $raw | ConvertFrom-Json
$before = @{}
foreach ($c in $nb.cells) { $before[$c.id] = ($c.source -join '') }

$cell = $nb.cells | Where-Object { $_.id -eq $CellId }
if (-not $cell) { throw "cell $CellId not found" }

# ipynb `source` is an array of lines, each keeping its trailing newline except
# the last one.
$text  = [System.IO.File]::ReadAllText((Resolve-Path $SourceFile)) -replace "`r`n", "`n"
$lines = $text -split "`n"
if ($lines[-1] -eq '') { $lines = $lines[0..($lines.Count-2)] }
$arr = @()
for ($i = 0; $i -lt $lines.Count; $i++) {
  if ($i -lt $lines.Count - 1) { $arr += ($lines[$i] + "`n") } else { $arr += $lines[$i] }
}
$cell.source = $arr
if ($cell.PSObject.Properties.Name -contains 'outputs') { $cell.outputs = @() }
if ($cell.PSObject.Properties.Name -contains 'execution_count') { $cell.execution_count = $null }

Copy-Item $Notebook "$Notebook.bak" -Force
$json = $nb | ConvertTo-Json -Depth 100
$utf8NoBom = New-Object System.Text.UTF8Encoding($false)
[System.IO.File]::WriteAllText((Resolve-Path $Notebook), $json, $utf8NoBom)

# verify
$nb2 = [System.IO.File]::ReadAllText((Resolve-Path $Notebook), [System.Text.Encoding]::UTF8) | ConvertFrom-Json
if ($nb2.cells.Count -ne $nb.cells.Count) { throw "cell count changed!" }
$changed = @()
foreach ($c in $nb2.cells) {
  $now = ($c.source -join '')
  if ($before[$c.id] -ne $now) { $changed += $c.id }
}
Write-Output ("cells: {0}   changed: {1}" -f $nb2.cells.Count, ($changed -join ', '))
if ($changed.Count -ne 1 -or $changed[0] -ne $CellId) {
  throw "unexpected cells changed -- restore from $Notebook.bak"
}
Write-Output "OK"
