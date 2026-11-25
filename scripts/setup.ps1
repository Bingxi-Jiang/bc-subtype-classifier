param(
    [string]$Python = "python",
    [switch]$SkipData
)
$ErrorActionPreference = "Stop"
Set-Location (Split-Path -Parent $PSScriptRoot)
function Assert-Exit([string]$Step) {
    if ($LASTEXITCODE -ne 0) { throw "$Step failed (exit $LASTEXITCODE)" }
}
if (-not (Test-Path .venv\Scripts\python.exe)) {
    & $Python -m venv .venv
    Assert-Exit "Create environment"
}
$ProjectPython = Join-Path (Get-Location) '.venv\Scripts\python.exe'
& $ProjectPython -m pip install --index-url https://pypi.org/simple -e '.[dev]'
Assert-Exit "Install project"
& $ProjectPython -m pytest -q
Assert-Exit "Tests"
if (-not $SkipData) {
    & $ProjectPython -m bc_subtype.cli download
    Assert-Exit "Download"
    & $ProjectPython -m bc_subtype.cli prepare
    Assert-Exit "Prepare"
}
Write-Host "Ready. Run: .venv\Scripts\bc-subtype.exe baseline --config configs/smoke.json --output-dir results/smoke"
