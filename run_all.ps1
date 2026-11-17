# =============================================================================
# run_all.ps1 -- reproduce the TRN -> APSIM sorghum tillering results end to end.
#
#   pwsh ./run_all.ps1                     # full run, both networks, includes APSIM
#   pwsh ./run_all.ps1 -Arm hand_curated   # one network only
#   pwsh ./run_all.ps1 -SkipApsim          # phases 00-06 only (no APSIM needed)
#   pwsh ./run_all.ps1 -Arm ?              # list the registered networks
#
# Works on Windows (PowerShell 5.1 or 7) and on macOS / Linux (PowerShell 7, `pwsh`).
# Everything is read from and written inside this folder.
#
# The networks ("arms") are defined once in scripts/arms.json and read from there
# by this script and by every Python script.
# =============================================================================
param(
    [string]$Arm = "all",
    [switch]$SkipApsim
)

# Continue (not Stop): R and APSIM write progress to stderr, which under 'Stop'
# PowerShell treats as fatal. We gate on $LASTEXITCODE instead.
$ErrorActionPreference = "Continue"
$ROOT = $PSScriptRoot

# --- the arm registry (single source of truth) ---
$Registry = Get-Content (Join-Path $ROOT "scripts/arms.json") -Raw | ConvertFrom-Json
$AllArms  = $Registry.arms
$ArmNames = $AllArms | ForEach-Object { $_.name }

if ($Arm -eq "?") {
    Write-Host "registered arms:"
    foreach ($a in $AllArms) { Write-Host ("  {0,-16}{1,-26}{2}" -f $a.name, $a.label, $a.network_dir) }
    exit 0
}
if ($Arm -ne "all" -and $ArmNames -notcontains $Arm) {
    throw "unknown arm '$Arm'. Choose 'all' or one of: $($ArmNames -join ', ')"
}

$Selected = if ($Arm -eq "all") { $AllArms } else { $AllArms | Where-Object name -eq $Arm }
$SelNames = @($Selected | ForEach-Object { $_.name })
$IsFull   = ($Arm -eq "all")

# --- locate the three toolchains (PATH first, then the usual install folders) ---
function Find-Tool($names, $globs, $what) {
    foreach ($n in $names) {
        $c = Get-Command $n -ErrorAction SilentlyContinue | Select-Object -First 1
        if ($c) { return $c.Source }
    }
    foreach ($g in $globs) {
        $c = Get-Item $g -ErrorAction SilentlyContinue |
             Sort-Object FullName -Descending | Select-Object -First 1
        if ($c) { return $c.FullName }
    }
    throw "$what not found. See README.md (Requirements)."
}

$Python  = Find-Tool @("python3", "python") @() "Python >= 3.10"
$Rscript = Find-Tool @("Rscript", "Rscript.exe") @(
               "C:\Program Files\R\*\bin\Rscript.exe",
               "/Library/Frameworks/R.framework/Resources/bin/Rscript",
               "/opt/homebrew/bin/Rscript", "/usr/local/bin/Rscript", "/usr/bin/Rscript"
           ) "Rscript (R >= 4.4)"
$Apsim   = $null
if (-not $SkipApsim) {
    $Apsim = Find-Tool @("Models", "Models.exe") @(
                 "C:\Program Files\APSIM*\bin\Models.exe",
                 "/Applications/APSIM*.app/Contents/Resources/bin/Models",
                 "/usr/local/lib/apsim/*/bin/Models"
             ) "APSIM Next Generation (Models / Models.exe)"
}

function Step($msg, [scriptblock]$block) {
    Write-Host "`n=== $msg ===" -ForegroundColor Cyan
    & $block
    if ($LASTEXITCODE -ne 0) { throw "step failed: $msg" }
}

Write-Host "root    $ROOT"
Write-Host "python  $Python"
Write-Host "Rscript $Rscript"
if ($Apsim) { Write-Host "APSIM   $Apsim" } else { Write-Host "APSIM   (skipped)" }
Write-Host "arm     $Arm  ->  $($SelNames -join ', ')"

# --- 0. environment check ---
# The sweep runs the committed model/model_baseline/nextStep.R, so base R is enough.
Step "00 environment check" {
    & $Python -c "import sys; assert sys.version_info >= (3, 10), 'Python >= 3.10 required, found ' + sys.version.split()[0]; import matplotlib, numpy, pandas; print('Python', sys.version.split()[0], '+ matplotlib/numpy/pandas OK')"
    if ($LASTEXITCODE -ne 0) { return }
    & $Rscript -e "cat('R', as.character(getRversion()), fill=TRUE)"
}
$pub = Join-Path $ROOT "00_reference/apsim_published/DensityManagement_NoTiller.db"
if (-not (Test-Path $pub)) {
    Write-Host "note: 00_reference/apsim_published/*.db not found -- the check against the" -ForegroundColor Yellow
    Write-Host "      published APSIM run (phase 08) will be skipped. See README.md to download it." -ForegroundColor Yellow
}

# --- 1. audit: gene coverage, and the Python cross-check of the curated model ---
Step "01 gene / genotype coverage"        { & $Python "$ROOT/scripts/audit_gene_coverage.py" }
Step "01b PSoup Python cross-check"       { & $Python "$ROOT/scripts/col1_reproduce_psoup.py" }

# --- 2. build the sucrose-extended Flash-P network from the frozen export ---
if ($Selected | Where-Object { $_.extend_sucrose }) {
    Step "02 extend Flash-P network with sucrose edges" { & $Python "$ROOT/scripts/extend_sucrose.py" }
}

# --- 3. sucrose sweeps: each network at SUC 0.5 / 1.0 / 1.5 ---
if ($Selected | Where-Object { $_.kind -eq "r_psoup" }) {
    Step "03a sucrose sweep - hand-curated (R)" { & $Rscript "$ROOT/scripts/run_psoup_sweep.R" }
}
$FlashPArms = @($Selected | Where-Object { $_.kind -eq "python_flashp" } | ForEach-Object { $_.name })
if ($FlashPArms) {
    Step "03b sucrose sweep - Flash-P (Python)"  { & $Python "$ROOT/scripts/run_sucrose_sweep.py" @FlashPArms }
}

# --- 4. how the solvers / networks compare on the same genotypes ---
if ($IsFull) {
    Step "04a pattern comparison" { & $Python "$ROOT/scripts/compare_patterns.py" }
}
if ($FlashPArms) {
    Step "04b solver comparison"  { & $Python "$ROOT/scripts/compare_solvers.py" @FlashPArms }
}

# --- 5. sucrose response -> APSIM tillerSd parameters, and the methods figures ---
Step "05a APSIM parameters"    { & $Python "$ROOT/scripts/build_apsim_parameters.py" }
Step "05b cultivar selection"  { & $Python "$ROOT/scripts/derive_apsim_params.py" }
if ($IsFull) {
    Step "05c methods figures"  { & $Python "$ROOT/scripts/make_2network_figures.py" --only methods }
}

# --- 6. stage a patched .apsimx per network from the pristine templates ---
# Arm-filtered: staging deletes the .db beside each .apsimx it writes.
Step "06 stage APSIM runs" { & $Python "$ROOT/scripts/build_apsim_runs.py" @SelNames }

# --- 7. run APSIM ---
if (-not $SkipApsim) {
    foreach ($a in $SelNames) {
        foreach ($f in @("Exploration_NoTiller.apsimx", "DensityManagement_NoTiller.apsimx")) {
            $p = Join-Path $ROOT "03_runs/apsim/$a/$f"
            # The templates predate recent APSIM releases; upgrade the staged COPY.
            Step "07 upgrade  $a/$f" { & $Apsim --upgrade "$p" }
            Step "07 APSIM    $a/$f" { & $Apsim "$p" }
        }
    }
}

$dbs = Get-ChildItem "$ROOT/03_runs/apsim" -Recurse -Filter *.db -ErrorAction SilentlyContinue |
       Where-Object { $_.Length -gt 0 }
if ($SkipApsim -and -not $dbs) {
    Write-Host "`nSTOPPING: -SkipApsim was given and step 06 cleared the .db files," -ForegroundColor Yellow
    Write-Host "so phases 08-09 have no APSIM output to read. Re-run without -SkipApsim." -ForegroundColor Yellow
    exit 0
}

# --- 8. compare the networks against each other and against the published run ---
if ($IsFull) {
    Step "08a APSIM arm comparison" { & $Python "$ROOT/scripts/compare_apsim_runs.py" }
    Step "08b crossover comparison" { & $Python "$ROOT/scripts/compare_crossover.py" }
}

# --- 9. result figures ---
if ($IsFull) {
    Step "09 result figures"  { & $Python "$ROOT/scripts/make_2network_figures.py" --only results }
    # --- 10. compare against the results shipped in expected_results/ ---
    Step "10 check against expected_results/" { & $Python "$ROOT/scripts/check_against_expected.py" }
} else {
    Step "09 figures 5-9 for $($SelNames -join ', ')" { & $Python "$ROOT/scripts/make_figures.py" @SelNames }
}

Write-Host "`nDONE. Outputs in 01_audit/, 03_runs/, 04_figures/." -ForegroundColor Green
