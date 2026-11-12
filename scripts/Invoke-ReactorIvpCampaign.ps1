[CmdletBinding()]
param(
    [ValidateSet('all', 'correctness', 'performance')]
    [string]$Suite = 'all',
    [switch]$PlanOnly,
    [switch]$FailOnAny,
    [ValidateRange(1, 32)]
    [int]$BuildJobs = 1
)

$ErrorActionPreference = 'Stop'
$repoRoot = (Resolve-Path (Join-Path $PSScriptRoot '..')).Path
$timestamp = (Get-Date).ToUniversalTime().ToString('yyyyMMdd-HHmmss-fffZ')
$runId = "reactor-ivp-$Suite-release-$timestamp"
$reportsRoot = Join-Path $repoRoot 'test_reports'
$runRoot = Join-Path $reportsRoot $runId
$suffix = 1
while (Test-Path -LiteralPath $runRoot) {
    $runRoot = Join-Path $reportsRoot "$runId-$suffix"
    $suffix++
}
$reportsDir = Join-Path $runRoot 'reports'
$technicalDir = Join-Path $runRoot 'technical'

function New-Step([string]$Id, [string]$SuiteName, [string]$Filter, [string]$Description, [string]$EnvName = '', [string]$EnvValue = '') {
    [pscustomobject]@{
        Id = $Id
        Suite = $SuiteName
        Description = $Description
        Filter = $Filter
        EnvironmentName = $EnvName
        EnvironmentValue = $EnvValue
        Arguments = @('test', '--release', '--lib', '--no-default-features', "-j$BuildJobs", $Filter, '--', '--nocapture', '--test-threads=1')
    }
}

$steps = @()
if ($Suite -in @('all', 'correctness')) {
    $steps += New-Step 'ivp-stories' 'correctness' 'SimpleReactorsIVP::reactor_ivp_story_tests' 'IVP correctness, bounded physical smoke, frontend/family parity and prepared lifecycle stories'
    $steps += New-Step 'ivp-multistep-frontend-parity' 'correctness' 'multistep_oxidation_expr_and_atomview_publish_matching_bounded_profiles' 'ExprLegacy/AtomView parity on the bounded two-step oxidation symbolic system'
    $steps += New-Step 'ivp-aot-stories' 'correctness' 'SimpleReactorsIVP::reactor_ivp_aot_story_tests' 'Gated LSODE2 AOT sparse/banded correctness stories' 'KITHE_RUN_IVP_AOT_TESTS' '1'
}
if ($Suite -in @('all', 'performance')) {
    $steps += [pscustomobject]@{
        Id = 'ivp-release-benchmark'; Suite = 'performance'
        Description = 'Criterion wall-time matrix plus detailed KiThe stages, RST telemetry, prepared lifecycle and AOT/Lambdify axes'
        Filter = 'reactor_ivp_baseline'; EnvironmentName = 'KITHE_RUN_IVP_AOT_BENCH'; EnvironmentValue = '1'
        Arguments = @('bench', '--bench', 'reactor_ivp_baseline', '--no-default-features', "-j$BuildJobs", '--', '--noplot', '--sample-size', '30', '--warm-up-time', '1', '--measurement-time', '2')
    }
}

function Format-Command($Step) {
    $quoted = foreach ($argument in $Step.Arguments) {
        if ($argument -match '[\s"]') { '"' + ($argument -replace '"', '\"') + '"' } else { $argument }
    }
    "cargo $($quoted -join ' ')"
}

Write-Host "ReactorIVP release campaign: $Suite; build jobs=$BuildJobs"
Write-Host "Run ID: $runId"
if ($PlanOnly) {
    foreach ($step in $steps) {
        $environment = if ($step.EnvironmentName) { "$($step.EnvironmentName)=$($step.EnvironmentValue)" } else { 'no extra environment' }
        Write-Host ("[{0}] {1} | {2} | {3}" -f $step.Id, $step.Suite, $environment, (Format-Command $step))
    }
    Write-Host 'Plan only: Cargo was not invoked and no report directory was created.'
    exit 0
}

New-Item -ItemType Directory -Path $reportsDir, $technicalDir -Force | Out-Null
$revision = (& git -C $repoRoot rev-parse HEAD 2>$null | Select-Object -First 1)
if (-not $revision) { $revision = 'unavailable' }
$dirty = [bool](& git -C $repoRoot status --porcelain 2>$null)
$rustc = (& rustc --version 2>$null | Select-Object -First 1)
$cargoVersion = (& cargo --version 2>$null | Select-Object -First 1)
$tcc = Get-Command tcc -ErrorAction SilentlyContinue
$records = [System.Collections.Generic.List[object]]::new()
$campaignClock = [System.Diagnostics.Stopwatch]::StartNew()

foreach ($step in $steps) {
    $logPath = Join-Path $technicalDir ("{0}.log" -f $step.Id)
    $timer = [System.Diagnostics.Stopwatch]::StartNew()
    $started = (Get-Date).ToUniversalTime()
    $status = 'failed'
    $reason = ''
    $exitCode = $null

    if ($step.Id -eq 'ivp-aot-stories' -and -not $tcc) {
        $status = 'not applicable'
        $reason = 'AOT correctness requires tcc on PATH; step skipped, remaining campaign continues.'
        'AOT correctness skipped: tcc was not found on PATH.' | Set-Content -LiteralPath $logPath -Encoding utf8
    }
    else {
        Write-Host ("START {0} -> {1}" -f $step.Id, $logPath)
        $oldEnvironmentValue = $null
        if ($step.EnvironmentName) {
            $oldEnvironmentValue = [Environment]::GetEnvironmentVariable($step.EnvironmentName, 'Process')
        }
        try {
            if ($step.EnvironmentName) { [Environment]::SetEnvironmentVariable($step.EnvironmentName, $step.EnvironmentValue, 'Process') }
            Push-Location $repoRoot
            try {
                $previousErrorActionPreference = $ErrorActionPreference
                $ErrorActionPreference = 'Continue'
                & cargo @($step.Arguments) *> $logPath
                $exitCode = $LASTEXITCODE
                $ErrorActionPreference = $previousErrorActionPreference
            }
            finally {
                if ($null -ne $previousErrorActionPreference) { $ErrorActionPreference = $previousErrorActionPreference }
                Pop-Location
            }
            if ($exitCode -eq 0) {
                $status = 'passed'
                if ($step.Suite -eq 'correctness') {
                    $logText = Get-Content -LiteralPath $logPath -Raw
                    if ($logText -notmatch '(?m)^test result: ok\.\s+[1-9][0-9]* passed') {
                        $status = 'failed'
                        $reason = 'Cargo succeeded but the filter executed no passing tests.'
                    }
                }
            }
            else { $reason = "Cargo exited with code $exitCode." }
        }
        catch {
            $reason = $_.Exception.Message
            try { $_ | Out-String | Set-Content -LiteralPath $logPath -Encoding utf8 } catch { }
        }
        finally {
            if ($step.EnvironmentName) { [Environment]::SetEnvironmentVariable($step.EnvironmentName, $oldEnvironmentValue, 'Process') }
        }
    }

    $timer.Stop()
    $record = [pscustomobject]@{
        id = $step.Id; suite = $step.Suite; description = $step.Description
        command = (Format-Command $step); status = $status; reason = $reason
        started_utc = $started.ToString('o'); duration_seconds = [math]::Round($timer.Elapsed.TotalSeconds, 3)
        exit_code = $exitCode; technical_log = "technical/$($step.Id).log"
    }
    $records.Add($record)
    Write-Host ("END   {0}: {1} ({2}s, exit={3})" -f $step.Id, $status, $record.duration_seconds, $exitCode)
    if ($reason) { Write-Host "      $reason" }
    if ($status -ne 'passed' -and (Test-Path -LiteralPath $logPath)) {
        Get-Content -LiteralPath $logPath -Tail 8 | ForEach-Object { Write-Host "        $_" }
    }
}
$campaignClock.Stop()

$criterionRows = @()
$criterionRoot = Join-Path $repoRoot 'target/criterion/reactor_ivp_h2_dissociation'
$benchmarkRecord = $records | Where-Object { $_.id -eq 'ivp-release-benchmark' } | Select-Object -First 1
if ($benchmarkRecord -and $benchmarkRecord.status -eq 'passed' -and (Test-Path -LiteralPath $criterionRoot)) {
    foreach ($estimateFile in Get-ChildItem -LiteralPath $criterionRoot -Recurse -Filter estimates.json | Where-Object { $_.Directory.Name -eq 'new' }) {
        $routeDir = $estimateFile.Directory.Parent
        $fixtureDir = $routeDir.Parent
        if (-not $fixtureDir -or $fixtureDir.FullName -eq $criterionRoot) { continue }
        try {
            $estimate = Get-Content -LiteralPath $estimateFile.FullName -Raw | ConvertFrom-Json
            $median = $estimate.median
            $criterionRows += [pscustomobject]@{
                fixture = $fixtureDir.Name; route = $routeDir.Name
                median_ms = [math]::Round($median.point_estimate / 1e6, 4)
                ci_low_ms = [math]::Round($median.confidence_interval.lower_bound / 1e6, 4)
                ci_high_ms = [math]::Round($median.confidence_interval.upper_bound / 1e6, 4)
            }
        }
        catch { Write-Warning "Could not parse Criterion estimate: $($estimateFile.FullName)" }
    }
}

$diagnosticTables = ''
$benchLog = Join-Path $technicalDir 'ivp-release-benchmark.log'
if (Test-Path -LiteralPath $benchLog) {
    $lines = Get-Content -LiteralPath $benchLog
    $markers = @(
        'Macro: instrumented wall time',
        'KiThe stage detail:',
        'RST native telemetry:',
        'Callback policy axis:',
        'Prepared RST lifecycle:',
        'LSODE2 execution axis:'
    )
    foreach ($marker in $markers) {
        $start = -1
        for ($i = 0; $i -lt $lines.Count; $i++) { if ($lines[$i].StartsWith($marker)) { $start = $i; break } }
        if ($start -lt 0) { continue }
        $tableLines = [System.Collections.Generic.List[string]]::new()
        $seenTable = $false
        for ($i = $start; $i -lt $lines.Count; $i++) {
            if ($lines[$i] -match '^\s*\|') { $seenTable = $true; $tableLines.Add($lines[$i].TrimEnd()) }
            elseif ($seenTable) { break }
        }
        if ($tableLines.Count -gt 0) { $diagnosticTables += "`n## $marker`n`n$($tableLines -join "`n")`n" }
    }
}

$manifest = [ordered]@{
    run_id = $runId; suite = $Suite; profile = 'release'
    created_utc = (Get-Date).ToUniversalTime().ToString('o')
    duration_seconds = [math]::Round($campaignClock.Elapsed.TotalSeconds, 3)
    git_revision = $revision; working_tree_dirty = $dirty
    host = [System.Runtime.InteropServices.RuntimeInformation]::OSDescription
    rustc = $rustc; cargo = $cargoVersion; build_jobs = $BuildJobs; tcc_available = [bool]$tcc
    steps = @($records); criterion_rows = @($criterionRows); diagnostic_tables_markdown = $diagnosticTables
}
$manifestPath = Join-Path $reportsDir 'campaign.json'
$summaryPath = Join-Path $reportsDir 'summary.md'
$manifest | ConvertTo-Json -Depth 8 | Set-Content -LiteralPath $manifestPath -Encoding utf8

$rendererLog = Join-Path $technicalDir 'report-renderer.log'
Push-Location $repoRoot
try {
    $previousErrorActionPreference = $ErrorActionPreference
    $ErrorActionPreference = 'Continue'
    & cargo run --quiet --example reactor_ivp_campaign_report --no-default-features "-j$BuildJobs" -- $manifestPath $summaryPath *> $rendererLog
    $rendererExit = $LASTEXITCODE
}
finally {
    if ($null -ne $previousErrorActionPreference) { $ErrorActionPreference = $previousErrorActionPreference }
    Pop-Location
}
if ($rendererExit -ne 0 -or -not (Test-Path -LiteralPath $summaryPath)) {
    @('# ReactorIVP campaign summary', '', 'The tabled report renderer failed; see `technical/report-renderer.log`.', '', 'Raw command logs and campaign JSON are preserved.') |
        Set-Content -LiteralPath $summaryPath -Encoding utf8
    Write-Warning "Summary renderer failed (exit $rendererExit); campaign logs are preserved."
}

Write-Host "Summary: $summaryPath"
Write-Host "Manifest: $manifestPath"
Write-Host "Technical logs: $technicalDir"
$failures = @($records | Where-Object { $_.status -in @('failed', 'interrupted') })
if ($FailOnAny -and $failures.Count -gt 0) { exit 1 }
exit 0
