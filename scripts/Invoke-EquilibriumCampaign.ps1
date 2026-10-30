[CmdletBinding()]
param(
    [ValidateSet('smoke', 'correctness', 'performance', 'all')]
    [string]$Suite = 'all',
    [switch]$PlanOnly,
    [switch]$FailOnAny,
    [ValidateRange(1, 32)]
    [int]$BuildJobs = 1
)

$ErrorActionPreference = 'Stop'
$repoRoot = (Resolve-Path (Join-Path $PSScriptRoot '..')).Path
$timestamp = (Get-Date).ToUniversalTime().ToString('yyyyMMdd-HHmmss-fffZ')
$profileLabel = if ($Suite -eq 'performance') { 'bench' } elseif ($Suite -eq 'all') { 'debug+bench' } else { 'debug' }
$baseRunId = "equilibrium-$Suite-$profileLabel-$timestamp"
$reportsRoot = Join-Path $repoRoot 'test_reports'
$runId = $baseRunId
$runRoot = Join-Path $reportsRoot $runId
$suffix = 1
while (Test-Path -LiteralPath $runRoot) {
    $runId = "$baseRunId-$suffix"
    $runRoot = Join-Path $reportsRoot $runId
    $suffix++
}
$reportsDir = Join-Path $runRoot 'reports'
$technicalDir = Join-Path $runRoot 'technical'

function New-TestStep([string]$Id, [string]$Filter, [string]$Description, [bool]$IncludeIgnored = $false) {
    $arguments = @('test', '--lib', '--no-default-features', "-j$BuildJobs", $Filter, '--', '--nocapture', '--test-threads=1')
    if ($IncludeIgnored) { $arguments += '--include-ignored' }
    [pscustomobject]@{
        Id = $Id; Suite = 'correctness'; Profile = 'debug'; Description = $Description
        FilePath = 'cargo'
        Arguments = $arguments
        Workload = $Filter; Telemetry = 'test-defined'
    }
}

function New-BenchStep([string]$Name, [string]$Features, [string]$Description) {
    $arguments = @('bench', '--bench', $Name, '--no-default-features', "-j$BuildJobs")
    if ($Features) { $arguments += @('--features', $Features) }
    $arguments += @('--', '--noplot')
    [pscustomobject]@{
        Id = "bench-$Name"; Suite = 'performance'; Profile = 'bench'; Description = $Description
        FilePath = 'cargo'; Arguments = $arguments
        Workload = $Name; Telemetry = $Features
    }
}

$correctnessSteps = @(
    (New-TestStep 'tp1907-backend-matrix' 'tp1907_upper_inventory_is_backend_invariant' 'TP-1907 upper-inventory backend acceptance matrix'),
    (New-TestStep 'tp1907-frontend-parity' 'tp1907_upper_inventory_facade_atom_view_and_expr_legacy_match_published_equilibrium' 'TP-1907 AtomView/ExprLegacy facade parity'),
    (New-TestStep 'stanjan15-small-constant' 'stanjan15_tiny_constant_and_public_ph_facade_are_frontend_equivalent' 'STANJAN-15 small-constant P,H frontend parity'),
    (New-TestStep 'elemental-ph-small-scale' 'elemental_ph_small_scale_backend_matrix' 'Elemental P,H small-scale backend characterization'),
    (New-TestStep 'tp1906-condensed-ph-range' 'tp1906_four_point_ph_lifecycle_preserves_pt_manifold_and_graphite_topology' 'TP-1906 condensed/mixed-phase P,H lifecycle'),
    (New-TestStep 'gui-ph-range-order' 'ph_range_solves_through_the_gui_request_facade_in_requested_order' 'GUI P,H range order and facade route'),
    (New-TestStep 'gui-ph-range-failure' 'failing_intermediate_ph_range_point_returns_typed_index_without_partial_outcome' 'GUI P,H range typed failure behavior'),
    (New-TestStep 'gui-ph-range-cancel' 'fixed_ph_range_cancellation_discards_late_worker_completion' 'GUI P,H range cancellation'),
    (New-TestStep 'nist-shomate-offline' 'local_nist_water_shomate_fixture_matches_the_published_500k_anchor_offline' 'Offline local NIST/Shomate fixture'),
    (New-TestStep 'thermo-content-manifest' 'bundled_thermochemistry_content_manifest_matches_checked_in_data' 'Thermochemistry content-hash manifest'),
    (New-TestStep 'rst-partial-telemetry' 'failed_rst_callback_keeps_partial_callback_telemetry' 'Failed RST callback partial telemetry'),
    (New-TestStep 'ph-range-wide-accuracy' 'i5_wide_ph_range_temperature_accuracy_characterization' 'Wide P,H range temperature accuracy characterization' $true)
)

$performanceSteps = @(
    (New-BenchStep 'equilibrium_baseline' 'equilibrium-instrumentation equilibrium-benchmark-prettyprint' 'Facade baseline with stage attribution'),
    (New-BenchStep 'equilibrium_rst_preparation' 'equilibrium-instrumentation' 'RST preparation and parameterized solve baseline'),
    (New-BenchStep 'equilibrium_prepared_session' 'equilibrium-instrumentation' 'Prepared-session reuse across targets'),
    (New-BenchStep 'equilibrium_ph_preparation' 'equilibrium-instrumentation' 'P,H trial preparation and workflow stages'),
    (New-BenchStep 'equilibrium_gibbs_interval' 'equilibrium-instrumentation' 'Gibbs interval reuse comparison'),
    (New-BenchStep 'equilibrium_native_callbacks' 'equilibrium-instrumentation' 'Native residual/Jacobian callback costs'),
    (New-BenchStep 'equilibrium_native_overhead' 'equilibrium-instrumentation' 'Native solver callback/orchestration overhead'),
    (New-BenchStep 'equilibrium_rst_workers' 'equilibrium-instrumentation' 'RST sequential/parallel worker crossover'),
    (New-BenchStep 'equilibrium_rst_policy_crossover' '' 'RST solver policy crossover')
)

$steps = switch ($Suite) {
    'smoke' { @($correctnessSteps | Where-Object { $_.Id -in @('thermo-content-manifest', 'nist-shomate-offline', 'rst-partial-telemetry') }) }
    'correctness' { $correctnessSteps }
    'performance' { $performanceSteps }
    'all' { @($correctnessSteps) + @($performanceSteps) }
}

function Format-Command($Step) {
    $quoted = foreach ($argument in $Step.Arguments) {
        if ($argument -match '[\s"]') { '"' + ($argument -replace '"', '\"') + '"' } else { $argument }
    }
    "$($Step.FilePath) $($quoted -join ' ')"
}

Write-Host "Equilibrium campaign: $Suite ($profileLabel), build jobs=$BuildJobs"
Write-Host "Run ID: $runId"
if ($PlanOnly) {
    foreach ($step in $steps) {
        Write-Host ("[{0}] {1} | {2}" -f $step.Id, $step.Profile, (Format-Command $step))
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
$records = [System.Collections.Generic.List[object]]::new()
$campaignClock = [System.Diagnostics.Stopwatch]::StartNew()

foreach ($step in $steps) {
    $logPath = Join-Path $technicalDir ("{0}.log" -f $step.Id)
    $started = (Get-Date).ToUniversalTime()
    $timer = [System.Diagnostics.Stopwatch]::StartNew()
    $status = 'failed'
    $reason = ''
    $exitCode = $null
    Write-Host ("START {0} [{1}] -> {2}" -f $step.Id, $step.Profile, $logPath)
    try {
        Push-Location $repoRoot
        try {
            $previousErrorActionPreference = $ErrorActionPreference
            $ErrorActionPreference = 'Continue'
            & $step.FilePath @($step.Arguments) *> $logPath
            $exitCode = $LASTEXITCODE
            $ErrorActionPreference = $previousErrorActionPreference
        }
        finally {
            if ($null -ne $previousErrorActionPreference) { $ErrorActionPreference = $previousErrorActionPreference }
            Pop-Location
        }
        if ($exitCode -eq 0) {
            $status = 'passed'
        }
        elseif ($exitCode -in @(130, 143, -1073741510)) {
            $status = 'interrupted'
            $reason = "Cargo process interrupted (exit $exitCode); later steps will continue."
        }
        elseif (Test-Path -LiteralPath $logPath) {
            $logText = Get-Content -LiteralPath $logPath -Raw
            if ($logText -match '(?i)unexpected argument|feature .* is not defined|does not contain this feature|no (test|bench|example) target named') {
                $status = 'not applicable'
                $reason = 'Command or feature is unsupported by this checkout; see technical log.'
            }
            else { $reason = "Cargo exited with code $exitCode." }
        }
        else { $reason = "Cargo exited with code $exitCode without creating a log." }
        if ($status -eq 'passed' -and $step.Suite -eq 'correctness') {
            $logText = Get-Content -LiteralPath $logPath -Raw
            if ($logText -notmatch '(?m)^test result: ok\.\s+[1-9][0-9]* passed') {
                $status = 'failed'
                $reason = 'Cargo succeeded but the filter did not execute any passing test; see technical log.'
            }
        }
    }
    catch {
        $reason = $_.Exception.Message
        try { $_ | Out-String | Set-Content -LiteralPath $logPath -Encoding utf8 } catch { }
    }
    $timer.Stop()
    $ended = (Get-Date).ToUniversalTime()
    $record = [pscustomobject]@{
        id = $step.Id; suite = $step.Suite; profile = $step.Profile
        description = $step.Description; workload = $step.Workload; telemetry = $step.Telemetry
        command = (Format-Command $step); status = $status; reason = $reason
        started_utc = $started.ToString('o'); ended_utc = $ended.ToString('o')
        duration_seconds = [math]::Round($timer.Elapsed.TotalSeconds, 3)
        exit_code = $exitCode; report_path = $null; technical_log = ("technical/{0}.log" -f $step.Id)
    }
    $records.Add($record)
    Write-Host ("END   {0}: {1} ({2}s, exit={3})" -f $step.Id, $status, $record.duration_seconds, $exitCode)
    if ($reason) { Write-Host "      $reason" }
    if ($status -ne 'passed' -and (Test-Path -LiteralPath $logPath)) {
        Write-Host '      Last log lines:'
        Get-Content -LiteralPath $logPath -Tail 8 | ForEach-Object { Write-Host "        $_" }
    }
}
$campaignClock.Stop()

$manifest = [ordered]@{
    run_id = $runId; suite = $Suite; profile = $profileLabel
    created_utc = (Get-Date).ToUniversalTime().ToString('o')
    duration_seconds = [math]::Round($campaignClock.Elapsed.TotalSeconds, 3)
    git_revision = $revision; working_tree_dirty = $dirty
    host = [System.Runtime.InteropServices.RuntimeInformation]::OSDescription
    rustc = $rustc; cargo = $cargoVersion; build_jobs = $BuildJobs
    command_policy = 'one Cargo process at a time; continue after step failures'
    steps = @($records)
}
$manifestPath = Join-Path $reportsDir 'campaign.json'
$summaryPath = Join-Path $reportsDir 'summary.md'
$manifest | ConvertTo-Json -Depth 8 | Set-Content -LiteralPath $manifestPath -Encoding utf8

$rendererLog = Join-Path $technicalDir 'report-renderer.log'
Push-Location $repoRoot
try {
    $previousErrorActionPreference = $ErrorActionPreference
    $ErrorActionPreference = 'Continue'
    & cargo run --quiet --example equilibrium_campaign_report --no-default-features "-j$BuildJobs" -- $manifestPath $summaryPath *> $rendererLog
    $rendererExit = $LASTEXITCODE
}
finally {
    if ($null -ne $previousErrorActionPreference) { $ErrorActionPreference = $previousErrorActionPreference }
    Pop-Location
}
if ($rendererExit -ne 0 -or -not (Test-Path -LiteralPath $summaryPath)) {
    @('# Equilibrium campaign summary', '', 'The `tabled` report renderer failed; see `technical/report-renderer.log`.', '', ('Manifest: `{0}`' -f 'campaign.json')) |
        Set-Content -LiteralPath $summaryPath -Encoding utf8
    Write-Warning "Summary renderer failed (exit $rendererExit); JSON manifest and raw logs are preserved."
}
else { Write-Host "Summary: $summaryPath" }
Write-Host "Manifest: $manifestPath"
Write-Host "Technical logs: $technicalDir"

$failures = @($records | Where-Object { $_.status -in @('failed', 'interrupted', 'not applicable') })
if ($FailOnAny -and $failures.Count -gt 0) { exit 1 }
exit 0
