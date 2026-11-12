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
$runId = "reactor-bvp-$Suite-release-$timestamp"
$runRoot = Join-Path (Join-Path $repoRoot 'test_reports') $runId
$suffix = 1
while (Test-Path -LiteralPath $runRoot) {
    $runRoot = Join-Path (Join-Path $repoRoot 'test_reports') "$runId-$suffix"
    $suffix++
}
$reportsDir = Join-Path $runRoot 'reports'
$technicalDir = Join-Path $runRoot 'technical'

function New-TestStep([string]$Id, [string]$Filter, [string]$Description, [string]$EnvName = '', [string]$EnvValue = '') {
    [pscustomobject]@{
        Id = $Id; Suite = 'correctness'; Description = $Description; Filter = $Filter
        EnvironmentName = $EnvName; EnvironmentValue = $EnvValue
        Arguments = @('test', '--release', '--lib', '--no-default-features', "-j$BuildJobs", $Filter, '--', '--nocapture', '--test-threads=1')
    }
}

$steps = @()
if ($Suite -in @('all', 'correctness')) {
    $steps += New-TestStep 'bvp-physical-stories' 'SimpleReactorsBVP::reactor_bvp_story_tests' 'Compact HMX physical endpoint, route lifecycle, parser and facade correctness stories'
    $steps += New-TestStep 'bvp-analytical-transport' 'bvp01_' 'BVP-01 analytical transport, tuned adaptive grids and CorrectionNewton physical reference'
    $steps += New-TestStep 'bvp-analytical-reaction-diffusion' 'bvp02' 'BVP-02 analytical reaction-diffusion oracle, grid strategies and CorrectionNewton'
    $steps += New-TestStep 'bvp-analytical-thermochemical' 'bvp03' 'BVP-03 thermochemical oracle, energy invariant, grid strategies and CorrectionNewton'
    $steps += New-TestStep 'bvp-analytical-arrhenius' 'bvp04' 'BVP-04 temperature-dependent Arrhenius reference and tuned grid strategies'
    $steps += New-TestStep 'bvp-literature-zz1962b' 'zz1962b' 'ZZ-1962-B independent branch, adaptive-grid and CorrectionNewton validation'
    $steps += New-TestStep 'bvp-symbolic-matrix' 'matrix_compact_hmx_' 'Physical HMX frontend/layout matrix and enlarged-grid correctness checks'
    $steps += New-TestStep 'bvp-sparse-callback-policies' 'story_large_bvp_callback_policies_match_and_report_dispatches' 'Largest analytic case: Sequential/Auto/Forced Parallel parity and callback dispatch telemetry for Damped and BVP_sci'
    $steps += New-TestStep 'bvp-kithe-preparation-parity' 'story_kithe_pre_rst_symbolic_policies_preserve_large_equation_system' '40-species/39-reaction KiThe setup parity across sequential, automatic, and parallel symbolic RHS assembly'
    $steps += New-TestStep 'bvp-frozen-naive-matrix' 'frozen_and_naive_routes_match_across_frontends_and_matrix_layouts' 'Frozen/Naive analytic profile correctness across ExprLegacy/AtomView and Dense/Sparse/Banded'
    $steps += New-TestStep 'bvp-sci-aot-stories' 'matrix_compact_hmx_sci_aot_matches_lambdify_across_frontends_and_layouts' 'Matched BVP_sci TCC AOT/Lambdify profile parity' 'KITHE_RUN_BVP_AOT_TESTS' '1'
    $steps += New-TestStep 'bvp-damped-aot-stories' 'matrix_compact_hmx_damped_aot_matches_lambdify_across_frontends_and_matrices' 'Matched Damped TCC AOT/Lambdify profile parity' 'KITHE_RUN_BVP_AOT_TESTS' '1'
}
if ($Suite -in @('all', 'performance')) {
    $steps += [pscustomobject]@{
        Id = 'bvp-release-benchmark'; Suite = 'performance'
        Description = 'Criterion release matrix plus detailed KiThe/API stages and RST BVP telemetry; optional matched AOT axis'
        Filter = 'reactor_bvp_baseline'; EnvironmentName = 'KITHE_RUN_BVP_AOT_BENCH'; EnvironmentValue = '1'
        Arguments = @('bench', '--bench', 'reactor_bvp_baseline', '--no-default-features', "-j$BuildJobs", '--', '--noplot', '--sample-size', '25', '--warm-up-time', '1', '--measurement-time', '2')
    }
}

function Format-Command($Step) {
    if ($Step.Suite -eq 'correctness') {
        if ($script:releaseTestBinary) {
            return '"{0}" {1} --nocapture --test-threads=1' -f $script:releaseTestBinary, $Step.Filter
        }
        return 'shared release lib-test executable (built once): {0} --nocapture --test-threads=1' -f $Step.Filter
    }
    $quoted = foreach ($argument in $Step.Arguments) {
        if ($argument -match '[\s"]') { '"' + ($argument -replace '"', '\"') + '"' } else { $argument }
    }
    "cargo $($quoted -join ' ')"
}

Write-Host "ReactorBVP release campaign: $Suite; build jobs=$BuildJobs"
Write-Host "Run ID: $runId"
if ($PlanOnly) {
    if ($Suite -in @('all', 'correctness')) {
        Write-Host ("[release-lib-test-build] cargo test --release --lib --no-default-features -j{0} --no-run" -f $BuildJobs)
    }
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
$releaseTestBinary = $null

if ($Suite -in @('all', 'correctness')) {
    $buildId = 'release-lib-test-build'
    $buildLog = Join-Path $technicalDir "$buildId.log"
    $buildTimer = [System.Diagnostics.Stopwatch]::StartNew()
    $buildStarted = (Get-Date).ToUniversalTime()
    $buildExitCode = $null
    $buildStatus = 'failed'
    $buildReason = ''
    $buildArguments = @('test', '--release', '--lib', '--no-default-features', "-j$BuildJobs", '--no-run')
    try {
        Push-Location $repoRoot
        try {
            $previousErrorActionPreference = $ErrorActionPreference
            $ErrorActionPreference = 'Continue'
            & cargo @buildArguments *> $buildLog
            $buildExitCode = $LASTEXITCODE
            $ErrorActionPreference = $previousErrorActionPreference
        }
        finally {
            if ($null -ne $previousErrorActionPreference) { $ErrorActionPreference = $previousErrorActionPreference }
            Pop-Location
        }
        if ($buildExitCode -eq 0) {
            $releaseTestBinary = Get-ChildItem -LiteralPath (Join-Path $repoRoot 'target/release/deps') -Filter 'KiThe-*.exe' |
                Sort-Object LastWriteTimeUtc -Descending | Select-Object -First 1 -ExpandProperty FullName
            if ($releaseTestBinary) { $buildStatus = 'passed' }
            else { $buildReason = 'Cargo succeeded but no KiThe release test executable was found.' }
        }
        else { $buildReason = "Cargo exited with code $buildExitCode." }
    }
    catch {
        $buildReason = $_.Exception.Message
        try { $_ | Out-String | Set-Content -LiteralPath $buildLog -Encoding utf8 } catch { }
    }
    $buildTimer.Stop()
    $records.Add([pscustomobject]@{
        id = $buildId; suite = 'correctness'; description = 'Compile the release lib test executable once for all correctness filters'
        command = "cargo $($buildArguments -join ' ')"; status = $buildStatus; reason = $buildReason
        started_utc = $buildStarted.ToString('o'); duration_seconds = [math]::Round($buildTimer.Elapsed.TotalSeconds, 3)
        exit_code = $buildExitCode; technical_log = "technical/$buildId.log"
    })
    Write-Host ("END   {0}: {1} ({2}s, exit={3})" -f $buildId, $buildStatus, [math]::Round($buildTimer.Elapsed.TotalSeconds, 3), $buildExitCode)
    if ($buildReason) { Write-Host "      $buildReason" }
}

foreach ($step in $steps) {
    $logPath = Join-Path $technicalDir ("{0}.log" -f $step.Id)
    $timer = [System.Diagnostics.Stopwatch]::StartNew()
    $started = (Get-Date).ToUniversalTime()
    $status = 'failed'; $reason = ''; $exitCode = $null

    if ($step.Suite -eq 'correctness' -and -not $releaseTestBinary) {
        $status = 'not run'
        $reason = 'The shared release test executable was not built successfully.'
        'Correctness step skipped because the shared release test build failed.' | Set-Content -LiteralPath $logPath -Encoding utf8
    }
    elseif ($step.EnvironmentName -and -not $tcc -and $step.Id -ne 'bvp-release-benchmark') {
        $status = 'not applicable'
        $reason = 'TCC is not on PATH; the optional AOT step was skipped and the rest of the campaign continues.'
        'AOT step skipped: tcc was not found on PATH.' | Set-Content -LiteralPath $logPath -Encoding utf8
    }
    else {
        Write-Host ("START {0} -> {1}" -f $step.Id, $logPath)
        $oldEnvironmentValue = $null
        if ($step.EnvironmentName) { $oldEnvironmentValue = [Environment]::GetEnvironmentVariable($step.EnvironmentName, 'Process') }
        try {
            if ($step.EnvironmentName -and ($tcc -or $step.Id -ne 'bvp-release-benchmark')) {
                [Environment]::SetEnvironmentVariable($step.EnvironmentName, $step.EnvironmentValue, 'Process')
            }
            elseif ($step.Id -eq 'bvp-release-benchmark') {
                [Environment]::SetEnvironmentVariable($step.EnvironmentName, $null, 'Process')
                if (-not $tcc) { $reason = 'TCC unavailable; running the full Lambdify benchmark while omitting the optional AOT axis.' }
            }
            Push-Location $repoRoot
            try {
                $previousErrorActionPreference = $ErrorActionPreference
                $ErrorActionPreference = 'Continue'
                if ($step.Suite -eq 'correctness') {
                    & $releaseTestBinary $step.Filter --nocapture --test-threads=1 *> $logPath
                }
                else {
                    & cargo @($step.Arguments) *> $logPath
                }
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
                        $status = 'failed'; $reason = 'Cargo succeeded but the filter executed no passing tests.'
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
    if ($status -eq 'failed' -and (Test-Path -LiteralPath $logPath)) { Get-Content -LiteralPath $logPath -Tail 8 | ForEach-Object { Write-Host "        $_" } }
}
$campaignClock.Stop()

$criterionRows = @()
$criterionRoot = Join-Path $repoRoot 'target/criterion'
$benchRecord = $records | Where-Object { $_.id -eq 'bvp-release-benchmark' -and $_.status -eq 'passed' } | Select-Object -First 1
if ($benchRecord -and (Test-Path -LiteralPath $criterionRoot)) {
    foreach ($estimateFile in Get-ChildItem -LiteralPath $criterionRoot -Recurse -Filter estimates.json | Where-Object { $_.Directory.Name -eq 'new' }) {
        try {
            $estimate = Get-Content -LiteralPath $estimateFile.FullName -Raw | ConvertFrom-Json
            $parts = $estimateFile.Directory.Parent.FullName.Substring($criterionRoot.Length).TrimStart([IO.Path]::DirectorySeparatorChar).Split([IO.Path]::DirectorySeparatorChar)
            if ($parts.Count -lt 2 -or $parts[0] -notin @('reactor_bvp_analytic_profiles', 'reactor_bvp_kithe_preparation')) { continue }
            $criterionRows += [pscustomobject]@{
                benchmark = $parts[0]
                fixture = if ($parts.Count -gt 1) { $parts[1] } else { 'unknown' }
                route = if ($parts.Count -gt 2) { $parts[2..($parts.Count - 1)] -join '/' } else { 'unknown' }
                median_ms = [math]::Round($estimate.median.point_estimate / 1e6, 4)
                ci_low_ms = [math]::Round($estimate.median.confidence_interval.lower_bound / 1e6, 4)
                ci_high_ms = [math]::Round($estimate.median.confidence_interval.upper_bound / 1e6, 4)
            }
        }
        catch { Write-Warning "Could not parse Criterion estimate: $($estimateFile.FullName)" }
    }
}

$diagnosticTables = ''
$benchLog = Join-Path $technicalDir 'bvp-release-benchmark.log'
if (Test-Path -LiteralPath $benchLog) {
    $lines = Get-Content -LiteralPath $benchLog
    foreach ($marker in @('Detailed route summary', 'Stage attribution:', 'Native BVP telemetry', 'KiThe-side preparation', 'Largest-case callback policy comparison', 'AOT-TCC versus Lambdify')) {
        $start = -1
        for ($i = 0; $i -lt $lines.Count; $i++) { if ($lines[$i].Contains($marker)) { $start = $i; break } }
        if ($start -lt 0) { continue }
        $tableLines = [System.Collections.Generic.List[string]]::new(); $seenTable = $false
        for ($i = $start; $i -lt $lines.Count; $i++) {
            if ($lines[$i] -match '^\s*\|') { $seenTable = $true; $tableLines.Add($lines[$i].TrimEnd()) }
            elseif ($seenTable) { break }
        }
        if ($tableLines.Count -gt 0) { $diagnosticTables += "`n### $marker`n`n$($tableLines -join "`n")`n" }
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
    $previousErrorActionPreference = $ErrorActionPreference; $ErrorActionPreference = 'Continue'
    & cargo run --quiet --example reactor_bvp_campaign_report --no-default-features "-j$BuildJobs" -- $manifestPath $summaryPath *> $rendererLog
    $rendererExit = $LASTEXITCODE
}
finally {
    if ($null -ne $previousErrorActionPreference) { $ErrorActionPreference = $previousErrorActionPreference }
    Pop-Location
}
if ($rendererExit -ne 0 -or -not (Test-Path -LiteralPath $summaryPath)) {
    @('# ReactorBVP campaign summary', '', 'The tabled report renderer failed; inspect `technical/report-renderer.log`.', '', 'Raw command logs and campaign JSON are retained.') | Set-Content -LiteralPath $summaryPath -Encoding utf8
    Write-Warning "Summary renderer failed (exit $rendererExit); campaign logs are preserved."
}

Write-Host "Summary: $summaryPath"
Write-Host "Manifest: $manifestPath"
Write-Host "Technical logs: $technicalDir"
$failures = @($records | Where-Object { $_.status -in @('failed', 'interrupted') })
if ($FailOnAny -and $failures.Count -gt 0) { exit 1 }
exit 0
