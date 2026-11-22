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
$runId = "reaction-network-ivp-$Suite-release-$timestamp"
$runRoot = Join-Path (Join-Path $repoRoot 'test_reports') $runId
$suffix = 1
while (Test-Path -LiteralPath $runRoot) {
    $runRoot = Join-Path (Join-Path $repoRoot 'test_reports') "$runId-$suffix"
    $suffix++
}
$reportsDir = Join-Path $runRoot 'reports'
$technicalDir = Join-Path $runRoot 'technical'

$steps = @()
if ($Suite -in @('all', 'correctness')) {
    $steps += [pscustomobject]@{
        Id = 'public-example'; Suite = 'correctness'
        Description = 'Compile the public-facade reaction-network example'
        Arguments = @('check', '--release', '--example', 'reaction_network_ivp', '--no-default-features', "-j$BuildJobs")
    }
    $steps += [pscustomobject]@{
        Id = 'reaction-network-stories'; Suite = 'correctness'
        Description = 'All reaction_network_ivp correctness and lifecycle tests'
        Filter = 'Kinetics::reaction_network_ivp'
    }
}
if ($Suite -in @('all', 'performance')) {
    $steps += [pscustomobject]@{
        Id = 'reaction-network-release-benchmark'; Suite = 'performance'
        Description = 'Cold/prepared Criterion baseline and separate Detailed stage-attribution matrix'
        Arguments = @('bench', '--bench', 'reaction_network_ivp_baseline', '--no-default-features', "-j$BuildJobs", '--', '--noplot', '--sample-size', '20', '--warm-up-time', '1', '--measurement-time', '2')
    }
}

function Format-Command($Step) {
    if ($Step.Suite -eq 'correctness') {
        if ($Step.Filter) { return "release lib-test executable $($Step.Filter) --nocapture --test-threads=1" }
    }
    $quoted = foreach ($argument in $Step.Arguments) {
        if ($argument -match '[\s"]') { '"' + ($argument -replace '"', '\"') + '"' } else { $argument }
    }
    "cargo $($quoted -join ' ')"
}

Write-Host "Reaction-network IVP release campaign: $Suite; build jobs=$BuildJobs"
Write-Host "Run ID: $runId"
if ($PlanOnly) {
    if ($Suite -in @('all', 'correctness')) {
        Write-Host "[release-lib-test-build] cargo test --release --lib --no-default-features -j$BuildJobs --no-run"
    }
    foreach ($step in $steps) {
        Write-Host ("[{0}] {1} | {2}" -f $step.Id, $step.Suite, (Format-Command $step))
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
$releaseTestBinary = $null

if ($Suite -in @('all', 'correctness')) {
    $buildId = 'release-lib-test-build'
    $buildLog = Join-Path $technicalDir "$buildId.log"
    $buildTimer = [System.Diagnostics.Stopwatch]::StartNew()
    $started = (Get-Date).ToUniversalTime()
    $buildArguments = @('test', '--release', '--lib', '--no-default-features', "-j$BuildJobs", '--no-run')
    $buildExit = $null
    $buildReason = ''
    try {
        Push-Location $repoRoot
        try {
            $oldPreference = $ErrorActionPreference
            $ErrorActionPreference = 'Continue'
            & cargo @buildArguments *> $buildLog
            $buildExit = $LASTEXITCODE
            $ErrorActionPreference = $oldPreference
        }
        finally {
            if ($null -ne $oldPreference) { $ErrorActionPreference = $oldPreference }
            Pop-Location
        }
        if ($buildExit -eq 0) {
            $releaseTestBinary = Get-ChildItem -LiteralPath (Join-Path $repoRoot 'target/release/deps') -Filter 'KiThe-*.exe' |
                Sort-Object LastWriteTimeUtc -Descending | Select-Object -First 1 -ExpandProperty FullName
            if (-not $releaseTestBinary) { $buildReason = 'Cargo succeeded but the release lib-test executable was not found.' }
        }
        else { $buildReason = "Cargo exited with code $buildExit." }
    }
    catch {
        $buildReason = $_.Exception.Message
        $_ | Out-String | Set-Content -LiteralPath $buildLog -Encoding utf8
    }
    $buildTimer.Stop()
    $buildStatus = if ($buildExit -eq 0 -and $releaseTestBinary) { 'passed' } else { 'failed' }
    $records.Add([pscustomobject]@{
        id = $buildId; suite = 'correctness'; description = 'Compile release lib tests once'
        command = "cargo $($buildArguments -join ' ')"; status = $buildStatus; reason = $buildReason
        started_utc = $started.ToString('o'); duration_seconds = [math]::Round($buildTimer.Elapsed.TotalSeconds, 3)
        exit_code = $buildExit; technical_log = "technical/$buildId.log"
    })
    Write-Host ("END   {0}: {1} ({2}s, exit={3})" -f $buildId, $buildStatus, [math]::Round($buildTimer.Elapsed.TotalSeconds, 3), $buildExit)
}

foreach ($step in $steps) {
    $logPath = Join-Path $technicalDir "$($step.Id).log"
    $timer = [System.Diagnostics.Stopwatch]::StartNew()
    $started = (Get-Date).ToUniversalTime()
    $status = 'failed'; $reason = ''; $exitCode = $null
    if ($step.Suite -eq 'correctness' -and $step.Filter -and -not $releaseTestBinary) {
        $status = 'not run'
        $reason = 'The shared release test executable was not built successfully.'
        $reason | Set-Content -LiteralPath $logPath -Encoding utf8
    }
    else {
        Write-Host "START $($step.Id) -> $logPath"
        try {
            Push-Location $repoRoot
            try {
                $oldPreference = $ErrorActionPreference
                $ErrorActionPreference = 'Continue'
                if ($step.Filter) {
                    & $releaseTestBinary $step.Filter --nocapture --test-threads=1 *> $logPath
                }
                else {
                    & cargo @($step.Arguments) *> $logPath
                }
                $exitCode = $LASTEXITCODE
                $ErrorActionPreference = $oldPreference
            }
            finally {
                if ($null -ne $oldPreference) { $ErrorActionPreference = $oldPreference }
                Pop-Location
            }
            if ($exitCode -eq 0) {
                $status = 'passed'
                if ($step.Filter -and (Get-Content -LiteralPath $logPath -Raw) -notmatch '(?m)^test result: ok\.\s+[1-9][0-9]* passed') {
                    $status = 'failed'
                    $reason = 'The filter ran no passing tests.'
                }
            }
            else { $reason = "Command exited with code $exitCode." }
        }
        catch {
            $reason = $_.Exception.Message
            $_ | Out-String | Set-Content -LiteralPath $logPath -Encoding utf8
        }
    }
    $timer.Stop()
    $records.Add([pscustomobject]@{
        id = $step.Id; suite = $step.Suite; description = $step.Description
        command = (Format-Command $step); status = $status; reason = $reason
        started_utc = $started.ToString('o'); duration_seconds = [math]::Round($timer.Elapsed.TotalSeconds, 3)
        exit_code = $exitCode; technical_log = "technical/$($step.Id).log"
    })
    Write-Host ("END   {0}: {1} ({2}s, exit={3})" -f $step.Id, $status, $records[$records.Count - 1].duration_seconds, $exitCode)
    if ($reason) { Write-Host "      $reason" }
    if ($status -eq 'failed' -and (Test-Path -LiteralPath $logPath)) {
        Get-Content -LiteralPath $logPath -Tail 8 | ForEach-Object { Write-Host "        $_" }
    }
}
$campaignClock.Stop()

$criterionRows = @()
$stageTable = @()
$configurationStageTable = @()
$criterionRoot = Join-Path $repoRoot 'target/criterion/reaction_network_ivp'
$benchRecord = $records | Where-Object { $_.id -eq 'reaction-network-release-benchmark' -and $_.status -eq 'passed' } | Select-Object -First 1
if ($benchRecord -and (Test-Path -LiteralPath $criterionRoot)) {
    foreach ($estimateFile in Get-ChildItem -LiteralPath $criterionRoot -Recurse -Filter estimates.json | Where-Object { $_.Directory.Name -eq 'new' }) {
        try {
            $estimate = Get-Content -LiteralPath $estimateFile.FullName -Raw | ConvertFrom-Json
            $criterionRows += [pscustomobject]@{
                benchmark = $estimateFile.Directory.Parent.FullName.Substring($criterionRoot.Length).TrimStart([IO.Path]::DirectorySeparatorChar)
                median_ms = [math]::Round($estimate.median.point_estimate / 1e6, 4)
                ci_low_ms = [math]::Round($estimate.median.confidence_interval.lower_bound / 1e6, 4)
                ci_high_ms = [math]::Round($estimate.median.confidence_interval.upper_bound / 1e6, 4)
            }
        }
        catch { Write-Warning "Could not parse Criterion estimate: $($estimateFile.FullName)" }
    }
}

$benchmarkLog = Join-Path $technicalDir 'reaction-network-release-benchmark.log'
if ($benchRecord -and (Test-Path -LiteralPath $benchmarkLog)) {
    $lines = Get-Content -LiteralPath $benchmarkLog
    for ($index = 0; $index -lt $lines.Count; $index++) {
        if (-not $lines[$index].StartsWith('Reaction-network IVP ') -or $lines[$index] -notmatch 'stages \(') { continue }
        $table = @()
        for ($tableIndex = $index + 1; $tableIndex -lt $lines.Count; $tableIndex++) {
            if ($lines[$tableIndex] -match '^\s*\|') { $table += $lines[$tableIndex].TrimEnd() }
            elseif ($table.Count -gt 0) { break }
        }
        if ($lines[$index] -match '40-state matrix/layout/callback') { $configurationStageTable = $table }
        else { $stageTable = $table }
    }
}

$manifest = [ordered]@{
    run_id = $runId; suite = $Suite; profile = 'release'
    created_utc = (Get-Date).ToUniversalTime().ToString('o')
    duration_seconds = [math]::Round($campaignClock.Elapsed.TotalSeconds, 3)
    git_revision = $revision; working_tree_dirty = $dirty
    host = [System.Runtime.InteropServices.RuntimeInformation]::OSDescription
    rustc = $rustc; cargo = $cargoVersion; build_jobs = $BuildJobs
    steps = @($records); criterion_rows = @($criterionRows)
}
$manifestPath = Join-Path $reportsDir 'campaign.json'
$manifest | ConvertTo-Json -Depth 8 | Set-Content -LiteralPath $manifestPath -Encoding utf8

$summary = [System.Collections.Generic.List[string]]::new()
$summary.Add('# Reaction-Network IVP Release Campaign')
$summary.Add('')
$summary.Add("Run: ``$runId`` | Profile: release | Duration: $($manifest.duration_seconds) s")
$summary.Add('')
$summary.Add('| Step | Suite | Status | Seconds | Exit |')
$summary.Add('|---|---|---:|---:|---:|')
foreach ($record in $records) {
    $summary.Add("| $($record.id) | $($record.suite) | $($record.status) | $($record.duration_seconds) | $($record.exit_code) |")
}
if ($criterionRows.Count -gt 0) {
    $summary.Add('')
    $summary.Add('## Criterion')
    $summary.Add('')
    $summary.Add('| Benchmark | Median ms | CI low ms | CI high ms |')
    $summary.Add('|---|---:|---:|---:|')
    foreach ($row in $criterionRows) {
        $summary.Add("| $($row.benchmark) | $($row.median_ms) | $($row.ci_low_ms) | $($row.ci_high_ms) |")
    }
}
if ($stageTable.Count -gt 0) {
    $summary.Add('')
    $summary.Add('## Stage Attribution')
    $summary.Add('')
    foreach ($line in $stageTable) { $summary.Add($line) }
}
if ($configurationStageTable.Count -gt 0) {
    $summary.Add('')
    $summary.Add('## 40-State Matrix And Callback Policy Stages')
    $summary.Add('')
    foreach ($line in $configurationStageTable) { $summary.Add($line) }
}
$summary.Add('')
$summary.Add('Raw command output is preserved in the per-step logs under `technical/`.')
$summaryPath = Join-Path $reportsDir 'summary.md'
$summary | Set-Content -LiteralPath $summaryPath -Encoding utf8

Write-Host "Summary: $summaryPath"
Write-Host "Manifest: $manifestPath"
Write-Host "Technical logs: $technicalDir"
$failures = @($records | Where-Object { $_.status -in @('failed', 'interrupted') })
if ($FailOnAny -and $failures.Count -gt 0) { exit 1 }
exit 0
