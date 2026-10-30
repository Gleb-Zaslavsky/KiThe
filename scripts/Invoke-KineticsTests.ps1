[CmdletBinding()]
param(
    [ValidateRange(1, 32)]
    [int]$BuildJobs = 1
)

$ErrorActionPreference = 'Stop'
$repoRoot = (Resolve-Path (Join-Path $PSScriptRoot '..')).Path

Push-Location $repoRoot
try {
    & cargo test --lib --no-default-features "-j$BuildJobs" 'Kinetics::' -- `
        --include-ignored --nocapture --test-threads=1 --skip experimental_kinetics
    $testExitCode = $LASTEXITCODE
}
finally {
    Pop-Location
}

exit $testExitCode
