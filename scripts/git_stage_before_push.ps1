$ErrorActionPreference = 'Stop'

$repoRoot = git rev-parse --show-toplevel
if ($LASTEXITCODE -ne 0 -or [string]::IsNullOrWhiteSpace($repoRoot)) {
    throw 'Current directory is not inside a Git repository.'
}

Set-Location -LiteralPath $repoRoot

$paths = @(
    'Cargo.toml'
    'scripts/'
    'examples/'
    'benches/'
    'src/'
)

git add -- @paths
if ($LASTEXITCODE -ne 0) {
    throw 'git add failed.'
}

Write-Host 'Staged changes in the requested paths:'
git diff --cached --stat -- @paths
if ($LASTEXITCODE -ne 0) {
    throw 'Unable to read the staged diff summary.'
}
