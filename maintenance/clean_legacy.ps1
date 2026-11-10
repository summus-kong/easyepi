<# Preview legacy artifacts by default; use -Apply only after reviewing the list. #>
param([switch]$Apply)
$ErrorActionPreference = 'Stop'
$projectRoot = [IO.Path]::GetFullPath((Join-Path $PSScriptRoot '..'))
if ((Split-Path $projectRoot -Leaf) -ne 'easyepiv1' -or -not (Test-Path -LiteralPath (Join-Path $projectRoot 'pyproject.toml'))) {
    throw 'Expected easyepiv1 project root'
}
$targets = @(Get-ChildItem -LiteralPath $projectRoot -Directory -Force | Where-Object {
    $_.Name -match '^(\.final-wheel.*|\.review-wheel|\.pytest.*|pytest_.*|build|easyepi\.egg-info|smoke_.*)$'
})
foreach ($relative in @('easyepi/python','easyepi/backends','easyepi/scripts','easyepi/.vscode',
                        'tests/_cli_tmp','tests/_cooc_tmp','tests/_pipe_tmp','tests/_pub_tmp','tests/_rep_tmp','tests/fig_check','data','scripts')) {
    $candidate = Join-Path $projectRoot $relative
    if (Test-Path -LiteralPath $candidate) { $targets += Get-Item -LiteralPath $candidate }
}
foreach ($folder in @('easyepi','tests','examples')) {
    $targets += Get-ChildItem -LiteralPath (Join-Path $projectRoot $folder) -Recurse -Directory -Force |
        Where-Object { $_.Name -eq '__pycache__' }
}
$targets = @($targets | Sort-Object FullName -Unique | Sort-Object { $_.FullName.Length } -Descending)
foreach ($target in $targets) {
    $absoluteTarget = [IO.Path]::GetFullPath($target.FullName)
    if (-not $absoluteTarget.StartsWith($projectRoot + [IO.Path]::DirectorySeparatorChar,[StringComparison]::OrdinalIgnoreCase)) {
        throw "Target escapes project: $absoluteTarget"
    }
    if ($target.Attributes -band [IO.FileAttributes]::ReparsePoint) { throw "Refusing linked target: $absoluteTarget" }
    Write-Output $absoluteTarget
    if ($Apply -and (Test-Path -LiteralPath $absoluteTarget)) {
        Remove-Item -LiteralPath $absoluteTarget -Recurse -Force
    }
}
if (-not $Apply) { Write-Output 'Preview only. Review targets, then run this script with -Apply.' }
