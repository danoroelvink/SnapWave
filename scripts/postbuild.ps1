[CmdletBinding()]
param(
    [Parameter(Mandatory = $true)]
    [ValidateNotNullOrEmpty()]
    [string]$TargetDir
)

$ErrorActionPreference = 'Stop'

function Resolve-FullPath([string]$p) {
    try { return (Resolve-Path -LiteralPath $p).Path }
    catch { return $p }
}

# Sanitize TargetDir coming from VS (strip quotes and whitespace)
$dst = $TargetDir.Trim() -replace '^[\"'']+|[\"'']+$',''

if ([string]::IsNullOrWhiteSpace($dst)) {
    throw "TargetDir is empty after sanitization."
}

if (!(Test-Path -LiteralPath $dst)) {
    throw "TargetDir does not exist: $dst"
}

# Normalize to full path
$dst = (Resolve-Path -LiteralPath $dst).Path

# Validate ONEAPI_ROOT
if (-not $env:ONEAPI_ROOT -or [string]::IsNullOrWhiteSpace($env:ONEAPI_ROOT)) {
    throw "ONEAPI_ROOT is not set. Please set it in your environment (or VS build environment)."
}

# oneAPI compiler root (stable indirection)
$compilerRoot = Join-Path $env:ONEAPI_ROOT 'compiler\latest\windows'
if (!(Test-Path -LiteralPath $compilerRoot)) {
    # fallback if your environment uses a slightly different layout
    $compilerRootAlt = Join-Path $env:ONEAPI_ROOT 'compiler\latest'
    if (Test-Path -LiteralPath $compilerRootAlt) {
        $compilerRoot = $compilerRootAlt
    } else {
        throw "Intel oneAPI compiler root not found under: `n- $compilerRoot`n- $compilerRootAlt`n(Is the Intel compiler component installed?)"
    }
}

$dllNames = @(
    'libifcoremd.dll',
    'libmmd.dll',
    'libiomp5md.dll',
    'svml_dispmd.dll'
)

# Prefer redist locations, fall back to bin locations.
$searchRoots = @(
    (Join-Path $compilerRoot 'redist'),
    (Join-Path $compilerRoot 'bin'),
    (Join-Path $compilerRoot 'bin\intel64')   # harmless if absent, helpful on some layouts
) | Where-Object { Test-Path -LiteralPath $_ }

if ($searchRoots.Count -eq 0) {
    throw "No searchable runtime folders found under compiler root: $compilerRoot"
}

Write-Host "postbuild.ps1: ONEAPI_ROOT=$($env:ONEAPI_ROOT)"
Write-Host "postbuild.ps1: compilerRoot=$compilerRoot"
Write-Host "postbuild.ps1: targetDir=$dst"
Write-Host "postbuild.ps1: searching roots:"
$searchRoots | ForEach-Object { Write-Host "  - $_" }

foreach ($name in $dllNames) {
    $hit = $null

    foreach ($sr in $searchRoots) {
        $hit = Get-ChildItem -Path $sr -Recurse -File -Filter $name -ErrorAction SilentlyContinue |
               Select-Object -First 1
        if ($hit) { break }
    }

    if (-not $hit) {
        $rootsText = ($searchRoots -join "`n  - ")
        throw "Required Intel runtime DLL not found: $name`nSearched:`n  - $rootsText"
    }

    Copy-Item -LiteralPath $hit.FullName -Destination $dst -Force
    Write-Host "postbuild.ps1: Copied $name from $($hit.FullName)"
}

Write-Host "postbuild.ps1: Completed successfully."
exit 0
