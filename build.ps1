param(
  [string]$BuildDir = (Join-Path $PSScriptRoot '.runs/build'),
  [string]$SourceDir = $PSScriptRoot,
  [string[]]$PrefixPath,
  [string]$InstallRoot,
  [string]$VisualStudioPath,
  [ValidateSet('Release', 'Debug')][string]$BuildType = 'Release',
  [ValidateSet('Dynamic')][string]$Runtime = 'Dynamic',
  [switch]$NoTests,
  [switch]$BuildProbe,
  [switch]$RunTests,
  [switch]$Clean
)
$ErrorActionPreference = 'Stop'
if ($RunTests -and $NoTests) { throw '-RunTests cannot be combined with -NoTests' }
if (-not $PSBoundParameters.ContainsKey('BuildDir')) {
  $buildSuffix = if ($BuildType -eq 'Debug') { '-debug' } else { '' }
  $BuildDir = Join-Path $PSScriptRoot ".runs/build$buildSuffix"
}
$runtimeLibrary = 'MultiThreaded$<$<CONFIG:Debug>:Debug>DLL'
if ($InstallRoot) { $InstallRoot = [System.IO.Path]::GetFullPath($InstallRoot) }
if ($Clean -and (Test-Path -LiteralPath "$BuildDir/CMakeCache.txt")) {
  $cache = Get-Content -LiteralPath "$BuildDir/CMakeCache.txt"
  $previousConfig = $cache | Where-Object { $_ -match '^CMAKE_BUILD_TYPE:STRING=' } | Select-Object -First 1
  $previousSource = $cache | Where-Object { $_ -match '^CMAKE_HOME_DIRECTORY:INTERNAL=' } | Select-Object -First 1
  if ($previousConfig -and $previousConfig.Split('=', 2)[1] -ne $BuildType) {
    throw 'Cleaning another build configuration is not allowed; use its own -BuildDir'
  }
  if ($previousSource -and [System.IO.Path]::GetFullPath($previousSource.Split('=', 2)[1]) -ne [System.IO.Path]::GetFullPath($SourceDir)) {
    throw 'Cleaning a build directory owned by another project is not allowed'
  }
}
$vswhere = Join-Path ${env:ProgramFiles(x86)} 'Microsoft Visual Studio/Installer/vswhere.exe'
$vsRoot = $VisualStudioPath
if (-not $vsRoot -and (Test-Path -LiteralPath $vswhere)) {
  $vsRoot = & $vswhere -latest -products '*' -requires Microsoft.VisualStudio.Component.VC.Tools.x86.x64 -property installationPath
}
if (-not $vsRoot) { throw 'Visual Studio C++ build tools not found; pass -VisualStudioPath or use CMake in an x64 developer shell' }
$compilerVersion = (Get-Content -LiteralPath "$vsRoot/VC/Auxiliary/Build/Microsoft.VCToolsVersion.default.txt" -Raw).Trim()
$compilerRoot = "$vsRoot/VC/Tools/MSVC/$compilerVersion"
$sdkRoot = (Get-ItemProperty 'HKLM:/SOFTWARE/Microsoft/Windows Kits/Installed Roots').KitsRoot10
$sdkVersion = (Get-ChildItem -LiteralPath "$sdkRoot/Include" -Directory | Where-Object { $_.Name -match '^10\.\d+\.\d+\.\d+$' } | Sort-Object { [version]$_.Name } -Descending | Select-Object -First 1).Name
if (-not $sdkVersion) { throw 'Windows 10/11 C++ SDK not found' }
$cmake = "$vsRoot/Common7/IDE/CommonExtensions/Microsoft/CMake/CMake/bin/cmake.exe"
$ninja = "$vsRoot/Common7/IDE/CommonExtensions/Microsoft/CMake/Ninja/ninja.exe"
$savedInclude=$env:INCLUDE; $savedLib=$env:LIB; $savedPath=$env:PATH
try {
  $env:INCLUDE="$compilerRoot/include;$sdkRoot/Include/$sdkVersion/ucrt;$sdkRoot/Include/$sdkVersion/shared;$sdkRoot/Include/$sdkVersion/um"
  $env:LIB="$compilerRoot/lib/x64;$sdkRoot/Lib/$sdkVersion/ucrt/x64;$sdkRoot/Lib/$sdkVersion/um/x64"
  $env:PATH="$compilerRoot/bin/Hostx64/x64;$sdkRoot/bin/$sdkVersion/x64;$env:PATH"
  if ($Runtime -eq 'Dynamic' -and $BuildType -eq 'Debug') {
    $debugCrtPath = "$vsRoot/VC/Redist/MSVC/$compilerVersion/debug_nonredist/x64"
    if (Test-Path -LiteralPath $debugCrtPath) {
      $debugCrt = Get-ChildItem -LiteralPath $debugCrtPath -Directory -Filter 'Microsoft.VC*.DebugCRT' | Select-Object -First 1
      if ($debugCrt) { $env:PATH="$($debugCrt.FullName);$env:PATH" }
    }
  }
  $configure = @('-S', $SourceDir, '-B', $BuildDir, '-G', 'Ninja', "-DCMAKE_MAKE_PROGRAM=$ninja", "-DCMAKE_C_COMPILER=$compilerRoot/bin/Hostx64/x64/cl.exe", "-DCMAKE_CXX_COMPILER=$compilerRoot/bin/Hostx64/x64/cl.exe", "-DCMAKE_BUILD_TYPE=$BuildType", "-DCMAKE_MSVC_RUNTIME_LIBRARY=$runtimeLibrary")
  if ($InstallRoot) { $configure += "-DIKO_INSTALL_ROOT=$InstallRoot" }
  if ($PrefixPath) {
    $resolvedPrefixes = ($PrefixPath | ForEach-Object { [System.IO.Path]::GetFullPath($_) }) -join ';'
    $configure += "-DCMAKE_PREFIX_PATH=$resolvedPrefixes"
  }
  $configure += "-DBUILD_TESTING=$(if ($NoTests) { 'OFF' } else { 'ON' })"
  $configure += "-DIKO_BUILD_PROBE=$(if ($BuildProbe -or -not $NoTests) { 'ON' } else { 'OFF' })"
  & $cmake @configure
  if ($LASTEXITCODE -ne 0) {throw 'Configure failed'}
  if ($Clean) { & $cmake --build $BuildDir --config $BuildType --clean-first } else { & $cmake --build $BuildDir --config $BuildType }
  if ($LASTEXITCODE -ne 0) {throw 'Build failed'}
  if ($RunTests) {
    $ctest = Join-Path (Split-Path -Parent $cmake) 'ctest.exe'
    & $ctest --test-dir $BuildDir -C $BuildType --output-on-failure
    if ($LASTEXITCODE -ne 0) { throw 'Tests failed' }
  }
  if ($InstallRoot) {
    & $cmake --install $BuildDir --config $BuildType --prefix $InstallRoot
    if ($LASTEXITCODE -ne 0) { throw 'Install failed' }
    Write-Host "iko installed to $InstallRoot"
  }
} finally { $env:INCLUDE=$savedInclude; $env:LIB=$savedLib; $env:PATH=$savedPath }
