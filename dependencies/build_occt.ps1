param(
  [string]$SourceDir = (Join-Path $PSScriptRoot '../.runs/deps/occt-8.0.1'),
  [string]$BuildDir = (Join-Path $PSScriptRoot '../.runs/occt-8.0.1-build'),
  [string]$InstallDir = (Join-Path $PSScriptRoot '../.runs/occt-8.0.1-install'),
  [string]$VisualStudioPath,
  [ValidateSet('Release', 'Debug')][string]$BuildType = 'Release',
  [ValidateSet('Dynamic')][string]$Runtime = 'Dynamic',
  [int]$Parallel = 8
)
$ErrorActionPreference = 'Stop'
if (-not $PSBoundParameters.ContainsKey('BuildDir')) {
  $buildSuffix = if ($BuildType -eq 'Debug') { '-debug' } else { '' }
  if ($Runtime -eq 'Dynamic') { $buildSuffix += '-dynamic' }
  $BuildDir = Join-Path $PSScriptRoot "../.runs/occt-8.0.1-build$buildSuffix"
}
if (-not $PSBoundParameters.ContainsKey('InstallDir') -and $Runtime -eq 'Dynamic') {
  $InstallDir = Join-Path $PSScriptRoot '../.runs/occt-8.0.1-install-dynamic'
}
$runtimeLibrary = 'MultiThreaded$<$<CONFIG:Debug>:Debug>'
if ($Runtime -eq 'Dynamic') { $runtimeLibrary += 'DLL' }
$revision = 'b8f597c677811d1f9f4d8a97f5ae2825c0353a42' # official V8_0_1
if (-not (Test-Path -LiteralPath "$SourceDir/CMakeLists.txt")) {
  throw 'Fetch official Open-Cascade-SAS/OCCT V8_0_1 before building with build_occt.ps1'
}
$occtGitArgs = @('-c', "safe.directory=$([System.IO.Path]::GetFullPath($SourceDir).Replace('\','/'))", '-C', $SourceDir)
$actual = & git @occtGitArgs rev-parse HEAD
if ($LASTEXITCODE -ne 0 -or $actual -ne $revision) { throw 'OCCT source must be official V8_0_1 (b8f597c677811d1f9f4d8a97f5ae2825c0353a42)' }
$projectionPatch = Join-Path $PSScriptRoot 'occt_projection-8.0.1.patch'
& git @occtGitArgs apply --reverse --check $projectionPatch 2>$null
if ($LASTEXITCODE -ne 0) {
  & git @occtGitArgs apply --check $projectionPatch
  if ($LASTEXITCODE -ne 0) { throw 'The pinned OCCT source cannot accept the SDK projection compatibility patch' }
  & git @occtGitArgs apply $projectionPatch
  if ($LASTEXITCODE -ne 0) { throw 'OCCT projection compatibility patch failed' }
  # Rebuild these template instantiations even when localized MSVC include
  # diagnostics did not register the changed .gxx dependency in an old build.
  foreach ($file in @('Approx_FitAndDivide_0.cxx','Approx_FitAndDivide2d_0.cxx')) {
    [System.IO.File]::SetLastWriteTimeUtc((Join-Path $SourceDir "src/ModelingData/TKGeomBase/Approx/$file"), [DateTime]::UtcNow)
  }
}
$vswhere = Join-Path ${env:ProgramFiles(x86)} 'Microsoft Visual Studio/Installer/vswhere.exe'
$vsRoot = $VisualStudioPath
if (-not $vsRoot) {
  $vsRoot = & $vswhere -latest -products '*' -requires Microsoft.VisualStudio.Component.VC.Tools.x86.x64 -property installationPath
}
if (-not $vsRoot) { throw 'Visual Studio C++ build tools not found' }
$compilerVersion = (Get-Content -LiteralPath "$vsRoot/VC/Auxiliary/Build/Microsoft.VCToolsVersion.default.txt" -Raw).Trim()
$compilerRoot = "$vsRoot/VC/Tools/MSVC/$compilerVersion"
$sdkRoot = (Get-ItemProperty 'HKLM:/SOFTWARE/Microsoft/Windows Kits/Installed Roots').KitsRoot10
$sdkVersion = (Get-ChildItem -LiteralPath "$sdkRoot/Include" -Directory | Where-Object { $_.Name -match '^10\.\d+\.\d+\.\d+$' } | Sort-Object { [version]$_.Name } -Descending | Select-Object -First 1).Name
$cmake = "$vsRoot/Common7/IDE/CommonExtensions/Microsoft/CMake/CMake/bin/cmake.exe"
$ninja = "$vsRoot/Common7/IDE/CommonExtensions/Microsoft/CMake/Ninja/ninja.exe"
$savedInclude=$env:INCLUDE; $savedLib=$env:LIB; $savedPath=$env:PATH
try {
  $env:INCLUDE="$compilerRoot/include;$sdkRoot/Include/$sdkVersion/ucrt;$sdkRoot/Include/$sdkVersion/shared;$sdkRoot/Include/$sdkVersion/um"
  $env:LIB="$compilerRoot/lib/x64;$sdkRoot/Lib/$sdkVersion/ucrt/x64;$sdkRoot/Lib/$sdkVersion/um/x64"
  $env:PATH="$compilerRoot/bin/Hostx64/x64;$sdkRoot/bin/$sdkVersion/x64;$env:PATH"
  $configure = @('-S', $SourceDir, '-B', $BuildDir, '-G', 'Ninja', "-DCMAKE_MAKE_PROGRAM=$ninja", "-DCMAKE_C_COMPILER=$compilerRoot/bin/Hostx64/x64/cl.exe", "-DCMAKE_CXX_COMPILER=$compilerRoot/bin/Hostx64/x64/cl.exe", "-DCMAKE_BUILD_TYPE=$BuildType", '-DCMAKE_POLICY_DEFAULT_CMP0091=NEW', "-DCMAKE_MSVC_RUNTIME_LIBRARY=$runtimeLibrary", '-DBUILD_LIBRARY_TYPE=Static', '-DCMAKE_C_FLAGS_RELEASE=/O2 /Ob2 /DNDEBUG', '-DCMAKE_CXX_FLAGS_RELEASE=/O2 /Ob2 /DNDEBUG', '-DCMAKE_C_FLAGS_DEBUG=/Zi /Ob0 /Od /RTC1', '-DCMAKE_CXX_FLAGS_DEBUG=/Zi /Ob0 /Od /RTC1', '-DBUILD_ADDITIONAL_TOOLKITS=TKDEIGES;TKDESTEP', '-DBUILD_CPP_STANDARD=C++17', '-DBUILD_USE_VCPKG=OFF', '-DBUILD_GTEST=OFF', '-DBUILD_MODULE_FoundationClasses=OFF', '-DBUILD_MODULE_ModelingData=OFF', '-DBUILD_MODULE_ModelingAlgorithms=OFF', '-DBUILD_MODULE_Visualization=OFF', '-DBUILD_MODULE_ApplicationFramework=OFF', '-DBUILD_MODULE_DataExchange=OFF', '-DBUILD_MODULE_Draw=OFF', '-DBUILD_DOC_Overview=OFF', '-DUSE_TBB=OFF', '-DUSE_FREETYPE=OFF', '-DBUILD_USE_PCH=OFF', "-DINSTALL_DIR=$([System.IO.Path]::GetFullPath($InstallDir))")
  # Upstream synchronizes both cache entries; set both for the install exports.
  $installPrefix = [System.IO.Path]::GetFullPath($InstallDir).Replace('\','/')
  $configure += "-DCMAKE_INSTALL_PREFIX=$installPrefix"
  & $cmake @configure
  if ($LASTEXITCODE -ne 0) { throw 'OCCT configure failed' }
  & $cmake --build $BuildDir --config $BuildType --parallel $Parallel
  if ($LASTEXITCODE -ne 0) { throw 'OCCT build failed' }
  & $cmake --install $BuildDir --config $BuildType
  if ($LASTEXITCODE -ne 0) { throw 'OCCT install failed' }
} finally { $env:INCLUDE=$savedInclude; $env:LIB=$savedLib; $env:PATH=$savedPath }
