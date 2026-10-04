cmake_minimum_required(VERSION 3.20)

set(package_prefix "${OUTPUT_DIR}/install")
set(module_dir "${package_prefix}/lib/cmake/iko_tessellation")

function(install_sdk)
  execute_process(COMMAND "${CMAKE_COMMAND}" --install "${SDK_BUILD_DIR}"
    --config "${CONFIG}" --prefix "${package_prefix}"
    RESULT_VARIABLE status OUTPUT_VARIABLE output ERROR_VARIABLE error)
  if(NOT status EQUAL 0)
    message(FATAL_ERROR "SDK installation failed:\n${output}\n${error}")
  endif()
endfunction()

function(configure_case name expected_error)
  execute_process(COMMAND "${CMAKE_COMMAND}"
    -S "${OUTPUT_DIR}/${name}" -B "${OUTPUT_DIR}/${name}-build"
    -G "${GENERATOR}" "-DCMAKE_MAKE_PROGRAM=${MAKE_PROGRAM}"
    "-DCMAKE_PREFIX_PATH=${package_prefix}"
    "-DCMAKE_CXX_COMPILER=${CXX_COMPILER}" "-DCMAKE_BUILD_TYPE=${CONFIG}"
    "-DCMAKE_RC_COMPILER=${RC_COMPILER}" "-DCMAKE_MT=${MT}"
    "-DCMAKE_TRY_COMPILE_PLATFORM_VARIABLES=CMAKE_RC_COMPILER;CMAKE_MT"
    -DCMAKE_DISABLE_FIND_PACKAGE_OpenCASCADE=ON ${ARGN}
    RESULT_VARIABLE status OUTPUT_VARIABLE output ERROR_VARIABLE error)
  if(expected_error STREQUAL "")
    if(NOT status EQUAL 0)
      message(FATAL_ERROR "Package case '${name}' failed:\n${output}\n${error}")
    endif()
  elseif(status EQUAL 0 OR NOT "${output}\n${error}" MATCHES "${expected_error}")
    message(FATAL_ERROR "Package case '${name}' did not reject the missing component:\n${output}\n${error}")
  endif()
endfunction()

install_sdk()
if(EXISTS "${package_prefix}/lib/cmake/iko/iko-config-version.cmake" OR
    EXISTS "${module_dir}/iko_tessellation-config-version.cmake")
  message(FATAL_ERROR "The installed iko SDK must not contain package version files")
endif()
file(READ "${SDK_SOURCE_DIR}/cmake/iko-config.cmake" common_source)
file(READ "${package_prefix}/lib/cmake/iko/iko-config.cmake" common_installed)
if(NOT common_source STREQUAL common_installed)
  message(FATAL_ERROR "The common iko entry point must be installed verbatim")
endif()

file(WRITE "${OUTPUT_DIR}/consumer/CMakeLists.txt" [=[
cmake_minimum_required(VERSION 3.20)
project(iko_package_check LANGUAGES CXX)
find_package(iko CONFIG REQUIRED COMPONENTS tessellation OPTIONAL_COMPONENTS unavailable)
if(NOT iko_tessellation_FOUND OR iko_unavailable_FOUND OR NOT TARGET iko::tessellation)
  message(FATAL_ERROR "Required and optional iko components were not resolved correctly")
endif()
get_target_property(sdk_include_dir iko::tessellation INTERFACE_INCLUDE_DIRECTORIES)
if(NOT EXISTS "${sdk_include_dir}/tessellation/tessellation.h")
  message(FATAL_ERROR "The installed module header is missing")
endif()
if(EXISTS "${sdk_include_dir}/tessellation/private")
  message(FATAL_ERROR "Private implementation headers must not be installed")
endif()
get_target_property(sdk_configs iko::tessellation IMPORTED_CONFIGURATIONS)
get_target_property(sdk_kind iko::tessellation TYPE)
get_target_property(sdk_dependencies iko::tessellation INTERFACE_LINK_LIBRARIES)
if(NOT sdk_kind STREQUAL "SHARED_LIBRARY" OR sdk_dependencies OR OpenCASCADE_FOUND OR TARGET TKDEIGES)
  message(FATAL_ERROR "The installed SDK must be a DLL with no OCCT link dependencies")
endif()
foreach(sdk_config IN LISTS sdk_configs)
  string(TOLOWER "${sdk_config}" sdk_config_dir)
  get_target_property(sdk_library iko::tessellation IMPORTED_LOCATION_${sdk_config})
  if(NOT EXISTS "${sdk_library}")
    message(FATAL_ERROR "The installed module library is missing: ${sdk_library}")
  endif()
  if(WIN32 AND NOT sdk_library MATCHES "/bin/${sdk_config_dir}/iko_tessellationd?\\.dll$")
    message(FATAL_ERROR "DLL must be installed under bin/debug or bin/release")
  endif()
  get_target_property(sdk_import_library iko::tessellation IMPORTED_IMPLIB_${sdk_config})
  if(WIN32 AND NOT EXISTS "${sdk_import_library}")
    message(FATAL_ERROR "The installed DLL import library is missing")
  endif()
  if(WIN32 AND NOT sdk_import_library MATCHES "/lib/${sdk_config_dir}/iko_tessellationd?\\.lib$")
    message(FATAL_ERROR "Import library must be installed under lib/debug or lib/release")
  endif()
endforeach()
find_package(iko QUIET CONFIG COMPONENTS unavailable)
if(iko_FOUND OR iko_unavailable_FOUND)
  message(FATAL_ERROR "An unavailable required component must reject the package")
endif()
set(iko_tessellation_FOUND FALSE)
find_package(iko CONFIG REQUIRED)
if(NOT iko_FOUND OR iko_tessellation_FOUND)
  message(FATAL_ERROR "An entry-only lookup must recover without loading a default module")
endif()
find_package(iko_tessellation CONFIG REQUIRED)
add_executable(consumer consumer.cpp)
target_link_libraries(consumer PRIVATE iko::tessellation)
add_custom_command(TARGET consumer POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E copy_if_different $<TARGET_FILE:iko::tessellation>
    $<TARGET_FILE_DIR:consumer>)
]=])
configure_file("${SDK_SOURCE_DIR}/tests/package_consumer.cpp"
  "${OUTPUT_DIR}/consumer/consumer.cpp" COPYONLY)

file(WRITE "${OUTPUT_DIR}/missing/CMakeLists.txt" [=[
cmake_minimum_required(VERSION 3.20)
project(iko_missing_component NONE)
find_package(iko CONFIG REQUIRED COMPONENTS unavailable)
]=])

# A separate prefix with just the common entry and an unrelated module proves
# that neither tessellation nor OCCT is required by the shared entry point.
set(entry_prefix "${OUTPUT_DIR}/entry-install")
file(MAKE_DIRECTORY "${entry_prefix}/lib/cmake/iko")
file(COPY "${SDK_SOURCE_DIR}/cmake/iko-config.cmake"
  DESTINATION "${entry_prefix}/lib/cmake/iko")
file(WRITE "${entry_prefix}/lib/cmake/iko_test_module/iko_test_module-config.cmake" [=[
set(iko_test_module_FOUND TRUE)
if(NOT TARGET iko::test_module)
  add_library(iko::test_module INTERFACE IMPORTED)
endif()
]=])
file(WRITE "${OUTPUT_DIR}/entry/CMakeLists.txt" [=[
cmake_minimum_required(VERSION 3.20)
project(iko_entry_check NONE)
find_package(iko CONFIG REQUIRED)
if(NOT iko_FOUND OR TARGET iko::tessellation OR OpenCASCADE_FOUND)
  message(FATAL_ERROR "The common entry point loaded a module or its dependencies")
endif()
find_package(iko CONFIG REQUIRED COMPONENTS test_module OPTIONAL_COMPONENTS tessellation unavailable)
if(NOT iko_test_module_FOUND OR NOT TARGET iko::test_module OR
    iko_tessellation_FOUND OR iko_unavailable_FOUND OR TARGET iko::tessellation)
  message(FATAL_ERROR "A separately installed module must work without a module list")
endif()
]=])

set(optional_source [=[
cmake_minimum_required(VERSION 3.20)
project(iko_optional_dependency NONE)
find_package(iko CONFIG REQUIRED OPTIONAL_COMPONENTS tessellation)
if(NOT iko_FOUND OR NOT iko_tessellation_FOUND OR NOT TARGET iko::tessellation OR OpenCASCADE_FOUND)
  message(FATAL_ERROR "An installed tessellation module must work without caller OCCT")
endif()
]=])
set(required_source [=[
cmake_minimum_required(VERSION 3.20)
project(iko_required_dependency NONE)
find_package(iko CONFIG REQUIRED COMPONENTS tessellation)
]=])
file(WRITE "${OUTPUT_DIR}/optional_dependency/CMakeLists.txt" "${optional_source}")
file(WRITE "${OUTPUT_DIR}/required_dependency/CMakeLists.txt" "${required_source}")
set(unpatched_occt "${OUTPUT_DIR}/unpatched-occt")
file(WRITE "${unpatched_occt}/OpenCASCADEConfig.cmake" [=[
set(OpenCASCADE_FOUND TRUE)
set(OpenCASCADE_BUILD_SHARED_LIBS FALSE)
set(OpenCASCADE_INCLUDE_DIR "${CMAKE_CURRENT_LIST_DIR}/include")
]=])
file(WRITE "${unpatched_occt}/OpenCASCADEConfigVersion.cmake" [=[
set(PACKAGE_VERSION "8.0.1")
set(PACKAGE_VERSION_COMPATIBLE TRUE)
set(PACKAGE_VERSION_EXACT TRUE)
]=])
file(WRITE "${unpatched_occt}/include/ProjLib_ComputeApproxOnPolarSurface.hxx"
  "// Fixture representing an unpatched OCCT header.\n")

configure_case(entry "" "-Diko_DIR=${entry_prefix}/lib/cmake/iko"
  -DCMAKE_DISABLE_FIND_PACKAGE_OpenCASCADE=ON)
configure_case(missing "Required iko component 'unavailable' was not found")
configure_case(optional_dependency "" -DCMAKE_DISABLE_FIND_PACKAGE_OpenCASCADE=ON)
configure_case(required_dependency ""
  -DCMAKE_DISABLE_FIND_PACKAGE_OpenCASCADE=ON)

# The patch remains a source-build requirement, never an installed dependency.
execute_process(COMMAND "${CMAKE_COMMAND}" -S "${SDK_SOURCE_DIR}"
  -B "${OUTPUT_DIR}/source-unpatched-build" -G "${GENERATOR}"
  "-DCMAKE_MAKE_PROGRAM=${MAKE_PROGRAM}" "-DCMAKE_C_COMPILER=${C_COMPILER}"
  "-DCMAKE_CXX_COMPILER=${CXX_COMPILER}" "-DCMAKE_RC_COMPILER=${RC_COMPILER}"
  "-DCMAKE_MT=${MT}" "-DCMAKE_BUILD_TYPE=${CONFIG}"
  "-DCMAKE_TRY_COMPILE_PLATFORM_VARIABLES=CMAKE_RC_COMPILER;CMAKE_MT"
  "-DOpenCASCADE_DIR=${unpatched_occt}"
  RESULT_VARIABLE status OUTPUT_VARIABLE output ERROR_VARIABLE error)
if(status EQUAL 0 OR NOT "${output}\n${error}" MATCHES "projection policy patch")
  message(FATAL_ERROR "Source builds must reject unpatched OCCT: ${output}\n${error}")
endif()

file(WRITE "${OUTPUT_DIR}/wrong_config/CMakeLists.txt" "${required_source}")
if(CONFIG STREQUAL "Debug")
  configure_case(wrong_config "configuration fallback is unsafe" -DCMAKE_BUILD_TYPE=Release)
else()
  configure_case(wrong_config "configuration fallback is unsafe" -DCMAKE_BUILD_TYPE=Debug)
endif()

if(WIN32)
  foreach(name IN ITEMS bad_crt bad_stl)
    file(WRITE "${OUTPUT_DIR}/${name}/CMakeLists.txt" [=[
cmake_minimum_required(VERSION 3.20)
project(iko_bad_abi LANGUAGES CXX)
find_package(iko CONFIG REQUIRED COMPONENTS tessellation)
add_library(check OBJECT check.cpp)
target_link_libraries(check PRIVATE iko::tessellation)
]=])
    file(WRITE "${OUTPUT_DIR}/${name}/check.cpp" "#include <tessellation/tessellation.h>\n")
  endforeach()
  configure_case(bad_crt "" "-DCMAKE_MSVC_RUNTIME_LIBRARY=MultiThreaded$<$<CONFIG:Debug>:Debug>")
  configure_case(bad_stl "" "-DCMAKE_CXX_FLAGS=/D_ITERATOR_DEBUG_LEVEL=1")
  foreach(name IN ITEMS bad_crt bad_stl)
    execute_process(COMMAND "${CMAKE_COMMAND}" --build "${OUTPUT_DIR}/${name}-build" --config "${CONFIG}"
      RESULT_VARIABLE status OUTPUT_VARIABLE output ERROR_VARIABLE error)
    if(name STREQUAL "bad_crt")
      set(expected "requires /MD")
    else()
      set(expected "requires matching Release/Debug and iterator debug level")
    endif()
    if(status EQUAL 0 OR NOT "${output}\n${error}" MATCHES "${expected}")
      message(FATAL_ERROR "Unsafe caller ABI was not rejected: ${output}\n${error}")
    endif()
  endforeach()
endif()

# Seed old, recorded files alongside files owned by another module/configuration
# and a user's unrecorded file. Reinstallation may remove only the obsolete one.
string(TOLOWER "${CONFIG}" config_lower)
set(record "${module_dir}/install-manifest-${config_lower}.txt")
if(NOT EXISTS "${record}")
  message(FATAL_ERROR "The module/configuration installation record is missing")
endif()
set(obsolete "include/tessellation/obsolete-recorded.h")
set(shared "include/tessellation/shared-recorded.h")
set(user_file "include/tessellation/user-note.txt")
set(other_module "include/test_module/other-owned.h")
set(other_export "lib/cmake/iko_tessellation/iko-tessellation-targets-retained.cmake")
set(obsolete_dll "bin/iko_tessellation-obsolete.dll")
if(CONFIG STREQUAL "Debug")
  set(other_library "${LIBDIR}/release/iko_tessellation.lib")
  set(other_dll "bin/release/iko_tessellation.dll")
  set(legacy_library "${LIBDIR}/iko_tessellationd.lib")
  set(legacy_dll "bin/iko_tessellationd.dll")
else()
  set(other_library "${LIBDIR}/debug/iko_tessellationd.lib")
  set(other_dll "bin/debug/iko_tessellationd.dll")
  set(legacy_library "${LIBDIR}/iko_tessellation.lib")
  set(legacy_dll "bin/iko_tessellation.dll")
endif()
foreach(path IN ITEMS "${legacy_library}" "${legacy_dll}" "bin/tessellate.exe")
  file(WRITE "${package_prefix}/${path}" "Recorded legacy flat-layout artifact.\n")
endforeach()
file(APPEND "${record}" "${legacy_library}\n${legacy_dll}\nbin/tessellate.exe\n")
foreach(path IN ITEMS "${obsolete}" "${shared}" "${user_file}" "${other_module}")
  file(WRITE "${package_prefix}/${path}" "Preserve unless owned and obsolete: ${path}\n")
endforeach()
file(WRITE "${package_prefix}/${obsolete_dll}" "Unowned DLL must be preserved.\n")
if(NOT EXISTS "${package_prefix}/${other_dll}")
  file(WRITE "${package_prefix}/${other_dll}" "Retained configuration DLL.\n")
endif()
file(WRITE "${package_prefix}/${other_export}" "# Retained configuration export.\n")
if(NOT EXISTS "${package_prefix}/${other_library}")
  file(WRITE "${package_prefix}/${other_library}" "Retained configuration library.\n")
endif()
file(WRITE "${module_dir}/install-manifest-retained.txt"
  "${shared}\n${other_export}\n${other_library}\n${other_dll}\nlib/cmake/iko/iko-config.cmake\n")
file(WRITE "${package_prefix}/lib/cmake/iko_test_module/install-manifest-release.txt"
  "${other_module}\nlib/cmake/iko/iko-config.cmake\n")
file(COPY "${entry_prefix}/lib/cmake/iko_test_module/iko_test_module-config.cmake"
  DESTINATION "${package_prefix}/lib/cmake/iko_test_module")
file(APPEND "${record}" "${obsolete}\n${shared}\n${other_export}\n${other_library}\n${other_dll}\n")
file(APPEND "${module_dir}/iko-tessellation-targets.cmake" "\n# Test changed common export.\n")
install_sdk()
if(EXISTS "${package_prefix}/${obsolete}")
  message(FATAL_ERROR "An obsolete file recorded for the current module/configuration was retained")
endif()
foreach(path IN ITEMS "${legacy_library}" "${legacy_dll}" "bin/tessellate.exe")
  if(EXISTS "${package_prefix}/${path}")
    message(FATAL_ERROR "Recorded obsolete flat-layout artifact was retained: ${path}")
  endif()
endforeach()
foreach(path IN ITEMS "${shared}" "${user_file}" "${other_module}" "${other_export}" "${other_library}" "${other_dll}" "${obsolete_dll}"
    "lib/cmake/iko_tessellation/install-manifest-retained.txt"
    "lib/cmake/iko_test_module/install-manifest-release.txt"
    "lib/cmake/iko_test_module/iko_test_module-config.cmake")
  if(NOT EXISTS "${package_prefix}/${path}")
    message(FATAL_ERROR "Installation removed a file belonging to another owner: ${path}")
  endif()
endforeach()
file(READ "${package_prefix}/lib/cmake/iko/iko-config.cmake" common_installed)
if(NOT common_source STREQUAL common_installed)
  message(FATAL_ERROR "Reinstallation changed or removed the shared entry point")
endif()
file(READ "${module_dir}/iko-tessellation-targets.cmake" module_export)
if(module_export MATCHES "Test changed common export")
  message(FATAL_ERROR "The current export was not updated during reinstallation")
endif()
configure_case(consumer "")
execute_process(COMMAND "${CMAKE_COMMAND}" --build "${OUTPUT_DIR}/consumer-build"
  --config "${CONFIG}" --verbose
  RESULT_VARIABLE status OUTPUT_VARIABLE output ERROR_VARIABLE error)
if(NOT status EQUAL 0 OR "${output}\n${error}" MATCHES "TK[A-Za-z0-9]+\\.lib")
  message(FATAL_ERROR "Installed consumer build failed or linked OCCT: ${output}\n${error}")
endif()
execute_process(COMMAND "${OUTPUT_DIR}/consumer-build/consumer${CMAKE_EXECUTABLE_SUFFIX}"
  "${SDK_SOURCE_DIR}/tests/data" "${OUTPUT_DIR}/consumer-output"
  RESULT_VARIABLE status OUTPUT_VARIABLE output ERROR_VARIABLE error)
if(NOT status EQUAL 0)
  message(FATAL_ERROR "Installed DLL consumer failed: ${status}\n${output}\n${error}")
endif()
message(STATUS "Installed DLL consumer, STL ownership and module/configuration installation ownership verified")
