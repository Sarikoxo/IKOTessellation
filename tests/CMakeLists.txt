function(iko_test_target name source)
  add_executable(${name} ${source})
  set_target_properties(${name} PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
  target_link_libraries(${name} PRIVATE iko::tessellation)
  target_include_directories(${name} PRIVATE ${PROJECT_SOURCE_DIR}/third_party/triangle)
  set_source_files_properties(${source} PROPERTIES OBJECT_DEPENDS "${sdk_headers}")
  if(MSVC)
    target_compile_options(${name} PRIVATE /fp:strict /utf-8)
  else()
    target_compile_options(${name} PRIVATE -ffp-contract=off -fno-fast-math)
  endif()
endfunction()

# In-tree tests/probes reuse hidden production objects in an uninstalled helper.
# Public API calls still use the DLL, including probe api_model operations.
add_library(iko_tessellation_test_support STATIC $<TARGET_OBJECTS:iko_tessellation_core>)
target_link_libraries(iko_tessellation_test_support PRIVATE TKDEIGES)
if(MSVC)
  target_link_libraries(iko_tessellation_test_support PRIVATE bcrypt ws2_32)
else()
  target_link_libraries(iko_tessellation_test_support PRIVATE m)
endif()

if(IKO_BUILD_PROBE)
  iko_test_target(tessellation_probe probe.cpp)
  target_link_libraries(tessellation_probe PRIVATE iko_tessellation_test_support)
  set_target_properties(tessellation_probe PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
endif()

if(BUILD_TESTING)
  iko_test_target(tessellation_topology_import topology_import_test.cpp)
  add_test(NAME sdk_topology_import COMMAND tessellation_topology_import
    ${CMAKE_CURRENT_SOURCE_DIR}/data)
  set_tests_properties(sdk_topology_import PROPERTIES TIMEOUT 180)
  iko_test_target(tessellation_topology topology_test.cpp)
  target_link_libraries(tessellation_topology PRIVATE iko_tessellation_test_support)
  add_test(NAME sdk_topology COMMAND tessellation_topology)
  set_tests_properties(sdk_topology PROPERTIES TIMEOUT 180)
  iko_test_target(tessellation_contract contract_test.cpp)
  target_link_libraries(tessellation_contract PRIVATE iko_tessellation_test_support)
  add_test(NAME sdk_contract COMMAND tessellation_contract
    ${CMAKE_CURRENT_SOURCE_DIR}/data ${CMAKE_CURRENT_BINARY_DIR}/output/contract)
  set_tests_properties(sdk_contract PROPERTIES TIMEOUT 180)
  add_test(NAME sdk_package COMMAND ${CMAKE_COMMAND}
    "-DSDK_BUILD_DIR=${CMAKE_BINARY_DIR}"
    "-DSDK_SOURCE_DIR=${PROJECT_SOURCE_DIR}"
    "-DOUTPUT_DIR=${CMAKE_CURRENT_BINARY_DIR}/output/package"
    "-DLIBDIR=${CMAKE_INSTALL_LIBDIR}"
    "-DCXX_COMPILER=${CMAKE_CXX_COMPILER}"
    "-DC_COMPILER=${CMAKE_C_COMPILER}"
    "-DRC_COMPILER=${CMAKE_RC_COMPILER}"
    "-DMT=${CMAKE_MT}"
    "-DCONFIG=$<CONFIG>"
    "-DGENERATOR=${CMAKE_GENERATOR}"
    "-DMAKE_PROGRAM=${CMAKE_MAKE_PROGRAM}"
    -P ${CMAKE_CURRENT_SOURCE_DIR}/verify_package.cmake)
  set_tests_properties(sdk_package PROPERTIES TIMEOUT 60)
  if(MSVC)
    get_filename_component(_iko_compiler_bin "${CMAKE_CXX_COMPILER}" DIRECTORY)
    add_test(NAME sdk_dll COMMAND ${CMAKE_COMMAND}
      "-DDUMPBIN=${_iko_compiler_bin}/dumpbin.exe"
      "-DDLL=$<TARGET_FILE:iko_tessellation>"
      "-DIMPLIB=$<TARGET_LINKER_FILE:iko_tessellation>"
      -P ${CMAKE_CURRENT_SOURCE_DIR}/verify_dll.cmake)
  endif()
  iko_test_target(tessellation_exact exact_test.cpp)
  add_test(NAME sdk_exact COMMAND tessellation_exact)
  if(IKO_BUILD_CLI)
    iko_test_target(tessellation_cli_test cli_test.cpp)
    add_test(NAME sdk_cli COMMAND tessellation_cli_test
      $<TARGET_FILE:tessellate> ${CMAKE_CURRENT_SOURCE_DIR}/data
      ${CMAKE_CURRENT_BINARY_DIR}/output/cli $<TARGET_FILE:iko_tessellation>)
    set_tests_properties(sdk_cli PROPERTIES TIMEOUT 240)
  endif()
endif()
