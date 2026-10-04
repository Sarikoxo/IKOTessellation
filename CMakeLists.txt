cmake_minimum_required(VERSION 3.20)
project(iko LANGUAGES C CXX)

set(IKO_INSTALL_ROOT "" CACHE PATH "Final shared installation directory for iko modules")
if(IKO_INSTALL_ROOT)
  get_filename_component(IKO_INSTALL_ROOT "${IKO_INSTALL_ROOT}" ABSOLUTE)
  set(CMAKE_INSTALL_PREFIX "${IKO_INSTALL_ROOT}" CACHE PATH "Installation prefix" FORCE)
endif()
include(GNUInstallDirs)
find_package(OpenCASCADE 8.0.1 EXACT CONFIG REQUIRED COMPONENTS
  FoundationClasses ModelingData ModelingAlgorithms Visualization ApplicationFramework DataExchange)
if(OpenCASCADE_BUILD_SHARED_LIBS)
  message(FATAL_ERROR "iko requires a static OCCT 8.0.1 installation")
endif()
include(cmake/check_occt_projection_policy.cmake)
iko_check_occt_projection_policy(_iko_occt_policy)
if(NOT _iko_occt_policy)
  message(FATAL_ERROR "OCCT 8.0.1 needs the SDK projection policy patch; build it with dependencies/build_occt.ps1")
endif()
unset(_iko_occt_policy)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(NOT DEFINED CMAKE_MSVC_RUNTIME_LIBRARY)
  set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>DLL")
endif()
if(MSVC AND NOT CMAKE_MSVC_RUNTIME_LIBRARY MATCHES "DLL$")
  message(FATAL_ERROR "iko DLL value types require the dynamic CRT: /MD Release, /MDd Debug")
endif()
option(IKO_BUILD_CLI "Build the tessellate command-line tool" ON)
option(IKO_BUILD_EXAMPLES "Build the iko example" ON)
include(CTest)
option(IKO_BUILD_PROBE "Build the machine-readable SDK test executable" ${BUILD_TESTING})

add_library(iko_tessellation_core OBJECT
  src/geometry.cpp src/parameterization.cpp src/mesh_generation.cpp
  src/content_hash.cpp src/options.cpp
  src/triangle.cpp src/trim.cpp src/pole.cpp src/trim_import.cpp
  src/bezier.cpp src/intersections.cpp src/repair.cpp
  src/combined_repair.cpp src/periodic_surface.cpp
  src/occt_import.cpp src/occt_adapt.cpp src/triangle_backend.c)
add_library(iko_tessellation SHARED src/mesh_export.cpp src/tessellation.cpp
  $<TARGET_OBJECTS:iko_tessellation_core>)
add_library(iko::tessellation ALIAS iko_tessellation)
set_target_properties(iko_tessellation PROPERTIES EXPORT_NAME tessellation DEBUG_POSTFIX d
  WINDOWS_EXPORT_ALL_SYMBOLS OFF CXX_VISIBILITY_PRESET hidden VISIBILITY_INLINES_HIDDEN ON)
set_target_properties(iko_tessellation_core PROPERTIES POSITION_INDEPENDENT_CODE ON
  C_VISIBILITY_PRESET hidden CXX_VISIBILITY_PRESET hidden VISIBILITY_INLINES_HIDDEN ON)
target_compile_definitions(iko_tessellation PRIVATE IKO_TESSELLATION_BUILD)
if(MSVC)
  target_compile_definitions(iko_tessellation PUBLIC
    IKO_TESSELLATION_MSC_VER=${MSVC_VERSION}
    IKO_TESSELLATION_ITERATOR_DEBUG_LEVEL=$<IF:$<CONFIG:Debug>,2,0>)
  target_compile_definitions(iko_tessellation_core PRIVATE
    IKO_TESSELLATION_ITERATOR_DEBUG_LEVEL=$<IF:$<CONFIG:Debug>,2,0>)
  # Normal archive extraction and function/data pruning; never /WHOLEARCHIVE.
  target_link_options(iko_tessellation PRIVATE /OPT:REF /OPT:ICF /INCREMENTAL:NO)
elseif(UNIX AND NOT APPLE)
  target_link_options(iko_tessellation PRIVATE "LINKER:--exclude-libs,ALL")
endif()
target_compile_features(iko_tessellation PUBLIC cxx_std_17)
target_include_directories(iko_tessellation PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
target_include_directories(iko_tessellation_core PRIVATE include third_party/triangle)
target_include_directories(iko_tessellation_core SYSTEM PRIVATE "${OpenCASCADE_INCLUDE_DIR}")
target_compile_definitions(iko_tessellation_core PRIVATE
  OCCT_STATIC_BUILD TRILIBRARY ANSI_DECLARATORS EXTERNAL_TEST NO_TIMER _CRT_SECURE_NO_WARNINGS)
# OCCT exports the static toolkit dependency closure of TKDEIGES.
target_link_libraries(iko_tessellation PRIVATE TKDEIGES)
target_link_libraries(iko_tessellation_core PRIVATE TKDEIGES)

file(GLOB_RECURSE sdk_headers CONFIGURE_DEPENDS include/tessellation/*.h)
file(GLOB_RECURSE sdk_sources CONFIGURE_DEPENDS src/*.cpp examples/*.cpp)
# Explicit dependencies also work with localized compiler diagnostics.
set_source_files_properties(${sdk_sources} PROPERTIES OBJECT_DEPENDS "${sdk_headers}")
set_source_files_properties(src/triangle_backend.c PROPERTIES OBJECT_DEPENDS
  "${CMAKE_CURRENT_SOURCE_DIR}/third_party/triangle/triangle.c;${CMAKE_CURRENT_SOURCE_DIR}/third_party/triangle/triangle.h")
if(MSVC)
  target_compile_options(iko_tessellation PRIVATE /fp:strict)
  target_compile_options(iko_tessellation_core PRIVATE /fp:strict)
  target_link_libraries(iko_tessellation PRIVATE bcrypt ws2_32)
else()
  target_compile_options(iko_tessellation PRIVATE -ffp-contract=off -fno-fast-math)
  target_compile_options(iko_tessellation_core PRIVATE -ffp-contract=off -fno-fast-math)
  target_link_libraries(iko_tessellation PRIVATE m)
endif()

if(IKO_BUILD_CLI)
  add_executable(tessellate src/tessellate.cpp src/options.cpp)
  target_link_libraries(tessellate PRIVATE iko::tessellation)
endif()
if(IKO_BUILD_EXAMPLES)
  add_executable(tessellation_example examples/example.cpp)
  target_link_libraries(tessellation_example PRIVATE iko::tessellation)
endif()

if(BUILD_TESTING OR IKO_BUILD_PROBE)
  add_subdirectory(tests)
endif()

install(CODE "list(LENGTH CMAKE_INSTALL_MANIFEST_FILES _iko_tessellation_manifest_start)")
install(TARGETS iko_tessellation EXPORT iko_tessellation_targets
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}/$<LOWER_CASE:$<CONFIG>>
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}/$<LOWER_CASE:$<CONFIG>>
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}/$<LOWER_CASE:$<CONFIG>>)
if(IKO_BUILD_CLI)
  install(TARGETS tessellate RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}/$<LOWER_CASE:$<CONFIG>>)
endif()
install(DIRECTORY include/tessellation DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
  PATTERN "private" EXCLUDE)
# Generate the standard export files, but install them individually below.
# The generator's bulk replacement rule must never run in the shared prefix.
install(CODE "message(FATAL_ERROR \"This internal component only generates exports; install the module normally\")"
  COMPONENT iko_tessellation_export_generation EXCLUDE_FROM_ALL)
install(EXPORT iko_tessellation_targets FILE iko-tessellation-targets.cmake NAMESPACE iko::
  DESTINATION lib/cmake/iko_tessellation
  COMPONENT iko_tessellation_export_generation EXCLUDE_FROM_ALL)
include(CMakePackageConfigHelpers)
configure_package_config_file(cmake/iko_tessellation-config.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/iko_tessellation-config.cmake
  INSTALL_DESTINATION lib/cmake/iko_tessellation)
install(FILES cmake/iko-config.cmake DESTINATION lib/cmake/iko)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/iko_tessellation-config.cmake
  DESTINATION lib/cmake/iko_tessellation)
install(FILES THIRD_PARTY_NOTICES.md DESTINATION ${CMAKE_INSTALL_DATADIR}/iko/tessellation)
install(FILES dependencies/occt_projection-8.0.1.patch
  DESTINATION ${CMAKE_INSTALL_DATADIR}/iko/tessellation/third_party/occt)
install(DIRECTORY third_party/ DESTINATION ${CMAKE_INSTALL_DATADIR}/iko/tessellation/third_party)
configure_file(cmake/install_tessellation.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/install_tessellation.in @ONLY)
file(GENERATE OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/install_tessellation-$<CONFIG>.cmake
  INPUT ${CMAKE_CURRENT_BINARY_DIR}/install_tessellation.in)
install(SCRIPT ${CMAKE_CURRENT_BINARY_DIR}/install_tessellation-$<CONFIG>.cmake)
