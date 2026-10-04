cmake_minimum_required(VERSION 3.20)
project(tessellation_consumer LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(NOT DEFINED CMAKE_MSVC_RUNTIME_LIBRARY)
  set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>DLL")
endif()
find_package(iko CONFIG REQUIRED COMPONENTS tessellation)
add_executable(tessellation_consumer example.cpp)
target_link_libraries(tessellation_consumer PRIVATE iko::tessellation)
add_custom_command(TARGET tessellation_consumer POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E copy_if_different $<TARGET_FILE:iko::tessellation>
    $<TARGET_FILE_DIR:tessellation_consumer>)
get_target_property(sdk_include_dir iko::tessellation INTERFACE_INCLUDE_DIRECTORIES)
# Keep installed-header changes visible even with localized MSVC diagnostics.
set_source_files_properties(example.cpp PROPERTIES OBJECT_DEPENDS
  "${sdk_include_dir}/tessellation/tessellation.h")
