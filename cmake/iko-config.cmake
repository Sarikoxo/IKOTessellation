# Shared entry point for every iko module. Install this file verbatim.
set(iko_FOUND TRUE)
unset(iko_NOT_FOUND_MESSAGE)

foreach(_iko_component IN LISTS iko_FIND_COMPONENTS)
  find_package(iko_${_iko_component} QUIET CONFIG
    PATHS "${CMAKE_CURRENT_LIST_DIR}/../iko_${_iko_component}" NO_DEFAULT_PATH)
  if(NOT iko_${_iko_component}_FOUND AND iko_FIND_REQUIRED_${_iko_component})
    set(iko_FOUND FALSE)
    string(APPEND iko_NOT_FOUND_MESSAGE
      "Required iko component '${_iko_component}' was not found. ")
    if(iko_${_iko_component}_NOT_FOUND_MESSAGE)
      string(APPEND iko_NOT_FOUND_MESSAGE "${iko_${_iko_component}_NOT_FOUND_MESSAGE} ")
    endif()
  endif()
endforeach()
unset(_iko_component)
