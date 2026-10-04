# The SDK's analytic projection compatibility requires this small patch on the
# exact official 8.0.1 source. Reject an unpatched installation with the same
# upstream version rather than silently changing complete meshes.
function(iko_check_occt_projection_policy result)
  set(_occt_projection_policy)
  set(_occt_policy_header "${OpenCASCADE_INCLUDE_DIR}/ProjLib_ComputeApproxOnPolarSurface.hxx")
  if(EXISTS "${_occt_policy_header}")
    file(STRINGS "${_occt_policy_header}" _occt_projection_policy
      REGEX "^#define OCCT_TESSELLATION_PROJECTION_POLICY 1$")
  endif()
  if(_occt_projection_policy)
    set(${result} TRUE PARENT_SCOPE)
  else()
    set(${result} FALSE PARENT_SCOPE)
  endif()
endfunction()
