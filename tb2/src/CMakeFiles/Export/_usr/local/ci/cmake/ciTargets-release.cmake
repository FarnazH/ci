#----------------------------------------------------------------
# Generated CMake target import file for configuration "Release".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "ci" for configuration "Release"
set_property(TARGET ci APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(ci PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_RELEASE "CXX"
  IMPORTED_LOCATION_RELEASE "/usr/local/ci/lib/libci.a"
  )

list(APPEND _IMPORT_CHECK_TARGETS ci )
list(APPEND _IMPORT_CHECK_FILES_FOR_ci "/usr/local/ci/lib/libci.a" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
