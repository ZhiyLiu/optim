# - Find a library installation or build tree.
# 
# The following variables are set if flvw-1.0 is found.  
# If flvw-1.0 is not found, flvw-1.0_FOUND is set to false.
#  flvw-1.0_FOUND         - Set to true when flvw-1.0 is found.
#  flvw-1.0_USE_FILE      - CMake file to use flvw-1.0.
#  flvw-1.0_MAJOR_VERSION - The flvw-1.0 major version number.
#  flvw-1.0_MINOR_VERSION - The flvw-1.0 minor version number 
#                       (odd non-release).
#  flvw-1.0_BUILD_VERSION - The flvw-1.0 patch level 
#                       (meaningless for odd minor).
#  flvw-1.0_INCLUDE_DIRS  - Include directories for flvw-1.0
#  flvw-1.0_LIBRARY_DIRS  - Link directories for flvw-1.0 libraries
#  flvw-1.0_LIBRARIES     - List of libraries to link against
#
# The following cache entries must be set by the user to locate flvw-1.0:
#  flvw-1.0_DIR  - The directory containing flvw-1.0Config.cmake.  
#             This is either the root of the build tree,
#             or the lib/flvw-1.0 directory.  This is the 
#             only cache entry.


# Construct consitent error messages for use below.
SET(flvw-1.0_DIR_DESCRIPTION "directory containing flvw-1.0Config.cmake.  This is either the root of the build tree, or PREFIX/lib/flvw-1.0 for an installation.")
SET(flvw-1.0_NOT_FOUND_MESSAGE "flvw-1.0 not found.  Set the flvw-1.0_DIR cmake cache entry to the ${flvw-1.0_DIR_DESCRIPTION}")

# Search only if the location is not already known.
IF(NOT flvw-1.0_DIR)
  # Get the system search path as a list.
  IF(UNIX)
    STRING(REGEX MATCHALL "[^:]+" flvw-1.0_DIR_SEARCH1 "$ENV{PATH}")
  ELSE(UNIX)
    STRING(REGEX REPLACE "\\\\" "/" flvw-1.0_DIR_SEARCH1 "$ENV{PATH}")
  ENDIF(UNIX)
  STRING(REGEX REPLACE "/;" ";" flvw-1.0_DIR_SEARCH2 "${flvw-1.0_DIR_SEARCH1}")

  # Construct a set of paths relative to the system search path.
  SET(flvw-1.0_DIR_SEARCH "")
  FOREACH(dir ${flvw-1.0_DIR_SEARCH2})
    SET(flvw-1.0_DIR_SEARCH ${flvw-1.0_DIR_SEARCH}
      ${dir}/../lib/lib64/flvw-1.0
      )
  ENDFOREACH(dir)

  #
  # Look for an installation or build tree.
  #
  FIND_PATH(flvw-1.0_DIR Useflvw-1.0.cmake
    # Look for an environment variable flvw-1.0_DIR.
    $ENV{flvw-1.0_DIR}

    # Look in places relative to the system executable search path.
    ${flvw-1.0_DIR_SEARCH}

    # Look in standard WIN install locations.
    "$ENV{ProgramFiles}/flvw-1.0"

    # Look in standard UNIX install locations.
    /usr/local/lib64/flvw-1.0
    /usr/lib64/flvw-1.0

    # Read from the CMakeSetup registry entries.  It is likely that
    # flvw-1.0 will have been recently built.
    [HKEY_CURRENT_USER\\Software\\Kitware\\CMakeSetup\\Settings\\StartPath;WhereBuild1]
    [HKEY_CURRENT_USER\\Software\\Kitware\\CMakeSetup\\Settings\\StartPath;WhereBuild2]
    [HKEY_CURRENT_USER\\Software\\Kitware\\CMakeSetup\\Settings\\StartPath;WhereBuild3]
    [HKEY_CURRENT_USER\\Software\\Kitware\\CMakeSetup\\Settings\\StartPath;WhereBuild4]
    [HKEY_CURRENT_USER\\Software\\Kitware\\CMakeSetup\\Settings\\StartPath;WhereBuild5]
    [HKEY_CURRENT_USER\\Software\\Kitware\\CMakeSetup\\Settings\\StartPath;WhereBuild6]
    [HKEY_CURRENT_USER\\Software\\Kitware\\CMakeSetup\\Settings\\StartPath;WhereBuild7]
    [HKEY_CURRENT_USER\\Software\\Kitware\\CMakeSetup\\Settings\\StartPath;WhereBuild8]
    [HKEY_CURRENT_USER\\Software\\Kitware\\CMakeSetup\\Settings\\StartPath;WhereBuild9]
    [HKEY_CURRENT_USER\\Software\\Kitware\\CMakeSetup\\Settings\\StartPath;WhereBuild10]

    # Help the user find it if we cannot.
    DOC "The ${flvw-1.0_DIR_DESCRIPTION}"
  )
ENDIF(NOT flvw-1.0_DIR)

# If flvw-1.0 was found, load the configuration file to get the rest of the
# settings.
IF(flvw-1.0_DIR)
  # Make sure the flvw-1.0Config.cmake file exists in the directory provided.
  IF(EXISTS ${flvw-1.0_DIR}/flvw-1.0Config.cmake)

    # We found flvw-1.0.  Load the settings.
    SET(flvw-1.0_FOUND 1)
    INCLUDE(${flvw-1.0_DIR}/flvw-1.0Config.cmake)

  ENDIF(EXISTS ${flvw-1.0_DIR}/flvw-1.0Config.cmake)
ELSE(flvw-1.0_DIR)
  # We did not find flvw-1.0.
  SET(flvw-1.0_FOUND 0)
ENDIF(flvw-1.0_DIR)

#-----------------------------------------------------------------------------
IF(NOT flvw-1.0_FOUND)
  # flvw-1.0 not found, explain to the user how to specify its location.
  IF(NOT flvw-1.0_FIND_QUIETLY)
    MESSAGE(FATAL_ERROR ${flvw-1.0_NOT_FOUND_MESSAGE})
  ELSE(NOT flvw-1.0_FIND_QUIETLY)
    IF(flvw-1.0_FIND_REQUIRED)
      MESSAGE(FATAL_ERROR ${flvw-1.0_NOT_FOUND_MESSAGE})
    ENDIF(flvw-1.0_FIND_REQUIRED)
  ENDIF(NOT flvw-1.0_FIND_QUIETLY)
ENDIF(NOT flvw-1.0_FOUND)
