# - Find a library installation or build tree.
# 
# The following variables are set if SRepInterpolation is found.  
# If SRepInterpolation is not found, SRepInterpolation_FOUND is set to false.
#  SRepInterpolation_FOUND         - Set to true when SRepInterpolation is found.
#  SRepInterpolation_USE_FILE      - CMake file to use SRepInterpolation.
#  SRepInterpolation_MAJOR_VERSION - The SRepInterpolation major version number.
#  SRepInterpolation_MINOR_VERSION - The SRepInterpolation minor version number 
#                       (odd non-release).
#  SRepInterpolation_BUILD_VERSION - The SRepInterpolation patch level 
#                       (meaningless for odd minor).
#  SRepInterpolation_INCLUDE_DIRS  - Include directories for SRepInterpolation
#  SRepInterpolation_LIBRARY_DIRS  - Link directories for SRepInterpolation libraries
#  SRepInterpolation_LIBRARIES     - List of libraries to link against
#
# The following cache entries must be set by the user to locate SRepInterpolation:
#  SRepInterpolation_DIR  - The directory containing SRepInterpolationConfig.cmake.  
#             This is either the root of the build tree,
#             or the lib/SRepInterpolation directory.  This is the 
#             only cache entry.


# Construct consitent error messages for use below.
SET(SRepInterpolation_DIR_DESCRIPTION "directory containing SRepInterpolationConfig.cmake.  This is either the root of the build tree, or PREFIX/lib/SRepInterpolation for an installation.")
SET(SRepInterpolation_NOT_FOUND_MESSAGE "SRepInterpolation not found.  Set the SRepInterpolation_DIR cmake cache entry to the ${SRepInterpolation_DIR_DESCRIPTION}")

# Search only if the location is not already known.
IF(NOT SRepInterpolation_DIR)
  # Get the system search path as a list.
  IF(UNIX)
    STRING(REGEX MATCHALL "[^:]+" SRepInterpolation_DIR_SEARCH1 "$ENV{PATH}")
  ELSE(UNIX)
    STRING(REGEX REPLACE "\\\\" "/" SRepInterpolation_DIR_SEARCH1 "$ENV{PATH}")
  ENDIF(UNIX)
  STRING(REGEX REPLACE "/;" ";" SRepInterpolation_DIR_SEARCH2 "${SRepInterpolation_DIR_SEARCH1}")

  # Construct a set of paths relative to the system search path.
  SET(SRepInterpolation_DIR_SEARCH "")
  FOREACH(dir ${SRepInterpolation_DIR_SEARCH2})
    SET(SRepInterpolation_DIR_SEARCH ${SRepInterpolation_DIR_SEARCH}
      ${dir}/../lib/lib64/SRepInterpolation
      )
  ENDFOREACH(dir)

  #
  # Look for an installation or build tree.
  #
  FIND_PATH(SRepInterpolation_DIR UseSRepInterpolation.cmake
    # Look for an environment variable SRepInterpolation_DIR.
    $ENV{SRepInterpolation_DIR}

    # Look in places relative to the system executable search path.
    ${SRepInterpolation_DIR_SEARCH}

    # Look in standard WIN install locations.
    "$ENV{ProgramFiles}/SRepInterpolation"

    # Look in standard UNIX install locations.
    /usr/local/lib64/SRepInterpolation
    /usr/lib64/SRepInterpolation

    # Read from the CMakeSetup registry entries.  It is likely that
    # SRepInterpolation will have been recently built.
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
    DOC "The ${SRepInterpolation_DIR_DESCRIPTION}"
  )
ENDIF(NOT SRepInterpolation_DIR)

# If SRepInterpolation was found, load the configuration file to get the rest of the
# settings.
IF(SRepInterpolation_DIR)
  # Make sure the SRepInterpolationConfig.cmake file exists in the directory provided.
  IF(EXISTS ${SRepInterpolation_DIR}/SRepInterpolationConfig.cmake)

    # We found SRepInterpolation.  Load the settings.
    SET(SRepInterpolation_FOUND 1)
    INCLUDE(${SRepInterpolation_DIR}/SRepInterpolationConfig.cmake)

  ENDIF(EXISTS ${SRepInterpolation_DIR}/SRepInterpolationConfig.cmake)
ELSE(SRepInterpolation_DIR)
  # We did not find SRepInterpolation.
  SET(SRepInterpolation_FOUND 0)
ENDIF(SRepInterpolation_DIR)

#-----------------------------------------------------------------------------
IF(NOT SRepInterpolation_FOUND)
  # SRepInterpolation not found, explain to the user how to specify its location.
  IF(NOT SRepInterpolation_FIND_QUIETLY)
    MESSAGE(FATAL_ERROR ${SRepInterpolation_NOT_FOUND_MESSAGE})
  ELSE(NOT SRepInterpolation_FIND_QUIETLY)
    IF(SRepInterpolation_FIND_REQUIRED)
      MESSAGE(FATAL_ERROR ${SRepInterpolation_NOT_FOUND_MESSAGE})
    ENDIF(SRepInterpolation_FIND_REQUIRED)
  ENDIF(NOT SRepInterpolation_FIND_QUIETLY)
ENDIF(NOT SRepInterpolation_FOUND)
