# - Find a library installation or build tree.
# 
# The following variables are set if calquantum is found.  
# If calquantum is not found, calquantum_FOUND is set to false.
#  calquantum_FOUND         - Set to true when calquantum is found.
#  calquantum_USE_FILE      - CMake file to use calquantum.
#  calquantum_MAJOR_VERSION - The calquantum major version number.
#  calquantum_MINOR_VERSION - The calquantum minor version number 
#                       (odd non-release).
#  calquantum_BUILD_VERSION - The calquantum patch level 
#                       (meaningless for odd minor).
#  calquantum_INCLUDE_DIRS  - Include directories for calquantum
#  calquantum_LIBRARY_DIRS  - Link directories for calquantum libraries
#  calquantum_LIBRARIES     - List of libraries to link against
#
# The following cache entries must be set by the user to locate calquantum:
#  calquantum_DIR  - The directory containing calquantumConfig.cmake.  
#             This is either the root of the build tree,
#             or the lib/calquantum directory.  This is the 
#             only cache entry.


# Construct consitent error messages for use below.
SET(calquantum_DIR_DESCRIPTION "directory containing calquantumConfig.cmake.  This is either the root of the build tree, or PREFIX/lib/calquantum for an installation.")
SET(calquantum_NOT_FOUND_MESSAGE "calquantum not found.  Set the calquantum_DIR cmake cache entry to the ${calquantum_DIR_DESCRIPTION}")

# Search only if the location is not already known.
IF(NOT calquantum_DIR)
  # Get the system search path as a list.
  IF(UNIX)
    STRING(REGEX MATCHALL "[^:]+" calquantum_DIR_SEARCH1 "$ENV{PATH}")
  ELSE(UNIX)
    STRING(REGEX REPLACE "\\\\" "/" calquantum_DIR_SEARCH1 "$ENV{PATH}")
  ENDIF(UNIX)
  STRING(REGEX REPLACE "/;" ";" calquantum_DIR_SEARCH2 "${calquantum_DIR_SEARCH1}")

  # Construct a set of paths relative to the system search path.
  SET(calquantum_DIR_SEARCH "")
  FOREACH(dir ${calquantum_DIR_SEARCH2})
    SET(calquantum_DIR_SEARCH ${calquantum_DIR_SEARCH}
      ${dir}/../lib/lib64/calquantum
      )
  ENDFOREACH(dir)

  #
  # Look for an installation or build tree.
  #
  FIND_PATH(calquantum_DIR Usecalquantum.cmake
    # Look for an environment variable calquantum_DIR.
    $ENV{calquantum_DIR}

    # Look in places relative to the system executable search path.
    ${calquantum_DIR_SEARCH}

    # Look in standard WIN install locations.
    "$ENV{ProgramFiles}/calquantum"

    # Look in standard UNIX install locations.
    /usr/local/lib64/calquantum
    /usr/lib64/calquantum

    # Read from the CMakeSetup registry entries.  It is likely that
    # calquantum will have been recently built.
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
    DOC "The ${calquantum_DIR_DESCRIPTION}"
  )
ENDIF(NOT calquantum_DIR)

# If calquantum was found, load the configuration file to get the rest of the
# settings.
IF(calquantum_DIR)
  # Make sure the calquantumConfig.cmake file exists in the directory provided.
  IF(EXISTS ${calquantum_DIR}/calquantumConfig.cmake)

    # We found calquantum.  Load the settings.
    SET(calquantum_FOUND 1)
    INCLUDE(${calquantum_DIR}/calquantumConfig.cmake)

  ENDIF(EXISTS ${calquantum_DIR}/calquantumConfig.cmake)
ELSE(calquantum_DIR)
  # We did not find calquantum.
  SET(calquantum_FOUND 0)
ENDIF(calquantum_DIR)

#-----------------------------------------------------------------------------
IF(NOT calquantum_FOUND)
  # calquantum not found, explain to the user how to specify its location.
  IF(NOT calquantum_FIND_QUIETLY)
    MESSAGE(FATAL_ERROR ${calquantum_NOT_FOUND_MESSAGE})
  ELSE(NOT calquantum_FIND_QUIETLY)
    IF(calquantum_FIND_REQUIRED)
      MESSAGE(FATAL_ERROR ${calquantum_NOT_FOUND_MESSAGE})
    ENDIF(calquantum_FIND_REQUIRED)
  ENDIF(NOT calquantum_FIND_QUIETLY)
ENDIF(NOT calquantum_FOUND)
