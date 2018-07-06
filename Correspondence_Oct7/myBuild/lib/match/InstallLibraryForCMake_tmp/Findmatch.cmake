# - Find a library installation or build tree.
# 
# The following variables are set if match is found.  
# If match is not found, match_FOUND is set to false.
#  match_FOUND         - Set to true when match is found.
#  match_USE_FILE      - CMake file to use match.
#  match_MAJOR_VERSION - The match major version number.
#  match_MINOR_VERSION - The match minor version number 
#                       (odd non-release).
#  match_BUILD_VERSION - The match patch level 
#                       (meaningless for odd minor).
#  match_INCLUDE_DIRS  - Include directories for match
#  match_LIBRARY_DIRS  - Link directories for match libraries
#  match_LIBRARIES     - List of libraries to link against
#
# The following cache entries must be set by the user to locate match:
#  match_DIR  - The directory containing matchConfig.cmake.  
#             This is either the root of the build tree,
#             or the lib/match directory.  This is the 
#             only cache entry.


# Construct consitent error messages for use below.
SET(match_DIR_DESCRIPTION "directory containing matchConfig.cmake.  This is either the root of the build tree, or PREFIX/lib/match for an installation.")
SET(match_NOT_FOUND_MESSAGE "match not found.  Set the match_DIR cmake cache entry to the ${match_DIR_DESCRIPTION}")

# Search only if the location is not already known.
IF(NOT match_DIR)
  # Get the system search path as a list.
  IF(UNIX)
    STRING(REGEX MATCHALL "[^:]+" match_DIR_SEARCH1 "$ENV{PATH}")
  ELSE(UNIX)
    STRING(REGEX REPLACE "\\\\" "/" match_DIR_SEARCH1 "$ENV{PATH}")
  ENDIF(UNIX)
  STRING(REGEX REPLACE "/;" ";" match_DIR_SEARCH2 "${match_DIR_SEARCH1}")

  # Construct a set of paths relative to the system search path.
  SET(match_DIR_SEARCH "")
  FOREACH(dir ${match_DIR_SEARCH2})
    SET(match_DIR_SEARCH ${match_DIR_SEARCH}
      ${dir}/../lib/lib64/match
      )
  ENDFOREACH(dir)

  #
  # Look for an installation or build tree.
  #
  FIND_PATH(match_DIR Usematch.cmake
    # Look for an environment variable match_DIR.
    $ENV{match_DIR}

    # Look in places relative to the system executable search path.
    ${match_DIR_SEARCH}

    # Look in standard WIN install locations.
    "$ENV{ProgramFiles}/match"

    # Look in standard UNIX install locations.
    /usr/local/lib64/match
    /usr/lib64/match

    # Read from the CMakeSetup registry entries.  It is likely that
    # match will have been recently built.
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
    DOC "The ${match_DIR_DESCRIPTION}"
  )
ENDIF(NOT match_DIR)

# If match was found, load the configuration file to get the rest of the
# settings.
IF(match_DIR)
  # Make sure the matchConfig.cmake file exists in the directory provided.
  IF(EXISTS ${match_DIR}/matchConfig.cmake)

    # We found match.  Load the settings.
    SET(match_FOUND 1)
    INCLUDE(${match_DIR}/matchConfig.cmake)

  ENDIF(EXISTS ${match_DIR}/matchConfig.cmake)
ELSE(match_DIR)
  # We did not find match.
  SET(match_FOUND 0)
ENDIF(match_DIR)

#-----------------------------------------------------------------------------
IF(NOT match_FOUND)
  # match not found, explain to the user how to specify its location.
  IF(NOT match_FIND_QUIETLY)
    MESSAGE(FATAL_ERROR ${match_NOT_FOUND_MESSAGE})
  ELSE(NOT match_FIND_QUIETLY)
    IF(match_FIND_REQUIRED)
      MESSAGE(FATAL_ERROR ${match_NOT_FOUND_MESSAGE})
    ENDIF(match_FIND_REQUIRED)
  ENDIF(NOT match_FIND_QUIETLY)
ENDIF(NOT match_FOUND)
