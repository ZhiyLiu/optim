# - Find a library installation or build tree.
# 
# The following variables are set if paul_code is found.  
# If paul_code is not found, paul_code_FOUND is set to false.
#  paul_code_FOUND         - Set to true when paul_code is found.
#  paul_code_USE_FILE      - CMake file to use paul_code.
#  paul_code_MAJOR_VERSION - The paul_code major version number.
#  paul_code_MINOR_VERSION - The paul_code minor version number 
#                       (odd non-release).
#  paul_code_BUILD_VERSION - The paul_code patch level 
#                       (meaningless for odd minor).
#  paul_code_INCLUDE_DIRS  - Include directories for paul_code
#  paul_code_LIBRARY_DIRS  - Link directories for paul_code libraries
#  paul_code_LIBRARIES     - List of libraries to link against
#
# The following cache entries must be set by the user to locate paul_code:
#  paul_code_DIR  - The directory containing paul_codeConfig.cmake.  
#             This is either the root of the build tree,
#             or the lib/paul_code directory.  This is the 
#             only cache entry.


# Construct consitent error messages for use below.
SET(paul_code_DIR_DESCRIPTION "directory containing paul_codeConfig.cmake.  This is either the root of the build tree, or PREFIX/lib/paul_code for an installation.")
SET(paul_code_NOT_FOUND_MESSAGE "paul_code not found.  Set the paul_code_DIR cmake cache entry to the ${paul_code_DIR_DESCRIPTION}")

# Search only if the location is not already known.
IF(NOT paul_code_DIR)
  # Get the system search path as a list.
  IF(UNIX)
    STRING(REGEX MATCHALL "[^:]+" paul_code_DIR_SEARCH1 "$ENV{PATH}")
  ELSE(UNIX)
    STRING(REGEX REPLACE "\\\\" "/" paul_code_DIR_SEARCH1 "$ENV{PATH}")
  ENDIF(UNIX)
  STRING(REGEX REPLACE "/;" ";" paul_code_DIR_SEARCH2 "${paul_code_DIR_SEARCH1}")

  # Construct a set of paths relative to the system search path.
  SET(paul_code_DIR_SEARCH "")
  FOREACH(dir ${paul_code_DIR_SEARCH2})
    SET(paul_code_DIR_SEARCH ${paul_code_DIR_SEARCH}
      ${dir}/../lib/lib64/paul_code
      )
  ENDFOREACH(dir)

  #
  # Look for an installation or build tree.
  #
  FIND_PATH(paul_code_DIR Usepaul_code.cmake
    # Look for an environment variable paul_code_DIR.
    $ENV{paul_code_DIR}

    # Look in places relative to the system executable search path.
    ${paul_code_DIR_SEARCH}

    # Look in standard WIN install locations.
    "$ENV{ProgramFiles}/paul_code"

    # Look in standard UNIX install locations.
    /usr/local/lib64/paul_code
    /usr/lib64/paul_code

    # Read from the CMakeSetup registry entries.  It is likely that
    # paul_code will have been recently built.
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
    DOC "The ${paul_code_DIR_DESCRIPTION}"
  )
ENDIF(NOT paul_code_DIR)

# If paul_code was found, load the configuration file to get the rest of the
# settings.
IF(paul_code_DIR)
  # Make sure the paul_codeConfig.cmake file exists in the directory provided.
  IF(EXISTS ${paul_code_DIR}/paul_codeConfig.cmake)

    # We found paul_code.  Load the settings.
    SET(paul_code_FOUND 1)
    INCLUDE(${paul_code_DIR}/paul_codeConfig.cmake)

  ENDIF(EXISTS ${paul_code_DIR}/paul_codeConfig.cmake)
ELSE(paul_code_DIR)
  # We did not find paul_code.
  SET(paul_code_FOUND 0)
ENDIF(paul_code_DIR)

#-----------------------------------------------------------------------------
IF(NOT paul_code_FOUND)
  # paul_code not found, explain to the user how to specify its location.
  IF(NOT paul_code_FIND_QUIETLY)
    MESSAGE(FATAL_ERROR ${paul_code_NOT_FOUND_MESSAGE})
  ELSE(NOT paul_code_FIND_QUIETLY)
    IF(paul_code_FIND_REQUIRED)
      MESSAGE(FATAL_ERROR ${paul_code_NOT_FOUND_MESSAGE})
    ENDIF(paul_code_FIND_REQUIRED)
  ENDIF(NOT paul_code_FIND_QUIETLY)
ENDIF(NOT paul_code_FOUND)
