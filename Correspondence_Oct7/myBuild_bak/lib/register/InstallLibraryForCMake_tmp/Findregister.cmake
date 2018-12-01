# - Find a library installation or build tree.
# 
# The following variables are set if register is found.  
# If register is not found, register_FOUND is set to false.
#  register_FOUND         - Set to true when register is found.
#  register_USE_FILE      - CMake file to use register.
#  register_MAJOR_VERSION - The register major version number.
#  register_MINOR_VERSION - The register minor version number 
#                       (odd non-release).
#  register_BUILD_VERSION - The register patch level 
#                       (meaningless for odd minor).
#  register_INCLUDE_DIRS  - Include directories for register
#  register_LIBRARY_DIRS  - Link directories for register libraries
#  register_LIBRARIES     - List of libraries to link against
#
# The following cache entries must be set by the user to locate register:
#  register_DIR  - The directory containing registerConfig.cmake.  
#             This is either the root of the build tree,
#             or the lib/register directory.  This is the 
#             only cache entry.


# Construct consitent error messages for use below.
SET(register_DIR_DESCRIPTION "directory containing registerConfig.cmake.  This is either the root of the build tree, or PREFIX/lib/register for an installation.")
SET(register_NOT_FOUND_MESSAGE "register not found.  Set the register_DIR cmake cache entry to the ${register_DIR_DESCRIPTION}")

# Search only if the location is not already known.
IF(NOT register_DIR)
  # Get the system search path as a list.
  IF(UNIX)
    STRING(REGEX MATCHALL "[^:]+" register_DIR_SEARCH1 "$ENV{PATH}")
  ELSE(UNIX)
    STRING(REGEX REPLACE "\\\\" "/" register_DIR_SEARCH1 "$ENV{PATH}")
  ENDIF(UNIX)
  STRING(REGEX REPLACE "/;" ";" register_DIR_SEARCH2 "${register_DIR_SEARCH1}")

  # Construct a set of paths relative to the system search path.
  SET(register_DIR_SEARCH "")
  FOREACH(dir ${register_DIR_SEARCH2})
    SET(register_DIR_SEARCH ${register_DIR_SEARCH}
      ${dir}/../lib/lib64/register
      )
  ENDFOREACH(dir)

  #
  # Look for an installation or build tree.
  #
  FIND_PATH(register_DIR Useregister.cmake
    # Look for an environment variable register_DIR.
    $ENV{register_DIR}

    # Look in places relative to the system executable search path.
    ${register_DIR_SEARCH}

    # Look in standard WIN install locations.
    "$ENV{ProgramFiles}/register"

    # Look in standard UNIX install locations.
    /usr/local/lib64/register
    /usr/lib64/register

    # Read from the CMakeSetup registry entries.  It is likely that
    # register will have been recently built.
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
    DOC "The ${register_DIR_DESCRIPTION}"
  )
ENDIF(NOT register_DIR)

# If register was found, load the configuration file to get the rest of the
# settings.
IF(register_DIR)
  # Make sure the registerConfig.cmake file exists in the directory provided.
  IF(EXISTS ${register_DIR}/registerConfig.cmake)

    # We found register.  Load the settings.
    SET(register_FOUND 1)
    INCLUDE(${register_DIR}/registerConfig.cmake)

  ENDIF(EXISTS ${register_DIR}/registerConfig.cmake)
ELSE(register_DIR)
  # We did not find register.
  SET(register_FOUND 0)
ENDIF(register_DIR)

#-----------------------------------------------------------------------------
IF(NOT register_FOUND)
  # register not found, explain to the user how to specify its location.
  IF(NOT register_FIND_QUIETLY)
    MESSAGE(FATAL_ERROR ${register_NOT_FOUND_MESSAGE})
  ELSE(NOT register_FIND_QUIETLY)
    IF(register_FIND_REQUIRED)
      MESSAGE(FATAL_ERROR ${register_NOT_FOUND_MESSAGE})
    ENDIF(register_FIND_REQUIRED)
  ENDIF(NOT register_FIND_QUIETLY)
ENDIF(NOT register_FOUND)
