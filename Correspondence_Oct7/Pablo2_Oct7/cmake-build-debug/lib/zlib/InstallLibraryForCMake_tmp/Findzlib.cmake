# - Find a library installation or build tree.
# 
# The following variables are set if zlib is found.  
# If zlib is not found, zlib_FOUND is set to false.
#  zlib_FOUND         - Set to true when zlib is found.
#  zlib_USE_FILE      - CMake file to use zlib.
#  zlib_MAJOR_VERSION - The zlib major version number.
#  zlib_MINOR_VERSION - The zlib minor version number 
#                       (odd non-release).
#  zlib_BUILD_VERSION - The zlib patch level 
#                       (meaningless for odd minor).
#  zlib_INCLUDE_DIRS  - Include directories for zlib
#  zlib_LIBRARY_DIRS  - Link directories for zlib libraries
#  zlib_LIBRARIES     - List of libraries to link against
#
# The following cache entries must be set by the user to locate zlib:
#  zlib_DIR  - The directory containing zlibConfig.cmake.  
#             This is either the root of the build tree,
#             or the lib/zlib directory.  This is the 
#             only cache entry.


# Construct consitent error messages for use below.
SET(zlib_DIR_DESCRIPTION "directory containing zlibConfig.cmake.  This is either the root of the build tree, or PREFIX/lib/zlib for an installation.")
SET(zlib_NOT_FOUND_MESSAGE "zlib not found.  Set the zlib_DIR cmake cache entry to the ${zlib_DIR_DESCRIPTION}")

# Search only if the location is not already known.
IF(NOT zlib_DIR)
  # Get the system search path as a list.
  IF(UNIX)
    STRING(REGEX MATCHALL "[^:]+" zlib_DIR_SEARCH1 "$ENV{PATH}")
  ELSE(UNIX)
    STRING(REGEX REPLACE "\\\\" "/" zlib_DIR_SEARCH1 "$ENV{PATH}")
  ENDIF(UNIX)
  STRING(REGEX REPLACE "/;" ";" zlib_DIR_SEARCH2 "${zlib_DIR_SEARCH1}")

  # Construct a set of paths relative to the system search path.
  SET(zlib_DIR_SEARCH "")
  FOREACH(dir ${zlib_DIR_SEARCH2})
    SET(zlib_DIR_SEARCH ${zlib_DIR_SEARCH}
      ${dir}/../lib/lib64/zlib
      )
  ENDFOREACH(dir)

  #
  # Look for an installation or build tree.
  #
  FIND_PATH(zlib_DIR Usezlib.cmake
    # Look for an environment variable zlib_DIR.
    $ENV{zlib_DIR}

    # Look in places relative to the system executable search path.
    ${zlib_DIR_SEARCH}

    # Look in standard WIN install locations.
    "$ENV{ProgramFiles}/zlib"

    # Look in standard UNIX install locations.
    /usr/local/lib64/zlib
    /usr/lib64/zlib

    # Read from the CMakeSetup registry entries.  It is likely that
    # zlib will have been recently built.
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
    DOC "The ${zlib_DIR_DESCRIPTION}"
  )
ENDIF(NOT zlib_DIR)

# If zlib was found, load the configuration file to get the rest of the
# settings.
IF(zlib_DIR)
  # Make sure the zlibConfig.cmake file exists in the directory provided.
  IF(EXISTS ${zlib_DIR}/zlibConfig.cmake)

    # We found zlib.  Load the settings.
    SET(zlib_FOUND 1)
    INCLUDE(${zlib_DIR}/zlibConfig.cmake)

  ENDIF(EXISTS ${zlib_DIR}/zlibConfig.cmake)
ELSE(zlib_DIR)
  # We did not find zlib.
  SET(zlib_FOUND 0)
ENDIF(zlib_DIR)

#-----------------------------------------------------------------------------
IF(NOT zlib_FOUND)
  # zlib not found, explain to the user how to specify its location.
  IF(NOT zlib_FIND_QUIETLY)
    MESSAGE(FATAL_ERROR ${zlib_NOT_FOUND_MESSAGE})
  ELSE(NOT zlib_FIND_QUIETLY)
    IF(zlib_FIND_REQUIRED)
      MESSAGE(FATAL_ERROR ${zlib_NOT_FOUND_MESSAGE})
    ENDIF(zlib_FIND_REQUIRED)
  ENDIF(NOT zlib_FIND_QUIETLY)
ENDIF(NOT zlib_FOUND)
