# - Find a library installation or build tree.
# 
# The following variables are set if planes is found.  
# If planes is not found, planes_FOUND is set to false.
#  planes_FOUND         - Set to true when planes is found.
#  planes_USE_FILE      - CMake file to use planes.
#  planes_MAJOR_VERSION - The planes major version number.
#  planes_MINOR_VERSION - The planes minor version number 
#                       (odd non-release).
#  planes_BUILD_VERSION - The planes patch level 
#                       (meaningless for odd minor).
#  planes_INCLUDE_DIRS  - Include directories for planes
#  planes_LIBRARY_DIRS  - Link directories for planes libraries
#  planes_LIBRARIES     - List of libraries to link against
#
# The following cache entries must be set by the user to locate planes:
#  planes_DIR  - The directory containing planesConfig.cmake.  
#             This is either the root of the build tree,
#             or the lib/planes directory.  This is the 
#             only cache entry.


# Construct consitent error messages for use below.
SET(planes_DIR_DESCRIPTION "directory containing planesConfig.cmake.  This is either the root of the build tree, or PREFIX/lib/planes for an installation.")
SET(planes_NOT_FOUND_MESSAGE "planes not found.  Set the planes_DIR cmake cache entry to the ${planes_DIR_DESCRIPTION}")

# Search only if the location is not already known.
IF(NOT planes_DIR)
  # Get the system search path as a list.
  IF(UNIX)
    STRING(REGEX MATCHALL "[^:]+" planes_DIR_SEARCH1 "$ENV{PATH}")
  ELSE(UNIX)
    STRING(REGEX REPLACE "\\\\" "/" planes_DIR_SEARCH1 "$ENV{PATH}")
  ENDIF(UNIX)
  STRING(REGEX REPLACE "/;" ";" planes_DIR_SEARCH2 "${planes_DIR_SEARCH1}")

  # Construct a set of paths relative to the system search path.
  SET(planes_DIR_SEARCH "")
  FOREACH(dir ${planes_DIR_SEARCH2})
    SET(planes_DIR_SEARCH ${planes_DIR_SEARCH}
      ${dir}/../lib/lib64/planes
      )
  ENDFOREACH(dir)

  #
  # Look for an installation or build tree.
  #
  FIND_PATH(planes_DIR Useplanes.cmake
    # Look for an environment variable planes_DIR.
    $ENV{planes_DIR}

    # Look in places relative to the system executable search path.
    ${planes_DIR_SEARCH}

    # Look in standard WIN install locations.
    "$ENV{ProgramFiles}/planes"

    # Look in standard UNIX install locations.
    /usr/local/lib64/planes
    /usr/lib64/planes

    # Read from the CMakeSetup registry entries.  It is likely that
    # planes will have been recently built.
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
    DOC "The ${planes_DIR_DESCRIPTION}"
  )
ENDIF(NOT planes_DIR)

# If planes was found, load the configuration file to get the rest of the
# settings.
IF(planes_DIR)
  # Make sure the planesConfig.cmake file exists in the directory provided.
  IF(EXISTS ${planes_DIR}/planesConfig.cmake)

    # We found planes.  Load the settings.
    SET(planes_FOUND 1)
    INCLUDE(${planes_DIR}/planesConfig.cmake)

  ENDIF(EXISTS ${planes_DIR}/planesConfig.cmake)
ELSE(planes_DIR)
  # We did not find planes.
  SET(planes_FOUND 0)
ENDIF(planes_DIR)

#-----------------------------------------------------------------------------
IF(NOT planes_FOUND)
  # planes not found, explain to the user how to specify its location.
  IF(NOT planes_FIND_QUIETLY)
    MESSAGE(FATAL_ERROR ${planes_NOT_FOUND_MESSAGE})
  ELSE(NOT planes_FIND_QUIETLY)
    IF(planes_FIND_REQUIRED)
      MESSAGE(FATAL_ERROR ${planes_NOT_FOUND_MESSAGE})
    ENDIF(planes_FIND_REQUIRED)
  ENDIF(NOT planes_FIND_QUIETLY)
ENDIF(NOT planes_FOUND)
