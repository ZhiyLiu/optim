# - Find a library installation or build tree.
# 
# The following variables are set if version is found.  
# If version is not found, version_FOUND is set to false.
#  version_FOUND         - Set to true when version is found.
#  version_USE_FILE      - CMake file to use version.
#  version_MAJOR_VERSION - The version major version number.
#  version_MINOR_VERSION - The version minor version number 
#                       (odd non-release).
#  version_BUILD_VERSION - The version patch level 
#                       (meaningless for odd minor).
#  version_INCLUDE_DIRS  - Include directories for version
#  version_LIBRARY_DIRS  - Link directories for version libraries
#  version_LIBRARIES     - List of libraries to link against
#
# The following cache entries must be set by the user to locate version:
#  version_DIR  - The directory containing versionConfig.cmake.  
#             This is either the root of the build tree,
#             or the lib/version directory.  This is the 
#             only cache entry.


# Construct consitent error messages for use below.
SET(version_DIR_DESCRIPTION "directory containing versionConfig.cmake.  This is either the root of the build tree, or PREFIX/lib/version for an installation.")
SET(version_NOT_FOUND_MESSAGE "version not found.  Set the version_DIR cmake cache entry to the ${version_DIR_DESCRIPTION}")

# Search only if the location is not already known.
IF(NOT version_DIR)
  # Get the system search path as a list.
  IF(UNIX)
    STRING(REGEX MATCHALL "[^:]+" version_DIR_SEARCH1 "$ENV{PATH}")
  ELSE(UNIX)
    STRING(REGEX REPLACE "\\\\" "/" version_DIR_SEARCH1 "$ENV{PATH}")
  ENDIF(UNIX)
  STRING(REGEX REPLACE "/;" ";" version_DIR_SEARCH2 "${version_DIR_SEARCH1}")

  # Construct a set of paths relative to the system search path.
  SET(version_DIR_SEARCH "")
  FOREACH(dir ${version_DIR_SEARCH2})
    SET(version_DIR_SEARCH ${version_DIR_SEARCH}
      ${dir}/../lib/lib64/version
      )
  ENDFOREACH(dir)

  #
  # Look for an installation or build tree.
  #
  FIND_PATH(version_DIR Useversion.cmake
    # Look for an environment variable version_DIR.
    $ENV{version_DIR}

    # Look in places relative to the system executable search path.
    ${version_DIR_SEARCH}

    # Look in standard WIN install locations.
    "$ENV{ProgramFiles}/version"

    # Look in standard UNIX install locations.
    /usr/local/lib64/version
    /usr/lib64/version

    # Read from the CMakeSetup registry entries.  It is likely that
    # version will have been recently built.
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
    DOC "The ${version_DIR_DESCRIPTION}"
  )
ENDIF(NOT version_DIR)

# If version was found, load the configuration file to get the rest of the
# settings.
IF(version_DIR)
  # Make sure the versionConfig.cmake file exists in the directory provided.
  IF(EXISTS ${version_DIR}/versionConfig.cmake)

    # We found version.  Load the settings.
    SET(version_FOUND 1)
    INCLUDE(${version_DIR}/versionConfig.cmake)

  ENDIF(EXISTS ${version_DIR}/versionConfig.cmake)
ELSE(version_DIR)
  # We did not find version.
  SET(version_FOUND 0)
ENDIF(version_DIR)

#-----------------------------------------------------------------------------
IF(NOT version_FOUND)
  # version not found, explain to the user how to specify its location.
  IF(NOT version_FIND_QUIETLY)
    MESSAGE(FATAL_ERROR ${version_NOT_FOUND_MESSAGE})
  ELSE(NOT version_FIND_QUIETLY)
    IF(version_FIND_REQUIRED)
      MESSAGE(FATAL_ERROR ${version_NOT_FOUND_MESSAGE})
    ENDIF(version_FIND_REQUIRED)
  ENDIF(NOT version_FIND_QUIETLY)
ENDIF(NOT version_FOUND)
