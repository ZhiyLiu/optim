# - Find a library installation or build tree.
# 
# The following variables are set if SRepIO is found.  
# If SRepIO is not found, SRepIO_FOUND is set to false.
#  SRepIO_FOUND         - Set to true when SRepIO is found.
#  SRepIO_USE_FILE      - CMake file to use SRepIO.
#  SRepIO_MAJOR_VERSION - The SRepIO major version number.
#  SRepIO_MINOR_VERSION - The SRepIO minor version number 
#                       (odd non-release).
#  SRepIO_BUILD_VERSION - The SRepIO patch level 
#                       (meaningless for odd minor).
#  SRepIO_INCLUDE_DIRS  - Include directories for SRepIO
#  SRepIO_LIBRARY_DIRS  - Link directories for SRepIO libraries
#  SRepIO_LIBRARIES     - List of libraries to link against
#
# The following cache entries must be set by the user to locate SRepIO:
#  SRepIO_DIR  - The directory containing SRepIOConfig.cmake.  
#             This is either the root of the build tree,
#             or the lib/SRepIO directory.  This is the 
#             only cache entry.


# Construct consitent error messages for use below.
SET(SRepIO_DIR_DESCRIPTION "directory containing SRepIOConfig.cmake.  This is either the root of the build tree, or PREFIX/lib/SRepIO for an installation.")
SET(SRepIO_NOT_FOUND_MESSAGE "SRepIO not found.  Set the SRepIO_DIR cmake cache entry to the ${SRepIO_DIR_DESCRIPTION}")

# Search only if the location is not already known.
IF(NOT SRepIO_DIR)
  # Get the system search path as a list.
  IF(UNIX)
    STRING(REGEX MATCHALL "[^:]+" SRepIO_DIR_SEARCH1 "$ENV{PATH}")
  ELSE(UNIX)
    STRING(REGEX REPLACE "\\\\" "/" SRepIO_DIR_SEARCH1 "$ENV{PATH}")
  ENDIF(UNIX)
  STRING(REGEX REPLACE "/;" ";" SRepIO_DIR_SEARCH2 "${SRepIO_DIR_SEARCH1}")

  # Construct a set of paths relative to the system search path.
  SET(SRepIO_DIR_SEARCH "")
  FOREACH(dir ${SRepIO_DIR_SEARCH2})
    SET(SRepIO_DIR_SEARCH ${SRepIO_DIR_SEARCH}
      ${dir}/../lib/lib64/SRepIO
      )
  ENDFOREACH(dir)

  #
  # Look for an installation or build tree.
  #
  FIND_PATH(SRepIO_DIR UseSRepIO.cmake
    # Look for an environment variable SRepIO_DIR.
    $ENV{SRepIO_DIR}

    # Look in places relative to the system executable search path.
    ${SRepIO_DIR_SEARCH}

    # Look in standard WIN install locations.
    "$ENV{ProgramFiles}/SRepIO"

    # Look in standard UNIX install locations.
    /usr/local/lib64/SRepIO
    /usr/lib64/SRepIO

    # Read from the CMakeSetup registry entries.  It is likely that
    # SRepIO will have been recently built.
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
    DOC "The ${SRepIO_DIR_DESCRIPTION}"
  )
ENDIF(NOT SRepIO_DIR)

# If SRepIO was found, load the configuration file to get the rest of the
# settings.
IF(SRepIO_DIR)
  # Make sure the SRepIOConfig.cmake file exists in the directory provided.
  IF(EXISTS ${SRepIO_DIR}/SRepIOConfig.cmake)

    # We found SRepIO.  Load the settings.
    SET(SRepIO_FOUND 1)
    INCLUDE(${SRepIO_DIR}/SRepIOConfig.cmake)

  ENDIF(EXISTS ${SRepIO_DIR}/SRepIOConfig.cmake)
ELSE(SRepIO_DIR)
  # We did not find SRepIO.
  SET(SRepIO_FOUND 0)
ENDIF(SRepIO_DIR)

#-----------------------------------------------------------------------------
IF(NOT SRepIO_FOUND)
  # SRepIO not found, explain to the user how to specify its location.
  IF(NOT SRepIO_FIND_QUIETLY)
    MESSAGE(FATAL_ERROR ${SRepIO_NOT_FOUND_MESSAGE})
  ELSE(NOT SRepIO_FIND_QUIETLY)
    IF(SRepIO_FIND_REQUIRED)
      MESSAGE(FATAL_ERROR ${SRepIO_NOT_FOUND_MESSAGE})
    ENDIF(SRepIO_FIND_REQUIRED)
  ENDIF(NOT SRepIO_FIND_QUIETLY)
ENDIF(NOT SRepIO_FOUND)
