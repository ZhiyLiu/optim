# - Find a library installation or build tree.
# 
# The following variables are set if SRepVisualization is found.  
# If SRepVisualization is not found, SRepVisualization_FOUND is set to false.
#  SRepVisualization_FOUND         - Set to true when SRepVisualization is found.
#  SRepVisualization_USE_FILE      - CMake file to use SRepVisualization.
#  SRepVisualization_MAJOR_VERSION - The SRepVisualization major version number.
#  SRepVisualization_MINOR_VERSION - The SRepVisualization minor version number 
#                       (odd non-release).
#  SRepVisualization_BUILD_VERSION - The SRepVisualization patch level 
#                       (meaningless for odd minor).
#  SRepVisualization_INCLUDE_DIRS  - Include directories for SRepVisualization
#  SRepVisualization_LIBRARY_DIRS  - Link directories for SRepVisualization libraries
#  SRepVisualization_LIBRARIES     - List of libraries to link against
#
# The following cache entries must be set by the user to locate SRepVisualization:
#  SRepVisualization_DIR  - The directory containing SRepVisualizationConfig.cmake.  
#             This is either the root of the build tree,
#             or the lib/SRepVisualization directory.  This is the 
#             only cache entry.


# Construct consitent error messages for use below.
SET(SRepVisualization_DIR_DESCRIPTION "directory containing SRepVisualizationConfig.cmake.  This is either the root of the build tree, or PREFIX/lib/SRepVisualization for an installation.")
SET(SRepVisualization_NOT_FOUND_MESSAGE "SRepVisualization not found.  Set the SRepVisualization_DIR cmake cache entry to the ${SRepVisualization_DIR_DESCRIPTION}")

# Search only if the location is not already known.
IF(NOT SRepVisualization_DIR)
  # Get the system search path as a list.
  IF(UNIX)
    STRING(REGEX MATCHALL "[^:]+" SRepVisualization_DIR_SEARCH1 "$ENV{PATH}")
  ELSE(UNIX)
    STRING(REGEX REPLACE "\\\\" "/" SRepVisualization_DIR_SEARCH1 "$ENV{PATH}")
  ENDIF(UNIX)
  STRING(REGEX REPLACE "/;" ";" SRepVisualization_DIR_SEARCH2 "${SRepVisualization_DIR_SEARCH1}")

  # Construct a set of paths relative to the system search path.
  SET(SRepVisualization_DIR_SEARCH "")
  FOREACH(dir ${SRepVisualization_DIR_SEARCH2})
    SET(SRepVisualization_DIR_SEARCH ${SRepVisualization_DIR_SEARCH}
      ${dir}/../lib/lib64/SRepVisualization
      )
  ENDFOREACH(dir)

  #
  # Look for an installation or build tree.
  #
  FIND_PATH(SRepVisualization_DIR UseSRepVisualization.cmake
    # Look for an environment variable SRepVisualization_DIR.
    $ENV{SRepVisualization_DIR}

    # Look in places relative to the system executable search path.
    ${SRepVisualization_DIR_SEARCH}

    # Look in standard WIN install locations.
    "$ENV{ProgramFiles}/SRepVisualization"

    # Look in standard UNIX install locations.
    /usr/local/lib64/SRepVisualization
    /usr/lib64/SRepVisualization

    # Read from the CMakeSetup registry entries.  It is likely that
    # SRepVisualization will have been recently built.
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
    DOC "The ${SRepVisualization_DIR_DESCRIPTION}"
  )
ENDIF(NOT SRepVisualization_DIR)

# If SRepVisualization was found, load the configuration file to get the rest of the
# settings.
IF(SRepVisualization_DIR)
  # Make sure the SRepVisualizationConfig.cmake file exists in the directory provided.
  IF(EXISTS ${SRepVisualization_DIR}/SRepVisualizationConfig.cmake)

    # We found SRepVisualization.  Load the settings.
    SET(SRepVisualization_FOUND 1)
    INCLUDE(${SRepVisualization_DIR}/SRepVisualizationConfig.cmake)

  ENDIF(EXISTS ${SRepVisualization_DIR}/SRepVisualizationConfig.cmake)
ELSE(SRepVisualization_DIR)
  # We did not find SRepVisualization.
  SET(SRepVisualization_FOUND 0)
ENDIF(SRepVisualization_DIR)

#-----------------------------------------------------------------------------
IF(NOT SRepVisualization_FOUND)
  # SRepVisualization not found, explain to the user how to specify its location.
  IF(NOT SRepVisualization_FIND_QUIETLY)
    MESSAGE(FATAL_ERROR ${SRepVisualization_NOT_FOUND_MESSAGE})
  ELSE(NOT SRepVisualization_FIND_QUIETLY)
    IF(SRepVisualization_FIND_REQUIRED)
      MESSAGE(FATAL_ERROR ${SRepVisualization_NOT_FOUND_MESSAGE})
    ENDIF(SRepVisualization_FIND_REQUIRED)
  ENDIF(NOT SRepVisualization_FIND_QUIETLY)
ENDIF(NOT SRepVisualization_FOUND)
