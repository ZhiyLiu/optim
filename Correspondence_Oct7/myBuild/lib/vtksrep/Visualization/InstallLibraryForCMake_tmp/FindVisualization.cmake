# - Find a library installation or build tree.
# 
# The following variables are set if Visualization is found.  
# If Visualization is not found, Visualization_FOUND is set to false.
#  Visualization_FOUND         - Set to true when Visualization is found.
#  Visualization_USE_FILE      - CMake file to use Visualization.
#  Visualization_MAJOR_VERSION - The Visualization major version number.
#  Visualization_MINOR_VERSION - The Visualization minor version number 
#                       (odd non-release).
#  Visualization_BUILD_VERSION - The Visualization patch level 
#                       (meaningless for odd minor).
#  Visualization_INCLUDE_DIRS  - Include directories for Visualization
#  Visualization_LIBRARY_DIRS  - Link directories for Visualization libraries
#  Visualization_LIBRARIES     - List of libraries to link against
#
# The following cache entries must be set by the user to locate Visualization:
#  Visualization_DIR  - The directory containing VisualizationConfig.cmake.  
#             This is either the root of the build tree,
#             or the lib/Visualization directory.  This is the 
#             only cache entry.


# Construct consitent error messages for use below.
SET(Visualization_DIR_DESCRIPTION "directory containing VisualizationConfig.cmake.  This is either the root of the build tree, or PREFIX/lib/Visualization for an installation.")
SET(Visualization_NOT_FOUND_MESSAGE "Visualization not found.  Set the Visualization_DIR cmake cache entry to the ${Visualization_DIR_DESCRIPTION}")

# Search only if the location is not already known.
IF(NOT Visualization_DIR)
  # Get the system search path as a list.
  IF(UNIX)
    STRING(REGEX MATCHALL "[^:]+" Visualization_DIR_SEARCH1 "$ENV{PATH}")
  ELSE(UNIX)
    STRING(REGEX REPLACE "\\\\" "/" Visualization_DIR_SEARCH1 "$ENV{PATH}")
  ENDIF(UNIX)
  STRING(REGEX REPLACE "/;" ";" Visualization_DIR_SEARCH2 "${Visualization_DIR_SEARCH1}")

  # Construct a set of paths relative to the system search path.
  SET(Visualization_DIR_SEARCH "")
  FOREACH(dir ${Visualization_DIR_SEARCH2})
    SET(Visualization_DIR_SEARCH ${Visualization_DIR_SEARCH}
      ${dir}/../lib/lib64/Visualization
      )
  ENDFOREACH(dir)

  #
  # Look for an installation or build tree.
  #
  FIND_PATH(Visualization_DIR UseVisualization.cmake
    # Look for an environment variable Visualization_DIR.
    $ENV{Visualization_DIR}

    # Look in places relative to the system executable search path.
    ${Visualization_DIR_SEARCH}

    # Look in standard WIN install locations.
    "$ENV{ProgramFiles}/Visualization"

    # Look in standard UNIX install locations.
    /usr/local/lib64/Visualization
    /usr/lib64/Visualization

    # Read from the CMakeSetup registry entries.  It is likely that
    # Visualization will have been recently built.
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
    DOC "The ${Visualization_DIR_DESCRIPTION}"
  )
ENDIF(NOT Visualization_DIR)

# If Visualization was found, load the configuration file to get the rest of the
# settings.
IF(Visualization_DIR)
  # Make sure the VisualizationConfig.cmake file exists in the directory provided.
  IF(EXISTS ${Visualization_DIR}/VisualizationConfig.cmake)

    # We found Visualization.  Load the settings.
    SET(Visualization_FOUND 1)
    INCLUDE(${Visualization_DIR}/VisualizationConfig.cmake)

  ENDIF(EXISTS ${Visualization_DIR}/VisualizationConfig.cmake)
ELSE(Visualization_DIR)
  # We did not find Visualization.
  SET(Visualization_FOUND 0)
ENDIF(Visualization_DIR)

#-----------------------------------------------------------------------------
IF(NOT Visualization_FOUND)
  # Visualization not found, explain to the user how to specify its location.
  IF(NOT Visualization_FIND_QUIETLY)
    MESSAGE(FATAL_ERROR ${Visualization_NOT_FOUND_MESSAGE})
  ELSE(NOT Visualization_FIND_QUIETLY)
    IF(Visualization_FIND_REQUIRED)
      MESSAGE(FATAL_ERROR ${Visualization_NOT_FOUND_MESSAGE})
    ENDIF(Visualization_FIND_REQUIRED)
  ENDIF(NOT Visualization_FIND_QUIETLY)
ENDIF(NOT Visualization_FOUND)
