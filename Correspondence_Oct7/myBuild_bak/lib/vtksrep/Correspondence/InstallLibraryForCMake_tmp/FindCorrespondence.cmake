# - Find a library installation or build tree.
# 
# The following variables are set if Correspondence is found.  
# If Correspondence is not found, Correspondence_FOUND is set to false.
#  Correspondence_FOUND         - Set to true when Correspondence is found.
#  Correspondence_USE_FILE      - CMake file to use Correspondence.
#  Correspondence_MAJOR_VERSION - The Correspondence major version number.
#  Correspondence_MINOR_VERSION - The Correspondence minor version number 
#                       (odd non-release).
#  Correspondence_BUILD_VERSION - The Correspondence patch level 
#                       (meaningless for odd minor).
#  Correspondence_INCLUDE_DIRS  - Include directories for Correspondence
#  Correspondence_LIBRARY_DIRS  - Link directories for Correspondence libraries
#  Correspondence_LIBRARIES     - List of libraries to link against
#
# The following cache entries must be set by the user to locate Correspondence:
#  Correspondence_DIR  - The directory containing CorrespondenceConfig.cmake.  
#             This is either the root of the build tree,
#             or the lib/Correspondence directory.  This is the 
#             only cache entry.


# Construct consitent error messages for use below.
SET(Correspondence_DIR_DESCRIPTION "directory containing CorrespondenceConfig.cmake.  This is either the root of the build tree, or PREFIX/lib/Correspondence for an installation.")
SET(Correspondence_NOT_FOUND_MESSAGE "Correspondence not found.  Set the Correspondence_DIR cmake cache entry to the ${Correspondence_DIR_DESCRIPTION}")

# Search only if the location is not already known.
IF(NOT Correspondence_DIR)
  # Get the system search path as a list.
  IF(UNIX)
    STRING(REGEX MATCHALL "[^:]+" Correspondence_DIR_SEARCH1 "$ENV{PATH}")
  ELSE(UNIX)
    STRING(REGEX REPLACE "\\\\" "/" Correspondence_DIR_SEARCH1 "$ENV{PATH}")
  ENDIF(UNIX)
  STRING(REGEX REPLACE "/;" ";" Correspondence_DIR_SEARCH2 "${Correspondence_DIR_SEARCH1}")

  # Construct a set of paths relative to the system search path.
  SET(Correspondence_DIR_SEARCH "")
  FOREACH(dir ${Correspondence_DIR_SEARCH2})
    SET(Correspondence_DIR_SEARCH ${Correspondence_DIR_SEARCH}
      ${dir}/../lib/lib64/Correspondence
      )
  ENDFOREACH(dir)

  #
  # Look for an installation or build tree.
  #
  FIND_PATH(Correspondence_DIR UseCorrespondence.cmake
    # Look for an environment variable Correspondence_DIR.
    $ENV{Correspondence_DIR}

    # Look in places relative to the system executable search path.
    ${Correspondence_DIR_SEARCH}

    # Look in standard WIN install locations.
    "$ENV{ProgramFiles}/Correspondence"

    # Look in standard UNIX install locations.
    /usr/local/lib64/Correspondence
    /usr/lib64/Correspondence

    # Read from the CMakeSetup registry entries.  It is likely that
    # Correspondence will have been recently built.
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
    DOC "The ${Correspondence_DIR_DESCRIPTION}"
  )
ENDIF(NOT Correspondence_DIR)

# If Correspondence was found, load the configuration file to get the rest of the
# settings.
IF(Correspondence_DIR)
  # Make sure the CorrespondenceConfig.cmake file exists in the directory provided.
  IF(EXISTS ${Correspondence_DIR}/CorrespondenceConfig.cmake)

    # We found Correspondence.  Load the settings.
    SET(Correspondence_FOUND 1)
    INCLUDE(${Correspondence_DIR}/CorrespondenceConfig.cmake)

  ENDIF(EXISTS ${Correspondence_DIR}/CorrespondenceConfig.cmake)
ELSE(Correspondence_DIR)
  # We did not find Correspondence.
  SET(Correspondence_FOUND 0)
ENDIF(Correspondence_DIR)

#-----------------------------------------------------------------------------
IF(NOT Correspondence_FOUND)
  # Correspondence not found, explain to the user how to specify its location.
  IF(NOT Correspondence_FIND_QUIETLY)
    MESSAGE(FATAL_ERROR ${Correspondence_NOT_FOUND_MESSAGE})
  ELSE(NOT Correspondence_FIND_QUIETLY)
    IF(Correspondence_FIND_REQUIRED)
      MESSAGE(FATAL_ERROR ${Correspondence_NOT_FOUND_MESSAGE})
    ENDIF(Correspondence_FIND_REQUIRED)
  ENDIF(NOT Correspondence_FIND_QUIETLY)
ENDIF(NOT Correspondence_FOUND)
