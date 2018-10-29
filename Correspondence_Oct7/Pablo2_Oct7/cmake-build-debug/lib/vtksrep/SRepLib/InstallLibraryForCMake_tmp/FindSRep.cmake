# - Find a library installation or build tree.
# 
# The following variables are set if SRep is found.  
# If SRep is not found, SRep_FOUND is set to false.
#  SRep_FOUND         - Set to true when SRep is found.
#  SRep_USE_FILE      - CMake file to use SRep.
#  SRep_MAJOR_VERSION - The SRep major version number.
#  SRep_MINOR_VERSION - The SRep minor version number 
#                       (odd non-release).
#  SRep_BUILD_VERSION - The SRep patch level 
#                       (meaningless for odd minor).
#  SRep_INCLUDE_DIRS  - Include directories for SRep
#  SRep_LIBRARY_DIRS  - Link directories for SRep libraries
#  SRep_LIBRARIES     - List of libraries to link against
#
# The following cache entries must be set by the user to locate SRep:
#  SRep_DIR  - The directory containing SRepConfig.cmake.  
#             This is either the root of the build tree,
#             or the lib/SRep directory.  This is the 
#             only cache entry.


# Construct consitent error messages for use below.
SET(SRep_DIR_DESCRIPTION "directory containing SRepConfig.cmake.  This is either the root of the build tree, or PREFIX/lib/SRep for an installation.")
SET(SRep_NOT_FOUND_MESSAGE "SRep not found.  Set the SRep_DIR cmake cache entry to the ${SRep_DIR_DESCRIPTION}")

# Search only if the location is not already known.
IF(NOT SRep_DIR)
  # Get the system search path as a list.
  IF(UNIX)
    STRING(REGEX MATCHALL "[^:]+" SRep_DIR_SEARCH1 "$ENV{PATH}")
  ELSE(UNIX)
    STRING(REGEX REPLACE "\\\\" "/" SRep_DIR_SEARCH1 "$ENV{PATH}")
  ENDIF(UNIX)
  STRING(REGEX REPLACE "/;" ";" SRep_DIR_SEARCH2 "${SRep_DIR_SEARCH1}")

  # Construct a set of paths relative to the system search path.
  SET(SRep_DIR_SEARCH "")
  FOREACH(dir ${SRep_DIR_SEARCH2})
    SET(SRep_DIR_SEARCH ${SRep_DIR_SEARCH}
      ${dir}/../lib/lib64/SRep
      )
  ENDFOREACH(dir)

  #
  # Look for an installation or build tree.
  #
  FIND_PATH(SRep_DIR UseSRep.cmake
    # Look for an environment variable SRep_DIR.
    $ENV{SRep_DIR}

    # Look in places relative to the system executable search path.
    ${SRep_DIR_SEARCH}

    # Look in standard WIN install locations.
    "$ENV{ProgramFiles}/SRep"

    # Look in standard UNIX install locations.
    /usr/local/lib64/SRep
    /usr/lib64/SRep

    # Read from the CMakeSetup registry entries.  It is likely that
    # SRep will have been recently built.
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
    DOC "The ${SRep_DIR_DESCRIPTION}"
  )
ENDIF(NOT SRep_DIR)

# If SRep was found, load the configuration file to get the rest of the
# settings.
IF(SRep_DIR)
  # Make sure the SRepConfig.cmake file exists in the directory provided.
  IF(EXISTS ${SRep_DIR}/SRepConfig.cmake)

    # We found SRep.  Load the settings.
    SET(SRep_FOUND 1)
    INCLUDE(${SRep_DIR}/SRepConfig.cmake)

  ENDIF(EXISTS ${SRep_DIR}/SRepConfig.cmake)
ELSE(SRep_DIR)
  # We did not find SRep.
  SET(SRep_FOUND 0)
ENDIF(SRep_DIR)

#-----------------------------------------------------------------------------
IF(NOT SRep_FOUND)
  # SRep not found, explain to the user how to specify its location.
  IF(NOT SRep_FIND_QUIETLY)
    MESSAGE(FATAL_ERROR ${SRep_NOT_FOUND_MESSAGE})
  ELSE(NOT SRep_FIND_QUIETLY)
    IF(SRep_FIND_REQUIRED)
      MESSAGE(FATAL_ERROR ${SRep_NOT_FOUND_MESSAGE})
    ENDIF(SRep_FIND_REQUIRED)
  ENDIF(NOT SRep_FIND_QUIETLY)
ENDIF(NOT SRep_FOUND)
