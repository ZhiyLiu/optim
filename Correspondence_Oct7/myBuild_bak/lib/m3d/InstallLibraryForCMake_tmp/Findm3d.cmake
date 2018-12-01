# - Find a library installation or build tree.
# 
# The following variables are set if m3d is found.  
# If m3d is not found, m3d_FOUND is set to false.
#  m3d_FOUND         - Set to true when m3d is found.
#  m3d_USE_FILE      - CMake file to use m3d.
#  m3d_MAJOR_VERSION - The m3d major version number.
#  m3d_MINOR_VERSION - The m3d minor version number 
#                       (odd non-release).
#  m3d_BUILD_VERSION - The m3d patch level 
#                       (meaningless for odd minor).
#  m3d_INCLUDE_DIRS  - Include directories for m3d
#  m3d_LIBRARY_DIRS  - Link directories for m3d libraries
#  m3d_LIBRARIES     - List of libraries to link against
#
# The following cache entries must be set by the user to locate m3d:
#  m3d_DIR  - The directory containing m3dConfig.cmake.  
#             This is either the root of the build tree,
#             or the lib/m3d directory.  This is the 
#             only cache entry.


# Construct consitent error messages for use below.
SET(m3d_DIR_DESCRIPTION "directory containing m3dConfig.cmake.  This is either the root of the build tree, or PREFIX/lib/m3d for an installation.")
SET(m3d_NOT_FOUND_MESSAGE "m3d not found.  Set the m3d_DIR cmake cache entry to the ${m3d_DIR_DESCRIPTION}")

# Search only if the location is not already known.
IF(NOT m3d_DIR)
  # Get the system search path as a list.
  IF(UNIX)
    STRING(REGEX MATCHALL "[^:]+" m3d_DIR_SEARCH1 "$ENV{PATH}")
  ELSE(UNIX)
    STRING(REGEX REPLACE "\\\\" "/" m3d_DIR_SEARCH1 "$ENV{PATH}")
  ENDIF(UNIX)
  STRING(REGEX REPLACE "/;" ";" m3d_DIR_SEARCH2 "${m3d_DIR_SEARCH1}")

  # Construct a set of paths relative to the system search path.
  SET(m3d_DIR_SEARCH "")
  FOREACH(dir ${m3d_DIR_SEARCH2})
    SET(m3d_DIR_SEARCH ${m3d_DIR_SEARCH}
      ${dir}/../lib/lib64/m3d
      )
  ENDFOREACH(dir)

  #
  # Look for an installation or build tree.
  #
  FIND_PATH(m3d_DIR Usem3d.cmake
    # Look for an environment variable m3d_DIR.
    $ENV{m3d_DIR}

    # Look in places relative to the system executable search path.
    ${m3d_DIR_SEARCH}

    # Look in standard WIN install locations.
    "$ENV{ProgramFiles}/m3d"

    # Look in standard UNIX install locations.
    /usr/local/lib64/m3d
    /usr/lib64/m3d

    # Read from the CMakeSetup registry entries.  It is likely that
    # m3d will have been recently built.
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
    DOC "The ${m3d_DIR_DESCRIPTION}"
  )
ENDIF(NOT m3d_DIR)

# If m3d was found, load the configuration file to get the rest of the
# settings.
IF(m3d_DIR)
  # Make sure the m3dConfig.cmake file exists in the directory provided.
  IF(EXISTS ${m3d_DIR}/m3dConfig.cmake)

    # We found m3d.  Load the settings.
    SET(m3d_FOUND 1)
    INCLUDE(${m3d_DIR}/m3dConfig.cmake)

  ENDIF(EXISTS ${m3d_DIR}/m3dConfig.cmake)
ELSE(m3d_DIR)
  # We did not find m3d.
  SET(m3d_FOUND 0)
ENDIF(m3d_DIR)

#-----------------------------------------------------------------------------
IF(NOT m3d_FOUND)
  # m3d not found, explain to the user how to specify its location.
  IF(NOT m3d_FIND_QUIETLY)
    MESSAGE(FATAL_ERROR ${m3d_NOT_FOUND_MESSAGE})
  ELSE(NOT m3d_FIND_QUIETLY)
    IF(m3d_FIND_REQUIRED)
      MESSAGE(FATAL_ERROR ${m3d_NOT_FOUND_MESSAGE})
    ENDIF(m3d_FIND_REQUIRED)
  ENDIF(NOT m3d_FIND_QUIETLY)
ENDIF(NOT m3d_FOUND)
