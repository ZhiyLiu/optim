# - Find a library installation or build tree.
# 
# The following variables are set if vtkpowercrust is found.  
# If vtkpowercrust is not found, vtkpowercrust_FOUND is set to false.
#  vtkpowercrust_FOUND         - Set to true when vtkpowercrust is found.
#  vtkpowercrust_USE_FILE      - CMake file to use vtkpowercrust.
#  vtkpowercrust_MAJOR_VERSION - The vtkpowercrust major version number.
#  vtkpowercrust_MINOR_VERSION - The vtkpowercrust minor version number 
#                       (odd non-release).
#  vtkpowercrust_BUILD_VERSION - The vtkpowercrust patch level 
#                       (meaningless for odd minor).
#  vtkpowercrust_INCLUDE_DIRS  - Include directories for vtkpowercrust
#  vtkpowercrust_LIBRARY_DIRS  - Link directories for vtkpowercrust libraries
#  vtkpowercrust_LIBRARIES     - List of libraries to link against
#
# The following cache entries must be set by the user to locate vtkpowercrust:
#  vtkpowercrust_DIR  - The directory containing vtkpowercrustConfig.cmake.  
#             This is either the root of the build tree,
#             or the lib/vtkpowercrust directory.  This is the 
#             only cache entry.


# Construct consitent error messages for use below.
SET(vtkpowercrust_DIR_DESCRIPTION "directory containing vtkpowercrustConfig.cmake.  This is either the root of the build tree, or PREFIX/lib/vtkpowercrust for an installation.")
SET(vtkpowercrust_NOT_FOUND_MESSAGE "vtkpowercrust not found.  Set the vtkpowercrust_DIR cmake cache entry to the ${vtkpowercrust_DIR_DESCRIPTION}")

# Search only if the location is not already known.
IF(NOT vtkpowercrust_DIR)
  # Get the system search path as a list.
  IF(UNIX)
    STRING(REGEX MATCHALL "[^:]+" vtkpowercrust_DIR_SEARCH1 "$ENV{PATH}")
  ELSE(UNIX)
    STRING(REGEX REPLACE "\\\\" "/" vtkpowercrust_DIR_SEARCH1 "$ENV{PATH}")
  ENDIF(UNIX)
  STRING(REGEX REPLACE "/;" ";" vtkpowercrust_DIR_SEARCH2 "${vtkpowercrust_DIR_SEARCH1}")

  # Construct a set of paths relative to the system search path.
  SET(vtkpowercrust_DIR_SEARCH "")
  FOREACH(dir ${vtkpowercrust_DIR_SEARCH2})
    SET(vtkpowercrust_DIR_SEARCH ${vtkpowercrust_DIR_SEARCH}
      ${dir}/../lib/lib64/vtkpowercrust
      )
  ENDFOREACH(dir)

  #
  # Look for an installation or build tree.
  #
  FIND_PATH(vtkpowercrust_DIR Usevtkpowercrust.cmake
    # Look for an environment variable vtkpowercrust_DIR.
    $ENV{vtkpowercrust_DIR}

    # Look in places relative to the system executable search path.
    ${vtkpowercrust_DIR_SEARCH}

    # Look in standard WIN install locations.
    "$ENV{ProgramFiles}/vtkpowercrust"

    # Look in standard UNIX install locations.
    /usr/local/lib64/vtkpowercrust
    /usr/lib64/vtkpowercrust

    # Read from the CMakeSetup registry entries.  It is likely that
    # vtkpowercrust will have been recently built.
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
    DOC "The ${vtkpowercrust_DIR_DESCRIPTION}"
  )
ENDIF(NOT vtkpowercrust_DIR)

# If vtkpowercrust was found, load the configuration file to get the rest of the
# settings.
IF(vtkpowercrust_DIR)
  # Make sure the vtkpowercrustConfig.cmake file exists in the directory provided.
  IF(EXISTS ${vtkpowercrust_DIR}/vtkpowercrustConfig.cmake)

    # We found vtkpowercrust.  Load the settings.
    SET(vtkpowercrust_FOUND 1)
    INCLUDE(${vtkpowercrust_DIR}/vtkpowercrustConfig.cmake)

  ENDIF(EXISTS ${vtkpowercrust_DIR}/vtkpowercrustConfig.cmake)
ELSE(vtkpowercrust_DIR)
  # We did not find vtkpowercrust.
  SET(vtkpowercrust_FOUND 0)
ENDIF(vtkpowercrust_DIR)

#-----------------------------------------------------------------------------
IF(NOT vtkpowercrust_FOUND)
  # vtkpowercrust not found, explain to the user how to specify its location.
  IF(NOT vtkpowercrust_FIND_QUIETLY)
    MESSAGE(FATAL_ERROR ${vtkpowercrust_NOT_FOUND_MESSAGE})
  ELSE(NOT vtkpowercrust_FIND_QUIETLY)
    IF(vtkpowercrust_FIND_REQUIRED)
      MESSAGE(FATAL_ERROR ${vtkpowercrust_NOT_FOUND_MESSAGE})
    ENDIF(vtkpowercrust_FIND_REQUIRED)
  ENDIF(NOT vtkpowercrust_FIND_QUIETLY)
ENDIF(NOT vtkpowercrust_FOUND)
