# - Find a library installation or build tree.
# 
# The following variables are set if seurat is found.  
# If seurat is not found, seurat_FOUND is set to false.
#  seurat_FOUND         - Set to true when seurat is found.
#  seurat_USE_FILE      - CMake file to use seurat.
#  seurat_MAJOR_VERSION - The seurat major version number.
#  seurat_MINOR_VERSION - The seurat minor version number 
#                       (odd non-release).
#  seurat_BUILD_VERSION - The seurat patch level 
#                       (meaningless for odd minor).
#  seurat_INCLUDE_DIRS  - Include directories for seurat
#  seurat_LIBRARY_DIRS  - Link directories for seurat libraries
#  seurat_LIBRARIES     - List of libraries to link against
#
# The following cache entries must be set by the user to locate seurat:
#  seurat_DIR  - The directory containing seuratConfig.cmake.  
#             This is either the root of the build tree,
#             or the lib/seurat directory.  This is the 
#             only cache entry.


# Construct consitent error messages for use below.
SET(seurat_DIR_DESCRIPTION "directory containing seuratConfig.cmake.  This is either the root of the build tree, or PREFIX/lib/seurat for an installation.")
SET(seurat_NOT_FOUND_MESSAGE "seurat not found.  Set the seurat_DIR cmake cache entry to the ${seurat_DIR_DESCRIPTION}")

# Search only if the location is not already known.
IF(NOT seurat_DIR)
  # Get the system search path as a list.
  IF(UNIX)
    STRING(REGEX MATCHALL "[^:]+" seurat_DIR_SEARCH1 "$ENV{PATH}")
  ELSE(UNIX)
    STRING(REGEX REPLACE "\\\\" "/" seurat_DIR_SEARCH1 "$ENV{PATH}")
  ENDIF(UNIX)
  STRING(REGEX REPLACE "/;" ";" seurat_DIR_SEARCH2 "${seurat_DIR_SEARCH1}")

  # Construct a set of paths relative to the system search path.
  SET(seurat_DIR_SEARCH "")
  FOREACH(dir ${seurat_DIR_SEARCH2})
    SET(seurat_DIR_SEARCH ${seurat_DIR_SEARCH}
      ${dir}/../lib/lib64/seurat
      )
  ENDFOREACH(dir)

  #
  # Look for an installation or build tree.
  #
  FIND_PATH(seurat_DIR Useseurat.cmake
    # Look for an environment variable seurat_DIR.
    $ENV{seurat_DIR}

    # Look in places relative to the system executable search path.
    ${seurat_DIR_SEARCH}

    # Look in standard WIN install locations.
    "$ENV{ProgramFiles}/seurat"

    # Look in standard UNIX install locations.
    /usr/local/lib64/seurat
    /usr/lib64/seurat

    # Read from the CMakeSetup registry entries.  It is likely that
    # seurat will have been recently built.
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
    DOC "The ${seurat_DIR_DESCRIPTION}"
  )
ENDIF(NOT seurat_DIR)

# If seurat was found, load the configuration file to get the rest of the
# settings.
IF(seurat_DIR)
  # Make sure the seuratConfig.cmake file exists in the directory provided.
  IF(EXISTS ${seurat_DIR}/seuratConfig.cmake)

    # We found seurat.  Load the settings.
    SET(seurat_FOUND 1)
    INCLUDE(${seurat_DIR}/seuratConfig.cmake)

  ENDIF(EXISTS ${seurat_DIR}/seuratConfig.cmake)
ELSE(seurat_DIR)
  # We did not find seurat.
  SET(seurat_FOUND 0)
ENDIF(seurat_DIR)

#-----------------------------------------------------------------------------
IF(NOT seurat_FOUND)
  # seurat not found, explain to the user how to specify its location.
  IF(NOT seurat_FIND_QUIETLY)
    MESSAGE(FATAL_ERROR ${seurat_NOT_FOUND_MESSAGE})
  ELSE(NOT seurat_FIND_QUIETLY)
    IF(seurat_FIND_REQUIRED)
      MESSAGE(FATAL_ERROR ${seurat_NOT_FOUND_MESSAGE})
    ENDIF(seurat_FIND_REQUIRED)
  ENDIF(NOT seurat_FIND_QUIETLY)
ENDIF(NOT seurat_FOUND)
