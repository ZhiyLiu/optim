# - Find a library installation or build tree.
# 
# The following variables are set if ImageIO is found.  
# If ImageIO is not found, ImageIO_FOUND is set to false.
#  ImageIO_FOUND         - Set to true when ImageIO is found.
#  ImageIO_USE_FILE      - CMake file to use ImageIO.
#  ImageIO_MAJOR_VERSION - The ImageIO major version number.
#  ImageIO_MINOR_VERSION - The ImageIO minor version number 
#                       (odd non-release).
#  ImageIO_BUILD_VERSION - The ImageIO patch level 
#                       (meaningless for odd minor).
#  ImageIO_INCLUDE_DIRS  - Include directories for ImageIO
#  ImageIO_LIBRARY_DIRS  - Link directories for ImageIO libraries
#  ImageIO_LIBRARIES     - List of libraries to link against
#
# The following cache entries must be set by the user to locate ImageIO:
#  ImageIO_DIR  - The directory containing ImageIOConfig.cmake.  
#             This is either the root of the build tree,
#             or the lib/ImageIO directory.  This is the 
#             only cache entry.


# Construct consitent error messages for use below.
SET(ImageIO_DIR_DESCRIPTION "directory containing ImageIOConfig.cmake.  This is either the root of the build tree, or PREFIX/lib/ImageIO for an installation.")
SET(ImageIO_NOT_FOUND_MESSAGE "ImageIO not found.  Set the ImageIO_DIR cmake cache entry to the ${ImageIO_DIR_DESCRIPTION}")

# Search only if the location is not already known.
IF(NOT ImageIO_DIR)
  # Get the system search path as a list.
  IF(UNIX)
    STRING(REGEX MATCHALL "[^:]+" ImageIO_DIR_SEARCH1 "$ENV{PATH}")
  ELSE(UNIX)
    STRING(REGEX REPLACE "\\\\" "/" ImageIO_DIR_SEARCH1 "$ENV{PATH}")
  ENDIF(UNIX)
  STRING(REGEX REPLACE "/;" ";" ImageIO_DIR_SEARCH2 "${ImageIO_DIR_SEARCH1}")

  # Construct a set of paths relative to the system search path.
  SET(ImageIO_DIR_SEARCH "")
  FOREACH(dir ${ImageIO_DIR_SEARCH2})
    SET(ImageIO_DIR_SEARCH ${ImageIO_DIR_SEARCH}
      ${dir}/../lib/lib64/ImageIO
      )
  ENDFOREACH(dir)

  #
  # Look for an installation or build tree.
  #
  FIND_PATH(ImageIO_DIR UseImageIO.cmake
    # Look for an environment variable ImageIO_DIR.
    $ENV{ImageIO_DIR}

    # Look in places relative to the system executable search path.
    ${ImageIO_DIR_SEARCH}

    # Look in standard WIN install locations.
    "$ENV{ProgramFiles}/ImageIO"

    # Look in standard UNIX install locations.
    /usr/local/lib64/ImageIO
    /usr/lib64/ImageIO

    # Read from the CMakeSetup registry entries.  It is likely that
    # ImageIO will have been recently built.
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
    DOC "The ${ImageIO_DIR_DESCRIPTION}"
  )
ENDIF(NOT ImageIO_DIR)

# If ImageIO was found, load the configuration file to get the rest of the
# settings.
IF(ImageIO_DIR)
  # Make sure the ImageIOConfig.cmake file exists in the directory provided.
  IF(EXISTS ${ImageIO_DIR}/ImageIOConfig.cmake)

    # We found ImageIO.  Load the settings.
    SET(ImageIO_FOUND 1)
    INCLUDE(${ImageIO_DIR}/ImageIOConfig.cmake)

  ENDIF(EXISTS ${ImageIO_DIR}/ImageIOConfig.cmake)
ELSE(ImageIO_DIR)
  # We did not find ImageIO.
  SET(ImageIO_FOUND 0)
ENDIF(ImageIO_DIR)

#-----------------------------------------------------------------------------
IF(NOT ImageIO_FOUND)
  # ImageIO not found, explain to the user how to specify its location.
  IF(NOT ImageIO_FIND_QUIETLY)
    MESSAGE(FATAL_ERROR ${ImageIO_NOT_FOUND_MESSAGE})
  ELSE(NOT ImageIO_FIND_QUIETLY)
    IF(ImageIO_FIND_REQUIRED)
      MESSAGE(FATAL_ERROR ${ImageIO_NOT_FOUND_MESSAGE})
    ENDIF(ImageIO_FIND_REQUIRED)
  ENDIF(NOT ImageIO_FIND_QUIETLY)
ENDIF(NOT ImageIO_FOUND)
