# - Find a library installation or build tree.
# 
# The following variables are set if gui is found.  
# If gui is not found, gui_FOUND is set to false.
#  gui_FOUND         - Set to true when gui is found.
#  gui_USE_FILE      - CMake file to use gui.
#  gui_MAJOR_VERSION - The gui major version number.
#  gui_MINOR_VERSION - The gui minor version number 
#                       (odd non-release).
#  gui_BUILD_VERSION - The gui patch level 
#                       (meaningless for odd minor).
#  gui_INCLUDE_DIRS  - Include directories for gui
#  gui_LIBRARY_DIRS  - Link directories for gui libraries
#  gui_LIBRARIES     - List of libraries to link against
#
# The following cache entries must be set by the user to locate gui:
#  gui_DIR  - The directory containing guiConfig.cmake.  
#             This is either the root of the build tree,
#             or the lib/gui directory.  This is the 
#             only cache entry.


# Construct consitent error messages for use below.
SET(gui_DIR_DESCRIPTION "directory containing guiConfig.cmake.  This is either the root of the build tree, or PREFIX/lib/gui for an installation.")
SET(gui_NOT_FOUND_MESSAGE "gui not found.  Set the gui_DIR cmake cache entry to the ${gui_DIR_DESCRIPTION}")

# Search only if the location is not already known.
IF(NOT gui_DIR)
  # Get the system search path as a list.
  IF(UNIX)
    STRING(REGEX MATCHALL "[^:]+" gui_DIR_SEARCH1 "$ENV{PATH}")
  ELSE(UNIX)
    STRING(REGEX REPLACE "\\\\" "/" gui_DIR_SEARCH1 "$ENV{PATH}")
  ENDIF(UNIX)
  STRING(REGEX REPLACE "/;" ";" gui_DIR_SEARCH2 "${gui_DIR_SEARCH1}")

  # Construct a set of paths relative to the system search path.
  SET(gui_DIR_SEARCH "")
  FOREACH(dir ${gui_DIR_SEARCH2})
    SET(gui_DIR_SEARCH ${gui_DIR_SEARCH}
      ${dir}/../lib/lib64/gui
      )
  ENDFOREACH(dir)

  #
  # Look for an installation or build tree.
  #
  FIND_PATH(gui_DIR Usegui.cmake
    # Look for an environment variable gui_DIR.
    $ENV{gui_DIR}

    # Look in places relative to the system executable search path.
    ${gui_DIR_SEARCH}

    # Look in standard WIN install locations.
    "$ENV{ProgramFiles}/gui"

    # Look in standard UNIX install locations.
    /usr/local/lib64/gui
    /usr/lib64/gui

    # Read from the CMakeSetup registry entries.  It is likely that
    # gui will have been recently built.
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
    DOC "The ${gui_DIR_DESCRIPTION}"
  )
ENDIF(NOT gui_DIR)

# If gui was found, load the configuration file to get the rest of the
# settings.
IF(gui_DIR)
  # Make sure the guiConfig.cmake file exists in the directory provided.
  IF(EXISTS ${gui_DIR}/guiConfig.cmake)

    # We found gui.  Load the settings.
    SET(gui_FOUND 1)
    INCLUDE(${gui_DIR}/guiConfig.cmake)

  ENDIF(EXISTS ${gui_DIR}/guiConfig.cmake)
ELSE(gui_DIR)
  # We did not find gui.
  SET(gui_FOUND 0)
ENDIF(gui_DIR)

#-----------------------------------------------------------------------------
IF(NOT gui_FOUND)
  # gui not found, explain to the user how to specify its location.
  IF(NOT gui_FIND_QUIETLY)
    MESSAGE(FATAL_ERROR ${gui_NOT_FOUND_MESSAGE})
  ELSE(NOT gui_FIND_QUIETLY)
    IF(gui_FIND_REQUIRED)
      MESSAGE(FATAL_ERROR ${gui_NOT_FOUND_MESSAGE})
    ENDIF(gui_FIND_REQUIRED)
  ENDIF(NOT gui_FIND_QUIETLY)
ENDIF(NOT gui_FOUND)
