# This is an implementation detail for using gui with the
# Findgui.cmake module.  Do not include directly by name.  
# This should be included only when Findgui.cmake sets 
# the gui_USE_FILE variable to point here.

IF(CREA_VERBOSE_CMAKE)
  MESSAGE(STATUS "Using gui")
ENDIF(CREA_VERBOSE_CMAKE)

# Load the compiler settings used for gui.
IF(gui_BUILD_SETTINGS_FILE)
  INCLUDE(CMakeImportBuildSettings)
  CMAKE_IMPORT_BUILD_SETTINGS(${gui_BUILD_SETTINGS_FILE})
ENDIF(gui_BUILD_SETTINGS_FILE)

# Add compiler flags needed to use gui.
SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${gui_REQUIRED_C_FLAGS}")
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${gui_REQUIRED_CXX_FLAGS}")
SET(CMAKE_LINK_FLAGS "${CMAKE_LINK_FLAGS} ${gui_REQUIRED_LINK_FLAGS}")

# Add include directories needed to use gui.
INCLUDE_DIRECTORIES(${gui_INCLUDE_DIRS})

# Add link directories needed to use gui.
LINK_DIRECTORIES(${gui_LIBRARY_DIRS})

# Set the version 
# Already done in bbtkConfigure.h
#ADD_DEFINITIONS( -Dgui_VERSION="\"${gui_VERSION}\"" )

# Additional use file 
IF (gui_HAS_ADDITIONAL_CONFIG_FILE)
  # Include it
  INCLUDE(${gui_DIR}/AdditionalUsegui.cmake)
ENDIF (gui_HAS_ADDITIONAL_CONFIG_FILE)
