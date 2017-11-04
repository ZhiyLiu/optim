# This is an implementation detail for using Visualization with the
# FindVisualization.cmake module.  Do not include directly by name.  
# This should be included only when FindVisualization.cmake sets 
# the Visualization_USE_FILE variable to point here.

IF(CREA_VERBOSE_CMAKE)
  MESSAGE(STATUS "Using Visualization")
ENDIF(CREA_VERBOSE_CMAKE)

# Load the compiler settings used for Visualization.
IF(Visualization_BUILD_SETTINGS_FILE)
  INCLUDE(CMakeImportBuildSettings)
  CMAKE_IMPORT_BUILD_SETTINGS(${Visualization_BUILD_SETTINGS_FILE})
ENDIF(Visualization_BUILD_SETTINGS_FILE)

# Add compiler flags needed to use Visualization.
SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${Visualization_REQUIRED_C_FLAGS}")
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${Visualization_REQUIRED_CXX_FLAGS}")
SET(CMAKE_LINK_FLAGS "${CMAKE_LINK_FLAGS} ${Visualization_REQUIRED_LINK_FLAGS}")

# Add include directories needed to use Visualization.
INCLUDE_DIRECTORIES(${Visualization_INCLUDE_DIRS})

# Add link directories needed to use Visualization.
LINK_DIRECTORIES(${Visualization_LIBRARY_DIRS})

# Set the version 
# Already done in bbtkConfigure.h
#ADD_DEFINITIONS( -DVisualization_VERSION="\"${Visualization_VERSION}\"" )

# Additional use file 
IF (Visualization_HAS_ADDITIONAL_CONFIG_FILE)
  # Include it
  INCLUDE(${Visualization_DIR}/AdditionalUseVisualization.cmake)
ENDIF (Visualization_HAS_ADDITIONAL_CONFIG_FILE)
