# This is an implementation detail for using SRepVisualization with the
# FindSRepVisualization.cmake module.  Do not include directly by name.  
# This should be included only when FindSRepVisualization.cmake sets 
# the SRepVisualization_USE_FILE variable to point here.

IF(CREA_VERBOSE_CMAKE)
  MESSAGE(STATUS "Using SRepVisualization")
ENDIF(CREA_VERBOSE_CMAKE)

# Load the compiler settings used for SRepVisualization.
IF(SRepVisualization_BUILD_SETTINGS_FILE)
  INCLUDE(CMakeImportBuildSettings)
  CMAKE_IMPORT_BUILD_SETTINGS(${SRepVisualization_BUILD_SETTINGS_FILE})
ENDIF(SRepVisualization_BUILD_SETTINGS_FILE)

# Add compiler flags needed to use SRepVisualization.
SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${SRepVisualization_REQUIRED_C_FLAGS}")
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SRepVisualization_REQUIRED_CXX_FLAGS}")
SET(CMAKE_LINK_FLAGS "${CMAKE_LINK_FLAGS} ${SRepVisualization_REQUIRED_LINK_FLAGS}")

# Add include directories needed to use SRepVisualization.
INCLUDE_DIRECTORIES(${SRepVisualization_INCLUDE_DIRS})

# Add link directories needed to use SRepVisualization.
LINK_DIRECTORIES(${SRepVisualization_LIBRARY_DIRS})

# Set the version 
# Already done in bbtkConfigure.h
#ADD_DEFINITIONS( -DSRepVisualization_VERSION="\"${SRepVisualization_VERSION}\"" )

# Additional use file 
IF (SRepVisualization_HAS_ADDITIONAL_CONFIG_FILE)
  # Include it
  INCLUDE(${SRepVisualization_DIR}/AdditionalUseSRepVisualization.cmake)
ENDIF (SRepVisualization_HAS_ADDITIONAL_CONFIG_FILE)
