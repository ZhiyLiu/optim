# This is an implementation detail for using SRepInterpolation with the
# FindSRepInterpolation.cmake module.  Do not include directly by name.  
# This should be included only when FindSRepInterpolation.cmake sets 
# the SRepInterpolation_USE_FILE variable to point here.

IF(CREA_VERBOSE_CMAKE)
  MESSAGE(STATUS "Using SRepInterpolation")
ENDIF(CREA_VERBOSE_CMAKE)

# Load the compiler settings used for SRepInterpolation.
IF(SRepInterpolation_BUILD_SETTINGS_FILE)
  INCLUDE(CMakeImportBuildSettings)
  CMAKE_IMPORT_BUILD_SETTINGS(${SRepInterpolation_BUILD_SETTINGS_FILE})
ENDIF(SRepInterpolation_BUILD_SETTINGS_FILE)

# Add compiler flags needed to use SRepInterpolation.
SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${SRepInterpolation_REQUIRED_C_FLAGS}")
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SRepInterpolation_REQUIRED_CXX_FLAGS}")
SET(CMAKE_LINK_FLAGS "${CMAKE_LINK_FLAGS} ${SRepInterpolation_REQUIRED_LINK_FLAGS}")

# Add include directories needed to use SRepInterpolation.
INCLUDE_DIRECTORIES(${SRepInterpolation_INCLUDE_DIRS})

# Add link directories needed to use SRepInterpolation.
LINK_DIRECTORIES(${SRepInterpolation_LIBRARY_DIRS})

# Set the version 
# Already done in bbtkConfigure.h
#ADD_DEFINITIONS( -DSRepInterpolation_VERSION="\"${SRepInterpolation_VERSION}\"" )

# Additional use file 
IF (SRepInterpolation_HAS_ADDITIONAL_CONFIG_FILE)
  # Include it
  INCLUDE(${SRepInterpolation_DIR}/AdditionalUseSRepInterpolation.cmake)
ENDIF (SRepInterpolation_HAS_ADDITIONAL_CONFIG_FILE)
