# This is an implementation detail for using planes with the
# Findplanes.cmake module.  Do not include directly by name.  
# This should be included only when Findplanes.cmake sets 
# the planes_USE_FILE variable to point here.

IF(CREA_VERBOSE_CMAKE)
  MESSAGE(STATUS "Using planes")
ENDIF(CREA_VERBOSE_CMAKE)

# Load the compiler settings used for planes.
IF(planes_BUILD_SETTINGS_FILE)
  INCLUDE(CMakeImportBuildSettings)
  CMAKE_IMPORT_BUILD_SETTINGS(${planes_BUILD_SETTINGS_FILE})
ENDIF(planes_BUILD_SETTINGS_FILE)

# Add compiler flags needed to use planes.
SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${planes_REQUIRED_C_FLAGS}")
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${planes_REQUIRED_CXX_FLAGS}")
SET(CMAKE_LINK_FLAGS "${CMAKE_LINK_FLAGS} ${planes_REQUIRED_LINK_FLAGS}")

# Add include directories needed to use planes.
INCLUDE_DIRECTORIES(${planes_INCLUDE_DIRS})

# Add link directories needed to use planes.
LINK_DIRECTORIES(${planes_LIBRARY_DIRS})

# Set the version 
# Already done in bbtkConfigure.h
#ADD_DEFINITIONS( -Dplanes_VERSION="\"${planes_VERSION}\"" )

# Additional use file 
IF (planes_HAS_ADDITIONAL_CONFIG_FILE)
  # Include it
  INCLUDE(${planes_DIR}/AdditionalUseplanes.cmake)
ENDIF (planes_HAS_ADDITIONAL_CONFIG_FILE)
