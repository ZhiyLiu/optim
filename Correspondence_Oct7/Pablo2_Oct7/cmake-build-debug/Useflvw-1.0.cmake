# This is an implementation detail for using flvw-1.0 with the
# Findflvw-1.0.cmake module.  Do not include directly by name.  
# This should be included only when Findflvw-1.0.cmake sets 
# the flvw-1.0_USE_FILE variable to point here.

IF(CREA_VERBOSE_CMAKE)
  MESSAGE(STATUS "Using flvw-1.0")
ENDIF(CREA_VERBOSE_CMAKE)

# Load the compiler settings used for flvw-1.0.
IF(flvw-1.0_BUILD_SETTINGS_FILE)
  INCLUDE(CMakeImportBuildSettings)
  CMAKE_IMPORT_BUILD_SETTINGS(${flvw-1.0_BUILD_SETTINGS_FILE})
ENDIF(flvw-1.0_BUILD_SETTINGS_FILE)

# Add compiler flags needed to use flvw-1.0.
SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${flvw-1.0_REQUIRED_C_FLAGS}")
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${flvw-1.0_REQUIRED_CXX_FLAGS}")
SET(CMAKE_LINK_FLAGS "${CMAKE_LINK_FLAGS} ${flvw-1.0_REQUIRED_LINK_FLAGS}")

# Add include directories needed to use flvw-1.0.
INCLUDE_DIRECTORIES(${flvw-1.0_INCLUDE_DIRS})

# Add link directories needed to use flvw-1.0.
LINK_DIRECTORIES(${flvw-1.0_LIBRARY_DIRS})

# Set the version 
# Already done in bbtkConfigure.h
#ADD_DEFINITIONS( -Dflvw-1.0_VERSION="\"${flvw-1.0_VERSION}\"" )

# Additional use file 
IF (flvw-1.0_HAS_ADDITIONAL_CONFIG_FILE)
  # Include it
  INCLUDE(${flvw-1.0_DIR}/AdditionalUseflvw-1.0.cmake)
ENDIF (flvw-1.0_HAS_ADDITIONAL_CONFIG_FILE)
