# This is an implementation detail for using m3d with the
# Findm3d.cmake module.  Do not include directly by name.  
# This should be included only when Findm3d.cmake sets 
# the m3d_USE_FILE variable to point here.

IF(CREA_VERBOSE_CMAKE)
  MESSAGE(STATUS "Using m3d")
ENDIF(CREA_VERBOSE_CMAKE)

# Load the compiler settings used for m3d.
IF(m3d_BUILD_SETTINGS_FILE)
  INCLUDE(CMakeImportBuildSettings)
  CMAKE_IMPORT_BUILD_SETTINGS(${m3d_BUILD_SETTINGS_FILE})
ENDIF(m3d_BUILD_SETTINGS_FILE)

# Add compiler flags needed to use m3d.
SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${m3d_REQUIRED_C_FLAGS}")
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${m3d_REQUIRED_CXX_FLAGS}")
SET(CMAKE_LINK_FLAGS "${CMAKE_LINK_FLAGS} ${m3d_REQUIRED_LINK_FLAGS}")

# Add include directories needed to use m3d.
INCLUDE_DIRECTORIES(${m3d_INCLUDE_DIRS})

# Add link directories needed to use m3d.
LINK_DIRECTORIES(${m3d_LIBRARY_DIRS})

# Set the version 
# Already done in bbtkConfigure.h
#ADD_DEFINITIONS( -Dm3d_VERSION="\"${m3d_VERSION}\"" )

# Additional use file 
IF (m3d_HAS_ADDITIONAL_CONFIG_FILE)
  # Include it
  INCLUDE(${m3d_DIR}/AdditionalUsem3d.cmake)
ENDIF (m3d_HAS_ADDITIONAL_CONFIG_FILE)
