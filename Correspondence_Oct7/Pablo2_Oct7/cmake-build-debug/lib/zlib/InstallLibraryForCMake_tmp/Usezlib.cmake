# This is an implementation detail for using zlib with the
# Findzlib.cmake module.  Do not include directly by name.  
# This should be included only when Findzlib.cmake sets 
# the zlib_USE_FILE variable to point here.

IF(CREA_VERBOSE_CMAKE)
  MESSAGE(STATUS "Using zlib")
ENDIF(CREA_VERBOSE_CMAKE)

# Load the compiler settings used for zlib.
IF(zlib_BUILD_SETTINGS_FILE)
  INCLUDE(CMakeImportBuildSettings)
  CMAKE_IMPORT_BUILD_SETTINGS(${zlib_BUILD_SETTINGS_FILE})
ENDIF(zlib_BUILD_SETTINGS_FILE)

# Add compiler flags needed to use zlib.
SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${zlib_REQUIRED_C_FLAGS}")
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${zlib_REQUIRED_CXX_FLAGS}")
SET(CMAKE_LINK_FLAGS "${CMAKE_LINK_FLAGS} ${zlib_REQUIRED_LINK_FLAGS}")

# Add include directories needed to use zlib.
INCLUDE_DIRECTORIES(${zlib_INCLUDE_DIRS})

# Add link directories needed to use zlib.
LINK_DIRECTORIES(${zlib_LIBRARY_DIRS})

# Set the version 
# Already done in bbtkConfigure.h
#ADD_DEFINITIONS( -Dzlib_VERSION="\"${zlib_VERSION}\"" )

# Additional use file 
IF (zlib_HAS_ADDITIONAL_CONFIG_FILE)
  # Include it
  INCLUDE(${zlib_DIR}/AdditionalUsezlib.cmake)
ENDIF (zlib_HAS_ADDITIONAL_CONFIG_FILE)
