# This is an implementation detail for using SRepIO with the
# FindSRepIO.cmake module.  Do not include directly by name.  
# This should be included only when FindSRepIO.cmake sets 
# the SRepIO_USE_FILE variable to point here.

IF(CREA_VERBOSE_CMAKE)
  MESSAGE(STATUS "Using SRepIO")
ENDIF(CREA_VERBOSE_CMAKE)

# Load the compiler settings used for SRepIO.
IF(SRepIO_BUILD_SETTINGS_FILE)
  INCLUDE(CMakeImportBuildSettings)
  CMAKE_IMPORT_BUILD_SETTINGS(${SRepIO_BUILD_SETTINGS_FILE})
ENDIF(SRepIO_BUILD_SETTINGS_FILE)

# Add compiler flags needed to use SRepIO.
SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${SRepIO_REQUIRED_C_FLAGS}")
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SRepIO_REQUIRED_CXX_FLAGS}")
SET(CMAKE_LINK_FLAGS "${CMAKE_LINK_FLAGS} ${SRepIO_REQUIRED_LINK_FLAGS}")

# Add include directories needed to use SRepIO.
INCLUDE_DIRECTORIES(${SRepIO_INCLUDE_DIRS})

# Add link directories needed to use SRepIO.
LINK_DIRECTORIES(${SRepIO_LIBRARY_DIRS})

# Set the version 
# Already done in bbtkConfigure.h
#ADD_DEFINITIONS( -DSRepIO_VERSION="\"${SRepIO_VERSION}\"" )

# Additional use file 
IF (SRepIO_HAS_ADDITIONAL_CONFIG_FILE)
  # Include it
  INCLUDE(${SRepIO_DIR}/AdditionalUseSRepIO.cmake)
ENDIF (SRepIO_HAS_ADDITIONAL_CONFIG_FILE)
