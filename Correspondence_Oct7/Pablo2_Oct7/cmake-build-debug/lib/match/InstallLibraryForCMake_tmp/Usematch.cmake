# This is an implementation detail for using match with the
# Findmatch.cmake module.  Do not include directly by name.  
# This should be included only when Findmatch.cmake sets 
# the match_USE_FILE variable to point here.

IF(CREA_VERBOSE_CMAKE)
  MESSAGE(STATUS "Using match")
ENDIF(CREA_VERBOSE_CMAKE)

# Load the compiler settings used for match.
IF(match_BUILD_SETTINGS_FILE)
  INCLUDE(CMakeImportBuildSettings)
  CMAKE_IMPORT_BUILD_SETTINGS(${match_BUILD_SETTINGS_FILE})
ENDIF(match_BUILD_SETTINGS_FILE)

# Add compiler flags needed to use match.
SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${match_REQUIRED_C_FLAGS}")
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${match_REQUIRED_CXX_FLAGS}")
SET(CMAKE_LINK_FLAGS "${CMAKE_LINK_FLAGS} ${match_REQUIRED_LINK_FLAGS}")

# Add include directories needed to use match.
INCLUDE_DIRECTORIES(${match_INCLUDE_DIRS})

# Add link directories needed to use match.
LINK_DIRECTORIES(${match_LIBRARY_DIRS})

# Set the version 
# Already done in bbtkConfigure.h
#ADD_DEFINITIONS( -Dmatch_VERSION="\"${match_VERSION}\"" )

# Additional use file 
IF (match_HAS_ADDITIONAL_CONFIG_FILE)
  # Include it
  INCLUDE(${match_DIR}/AdditionalUsematch.cmake)
ENDIF (match_HAS_ADDITIONAL_CONFIG_FILE)
