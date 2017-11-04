# This is an implementation detail for using paul_code with the
# Findpaul_code.cmake module.  Do not include directly by name.  
# This should be included only when Findpaul_code.cmake sets 
# the paul_code_USE_FILE variable to point here.

IF(CREA_VERBOSE_CMAKE)
  MESSAGE(STATUS "Using paul_code")
ENDIF(CREA_VERBOSE_CMAKE)

# Load the compiler settings used for paul_code.
IF(paul_code_BUILD_SETTINGS_FILE)
  INCLUDE(CMakeImportBuildSettings)
  CMAKE_IMPORT_BUILD_SETTINGS(${paul_code_BUILD_SETTINGS_FILE})
ENDIF(paul_code_BUILD_SETTINGS_FILE)

# Add compiler flags needed to use paul_code.
SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${paul_code_REQUIRED_C_FLAGS}")
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${paul_code_REQUIRED_CXX_FLAGS}")
SET(CMAKE_LINK_FLAGS "${CMAKE_LINK_FLAGS} ${paul_code_REQUIRED_LINK_FLAGS}")

# Add include directories needed to use paul_code.
INCLUDE_DIRECTORIES(${paul_code_INCLUDE_DIRS})

# Add link directories needed to use paul_code.
LINK_DIRECTORIES(${paul_code_LIBRARY_DIRS})

# Set the version 
# Already done in bbtkConfigure.h
#ADD_DEFINITIONS( -Dpaul_code_VERSION="\"${paul_code_VERSION}\"" )

# Additional use file 
IF (paul_code_HAS_ADDITIONAL_CONFIG_FILE)
  # Include it
  INCLUDE(${paul_code_DIR}/AdditionalUsepaul_code.cmake)
ENDIF (paul_code_HAS_ADDITIONAL_CONFIG_FILE)
