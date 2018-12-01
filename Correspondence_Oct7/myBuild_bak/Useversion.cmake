# This is an implementation detail for using version with the
# Findversion.cmake module.  Do not include directly by name.  
# This should be included only when Findversion.cmake sets 
# the version_USE_FILE variable to point here.

IF(CREA_VERBOSE_CMAKE)
  MESSAGE(STATUS "Using version")
ENDIF(CREA_VERBOSE_CMAKE)

# Load the compiler settings used for version.
IF(version_BUILD_SETTINGS_FILE)
  INCLUDE(CMakeImportBuildSettings)
  CMAKE_IMPORT_BUILD_SETTINGS(${version_BUILD_SETTINGS_FILE})
ENDIF(version_BUILD_SETTINGS_FILE)

# Add compiler flags needed to use version.
SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${version_REQUIRED_C_FLAGS}")
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${version_REQUIRED_CXX_FLAGS}")
SET(CMAKE_LINK_FLAGS "${CMAKE_LINK_FLAGS} ${version_REQUIRED_LINK_FLAGS}")

# Add include directories needed to use version.
INCLUDE_DIRECTORIES(${version_INCLUDE_DIRS})

# Add link directories needed to use version.
LINK_DIRECTORIES(${version_LIBRARY_DIRS})

# Set the version 
# Already done in bbtkConfigure.h
#ADD_DEFINITIONS( -Dversion_VERSION="\"${version_VERSION}\"" )

# Additional use file 
IF (version_HAS_ADDITIONAL_CONFIG_FILE)
  # Include it
  INCLUDE(${version_DIR}/AdditionalUseversion.cmake)
ENDIF (version_HAS_ADDITIONAL_CONFIG_FILE)
