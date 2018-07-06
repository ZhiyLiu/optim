# This is an implementation detail for using SRep with the
# FindSRep.cmake module.  Do not include directly by name.  
# This should be included only when FindSRep.cmake sets 
# the SRep_USE_FILE variable to point here.

IF(CREA_VERBOSE_CMAKE)
  MESSAGE(STATUS "Using SRep")
ENDIF(CREA_VERBOSE_CMAKE)

# Load the compiler settings used for SRep.
IF(SRep_BUILD_SETTINGS_FILE)
  INCLUDE(CMakeImportBuildSettings)
  CMAKE_IMPORT_BUILD_SETTINGS(${SRep_BUILD_SETTINGS_FILE})
ENDIF(SRep_BUILD_SETTINGS_FILE)

# Add compiler flags needed to use SRep.
SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${SRep_REQUIRED_C_FLAGS}")
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SRep_REQUIRED_CXX_FLAGS}")
SET(CMAKE_LINK_FLAGS "${CMAKE_LINK_FLAGS} ${SRep_REQUIRED_LINK_FLAGS}")

# Add include directories needed to use SRep.
INCLUDE_DIRECTORIES(${SRep_INCLUDE_DIRS})

# Add link directories needed to use SRep.
LINK_DIRECTORIES(${SRep_LIBRARY_DIRS})

# Set the version 
# Already done in bbtkConfigure.h
#ADD_DEFINITIONS( -DSRep_VERSION="\"${SRep_VERSION}\"" )

# Additional use file 
IF (SRep_HAS_ADDITIONAL_CONFIG_FILE)
  # Include it
  INCLUDE(${SRep_DIR}/AdditionalUseSRep.cmake)
ENDIF (SRep_HAS_ADDITIONAL_CONFIG_FILE)
