# This is an implementation detail for using Correspondence with the
# FindCorrespondence.cmake module.  Do not include directly by name.  
# This should be included only when FindCorrespondence.cmake sets 
# the Correspondence_USE_FILE variable to point here.

IF(CREA_VERBOSE_CMAKE)
  MESSAGE(STATUS "Using Correspondence")
ENDIF(CREA_VERBOSE_CMAKE)

# Load the compiler settings used for Correspondence.
IF(Correspondence_BUILD_SETTINGS_FILE)
  INCLUDE(CMakeImportBuildSettings)
  CMAKE_IMPORT_BUILD_SETTINGS(${Correspondence_BUILD_SETTINGS_FILE})
ENDIF(Correspondence_BUILD_SETTINGS_FILE)

# Add compiler flags needed to use Correspondence.
SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${Correspondence_REQUIRED_C_FLAGS}")
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${Correspondence_REQUIRED_CXX_FLAGS}")
SET(CMAKE_LINK_FLAGS "${CMAKE_LINK_FLAGS} ${Correspondence_REQUIRED_LINK_FLAGS}")

# Add include directories needed to use Correspondence.
INCLUDE_DIRECTORIES(${Correspondence_INCLUDE_DIRS})

# Add link directories needed to use Correspondence.
LINK_DIRECTORIES(${Correspondence_LIBRARY_DIRS})

# Set the version 
# Already done in bbtkConfigure.h
#ADD_DEFINITIONS( -DCorrespondence_VERSION="\"${Correspondence_VERSION}\"" )

# Additional use file 
IF (Correspondence_HAS_ADDITIONAL_CONFIG_FILE)
  # Include it
  INCLUDE(${Correspondence_DIR}/AdditionalUseCorrespondence.cmake)
ENDIF (Correspondence_HAS_ADDITIONAL_CONFIG_FILE)
