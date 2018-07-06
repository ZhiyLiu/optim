# This is an implementation detail for using calquantum with the
# Findcalquantum.cmake module.  Do not include directly by name.  
# This should be included only when Findcalquantum.cmake sets 
# the calquantum_USE_FILE variable to point here.

IF(CREA_VERBOSE_CMAKE)
  MESSAGE(STATUS "Using calquantum")
ENDIF(CREA_VERBOSE_CMAKE)

# Load the compiler settings used for calquantum.
IF(calquantum_BUILD_SETTINGS_FILE)
  INCLUDE(CMakeImportBuildSettings)
  CMAKE_IMPORT_BUILD_SETTINGS(${calquantum_BUILD_SETTINGS_FILE})
ENDIF(calquantum_BUILD_SETTINGS_FILE)

# Add compiler flags needed to use calquantum.
SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${calquantum_REQUIRED_C_FLAGS}")
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${calquantum_REQUIRED_CXX_FLAGS}")
SET(CMAKE_LINK_FLAGS "${CMAKE_LINK_FLAGS} ${calquantum_REQUIRED_LINK_FLAGS}")

# Add include directories needed to use calquantum.
INCLUDE_DIRECTORIES(${calquantum_INCLUDE_DIRS})

# Add link directories needed to use calquantum.
LINK_DIRECTORIES(${calquantum_LIBRARY_DIRS})

# Set the version 
# Already done in bbtkConfigure.h
#ADD_DEFINITIONS( -Dcalquantum_VERSION="\"${calquantum_VERSION}\"" )

# Additional use file 
IF (calquantum_HAS_ADDITIONAL_CONFIG_FILE)
  # Include it
  INCLUDE(${calquantum_DIR}/AdditionalUsecalquantum.cmake)
ENDIF (calquantum_HAS_ADDITIONAL_CONFIG_FILE)
