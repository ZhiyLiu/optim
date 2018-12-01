# This is an implementation detail for using register with the
# Findregister.cmake module.  Do not include directly by name.  
# This should be included only when Findregister.cmake sets 
# the register_USE_FILE variable to point here.

IF(CREA_VERBOSE_CMAKE)
  MESSAGE(STATUS "Using register")
ENDIF(CREA_VERBOSE_CMAKE)

# Load the compiler settings used for register.
IF(register_BUILD_SETTINGS_FILE)
  INCLUDE(CMakeImportBuildSettings)
  CMAKE_IMPORT_BUILD_SETTINGS(${register_BUILD_SETTINGS_FILE})
ENDIF(register_BUILD_SETTINGS_FILE)

# Add compiler flags needed to use register.
SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${register_REQUIRED_C_FLAGS}")
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${register_REQUIRED_CXX_FLAGS}")
SET(CMAKE_LINK_FLAGS "${CMAKE_LINK_FLAGS} ${register_REQUIRED_LINK_FLAGS}")

# Add include directories needed to use register.
INCLUDE_DIRECTORIES(${register_INCLUDE_DIRS})

# Add link directories needed to use register.
LINK_DIRECTORIES(${register_LIBRARY_DIRS})

# Set the version 
# Already done in bbtkConfigure.h
#ADD_DEFINITIONS( -Dregister_VERSION="\"${register_VERSION}\"" )

# Additional use file 
IF (register_HAS_ADDITIONAL_CONFIG_FILE)
  # Include it
  INCLUDE(${register_DIR}/AdditionalUseregister.cmake)
ENDIF (register_HAS_ADDITIONAL_CONFIG_FILE)
