# This is an implementation detail for using vtkpowercrust with the
# Findvtkpowercrust.cmake module.  Do not include directly by name.  
# This should be included only when Findvtkpowercrust.cmake sets 
# the vtkpowercrust_USE_FILE variable to point here.

IF(CREA_VERBOSE_CMAKE)
  MESSAGE(STATUS "Using vtkpowercrust")
ENDIF(CREA_VERBOSE_CMAKE)

# Load the compiler settings used for vtkpowercrust.
IF(vtkpowercrust_BUILD_SETTINGS_FILE)
  INCLUDE(CMakeImportBuildSettings)
  CMAKE_IMPORT_BUILD_SETTINGS(${vtkpowercrust_BUILD_SETTINGS_FILE})
ENDIF(vtkpowercrust_BUILD_SETTINGS_FILE)

# Add compiler flags needed to use vtkpowercrust.
SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${vtkpowercrust_REQUIRED_C_FLAGS}")
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${vtkpowercrust_REQUIRED_CXX_FLAGS}")
SET(CMAKE_LINK_FLAGS "${CMAKE_LINK_FLAGS} ${vtkpowercrust_REQUIRED_LINK_FLAGS}")

# Add include directories needed to use vtkpowercrust.
INCLUDE_DIRECTORIES(${vtkpowercrust_INCLUDE_DIRS})

# Add link directories needed to use vtkpowercrust.
LINK_DIRECTORIES(${vtkpowercrust_LIBRARY_DIRS})

# Set the version 
# Already done in bbtkConfigure.h
#ADD_DEFINITIONS( -Dvtkpowercrust_VERSION="\"${vtkpowercrust_VERSION}\"" )

# Additional use file 
IF (vtkpowercrust_HAS_ADDITIONAL_CONFIG_FILE)
  # Include it
  INCLUDE(${vtkpowercrust_DIR}/AdditionalUsevtkpowercrust.cmake)
ENDIF (vtkpowercrust_HAS_ADDITIONAL_CONFIG_FILE)
