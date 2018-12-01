# This is an implementation detail for using ImageIO with the
# FindImageIO.cmake module.  Do not include directly by name.  
# This should be included only when FindImageIO.cmake sets 
# the ImageIO_USE_FILE variable to point here.

IF(CREA_VERBOSE_CMAKE)
  MESSAGE(STATUS "Using ImageIO")
ENDIF(CREA_VERBOSE_CMAKE)

# Load the compiler settings used for ImageIO.
IF(ImageIO_BUILD_SETTINGS_FILE)
  INCLUDE(CMakeImportBuildSettings)
  CMAKE_IMPORT_BUILD_SETTINGS(${ImageIO_BUILD_SETTINGS_FILE})
ENDIF(ImageIO_BUILD_SETTINGS_FILE)

# Add compiler flags needed to use ImageIO.
SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${ImageIO_REQUIRED_C_FLAGS}")
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${ImageIO_REQUIRED_CXX_FLAGS}")
SET(CMAKE_LINK_FLAGS "${CMAKE_LINK_FLAGS} ${ImageIO_REQUIRED_LINK_FLAGS}")

# Add include directories needed to use ImageIO.
INCLUDE_DIRECTORIES(${ImageIO_INCLUDE_DIRS})

# Add link directories needed to use ImageIO.
LINK_DIRECTORIES(${ImageIO_LIBRARY_DIRS})

# Set the version 
# Already done in bbtkConfigure.h
#ADD_DEFINITIONS( -DImageIO_VERSION="\"${ImageIO_VERSION}\"" )

# Additional use file 
IF (ImageIO_HAS_ADDITIONAL_CONFIG_FILE)
  # Include it
  INCLUDE(${ImageIO_DIR}/AdditionalUseImageIO.cmake)
ENDIF (ImageIO_HAS_ADDITIONAL_CONFIG_FILE)
