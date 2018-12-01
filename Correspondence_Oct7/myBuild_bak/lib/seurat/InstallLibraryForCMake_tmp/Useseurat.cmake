# This is an implementation detail for using seurat with the
# Findseurat.cmake module.  Do not include directly by name.  
# This should be included only when Findseurat.cmake sets 
# the seurat_USE_FILE variable to point here.

IF(CREA_VERBOSE_CMAKE)
  MESSAGE(STATUS "Using seurat")
ENDIF(CREA_VERBOSE_CMAKE)

# Load the compiler settings used for seurat.
IF(seurat_BUILD_SETTINGS_FILE)
  INCLUDE(CMakeImportBuildSettings)
  CMAKE_IMPORT_BUILD_SETTINGS(${seurat_BUILD_SETTINGS_FILE})
ENDIF(seurat_BUILD_SETTINGS_FILE)

# Add compiler flags needed to use seurat.
SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${seurat_REQUIRED_C_FLAGS}")
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${seurat_REQUIRED_CXX_FLAGS}")
SET(CMAKE_LINK_FLAGS "${CMAKE_LINK_FLAGS} ${seurat_REQUIRED_LINK_FLAGS}")

# Add include directories needed to use seurat.
INCLUDE_DIRECTORIES(${seurat_INCLUDE_DIRS})

# Add link directories needed to use seurat.
LINK_DIRECTORIES(${seurat_LIBRARY_DIRS})

# Set the version 
# Already done in bbtkConfigure.h
#ADD_DEFINITIONS( -Dseurat_VERSION="\"${seurat_VERSION}\"" )

# Additional use file 
IF (seurat_HAS_ADDITIONAL_CONFIG_FILE)
  # Include it
  INCLUDE(${seurat_DIR}/AdditionalUseseurat.cmake)
ENDIF (seurat_HAS_ADDITIONAL_CONFIG_FILE)
