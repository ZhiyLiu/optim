#-----------------------------------------------------------------------------
#
# seuratConfig.cmake - CMake configuration file for external projects.
# This file was automatically generated by the cmake macro 
# CREA_INSTALL_LIBRARY_FOR_CMAKE of the package CREA
#
# This file is configured by cmake and used by the
# Useseurat.cmake module to load the lib settings 
# for an external project.

# Build tree config ?
SET(CILC_BUILD_TREE_CONFIGURATION TRUE)


IF(UNIX)
SET(GOTO_INSTALL_PREFIX /../../..)
ENDIF(UNIX)


# The seurat include file *RELATIVE* directories.
SET(CILC_RELATIVE_INCLUDE_PATHS "lib/seurat")
# Compute the prefix for include and library paths
IF(CILC_BUILD_TREE_CONFIGURATION)
  # In build tree 
  # the include paths are relative to the source tree *AND* the binary tree 
  # for generated files 
  SET(CILC_INCLUDE_PATH_PREFIX /playpen/workspace/newuoa/Correspondence_Oct7/Pablo2_Oct7)
  # Build the *ABSOLUTE* directories
  FOREACH(path ${CILC_RELATIVE_INCLUDE_PATHS})
    SET(seurat_INCLUDE_DIRS
      ${seurat_INCLUDE_DIRS}
      ${CILC_INCLUDE_PATH_PREFIX}/${path}
      )
  ENDFOREACH(path ${CILC_RELATIVE_INCLUDE_PATHS})
  SET(CILC_INCLUDE_PATH_PREFIX /playpen/workspace/newuoa/Correspondence_Oct7/myBuild)
  # Build the *ABSOLUTE* directories
  FOREACH(path ${CILC_RELATIVE_INCLUDE_PATHS})
    SET(seurat_INCLUDE_DIRS
      ${seurat_INCLUDE_DIRS}
      ${CILC_INCLUDE_PATH_PREFIX}/${path}
      )
  ENDFOREACH(path ${CILC_RELATIVE_INCLUDE_PATHS})
ELSE(CILC_BUILD_TREE_CONFIGURATION)
  # In install tree 
  # the include paths are relative to install prefix 
  # On unix , GOTO_INSTALL_PREFIX allows to get back to the 
  # installation prefix from  seurat_DIR
  SET(CILC_INCLUDE_PATH_PREFIX ${seurat_DIR}${GOTO_INSTALL_PREFIX})
  # Build the *ABSOLUTE* directories
  FOREACH(path ${CILC_RELATIVE_INCLUDE_PATHS})
    SET(seurat_INCLUDE_DIRS
      ${seurat_INCLUDE_DIRS}
      ${CILC_INCLUDE_PATH_PREFIX}/${path}
      )
  ENDFOREACH(path ${CILC_RELATIVE_INCLUDE_PATHS})
ENDIF(CILC_BUILD_TREE_CONFIGURATION)



# Compute the prefix for library paths
IF(CILC_BUILD_TREE_CONFIGURATION)
  # In build tree
  # the library paths are relative to the binary tree 
  SET(CILC_LIBRARY_PATH_PREFIX /playpen/workspace/newuoa/Correspondence_Oct7/myBuild)
ELSE(CILC_BUILD_TREE_CONFIGURATION)
  # In install tree 
  # the library paths are relative to install prefix
  SET(CILC_LIBRARY_PATH_PREFIX ${seurat_DIR}${GOTO_INSTALL_PREFIX})
ENDIF(CILC_BUILD_TREE_CONFIGURATION)
# The seurat library file *RELATIVE* directories.
SET(CILC_RELATIVE_LIBRARY_PATHS ".")
# Build the *ABSOLUTE* directories
FOREACH(path ${CILC_RELATIVE_LIBRARY_PATHS})
  SET(seurat_LIBRARY_DIRS
    ${seurat_LIBRARY_DIRS}
    ${CILC_LIBRARY_PATH_PREFIX}/${path}
    )
ENDFOREACH(path ${CILC_RELATIVE_LIBRARY_PATHS})

# Set the "prefix path"
SET(seurat_INSTALL_PREFIX ${CILC_LIBRARY_PATH_PREFIX})

# The C and C++ flags added by seurat to the cmake-configured flags.
SET(seurat_REQUIRED_C_FLAGS "")
SET(seurat_REQUIRED_CXX_FLAGS "")
SET(seurat_REQUIRED_LINK_FLAGS "")

# The seurat version 
SET(seurat_MAJOR_VERSION )
SET(seurat_MINOR_VERSION )
SET(seurat_BUILD_VERSION )
SET(seurat_VERSION ..)

# The location of the Useseurat.cmake file.
SET(seurat_USE_FILE "${seurat_DIR}/Useseurat.cmake")

# The build settings file.
SET(seurat_BUILD_SETTINGS_FILE 
  "${seurat_DIR}/seuratBuildSettings.cmake")

# A list of all libraries for seurat.  Those listed here should
# automatically pull in their dependencies.
SET(seurat_LIBRARIES seurat)

# Messages
IF(CREA_VERBOSE_CMAKE)
  MESSAGE(STATUS "=======================================")
  MESSAGE(STATUS "Looking for seurat... found:")
  MESSAGE(STATUS "* seurat_DIR          = ${seurat_DIR}")
  MESSAGE(STATUS "* seurat_VERSION      = ${seurat_VERSION}")
  MESSAGE(STATUS "* seurat_USE_FILE     = ${seurat_USE_FILE}")

  MESSAGE(STATUS "* seurat_INCLUDE_DIRS = ${seurat_INCLUDE_DIRS}")
  MESSAGE(STATUS "* seurat_LIBRARY_DIRS = ${seurat_LIBRARY_DIRS}")
  MESSAGE(STATUS "* seurat_LIBRARIES    = ${seurat_LIBRARIES}")
ENDIF(CREA_VERBOSE_CMAKE)

# Does the library has an additional config file (user provided) ?
SET(seurat_HAS_ADDITIONAL_CONFIG_FILE FALSE)

IF (seurat_HAS_ADDITIONAL_CONFIG_FILE)
  IF(CREA_VERBOSE_CMAKE)
    MESSAGE(STATUS "Reading seurat additional configuration file")
  ENDIF(CREA_VERBOSE_CMAKE)
  # Include it
  INCLUDE(${seurat_DIR}/AdditionalseuratConfig.cmake)
ENDIF (seurat_HAS_ADDITIONAL_CONFIG_FILE)
