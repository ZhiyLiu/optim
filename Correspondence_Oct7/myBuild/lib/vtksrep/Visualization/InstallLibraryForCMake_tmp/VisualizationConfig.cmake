#-----------------------------------------------------------------------------
#
# VisualizationConfig.cmake - CMake configuration file for external projects.
# This file was automatically generated by the cmake macro 
# CREA_INSTALL_LIBRARY_FOR_CMAKE of the package CREA
#
# This file is configured by cmake and used by the
# UseVisualization.cmake module to load the lib settings 
# for an external project.

# Build tree config ?
SET(CILC_BUILD_TREE_CONFIGURATION FALSE)


IF(UNIX)
SET(GOTO_INSTALL_PREFIX /../../..)
ENDIF(UNIX)


# The Visualization include file *RELATIVE* directories.
SET(CILC_RELATIVE_INCLUDE_PATHS "include/Visualization")
# Compute the prefix for include and library paths
IF(CILC_BUILD_TREE_CONFIGURATION)
  # In build tree 
  # the include paths are relative to the source tree *AND* the binary tree 
  # for generated files 
  SET(CILC_INCLUDE_PATH_PREFIX /playpen/workspace/newuoa/Correspondence_Oct7/Pablo2_Oct7)
  # Build the *ABSOLUTE* directories
  FOREACH(path ${CILC_RELATIVE_INCLUDE_PATHS})
    SET(Visualization_INCLUDE_DIRS
      ${Visualization_INCLUDE_DIRS}
      ${CILC_INCLUDE_PATH_PREFIX}/${path}
      )
  ENDFOREACH(path ${CILC_RELATIVE_INCLUDE_PATHS})
  SET(CILC_INCLUDE_PATH_PREFIX /playpen/workspace/newuoa/Correspondence_Oct7/myBuild)
  # Build the *ABSOLUTE* directories
  FOREACH(path ${CILC_RELATIVE_INCLUDE_PATHS})
    SET(Visualization_INCLUDE_DIRS
      ${Visualization_INCLUDE_DIRS}
      ${CILC_INCLUDE_PATH_PREFIX}/${path}
      )
  ENDFOREACH(path ${CILC_RELATIVE_INCLUDE_PATHS})
ELSE(CILC_BUILD_TREE_CONFIGURATION)
  # In install tree 
  # the include paths are relative to install prefix 
  # On unix , GOTO_INSTALL_PREFIX allows to get back to the 
  # installation prefix from  Visualization_DIR
  SET(CILC_INCLUDE_PATH_PREFIX ${Visualization_DIR}${GOTO_INSTALL_PREFIX})
  # Build the *ABSOLUTE* directories
  FOREACH(path ${CILC_RELATIVE_INCLUDE_PATHS})
    SET(Visualization_INCLUDE_DIRS
      ${Visualization_INCLUDE_DIRS}
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
  SET(CILC_LIBRARY_PATH_PREFIX ${Visualization_DIR}${GOTO_INSTALL_PREFIX})
ENDIF(CILC_BUILD_TREE_CONFIGURATION)
# The Visualization library file *RELATIVE* directories.
SET(CILC_RELATIVE_LIBRARY_PATHS "lib64")
# Build the *ABSOLUTE* directories
FOREACH(path ${CILC_RELATIVE_LIBRARY_PATHS})
  SET(Visualization_LIBRARY_DIRS
    ${Visualization_LIBRARY_DIRS}
    ${CILC_LIBRARY_PATH_PREFIX}/${path}
    )
ENDFOREACH(path ${CILC_RELATIVE_LIBRARY_PATHS})

# Set the "prefix path"
SET(Visualization_INSTALL_PREFIX ${CILC_LIBRARY_PATH_PREFIX})

# The C and C++ flags added by Visualization to the cmake-configured flags.
SET(Visualization_REQUIRED_C_FLAGS "")
SET(Visualization_REQUIRED_CXX_FLAGS "")
SET(Visualization_REQUIRED_LINK_FLAGS "")

# The Visualization version 
SET(Visualization_MAJOR_VERSION )
SET(Visualization_MINOR_VERSION )
SET(Visualization_BUILD_VERSION )
SET(Visualization_VERSION ..)

# The location of the UseVisualization.cmake file.
SET(Visualization_USE_FILE "${Visualization_DIR}/UseVisualization.cmake")

# The build settings file.
SET(Visualization_BUILD_SETTINGS_FILE 
  "${Visualization_DIR}/VisualizationBuildSettings.cmake")

# A list of all libraries for Visualization.  Those listed here should
# automatically pull in their dependencies.
SET(Visualization_LIBRARIES Visualization)

# Messages
IF(CREA_VERBOSE_CMAKE)
  MESSAGE(STATUS "=======================================")
  MESSAGE(STATUS "Looking for Visualization... found:")
  MESSAGE(STATUS "* Visualization_DIR          = ${Visualization_DIR}")
  MESSAGE(STATUS "* Visualization_VERSION      = ${Visualization_VERSION}")
  MESSAGE(STATUS "* Visualization_USE_FILE     = ${Visualization_USE_FILE}")

  MESSAGE(STATUS "* Visualization_INCLUDE_DIRS = ${Visualization_INCLUDE_DIRS}")
  MESSAGE(STATUS "* Visualization_LIBRARY_DIRS = ${Visualization_LIBRARY_DIRS}")
  MESSAGE(STATUS "* Visualization_LIBRARIES    = ${Visualization_LIBRARIES}")
ENDIF(CREA_VERBOSE_CMAKE)

# Does the library has an additional config file (user provided) ?
SET(Visualization_HAS_ADDITIONAL_CONFIG_FILE FALSE)

IF (Visualization_HAS_ADDITIONAL_CONFIG_FILE)
  IF(CREA_VERBOSE_CMAKE)
    MESSAGE(STATUS "Reading Visualization additional configuration file")
  ENDIF(CREA_VERBOSE_CMAKE)
  # Include it
  INCLUDE(${Visualization_DIR}/AdditionalVisualizationConfig.cmake)
ENDIF (Visualization_HAS_ADDITIONAL_CONFIG_FILE)
