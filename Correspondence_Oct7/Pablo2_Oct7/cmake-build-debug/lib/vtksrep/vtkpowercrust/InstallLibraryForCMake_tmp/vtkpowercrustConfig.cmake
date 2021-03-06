#-----------------------------------------------------------------------------
#
# vtkpowercrustConfig.cmake - CMake configuration file for external projects.
# This file was automatically generated by the cmake macro 
# CREA_INSTALL_LIBRARY_FOR_CMAKE of the package CREA
#
# This file is configured by cmake and used by the
# Usevtkpowercrust.cmake module to load the lib settings 
# for an external project.

# Build tree config ?
SET(CILC_BUILD_TREE_CONFIGURATION FALSE)


IF(UNIX)
SET(GOTO_INSTALL_PREFIX /../../..)
ENDIF(UNIX)


# The vtkpowercrust include file *RELATIVE* directories.
SET(CILC_RELATIVE_INCLUDE_PATHS "include/vtkpowercrust")
# Compute the prefix for include and library paths
IF(CILC_BUILD_TREE_CONFIGURATION)
  # In build tree 
  # the include paths are relative to the source tree *AND* the binary tree 
  # for generated files 
  SET(CILC_INCLUDE_PATH_PREFIX /playpen/workspace/newuoa/Correspondence_Oct7/Pablo2_Oct7)
  # Build the *ABSOLUTE* directories
  FOREACH(path ${CILC_RELATIVE_INCLUDE_PATHS})
    SET(vtkpowercrust_INCLUDE_DIRS
      ${vtkpowercrust_INCLUDE_DIRS}
      ${CILC_INCLUDE_PATH_PREFIX}/${path}
      )
  ENDFOREACH(path ${CILC_RELATIVE_INCLUDE_PATHS})
  SET(CILC_INCLUDE_PATH_PREFIX /playpen/workspace/newuoa/Correspondence_Oct7/Pablo2_Oct7/cmake-build-debug)
  # Build the *ABSOLUTE* directories
  FOREACH(path ${CILC_RELATIVE_INCLUDE_PATHS})
    SET(vtkpowercrust_INCLUDE_DIRS
      ${vtkpowercrust_INCLUDE_DIRS}
      ${CILC_INCLUDE_PATH_PREFIX}/${path}
      )
  ENDFOREACH(path ${CILC_RELATIVE_INCLUDE_PATHS})
ELSE(CILC_BUILD_TREE_CONFIGURATION)
  # In install tree 
  # the include paths are relative to install prefix 
  # On unix , GOTO_INSTALL_PREFIX allows to get back to the 
  # installation prefix from  vtkpowercrust_DIR
  SET(CILC_INCLUDE_PATH_PREFIX ${vtkpowercrust_DIR}${GOTO_INSTALL_PREFIX})
  # Build the *ABSOLUTE* directories
  FOREACH(path ${CILC_RELATIVE_INCLUDE_PATHS})
    SET(vtkpowercrust_INCLUDE_DIRS
      ${vtkpowercrust_INCLUDE_DIRS}
      ${CILC_INCLUDE_PATH_PREFIX}/${path}
      )
  ENDFOREACH(path ${CILC_RELATIVE_INCLUDE_PATHS})
ENDIF(CILC_BUILD_TREE_CONFIGURATION)



# Compute the prefix for library paths
IF(CILC_BUILD_TREE_CONFIGURATION)
  # In build tree
  # the library paths are relative to the binary tree 
  SET(CILC_LIBRARY_PATH_PREFIX /playpen/workspace/newuoa/Correspondence_Oct7/Pablo2_Oct7/cmake-build-debug)
ELSE(CILC_BUILD_TREE_CONFIGURATION)
  # In install tree 
  # the library paths are relative to install prefix
  SET(CILC_LIBRARY_PATH_PREFIX ${vtkpowercrust_DIR}${GOTO_INSTALL_PREFIX})
ENDIF(CILC_BUILD_TREE_CONFIGURATION)
# The vtkpowercrust library file *RELATIVE* directories.
SET(CILC_RELATIVE_LIBRARY_PATHS "lib64")
# Build the *ABSOLUTE* directories
FOREACH(path ${CILC_RELATIVE_LIBRARY_PATHS})
  SET(vtkpowercrust_LIBRARY_DIRS
    ${vtkpowercrust_LIBRARY_DIRS}
    ${CILC_LIBRARY_PATH_PREFIX}/${path}
    )
ENDFOREACH(path ${CILC_RELATIVE_LIBRARY_PATHS})

# Set the "prefix path"
SET(vtkpowercrust_INSTALL_PREFIX ${CILC_LIBRARY_PATH_PREFIX})

# The C and C++ flags added by vtkpowercrust to the cmake-configured flags.
SET(vtkpowercrust_REQUIRED_C_FLAGS "")
SET(vtkpowercrust_REQUIRED_CXX_FLAGS "")
SET(vtkpowercrust_REQUIRED_LINK_FLAGS "")

# The vtkpowercrust version 
SET(vtkpowercrust_MAJOR_VERSION )
SET(vtkpowercrust_MINOR_VERSION )
SET(vtkpowercrust_BUILD_VERSION )
SET(vtkpowercrust_VERSION ..)

# The location of the Usevtkpowercrust.cmake file.
SET(vtkpowercrust_USE_FILE "${vtkpowercrust_DIR}/Usevtkpowercrust.cmake")

# The build settings file.
SET(vtkpowercrust_BUILD_SETTINGS_FILE 
  "${vtkpowercrust_DIR}/vtkpowercrustBuildSettings.cmake")

# A list of all libraries for vtkpowercrust.  Those listed here should
# automatically pull in their dependencies.
SET(vtkpowercrust_LIBRARIES vtkpowercrust)

# Messages
IF(CREA_VERBOSE_CMAKE)
  MESSAGE(STATUS "=======================================")
  MESSAGE(STATUS "Looking for vtkpowercrust... found:")
  MESSAGE(STATUS "* vtkpowercrust_DIR          = ${vtkpowercrust_DIR}")
  MESSAGE(STATUS "* vtkpowercrust_VERSION      = ${vtkpowercrust_VERSION}")
  MESSAGE(STATUS "* vtkpowercrust_USE_FILE     = ${vtkpowercrust_USE_FILE}")

  MESSAGE(STATUS "* vtkpowercrust_INCLUDE_DIRS = ${vtkpowercrust_INCLUDE_DIRS}")
  MESSAGE(STATUS "* vtkpowercrust_LIBRARY_DIRS = ${vtkpowercrust_LIBRARY_DIRS}")
  MESSAGE(STATUS "* vtkpowercrust_LIBRARIES    = ${vtkpowercrust_LIBRARIES}")
ENDIF(CREA_VERBOSE_CMAKE)

# Does the library has an additional config file (user provided) ?
SET(vtkpowercrust_HAS_ADDITIONAL_CONFIG_FILE FALSE)

IF (vtkpowercrust_HAS_ADDITIONAL_CONFIG_FILE)
  IF(CREA_VERBOSE_CMAKE)
    MESSAGE(STATUS "Reading vtkpowercrust additional configuration file")
  ENDIF(CREA_VERBOSE_CMAKE)
  # Include it
  INCLUDE(${vtkpowercrust_DIR}/AdditionalvtkpowercrustConfig.cmake)
ENDIF (vtkpowercrust_HAS_ADDITIONAL_CONFIG_FILE)
