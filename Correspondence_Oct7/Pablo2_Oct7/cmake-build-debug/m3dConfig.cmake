#-----------------------------------------------------------------------------
#
# m3dConfig.cmake - CMake configuration file for external projects.
# This file was automatically generated by the cmake macro 
# CREA_INSTALL_LIBRARY_FOR_CMAKE of the package CREA
#
# This file is configured by cmake and used by the
# Usem3d.cmake module to load the lib settings 
# for an external project.

# Build tree config ?
SET(CILC_BUILD_TREE_CONFIGURATION TRUE)


IF(UNIX)
SET(GOTO_INSTALL_PREFIX /../../..)
ENDIF(UNIX)


# The m3d include file *RELATIVE* directories.
SET(CILC_RELATIVE_INCLUDE_PATHS "lib/m3d")
# Compute the prefix for include and library paths
IF(CILC_BUILD_TREE_CONFIGURATION)
  # In build tree 
  # the include paths are relative to the source tree *AND* the binary tree 
  # for generated files 
  SET(CILC_INCLUDE_PATH_PREFIX /playpen/workspace/newuoa/Correspondence_Oct7/Pablo2_Oct7)
  # Build the *ABSOLUTE* directories
  FOREACH(path ${CILC_RELATIVE_INCLUDE_PATHS})
    SET(m3d_INCLUDE_DIRS
      ${m3d_INCLUDE_DIRS}
      ${CILC_INCLUDE_PATH_PREFIX}/${path}
      )
  ENDFOREACH(path ${CILC_RELATIVE_INCLUDE_PATHS})
  SET(CILC_INCLUDE_PATH_PREFIX /playpen/workspace/newuoa/Correspondence_Oct7/Pablo2_Oct7/cmake-build-debug)
  # Build the *ABSOLUTE* directories
  FOREACH(path ${CILC_RELATIVE_INCLUDE_PATHS})
    SET(m3d_INCLUDE_DIRS
      ${m3d_INCLUDE_DIRS}
      ${CILC_INCLUDE_PATH_PREFIX}/${path}
      )
  ENDFOREACH(path ${CILC_RELATIVE_INCLUDE_PATHS})
ELSE(CILC_BUILD_TREE_CONFIGURATION)
  # In install tree 
  # the include paths are relative to install prefix 
  # On unix , GOTO_INSTALL_PREFIX allows to get back to the 
  # installation prefix from  m3d_DIR
  SET(CILC_INCLUDE_PATH_PREFIX ${m3d_DIR}${GOTO_INSTALL_PREFIX})
  # Build the *ABSOLUTE* directories
  FOREACH(path ${CILC_RELATIVE_INCLUDE_PATHS})
    SET(m3d_INCLUDE_DIRS
      ${m3d_INCLUDE_DIRS}
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
  SET(CILC_LIBRARY_PATH_PREFIX ${m3d_DIR}${GOTO_INSTALL_PREFIX})
ENDIF(CILC_BUILD_TREE_CONFIGURATION)
# The m3d library file *RELATIVE* directories.
SET(CILC_RELATIVE_LIBRARY_PATHS ".")
# Build the *ABSOLUTE* directories
FOREACH(path ${CILC_RELATIVE_LIBRARY_PATHS})
  SET(m3d_LIBRARY_DIRS
    ${m3d_LIBRARY_DIRS}
    ${CILC_LIBRARY_PATH_PREFIX}/${path}
    )
ENDFOREACH(path ${CILC_RELATIVE_LIBRARY_PATHS})

# Set the "prefix path"
SET(m3d_INSTALL_PREFIX ${CILC_LIBRARY_PATH_PREFIX})

# The C and C++ flags added by m3d to the cmake-configured flags.
SET(m3d_REQUIRED_C_FLAGS "")
SET(m3d_REQUIRED_CXX_FLAGS "")
SET(m3d_REQUIRED_LINK_FLAGS "")

# The m3d version 
SET(m3d_MAJOR_VERSION )
SET(m3d_MINOR_VERSION )
SET(m3d_BUILD_VERSION )
SET(m3d_VERSION ..)

# The location of the Usem3d.cmake file.
SET(m3d_USE_FILE "${m3d_DIR}/Usem3d.cmake")

# The build settings file.
SET(m3d_BUILD_SETTINGS_FILE 
  "${m3d_DIR}/m3dBuildSettings.cmake")

# A list of all libraries for m3d.  Those listed here should
# automatically pull in their dependencies.
SET(m3d_LIBRARIES m3d)

# Messages
IF(CREA_VERBOSE_CMAKE)
  MESSAGE(STATUS "=======================================")
  MESSAGE(STATUS "Looking for m3d... found:")
  MESSAGE(STATUS "* m3d_DIR          = ${m3d_DIR}")
  MESSAGE(STATUS "* m3d_VERSION      = ${m3d_VERSION}")
  MESSAGE(STATUS "* m3d_USE_FILE     = ${m3d_USE_FILE}")

  MESSAGE(STATUS "* m3d_INCLUDE_DIRS = ${m3d_INCLUDE_DIRS}")
  MESSAGE(STATUS "* m3d_LIBRARY_DIRS = ${m3d_LIBRARY_DIRS}")
  MESSAGE(STATUS "* m3d_LIBRARIES    = ${m3d_LIBRARIES}")
ENDIF(CREA_VERBOSE_CMAKE)

# Does the library has an additional config file (user provided) ?
SET(m3d_HAS_ADDITIONAL_CONFIG_FILE FALSE)

IF (m3d_HAS_ADDITIONAL_CONFIG_FILE)
  IF(CREA_VERBOSE_CMAKE)
    MESSAGE(STATUS "Reading m3d additional configuration file")
  ENDIF(CREA_VERBOSE_CMAKE)
  # Include it
  INCLUDE(${m3d_DIR}/Additionalm3dConfig.cmake)
ENDIF (m3d_HAS_ADDITIONAL_CONFIG_FILE)
