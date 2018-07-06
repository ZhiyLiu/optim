#-----------------------------------------------------------------------------
#
# registerConfig.cmake - CMake configuration file for external projects.
# This file was automatically generated by the cmake macro 
# CREA_INSTALL_LIBRARY_FOR_CMAKE of the package CREA
#
# This file is configured by cmake and used by the
# Useregister.cmake module to load the lib settings 
# for an external project.

# Build tree config ?
SET(CILC_BUILD_TREE_CONFIGURATION FALSE)


IF(UNIX)
SET(GOTO_INSTALL_PREFIX /../../..)
ENDIF(UNIX)


# The register include file *RELATIVE* directories.
SET(CILC_RELATIVE_INCLUDE_PATHS "include/register")
# Compute the prefix for include and library paths
IF(CILC_BUILD_TREE_CONFIGURATION)
  # In build tree 
  # the include paths are relative to the source tree *AND* the binary tree 
  # for generated files 
  SET(CILC_INCLUDE_PATH_PREFIX /playpen/workspace/newuoa/Correspondence_Oct7/Pablo2_Oct7)
  # Build the *ABSOLUTE* directories
  FOREACH(path ${CILC_RELATIVE_INCLUDE_PATHS})
    SET(register_INCLUDE_DIRS
      ${register_INCLUDE_DIRS}
      ${CILC_INCLUDE_PATH_PREFIX}/${path}
      )
  ENDFOREACH(path ${CILC_RELATIVE_INCLUDE_PATHS})
  SET(CILC_INCLUDE_PATH_PREFIX /playpen/workspace/newuoa/Correspondence_Oct7/myBuild)
  # Build the *ABSOLUTE* directories
  FOREACH(path ${CILC_RELATIVE_INCLUDE_PATHS})
    SET(register_INCLUDE_DIRS
      ${register_INCLUDE_DIRS}
      ${CILC_INCLUDE_PATH_PREFIX}/${path}
      )
  ENDFOREACH(path ${CILC_RELATIVE_INCLUDE_PATHS})
ELSE(CILC_BUILD_TREE_CONFIGURATION)
  # In install tree 
  # the include paths are relative to install prefix 
  # On unix , GOTO_INSTALL_PREFIX allows to get back to the 
  # installation prefix from  register_DIR
  SET(CILC_INCLUDE_PATH_PREFIX ${register_DIR}${GOTO_INSTALL_PREFIX})
  # Build the *ABSOLUTE* directories
  FOREACH(path ${CILC_RELATIVE_INCLUDE_PATHS})
    SET(register_INCLUDE_DIRS
      ${register_INCLUDE_DIRS}
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
  SET(CILC_LIBRARY_PATH_PREFIX ${register_DIR}${GOTO_INSTALL_PREFIX})
ENDIF(CILC_BUILD_TREE_CONFIGURATION)
# The register library file *RELATIVE* directories.
SET(CILC_RELATIVE_LIBRARY_PATHS "lib64")
# Build the *ABSOLUTE* directories
FOREACH(path ${CILC_RELATIVE_LIBRARY_PATHS})
  SET(register_LIBRARY_DIRS
    ${register_LIBRARY_DIRS}
    ${CILC_LIBRARY_PATH_PREFIX}/${path}
    )
ENDFOREACH(path ${CILC_RELATIVE_LIBRARY_PATHS})

# Set the "prefix path"
SET(register_INSTALL_PREFIX ${CILC_LIBRARY_PATH_PREFIX})

# The C and C++ flags added by register to the cmake-configured flags.
SET(register_REQUIRED_C_FLAGS "")
SET(register_REQUIRED_CXX_FLAGS "")
SET(register_REQUIRED_LINK_FLAGS "")

# The register version 
SET(register_MAJOR_VERSION )
SET(register_MINOR_VERSION )
SET(register_BUILD_VERSION )
SET(register_VERSION ..)

# The location of the Useregister.cmake file.
SET(register_USE_FILE "${register_DIR}/Useregister.cmake")

# The build settings file.
SET(register_BUILD_SETTINGS_FILE 
  "${register_DIR}/registerBuildSettings.cmake")

# A list of all libraries for register.  Those listed here should
# automatically pull in their dependencies.
SET(register_LIBRARIES register)

# Messages
IF(CREA_VERBOSE_CMAKE)
  MESSAGE(STATUS "=======================================")
  MESSAGE(STATUS "Looking for register... found:")
  MESSAGE(STATUS "* register_DIR          = ${register_DIR}")
  MESSAGE(STATUS "* register_VERSION      = ${register_VERSION}")
  MESSAGE(STATUS "* register_USE_FILE     = ${register_USE_FILE}")

  MESSAGE(STATUS "* register_INCLUDE_DIRS = ${register_INCLUDE_DIRS}")
  MESSAGE(STATUS "* register_LIBRARY_DIRS = ${register_LIBRARY_DIRS}")
  MESSAGE(STATUS "* register_LIBRARIES    = ${register_LIBRARIES}")
ENDIF(CREA_VERBOSE_CMAKE)

# Does the library has an additional config file (user provided) ?
SET(register_HAS_ADDITIONAL_CONFIG_FILE FALSE)

IF (register_HAS_ADDITIONAL_CONFIG_FILE)
  IF(CREA_VERBOSE_CMAKE)
    MESSAGE(STATUS "Reading register additional configuration file")
  ENDIF(CREA_VERBOSE_CMAKE)
  # Include it
  INCLUDE(${register_DIR}/AdditionalregisterConfig.cmake)
ENDIF (register_HAS_ADDITIONAL_CONFIG_FILE)
