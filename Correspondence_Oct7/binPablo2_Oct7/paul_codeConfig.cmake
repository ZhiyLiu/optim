#-----------------------------------------------------------------------------
#
# paul_codeConfig.cmake - CMake configuration file for external projects.
# This file was automatically generated by the cmake macro 
# CREA_INSTALL_LIBRARY_FOR_CMAKE of the package CREA
#
# This file is configured by cmake and used by the
# Usepaul_code.cmake module to load the lib settings 
# for an external project.

# Build tree config ?
SET(CILC_BUILD_TREE_CONFIGURATION TRUE)


IF(UNIX)
SET(GOTO_INSTALL_PREFIX /../../..)
ENDIF(UNIX)


# The paul_code include file *RELATIVE* directories.
SET(CILC_RELATIVE_INCLUDE_PATHS "lib/paul_code")
# Compute the prefix for include and library paths
IF(CILC_BUILD_TREE_CONFIGURATION)
  # In build tree 
  # the include paths are relative to the source tree *AND* the binary tree 
  # for generated files 
  SET(CILC_INCLUDE_PATH_PREFIX /work/ltu/Correspondence_Oct7/Pablo2_Oct7)
  # Build the *ABSOLUTE* directories
  FOREACH(path ${CILC_RELATIVE_INCLUDE_PATHS})
    SET(paul_code_INCLUDE_DIRS
      ${paul_code_INCLUDE_DIRS}
      ${CILC_INCLUDE_PATH_PREFIX}/${path}
      )
  ENDFOREACH(path ${CILC_RELATIVE_INCLUDE_PATHS})
  SET(CILC_INCLUDE_PATH_PREFIX /work/ltu/Correspondence_Oct7/binPablo2_Oct7)
  # Build the *ABSOLUTE* directories
  FOREACH(path ${CILC_RELATIVE_INCLUDE_PATHS})
    SET(paul_code_INCLUDE_DIRS
      ${paul_code_INCLUDE_DIRS}
      ${CILC_INCLUDE_PATH_PREFIX}/${path}
      )
  ENDFOREACH(path ${CILC_RELATIVE_INCLUDE_PATHS})
ELSE(CILC_BUILD_TREE_CONFIGURATION)
  # In install tree 
  # the include paths are relative to install prefix 
  # On unix , GOTO_INSTALL_PREFIX allows to get back to the 
  # installation prefix from  paul_code_DIR
  SET(CILC_INCLUDE_PATH_PREFIX ${paul_code_DIR}${GOTO_INSTALL_PREFIX})
  # Build the *ABSOLUTE* directories
  FOREACH(path ${CILC_RELATIVE_INCLUDE_PATHS})
    SET(paul_code_INCLUDE_DIRS
      ${paul_code_INCLUDE_DIRS}
      ${CILC_INCLUDE_PATH_PREFIX}/${path}
      )
  ENDFOREACH(path ${CILC_RELATIVE_INCLUDE_PATHS})
ENDIF(CILC_BUILD_TREE_CONFIGURATION)



# Compute the prefix for library paths
IF(CILC_BUILD_TREE_CONFIGURATION)
  # In build tree
  # the library paths are relative to the binary tree 
  SET(CILC_LIBRARY_PATH_PREFIX /work/ltu/Correspondence_Oct7/binPablo2_Oct7)
ELSE(CILC_BUILD_TREE_CONFIGURATION)
  # In install tree 
  # the library paths are relative to install prefix
  SET(CILC_LIBRARY_PATH_PREFIX ${paul_code_DIR}${GOTO_INSTALL_PREFIX})
ENDIF(CILC_BUILD_TREE_CONFIGURATION)
# The paul_code library file *RELATIVE* directories.
SET(CILC_RELATIVE_LIBRARY_PATHS ".")
# Build the *ABSOLUTE* directories
FOREACH(path ${CILC_RELATIVE_LIBRARY_PATHS})
  SET(paul_code_LIBRARY_DIRS
    ${paul_code_LIBRARY_DIRS}
    ${CILC_LIBRARY_PATH_PREFIX}/${path}
    )
ENDFOREACH(path ${CILC_RELATIVE_LIBRARY_PATHS})

# Set the "prefix path"
SET(paul_code_INSTALL_PREFIX ${CILC_LIBRARY_PATH_PREFIX})

# The C and C++ flags added by paul_code to the cmake-configured flags.
SET(paul_code_REQUIRED_C_FLAGS "")
SET(paul_code_REQUIRED_CXX_FLAGS "")
SET(paul_code_REQUIRED_LINK_FLAGS "")

# The paul_code version 
SET(paul_code_MAJOR_VERSION )
SET(paul_code_MINOR_VERSION )
SET(paul_code_BUILD_VERSION )
SET(paul_code_VERSION ..)

# The location of the Usepaul_code.cmake file.
SET(paul_code_USE_FILE "${paul_code_DIR}/Usepaul_code.cmake")

# The build settings file.
SET(paul_code_BUILD_SETTINGS_FILE 
  "${paul_code_DIR}/paul_codeBuildSettings.cmake")

# A list of all libraries for paul_code.  Those listed here should
# automatically pull in their dependencies.
SET(paul_code_LIBRARIES paul_code)

# Messages
IF(CREA_VERBOSE_CMAKE)
  MESSAGE(STATUS "=======================================")
  MESSAGE(STATUS "Looking for paul_code... found:")
  MESSAGE(STATUS "* paul_code_DIR          = ${paul_code_DIR}")
  MESSAGE(STATUS "* paul_code_VERSION      = ${paul_code_VERSION}")
  MESSAGE(STATUS "* paul_code_USE_FILE     = ${paul_code_USE_FILE}")

  MESSAGE(STATUS "* paul_code_INCLUDE_DIRS = ${paul_code_INCLUDE_DIRS}")
  MESSAGE(STATUS "* paul_code_LIBRARY_DIRS = ${paul_code_LIBRARY_DIRS}")
  MESSAGE(STATUS "* paul_code_LIBRARIES    = ${paul_code_LIBRARIES}")
ENDIF(CREA_VERBOSE_CMAKE)

# Does the library has an additional config file (user provided) ?
SET(paul_code_HAS_ADDITIONAL_CONFIG_FILE FALSE)

IF (paul_code_HAS_ADDITIONAL_CONFIG_FILE)
  IF(CREA_VERBOSE_CMAKE)
    MESSAGE(STATUS "Reading paul_code additional configuration file")
  ENDIF(CREA_VERBOSE_CMAKE)
  # Include it
  INCLUDE(${paul_code_DIR}/Additionalpaul_codeConfig.cmake)
ENDIF (paul_code_HAS_ADDITIONAL_CONFIG_FILE)
