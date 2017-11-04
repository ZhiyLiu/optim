#-----------------------------------------------------------------------------
#
# zlibConfig.cmake - CMake configuration file for external projects.
# This file was automatically generated by the cmake macro 
# CREA_INSTALL_LIBRARY_FOR_CMAKE of the package CREA
#
# This file is configured by cmake and used by the
# Usezlib.cmake module to load the lib settings 
# for an external project.

# Build tree config ?
SET(CILC_BUILD_TREE_CONFIGURATION TRUE)


IF(UNIX)
SET(GOTO_INSTALL_PREFIX /../../..)
ENDIF(UNIX)


# The zlib include file *RELATIVE* directories.
SET(CILC_RELATIVE_INCLUDE_PATHS "lib/zlib")
# Compute the prefix for include and library paths
IF(CILC_BUILD_TREE_CONFIGURATION)
  # In build tree 
  # the include paths are relative to the source tree *AND* the binary tree 
  # for generated files 
  SET(CILC_INCLUDE_PATH_PREFIX /work/ltu/Correspondence_Oct7/Pablo2_Oct7)
  # Build the *ABSOLUTE* directories
  FOREACH(path ${CILC_RELATIVE_INCLUDE_PATHS})
    SET(zlib_INCLUDE_DIRS
      ${zlib_INCLUDE_DIRS}
      ${CILC_INCLUDE_PATH_PREFIX}/${path}
      )
  ENDFOREACH(path ${CILC_RELATIVE_INCLUDE_PATHS})
  SET(CILC_INCLUDE_PATH_PREFIX /work/ltu/Correspondence_Oct7/binPablo2_Oct7)
  # Build the *ABSOLUTE* directories
  FOREACH(path ${CILC_RELATIVE_INCLUDE_PATHS})
    SET(zlib_INCLUDE_DIRS
      ${zlib_INCLUDE_DIRS}
      ${CILC_INCLUDE_PATH_PREFIX}/${path}
      )
  ENDFOREACH(path ${CILC_RELATIVE_INCLUDE_PATHS})
ELSE(CILC_BUILD_TREE_CONFIGURATION)
  # In install tree 
  # the include paths are relative to install prefix 
  # On unix , GOTO_INSTALL_PREFIX allows to get back to the 
  # installation prefix from  zlib_DIR
  SET(CILC_INCLUDE_PATH_PREFIX ${zlib_DIR}${GOTO_INSTALL_PREFIX})
  # Build the *ABSOLUTE* directories
  FOREACH(path ${CILC_RELATIVE_INCLUDE_PATHS})
    SET(zlib_INCLUDE_DIRS
      ${zlib_INCLUDE_DIRS}
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
  SET(CILC_LIBRARY_PATH_PREFIX ${zlib_DIR}${GOTO_INSTALL_PREFIX})
ENDIF(CILC_BUILD_TREE_CONFIGURATION)
# The zlib library file *RELATIVE* directories.
SET(CILC_RELATIVE_LIBRARY_PATHS ".")
# Build the *ABSOLUTE* directories
FOREACH(path ${CILC_RELATIVE_LIBRARY_PATHS})
  SET(zlib_LIBRARY_DIRS
    ${zlib_LIBRARY_DIRS}
    ${CILC_LIBRARY_PATH_PREFIX}/${path}
    )
ENDFOREACH(path ${CILC_RELATIVE_LIBRARY_PATHS})

# Set the "prefix path"
SET(zlib_INSTALL_PREFIX ${CILC_LIBRARY_PATH_PREFIX})

# The C and C++ flags added by zlib to the cmake-configured flags.
SET(zlib_REQUIRED_C_FLAGS "")
SET(zlib_REQUIRED_CXX_FLAGS "")
SET(zlib_REQUIRED_LINK_FLAGS "")

# The zlib version 
SET(zlib_MAJOR_VERSION )
SET(zlib_MINOR_VERSION )
SET(zlib_BUILD_VERSION )
SET(zlib_VERSION ..)

# The location of the Usezlib.cmake file.
SET(zlib_USE_FILE "${zlib_DIR}/Usezlib.cmake")

# The build settings file.
SET(zlib_BUILD_SETTINGS_FILE 
  "${zlib_DIR}/zlibBuildSettings.cmake")

# A list of all libraries for zlib.  Those listed here should
# automatically pull in their dependencies.
SET(zlib_LIBRARIES zlib)

# Messages
IF(CREA_VERBOSE_CMAKE)
  MESSAGE(STATUS "=======================================")
  MESSAGE(STATUS "Looking for zlib... found:")
  MESSAGE(STATUS "* zlib_DIR          = ${zlib_DIR}")
  MESSAGE(STATUS "* zlib_VERSION      = ${zlib_VERSION}")
  MESSAGE(STATUS "* zlib_USE_FILE     = ${zlib_USE_FILE}")

  MESSAGE(STATUS "* zlib_INCLUDE_DIRS = ${zlib_INCLUDE_DIRS}")
  MESSAGE(STATUS "* zlib_LIBRARY_DIRS = ${zlib_LIBRARY_DIRS}")
  MESSAGE(STATUS "* zlib_LIBRARIES    = ${zlib_LIBRARIES}")
ENDIF(CREA_VERBOSE_CMAKE)

# Does the library has an additional config file (user provided) ?
SET(zlib_HAS_ADDITIONAL_CONFIG_FILE FALSE)

IF (zlib_HAS_ADDITIONAL_CONFIG_FILE)
  IF(CREA_VERBOSE_CMAKE)
    MESSAGE(STATUS "Reading zlib additional configuration file")
  ENDIF(CREA_VERBOSE_CMAKE)
  # Include it
  INCLUDE(${zlib_DIR}/AdditionalzlibConfig.cmake)
ENDIF (zlib_HAS_ADDITIONAL_CONFIG_FILE)