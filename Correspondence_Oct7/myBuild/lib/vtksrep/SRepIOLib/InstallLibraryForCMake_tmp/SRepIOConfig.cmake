#-----------------------------------------------------------------------------
#
# SRepIOConfig.cmake - CMake configuration file for external projects.
# This file was automatically generated by the cmake macro 
# CREA_INSTALL_LIBRARY_FOR_CMAKE of the package CREA
#
# This file is configured by cmake and used by the
# UseSRepIO.cmake module to load the lib settings 
# for an external project.

# Build tree config ?
SET(CILC_BUILD_TREE_CONFIGURATION FALSE)


IF(UNIX)
SET(GOTO_INSTALL_PREFIX /../../..)
ENDIF(UNIX)


# The SRepIO include file *RELATIVE* directories.
SET(CILC_RELATIVE_INCLUDE_PATHS "include/SRepIO")
# Compute the prefix for include and library paths
IF(CILC_BUILD_TREE_CONFIGURATION)
  # In build tree 
  # the include paths are relative to the source tree *AND* the binary tree 
  # for generated files 
  SET(CILC_INCLUDE_PATH_PREFIX /playpen/workspace/newuoa/Correspondence_Oct7/Pablo2_Oct7)
  # Build the *ABSOLUTE* directories
  FOREACH(path ${CILC_RELATIVE_INCLUDE_PATHS})
    SET(SRepIO_INCLUDE_DIRS
      ${SRepIO_INCLUDE_DIRS}
      ${CILC_INCLUDE_PATH_PREFIX}/${path}
      )
  ENDFOREACH(path ${CILC_RELATIVE_INCLUDE_PATHS})
  SET(CILC_INCLUDE_PATH_PREFIX /playpen/workspace/newuoa/Correspondence_Oct7/myBuild)
  # Build the *ABSOLUTE* directories
  FOREACH(path ${CILC_RELATIVE_INCLUDE_PATHS})
    SET(SRepIO_INCLUDE_DIRS
      ${SRepIO_INCLUDE_DIRS}
      ${CILC_INCLUDE_PATH_PREFIX}/${path}
      )
  ENDFOREACH(path ${CILC_RELATIVE_INCLUDE_PATHS})
ELSE(CILC_BUILD_TREE_CONFIGURATION)
  # In install tree 
  # the include paths are relative to install prefix 
  # On unix , GOTO_INSTALL_PREFIX allows to get back to the 
  # installation prefix from  SRepIO_DIR
  SET(CILC_INCLUDE_PATH_PREFIX ${SRepIO_DIR}${GOTO_INSTALL_PREFIX})
  # Build the *ABSOLUTE* directories
  FOREACH(path ${CILC_RELATIVE_INCLUDE_PATHS})
    SET(SRepIO_INCLUDE_DIRS
      ${SRepIO_INCLUDE_DIRS}
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
  SET(CILC_LIBRARY_PATH_PREFIX ${SRepIO_DIR}${GOTO_INSTALL_PREFIX})
ENDIF(CILC_BUILD_TREE_CONFIGURATION)
# The SRepIO library file *RELATIVE* directories.
SET(CILC_RELATIVE_LIBRARY_PATHS "lib64")
# Build the *ABSOLUTE* directories
FOREACH(path ${CILC_RELATIVE_LIBRARY_PATHS})
  SET(SRepIO_LIBRARY_DIRS
    ${SRepIO_LIBRARY_DIRS}
    ${CILC_LIBRARY_PATH_PREFIX}/${path}
    )
ENDFOREACH(path ${CILC_RELATIVE_LIBRARY_PATHS})

# Set the "prefix path"
SET(SRepIO_INSTALL_PREFIX ${CILC_LIBRARY_PATH_PREFIX})

# The C and C++ flags added by SRepIO to the cmake-configured flags.
SET(SRepIO_REQUIRED_C_FLAGS "")
SET(SRepIO_REQUIRED_CXX_FLAGS "")
SET(SRepIO_REQUIRED_LINK_FLAGS "")

# The SRepIO version 
SET(SRepIO_MAJOR_VERSION )
SET(SRepIO_MINOR_VERSION )
SET(SRepIO_BUILD_VERSION )
SET(SRepIO_VERSION ..)

# The location of the UseSRepIO.cmake file.
SET(SRepIO_USE_FILE "${SRepIO_DIR}/UseSRepIO.cmake")

# The build settings file.
SET(SRepIO_BUILD_SETTINGS_FILE 
  "${SRepIO_DIR}/SRepIOBuildSettings.cmake")

# A list of all libraries for SRepIO.  Those listed here should
# automatically pull in their dependencies.
SET(SRepIO_LIBRARIES SRepIO)

# Messages
IF(CREA_VERBOSE_CMAKE)
  MESSAGE(STATUS "=======================================")
  MESSAGE(STATUS "Looking for SRepIO... found:")
  MESSAGE(STATUS "* SRepIO_DIR          = ${SRepIO_DIR}")
  MESSAGE(STATUS "* SRepIO_VERSION      = ${SRepIO_VERSION}")
  MESSAGE(STATUS "* SRepIO_USE_FILE     = ${SRepIO_USE_FILE}")

  MESSAGE(STATUS "* SRepIO_INCLUDE_DIRS = ${SRepIO_INCLUDE_DIRS}")
  MESSAGE(STATUS "* SRepIO_LIBRARY_DIRS = ${SRepIO_LIBRARY_DIRS}")
  MESSAGE(STATUS "* SRepIO_LIBRARIES    = ${SRepIO_LIBRARIES}")
ENDIF(CREA_VERBOSE_CMAKE)

# Does the library has an additional config file (user provided) ?
SET(SRepIO_HAS_ADDITIONAL_CONFIG_FILE FALSE)

IF (SRepIO_HAS_ADDITIONAL_CONFIG_FILE)
  IF(CREA_VERBOSE_CMAKE)
    MESSAGE(STATUS "Reading SRepIO additional configuration file")
  ENDIF(CREA_VERBOSE_CMAKE)
  # Include it
  INCLUDE(${SRepIO_DIR}/AdditionalSRepIOConfig.cmake)
ENDIF (SRepIO_HAS_ADDITIONAL_CONFIG_FILE)
