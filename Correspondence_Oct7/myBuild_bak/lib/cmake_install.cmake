# Install script for directory: /playpen/workspace/newuoa/Correspondence_Oct7/Pablo2_Oct7/lib

# Set the install prefix
IF(NOT DEFINED CMAKE_INSTALL_PREFIX)
  SET(CMAKE_INSTALL_PREFIX "/usr/local")
ENDIF(NOT DEFINED CMAKE_INSTALL_PREFIX)
STRING(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
IF(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  IF(BUILD_TYPE)
    STRING(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  ELSE(BUILD_TYPE)
    SET(CMAKE_INSTALL_CONFIG_NAME "")
  ENDIF(BUILD_TYPE)
  MESSAGE(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
ENDIF(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)

# Set the component getting installed.
IF(NOT CMAKE_INSTALL_COMPONENT)
  IF(COMPONENT)
    MESSAGE(STATUS "Install component: \"${COMPONENT}\"")
    SET(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  ELSE(COMPONENT)
    SET(CMAKE_INSTALL_COMPONENT)
  ENDIF(COMPONENT)
ENDIF(NOT CMAKE_INSTALL_COMPONENT)

# Install shared libraries without execute permission?
IF(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  SET(CMAKE_INSTALL_SO_NO_EXE "1")
ENDIF(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)

IF(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  INCLUDE("/playpen/workspace/newuoa/Correspondence_Oct7/myBuild/lib/gui/cmake_install.cmake")
  INCLUDE("/playpen/workspace/newuoa/Correspondence_Oct7/myBuild/lib/ImageIO/cmake_install.cmake")
  INCLUDE("/playpen/workspace/newuoa/Correspondence_Oct7/myBuild/lib/m3d/cmake_install.cmake")
  INCLUDE("/playpen/workspace/newuoa/Correspondence_Oct7/myBuild/lib/match/cmake_install.cmake")
  INCLUDE("/playpen/workspace/newuoa/Correspondence_Oct7/myBuild/lib/paul_code/cmake_install.cmake")
  INCLUDE("/playpen/workspace/newuoa/Correspondence_Oct7/myBuild/lib/planes/cmake_install.cmake")
  INCLUDE("/playpen/workspace/newuoa/Correspondence_Oct7/myBuild/lib/register/cmake_install.cmake")
  INCLUDE("/playpen/workspace/newuoa/Correspondence_Oct7/myBuild/lib/seurat/cmake_install.cmake")
  INCLUDE("/playpen/workspace/newuoa/Correspondence_Oct7/myBuild/lib/version/cmake_install.cmake")
  INCLUDE("/playpen/workspace/newuoa/Correspondence_Oct7/myBuild/lib/zlib/cmake_install.cmake")
  INCLUDE("/playpen/workspace/newuoa/Correspondence_Oct7/myBuild/lib/flvw-1.0/cmake_install.cmake")
  INCLUDE("/playpen/workspace/newuoa/Correspondence_Oct7/myBuild/lib/vtksrep/cmake_install.cmake")

ENDIF(NOT CMAKE_INSTALL_LOCAL_ONLY)

