# Install script for directory: /work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib

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
    SET(CMAKE_INSTALL_CONFIG_NAME "relwithdebinfo")
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
  SET(CMAKE_INSTALL_SO_NO_EXE "0")
ENDIF(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)

IF(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  INCLUDE("/work/ltu/Correspondence_Oct7/binPablo2_Oct7/lib/gui/cmake_install.cmake")
  INCLUDE("/work/ltu/Correspondence_Oct7/binPablo2_Oct7/lib/ImageIO/cmake_install.cmake")
  INCLUDE("/work/ltu/Correspondence_Oct7/binPablo2_Oct7/lib/m3d/cmake_install.cmake")
  INCLUDE("/work/ltu/Correspondence_Oct7/binPablo2_Oct7/lib/match/cmake_install.cmake")
  INCLUDE("/work/ltu/Correspondence_Oct7/binPablo2_Oct7/lib/paul_code/cmake_install.cmake")
  INCLUDE("/work/ltu/Correspondence_Oct7/binPablo2_Oct7/lib/planes/cmake_install.cmake")
  INCLUDE("/work/ltu/Correspondence_Oct7/binPablo2_Oct7/lib/register/cmake_install.cmake")
  INCLUDE("/work/ltu/Correspondence_Oct7/binPablo2_Oct7/lib/seurat/cmake_install.cmake")
  INCLUDE("/work/ltu/Correspondence_Oct7/binPablo2_Oct7/lib/version/cmake_install.cmake")
  INCLUDE("/work/ltu/Correspondence_Oct7/binPablo2_Oct7/lib/zlib/cmake_install.cmake")
  INCLUDE("/work/ltu/Correspondence_Oct7/binPablo2_Oct7/lib/flvw-1.0/cmake_install.cmake")
  INCLUDE("/work/ltu/Correspondence_Oct7/binPablo2_Oct7/lib/vtksrep/cmake_install.cmake")

ENDIF(NOT CMAKE_INSTALL_LOCAL_ONLY)

