# Install script for directory: /work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/vtksrep

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
  INCLUDE("/work/ltu/Correspondence_Oct7/binPablo2_Oct7/lib/vtksrep/SRepInterpolationLib/cmake_install.cmake")
  INCLUDE("/work/ltu/Correspondence_Oct7/binPablo2_Oct7/lib/vtksrep/SRepLib/cmake_install.cmake")
  INCLUDE("/work/ltu/Correspondence_Oct7/binPablo2_Oct7/lib/vtksrep/SRepIOLib/cmake_install.cmake")
  INCLUDE("/work/ltu/Correspondence_Oct7/binPablo2_Oct7/lib/vtksrep/SRepVisualizationLib/cmake_install.cmake")
  INCLUDE("/work/ltu/Correspondence_Oct7/binPablo2_Oct7/lib/vtksrep/calquantum/cmake_install.cmake")
  INCLUDE("/work/ltu/Correspondence_Oct7/binPablo2_Oct7/lib/vtksrep/Correspondence/cmake_install.cmake")
  INCLUDE("/work/ltu/Correspondence_Oct7/binPablo2_Oct7/lib/vtksrep/Visualization/cmake_install.cmake")
  INCLUDE("/work/ltu/Correspondence_Oct7/binPablo2_Oct7/lib/vtksrep/vtkpowercrust/cmake_install.cmake")

ENDIF(NOT CMAKE_INSTALL_LOCAL_ONLY)

