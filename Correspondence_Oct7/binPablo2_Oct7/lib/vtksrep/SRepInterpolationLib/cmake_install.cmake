# Install script for directory: /work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/vtksrep/SRepInterpolationLib

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

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib64/SRepInterpolation" TYPE FILE FILES "/work/ltu/Correspondence_Oct7/binPablo2_Oct7/lib/vtksrep/SRepInterpolationLib/InstallLibraryForCMake_tmp/UseSRepInterpolation.cmake")
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib64/SRepInterpolation" TYPE FILE FILES "/work/ltu/Correspondence_Oct7/binPablo2_Oct7/lib/vtksrep/SRepInterpolationLib/InstallLibraryForCMake_tmp/SRepInterpolationConfig.cmake")
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib64/SRepInterpolation" TYPE FILE FILES "/work/ltu/Correspondence_Oct7/binPablo2_Oct7/lib/vtksrep/SRepInterpolationLib/InstallLibraryForCMake_tmp/SRepInterpolationBuildSettings.cmake")
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  list(APPEND CPACK_ABSOLUTE_DESTINATION_FILES
   "/usr/local/share/cmake/Modules/FindSRepInterpolation.cmake")
FILE(INSTALL DESTINATION "/usr/local/share/cmake/Modules" TYPE FILE FILES "/work/ltu/Correspondence_Oct7/binPablo2_Oct7/lib/vtksrep/SRepInterpolationLib/InstallLibraryForCMake_tmp/FindSRepInterpolation.cmake")
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/SRepInterpolation" TYPE FILE FILES
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/vtksrep/SRepInterpolationLib/vtksrepinterpolatecrestspokes.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/vtksrep/SRepInterpolationLib/vtkquadmeshtotriangularmesh.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/vtksrep/SRepInterpolationLib/vtksrepinterpolator.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/vtksrep/SRepInterpolationLib/vtksrepinterpolatemedialsheet.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/vtksrep/SRepInterpolationLib/vtkinterpolatecurve.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/vtksrep/SRepInterpolationLib/vtksrepinterpolatemedialcrestcurve.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/vtksrep/SRepInterpolationLib/minimizecurvaturefunction.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/vtksrep/SRepInterpolationLib/vtksrepinterpolatecrestspokesquartic.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/vtksrep/SRepInterpolationLib/vtksrepinterpolatemedialspokes.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/vtksrep/SRepInterpolationLib/vtksrepinterpolatemedialspokeshermite.h"
    "/work/ltu/Correspondence_Oct7/binPablo2_Oct7/lib/vtksrep/SRepInterpolationLib/SRepInterpolation_EXPORT.h"
    )
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib64" TYPE STATIC_LIBRARY FILES "/work/ltu/Correspondence_Oct7/binPablo2_Oct7/libSRepInterpolation.a")
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

