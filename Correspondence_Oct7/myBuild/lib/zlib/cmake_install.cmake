# Install script for directory: /playpen/workspace/newuoa/Correspondence_Oct7/Pablo2_Oct7/lib/zlib

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

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib64/zlib" TYPE FILE FILES "/playpen/workspace/newuoa/Correspondence_Oct7/myBuild/lib/zlib/InstallLibraryForCMake_tmp/Usezlib.cmake")
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib64/zlib" TYPE FILE FILES "/playpen/workspace/newuoa/Correspondence_Oct7/myBuild/lib/zlib/InstallLibraryForCMake_tmp/zlibConfig.cmake")
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib64/zlib" TYPE FILE FILES "/playpen/workspace/newuoa/Correspondence_Oct7/myBuild/lib/zlib/InstallLibraryForCMake_tmp/zlibBuildSettings.cmake")
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/usr/local/share/cmake/Modules/Findzlib.cmake")
  IF (CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  ENDIF (CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
  IF (CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  ENDIF (CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
FILE(INSTALL DESTINATION "/usr/local/share/cmake/Modules" TYPE FILE FILES "/playpen/workspace/newuoa/Correspondence_Oct7/myBuild/lib/zlib/InstallLibraryForCMake_tmp/Findzlib.cmake")
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/zlib" TYPE FILE FILES
    "/playpen/workspace/newuoa/Correspondence_Oct7/Pablo2_Oct7/lib/zlib/include/inffixed.h"
    "/playpen/workspace/newuoa/Correspondence_Oct7/Pablo2_Oct7/lib/zlib/include/zconf.in.h"
    "/playpen/workspace/newuoa/Correspondence_Oct7/Pablo2_Oct7/lib/zlib/include/trees.h"
    "/playpen/workspace/newuoa/Correspondence_Oct7/Pablo2_Oct7/lib/zlib/include/zlib.h"
    "/playpen/workspace/newuoa/Correspondence_Oct7/Pablo2_Oct7/lib/zlib/include/deflate.h"
    "/playpen/workspace/newuoa/Correspondence_Oct7/Pablo2_Oct7/lib/zlib/include/zconf.h"
    "/playpen/workspace/newuoa/Correspondence_Oct7/Pablo2_Oct7/lib/zlib/include/inflate.h"
    "/playpen/workspace/newuoa/Correspondence_Oct7/Pablo2_Oct7/lib/zlib/include/zutil.h"
    "/playpen/workspace/newuoa/Correspondence_Oct7/Pablo2_Oct7/lib/zlib/include/inftrees.h"
    "/playpen/workspace/newuoa/Correspondence_Oct7/Pablo2_Oct7/lib/zlib/include/crc32.h"
    "/playpen/workspace/newuoa/Correspondence_Oct7/Pablo2_Oct7/lib/zlib/include/inffast.h"
    "/playpen/workspace/newuoa/Correspondence_Oct7/myBuild/lib/zlib/zlib_EXPORT.h"
    )
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib64" TYPE STATIC_LIBRARY FILES "/playpen/workspace/newuoa/Correspondence_Oct7/myBuild/libzlib.a")
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

