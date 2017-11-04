# Install script for directory: /work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/paul_code

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
  FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib64/paul_code" TYPE FILE FILES "/work/ltu/Correspondence_Oct7/binPablo2_Oct7/lib/paul_code/InstallLibraryForCMake_tmp/Usepaul_code.cmake")
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib64/paul_code" TYPE FILE FILES "/work/ltu/Correspondence_Oct7/binPablo2_Oct7/lib/paul_code/InstallLibraryForCMake_tmp/paul_codeConfig.cmake")
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib64/paul_code" TYPE FILE FILES "/work/ltu/Correspondence_Oct7/binPablo2_Oct7/lib/paul_code/InstallLibraryForCMake_tmp/paul_codeBuildSettings.cmake")
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  list(APPEND CPACK_ABSOLUTE_DESTINATION_FILES
   "/usr/local/share/cmake/Modules/Findpaul_code.cmake")
FILE(INSTALL DESTINATION "/usr/local/share/cmake/Modules" TYPE FILE FILES "/work/ltu/Correspondence_Oct7/binPablo2_Oct7/lib/paul_code/InstallLibraryForCMake_tmp/Findpaul_code.cmake")
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/paul_code" TYPE FILE FILES
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/paul_code/include/libs.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/paul_code/include/Registry.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/paul_code/include/optima.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/paul_code/include/matrix.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/paul_code/include/LogManager.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/paul_code/include/ConjugateGradientMethod.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/paul_code/include/problems.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/paul_code/include/SimplexMethod.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/paul_code/include/BrentLinearMethod.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/paul_code/include/support.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/paul_code/include/f2c.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/paul_code/include/Classes.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/paul_code/include/EvolutionaryStrategy.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/paul_code/include/Solution.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/paul_code/include/tsp.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/paul_code/include/simplex.h"
    "/work/ltu/Correspondence_Oct7/binPablo2_Oct7/lib/paul_code/paul_code_EXPORT.h"
    )
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib64" TYPE STATIC_LIBRARY FILES "/work/ltu/Correspondence_Oct7/binPablo2_Oct7/libpaul_code.a")
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

