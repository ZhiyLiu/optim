# Install script for directory: /work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/match

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
  FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib64/match" TYPE FILE FILES "/work/ltu/Correspondence_Oct7/binPablo2_Oct7/lib/match/InstallLibraryForCMake_tmp/Usematch.cmake")
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib64/match" TYPE FILE FILES "/work/ltu/Correspondence_Oct7/binPablo2_Oct7/lib/match/InstallLibraryForCMake_tmp/matchConfig.cmake")
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib64/match" TYPE FILE FILES "/work/ltu/Correspondence_Oct7/binPablo2_Oct7/lib/match/InstallLibraryForCMake_tmp/matchBuildSettings.cmake")
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  list(APPEND CPACK_ABSOLUTE_DESTINATION_FILES
   "/usr/local/share/cmake/Modules/Findmatch.cmake")
FILE(INSTALL DESTINATION "/usr/local/share/cmake/Modules" TYPE FILE FILES "/work/ltu/Correspondence_Oct7/binPablo2_Oct7/lib/match/InstallLibraryForCMake_tmp/Findmatch.cmake")
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/match" TYPE FILE FILES
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/match/include/TemplateProfiles.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/match/include/Mask.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/match/include/Tuning.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/match/include/ImageDistanceMapBkp1.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/match/include/MatchBkp2.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/match/include/MatchBkp3.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/match/include/SimpleMaskFile.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/match/include/ImageDistanceMap.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/match/include/MethodOfMoments.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/match/include/ccl.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/match/include/Curvature.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/match/include/gpTuning.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/match/include/SurfacePatchEnsemble.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/match/include/bpTuningBkp2.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/match/include/Danielsson.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/match/include/ObjectRelativeSampling.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/match/include/bpTuningBkp1.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/match/include/MatchBkp1.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/match/include/MatchBkp4.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/match/include/bpTuning.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/match/include/MaskFile.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/match/include/DistanceVectorList.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/match/include/MatchUtility.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/match/include/DistanceMap3D.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/match/include/DQFMatch.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/match/include/deltas.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/match/include/Match.h"
    "/work/ltu/Correspondence_Oct7/binPablo2_Oct7/lib/match/match_EXPORT.h"
    )
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib64" TYPE STATIC_LIBRARY FILES "/work/ltu/Correspondence_Oct7/binPablo2_Oct7/libmatch.a")
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

