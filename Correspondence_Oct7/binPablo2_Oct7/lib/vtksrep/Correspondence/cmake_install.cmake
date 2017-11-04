# Install script for directory: /work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/vtksrep/Correspondence

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
  FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib64/Correspondence" TYPE FILE FILES "/work/ltu/Correspondence_Oct7/binPablo2_Oct7/lib/vtksrep/Correspondence/InstallLibraryForCMake_tmp/UseCorrespondence.cmake")
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib64/Correspondence" TYPE FILE FILES "/work/ltu/Correspondence_Oct7/binPablo2_Oct7/lib/vtksrep/Correspondence/InstallLibraryForCMake_tmp/CorrespondenceConfig.cmake")
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib64/Correspondence" TYPE FILE FILES "/work/ltu/Correspondence_Oct7/binPablo2_Oct7/lib/vtksrep/Correspondence/InstallLibraryForCMake_tmp/CorrespondenceBuildSettings.cmake")
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  list(APPEND CPACK_ABSOLUTE_DESTINATION_FILES
   "/usr/local/share/cmake/Modules/FindCorrespondence.cmake")
FILE(INSTALL DESTINATION "/usr/local/share/cmake/Modules" TYPE FILE FILES "/work/ltu/Correspondence_Oct7/binPablo2_Oct7/lib/vtksrep/Correspondence/InstallLibraryForCMake_tmp/FindCorrespondence.cmake")
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/Correspondence" TYPE FILE FILES
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/vtksrep/Correspondence/visualization.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/vtksrep/Correspondence/macro_const_values.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/vtksrep/Correspondence/procrustes.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/vtksrep/Correspondence/thinplatesplinepdmtosrep.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/vtksrep/Correspondence/slidestandardspokes.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/vtksrep/Correspondence/itkMesh3DProcrustesAlignFilter.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/vtksrep/Correspondence/oneplusonecostfunction.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/vtksrep/Correspondence/samplesreppoints.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/vtksrep/Correspondence/slidecrestspokes.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/vtksrep/Correspondence/cost_function.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/vtksrep/Correspondence/quadfigattribution.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/vtksrep/Correspondence/optimizationusingnewuoa.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/vtksrep/Correspondence/regularityentropy.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/vtksrep/Correspondence/movespokes.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/vtksrep/Correspondence/newuoa.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/vtksrep/Correspondence/weightedprocrustes.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/vtksrep/Correspondence/srepprocrustes.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/vtksrep/Correspondence/srepCorrespondenceEvaluator.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/vtksrep/Correspondence/alignsrep.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/vtksrep/Correspondence/calcrestregularityfeatures.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/vtksrep/Correspondence/itkCorrespondenceEvaluator.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/vtksrep/Correspondence/thinplatesplinesrep.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/vtksrep/Correspondence/calregularityfeatures.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/vtksrep/Correspondence/correspondenceevaluation.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/vtksrep/Correspondence/uvmap.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/vtksrep/Correspondence/matlabengine.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/vtksrep/Correspondence/toolsfunc.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/vtksrep/Correspondence/shiftedsrep.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/vtksrep/Correspondence/objectivefunction.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/vtksrep/Correspondence/srepCorrespondenceEvaluator.txx"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/vtksrep/Correspondence/itkCorrespondenceEvaluator.txx"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/vtksrep/Correspondence/itkMesh3DProcrustesAlignFilter.txx"
    "/work/ltu/Correspondence_Oct7/binPablo2_Oct7/lib/vtksrep/Correspondence/Correspondence_EXPORT.h"
    )
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib64" TYPE STATIC_LIBRARY FILES "/work/ltu/Correspondence_Oct7/binPablo2_Oct7/libCorrespondence.a")
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

