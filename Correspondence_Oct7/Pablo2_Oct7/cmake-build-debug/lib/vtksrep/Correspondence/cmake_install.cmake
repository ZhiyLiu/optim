# Install script for directory: /playpen/workspace/newuoa/Correspondence_Oct7/Pablo2_Oct7/lib/vtksrep/Correspondence

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
    SET(CMAKE_INSTALL_CONFIG_NAME "Debug")
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
  FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib64/Correspondence" TYPE FILE FILES "/playpen/workspace/newuoa/Correspondence_Oct7/Pablo2_Oct7/cmake-build-debug/lib/vtksrep/Correspondence/InstallLibraryForCMake_tmp/UseCorrespondence.cmake")
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib64/Correspondence" TYPE FILE FILES "/playpen/workspace/newuoa/Correspondence_Oct7/Pablo2_Oct7/cmake-build-debug/lib/vtksrep/Correspondence/InstallLibraryForCMake_tmp/CorrespondenceConfig.cmake")
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib64/Correspondence" TYPE FILE FILES "/playpen/workspace/newuoa/Correspondence_Oct7/Pablo2_Oct7/cmake-build-debug/lib/vtksrep/Correspondence/InstallLibraryForCMake_tmp/CorrespondenceBuildSettings.cmake")
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/usr/local/share/cmake/Modules/FindCorrespondence.cmake")
  IF (CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  ENDIF (CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
  IF (CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  ENDIF (CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
FILE(INSTALL DESTINATION "/usr/local/share/cmake/Modules" TYPE FILE FILES "/playpen/workspace/newuoa/Correspondence_Oct7/Pablo2_Oct7/cmake-build-debug/lib/vtksrep/Correspondence/InstallLibraryForCMake_tmp/FindCorrespondence.cmake")
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/Correspondence" TYPE FILE FILES
    "/playpen/workspace/newuoa/Correspondence_Oct7/Pablo2_Oct7/lib/vtksrep/Correspondence/visualization.h"
    "/playpen/workspace/newuoa/Correspondence_Oct7/Pablo2_Oct7/lib/vtksrep/Correspondence/samplesreppoints.h"
    "/playpen/workspace/newuoa/Correspondence_Oct7/Pablo2_Oct7/lib/vtksrep/Correspondence/matlabengine.h"
    "/playpen/workspace/newuoa/Correspondence_Oct7/Pablo2_Oct7/lib/vtksrep/Correspondence/correspondenceevaluation.h"
    "/playpen/workspace/newuoa/Correspondence_Oct7/Pablo2_Oct7/lib/vtksrep/Correspondence/thinplatesplinepdmtosrep.h"
    "/playpen/workspace/newuoa/Correspondence_Oct7/Pablo2_Oct7/lib/vtksrep/Correspondence/M3DSpokeAngleOptimizer.h"
    "/playpen/workspace/newuoa/Correspondence_Oct7/Pablo2_Oct7/lib/vtksrep/Correspondence/srepCorrespondenceEvaluator.h"
    "/playpen/workspace/newuoa/Correspondence_Oct7/Pablo2_Oct7/lib/vtksrep/Correspondence/weightedprocrustes.h"
    "/playpen/workspace/newuoa/Correspondence_Oct7/Pablo2_Oct7/lib/vtksrep/Correspondence/calregularityfeatures.h"
    "/playpen/workspace/newuoa/Correspondence_Oct7/Pablo2_Oct7/lib/vtksrep/Correspondence/thinplatesplinesrep.h"
    "/playpen/workspace/newuoa/Correspondence_Oct7/Pablo2_Oct7/lib/vtksrep/Correspondence/cost_function.h"
    "/playpen/workspace/newuoa/Correspondence_Oct7/Pablo2_Oct7/lib/vtksrep/Correspondence/RegularityEntropy.h"
    "/playpen/workspace/newuoa/Correspondence_Oct7/Pablo2_Oct7/lib/vtksrep/Correspondence/itkMesh3DProcrustesAlignFilter.h"
    "/playpen/workspace/newuoa/Correspondence_Oct7/Pablo2_Oct7/lib/vtksrep/Correspondence/NormalMatchComputer.h"
    "/playpen/workspace/newuoa/Correspondence_Oct7/Pablo2_Oct7/lib/vtksrep/Correspondence/slidecrestspokes.h"
    "/playpen/workspace/newuoa/Correspondence_Oct7/Pablo2_Oct7/lib/vtksrep/Correspondence/srepprocrustes.h"
    "/playpen/workspace/newuoa/Correspondence_Oct7/Pablo2_Oct7/lib/vtksrep/Correspondence/calcrestregularityfeatures.h"
    "/playpen/workspace/newuoa/Correspondence_Oct7/Pablo2_Oct7/lib/vtksrep/Correspondence/uvmap.h"
    "/playpen/workspace/newuoa/Correspondence_Oct7/Pablo2_Oct7/lib/vtksrep/Correspondence/quaternion_utils.h"
    "/playpen/workspace/newuoa/Correspondence_Oct7/Pablo2_Oct7/lib/vtksrep/Correspondence/slidestandardspokes.h"
    "/playpen/workspace/newuoa/Correspondence_Oct7/Pablo2_Oct7/lib/vtksrep/Correspondence/toolsfunc.h"
    "/playpen/workspace/newuoa/Correspondence_Oct7/Pablo2_Oct7/lib/vtksrep/Correspondence/quadfigattribution.h"
    "/playpen/workspace/newuoa/Correspondence_Oct7/Pablo2_Oct7/lib/vtksrep/Correspondence/M3DSpokeLengthOptimizer2D.h"
    "/playpen/workspace/newuoa/Correspondence_Oct7/Pablo2_Oct7/lib/vtksrep/Correspondence/itkCorrespondenceEvaluator.h"
    "/playpen/workspace/newuoa/Correspondence_Oct7/Pablo2_Oct7/lib/vtksrep/Correspondence/procrustes.h"
    "/playpen/workspace/newuoa/Correspondence_Oct7/Pablo2_Oct7/lib/vtksrep/Correspondence/SimilarityComputer2D.h"
    "/playpen/workspace/newuoa/Correspondence_Oct7/Pablo2_Oct7/lib/vtksrep/Correspondence/regularityentropy.h"
    "/playpen/workspace/newuoa/Correspondence_Oct7/Pablo2_Oct7/lib/vtksrep/Correspondence/objectivefunction.h"
    "/playpen/workspace/newuoa/Correspondence_Oct7/Pablo2_Oct7/lib/vtksrep/Correspondence/SimilarityComputer.h"
    "/playpen/workspace/newuoa/Correspondence_Oct7/Pablo2_Oct7/lib/vtksrep/Correspondence/macro_const_values.h"
    "/playpen/workspace/newuoa/Correspondence_Oct7/Pablo2_Oct7/lib/vtksrep/Correspondence/oneplusonecostfunction.h"
    "/playpen/workspace/newuoa/Correspondence_Oct7/Pablo2_Oct7/lib/vtksrep/Correspondence/movespokes.h"
    "/playpen/workspace/newuoa/Correspondence_Oct7/Pablo2_Oct7/lib/vtksrep/Correspondence/quaternion_io.h"
    "/playpen/workspace/newuoa/Correspondence_Oct7/Pablo2_Oct7/lib/vtksrep/Correspondence/newuoa.h"
    "/playpen/workspace/newuoa/Correspondence_Oct7/Pablo2_Oct7/lib/vtksrep/Correspondence/optimizationusingnewuoa.h"
    "/playpen/workspace/newuoa/Correspondence_Oct7/Pablo2_Oct7/lib/vtksrep/Correspondence/alignsrep.h"
    "/playpen/workspace/newuoa/Correspondence_Oct7/Pablo2_Oct7/lib/vtksrep/Correspondence/shiftedsrep.h"
    "/playpen/workspace/newuoa/Correspondence_Oct7/Pablo2_Oct7/lib/vtksrep/Correspondence/M3DSpokeLengthOptimizer.h"
    "/playpen/workspace/newuoa/Correspondence_Oct7/Pablo2_Oct7/lib/vtksrep/Correspondence/quaternion.h"
    "/playpen/workspace/newuoa/Correspondence_Oct7/Pablo2_Oct7/lib/vtksrep/Correspondence/srepCorrespondenceEvaluator.txx"
    "/playpen/workspace/newuoa/Correspondence_Oct7/Pablo2_Oct7/lib/vtksrep/Correspondence/itkMesh3DProcrustesAlignFilter.txx"
    "/playpen/workspace/newuoa/Correspondence_Oct7/Pablo2_Oct7/lib/vtksrep/Correspondence/itkCorrespondenceEvaluator.txx"
    "/playpen/workspace/newuoa/Correspondence_Oct7/Pablo2_Oct7/cmake-build-debug/lib/vtksrep/Correspondence/Correspondence_EXPORT.h"
    )
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib64" TYPE STATIC_LIBRARY FILES "/playpen/workspace/newuoa/Correspondence_Oct7/Pablo2_Oct7/cmake-build-debug/libCorrespondence.a")
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

