# Install script for directory: /work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/register

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
  FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib64/register" TYPE FILE FILES "/work/ltu/Correspondence_Oct7/binPablo2_Oct7/lib/register/InstallLibraryForCMake_tmp/Useregister.cmake")
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib64/register" TYPE FILE FILES "/work/ltu/Correspondence_Oct7/binPablo2_Oct7/lib/register/InstallLibraryForCMake_tmp/registerConfig.cmake")
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib64/register" TYPE FILE FILES "/work/ltu/Correspondence_Oct7/binPablo2_Oct7/lib/register/InstallLibraryForCMake_tmp/registerBuildSettings.cmake")
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  list(APPEND CPACK_ABSOLUTE_DESTINATION_FILES
   "/usr/local/share/cmake/Modules/Findregister.cmake")
FILE(INSTALL DESTINATION "/usr/local/share/cmake/Modules" TYPE FILE FILES "/work/ltu/Correspondence_Oct7/binPablo2_Oct7/lib/register/InstallLibraryForCMake_tmp/Findregister.cmake")
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/register" TYPE FILE FILES
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/register/include/Trackball.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/register/include/M3DSubfigureOptimizer.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/register/include/M3DMainFigResiduePGAProblem.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/register/include/M3DSpokeOptimizer.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/register/include/OctTreeTauBand.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/register/include/M3DSpokeProblem.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/register/include/MakeOptimizationVideo.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/register/include/M3DMainFigureOptimizer.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/register/include/M3DFigureElongater.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/register/include/M3DBoundingSphere.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/register/include/M3DRegistrationPGAProblem.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/register/include/M3DSubfigureTransformation.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/register/include/FitUnlabeledPoints.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/register/include/M3DVoxelOptimizer.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/register/include/M3DRegistrationSimilarityProblem.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/register/include/DistanceToPointSetFunction.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/register/include/M3DDeformationProblem.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/register/include/M3DtoPovray.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/register/include/M3DSubfigureProblem.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/register/include/M3DSimilarityCPNSOptimizer.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/register/include/M3DSimilarityCPNSProblem.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/register/include/FunctionExplorer.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/register/include/M3DDeformationOptimizer.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/register/include/M3DRegistrationPGAOptimizer.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/register/include/ImagePlanes.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/register/include/P3DControlBkp1.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/register/include/M3DMainFigureProblem.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/register/include/M3DSimilarityElongationProblem.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/register/include/P3DUndoList.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/register/include/M3DPrimitiveCorrector.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/register/include/M3DRegistrationOptimizer.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/register/include/PCA.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/register/include/P3DControl.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/register/include/M3DSimilarityPGAOptimizer.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/register/include/M3DMainFigResidueCPNSOptimizer.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/register/include/M3DRegistrationSimilarityOptimizer.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/register/include/M3DSimilarityPGAProblem.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/register/include/M3DRegistrationProblem.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/register/include/M3DMainFigResidueCPNSProblem.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/register/include/M3DTubeSimilarityPGAProblem.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/register/include/M3DAdaptiveRegistrationPGAProblem.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/register/include/LandmarkDeformation.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/register/include/M3DAdaptiveRegistrationPGAOptimizer.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/register/include/OptimizerBase.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/register/include/M3DSRepOptimizer.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/register/include/M3DMainFigResiduePGAOptimizer.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/register/include/Anastruct.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/register/include/M3DSimilarityElongationOptimizer.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/register/include/M3DSRepProblem.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/register/include/M3DDeformationPGAProblem.h"
    "/work/ltu/Correspondence_Oct7/binPablo2_Oct7/lib/register/register_EXPORT.h"
    )
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib64" TYPE STATIC_LIBRARY FILES "/work/ltu/Correspondence_Oct7/binPablo2_Oct7/libregister.a")
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

