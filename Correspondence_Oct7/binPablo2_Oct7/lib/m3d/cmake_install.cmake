# Install script for directory: /work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/m3d

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
  FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib64/m3d" TYPE FILE FILES "/work/ltu/Correspondence_Oct7/binPablo2_Oct7/lib/m3d/InstallLibraryForCMake_tmp/Usem3d.cmake")
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib64/m3d" TYPE FILE FILES "/work/ltu/Correspondence_Oct7/binPablo2_Oct7/lib/m3d/InstallLibraryForCMake_tmp/m3dConfig.cmake")
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib64/m3d" TYPE FILE FILES "/work/ltu/Correspondence_Oct7/binPablo2_Oct7/lib/m3d/InstallLibraryForCMake_tmp/m3dBuildSettings.cmake")
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  list(APPEND CPACK_ABSOLUTE_DESTINATION_FILES
   "/usr/local/share/cmake/Modules/Findm3d.cmake")
FILE(INSTALL DESTINATION "/usr/local/share/cmake/Modules" TYPE FILE FILES "/work/ltu/Correspondence_Oct7/binPablo2_Oct7/lib/m3d/InstallLibraryForCMake_tmp/Findm3d.cmake")
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/m3d" TYPE FILE FILES
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/m3d/include/WorldSystem.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/m3d/include/M3DPrimitiveRenderer.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/m3d/include/M3DEndPrimitive.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/m3d/include/GeodesicDistanceFunction.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/m3d/include/ControlParms.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/m3d/include/Quat.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/m3d/include/Image3D.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/m3d/include/M3DQuadInterpolater.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/m3d/include/M3DFigureTreeNode.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/m3d/include/M3DSpoke.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/m3d/include/M3DQuadEndPrimitive.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/m3d/include/M3DPNSTransform.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/m3d/include/pablo_version.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/m3d/include/M3DTubeEndPrimitive.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/m3d/include/M3DObjectRenderer.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/m3d/include/Bezier2D.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/m3d/include/M3DQuadPrimitive.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/m3d/include/M3DDisplayGlobals.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/m3d/include/M3DPGAPrimitiveStats.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/m3d/include/ScaleTransShapeSpace.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/m3d/include/M3DAtomPredictorTube.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/m3d/include/M3DCPNSStats.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/m3d/include/GaussianBlur3D.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/m3d/include/control_parms_defaults.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/m3d/include/ControlFile.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/m3d/include/VectorND.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/m3d/include/M3DInterpolater.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/m3d/include/M3DTubePrimitive.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/m3d/include/M3DPGA.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/m3d/include/Vector3D.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/m3d/include/M3DTubeFigure.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/m3d/include/ControlParmsAccess.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/m3d/include/SimilarityTransform3D.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/m3d/include/M3DQuadInterpolator.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/m3d/include/M3DFigureStats.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/m3d/include/Bezier1D.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/m3d/include/M3DSpoke_bkp_1.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/m3d/include/IntervalTimer.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/m3d/include/M3DObject.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/m3d/include/SubdivBoundary.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/m3d/include/M3DFigure.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/m3d/include/utility.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/m3d/include/M3DQuadFigureRenderer.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/m3d/include/PluncMatrixFile.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/m3d/include/RigidShapeSpace.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/m3d/include/NormalAlignedShapeSpace.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/m3d/include/M3DQuadFigure.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/m3d/include/GeodesicSym.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/m3d/include/M3DFigureRenderer.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/m3d/include/DQFImage.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/m3d/include/M3DAtomPredictor.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/m3d/include/snapshot.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/m3d/include/Geodesic.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/m3d/include/Mathdefs.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/m3d/include/M3DPrimitive.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/m3d/include/M3DPGAStats.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/m3d/include/M3DObjectFile.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/m3d/include/RotationAndScaleShapeSpace.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/m3d/include/RAWImageFile.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/m3d/include/InterfiguralConstraints.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/m3d/include/Vector2D.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/m3d/include/M3DSpoke_bkp_0.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/m3d/include/ImageResample3D.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/m3d/include/BYU.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/m3d/include/TileSet.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/m3d/include/SimTransShapeSpace.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/m3d/include/M3DAtomPredictorQuad.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/m3d/include/Hermite.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/m3d/include/M3DFigurePredictor.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/m3d/include/M3DTubeFigureRenderer.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/m3d/include/Matrix4D.h"
    "/work/ltu/Correspondence_Oct7/binPablo2_Oct7/lib/m3d/m3d_EXPORT.h"
    )
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib64" TYPE STATIC_LIBRARY FILES "/work/ltu/Correspondence_Oct7/binPablo2_Oct7/libm3d.a")
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

