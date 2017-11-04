# Install script for directory: /work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/seurat

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
  FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib64/seurat" TYPE FILE FILES "/work/ltu/Correspondence_Oct7/binPablo2_Oct7/lib/seurat/InstallLibraryForCMake_tmp/Useseurat.cmake")
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib64/seurat" TYPE FILE FILES "/work/ltu/Correspondence_Oct7/binPablo2_Oct7/lib/seurat/InstallLibraryForCMake_tmp/seuratConfig.cmake")
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib64/seurat" TYPE FILE FILES "/work/ltu/Correspondence_Oct7/binPablo2_Oct7/lib/seurat/InstallLibraryForCMake_tmp/seuratBuildSettings.cmake")
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  list(APPEND CPACK_ABSOLUTE_DESTINATION_FILES
   "/usr/local/share/cmake/Modules/Findseurat.cmake")
FILE(INSTALL DESTINATION "/usr/local/share/cmake/Modules" TYPE FILE FILES "/work/ltu/Correspondence_Oct7/binPablo2_Oct7/lib/seurat/InstallLibraryForCMake_tmp/Findseurat.cmake")
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/seurat" TYPE FILE FILES
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/seurat/include/Mesh.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/seurat/include/MyQueue.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/seurat/include/CCSubdivsurf.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/seurat/include/renderDefinitions.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/seurat/include/Intersection.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/seurat/include/Zerofinder.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/seurat/include/Plist_subdivcomp.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/seurat/include/HanInterpolation.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/seurat/include/PseudoSet.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/seurat/include/QuadMesh.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/seurat/include/Pointlist_serverB.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/seurat/include/MyList.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/seurat/include/Tritri.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/seurat/include/Conjgrad2.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/seurat/include/LinAlg.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/seurat/include/TubeMesh.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/seurat/include/Shapeheader.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/seurat/include/Xferlist.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/seurat/include/Samplestats.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/seurat/include/Pointlist_server2.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/seurat/include/Shapedepend.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/seurat/include/M3DObjectSurfaceRenderer.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/seurat/include/TileSetRenderer.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/seurat/include/M3DBlendedRenderer.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/seurat/include/SurfaceColorMap.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/seurat/include/Diatomgrid.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/seurat/include/SelectedPartialFigures.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/seurat/include/Subdivsurf.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/seurat/include/M3DObjectSurfaceVisualizer.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/seurat/include/Diatom.h"
    "/work/ltu/Correspondence_Oct7/binPablo2_Oct7/lib/seurat/seurat_EXPORT.h"
    )
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib64" TYPE STATIC_LIBRARY FILES "/work/ltu/Correspondence_Oct7/binPablo2_Oct7/libseurat.a")
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

