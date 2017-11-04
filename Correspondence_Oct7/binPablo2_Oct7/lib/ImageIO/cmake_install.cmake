# Install script for directory: /work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/ImageIO

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
  FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib64/ImageIO" TYPE FILE FILES "/work/ltu/Correspondence_Oct7/binPablo2_Oct7/lib/ImageIO/InstallLibraryForCMake_tmp/UseImageIO.cmake")
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib64/ImageIO" TYPE FILE FILES "/work/ltu/Correspondence_Oct7/binPablo2_Oct7/lib/ImageIO/InstallLibraryForCMake_tmp/ImageIOConfig.cmake")
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib64/ImageIO" TYPE FILE FILES "/work/ltu/Correspondence_Oct7/binPablo2_Oct7/lib/ImageIO/InstallLibraryForCMake_tmp/ImageIOBuildSettings.cmake")
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  list(APPEND CPACK_ABSOLUTE_DESTINATION_FILES
   "/usr/local/share/cmake/Modules/FindImageIO.cmake")
FILE(INSTALL DESTINATION "/usr/local/share/cmake/Modules" TYPE FILE FILES "/work/ltu/Correspondence_Oct7/binPablo2_Oct7/lib/ImageIO/InstallLibraryForCMake_tmp/FindImageIO.cmake")
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/ImageIO" TYPE FILE FILES
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/ImageIO/utilities/include/BinaryIO.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/ImageIO/utilities/include/Debugging.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/ImageIO/utilities/include/BasicException.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/ImageIO/Analyze/include/mayo_analyze.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/ImageIO/Analyze/include/analyze_io.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/ImageIO/Analyze/include/Analyze.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/ImageIO/Gipl/include/ipmatrix.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/ImageIO/Gipl/include/ipimage.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/ImageIO/Gipl/include/giplio.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/ImageIO/Gipl/include/gipl_header.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/ImageIO/Gipl/include/misc.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/ImageIO/Meta/include/metaImageTypes.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/ImageIO/Meta/include/metaUtils.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/ImageIO/Meta/include/metaEvent.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/ImageIO/Meta/include/metaImageUtils.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/ImageIO/Meta/include/metaTypes.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/ImageIO/Meta/include/metaImage.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/ImageIO/Meta/include/metaObject.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/ImageIO/Meta/include/metaIOConfig.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/ImageIO/Meta/include/localMetaConfiguration.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/ImageIO/PlanIm/include/plan_xdr.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/ImageIO/PlanIm/include/plan_config.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/ImageIO/PlanIm/include/strings.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/ImageIO/PlanIm/include/gen.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/ImageIO/PlanIm/include/types.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/ImageIO/PlanIm/include/xdr_ll_planio.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/ImageIO/PlanIm/include/xdr.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/ImageIO/PlanIm/include/plan_sys.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/ImageIO/PlanIm/include/unistd.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/ImageIO/PlanIm/include/getopt.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/ImageIO/PlanIm/include/plan_strings.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/ImageIO/PlanIm/include/rpc_winnt/types.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/ImageIO/PlanIm/include/rpc_winnt/xdr.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/ImageIO/PlanIm/include/libplan.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/ImageIO/PlanIm/include/plan_im.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/ImageIO/PlanIm/include/libplanio.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/ImageIO/PlanIm/include/libbrachy.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/ImageIO/PlanIm/include/extbm.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/ImageIO/PlanIm/include/plan_file_io.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/ImageIO/PlanIm/include/libmisc.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/ImageIO/include/const.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/ImageIO/include/ImageStruct.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/ImageIO/include/AllImageIO.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/ImageIO/include/macros.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/ImageIO/include/ImageIO.h"
    "/work/ltu/Correspondence_Oct7/binPablo2_Oct7/lib/ImageIO/ImageIO_EXPORT.h"
    )
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib64" TYPE STATIC_LIBRARY FILES "/work/ltu/Correspondence_Oct7/binPablo2_Oct7/libImageIO.a")
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

