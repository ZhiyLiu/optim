# Install script for directory: /work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/gui

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
  FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib64/gui" TYPE FILE FILES "/work/ltu/Correspondence_Oct7/binPablo2_Oct7/lib/gui/InstallLibraryForCMake_tmp/Usegui.cmake")
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib64/gui" TYPE FILE FILES "/work/ltu/Correspondence_Oct7/binPablo2_Oct7/lib/gui/InstallLibraryForCMake_tmp/guiConfig.cmake")
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib64/gui" TYPE FILE FILES "/work/ltu/Correspondence_Oct7/binPablo2_Oct7/lib/gui/InstallLibraryForCMake_tmp/guiBuildSettings.cmake")
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  list(APPEND CPACK_ABSOLUTE_DESTINATION_FILES
   "/usr/local/share/cmake/Modules/Findgui.cmake")
FILE(INSTALL DESTINATION "/usr/local/share/cmake/Modules" TYPE FILE FILES "/work/ltu/Correspondence_Oct7/binPablo2_Oct7/lib/gui/InstallLibraryForCMake_tmp/Findgui.cmake")
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/gui" TYPE FILE FILES
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/gui/src/button_labels.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/gui/src/P3DUserInterface.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/gui/src/make_windows.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/gui/src/LogManagerChart.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/gui/src/OptVisualizerCallback.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/gui/src/Plot3DWindowWrapper.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/gui/src/menu.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/gui/src/P3DUserInterfaceCallback.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/gui/src/OptVisualizerUIDefs.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/gui/src/OptVisualizerUI.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/gui/src/Fl_Aspect_Ratio_Group.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/gui/src/movable_Fl_Window.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/gui/src/P3DDisplayGlobals.h"
    "/work/ltu/Correspondence_Oct7/Pablo2_Oct7/lib/gui/src/P3DView.h"
    "/work/ltu/Correspondence_Oct7/binPablo2_Oct7/lib/gui/gui_EXPORT.h"
    )
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib64" TYPE STATIC_LIBRARY FILES "/work/ltu/Correspondence_Oct7/binPablo2_Oct7/libgui.a")
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

