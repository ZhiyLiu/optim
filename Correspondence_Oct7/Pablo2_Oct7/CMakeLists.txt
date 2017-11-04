# ---------------------------------
#
# Unless you are an experimented cmake user,
# have a look at :
# http://www.creatis.insa-lyon.fr/site/fr/CreatoosGettingStarted-v2.0.X
# before starting.
#
# Generated with creaNewProject, version of : 18/03/2010
# ---------------------------------

CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

#MARK_AS_ADVANCED( FORCE CMAKE_BACKWARDS_COMPATIBILITY )

# for CMake 2.6 corrected behaviour (see "cmake --help-policy CMP0003")
IF(COMMAND cmake_policy AND ${CMAKE_MAJOR_VERSION} EQUAL 2 AND ${CMAKE_MINOR_VERSION} GREATER 4)
  CMAKE_POLICY(SET CMP0003 NEW)
  CMAKE_POLICY(SET CMP0005 NEW)
  CMAKE_POLICY(SET CMP0011 NEW)
ENDIF(COMMAND cmake_policy AND ${CMAKE_MAJOR_VERSION} EQUAL 2 AND ${CMAKE_MINOR_VERSION} GREATER 4)


###################################
PROJECT( Pablo2 )
###################################

#==================================
# The project version
SET(PROJECT_MAJOR_VERSION 2)
SET(PROJECT_MINOR_VERSION 1)
SET(PROJECT_BUILD_VERSION 0)

SET(PROJECT_VERSION "${PROJECT_MAJOR_VERSION}.${PROJECT_MINOR_VERSION}.${PROJECT_BUILD_VERSION}")
SET(PROJECT_VERSION_DATE "11/02/2010")
#ADD_DEFINITIONS( -DPROJECT_VERSION=${PROJECT_VERSION} )
#ADD_DEFINITIONS( -DPROJECT_VERSION_DATE=${PROJECT_VERSION_DATE} )


# USER! : Add here as many 'FIND_PACKAGE' blocks as *your* aplication requires
# -----

# WARNING : we are talking here about 'bbtk packages', not about libraries!
# if your project just requires the GDCM library, and not bbtk/packages/gdcm,
# do not add FIND_PACKAGE(gdcm REQUIRED).
# Just set USE_GDCM to ON, after " Libraries/tools used

INCLUDE(${PROJECT_SOURCE_DIR}/cmake/UserRootCMakeLists.cmake)

#INCLUDE(${PROJECT_SOURCE_DIR}/cmake/FindMATLAB.cmake)

#set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${PROJECT_SOURCE_DIR}/cmake")

#IF(MATLAB_FOUND)
#    find_package(MATLAB REQUIRED)
#    include_directories(${MATLAB_INCLUDE_DIR})
#ENDIF(MATLAB_FOUND)

#INCLUDE_DIRECTORIES("/NIRAL/tools/Matlab2011a/extern/include")


#==================================
# Do not allow to build inside the source tree
CREA_PREVENT_IN_SOURCE_BUILD()
#==================================

#==================================
# Libraries/tools used

# USER! : Note : Set USE_CREA to ON
# -----

#  if you need to LINK against crea
FIND_PACKAGE(OpenGL REQUIRED)


IF(NOT WIN32)
	find_package(BLAS REQUIRED )
	find_package(LAPACK REQUIRED)
	FIND_PACKAGE(FLTK REQUIRED)
    INCLUDE_DIRECTORIES(${FLTK_INCLUDE_DIR}})
ELSE(NOT WIN32)

	SET(FLTK_LIBRARIES
		fltk
		fltkd
		fltkforms
		fltkformsd
		fltkgl
		fltkgld
		fltkimages
		fltkimagesd
		fltkjpeg
		fltkjpegd
		fltkpng
		fltkpngd
		fltkz
		fltkzd
		)
	INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR}/include_ext)
	LINK_DIRECTORIES(${PROJECT_SOURCE_DIR}/include_ext/lib_ext/FLTK)

	LINK_DIRECTORIES(${PROJECT_SOURCE_DIR}/include_ext/lib_ext/LAPACK)

    INCLUDE_DIRECTORIES(${FLTK_INCLUDE_DIR}})

	SET(BLAS_LIBRARIES BLASd_nowrap)
	SET(LAPACK_LIBRARIES
		#clapack_debug
		#clapackMD
		clapack
		clapackd_nowrap
		)

ENDIF(NOT WIN32)
#==================================

#==================================

#==================================
# Where to put executables and libs
SET(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR})
SET(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR})
MARK_AS_ADVANCED(
  CMAKE_BACKWARDS_COMPATIBILITY
  EXECUTABLE_OUTPUT_PATH
  LIBRARY_OUTPUT_PATH
  )
#==================================

IF(WIN32)
  INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR}/win32)
  LINK_DIRECTORIES(${PROJECT_SOURCE_DIR}/win32)
ENDIF(WIN32)

#==================================

#==================================

  SET(LIBNAME $(Pablo2))


add_definitions(-DBINARY -DLM_METHOD_OBJ -DNO_MATCH_PENALTY)

IF(WIN32)
	add_definitions(-DPLAN_WINNT)
ENDIF(WIN32)


  INCLUDE_DIRECTORIES(
    ${PROJECT_BINARY_DIR}

    ${PROJECT_SOURCE_DIR}/lib/include
    ${PROJECT_SOURCE_DIR}/lib/gui/src
    ${PROJECT_SOURCE_DIR}/lib/ImageIO/include
    ${PROJECT_SOURCE_DIR}/lib/ImageIO/PlanIm/include
    ${PROJECT_SOURCE_DIR}/lib/ImageIO/utilities/include
    ${PROJECT_SOURCE_DIR}/lib/ImageIO/Analyze/include/
    ${PROJECT_SOURCE_DIR}/lib/ImageIO/Gipl/include
    ${PROJECT_SOURCE_DIR}/lib/ImageIO/Meta/include/
    ${PROJECT_SOURCE_DIR}/lib/m3d/include
    ${PROJECT_SOURCE_DIR}/lib/match/include
    ${PROJECT_SOURCE_DIR}/lib/paul_code/include
    ${PROJECT_SOURCE_DIR}/lib/planes/include
    ${PROJECT_SOURCE_DIR}/lib/register/include
    ${PROJECT_SOURCE_DIR}/lib/seurat/include
    ${PROJECT_SOURCE_DIR}/lib/seurat/src
    ${PROJECT_SOURCE_DIR}/lib/version/include
    ${PROJECT_SOURCE_DIR}/lib/zlib/include
    ${PROJECT_SOURCE_DIR}/lib/flvw-1.0/FL
    ${PROJECT_SOURCE_DIR}/lib/flvw-1.0/

    ${PROJECT_SOURCE_DIR}/lib/vtksrep/SRepLib
    ${PROJECT_SOURCE_DIR}/lib/vtksrep/SRepIOLib
    ${PROJECT_SOURCE_DIR}/lib/vtksrep/SRepInterpolationLib
    ${PROJECT_SOURCE_DIR}/lib/vtksrep/SRepVisualizationLib
    ${PROJECT_SOURCE_DIR}/lib/vtksrep/calquantum

    ${PROJECT_SOURCE_DIR}/lib/vtksrep/Correspondence
    ${PROJECT_SOURCE_DIR}/lib/vtksrep/Visualization
    ${PROJECT_SOURCE_DIR}/lib/vtksrep/vtkpowercrust



  )

 # USER! : Add here the directories holding the .h you need

 # INCLUDE_DIRECTORIES(appli/_YOUR_OWN_APPLI_SOURCE_DIRECTORY_)
 # INCLUDE_DIRECTORIES(lib/_YOUR_OWN_LIBRARY_SOURCE_DIRECTORY_)

#==================================



FIND_PACKAGE(VTK)
IF(VTK_FOUND)
    INCLUDE(${VTK_USE_FILE})
ENDIF(VTK_FOUND)

FIND_PACKAGE(ITK)
IF(ITK_FOUND)

    IF(ITK_BUILD_SETTINGS_FILE AND NOT SKIP_ITK_BUILD_SETTINGS_FILE)
      INCLUDE(${CMAKE_ROOT}/Modules/CMakeImportBuildSettings.cmake)
      CMAKE_IMPORT_BUILD_SETTINGS(${ITK_BUILD_SETTINGS_FILE})
    ENDIF(ITK_BUILD_SETTINGS_FILE AND NOT SKIP_ITK_BUILD_SETTINGS_FILE)

    SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${ITK_REQUIRED_C_FLAGS}")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${ITK_REQUIRED_CXX_FLAGS}")
    SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${ITK_REQUIRED_LINK_FLAGS}")
    SET(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} ${ITK_REQUIRED_LINK_FLAGS}")
    SET(CMAKE_MODULE_LINKER_FLAGS "${CMAKE_MODULE_LINKER_FLAGS} ${ITK_REQUIRED_LINK_FLAGS}")

    INCLUDE_DIRECTORIES(${ITK_INCLUDE_DIRS})
    LINK_DIRECTORIES(${ITK_LIBRARY_DIRS})

ENDIF(ITK_FOUND)



#----------------------------------------------------------------------------


#==================================
# Where to generate documentation

SET(BBTK_DOC_BUILD_PATH
     ${PROJECT_BINARY_DIR}/doc
   )
#==================================
# USER! : Hands off!
# ----

OPTION(BUILD_TEST OFF)

# ADD SUB DIRECTORIES
# ----------------
ADD_SUBDIRECTORY(lib)
ADD_SUBDIRECTORY(appli)
ADD_SUBDIRECTORY(doc)
ADD_SUBDIRECTORY(data)
ADD_SUBDIRECTORY(install)
#==================================



