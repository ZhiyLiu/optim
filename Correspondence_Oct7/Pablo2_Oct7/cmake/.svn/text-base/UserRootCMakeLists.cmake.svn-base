
SET(CREA_CMAKE_DIR ${PROJECT_SOURCE_DIR}/cmake)

IF (WIN32)
    SET(CMAKE_CREA_LIB_PATH bin)
ELSE (WIN32)
    if( NOT APPLE )
      # check 64 bit
      if( ${CMAKE_SIZEOF_VOID_P} EQUAL 4 )
	 set( HAVE_64_BIT 0 )
	 SET(CMAKE_CREA_LIB_PATH lib)
      else( ${CMAKE_SIZEOF_VOID_P}EQUAL 4 )
	 set( HAVE_64_BIT 1 )
	 SET(CMAKE_CREA_LIB_PATH lib64)
      endif( ${CMAKE_SIZEOF_VOID_P} EQUAL 4 )
    else ( NOT APPLE )
	 SET(CMAKE_CREA_LIB_PATH lib)

      if( ${CMAKE_SIZEOF_VOID_P} EQUAL 4 )
	 message("EED crea definitions  --------APPLE------------ 64 0")
      else( ${CMAKE_SIZEOF_VOID_P}EQUAL 4 )
	 message("EED crea definitions  --------APPLE------------ 64 1")
      endif( ${CMAKE_SIZEOF_VOID_P} EQUAL 4 )


    endif( NOT APPLE )
ENDIF(WIN32)

INCLUDE(${CREA_CMAKE_DIR}/UserMacros.cmake)
INCLUDE(${CREA_CMAKE_DIR}/UserConfig.cmake)
INCLUDE(${CREA_CMAKE_DIR}/UserSetDeducedPaths.cmake)
INCLUDE(${CREA_CMAKE_DIR}/UserBuildAllOption.cmake)
INCLUDE(${CREA_CMAKE_DIR}/UserDefineOptions.cmake)
