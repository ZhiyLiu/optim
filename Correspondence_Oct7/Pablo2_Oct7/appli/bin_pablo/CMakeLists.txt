#----------------------------------------------------------------------------
# USER! : SET THE NAME OF YOUR EXECUTABLE
# Replace 'MyExe' by the name you want to give your executable.
# (a good policy is to give the executable the same name that the directory)

#########################
SET ( EXE_NAME   bin_pablo  )
#########################

#----------------------------------------------------------------------------

#----------------------------------------------------------------------------
# EXECUTABLE SOURCES (TO BE COMPILED)
# EITHER LIST ALL .cxx, *.cpp, *.cc IN CURRENT DIR USING NEXT LINE:

#FILE(GLOB ${EXE_NAME}_SOURCES *.cxx *.cpp *.cc)

# OR MANUALLY LIST YOUR FILES WITH NEXT COMMAND (WITHOUT EXTENSION)
  SET ( ${EXE_NAME}_SOURCES
      main.cpp
    )
#----------------------------------------------------------------------------

INCLUDE_DIRECTORIES (

# USER! : Add here the directories holding th extra .h files you need
# e.g.
# ../../lib/<my_library_I_just_created>

)

#----------------------------------------------------------------------------
# DEPENDENCIES (LIBRARIES TO LINK WITH)
SET ( ${EXE_NAME}_LINK_LIBRARIES
  #    ${crea_LIBRARIES}
  #    ${WXWIDGETS_LIBRARIES}
  #    ${KWWidgets_LIBRARIES}
  #    ${VTK_LIBRARIES}
  #    ${ITK_LIBRARIES}
  #    ${GDCM_LIBRARIES}
  #    ${BOOST_LIBRARIES}
  ${FLTK_LIBRARIES}
  flvw-1.0
  ${CORELIBS}
  ${OPENGL_LIBRARIES}
  gui
  ImageIO
  m3d
  #match
  paul_code
  planes
  register
  seurat
  version
  zlib
  
  

  
  
  # USER! : Add here those agmonst the various (?) PROJECT LIBRARIES
  # you need for the current executable
  # (If you created only one Library, don't forget it !...) 
  
  )
#----------------------------------------------------------------------------

#----------------------------------------------------------------------------
# USER! : UNCOMMENT NEXT LINE IF YOU WANT A CONSOLE ON WINDOWS
# NB : YOUR MAIN MUST BE ADAPTED ALSO
#      SEE THE MACRO CREA_WXMAIN_WITH_CONSOLE IN creaWx.h
#SET(${EXE_NAME}_CONSOLE TRUE)
#----------------------------------------------------------------------------

#----------------------------------------------------------------------------
# CREATES AND INSTALLS THE EXE
CREA_ADD_EXECUTABLE( ${EXE_NAME} )
#add_executable(optimizeSrep OptimizeSrep.cpp)
#target_link_libraries(optimizeSrep ${${EXE_NAME}_LINK_LIBRARIES})

IF(WIN32)
	SET(IGNORE_LIBRARY_FLAGS
	   "/NODEFAULTLIB:libcd.lib"
	)
	SET_TARGET_PROPERTIES(${EXE_NAME} PROPERTIES
		   LINK_FLAGS_DEBUG
		   ${IGNORE_LIBRARY_FLAGS})
ENDIF(WIN32)
#----------------------------------------------------------------------------


#----------------------------------------------------------------------------
# USER! : SET THE NAME OF YOUR EXECUTABLE
# Replace 'MyExe' by the name you want to give your executable.
# (a good policy is to give the executable the same name that the directory)

OPTION(BUILD_ANTIALIAS OFF)

IF(BUILD_ANTIALIAS)


CREA_CPDIR(${PROJECT_SOURCE_DIR}/appli/bin_pablo/antiAlias ${PROJECT_BINARY_DIR})

#########################
SET ( EXE_NAME   Image2SignedDistanceMap  )
#########################

#----------------------------------------------------------------------------

#----------------------------------------------------------------------------
# EXECUTABLE SOURCES (TO BE COMPILED)
# EITHER LIST ALL .cxx, *.cpp, *.cc IN CURRENT DIR USING NEXT LINE:

#FILE(GLOB ${EXE_NAME}_SOURCES *.cxx *.cpp *.cc)

# OR MANUALLY LIST YOUR FILES WITH NEXT COMMAND (WITHOUT EXTENSION)
  SET ( ${EXE_NAME}_SOURCES
     Image2Distance.cxx
    )
#----------------------------------------------------------------------------

INCLUDE_DIRECTORIES (

# USER! : Add here the directories holding th extra .h files you need
# e.g.
# ../../lib/<my_library_I_just_created>

)

#----------------------------------------------------------------------------
# DEPENDENCIES (LIBRARIES TO LINK WITH)
SET ( ${EXE_NAME}_LINK_LIBRARIES
  #    ${crea_LIBRARIES}
  #    ${WXWIDGETS_LIBRARIES}
  #    ${KWWidgets_LIBRARIES}
  #    ${VTK_LIBRARIES}
      ${ITK_LIBRARIES}
  #    ${GDCM_LIBRARIES}
  #    ${BOOST_LIBRARIES}
  
  
  # USER! : Add here those agmonst the various (?) PROJECT LIBRARIES
  # you need for the current executable
  # (If you created only one Library, don't forget it !...) 
  
  )
#----------------------------------------------------------------------------

#----------------------------------------------------------------------------
# USER! : UNCOMMENT NEXT LINE IF YOU WANT A CONSOLE ON WINDOWS
# NB : YOUR MAIN MUST BE ADAPTED ALSO
#      SEE THE MACRO CREA_WXMAIN_WITH_CONSOLE IN creaWx.h
#SET(${EXE_NAME}_CONSOLE TRUE)
#----------------------------------------------------------------------------

#-------------BUILD_IMAGE2DISTANCEMAP---------------------------------------------------------------
# CREATES AND INSTALLS THE EXE
CREA_ADD_EXECUTABLE( ${EXE_NAME} )
ENDIF(BUILD_ANTIALIAS)

### ----------------Build NEWUOA optimizer



