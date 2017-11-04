IF(BUILD_TEST)



#----------------------------------------------------------------------------
# USER! : SET THE NAME OF YOUR EXECUTABLE
# Replace 'MyExe' by the name you want to give your executable.
# (a good policy is to give the executable the same name that the directory)

#########################
SET ( EXE_NAME   brainDeflation  )
#########################

#----------------------------------------------------------------------------

#----------------------------------------------------------------------------
# EXECUTABLE SOURCES (TO BE COMPILED)
# EITHER LIST ALL .cxx, *.cpp, *.cc IN CURRENT DIR USING NEXT LINE:

#FILE(GLOB ${EXE_NAME}_SOURCES *.cxx *.cpp *.cc)

# OR MANUALLY LIST YOUR FILES WITH NEXT COMMAND (WITHOUT EXTENSION)
  SET ( ${EXE_NAME}_SOURCES
      brainDeflation.cxx
    )
#----------------------------------------------------------------------------

INCLUDE_DIRECTORIES (

# USER! : Add here the directories holding the extra .h files you need
# e.g.
 ../../lib/mySampleLib

)

#----------------------------------------------------------------------------
# DEPENDENCIES (LIBRARIES TO LINK WITH)
SET ( ${EXE_NAME}_LINK_LIBRARIES
  #    ${crea_LIBRARIES}
  #    ${WXWIDGETS_LIBRARIES}
  #    ${KWWidgets_LIBRARIES}
      vtkRendering
      vtkFiltering
      vtkCommon
      vtkHybrid
      itkvnl
  #    ${ITK_LIBRARIES}
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

#----------------------------------------------------------------------------
# CREATES AND INSTALLS THE EXE
#CREA_ADD_EXECUTABLE( ${EXE_NAME} )
#----------------------------------------------------------------------------


#########################
SET ( EXE_NAME   brainDeflationMultiple  )
#########################

#----------------------------------------------------------------------------

#----------------------------------------------------------------------------
# EXECUTABLE SOURCES (TO BE COMPILED)
# EITHER LIST ALL .cxx, *.cpp, *.cc IN CURRENT DIR USING NEXT LINE:

#FILE(GLOB ${EXE_NAME}_SOURCES *.cxx *.cpp *.cc)

# OR MANUALLY LIST YOUR FILES WITH NEXT COMMAND (WITHOUT EXTENSION)
  SET ( ${EXE_NAME}_SOURCES
      brainDeflationMultiple.cxx
    )
#----------------------------------------------------------------------------

INCLUDE_DIRECTORIES (

# USER! : Add here the directories holding the extra .h files you need
# e.g.
 ../../lib/mySampleLib

)

#----------------------------------------------------------------------------
# DEPENDENCIES (LIBRARIES TO LINK WITH)
SET ( ${EXE_NAME}_LINK_LIBRARIES
  #    ${crea_LIBRARIES}
  #    ${WXWIDGETS_LIBRARIES}
  #    ${KWWidgets_LIBRARIES}
      vtkRendering
      vtkFiltering
      vtkCommon
      vtkHybrid
      itkvnl    
  #    ${ITK_LIBRARIES}
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

#----------------------------------------------------------------------------
# CREATES AND INSTALLS THE EXE
#CREA_ADD_EXECUTABLE( ${EXE_NAME} )





#########################
SET ( EXE_NAME   brainDeflationCortex  )
#########################

#----------------------------------------------------------------------------

#----------------------------------------------------------------------------
# EXECUTABLE SOURCES (TO BE COMPILED)
# EITHER LIST ALL .cxx, *.cpp, *.cc IN CURRENT DIR USING NEXT LINE:

#FILE(GLOB ${EXE_NAME}_SOURCES *.cxx *.cpp *.cc)

# OR MANUALLY LIST YOUR FILES WITH NEXT COMMAND (WITHOUT EXTENSION)
  SET ( ${EXE_NAME}_SOURCES
      brainDeflationCortex.cxx
    )
#----------------------------------------------------------------------------

INCLUDE_DIRECTORIES (

# USER! : Add here the directories holding the extra .h files you need
# e.g.
 ../../lib/mySampleLib

)

#----------------------------------------------------------------------------
# DEPENDENCIES (LIBRARIES TO LINK WITH)
SET ( ${EXE_NAME}_LINK_LIBRARIES
  #    ${crea_LIBRARIES}
  #    ${WXWIDGETS_LIBRARIES}
  #    ${KWWidgets_LIBRARIES}
      vtkRendering
      vtkFiltering      
      vtkCommon 
      vtkHybrid
      itkvnl
  #    ${ITK_LIBRARIES}
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

#----------------------------------------------------------------------------
# CREATES AND INSTALLS THE EXE
CREA_ADD_EXECUTABLE( ${EXE_NAME} )



#########################
SET ( EXE_NAME   brainExtractCorticalSheet  )
#########################

# OR MANUALLY LIST YOUR FILES WITH NEXT COMMAND (WITHOUT EXTENSION)
  SET ( ${EXE_NAME}_SOURCES
      brainExtractCorticalSheet.cxx
    )
#----------------------------------------------------------------------------

INCLUDE_DIRECTORIES (

# USER! : Add here the directories holding the extra .h files you need
# e.g.
 ../../lib/mySampleLib

)

#----------------------------------------------------------------------------
# DEPENDENCIES (LIBRARIES TO LINK WITH)
SET ( ${EXE_NAME}_LINK_LIBRARIES
  #    ${crea_LIBRARIES}
  #    ${WXWIDGETS_LIBRARIES}
  #    ${KWWidgets_LIBRARIES}
      vtkRendering
      vtkFiltering
      vtkCommon
      itkvnl
  #    ${ITK_LIBRARIES}
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

#----------------------------------------------------------------------------
# CREATES AND INSTALLS THE EXE
#CREA_ADD_EXECUTABLE( ${EXE_NAME} )



#########################
SET ( EXE_NAME   brainQuadMesh  )
#########################

# OR MANUALLY LIST YOUR FILES WITH NEXT COMMAND (WITHOUT EXTENSION)
  SET ( ${EXE_NAME}_SOURCES
      brainQuadMesh.cxx
    )
#----------------------------------------------------------------------------

INCLUDE_DIRECTORIES (

# USER! : Add here the directories holding the extra .h files you need
# e.g.
 ../../lib/mySampleLib

)

#----------------------------------------------------------------------------
# DEPENDENCIES (LIBRARIES TO LINK WITH)
SET ( ${EXE_NAME}_LINK_LIBRARIES
  #    ${crea_LIBRARIES}
  #    ${WXWIDGETS_LIBRARIES}
  #    ${KWWidgets_LIBRARIES}
      vtkRendering
      vtkFiltering
      vtkHybrid
      vtkCommon
      itkvnl

      SRep
      SRepInterpolation
      SRepVisualization
  #    ${ITK_LIBRARIES}
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

#----------------------------------------------------------------------------
# CREATES AND INSTALLS THE EXE
CREA_ADD_EXECUTABLE( ${EXE_NAME} )






#########################
SET ( EXE_NAME   testStereographic  )
#########################

# OR MANUALLY LIST YOUR FILES WITH NEXT COMMAND (WITHOUT EXTENSION)
  SET ( ${EXE_NAME}_SOURCES
      testStereographic.cxx
    )
#----------------------------------------------------------------------------

INCLUDE_DIRECTORIES (

# USER! : Add here the directories holding the extra .h files you need
# e.g.
 ../../lib/mySampleLib

)

#----------------------------------------------------------------------------
# DEPENDENCIES (LIBRARIES TO LINK WITH)
SET ( ${EXE_NAME}_LINK_LIBRARIES
  #    ${crea_LIBRARIES}
  #    ${WXWIDGETS_LIBRARIES}
  #    ${KWWidgets_LIBRARIES}
      vtkRendering
      vtkFiltering
      vtkHybrid
      vtkCommon
      itkvnl

      SRep
      SRepInterpolation
      SRepVisualization
  #    ${ITK_LIBRARIES}
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

#----------------------------------------------------------------------------
# CREATES AND INSTALLS THE EXE
CREA_ADD_EXECUTABLE( ${EXE_NAME} )


ENDIF(BUILD_TEST)
