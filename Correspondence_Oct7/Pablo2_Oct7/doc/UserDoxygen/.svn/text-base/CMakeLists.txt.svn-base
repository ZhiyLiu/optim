MESSAGE(STATUS "")
MESSAGE(STATUS "=======================================")
MESSAGE(STATUS "Configuring user doxygen doc")


#SET THE .h FILES INHERITED FROM OTHER PROJECTS, THIS WILL GENERATE COMPLETE DIAGRAMS OF HERITAGE IN DOXYGEN
SET(EXTERNALINHERITED 
   )


GET_PROPERTY(INC_DIRS DIRECTORY PROPERTY INCLUDE_DIRECTORIES)

#message(status "*********** EXTERNALINHERITED ${EXTERNALINHERITED}")
foreach(fileinherited ${EXTERNALINHERITED})
	FIND_FILE(F0_${fileinherited} ${fileinherited}
		  HINTS ${INC_DIRS})

	SET(ALL_INHERITED "${F0_${fileinherited}} ${ALL_INHERITED}")	
endforeach(fileinherited)

#message(status "*********** ALL_INHERITED ${ALL_INHERITED}")



SET(INPUT 
  ${CMAKE_CURRENT_BINARY_DIR}/UserDocumentationLinks.txt
  ${PROJECT_SOURCE_DIR}/lib
  ${PROJECT_SOURCE_DIR}/appli
  ${ALL_INHERITED}
  )

INCLUDE(../../cmake/UserBuildDoxygenDoc.cmake)

USER_BUILD_DOXYGEN_DOC(
  "${PROJECT_NAME}" 
  "${INPUT}" 
  "${PROJECT_NAME}" 
  "${DOXYGEN_DOC_PREDEFINED}"
  )

MESSAGE(STATUS "=======================================")
