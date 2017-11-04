MESSAGE(STATUS "")
MESSAGE(STATUS "=======================================")
MESSAGE(STATUS "Configuring User html doc")

# Configure & install images 
USER_DOC_INSTALL_IMAGES(template_WebSite)

# Configure html pages for build tree and install 
FILE(GLOB HTML_PAGES "." "*.html")
FOREACH(page ${HTML_PAGES})  
  GET_FILENAME_COMPONENT(filename "${page}" NAME)
  CONFIGURE_FILE(
    ${page}
    ${${PROJECT_NAME}_DOC_BUILD_PATH}/template_WebSite/${filename}
    @ONLY
    )
  # Install
  INSTALL(
    FILES ${${PROJECT_NAME}_DOC_BUILD_PATH}/template_WebSite/${filename}
    DESTINATION ${${PROJECT_NAME}_DOC_INSTALL_PATH}/template_WebSite
    )
ENDFOREACH(page)

	   
MESSAGE(STATUS "=======================================")
