MESSAGE(STATUS "")
MESSAGE(STATUS "=======================================")
MESSAGE(STATUS "Configuring user documentation ")

SET(ARTICLE_LINKS "<h1> Conference Articles </h1><br>")

#The name of the conference should be the same name as the folder where the *.tex file is located
SET(CONFERENCE_NAME "TemplateIEEE")
#Filename of the *.tex file without extension, the *.bib file must be named the same
SET(ARTICLE_NAME "bare_conf")

USER_BUILD_LATEX_DOC("1" ${CONFERENCE_NAME} ${ARTICLE_NAME})


#SET(ARTICLE_LINKS "<h1> Journal Articles </h1><br>")

#The name of the conference should be the same name as the folder where the *.tex file is located
#SET(CONFERENCE_NAME "TemplateIEEE")
#Filename of the *.tex file without extension, the *.bib file must be named the same
#SET(ARTICLE_NAME "bare_conf")

#USER_BUILD_LATEX_DOC("1" ${CONFERENCE_NAME} ${ARTICLE_NAME})


#SET(ARTICLE_LINKS "<h1> Posters </h1><br>")
#The name of the conference should be the same name as the folder where the *.tex file is located
#SET(CONFERENCE_NAME "POSTERCONF")
#Filename of the *.tex file without extension, the *.bib file must be named the same
#SET(ARTICLE_NAME "MYPOSTER")


CONFIGURE_FILE(
    ${CMAKE_CURRENT_SOURCE_DIR}/UserDocumentationLinks.txt.in
    ${CMAKE_CURRENT_BINARY_DIR}/../UserDoxygen/UserDocumentationLinks.txt    
    ) 



MESSAGE(STATUS "=======================================")
