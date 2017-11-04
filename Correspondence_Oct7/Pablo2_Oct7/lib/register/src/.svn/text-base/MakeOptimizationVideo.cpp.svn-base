#include <iostream>
using namespace std ;

#include <fstream>
#include <string> 
#include <sstream>

#include "MakeOptimizationVideo.h"
#include "M3DObjectFile.h"

/*
	Author: Dibyendusekhar Goswami, May 2011

	MakeOptimizationVideo:	A class (namespace) to wrap a function that can write out the M3DObject that is passed to 
							it in a sequence. 

	So, this can be used to record the sequence of models during the Pablo optimization process 
	(I call this sequence a video).
	These models can be later loaded into Pablo (GUI) using [ Tools -> Model/Tile Slide Show ]
	and a video (of .ppm images) can be made by checking the [ Tools -> Make Movie ] checkbox

*/

namespace MakeOptimizationVideo {

	bool makeVideoOn = 0 ;

	char * videoFolder = NULL ;

	int frameCount = 0 ;

	// dibyendu 
	// This function writes out the M3D object that is passed to it in the folder videoFolder 
	// specified by -videoFolder option in Pablo input

	void writeCurrentModel( M3DObject* const currObject ) {	

		if( ! makeVideoOn || ! videoFolder ) {
			cout << "Video option not specified in -videoFolder or -vf in Pablo input" << endl ;
			return ;
		}

		if( currObject == NULL ) {
			cout << "NULL model passed to writeCurrentModel()" << endl ;
			return ;
		}

		frameCount ++ ;

		char outputFileName[100] ;

		// if the last character in videoFolder is a /, then just append "n.m3d" to it, else apppend "

		if( videoFolder[strlen(videoFolder)-1] == '/' || videoFolder[strlen(videoFolder)-1] == '\\' )
			sprintf( outputFileName, "%s%d%s", videoFolder, frameCount, ".m3d" ) ;
		else
			sprintf( outputFileName, "%s%s%d%s", videoFolder, "/", frameCount, ".m3d" ) ;			

		M3DObjectFile mFile ;

#define _DIBYENDU_CHECK_DILATION

#ifdef _DIBYENDU_CHECK_DILATION

		// check if the current object is dilated
		// if it is dilated, then copy it into a new object, erode the new object by the dilation amount, and write the new object (in original size)

		if( currObject->isObjectDilated() ) {

			M3DObject * objToWrite = currObject->assign() ;						// copy into a new object

			objToWrite->erode( currObject->getDilationFactorInModelUnits() ) ;	// erode the new object
			
			mFile.write( outputFileName, *objToWrite, 0, 0, -1, 0, 0, 0 ) ;

			if( objToWrite ) {
				delete objToWrite ;											// delete the new object
				objToWrite = NULL ;
			}
		}	
		else
			mFile.write( outputFileName, *currObject, 0, 0, -1, 0, 0, 0 ) ;

#else

		mFile.write( outputFileName, *currObject, 0, 0, -1, 0, 0, 0 ) ;

	#endif // _DIBYENDU_CHECK_DILATION

	}

}	// end of namespace MakeOptimizationVideo
