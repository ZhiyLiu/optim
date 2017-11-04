#ifndef _MAKE_OPTIMIZATION_VIDEO
#define _MAKE_OPTIMIZATION_VIDEO

/*
	Author: Dibyendusekhar Goswami, May 2011

	MakeOptimizationVideo:	A class (namespace) to wrap a function that can write out the M3DObject that is passed to 
							it in a sequence. 

	So, this can be used to record the sequence of models during the Pablo optimization process 
	(I call this sequence a video).
	These models can be later loaded into Pablo (GUI) using [ Tools -> Model/Tile Slide Show ]
	and a video (of .ppm images) can be made by checking the [ Tools -> Make Movie ] checkbox

*/


#include "M3DObject.h"

namespace MakeOptimizationVideo {

// public :

	extern bool makeVideoOn ;

	extern char * videoFolder ;

	extern int frameCount ;

	void turnVideoOn( char * _videoFolder ) ;

	void turnVideoOff() ;

	void setVideoFolderName( char * _videoFolder ) ;

	// dibyendu 
	// This function writes out the M3D object that is passed to it in the folder videoFolder 
	// specified by -videoFolder option in Pablo input
	
	void writeCurrentModel( M3DObject* const currObject ) ;
}

#endif  //_MAKE_OPTIMIZATION_VIDEO 

