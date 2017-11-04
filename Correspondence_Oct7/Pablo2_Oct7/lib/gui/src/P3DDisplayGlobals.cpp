// -- P3DDisplayGlobals.cpp --
// This file contains declarations of globals used in the display classes

#include "P3DDisplayGlobals.h"


float gLastScale;

unsigned int gMasterDisplayList = 0;
bool gMasterDispListIsValid = false;

unsigned int gSurfaceDisplayList = 0;

unsigned int gAuxSurfaceDisplayList = 0;

unsigned int gModelDisplayList = 0;

unsigned int gImagePlanesDisplayList = 0;

unsigned int gBoundaryCurveDisplayList = 0;

#ifdef BINARY
unsigned int gDistanceVectorDisplayList = 0;
#endif

bool displayListsInit = false;

const double VIEW_HEIGHT = 2.0;
const double VIEW_WIDTH = 2.0;

const double STD_AXES_LENGTH = 1.0;

// Affects the number of atoms that can be simultaneously selected by a mouse drag
const int BUFSIZE = 1024;


