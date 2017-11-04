#ifndef P3D_DISPLAY_GLOBALS_H
#define P3D_DISPLAY_GLOBALS_H

#include "M3DDisplayGlobals.h"

// Saves the last scale at which the model was drawn
extern float gLastScale;

// Master display list.  All other display lists are
//   nested in this one.
extern unsigned int gMasterDisplayList;
extern bool gMasterDispListIsValid;

// Global display list for the figure surfaces
extern unsigned int gSurfaceDisplayList;

// Global display list for the auxiliary surfaces
extern unsigned int gAuxSurfaceDisplayList;

// Global display list for the medial model
extern unsigned int gModelDisplayList;

// Global display list for the image planes
extern unsigned int gImagePlanesDisplayList;

// Global display list for the boundary curves
extern unsigned int gBoundaryCurveDisplayList;

#ifdef BINARY
// Global display list for the distance vectors
extern unsigned int gDistanceVectorDisplayList;
#endif

// Flag assuring the display lists are generated once
extern bool displayListsInit;

// The height and width of the view window, in model units
extern const double VIEW_HEIGHT;
extern const double VIEW_WIDTH;

// Length of the standard axes
extern const double STD_AXES_LENGTH;

// Size of the glSelectBuffer
extern const int BUFSIZE;

#endif

