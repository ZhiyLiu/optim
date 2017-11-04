#ifndef M3D_DISPLAY_GLOBALS_H
#define M3D_DISPLAY_GLOBALS_H

#ifdef WIN32
    #include <windows.h>
#endif
#ifdef __APPLE__
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#else
#include <GL/gl.h>
#include <GL/glu.h>
#endif
#ifndef GLUquadric
typedef GLUquadricObj GLUquadric;
#endif

// ID signifying to ignore selection (used with GLselect)
extern const int IGNORE_SELECTION_ID;

// The sizes of selected and unselected primitive points
extern const double LARGE_POINT_SIZE;
extern const double SMALL_POINT_SIZE;

#endif

