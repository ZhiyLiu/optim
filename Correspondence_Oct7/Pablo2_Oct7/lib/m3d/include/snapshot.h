#ifndef SNAPSHOT_H
#define SNAPSHOT_H


#include <float.h>
#include "Image3D.h"

#include <limits.h>

/*  Image snapshot functions.

	Note: These should not be used for routine image output.  They are
	designed only to be used for debugging operations.

	The first 3 functions produce a Raw3 image, given the file name, an
	array of pixels, and the dimensions of the array.  The intensity
	range may be provided, if known, but is not needed.

	The last function is similar, but takes an Image3D object.
*/


void snapshot(const char * filename, float * pixels, int x, int y, int z,
	float min = FLT_MAX, float max = FLT_MIN, float xSpacing = 0.0f,
	float ySpacing = 0.0f, float zSpacing = 0.0f);

void snapshot(const char * filename, GreyValue * pixels, int x, int y, int z,
	GreyValue min = MAX_GREY_VALUE, GreyValue max = MIN_GREY_VALUE,
	float xSpacing = 0.0f, float ySpacing = 0.0f, float zSpacing = 0.0f);

void snapshot(const char * filename, unsigned char * pixels, int x, int y, int z,
	unsigned char min = CHAR_MAX, unsigned char max = CHAR_MIN,
	float xSpacing = 0.0f, float ySpacing = 0.0f, float zSpacing = 0.0f);


void snapshot(const char * filename, Image3D & image,
	GreyValue min = MAX_GREY_VALUE, GreyValue max = MIN_GREY_VALUE);


#endif

