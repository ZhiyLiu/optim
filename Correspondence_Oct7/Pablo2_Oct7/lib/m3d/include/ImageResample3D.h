#ifndef IMAGE_RESAMPLE_3D
#define IMAGE_RESAMPLE_3D

#include "Image3D.h"


class ImageResample3D
{
public:

    ImageResample3D() {}

    // Uses the smallest dimension of a voxel as the isotropic
    // size for new voxels.  Attraction and repulsion masks are
    // only used with stacked images.
    void isotropicSample(Image3D & image, double userSpacing = 0.0,
		GreyValue attractionMask = 1, GreyValue repulsionMask = 0);

#ifdef BINARY

	// Pad image out so that all extents are maxExtent cm long.
	// If maxExtent is omitted, use the largest extent of any
	// dimension as the extent for the new image.  Return
	// 1 if OK, 0 if any dimension's extent > maxExtent
	// or cannot alloc memory.
	int cubeSample(Image3D & image, float maxExtent = 0);

#endif

	bool superSample(Image3D & image, int outGridDims);

};



#endif

