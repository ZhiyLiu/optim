#include <memory.h>
#include <math.h>
#include "ImagePlanes.h"

//#define DEBUG


#ifdef DEBUG
#include <iostream>
using namespace std;
#endif


void ImagePlanes::init(Image3D * im, int size, bool smoothness)
{
    image = im;
	yFlip = false;
	zFlip = false;

	if (size > 0 && size != imagePlaneSize) {
	    imagePlaneSize = size;
		allocatePlaneMemory();

		if (image != NULL) {
			// Start with the central slice of each axis
			xCutPlaneWorldPos = image->imageXToWorld(image->getXDim()/2);
			yCutPlaneWorldPos = image->imageYToWorld(image->getYDim()/2);
			zCutPlaneWorldPos = image->imageZToWorld(image->getZDim()/2);
		}
	}
	else
	    imagePlaneSize = size;

	smooth = smoothness;

	if (im == NULL) {
		imagePlaneSize = 0;
		return;
	}

#ifndef UNFLIPPED
	if (image->getYSpacing() < 0.0)
		yFlip = true;
	if (image->getZSpacing() < 0.0)
		zFlip = true;
#endif
	
	if (im->getIsImageStacked() == true)
		lut_max = MAX_STACKED_GREY_VALUE;
	else {
	    GreyValue min;
		image->range(min, lut_max);
		if (lut_max < min)
		    lut_max = MAX_GREY_VALUE;
	}
#ifdef DEBUG
    cout << "ImagePlanes::lut_max = " << lut_max << endl;
#endif

	if (colorLUT == NULL) {
		unsigned int index;
   		colorLUT = new unsigned char[MAX_GREY_VALUE - MIN_GREY_VALUE + 1];
		for (index = MIN_GREY_VALUE; index <= (unsigned int) MAX_GREY_VALUE; index++)
			colorLUT[index] = (unsigned char)(((float) index / (float) lut_max) * 255.0f);
	}
}

ImagePlanes::ImagePlanes()
{
	colorLUT = NULL;
    xCutPlane = NULL;
    yCutPlane = NULL;
    zCutPlane = NULL;
    arbCutPlane = NULL;
	imagePlaneSize = 0;

	init(NULL, 0, true);
}

ImagePlanes::ImagePlanes(Image3D * im, int size)
{
	colorLUT = NULL;
    xCutPlane = NULL;
    yCutPlane = NULL;
    zCutPlane = NULL;
    arbCutPlane = NULL;
	imagePlaneSize = 0;

	init(im, size, true);
}

// If external is true, any previous image owned by this class will not be
// deleted.  This capability is provided by use by programs other than Pablo.
void ImagePlanes::setImagePtr(Image3D * im, int size, bool external)
{
#ifdef DEBUG
    cout << "ImagePlanes::setImagePtr(im = 0x" << hex << im << dec << ")\n";
#endif
	if (im == NULL || size != imagePlaneSize)
		freePlaneMemory();

	if (image != NULL && image != im && ! external) 
	{
		delete image;
		image = NULL ;
	}

	init(im, size, smooth);
}

ImagePlanes::~ImagePlanes()
{
	freePlaneMemory();

	if( colorLUT != NULL )
		delete [] colorLUT;
	colorLUT = NULL ;

	// dibyendu : there is no need to delete the image because it was never assigned a new memory block
	// the image was only copied from another place

	//if( image != NULL )
	//	delete image;
	//image = NULL ;
}

void ImagePlanes::setSmoothing(bool smoothImages) {
	smooth = smoothImages;
}

unsigned char * ImagePlanes::getXCutPlane(bool final)
{
#ifdef DEBUG
//	std::cout << "X(" << final << ", " << smooth << ")\n";
#endif
	if (final && smooth)
		fillXCutPlane();
	else
		fastFillXCutPlane();

    return xCutPlane;
}

unsigned char * ImagePlanes::getYCutPlane(bool final)
{
#ifdef DEBUG
//	std::cout << "Y(" << final << ", " << smooth << ")\n";
#endif
	if (final && smooth)
		fillYCutPlane();
	else
		fastFillYCutPlane();

    return yCutPlane;
}

unsigned char * ImagePlanes::getZCutPlane(bool final)
{
#ifdef DEBUG
//	std::cout << "Z(" << final << ", " << smooth << ")\n";
#endif
	if (final && smooth)
		fillZCutPlane();
	else
		fastFillZCutPlane();

    return zCutPlane;
}
double ImagePlanes::getXCutPlaneModelPos() const
{
	if (image == NULL)
		return -1.0;

	return image->worldXToModel(xCutPlaneWorldPos);
}

double ImagePlanes::getYCutPlaneModelPos() const
{
	if (image == NULL)
		return -1.0;

	return image->worldYToModel(yCutPlaneWorldPos);
}

double ImagePlanes::getZCutPlaneModelPos() const
{
	if (image == NULL)
		return -1.0;

	return image->worldZToModel(zCutPlaneWorldPos);
}

/*  This function computes the model-space positions of the three slices.
	It combines the three functions above.
*/
bool ImagePlanes::getCutPlaneModelPositions(double & x, double & y,
    double & z) const
{
	if (image == NULL)
		return false;

	x = image->worldXToModel(xCutPlaneWorldPos);
	y = image->worldYToModel(yCutPlaneWorldPos);
	z = image->worldZToModel(zCutPlaneWorldPos);
	return true;
}

int ImagePlanes::getXCutPlaneIndex() const
{
    int xIndex;

	if (image == NULL)
		return -1;

	xIndex = image->worldXToImage(xCutPlaneWorldPos);

    // Set image to black, if we are out of bounds
    if (xIndex < 0 || xIndex >= image->getXDim())
		return -1;

	return xIndex;
}

int ImagePlanes::getYCutPlaneIndex() const
{
    int yIndex;

	if (image == NULL)
		return -1;

	yIndex = image->worldYToImage(yCutPlaneWorldPos);

    // Set image to black, if we are out of bounds
    if (yIndex < 0 || yIndex >= image->getYDim())
		return -1;

	return yIndex;
}

int ImagePlanes::getZCutPlaneIndex() const
{
    int zIndex;

	if (image == NULL)
		return -1;

	zIndex = image->worldZToImage(zCutPlaneWorldPos);

    // Set image to black, if we are out of bounds
    if (zIndex < 0 || zIndex >= image->getZDim())
		return -1;

	return zIndex;
}

unsigned char * ImagePlanes::getArbCutPlane(bool final)
{
	if (final && smooth)
        fillArbCutPlane();
	else
		fastFillArbCutPlane();

    return arbCutPlane;
}

void ImagePlanes::setXCutPlaneSlice(int sliceIndex)
{
	if (image == NULL)
		return;

	xCutPlaneWorldPos = image->imageXToWorld(sliceIndex);
}

void ImagePlanes::setYCutPlaneSlice(int sliceIndex)
{
	if (image == NULL)
		return;

	yCutPlaneWorldPos = image->imageYToWorld(sliceIndex);
}

void ImagePlanes::setZCutPlaneSlice(int sliceIndex)
{
	if (image == NULL)
		return;

	zCutPlaneWorldPos = image->imageZToWorld(sliceIndex);
}

void ImagePlanes::setArbCutPlaneWorldPos(const Vector3D & pos)
{
    if (arbCutPlaneWorldPos != pos)
        arbCutPlaneWorldPos = pos;
}

void ImagePlanes::setArbCutPlaneXDir(const Vector3D & xdir)
{
    if (arbCutPlaneXDir != xdir)
        arbCutPlaneXDir = xdir;
}

void ImagePlanes::setArbCutPlaneYDir(const Vector3D & ydir)
{
    if (arbCutPlaneYDir != ydir)
        arbCutPlaneYDir = ydir;
}

void ImagePlanes::setArbCutPlaneExtent(double extent)
{
    if (arbCutPlaneExtent != extent)
        arbCutPlaneExtent = extent;
}

// Min and max should be in the range [0, 1] (relative intensities)
void ImagePlanes::setIntensityWindow(double min, double max)
{
    int index;

#ifdef DEBUG
    cout << "Called ImagePlanes::setIntensityWindow() with min = " << min
        << " and max = " << max << endl;
#endif
	for (index = MIN_GREY_VALUE; index <= MAX_GREY_VALUE; index++) {
        double level = (double) index / (double) lut_max;
        if (level <= min)
            colorLUT[index] = 0;
        else if (level >= max)
            colorLUT[index] = 255;
        else
            colorLUT[index] = (char) (255.0 * (level - min) / (max - min));
    }
}

// Saggital plane, interpolated
void ImagePlanes::fillXCutPlane()
{
    int xIndex;
    double yIndex, zIndex;
    int yDim, zDim;
    int i, j;
    GreyValue * voxelPtr;
    GreyValue voxel;
    unsigned char * xCutPlanePtr;

    if (image == NULL) {
        freePlaneMemory();
        return;
    }

    if (xCutPlane == NULL)
        return;

    voxelPtr = image->getVoxels();
    if (voxelPtr == NULL) {
        delete [] xCutPlane;
        xCutPlane = NULL;
        return;
    }

	xIndex = image->worldXToImage(xCutPlaneWorldPos);
#ifdef DEBUG
	cout << "ImagePlanes::fillXCutPlane(): Filling X slice "
		<< xIndex << ", " << xCutPlaneWorldPos << " cm" << endl;
#endif

    // Set image to black, if we are out of bounds
    if (xIndex < 0 || xIndex >= image->getXDim()) {
        memset((void*) xCutPlane, 0, imagePlaneSize*imagePlaneSize);
        return;
    }

    yDim = image->getYDim() - 1;
    zDim = image->getZDim() - 1;

    xCutPlanePtr = xCutPlane;
    for (i = 0; i < imagePlaneSize; i++) {
		zIndex = ((double) (i*zDim))/(imagePlaneSize - 1);
        for (j = 0; j < imagePlaneSize; j++) {
			yIndex = ((double) (j*yDim))/(imagePlaneSize - 1);
			// For AE2, masking is done in getInterpolatedVoxelValue()
			voxel = (GreyValue) (0.5 + image->getInterpolatedVoxelValue(
                (double) xIndex, yIndex, zIndex));
            (*xCutPlanePtr) = colorLUT[voxel];
            xCutPlanePtr++;
        }
    }
}

// Note that the fast-fill functions address the voxel array in Image3D
// directly.  Therefore, the pointer arithmetic for them must account
// for the directions of the axes.  The Y and Z axes change direction
// when their corresponding spacing values are negative.  So flipping
// the voxel access in the Y and Z directions must be supported.

// Saggital plane
void ImagePlanes::fastFillXCutPlane()
{
    int xIndex;
    int xDim, yDim, zDim;
    int yIncrement, zIncrement;
    int zLeap;
    int yModulo, zModulo;
    int yCounter, zCounter;
    int j, k;

    GreyValue * voxelPtr;
    GreyValue * oldPtr;
    unsigned char * xCutPlanePtr;

    if (image == NULL) {
        freePlaneMemory();
        return;
    }

    if (xCutPlane == NULL)
        return;

    voxelPtr = image->getVoxels();
    if (voxelPtr == NULL) {
        delete [] xCutPlane;
        xCutPlane = NULL;
        return;
    }

    xDim = image->getXDim();
    yDim = image->getYDim();
    zDim = image->getZDim();

	xIndex = image->worldXToImage(xCutPlaneWorldPos);
	// Note that X cannot be flipped
#ifdef DEBUG
	cout << "ImagePlanes::fastFillXCutPlane() Filling X slice "
		<< xIndex << ", " << xCutPlaneWorldPos << " cm" << endl;
#endif

    // Set image to black, if we are out of bounds
    if (xIndex < 0 || xIndex >= xDim) {
        memset((void*) xCutPlane, 0, imagePlaneSize*imagePlaneSize);
        return;
    }

    zLeap = xDim * yDim;

    yIncrement = (yDim / imagePlaneSize) * xDim;
    yModulo = yDim % imagePlaneSize;    // Remainder

    zIncrement = (zDim / imagePlaneSize) * zLeap;
    zModulo = zDim % imagePlaneSize;

    xCutPlanePtr = xCutPlane;
    yCounter = 0;
    zCounter = 0;

	if (yFlip) {
        if (zFlip) {
            yIncrement = -yIncrement;
            voxelPtr += (xDim*yDim*zDim - 1) - (xDim - xIndex - 1);
            xDim = -xDim;
            zLeap = -zLeap;
        }
        else {
            yIncrement = -yIncrement;
            voxelPtr += (xDim*yDim - 1) - (xDim - xIndex - 1);
            xDim = -xDim;
        }
	}
	else if (zFlip) {
        zIncrement = -zIncrement;
        voxelPtr += xDim*yDim*(zDim - 1) + xIndex;
        zLeap = -zLeap;
	}
	else
        voxelPtr += xIndex;

    for (k = 0; k < imagePlaneSize; k++) {		// Loop in Z direction
        oldPtr = voxelPtr;						// Save value for reuse
        for (j = 0; j < imagePlaneSize; j++) {		// Loop in Y direction
#ifdef AE2_BUILD
            (*xCutPlanePtr) = colorLUT[IMAGE_BITS & *voxelPtr];
#else
            (*xCutPlanePtr) = colorLUT[*voxelPtr];
#endif
            voxelPtr += yIncrement;		// Next input voxel
            yCounter += yModulo;
            if (yCounter >= imagePlaneSize) {
                yCounter -= imagePlaneSize;
                voxelPtr += xDim;
            }
            xCutPlanePtr++;	// Next output pixel
        }
        voxelPtr = oldPtr;		// Restore starting value
        voxelPtr += zIncrement;		// Compute start of next plane	
        zCounter += zModulo;
        if (zCounter >= imagePlaneSize) {
            zCounter -= imagePlaneSize;
            voxelPtr += zLeap;
        }
    }
}

// Coronal plane, interpolated
void ImagePlanes::fillYCutPlane()
{
    int yIndex;
    double xIndex, zIndex;
    int xDim, zDim;
    int i, j;
    GreyValue voxel;
    GreyValue * voxelPtr;
    unsigned char * yCutPlanePtr;

    if (image == NULL) {
        freePlaneMemory();
        return;
    }

    if (yCutPlane == NULL)
        return;

    voxelPtr = image->getVoxels();
    if (voxelPtr == NULL) {
        delete [] yCutPlane;
        yCutPlane = NULL;
        return;
    }

	yIndex = image->worldYToImage(yCutPlaneWorldPos);
#ifdef DEBUG
	cout << "ImagePlanes::fillYCutPlane() Filling Y slice "
		<< yIndex << ", " << yCutPlaneWorldPos << " cm" << endl;
#endif

    // Set image to black, if we are out of bounds
    if (yIndex < 0 || yIndex >= image->getYDim()) {
        memset((void*) yCutPlane, 0, imagePlaneSize*imagePlaneSize);
        return;
    }

    xDim = image->getXDim() - 1;
    zDim = image->getZDim() - 1;

    yCutPlanePtr = yCutPlane;
    for (i = 0; i < imagePlaneSize; i++) {
		zIndex = ((double) (i*zDim))/(imagePlaneSize - 1);
        for (j = 0; j < imagePlaneSize; j++) {
			xIndex = ((double) (j*xDim))/(imagePlaneSize - 1);
			// For AE2, masking is done in getInterpolatedVoxelValue()
			voxel = (GreyValue) (0.5 + image->getInterpolatedVoxelValue(
                xIndex, (double) yIndex, zIndex));
            (*yCutPlanePtr) = colorLUT[voxel];
            yCutPlanePtr++;
        }
    }
}

// Coronal plane
void ImagePlanes::fastFillYCutPlane()
{
    int yIndex;
    int xDim, yDim, zDim;
    int xIncrement, zIncrement;
    int zLeap;
    int xModulo, zModulo;
    int xCounter, zCounter;
    int i, k;

    GreyValue * voxelPtr;
    GreyValue * oldPtr;
    unsigned char * yCutPlanePtr;

    if (image == NULL) {
        freePlaneMemory();
        return;
    }

    if (yCutPlane == NULL)
        return;

    voxelPtr = image->getVoxels();
    if (voxelPtr == NULL) {
        delete [] yCutPlane;
        yCutPlane = NULL;
        return;
    }

    xDim = image->getXDim();
    yDim = image->getYDim();
    zDim = image->getZDim();

	yIndex = image->worldYToImage(yCutPlaneWorldPos);
	if (yFlip) {
		// Voxel access using Image3D::getVoxelValue() is flipped in Y.
		// This function accesses the voxel array directly, so we must
		// "unflip" here.
		yIndex = yDim - 1 - yIndex;
	}
#ifdef DEBUG
	cout << "ImagePlanes::fastFillYCutPlane() Filling Y slice "
		<< yIndex << ", " << yCutPlaneWorldPos << " cm" << endl;
#endif

    // Set image to black, if we are out of bounds
    if (yIndex < 0 || yIndex >= yDim) {
        memset((void*) yCutPlane, 0, imagePlaneSize*imagePlaneSize);
        return;
    }

    zLeap = xDim * yDim;

    xIncrement = (xDim / imagePlaneSize);
    xModulo = xDim % imagePlaneSize;    // Remainder

    zIncrement = (zDim / imagePlaneSize) * zLeap;
    zModulo = zDim % imagePlaneSize;

	if (zFlip) {
        zIncrement = -zIncrement;
        voxelPtr += xDim*yDim*(zDim - 1);
        zLeap = -zLeap;
	}
    voxelPtr += yIndex * xDim;

    yCutPlanePtr = yCutPlane;
    xCounter = 0;
    zCounter = 0;

    for (k = 0; k < imagePlaneSize; k++) {		// Loop in Z direction
        oldPtr = voxelPtr;						// Save value for reuse
        for (i = 0; i < imagePlaneSize; i++) {		// Loop in X direction
#ifdef AE2_BUILD
            (*yCutPlanePtr) = colorLUT[IMAGE_BITS & *voxelPtr];
#else
            (*yCutPlanePtr) = colorLUT[*voxelPtr];
#endif
            voxelPtr += xIncrement;		// Next input voxel
            xCounter += xModulo;
            if (xCounter >= imagePlaneSize) {
                xCounter -= imagePlaneSize;
                voxelPtr++;
            }
            yCutPlanePtr++;    // Next output pixel
        }
        voxelPtr = oldPtr;		// Restore starting value
        voxelPtr += zIncrement;		// Compute start of next row
        zCounter += zModulo;
        if (zCounter >= imagePlaneSize) {
            zCounter -= imagePlaneSize;
            voxelPtr += zLeap;
        }
    }
}

// Axial plane, interpolated
void ImagePlanes::fillZCutPlane()
{
    int zIndex;
    double xIndex, yIndex;
    int xDim, yDim;
    int i, j;
    GreyValue voxel;
    GreyValue * voxelPtr;
    unsigned char * zCutPlanePtr;

    if (image == NULL) {
        freePlaneMemory();
        return;
    }

    if (zCutPlane == NULL)
        return;

    voxelPtr = image->getVoxels();
    if (voxelPtr == NULL) {
        delete [] zCutPlane;
        zCutPlane = NULL;
        return;
    }

	zIndex = image->worldZToImage(zCutPlaneWorldPos);
#ifdef DEBUG
	cout << "ImagePlanes::fillZCutPlane(): Filling Z slice "
		<< zIndex << ", " << zCutPlaneWorldPos << " cm" << endl;
#endif

    // Set image to black, if we are out of bounds
    if (zIndex < 0 || zIndex >= image->getZDim()) {
        memset((void*) zCutPlane, 0, imagePlaneSize*imagePlaneSize);
        return;
    }

    xDim = image->getXDim() - 1;
    yDim = image->getYDim() - 1;

    zCutPlanePtr = zCutPlane;
    for (i = 0; i < imagePlaneSize; i++) {
		yIndex = ((double) (i*yDim))/(imagePlaneSize - 1);
        for (j = 0; j < imagePlaneSize; j++) {
			xIndex = ((double) (j*xDim))/(imagePlaneSize - 1);
			// For AE2, masking is done in getInterpolatedVoxelValue()
			voxel = (GreyValue) (0.5 + image->getInterpolatedVoxelValue(
                xIndex, yIndex, (double) zIndex));
            (*zCutPlanePtr) = colorLUT[voxel];
            zCutPlanePtr++;
        }
    }
}

// Axial plane
void ImagePlanes::fastFillZCutPlane()
{
    int zIndex;
    int xDim, yDim, zDim;
    int xIncrement, yIncrement;
    int xModulo, yModulo;
    int xCounter, yCounter;
    int i, j;

    GreyValue * voxelPtr;
    GreyValue * oldPtr;
    unsigned char * zCutPlanePtr;

    if (image == NULL) {
        freePlaneMemory();
        return;
    }

    if (zCutPlane == NULL)
        return;

    voxelPtr = image->getVoxels();
    if (voxelPtr == NULL) {
        delete [] zCutPlane;
        zCutPlane = NULL;
        return;
    }

    xDim = image->getXDim();
    yDim = image->getYDim();
    zDim = image->getZDim();

	zIndex = image->worldZToImage(zCutPlaneWorldPos);
	if (zFlip) {
		// Voxel access using Image3D::getVoxelValue() is flipped in Z.
		// This function accesses the voxel array directly, so we must
		// "unflip" here.
		zIndex = zDim - 1 - zIndex;
	}
#ifdef DEBUG
	cout << "ImagePlanes::fastFillZCutPlane(): Filling Z slice "
		<< zIndex << ", " << zCutPlaneWorldPos << " cm" << endl;
#endif

    // Set image to black, if we are out of bounds
    if (zIndex < 0 || zIndex >= zDim) {
        memset((void*) zCutPlane, 0, imagePlaneSize*imagePlaneSize);
        return;
    }

    xIncrement = (xDim / imagePlaneSize);
    xModulo = xDim % imagePlaneSize;    // Remainder

    yIncrement = (yDim / imagePlaneSize) * xDim;
    yModulo = yDim % imagePlaneSize;

	if (yFlip) {
        yIncrement = -yIncrement;
        voxelPtr += xDim*yDim*zIndex + xDim*(yDim - 1);
		xDim = -xDim;
	}
	else
		voxelPtr += zIndex * xDim * yDim;

    zCutPlanePtr = zCutPlane;
    xCounter = 0;
    yCounter = 0;

    for (j = 0; j < imagePlaneSize; j++) {		// Loop in Y direction
        oldPtr = voxelPtr;						// Save value for reuse
        for (i = 0; i < imagePlaneSize; i++) {		// Loop in X direction
#ifdef AE2_BUILD
            (*zCutPlanePtr) = colorLUT[IMAGE_BITS & *voxelPtr];
#else
            (*zCutPlanePtr) = colorLUT[*voxelPtr];
#endif
            voxelPtr += xIncrement;		// Next input voxel
            xCounter += xModulo;
            if (xCounter >= imagePlaneSize) {
                xCounter -= imagePlaneSize;
                voxelPtr++;
            }
            zCutPlanePtr++;		// Next output pixel
        }
        voxelPtr = oldPtr;		// Restore starting value
        voxelPtr += yIncrement;		// Compute start of next row
        yCounter += yModulo;
        if (yCounter >= imagePlaneSize) {
            yCounter -= imagePlaneSize;
            voxelPtr += xDim;
        }
    }
}

void ImagePlanes::fillArbCutPlane()
{

}

void ImagePlanes::fastFillArbCutPlane()
{

}

void ImagePlanes::allocatePlaneMemory()
{
    int totalSize;

#ifdef DEBUG
    cout << "ImagePlanes::allocatePlaneMemory()\n";
#endif
    totalSize = imagePlaneSize * imagePlaneSize;

    if (xCutPlane != NULL)
        delete [] xCutPlane;
    if (yCutPlane != NULL)
        delete [] yCutPlane;
    if (zCutPlane != NULL)
        delete [] zCutPlane;
    if (arbCutPlane != NULL)
        delete [] arbCutPlane;

    if (totalSize <= 0) {
        xCutPlane = NULL;
        yCutPlane = NULL;
        zCutPlane = NULL;
        arbCutPlane = NULL;
    }
    else {
        xCutPlane = new unsigned char[totalSize];
        yCutPlane = new unsigned char[totalSize];
        zCutPlane = new unsigned char[totalSize];
        arbCutPlane = new unsigned char[totalSize];
    }
}

void ImagePlanes::freePlaneMemory()
{
#ifdef DEBUG
    cout << "ImagePlanes::freePlaneMemory()\n";
#endif
    if (xCutPlane != NULL) {
        delete [] xCutPlane;
        xCutPlane = NULL;
    }
    if (yCutPlane != NULL) {
        delete [] yCutPlane;
        yCutPlane = NULL;
    }
    if (zCutPlane != NULL) {
        delete [] zCutPlane;
        zCutPlane = NULL;
    }
    if (arbCutPlane != NULL) {
        delete [] arbCutPlane;
        arbCutPlane = NULL;
    }
}


