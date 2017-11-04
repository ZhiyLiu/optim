#ifndef DQF_IMAGE_H
#define DQF_IMAGE_H

#include "Image3D.h"


/*  Class DQFImage.

    This class provides a way of storing an ROI (region of interest) of
	a spectral image based on another Image3D object.  The dimensions
	of the larger, separate image must be specified, and the usual call
	to specify the world coordinate system for the larger image should
	be made.  However, in addition, the spectral width of each voxel,
	the number of grey values per voxel must be provided.  Also, one of
	the roi() functions provided must be called to specify the region of
	interest.  The data stored will actually only be that for the ROI.
	The width of each voxel will be accounted for by extending the X
	dimension of the base image by a factor of the width.

    In order for these images to be transmitted to and from disk, they
	must be given the modality DQF and they must be accessed using this
	class.  Otherwise they will not make sense.
*/

class DQFImage : public Image3D
{
public:

    DQFImage();
	// Specify ROI's dimensions and spectral width (items/voxel).
	// Optionally specify the dimensions of the separate image, which are
	// considered invariant by this class.
    DQFImage(int rX, int rY, int rZ, int w, int iX = 0, int iY = 0,
		int iZ = 0);
    ~DQFImage() { }

	// Specify the ROI dimensions and spectral width (items/voxel),
	// followed by the larger image's dimensions, which are considered
	// invariant by this class.  Use this function if the default
	// constructor was used.
	void setDims(int rX, int rY, int rZ, int w, int iX = 0, int iY = 0,
		int iZ = 0);

	// Specify the relationship between actual DQF data and grey values
	// stored in this class for I/O.  The conversions are:
	//    greyEntry = round(doubleDQF*scale) - shift
	//    doubleDQF = (greyEntry + shift)/scale
	// This should only be called once before setting any voxels.  If
	// not called a scale of 1.0 and shift of 0 are assumed.
	void setMapping(double scale, int shift = 0) {
		multiplier = scale;
		delta = shift;
	}
	void mapping(double & scale, int & shift) {
		scale = multiplier;
		shift = delta;
	}

	// Set the ROI after construction.  This should  be called at least
	// once immediately after construction.  Thereafter, it may be called
	// at any time to change the ROI.  The two points must both be inside
	// the separate image.  Note that changing the ROI does not causes the
	// stored intensities to move to other voxels.
	bool roi(int x1, int y1, int z1, int x2, int y2, int z2);
	// Note: the points are integers, not doubles, in the function below.
	bool roi(Vector3D & r1, Vector3D & r2);

	// Number of slices per Cartesian direction in the ROI
    int getROIXDim() const { return aDim; }
    int getROIYDim() const { return bDim; }
    int getROIZDim() const { return cDim; }
	Vector3D getROIDims() const {
		return Vector3D((double) aDim, (double) bDim, (double) cDim);
	}

	// Number of items per voxel
	int getWidth() const { return sWidth; }

	// Location of the current ROI on the separate image.  The vectors
	// returned contain integers.
	void getROI(Vector3D & r1, Vector3D & r2);
	Vector3D getROIOrigin() const {
		return Vector3D((double) xOrigin, (double) yOrigin, (double) zOrigin);
	}
	Vector3D getROIBound() const {
		return Vector3D((double) xBound, (double) yBound, (double) zBound);
	}

	// Set DQF values of voxels relative to the ROI origin
    void setVoxel(int x, int y, int z, double * val);

	// Set scaled and shifted voxels relative to the ROI origin
    void setVoxel(int x, int y, int z, GreyValue * val);

	// The voxels array must contain the entire spectral ROI stored as
	// scaled and shifted grey values.  Argument w must be spectral width
	// of each voxel.
    void setVoxels(GreyValue * voxels, int rX, int rY, int rZ, int w);
	// The getVoxels() function of the base class may be used, but the
	// user must be careful to properly index the array.  Note the
	// cautionary comment on the declaration in Image3D.h.

	// The number of voxels in the ROI
	int getVoxelCount() const { return nVoxels; }

	// Obtain DQF values of voxels by slice index of the ROI.
	// The array must have length equal to the spectral width.
    void getVoxelValue(int x, int y, int z, double * voxel);

	// Obtain scaled and shifted voxels by slice index of the ROI.
	// The array must have length equal to the spectral width.
    void getVoxelValue(int x, int y, int z, GreyValue * voxel);

	// These two functions are the same as the two above, except that they
	// use indexes based on the larger, separate image.  That image's
	// dimensions must have been provided (above).
    void getImageVoxelValue(int x, int y, int z, double * voxel);
    void getImageVoxelValue(int x, int y, int z, GreyValue * voxel);

	void getInterpolatedImageVoxelValue(double x, double y, double z, double * voxel);

	// These functions return the dimensions of the separate image
	int getXDim() const { return xDim; }
    int getYDim() const { return yDim; }
    int getZDim() const { return zDim; }

	// These functions are intentionally disabled
	int imageIndex(int x, int y, int z) const { return -1; }
	int imageIndex(int slices[3]) const { return -1; }
	bool replaceVoxels(GreyValue *) { return false; }

	// Set the image intensities to val.  If reset is true, the intensity
        // mapping variables will be reset.
	void clear(GreyValue val = 0, bool reset = true);

	void makeBinary(GreyValue, GreyValue) { }
    void getInterpolatedVoxelValue(double x, double y, double z,
		double * voxel) { }
    GreyValue getWindowedVoxelValue(int, int, int) { return 0; }
    double getWindowedInterpolatedVoxelValue(double, double, double)
	{
		return 0.0;
	}
    void intensityWindow(double, double) { }
    void setModality(modality_t modality) { }	// The modality is DQF
	void setIsImageStacked(bool) { }	// DQF images are not stacked
	void setStackedMask(GreyValue) { }
	void pushImageIsStacked(bool, GreyValue) { }

protected:

	int nVoxels;	// Number of spectral voxels stored
	int sWidth;		// Number of items per spectral voxel
	int xDim, yDim, zDim;		// Separate image's dimensions
	int aDim, bDim, cDim;		// ROI dimensions
	int xOrigin, yOrigin, zOrigin;	// ROI origin (first voxel)
	int xBound, yBound, zBound;     // ROI bound (last voxel)
	double multiplier;
	int delta;

private:

	DQFImage(DQFImage & image);	// Not implemented
	DQFImage operator=(DQFImage & image);	// Not implemented
};


void dqfTest();


#endif


