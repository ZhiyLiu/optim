#include <iostream>

#include "ImageResample3D.h"
#include "ControlParms.h"


#include <string.h>

using namespace std;


void ImageResample3D::isotropicSample(Image3D & image, double userSpacing,
	GreyValue attractionMask, GreyValue repulsionMask)
{
    GreyValue * newData;
    double spacing, maxSpacing;
    Vector3D scale;
    double xImageSpacing,
           yImageSpacing,
           zImageSpacing;
    int xDim, yDim, zDim;
    int index, jIndex, kIndex;
	int jStart, jIncr, kStart, kIncr;
    int i, j, k;
    double x, y, z;
	int sign[3];

	if (userSpacing <= 0.0) {
		// Find minimum voxel dimension (spacing is in world coordinates)
		spacing = image.getAbsXSpacing();
		if (spacing > image.getAbsYSpacing())
			spacing = image.getAbsYSpacing();
		if (spacing > image.getAbsZSpacing())
			spacing = image.getAbsZSpacing();
		if (globalVerbosity >= 0)
			cout << "Resampling at voxel size: " << spacing << endl;
	}
	else {
		spacing = userSpacing;
#ifdef BINARY
		cout << "Resampling at USER'S voxel size: " << spacing << endl;
#endif
	}

    // Find maximum voxel dimension for isotropicity check
    maxSpacing = image.getAbsXSpacing();
    if (maxSpacing < image.getAbsYSpacing())
        maxSpacing = image.getAbsYSpacing();
    if (maxSpacing < image.getAbsZSpacing())
        maxSpacing = image.getAbsZSpacing();

	if (globalVerbosity >= 1) {
		cout << "Old spacing: ";
		cout << image.getAbsXSpacing() << " " << image.getAbsYSpacing() << " "
			<< image.getAbsZSpacing() << '\n';

		if (spacing != maxSpacing)
			cout << "New spacing: " << spacing << endl;
	}

	if (spacing == maxSpacing) {
		cout << "Skipping resampling: image is already isotropic" << endl;
		return;
	}

    // Put spacing into the original image's coordinates
    xImageSpacing = image.worldXToImageDistance(spacing);
    yImageSpacing = image.worldYToImageDistance(spacing);
    zImageSpacing = image.worldZToImageDistance(spacing);

	if (globalVerbosity >= 1) {
		cout << "Old image dimensions: ";
		cout << image.getXDim() << " " << image.getYDim() << " " << image.getZDim() << endl;
	}

    // Get new number of voxels in each dimension
#ifdef BINARY										/* AGG: this block probably should be discarded in favor of the #else block */
    xDim = (int) (0.5 + (((double) image.getXDim())*image.getAbsXSpacing()/spacing));
    yDim = (int) (0.5 + (((double) image.getYDim())*image.getAbsYSpacing()/spacing));
    zDim = (int) (0.5 + (((double) image.getZDim())*image.getAbsZSpacing()/spacing));
#else
    xDim = (int) (((double) image.getXDim())*image.getAbsXSpacing()/spacing);
    yDim = (int) (((double) image.getYDim())*image.getAbsYSpacing()/spacing);
    zDim = (int) (((double) image.getZDim())*image.getAbsZSpacing()/spacing);
#endif
	if (globalVerbosity >= 1) {
		cout << "New image dimensions: ";
		cout << xDim << " " << yDim << " " << zDim << endl;
	}

    newData = new GreyValue[xDim * yDim * zDim];
    if (newData == NULL) {
        cout << "Error: Could not allocate space for new image" << endl;
        return;
    }

    // Resample
	// Notice that the assignment of voxels must be done in a way
	// that maintains any negative signs in the original spacings.
	// This is because the corresponding indexing look up tables in
	// the Image3D object will be flipped.
    image.getSign(sign);
	if (sign[1] < 0) {
		jStart = yDim - 1;
		jIncr = -1;
	}
	else {
		jStart = 0;
		jIncr = 1;
	}
	if (sign[2] < 0) {
		kStart = zDim - 1;
		kIncr = -1;
	}
	else {
		kStart = 0;
		kIncr = 1;
	}
	jIndex = jStart;
	kIndex = kStart;
    for (k = 0; k < zDim; k++) {
		if (globalVerbosity >= 1 && (k & 0100))
			cout << '.' << flush;
        z = ((double) k) * zImageSpacing;

        for (j = 0; j < yDim; j++) {
            y = ((double) j) * yImageSpacing;

            for (i = 0; i < xDim; i++) {
                x = ((double) i) * xImageSpacing;

				index = i + (jIndex + kIndex*yDim)*xDim;
#ifdef BINARY							/* AGG: This block compares using <, while the else block uses <=; discard this block? */
				double d = image.getWindowedInterpolatedVoxelValue(x, y, z);
				newData[index] = (d > 0.5);
#else
                newData[index] = (GreyValue)
					(0.5 + image.getWindowedInterpolatedVoxelValue(x, y, z));
#endif
            }
			jIndex += jIncr;
        }
		jIndex = jStart;
		kIndex += kIncr;
    }

#ifdef BINARY
	jIndex = jStart;
	kIndex = kStart;
	if (image.getIsImageStacked()) {
		// When resampling a stacked image,
		// assume we just processed the attraction plane(s) as a binary.
		// Now shift it over and process the retraction plane(s) as well.
		image.pushImageIsStacked(true, repulsionMask);

        for (k = 0; k < zDim; k++) {
			if (globalVerbosity >= 1 && (k & 0100))
				cout << '.' << flush;
			z = ((double) k) * zImageSpacing;

            for (j = 0; j < yDim; j++) {
				y = ((double) j) * yImageSpacing;

				for (i = 0; i < xDim; i++) {
					x = ((double) i) * xImageSpacing;

					index = i + (jIndex + kIndex*yDim)*xDim;
					newData[index] *= attractionMask;
					double d = image.getWindowedInterpolatedVoxelValue(x, y, z);
					newData[index] += repulsionMask * (d > 0.5);
				}
				jIndex += jIncr;
			}
			jIndex = jStart;
			kIndex += kIncr;
		}
		image.popImageIsStacked();
	}
#endif

    // Set the new voxel data (Image3D takes care of deleting old data)
    image.setVoxels(newData, xDim, yDim, zDim);
		    //	AGG: The following line is ignoring the origin.  Should it?  It also needs the signs.
	Vector3D origin = image.getWorldOrigin();
    image.setSpacingAndOrigin(spacing, sign[1]*spacing, sign[2]*spacing, &origin);
	if (globalVerbosity >= 1)
		cout << " done" << endl;
}

#ifdef BINARY

int ImageResample3D::cubeSample(Image3D & image, float maxExtent)
{
	float xExtent, yExtent, zExtent, realMaxExtent;

	xExtent = image.getXExtent();
	yExtent = image.getYExtent();
	zExtent = image.getZExtent();

	realMaxExtent = xExtent;
	if (yExtent < realMaxExtent)
		realMaxExtent = yExtent;
	if (zExtent < realMaxExtent)
		realMaxExtent = zExtent;

	if (globalVerbosity >= 1)
		cout << " maxExtent is " << realMaxExtent <<
			" (world coords), from image header" << endl;

	if (maxExtent != 0) {		// verify user's maxExtent
		if (realMaxExtent > maxExtent) {
			cout << "Error: specified extent (" << maxExtent <<
				") is shorter than a dimension's extent (" << realMaxExtent
				<< ")" << endl;
			return 0;
		}
		if (globalVerbosity >= 1)
			cout << " user overrides maxExtent with " << maxExtent
				<< " (world coords)" << endl;
	}
	else
		maxExtent = realMaxExtent;

    int newxDim = (int) (.5 + (((double) image.getXDim()) * maxExtent / xExtent));
    int newyDim = (int) (.5 + (((double) image.getYDim()) * maxExtent / yExtent));
    int newzDim = (int) (.5 + (((double) image.getZDim()) * maxExtent / zExtent));

	GreyValue *newData = new GreyValue[newxDim * newyDim * newzDim];
    if (! newData) {
        cout << "Error: Could not allocate space for new image." << endl;
        return 0;
    }

    // pad with 0 intensity voxels
	memset(newData, 0, sizeof(GreyValue) * newxDim * newyDim * newzDim);

	int oldxDim = image.getXDim();
	int oldyDim = image.getYDim();
	int oldzDim = image.getZDim();
	int oldxyDim = oldxDim * oldyDim;

	int newxyDim = newxDim * newyDim;

	if (oldxDim >= newxDim && oldyDim >= newyDim && oldzDim >= newzDim)
		return 1;		// nothing to do

	if (globalVerbosity >= 1)
		cout << " cubing extents by padding dims (" << oldxDim << " " << oldyDim
			<< " " << oldzDim << ") up to (" << newxDim << " " << newyDim
			<< " " << newzDim << ")" << endl;

	GreyValue *oldData = image.getVoxels();
    for (int z = 0; z < oldzDim; z++) {
        for (int y = 0; y < oldyDim; y++) {
			// fast-copy each scan line as a block because it's guaranteed contiguous
			memcpy(newData + z*newxyDim + y*newxDim,
				   oldData + z*oldxyDim + y*oldxDim,
				   oldxDim * sizeof(GreyValue));
        }
    }

    // Set the new voxel data (Image3D takes care of deleting old data)
    image.setVoxels(newData, newxDim, newyDim, newzDim);
	return 1;
}

#endif

bool ImageResample3D::superSample(Image3D & image, int outGridDims)
{
    GreyValue * newData;
    double xspacing, yspacing, zspacing;
    Vector3D worldToModelScale;
    int xDim, yDim, zDim, mdim;
    int index, jIndex, kIndex;
	int jStart, jIncr, kStart, kIncr;
    int i, j, k;
    double x, y, z;
	int sign[3];
	double stepsize, voxval;


    if (outGridDims<= 0)
        return false;

//	report_intensity_range(image);
	if (globalVerbosity >= 0)
		cout << "Super-sampling to dimensions of : " << outGridDims << endl;

	xDim = image.getXDim();
	yDim = image.getYDim();
	zDim = image.getZDim();

	if (xDim > yDim) {
	    if (xDim > zDim)
			mdim = xDim;
		else
			mdim = zDim;
	}
	else if (yDim > zDim)
		mdim = yDim;
	else
		mdim = zDim;

	stepsize = (double) mdim/outGridDims;

	xDim = (int) (0.5 + ((double)xDim/stepsize));
	yDim = (int) (0.5 + ((double)yDim/stepsize));
	zDim = (int) (0.5 + ((double)zDim/stepsize));

    newData = new GreyValue[xDim * yDim * zDim];
    if (newData == NULL) {
        cerr << "Error: Could not allocate space for resampled image." << endl;
        return false;
    }

    // Resample
	// Notice that the assignment of voxels must be done in a way
	// that maintains any negative signs in the original spacings.
	// This is because the corresponding indexing look up tables in
	// the Image3D object will be flipped.
    image.getSign(sign);
	if (sign[1] < 0) {
		jStart = yDim - 1;
		jIncr = -1;
	}
	else {
		jStart = 0;
		jIncr = 1;
	}
	if (sign[2] < 0) {
		kStart = zDim - 1;
		kIncr = -1;
	}
	else {
		kStart = 0;
		kIncr = 1;
	}
	jIndex = jStart;
	kIndex = kStart;

	z = -stepsize; 
    for (k = 0; k < zDim; k++) {
		if (globalVerbosity >= 1 && (k & 0100))
			cout << '.' << flush;
        z += stepsize;

		y = -stepsize;
        for (j = 0; j < yDim; j++) {
            y += stepsize;

			x = -stepsize;
            for (i = 0; i < xDim; i++) {
                x += stepsize;

				voxval = image.getVoxelValue((int) x, (int) y, (int) z);
				index = i + (jIndex + kIndex*yDim)*xDim;
                newData[index] = (GreyValue) (0.5 + voxval);
            }
			jIndex += jIncr;
        }
		jIndex = jStart;
		kIndex += kIncr;
    }

    // Set the new voxel data
    image.setVoxels(newData, xDim, yDim, zDim);

	xspacing = image.getXSpacing()*stepsize;
	yspacing = image.getYSpacing()*stepsize;
	zspacing = image.getZSpacing()*stepsize;

	Vector3D origin = image.getWorldOrigin();
	image.setSpacingAndOrigin(xspacing, yspacing, zspacing, &origin);
	if (globalVerbosity >= 1)
		cout << " done" << endl;

//	report_intensity_range(image);
	return true;
}


