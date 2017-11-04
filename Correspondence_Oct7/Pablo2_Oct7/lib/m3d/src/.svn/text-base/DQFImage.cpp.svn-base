#include <iostream>
#include "DQFImage.h"


using namespace std;


DQFImage::DQFImage() : Image3D()
{
	xDim = 0;
	yDim = 0;
	zDim = 0;
	xOrigin = 0;
	yOrigin = 0;
	zOrigin = 0;
	xBound = 0;    // The bound is only stored for speed in voxel access
	yBound = 0;
	zBound = 0;
	aDim = 0;
	bDim = 0;
	cDim = 0;
	nVoxels = 0;
	sWidth = 0;
	multiplier = 1.0;
	delta = 0;
	Image3D::setModality(DQF);
}

DQFImage::DQFImage(int rX, int rY, int rZ, int w, int iX, int iY, int iZ)
		: Image3D()
{
	setDims(rX, rY, rZ, w, iX, iY, iZ);
	xOrigin = 0;
	yOrigin = 0;
	zOrigin = 0;
	xBound = 0;
	yBound = 0;
	zBound = 0;
	multiplier = 1.0;
	delta = 0;
	Image3D::setModality(DQF);
}

void DQFImage::setDims(int rX, int rY, int rZ, int w,
    int iX, int iY, int iZ)
{
	if (Image3D::getVoxelCount() != 0) {
		cout << "Error: Dimensions of ROI can only be specified once" << endl;
		return;
	}
	xDim = iX;
	yDim = iY;
	zDim = iZ;
	aDim = rX;
	bDim = rY;
	cDim = rZ;
	sWidth = w;
	nVoxels = rX*rY*rZ;
	GreyValue * voxels = new GreyValue[nVoxels*w];
	Image3D::setVoxels(voxels, rX*w, rY, rZ);    // Construct other arrays
	// DQF images are considered to have the full grey range, even if
	// they do not.  This is because they must be written to disk using
	// this range, so that scaling and shifting do not occur when they are
	// read.
	Image3D::setRange(MIN_GREY_VALUE, MAX_GREY_VALUE);
}

bool DQFImage::roi(int x1, int y1, int z1, int x2, int y2, int z2)
{
	int n;

	n = (x2 - x1 + 1)*(y2 - y1 + 1)*(z2 - z1 + 1);
	if (nVoxels != 0 && n != nVoxels) {
		cout << "Error: DQF ROI's size is incorrect" << endl;
		return false;
	}

	if (x1 < 0 || y1 < 0 || z1 < 0 || x2 >= xDim || y2 >= yDim || z2 >= zDim) {
		cout << "Error: DQF ROI is outside of the separate image" << endl;
		return false;
	}

	xOrigin = x1;
	yOrigin = y1;
	zOrigin = z1;
	xBound = x2;
	yBound = y2;
	zBound = z2;
	aDim = x2 - x1 + 1;
	bDim = y2 - y1 + 1;
	cDim = z2 - z2 + 1;
	return true;
}

bool DQFImage::roi(Vector3D & r1, Vector3D & r2)
{
	return roi((int) r1.getX(), (int) r1.getY(), (int) r1.getZ(),
		(int) r2.getX(), (int) r2.getY(), (int) r2.getZ());
}

void DQFImage::getROI(Vector3D & r1, Vector3D & r2)    // Query the ROI
{
	r1.set((double) xOrigin, (double) yOrigin, (double) zOrigin);
	r2.set((double) xBound, (double) yBound, (double) zBound);
}

void DQFImage::clear(GreyValue val, bool reset)
{
    for (int i = 0; i < nVoxels*sWidth; i++)
	    voxelArray[i] = val;

    if (! reset)
	return;
    scale_factor = 1.0;
    intens_shift = 0;

    minIntens = 1;
    maxIntens = 0;
    intensRange = 0;
}

void DQFImage::setVoxels(GreyValue * voxels, int rX, int rY, int rZ, int w)
{
	// Alternate way to construct the base image
	Image3D::setVoxels(voxels, rX*w, rY, rZ, true);
	nVoxels = rX*rY*rZ;    // Prevents constructor call in roi() above
}

// The coordinates are ROI-relative
void DQFImage::setVoxel(int x, int y, int z, GreyValue * val) {
	for (int n = 0; n < sWidth; n++)
		Image3D::setVoxel(x*sWidth + n, y, z, val[n]);
}

void DQFImage::setVoxel(int x, int y, int z, double * val)
{
	for (int n = 0; n < sWidth; n++) {
		int voxel = (int) (val[n]*multiplier + (val[n] < 0 ? -0.5 : 0.5));
		Image3D::setVoxel(x*sWidth + n, y, z, (GreyValue) (voxel - delta));
	}
}

// The coordinates are ROI-relative
void DQFImage::getVoxelValue(int x, int y, int z, double * voxel)
{
	for (int n = 0; n < sWidth; n++) {
		int val;
		val = (int) Image3D::getVoxelValue(x*sWidth + n, y, z);
		voxel[n] = (val + delta)/multiplier;
	}
}

// The coordinates are ROI-relative
void DQFImage::getVoxelValue(int x, int y, int z, GreyValue * voxel)
{
	for (int n = 0; n < sWidth; n++)
		voxel[n] = Image3D::getVoxelValue(x*sWidth + n, y, z);
}

void DQFImage::getImageVoxelValue(int x, int y, int z, GreyValue * voxel)
{
	if (x < xOrigin || y < yOrigin || z < zOrigin
			|| x > xBound || y > yBound || z > zBound)
	{
		for (int n = 0; n < sWidth; n++)
			voxel[n] = 0;
	}
	else {
		x -= xOrigin;
		y -= yOrigin;
		z -= zOrigin;
		for (int n = 0; n < sWidth; n++)
			voxel[n] = Image3D::getVoxelValue(x*sWidth + n, y, z);
	}
}

void DQFImage::getImageVoxelValue(int x, int y, int z, double * voxel)
{
	if (x < xOrigin || y < yOrigin || z < zOrigin
			|| x > xBound || y > yBound || z > zBound)
	{
		for (int n = 0; n < sWidth; n++)
			voxel[n] = 0;
	}
	else {
		x -= xOrigin;
		y -= yOrigin;
		z -= zOrigin;
		for (int n = 0; n < sWidth; n++) {
			int val;
			val = (int) Image3D::getVoxelValue(x*sWidth + n, y, z);
			voxel[n] = (val + delta)/multiplier;
		}
	}
}

//I don't know about error checking this in case the neighbor voxels aren't in bounds...
void DQFImage::getInterpolatedImageVoxelValue(double x, double y, double z, double * voxel)
{
    if (x < xOrigin || y < yOrigin || z < zOrigin
			|| x > xBound || y > yBound || z > zBound)
	{
		for (int n = 0; n < sWidth; n++)
			voxel[n] = 0;
	}
	else {
		x -= xOrigin;
		y -= yOrigin;
		z -= zOrigin;
		
		int xIndex,
			yIndex,
			zIndex;
		
		double xt, yt, zt;
		
		double v[8];
		
		xIndex = (int) x;
		yIndex = (int) y;
		zIndex = (int) z;
		
		xt = (double) x - (double) xIndex;
		yt = (double) y - (double) yIndex;
		zt = (double) z - (double) zIndex;
		
		for (int n = 0; n < sWidth; n++)
		{
			v[0] = (double) Image3D::getVoxelValue(xIndex*sWidth + n, yIndex, zIndex);
			v[0] = (v[0] + delta)/multiplier;

			v[1] = (double) Image3D::getVoxelValue((xIndex + 1)*sWidth + n, yIndex, zIndex);
			v[1] = (v[1] + delta)/multiplier;

			v[2] = (double) Image3D::getVoxelValue(xIndex*sWidth + n, yIndex + 1, zIndex);
			v[2] = (v[2] + delta)/multiplier;

			v[3] = (double) Image3D::getVoxelValue((xIndex + 1)*sWidth + n, yIndex + 1, zIndex);
			v[3] = (v[3] + delta)/multiplier;

			v[4] = (double) Image3D::getVoxelValue(xIndex*sWidth + n, yIndex, zIndex + 1);
			v[4] = (v[4] + delta)/multiplier;

			v[5] = (double) Image3D::getVoxelValue((xIndex + 1)*sWidth + n, yIndex, zIndex + 1);
			v[5] = (v[5] + delta)/multiplier;

			v[6] = (double) Image3D::getVoxelValue(xIndex*sWidth + n, yIndex + 1, zIndex + 1);
			v[6] = (v[6] + delta)/multiplier;

			v[7] = (double) Image3D::getVoxelValue((xIndex + 1)*sWidth + n, yIndex + 1, zIndex + 1);
			v[7] = (v[7] + delta)/multiplier;
			
			// This works on voxel centers:  if xt, yt and zt are 0, no
			// interpolation occurs.
			voxel[n] = (1 - zt) * ((1 - yt) * (v[0] * (1 - xt) + v[1] * xt)
				+ yt * (v[2] * (1 - xt) + v[3] * xt))
				+ zt * ((1 - yt) * (v[4] * (1 - xt) + v[5] * xt)
				+ yt * (v[6] * (1 - xt) + v[7] * xt));
			
		}
		
	}
    
}
