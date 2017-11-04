#ifndef IMAGE_DISTANCE_MAP_3D_H
#define IMAGE_DISTANCE_MAP_3D_H

#include "DistanceMap3D.h"
#include "Image3D.h"

//xiaojie

#include "Curvature.h"


/*  Class ImageDistanceMap.

This class produces 3D distance maps from a binary image containing
a single figure.

To use this class, the image to which the model is being fit must be
provided, either using a constructor of by the initialize function.
Then createMap() may be called.

Thereafter, the various access functions can be used to get values
from the map.  All distances computed are positive, whether inside
or outside the object.  Thus, on the boundary of the object, the
distances equal zero.  None of the getDistance*() functions are
interpolating.

Distances are represented as fixed point 12.3 integers, that is,
the bottom 3 bits are to the right of the binary point and the top 12
bits are to the left. Another way to think of this is that the stored
value is 8 times the actual voxel distance, eg, a value of 8 represents
1 voxel, 4 means 1/2 voxel, -12 means -1.5 voxels.

The map produced will be the same size as the input image.
*/



class ImageDistanceMap
{
public:

	ImageDistanceMap();
	ImageDistanceMap(Image3D * im);
	virtual ~ImageDistanceMap();

	void initialize();
	void initialize(Image3D * im);

	// Distance map creation functions
	bool createMap();
	bool status() { return distIm != NULL; }	// True if a good map exists

	// Xiaojie, Dibyendu
	// Create a distance map gradient from three gradient images

	void createGradDistMapFromGradXYZ( Image3D* gradX, Image3D* gradY, Image3D* gradZ ) ;


	// Returns array of 3 axial voxel lengths
	const int * getMapSize() { return mapSize; }
	const double * getMapSpacing() { return mapSpacing; }
	const double * getModelToImageScale() { return modelToImageScale; }

	// Return distance to surface in model coordinates at an image coordinate;
	// use city-block distance if off-the-image.
	float getDistance(int i, int j, int k);
	//xiaojie
	Vector3D getGradDistance(int i, int j, int k);

	// These functions return the distance to the surface in model coordinates
	// at a model coordinate. If the point is outside the unit cube, a city block
	// approximation to the distance is returned.
	float getDistance(double x, double y, double z) {
		return get_distance(x, y, z);
	}
	float getDistance(Vector3D point) {
		return getDistance(point.getX(), point.getY(), point.getZ());
	}

	//xiaojie
	Vector3D getGradDistance(double x, double y, double z) {
		return get_graddistance(x, y, z);
	}

	//xiaojie
	Vector3D getGradDistance(Vector3D point) {
		return getGradDistance(point.getX(), point.getY(), point.getZ());
	}

	//xiaojie

	// function to get the curvature at a point on the distance map (level surface)

	Curvature getCurvature(Vector3D point) {
		return getCurvature(point.getX(), point.getY(), point.getZ());
	}

	// function to get the curvature at a point on the distance map (level surface)

	Curvature getCurvature(double x, double y, double z);

	//xiaojie

	// function to get the gradient of curvature at a point on the distance map (level surface)

	bool getGradCurvature(Vector3D point, Vector3D& gradKappa1, Vector3D& gradKappa2) {
		return getGradCurvature(point.getX(), point.getY(), point.getZ(), gradKappa1, gradKappa2);
	}

	// function to get the gradient of curvature at a point on the distance map (level surface)

	bool getGradCurvature(double x, double y, double z, Vector3D& gradKappa1, Vector3D& gradKappa2);

	//xiaojie
	void setgradXDistImfromImage3D(Image3D * gradXdistImage);
    void setgradYDistImfromImage3D(Image3D * gradYdistImage);
    void setgradZDistImfromImage3D(Image3D * gradZdistImage);

// Return the distance to the surface in world coordinates at a model
// coordinate, or INFINITE_DISTANCE, if the argument lies outside the
// unit cube.
float getWorldDistance(double x, double y, double z);

operator Image3D *();
void fromImage3D(Image3D * distMap);
void printMinMaxDistance();

private:

	bool createRawSignedImageDistanceMaps(char * bdryIm);

	// is voxel @ (x,y,z) on object boundary?
	// CAREFUL: (x,y,z) is index into voxel array here: read the source!
	int zeroNeighbors(int x, int y, int z);

	float bilerp_in_z(float x, float y, int z);
	float trilerp(float x, float y, float z);
	float get_distance(double x, double y, double z);
	//xiaojie
	Vector3D bilerp_in_z_grad(float x, float y, int z);
	Vector3D trilerp_grad(float x, float y, float z);
	Vector3D get_graddistance(double x, double y, double z);

	Image3D * image;
	short * distIm;

	//xiaojie
	//Assume they are normalized?
	Vector3D * gradDistIm;

	//short * gradDistIm_X ;
	//short * gradDistIm_X ;
	//short * gradDistIm_X ;

	int mapSize[3];
	int size;
	double scale;
	double mapSpacing[3];
	double modelToImageScale[3];

};


#endif

