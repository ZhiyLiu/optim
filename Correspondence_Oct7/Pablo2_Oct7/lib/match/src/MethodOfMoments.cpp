#include "Mathdefs.h"
#include "Vector3D.h"
#include "Quat.h"
#include "SimilarityTransform3D.h"
#include "M3DObject.h"
#include "Image3D.h"
#include "matrix.h"

#define D_XFERLIST
#define D_POINTLIST_SERVER2
#include "Shapedepend.h"

#include "MethodOfMoments.h"
#include "Tuning.h"


//#define BINARY_PABLO_REPORT
//#define DEBUG

using namespace std;


#ifdef BINARY
void MethodOfMoments::objectMoments(M3DObject * object, double & volume, Vector3D & center, Quat & rotation)
{
	// For BINARY_PABLO, this stuff only works on single figure constructs
	M3DFigure * fig = object->getFigurePtr((int) tuningWt(BpFigureId));
	MethodOfMoments::figureMoments(fig, volume, center, rotation);
}
#endif

void MethodOfMoments::figureMoments(M3DFigure * figure, double & volume, Vector3D & center, Quat & rotation)
{
	/*  This was largely swiped/adapted from Tom's DamFit Method of Moments
		Based on Stokes Theorem/Divergence Theorem

		The idea:

		  Integral (Divergence of some vector function times volume)
		 over space

        = Integral (Function times area)
         over boundary

		=   Sum    (Function times area)     (Discrete representation)
		 over tiles

	    If we make the the divergence of some function:

		1. (grad F) = unity, we get the VOLUME.
		   Thus F MAY be [x 0 0], so dF/dx = [1 0 0] --> Tom's choice
		   or F MAY be [x/3 y/3 z/3], so dF/dx = [1/3 1/3 1/3] --> Sarang's choice, implemented here

		2. (grad F) = x, we get the SUMMED X coord
		   Thus F must be [x^2/2 0 0], so dF/dx = [x 0 0]

		3. (grad F) = [x y z] -> F = 0.5 * [x^2 y^2 z^2]
		   Then, divide through by the volume at the end to get the average (mu) instead
		   of the sum [x y z]

		4. (grad F) = [X-mu][X-mu]^t -> F = [X-mu][X-mu]^t[X-mu]
		   Then, divide through by the volume * 3 (?) at the end to get the COVARIANCE
		   COV = 1/N Sum( [X-mu][X-mu]^t )

		F gets applied to a POINT, which is taken to be the midpoint of the triangle tile.

		Area of a triangle is the (unnormalized CROSS PRODUCT of two edges)/2.
	*/

	Xferlist * xferList;
    ThallCode::Pointlist_server2 pList;
    int level = 4; // A high subdivision level allows crude approx.
    int nTiles;
	Bpoint * bPointList;
    int which_tile,index;

    Vector3D p1, p2, p3, midpoint, norm, testNorm;

    // Produce boundary tiles
    xferList = convertM3DtoXfer(figure);
    pList.init(xferList);
	delete [] (xferList->atomlist);
	delete xferList;
    pList.ComputeSubdivBoundaryTiles(level);
    // List of quad-tiles, each stored as nTiles sets of 4 sequential Bpoints.
    pList.subdivtileinfo(&nTiles, &bPointList);

	// Initialize counters
    center.set(0.0, 0.0, 0.0);
	volume = 0.0;

    // Traverse and compute the volume and CoG...
    for(which_tile = 0; which_tile < nTiles; which_tile ++)
    {
        index = which_tile * 4;

        // First triangle of tile
        p1.set(bPointList[index].pnt);
        p2.set(bPointList[index + 1].pnt);
        p3.set(bPointList[index + 2].pnt);

        // Norm here is the triangle normal times the area (not unit)
        norm = 0.5 * (p2 - p1).cross(p3 - p1);
        midpoint = (p1 + p2 + p3) / 3.0;

		// For COG, 
		// F that we're summing over is X^2/2 times Area
		// X = midpoint, Area = norm
        center += 0.5 * midpoint.vprod(midpoint.vprod(norm));

		// For Volume
		// F that we're summing over is X times Area
		// X = midpoint, Area = norm
        volume += (midpoint/3)*norm;

        // Second triangle of tile (yes, I am repeating code)
        p2.set(bPointList[index + 2].pnt);
        p3.set(bPointList[index + 3].pnt);

        // Norm here is the triangle normal times the area (not unit)
        norm = 0.5 * (p2 - p1).cross(p3 - p1);
        midpoint = (p1 + p2 + p3) / 3.0;
        center += 0.5 * midpoint.vprod(midpoint.vprod(norm));  // [x^2/2 y^2/2 z^2/2] * area
        volume += (midpoint/3)*norm;
    }
	// Divide center through by volume because we want AVERAGE x,y,z
	// not SUM x,y,z
	center /= volume;


    // Traverse again to compute the covariance...
    Vector3D diff;
	Matrix covariance(3,3);
	Vector tempVector(3);
	
    for(which_tile = 0; which_tile < nTiles; which_tile++)
    {
        index = which_tile * 4;

        // First triangle of tile
        p1.set(bPointList[index].pnt);
        p2.set(bPointList[index + 1].pnt);
        p3.set(bPointList[index + 2].pnt);

		// Norm here is the triangle normal times the area (not unit)
        norm = 0.5 * (p2 - p1).cross(p3 - p1);
        midpoint = (p1 + p2 + p3) / 3.0;

		diff = midpoint - center;  // Diff is [X-mu]

        // Need to convert to Paul vectors for matrix computation
        tempVector(0) = diff.getX();
        tempVector(1) = diff.getY();
        tempVector(2) = diff.getZ();

		//               [X-mu]   *    [X-mu]^trnsp *  [X-mu] * dA
        covariance += (tempVector * tempVector.t()) * (diff   * norm);

        // Second triangle of tile (again, repeating code)
        p1.set(bPointList[index].pnt);
        p2.set(bPointList[index + 1].pnt);
        p3.set(bPointList[index + 2].pnt);

		// Norm here is the triangle normal times the area (not unit)
        norm = 0.5 * (p2 - p1).cross(p3 - p1);
        midpoint = (p1 + p2 + p3) / 3.0;

		diff = midpoint - center;

        // Need to convert to Paul vectors for matrix computation
        tempVector(0) = diff.getX();
        tempVector(1) = diff.getY();
        tempVector(2) = diff.getZ();

		//               [X-mu]   *    [X-mu]^trnsp *  [X-mu] * dA
        covariance += (tempVector * tempVector.t()) * (diff * norm);
    }
	// We don't need to do this uniform scale for rotation only,
	// but the scaling should make the figure and voxel covariances come out the same.

	// Old Tom suggests that this number should be 5, new Tom suggests that this number should be 3(?!)
	//covariance /= (5.0 * volume);

/*	cout << "Figure Covariance: " << endl;
	covariance.print();
	cout << "Det: " << covariance.det() << endl;
*/
    covarianceToQuat(covariance, rotation);

#ifdef BINARY_PABLO_REPORT
	cout << "-------------------------- FIGURE MOMENTS ---" << endl;
	cout << "Centroid: ";
	center.print();
	cout << "Volume: ";
	cout << volume << endl;
	cout << "Rotation: ";
	rotation.print();
	cout << endl;
#endif

	// Everything is still in unit cube space now

}


void MethodOfMoments::imageMoments(Image3D * image, double & volume, Vector3D & center, Quat & rotation)
{
	/*  Determines 'Moments' for segmented binary image objects  
		Basically just average voxel position = mean and sum of "on's" = volume

		This was largely swiped/adapted from Eric's DamFit Code
		Doesn't work with bitplanes anymore, but trivial to add back in
	*/

#ifndef BINARY
    GreyValue * p = image->getVoxels();   // This is the first intensity
#endif
    long xsum=0, ysum=0, zsum=0, voxelCount=0;
	Vector3D voxelCenter, scale;

	int x, y, z;
    for (z = 0; z < image->getZDim(); z++) {
        for (y = 0; y < image->getYDim(); y++) {
#ifdef BINARY
			for (x = 0; x < image->getXDim(); x++) {
				if (image->getVoxelValue(x,y,z) > 0) {
#else
            for (x = 0; x < image->getXDim(); x++, p++) {
                if (*p > 0) {
#endif
					xsum += x;
                    ysum += y;
                    zsum += z;
                    voxelCount ++;
                }
            }
        }
    }

    // express volume in model coordinate space.
	// volume should be positive, even if spacing is negative.
	volume = voxelCount * image->getAbsXSpacing() * image->getAbsYSpacing() * image->getAbsZSpacing();
    // put volume into world spacing (3 dimensions)
	volume /= pow(image->getModelToWorldScale(), 3);	// AGG: This introduces a  half-voxel error 

	// find voxel center
	voxelCenter.set(xsum, ysum, zsum);
	voxelCenter /= voxelCount;
	voxelCenter += .5;	// offset from the voxel's corner to the voxel's center


    // express center in model coordinate space.
	scale.set(image->getAbsXSpacing(), image->getAbsYSpacing(), image->getAbsZSpacing());
	center = voxelCenter.vprod(scale);
	// put center into world spacing
	center /= image->getModelToWorldScale();

    // Covariance calculation
	Matrix covariance(3,3);
	double xDiff, yDiff, zDiff;
	covariance.setAll(0.0);

#ifndef BINARY
    p = image->getVoxels();  // Start at first voxel again
#endif
    for (z = 0; z < image->getZDim(); z++) {
        for (y = 0; y < image->getYDim(); y++) {
#ifdef BINARY
            for (x = 0; x < image->getXDim(); x++) {
	      if (image->getVoxelValue(x,y,z) > 0) {
#else
            for (x = 0; x < image->getXDim(); x++, p++) {
                if (*p > 0) {
#endif
                    xDiff = x - voxelCenter.getX();
                    yDiff = y - voxelCenter.getY();
                    zDiff = z - voxelCenter.getZ();

                    covariance(0,0) += xDiff * xDiff;
                    covariance(0,1) += xDiff * yDiff;
                    covariance(0,2) += xDiff * zDiff;
                    covariance(1,1) += yDiff * yDiff;
                    covariance(1,2) += yDiff * zDiff;
                    covariance(2,2) += zDiff * zDiff;
				}
            }
        }
    }


    // True covariance matrix needs to be scaled, but we don't care about
    // uniform scaling, just rotations. We still need any anisotropic scaling
    covariance(0, 0) *= scale.getX() * scale.getX();
    covariance(0, 1) *= scale.getX() * scale.getY();
    covariance(0, 2) *= scale.getX() * scale.getZ();
    covariance(1, 1) *= scale.getY() * scale.getY();
    covariance(1, 2) *= scale.getY() * scale.getZ();
    covariance(2, 2) *= scale.getZ() * scale.getZ();
    covariance(1, 0) = covariance(0, 1);
    covariance(2, 0) = covariance(0, 2);
    covariance(2, 1) = covariance(1, 2);

	covariance /= voxelCount;  // This should be right for scale 1/N

/*	cout << "Image Covariance: " << endl;
	covariance.print();
	cout << "Det: " << covariance.det() << endl;
*/
	covarianceToQuat(covariance, rotation);

#ifdef BINARY_PABLO_REPORT
	cout << "-------------------------- VOXEL MOMENTS ---" << endl;
	cout << "Centroid: ";
	center.print();
	cout << "Volume: ";
    cout << volume << endl;
	cout << "Rotation: ";
	rotation.print();
	cout << endl;
#endif

}


void MethodOfMoments::covarianceToQuat(Matrix covariance, Quat & quat)
// Again, swiped from Tom's DamFit code
{

	// Eigenvalue decomposition
    Matrix rotMat(3, 3);
    Vector diag(3);
	int ret = covariance.factorEV(diag, rotMat, SPD);
    if(ret != 0)
        cout << "Error: factorEV failed" << endl;

	if(rotMat.det() < 0)
        rotMat = -rotMat;

    // Convert rotation matrix into a quaternion (transposed because we need column-major)
    double tempMat[4][4];
    for(int x = 0; x < 3; x++){
       for(int y = 0; y < 3; y++)
           tempMat[x][y] = rotMat(y, x);
    }
	quat.matrixToQuat(tempMat);
}

void MethodOfMoments::leastAngleRotation(Quat rotation1, Quat rotation2, Quat & returnQuat)
// Again, swiped from Tom's DamFit code
{

	Quat rotation, testRotation, bestRotation;
	Quat axis1Rotation, axis2Rotation, axis3Rotation;

	rotation = rotation1 * rotation2.conj();

    // Find best (smallest) rotation modulo 180 degree rotations about major axes
    // Note: If angle of rotation, theta, is within [-pi, pi], then quaternion
    // with largest w = cos(theta/2) will be smallest. To make sure we are within
    // [-pi, pi], we make sure that we choose the quaternion to have w > 0.
    if(rotation.getW() < 0)
        rotation = -rotation;
    double testW = rotation.getW();
    bestRotation = rotation;

	// Go backwards 2 steps to get the rotation matrix back out of one of the quats
	// to try different flips
	Matrix rotMat(3,3);
	double tempMat[4][4];
	rotation1.buildRotMatrix(*tempMat);
	for (int i=0;i<3;i++){
		for (int j=0;j<3;j++){
			rotMat(i,j) = tempMat[i][j];
		}
	}
	if(rotMat.det() < 0)
        rotMat = -rotMat;

	//rotMat.print();

	// Recover 180 degree rotations
	// Isn't there an easier way to get the quaternions for these?
    axis1Rotation.setAxisAngle(Vector3D(rotMat(0, 0), rotMat(1, 0), rotMat(2, 0)), R_PI);
    axis2Rotation.setAxisAngle(Vector3D(rotMat(0, 1), rotMat(1, 1), rotMat(2, 1)), R_PI);
    axis3Rotation.setAxisAngle(Vector3D(rotMat(0, 2), rotMat(1, 2), rotMat(2, 2)), R_PI);

    testRotation = axis1Rotation * rotation;
    if(testRotation.getW() < 0)
       testRotation = -testRotation;
    if(testRotation.getW() > testW)
    {
#ifdef DEBUG
		cout << "Changed rotation in axis1" << endl;
		cout << "Old:" << endl;
		rotation.print();
		cout << "New:" << endl;
		testRotation.print();
#endif
		bestRotation = testRotation;
		testW = bestRotation.getW();
    }

    testRotation = axis2Rotation * rotation;
    if(testRotation.getW() < 0)
       testRotation = -testRotation;
    if(testRotation.getW() > testW)
    {
#ifdef DEBUG
		cout << "Changed rotation in axis2" << endl;
		cout << "Old:" << endl;
		rotation.print();
		cout << "New:" << endl;
		testRotation.print();
#endif
		bestRotation = testRotation;
		testW = bestRotation.getW();
    }

    testRotation = axis3Rotation * rotation;
    if(testRotation.getW() < 0)
       testRotation = -testRotation;
    if(testRotation.getW() > testW)
    {
#ifdef DEBUG
		cout << "Changed rotation in axis3" << endl;
		cout << "Old:" << endl;
		rotation.print();
		cout << "New:" << endl;
		testRotation.print();
#endif
		bestRotation = testRotation;
    }

	returnQuat = Quat(bestRotation);
}

#ifdef BINARY
void MethodOfMoments::positionObjectOnImage(M3DObject * object, Image3D * image, bool verbose)
{

	// We assume the optimizer selects all appropriate atoms before we start

	Vector3D figureCenter, imageCenter, initialTranslation;
	double figureVolume, imageVolume, initialScale;
	Quat figureRotation, imageRotation, initialRotation;

	MethodOfMoments::objectMoments(object, figureVolume, figureCenter, figureRotation);
	MethodOfMoments::imageMoments(image, imageVolume, imageCenter, imageRotation);

	// Translation is just the diff of the two center's
	initialTranslation = imageCenter - figureCenter;

	// Note here that the SCALE is the CUBE ROOT of the RATIO of VOLUMES
	initialScale = pow(imageVolume/figureVolume,1./3.);

	// Determine the smallest angle rotation we can make to align the figure and image cov's
	MethodOfMoments::leastAngleRotation(imageRotation, figureRotation, initialRotation);

	if (verbose)
	{
		cout << "-------------------------- INITIALIZATION PARAMS ---" << endl;
		cout << "Translation: ";
		initialTranslation.print();
		cout << "Scale: ";
		cout << initialScale << endl;
		cout << "Rotation: ";
		initialRotation.print();
	}
	M3DObject * ref_object = object->clone();

	// Setup the similarity xform & apply it
    SimilarityTransform3D initialTransform;
    initialTransform.setTranslation(initialTranslation);
    initialTransform.setScale(initialScale);
    initialTransform.setRotation(initialRotation);
    object->applySimilarityAboutCOG(&initialTransform);
}
#endif

