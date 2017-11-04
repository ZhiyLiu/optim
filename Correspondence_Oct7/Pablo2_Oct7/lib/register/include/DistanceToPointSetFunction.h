#ifndef _DISTANCE_TO_POINT_SET_FUNCTION_H
#define _DISTANCE_TO_POINT_SET_FUNCTION_H

// MSVC6: disable compiler warning about debugging symbol names too long
#pragma warning (disable:4786)


#include "optima.h"
#include "M3DObject.h"
#include "Image3D.h"
#include <vector>
#include <map>


class SimTransShapeSpace;


/*  DistanceToPointSetFunction

	This class is used for "contour-based" initialization.  It
	contains a function to compute the distance from an object to
	a set of points (in model coordinates).

	This class is used to compute the contour penalty in Match.
	Once the object is constructed and all necessary parameters
	are set, the computeDistanceSqr() function can be called
	with a target object to compute the penalty.

	This class is also used by class FitUnlabeledPoints to fit a
	model to a set of contours.  This is the so-called
	"three-contour" method of alignment.  In this mode, it will
	compute a similarity transform between the model and set of
	contours, and optionally do a PGA fit of the model to the
	contours.

	Originally this class was part of a stand-alone program checked
	into SVN, .../branches/levy/pipeline/fit_unlabled_pts.
*/


class DistanceToPointSetFunction : public Function
{

public:

	enum SimTransType { Scale, Translation, Both, Full, Fail };
	enum stage { None, SimT, PGA };

    // Constructors

	// This constructor reads the specified text file of boundary
    // points (one x,y,z triple per line).  A subdivision level for
	// sampling the object and the number of the figure to match
	// may be provided.
	//
	// If an image is provided, the points will be converted from
	// world to model coordinates.  Otherwise, they should be in
	// model coordinates as provided.
	//
    // If the points are known to come from contours that are a
    // known number of slices in from the edge of the object,
    // then img can be used to convert slices to model coordinates,
    // and padding is the number of slices to adjust for.
	//
	// With no arguments, this constructor serves as the default
	// constructor.
    DistanceToPointSetFunction(const char * pointsFileName = NULL,
		SimTransType simType = Full, int level = 2, int figureId = 0,
		Image3D * img = NULL, int padding = 0);

	// This constructor is like the one above, except that the points
	// are supplied in a Vector3D array of length numPoints.
	DistanceToPointSetFunction(Vector3D * pointSet, int numPoints,
		SimTransType simType = Full, int level = 2, int figureId = 0,
		Image3D * img = NULL, int padding = 0);

	~DistanceToPointSetFunction();

	// Deferred specification of the point set (contours/surface points).
    bool loadPointSet(const char * pointsFileName);
    bool addPointSet(Vector3D * pointSet, int numPoints);

	// Functions to set other parameters after construction
    void setLevel(int level) { subdivLevel = level; }
    void setSizeWeight(double weight) { sizeWeight = weight; }

    // Used to match the top of the object to the top contour
    // 0.04 is a good value (corresponds to matching 20% regions)
    //
    // A value of 0 could only match a boundary point to a contour
    // with the exact same Z value.  That seems non-generic, and the
    // optimizer probably won't find a gradient.
    //
    // A value of 1 allows any boundary point to match any contour.
    // This has been known to cause some serious scaling errors.
	// If not specified, a default value of 0.04 is used.
    void setClosenessLimit(double limit) {
		sqrtClosenessLimit = sqrt(limit);
	}

	// Specify the model to be used for computing a similarity
	// transform
	void setModel(M3DObject * obj) {
		referenceObject = obj;
		referenceObject->select();
	}

	// Specify how the model will be fit to the point set
	void setStage(stage s) { optimizationStep = s; }

	// Access functions for use after construction
    int getNumPoints() { return points.size(); }
    int parameterCount();
    SimTransShapeSpace * getSpace() { return space; }

    // Compute the square distance from the target object to the
	// contours
    double computeDistanceSqr(M3DObject* targetObject);

	// Use only when doing a PGA fit of a model to a point set
    void setPGAOrder(int order) { pgaOrder = order; }
	void setPGACoefs(const Vector & pg);

    void setPostSimTrans(const Vector & simtrans_v) {
		postSimTransCoefs = simtrans_v;
	}

    M3DObject * createTargetObject(const Vector & v);

protected:

	void adjustPointSet();  // Convert the specified points from world to model coordinates
	void normalizePointSet();
    virtual double evaluate(const Vector & v);

    int subdivLevel;	// For constructing bpoints		
    int figure_id;		// For constructing bpoints
    // Only consider distances to points within this limit		
	double sqrtClosenessLimit;
    int slicePadding;	// Extreme z-values of the points came from this many
						// slices inside from the end of the object

    double sizeWeight; // Weight for extra objective function term to encourage
					// the model to span the same number of axial slices as
					// the contours do (after padding adjustment)

    M3DObject * referenceObject;
    Image3D * image;
    SimTransShapeSpace * space;
	stage optimizationStep;
	Vector postSimTransCoefs;
    int pgaOrder;	// Which pgSet to use in the PGA deformation				
	std::vector<double> pgaCoefs;

    std::vector<Vector3D> points;  
	std::vector<Vector3D> normalizedPoints;  
	double minX, maxX, minY, maxY, minZ, maxZ;
    std::vector<double> zSize;
};



#endif /* _DISTANCE_TO_POINT_SET_FUNCTION_H */

