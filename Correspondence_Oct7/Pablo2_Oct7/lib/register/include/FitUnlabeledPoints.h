#ifndef _FIT_UNLABELED_POINTS_H_
#define _FIT_UNLABELED_POINTS_H_


class Image3D;
class M3DObject;
class SimilarityTransform3D;

/*  Class FitUnlabeledPoints

	Fit align a model to a set of surface points, typically comprising
	contours from an anastruct.  The fit may simply produce a similarity
	transform, or may also include PGA deformation.  In the later case,
	the fitted model may be obtained.

    There are three kinds of fits that can be done: similarity transform,
    PGA, or alternating between them both.  The type argument of the
    fit() command is used to specify the kind of fit.

    If order is 0 or more, fitting will alternate
    between optimizing a similarity transform and optimizing the PGA
    coefficients, until a reasonable stopping point is reached.  If
    order is negative, only a similarity transform fit will be done.
    Note that the usual correspondence between figureId and order should
    be observed here.

*/
class FitUnlabeledPoints
{

public:

	// The order in this enum matters: PGA methods must be clumped at the end.
	// See FitUnlabeledPoints.cpp for details.
    enum PointFitType { Scale, Trans, ScaleTrans, Similarity, PGA, ScalePGA, SimPGA };

	// Constructor.  Note that the usual correspondence between figureId
	// and order should be observed here.
    FitUnlabeledPoints(int level = 2, int figureId = 0, int order = -1,
		double closeness = 0.04, double outerLimit = 1.0e-5,
		double innerLimit = 1.0e-6);
    ~FitUnlabeledPoints();

	// Perform a fit.  If type is SimPGA and order is 0 or more,
	// fitting will alternate between optimizing a similarity transform
	// and optimizing the PGA coefficients, until a reasonable stopping
	// point is reached.  If type is PGA, only a PGA fit will be done.
	// If type is SimPGA and order is negative, only a similarity
	// transform fit will be done.  If type is ScalePGA  and order is
	// negative, only a scaling fit will be done.
	//
	// If no image is provided, the points are assumed to be in model
	// coordinates.  If an image is provided, they must be in world
	// coordinates.
    bool fit(M3DObject * model, double * x, double * y, double * z,
		int numPoints, Image3D * image = NULL,
		PointFitType type = SimPGA);

	// This works like the preceding function, except that the point list
	// is provided via a file.
    bool fit(M3DObject * model, const char * filename, Image3D * image = NULL,
		PointFitType type = SimPGA);

    // Functions for changing parameters
    void setCloseness(double closeness) { closenessLimit = closeness; }
    void setModelParameters(int figureId, int level = 2, int order = -1) {
	    figure_id = figureId;
        surf_level = level;
	    pg_order = order;
    }
    void setOuterLimit(double outerLimit) { convergence_limit = outerLimit; }
    void setInnerLimit(double innerLimit) { fit_limit = innerLimit; }

	// Obtain results of the similarity transformation fit
    void result(SimilarityTransform3D & outputXform);

	// Obtain results of the PGA fit.  Note that ownership of the
	// object passes to the caller.  When done with it, it should
	// be deleted.
    M3DObject * result() {
		M3DObject * obj = bestObj;
		bestObj = NULL;
		return obj;
	}
    void result(std::vector<double> & pgaCoefs);

protected:

    double closenessLimit;
    int figure_id;
    int pg_order;
    double convergence_limit;
	double fit_limit;
    int surf_level;
	SimilarityTransform3D * xform;
	M3DObject * bestObj;
	Vector * pga_coefs;

private:

	bool runPointsFit(DistanceToPointSetFunction & f,
		M3DObject * model, Image3D * image, PointFitType type);
};


#endif	/* _FIT_UNLABELED_POINTS_H_ */

