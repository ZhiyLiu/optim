#ifndef M3D_FIGUREE_PREDICTOR_H
#define M3D_FIGUREE_PREDICTOR_H

#include "M3DFigure.h"

struct SimilarityTransform
{
	Vector3D	COR;			// center of rotation, defaulted as the COG
	Matrix		cov;			// covariance matrix

	Matrix		rotM;			// rotation matrix
	Quat		rotQ;			// rotation quaternion
	double		scale;			// scale
	Vector3D	trans;			// translation
	Vector3D	transTemp;
};

class M3DFigurePredictor
{
public:
	// Either the newFigure or the newPrimitives hold the predicting "link"
	// atoms specified by linkCount and linkIds[].
	// Note: the predicted figure will be stored back in *figure, instead
	// of *newFigure
    M3DFigurePredictor(int linkCount, int *linkIds, M3DFigure *figure,
		M3DFigure *newFigure);
    M3DFigurePredictor(int linkCount, int *linkIds, M3DFigure *figure,
		M3DPrimitive **newPrimitives);

	// Predict the arbitrarily cut end of the bones by the rest of the atoms.
	// Only the "link" atoms specified by realLinkCoiunt/realLinkIds are 
	// used inthe estimation of the similarity transformation.
    M3DFigurePredictor(int linkCount, int *linkIds, int realLinkCount,
		int *realLinkIds, M3DFigure *figure, M3DFigure *newFigure);

    ~M3DFigurePredictor() {}

	const SimilarityTransform * bestTransform() { return &bestSim; }

private:

	// Estimate the best similarity transformation given two sets of points
	// []x and []y.
	void estimateSimilarityTransformation(int n, Vector3D *x, Vector3D *y, 
		SimilarityTransform *sTrans, Vector3D *COG=NULL);

	// Main prediction function
	// linkCount/linkIds	the # and indices of the "link" atoms
	// figure				the entire figure, with the atoms before transformation
	// newLinkAtoms			all the "link" atoms after the transformation
	// the predicted figure will be stored back to *figure
	void predictFigureByLinkAtoms(int linkCount, int *linkIds, M3DFigure *figure,
		M3DPrimitive **newLinkAtoms);

	// Prediction function for the arbitrarily cut end of the bones
	void predictFigureByLinkAtoms(int linkCount, int *linkIds, int realLinkCount,
		int *realLinkIds, M3DFigure *figure, M3DPrimitive **newLinkAtoms);

	SimilarityTransform bestSim;
};

#endif

