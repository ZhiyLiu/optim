#include <math.h>
#include <iostream>
#include "M3DSimilarityElongationProblem.h"
#include "M3DFigureElongater.h"
#include "LogManager.h"
#include "Tuning.h"

//#define DEBUG_SAVE_INTERMEDIATE_MODELS
//#include "P3DControl.h"


extern int globalVerbosity;

using namespace std;


#ifdef BINARY
extern double computeLandmarkPenalty(M3DObject * targetObject, int figureId,
	M3DObject * liObject);
#endif


M3DSimilarityElongationProblem::M3DSimilarityElongationProblem()
	: M3DMainFigureProblem()
{
	tubePhiMode	= true;
}


void M3DSimilarityElongationProblem::initialize(Match * m, M3DObject * referenceObject,
	int treeIndex)
{
	M3DMainFigureProblem::initialize(m, referenceObject, treeIndex );
    M3DFigureTreeNode* treeNode = refObject->getFigureTreeRoot(treeIndex);
	int figureId = treeNode->getFigureId();
	if( dynamic_cast<M3DTubeFigure*>(refObject->getFigurePtr(figureId)) != NULL ) {
		tubePhiMode	= true;
		setParameterCount(parameterCount()+1);
	}
	else {
		tubePhiMode	= false;
		// When initialized outside, tubePhi is already taken into account.
	}
}


double M3DSimilarityElongationProblem::evaluate(const Vector & x)
{
	double totalPenalty;
    SimilarityTransform3D transformMatrix;
    M3DObject * targetObject;
	double geomPenalty;
    int figureId;

    evaluationCost += 1;

    if (match == NULL)
        return MAIN_FIGURE_MAX_RETURN;

#ifdef BINARY
	if (globalVerbosity > 0) {
		cout << "   x:";
		for (int i = 0; i < x.size(); i++)
			cout << ' ' << x(i);
		cout << '\n';
	}
#endif	/* BINARY */

    targetObject = createTargetObject(x, figureId);


#ifdef BINARY
	// Object for match to compare against the target object
	totalPenalty = match->computeFigureMatchAndPenalties(targetObject,
		figureId, refObject, std::vector<double>(), geomPenalty);
    lastPenalty = totalPenalty;		// Save for future queries
#else	/* BINARY */
	if (match->getType() == BINARY_IMAGE_MATCH)
		return 0.0;

	totalPenalty = match->computeFigureMatchAndPenalties(targetObject,
		figureId, NULL, std::vector<double>(), geomPenalty);
	lastPenalty = geomPenalty;		// Save for future queries
#endif	/* ! BINARY */

    delete targetObject;

#ifdef OPTIMIZATION_VISUALIZER
	copyResultsToVisualizer(&x, totalPenalty);
#endif

    return totalPenalty;
}

Components M3DSimilarityElongationProblem::evaluateTerms(const Vector & x)
{
	Components totalPenalty;
    SimilarityTransform3D transformMatrix;
    M3DObject * targetObject;
	double geomPenalty;
    int figureId;

    evaluationCost += 1;

    if (match == NULL)
        return Components(MAIN_FIGURE_MAX_RETURN);

#ifdef BINARY
	if (globalVerbosity > 0) {
		cout << "   x:";
		for (int i = 0; i < x.size(); i++)
			cout << ' ' << x(i);
		cout << '\n';
	}
#endif	/* BINARY */

    targetObject = createTargetObject(x, figureId);



#ifdef BINARY
	// Object for match to compare against the target object
	totalPenalty = match->computeFigureMatchAndPenaltiesAsComponents(targetObject,
		figureId, refObject, std::vector<double>(), geomPenalty);
    lastPenalty = totalPenalty.rough + totalPenalty.smooth;		// Save for future queries
#else	/* BINARY */
	if (match->getType() == BINARY_IMAGE_MATCH)
		return Components(0.0);

	totalPenalty = match->computeFigureMatchAndPenaltiesAsComponents(targetObject,
		figureId, NULL, std::vector<double>(), geomPenalty);
	lastPenalty = geomPenalty;		// Save for future queries
#endif	/* ! BINARY */

    delete targetObject;

#ifdef OPTIMIZATION_VISUALIZER
	copyResultsToVisualizer(&x, totalPenalty.rough + totalPenalty.smooth);
#endif

    return totalPenalty;
}


// For Binary Pablo, this returns the reference object, which has been transformed by
// MOM and by the SimTranform (as indicated by x).
// The predict argument is only used in Adaptive Pablo
M3DObject * M3DSimilarityElongationProblem::createTargetObject(const Vector & x,
	int & figureId, bool predict)
{
    Vector3D trans;
    Quat finalRot;
    Vector3D axis;
    double angle;
    double elongation	= 0.0;
    M3DFigureElongater elongater;
    M3DFigure * figure;
    M3DObject * targetObject;
#ifndef BINARY
    M3DFigureTreeNode * treeNode;
#endif

    if (refObject == NULL)
        return NULL;

#ifndef BINARY
    treeNode = refObject->getFigureTreeRoot(tree_index);
    if (treeNode == NULL)
        return NULL;

    figureId = treeNode->getFigureId();
#else
	figureId = (int) tuningWt(BpFigureId); // MultiObject
#endif

    targetObject = refObject->assign();

    figure = targetObject->getFigurePtr(figureId);
    if (figure == NULL) {
        delete targetObject;
        return NULL;
    }
    figure->select();

    // Apply initial transformation
#ifdef UNFLIPPED
	Vector3D center(0.0, 1.0, 0.0);
#else
	Vector3D center(0.0, 0.0, 0.0);
#endif
#ifndef BINARY
	figure->applySimilarity(mainFigProblemSimTransform, center);
#else
	figure->applySimilarity(GLOBALMOMTransform, Vector3D(0.0, 0.0, 0.0));
#endif

    // Apply similarity transform
	SimilarityTransform3D trf;
	trf.setTranslation( Vector3D(x(0), x(1), x(2)) 
		* tuningWt(FigureTranslationFactor) );

    axis.set(x(3), x(4), x(5));
    angle = tuningWt(FigureRotationFactor) * axis.normalize();
    finalRot.setAxisAngle(axis, angle);
	trf.setRotation( finalRot );

	if( tubePhiMode ) {
		trf.setTubePhi( x(6) * tuningWt(FigureTubePhiFactor) );
		if (parameterCount() > 7) {
			trf.setScale( exp(tuningWt(FigureScaleFactor) * x(7)) );
		}
	}
	else {
		if (parameterCount() > 6) {
			trf.setScale( exp(tuningWt(FigureScaleFactor) * x(6)) );
		}
	}

	figure->applySimilarityAboutCOG( &trf );

	if (parameterCount() > 7 && !tubePhiMode) {
		// Elongate the figure
		elongation = exp(tuningWt(FigureElongationFactor) * x(7));
		elongater.initialize(dynamic_cast<M3DQuadFigure*>( figure));
		elongater.elongate(elongation);
	}

	if( globalVerbosity >= 1 ) {
	    static int iterationCounter = 0;
    	std::cout << iterationCounter++ << ": T:";
	    trf.getTranslation().print();
    	cout << " R:";
		trf.getRotation().print();
    	cout << " S: " << trf.getScale() << " Phi: " << trf.getTubePhi() << " E: " << elongation << endl;
	}

#ifdef DEBUG_SAVE_INTERMEDIATE_MODELS
	char tempName[50];
	static int counterAtomIterations=0;
	static P3DControl * gControl	= NULL;
	if( gControl == NULL ) {
		gControl= new P3DControl(0);
	}

	M3DObject * temp = targetObject->assign();
	gControl->newModel(temp);
	sprintf(tempName,"figure_iter%03d.m3d",counterAtomIterations++);
	gControl->write(tempName);
#endif


    return targetObject;
}

double M3DSimilarityElongationProblemBoundsFunction::bound( int i )
{
	assert(i >= 0);
	Vector bounds;
	// tubes don't have support for elongation.
	if( problem->tubePhiMode ) {
		bounds	= Vector(8,
			tuningWt(FigureTranslationBound), tuningWt(FigureTranslationBound), tuningWt(FigureTranslationBound),
			tuningWt(FigureRotationBound)*R_PI/180.0, tuningWt(FigureRotationBound)*R_PI/180.0, tuningWt(FigureRotationBound)*R_PI/180.0,
			tuningWt(FigureTubePhiBound)*R_PI/180.0,
			log(tuningWt(FigureScaleBound)));
	}
	else {
		bounds	= Vector(8,
			tuningWt(FigureTranslationBound), tuningWt(FigureTranslationBound), tuningWt(FigureTranslationBound),
			tuningWt(FigureRotationBound)*R_PI/180.0, tuningWt(FigureRotationBound)*R_PI/180.0, tuningWt(FigureRotationBound)*R_PI/180.0,
			log(tuningWt(FigureScaleBound)),
			log(tuningWt(FigureElongationBound)));
	}
	return bounds( (i < bounds.size()) ? i : (bounds.size() - 1) );
}

double M3DSimilarityElongationProblemBoundsFunction::factor(int i )
{
	assert(i >= 0);
	Vector factors;
	// tubes don't have support for elongation.
	if( problem->tubePhiMode ) {
		factors	= Vector(8,
			tuningWt(FigureTranslationFactor), tuningWt(FigureTranslationFactor), tuningWt(FigureTranslationFactor),
			tuningWt(FigureRotationFactor), tuningWt(FigureRotationFactor), tuningWt(FigureRotationFactor),
    	    tuningWt(FigureTubePhiFactor),
			tuningWt(FigureScaleFactor));
	}
	else {
		factors	= Vector(8,
			tuningWt(FigureTranslationFactor), tuningWt(FigureTranslationFactor), tuningWt(FigureTranslationFactor),
			tuningWt(FigureRotationFactor), tuningWt(FigureRotationFactor), tuningWt(FigureRotationFactor),
			tuningWt(FigureScaleFactor), 
			tuningWt(FigureElongationFactor));
	}
	return factors( (i < factors.size()) ? i : (factors.size() - 1) );
}

