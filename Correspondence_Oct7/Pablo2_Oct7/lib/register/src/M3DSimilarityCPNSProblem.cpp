#include <math.h>
#include <iostream>
#include "M3DSimilarityCPNSProblem.h"
#include "Tuning.h"

#ifdef USE_KNOWN_QUADRATIC

double DEBUG_evaluate(const Vector & x);

#endif

extern int globalVerbosity;

using namespace std;

M3DSimilarityCPNSProblem::M3DSimilarityCPNSProblem() : M3DMainFigureProblem()
{
#ifdef BINARY
	similarityObject = NULL;
#endif

	cout << "Problem" << endl;
}

M3DSimilarityCPNSProblem::~M3DSimilarityCPNSProblem()
{
#ifdef BINARY
	delete similarityObject;
#endif
}

void M3DSimilarityCPNSProblem::initialize(Match *m, M3DObject * referenceObject, int treeIndex)
{
	M3DMainFigureProblem::initialize(m, referenceObject, treeIndex);
#ifdef BINARY
	initializeCPNS();
	
	if (similarityObject != NULL) 
	{
		delete similarityObject;
		similarityObject = NULL;
	}
#endif

}

void M3DSimilarityCPNSProblem::setCPNS(M3DCPNSStats * cpnsStat)
{
	cout << "setCPNS" << endl;
	if (cpnsStat == NULL)
		return;
	else
		cpns = cpnsStat;
}

int M3DSimilarityCPNSProblem::parameterCount() const
{
	int maxModes, numModes;

	if (refObject == NULL)
		return 0;

	M3DFigureTreeNode* treeNode = refObject->getFigureTreeRoot(tree_index);
	if (treeNode == NULL)
		return 0;

	int figureId = treeNode->getFigureId();

	maxModes = (int) tuningWt(FigureMaxPGAModes);
	numModes = cpns->getNEigenModes();

	if (maxModes < numModes)
		numModes = maxModes;

	if (0 == (int) tuningWt(FigureMahalanobisOnly))
		numModes += 7;

	return numModes;

}

double M3DSimilarityCPNSProblem::evaluate(const Vector & x)
{
	double geomPenalty;
	double totalPenalty;

	M3DObject * targetObject;
	int figureId;

	evaluationCost += 1;

	if (match == NULL)
		return MAIN_FIGURE_MAX_RETURN;
	
	if (globalVerbosity > 1) 
	{
		cout << "  x:";
		for (int i = 0; i < x.size(); i++)
			cout << ' ' << x(i);
		cout << '\n';
	}
	
	targetObject = createTargetObject(x, figureId);
	std::vector<double> cpns = cpnsCoefficients(x);

#ifdef BINARY
	totalPenalty = match->computeFigureMatchAndPenalties(targetObject, figureId, similarityObject, cpns, geomPenalty);
	lastPenalty = totalPenalty;
#else
	totalPenalty = match->computeFigureMatchAndPenalties(targetObject, figureId, NULL, cpns, geomPenalty);
	lastPenalty = geomPenalty;
#endif

	delete targetObject;

	return totalPenalty;

}

Components M3DSimilarityCPNSProblem::evaluateTerms(const Vector &x)
{
	double geomPenalty;
	Components totalPenalty;
	M3DObject * targetObject;
	int figureId;

	evaluationCost += 1;

	if (match == NULL)
		return Components(MAIN_FIGURE_MAX_RETURN);

	if (globalVerbosity > 1)
	{
		cout << "  x:";
		for (int i = 0; i < x.size(); i++)
			cout << ' ' << x(i);
		cout << '\n';
	}
	

	targetObject = createTargetObject(x, figureId);
	std::vector<double> coeffs = cpnsCoefficients(x);

#ifdef BINARY
	totalPenalty = match->computeFigureMatchAndPenaltiesAsComponents(targetObject, figureId, similarityObject, coeffs, geomPenalty);
	lastPenalty = totalPenalty.rough + totalPenalty.smooth;
#else
	totalPenalty = match->computeFigureMatchAndPenaltiesAsComponents(targetObject, figureId, NULL, cpns, geomPenalty);
	lastPenalty = geomPenalty;
#endif

	delete targetObject;
	return totalPenalty;
}

M3DObject * M3DSimilarityCPNSProblem::createTargetObject(const Vector &x, int & figureId, bool predict)
{
	int primitiveCount;
	M3DObject *targetObject;
	M3DFigureTreeNode *treeNode;
	M3DFigure *figure;

	if (refObject == NULL)
		return NULL;

	treeNode = refObject->getFigureTreeRoot(tree_index);
	if (treeNode == NULL)
		return NULL;

#ifdef BINARY
	figureId = (int) tuningWt(BpFigureId);
#else
	figureId = treeNode->getFigureId();
#endif

	targetObject = refObject->assign();

	figure = targetObject->getFigurePtr(figureId);
	if (figure == NULL)
	{
		delete targetObject;
		return NULL;
	}
	figure->select();

	primitiveCount = figure->getPrimitiveCount();

	vector<double> vals = cpnsCoefficients(x);
	double * coeffs = new double[vals.size()];

	for(int i = 0; i < vals.size(); i++)
		coeffs[i]=vals[i];

	targetObject = cpns->eigenmodeDeformMean(coeffs, refObject);

#ifdef UNFLIPPED
	Vector3D center(0.0, 1.0, 0.0);
#else
	Vector3D center(0.0, 0.0, 0.0);
#endif
#ifdef BINARY
	int useMOM = (tuningWt(BpDoMethodOfMoments) && 
		! tuningWt(BpFigureResetModelBefore));	
    SimilarityTransform3D * xform = useMOM ? &GLOBALMOMTransform 
		: targetObject->getTransformation();
	figure->applySimilarity(*xform, center);

	similarityObject = targetObject->assign();
#endif

	if ((int) tuningWt(FigureMahalanobisOnly) == 0)
	{
		Quat finalRot;
		Vector3D axis;
		axis.set(x(3), x(4), x(5));
		double angle = tuningWt(FigureRotationFactor) * axis.normalize();
		finalRot.setAxisAngle(axis, angle);

		// Apply similarity transform.
		SimilarityTransform3D trf;
		trf.setTranslation( Vector3D(x(0), x(1), x(2)) 
			* tuningWt(FigureTranslationFactor) );
		trf.setRotation( finalRot );
		trf.setScale( exp(tuningWt(FigureScaleFactor) * x(6)) );

		figure->applySimilarityAboutCOG( &trf );
	}

	return targetObject;


	
}

std::vector<double> M3DSimilarityCPNSProblem::cpnsCoefficients( const Vector &x )
{
	vector<double> coeffs;

	int start;

	if (0 != (int) tuningWt(FigureMahalanobisOnly))
		start = 0;
	else
		start = 7;
	for (int i = start; i < x.size(); i++)
	{
		coeffs.push_back(x(i)*tuningWt(FigurePGAFactor));
	}

	return coeffs;
}


void M3DSimilarityCPNSProblem::initializeCPNS()
{
	return;
}

M3DSimilarityCPNSProblemBoundsFunction::M3DSimilarityCPNSProblemBoundsFunction()
{
	if ((int) tuningWt(FigureMahalanobisOnly) == 0) {
		bounds	= Vector(8,
			tuningWt(FigureTranslationBound), tuningWt(FigureTranslationBound), tuningWt(FigureTranslationBound),
			tuningWt(FigureRotationBound)*R_PI/180.0, tuningWt(FigureRotationBound)*R_PI/180.0, tuningWt(FigureRotationBound)*R_PI/180.0,
			log(tuningWt(FigureScaleBound)),
			tuningWt(FigurePGABound));
		factors	= Vector(8,
			tuningWt(FigureTranslationFactor),
			tuningWt(FigureTranslationFactor),
			tuningWt(FigureTranslationFactor),
			tuningWt(FigureRotationFactor),
			tuningWt(FigureRotationFactor),
			tuningWt(FigureRotationFactor),
			tuningWt(FigureScaleFactor),
			tuningWt(FigurePGAFactor));
	}
	else {
		bounds	= Vector(1, 3.0);
		factors	= Vector(1,
			tuningWt(FigurePGAFactor));
	}
}
