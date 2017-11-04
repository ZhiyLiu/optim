#include <math.h>
#include <iostream>
#include "M3DSimilarityPGAProblem.h"
#include "Tuning.h"

//#define DEBUG_SAVE_INTERMEDIATE_MODELS
//#include "P3DControl.h"

//#define USE_KNOWN_QUADRATIC


#ifdef USE_KNOWN_QUADRATIC

/* This is an aid to debugging an optimization method (ConjugateGradientMethod)
 * when this symbol is defined the objective function will be a known quadratic
 */
double DEBUG_evaluate(const Vector & x);

#endif // USE_KNOWN_QUADRATIC


extern int globalVerbosity;


using namespace std;


M3DSimilarityPGAProblem::M3DSimilarityPGAProblem() : M3DMainFigureProblem()
{
#ifdef BINARY
	similarityObject = NULL;
#endif
}

M3DSimilarityPGAProblem::~M3DSimilarityPGAProblem()
{
#ifdef BINARY
	delete similarityObject;
#endif
}

void M3DSimilarityPGAProblem::initialize(Match * m, M3DObject * referenceObject,
	int treeIndex)
{
    M3DMainFigureProblem::initialize(m, referenceObject, treeIndex);
#ifdef BINARY
    initializePGA();

	if (similarityObject != NULL) {
		delete similarityObject;
		similarityObject = NULL;
	}
#endif
}

void M3DSimilarityPGAProblem::setPGA(M3DPGAStats * pgaStat)
{
	if (pgaStat == NULL)
		return;
	else
		pga = pgaStat;
}

// Return # of modes (up to FMAHAMAX) + count of fixed components
// (3*trans, 3*rot, scale).
//		AGG: This function probably should take the figureId as an argument.
int M3DSimilarityPGAProblem::parameterCount() const
{
	int maxModes, numModes;

    if (refObject == NULL) 
        return 0;

    M3DFigureTreeNode* treeNode = refObject->getFigureTreeRoot(tree_index);
    if (treeNode == NULL)
        return 0;

	// Need to know figure id to get the right PG set
	int figureId = treeNode->getFigureId();
	int order = pga->findPGSetId(figureId);

	PGData * pgPtr = pga->getPGDataPtr(order);

	maxModes = (int) tuningWt(FigureMaxPGAModes);
	numModes = pgPtr->numPGs;

	if (maxModes < numModes)
		numModes = maxModes;

	if (0 == (int) tuningWt(FigureMahalanobisOnly))
		numModes += 7;

	return numModes;
}

double M3DSimilarityPGAProblem::evaluate(const Vector & x)
{

#ifdef USE_KNOWN_QUADRATIC
	return DEBUG_evaluate(x);
#endif

	double geomPenalty;
	double totalPenalty;
    M3DObject * targetObject;
    int figureId;

    evaluationCost += 1;

    if (match == NULL)
        return MAIN_FIGURE_MAX_RETURN;

	if (globalVerbosity > 1) {
		cout << "   x:";
		for (int i = 0; i < x.size(); i++)
			cout << ' ' << x(i);
		cout << '\n';
	}

	// In Binary Pablo, this also computes the similarityObject
    targetObject = createTargetObject(x, figureId);
    std::vector<double> pga	= pgaCoefficients(x);

#ifdef BINARY
	totalPenalty = match->computeFigureMatchAndPenalties(targetObject, figureId,
		similarityObject, pga, geomPenalty);
    lastPenalty = totalPenalty;		// Save for future queries
#else	/* BINARY */
	totalPenalty = match->computeFigureMatchAndPenalties(targetObject, figureId,
		NULL, pga, geomPenalty);
    lastPenalty = geomPenalty;		// Save for future queries
#endif	/* BINARY */

    delete targetObject;

#ifdef OPTIMIZATION_VISUALIZER
	copyResultsToVisualizer(&x, totalPenalty);
#endif

    return totalPenalty;
}

Components M3DSimilarityPGAProblem::evaluateTerms(const Vector & x)
{

	double geomPenalty;
	Components totalPenalty;
    M3DObject * targetObject;
    int figureId;

    evaluationCost += 1;

    if (match == NULL)
        return Components(MAIN_FIGURE_MAX_RETURN);

	if (globalVerbosity > 1) {
		cout << "   x:";
		for (int i = 0; i < x.size(); i++)
			cout << ' ' << x(i);
		cout << '\n';
	}

	// In Binary Pablo, this also computes the similarityObject
    targetObject = createTargetObject(x, figureId);
    std::vector<double> pga	= pgaCoefficients(x);

#ifdef BINARY
	totalPenalty = match->computeFigureMatchAndPenaltiesAsComponents(targetObject, figureId,
		similarityObject, pga, geomPenalty);
    lastPenalty = totalPenalty.rough + totalPenalty.smooth;		// Save for future queries
#else	/* BINARY */
	totalPenalty = match->computeFigureMatchAndPenaltiesAsComponents(targetObject, figureId,
		NULL, pga, geomPenalty);
    lastPenalty = geomPenalty;		// Save for future queries
#endif	/* BINARY */

#ifdef DEBUG_SAVE_INTERMEDIATE_MODELS
	char tempName[50];
	static int counterAtomIterations=0;
	static P3DControl * gControl	= NULL;
	if( gControl == NULL ) {
		gControl= new P3DControl(0);
	}

	M3DObject * temp = targetObject->assign();
	gControl->newModel(temp);
	sprintf(tempName,"./figure_iter%03d.m3d",counterAtomIterations++);
	gControl->write(tempName);
#endif

    delete targetObject;

#ifdef OPTIMIZATION_VISUALIZER
	copyResultsToVisualizer(&x, totalPenalty.rough + totalPenalty.smooth);
#endif

    return totalPenalty;
}


// This function returns a copy of the model after the PGA, initial, and
// similarity transforms are applied.  In Binary Pablo, as a side effect,
// it also sets similarityObject to contain a copy of the model after the
// PGA and initial transforms are applied, but without the similarity
// transform applied.
M3DObject * M3DSimilarityPGAProblem::createTargetObject(const Vector & x,
	int & figureId, bool predict )
{
    int primitiveCount;
    M3DObject * targetObject;
    M3DFigureTreeNode * treeNode;
    M3DFigure * figure;

    if (refObject == NULL || x.size() != parameterCount())
        return NULL;

    treeNode = refObject->getFigureTreeRoot(tree_index);
    if (treeNode == NULL)
        return NULL;

#ifdef BINARY
	figureId = (int) tuningWt(BpFigureId);	// MultiObject
#else
    figureId = treeNode->getFigureId();
#endif

	//JJ: It does not matter what object is assigned
	//    to 'targetObject'.
    targetObject = refObject->assign();

    figure = targetObject->getFigurePtr(figureId);
    if (figure == NULL) {
        delete targetObject;
        return NULL;
    }
    figure->select();

    primitiveCount = figure->getPrimitiveCount();

	vector<double> vals	= pgaCoefficients(x);

	// Apply single figure statistics
	int order = pga->findPGSetId(figureId);
	pga->doPGADeform(targetObject, vals, order, false);

    // Apply initial transformation
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
#else
	figure->applySimilarity(mainFigProblemSimTransform, center);

	// Apply hand placement transform.
	SimilarityTransform3D handTransform = pga->getHandPlacementTransform();
	figure->applySimilarity(handTransform, center);
#endif

	if ((int) tuningWt(FigureMahalanobisOnly) == 0) {
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

std::vector<double> M3DSimilarityPGAProblem::pgaCoefficients( const Vector&x )
{
	vector<double> coeffs;
	int start;
	if (0 != (int) tuningWt(FigureMahalanobisOnly))
		start = 0;
	else
		start = 7;
	for (int i = start; i < x.size(); i++) {
		coeffs.push_back(x(i)*tuningWt(FigurePGAFactor));
	}
	return coeffs;
}


M3DSimilarityPGAProblemBoundsFunction::M3DSimilarityPGAProblemBoundsFunction()
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


#ifdef USE_KNOWN_QUADRATIC

/* This is an aid to debugging an optimization method (ConjugateGradientMethod)
 * when this symbol is defined the objective function will be the following quadratic
 *
 * f(X) = (1/2) x A x' - b x
 *
 * where A =   {55, 62, 50, 41, 34, 
 *              62, 92, 53, 42, 38, 
 *							50, 53, 55, 37, 39, 
 *							41, 42, 37, 37, 23, 
 *		  				34, 38, 39, 23, 31};
 *
 * and B = {0.125, 0.125, 0.125, 0.125, 0.125}
 *
 * When Mathematica optimizes this equation, using the conjugate gradient method
 * the sequence of X vectors is:
 *
   {{0.000564023, 0.000564023, 0.000564023, 0.000564023, 0.000564023}, 
	 {0.000094392, -0.00173169, 0.000419029, 0.00261033, 0.00321902}, 
	 {-0.00573557, -0.000126562, -0.0134987, 0.0130744, 0.017357}, 
	 {-0.00607567, -0.000128164, -0.0139304, 0.0124667, 0.01846}, 
	 {-0.00611456, 0.0000709823, -0.0140622, 0.0121322, 0.0193457}, 
	 {-0.00333643, -0.000744002, -0.0159763, 0.0112274, 0.0203819}, 
	 {-0.00333019, -0.000723961, -0.0159862, 0.0112256, 0.0203771}, 
	 {-0.00333134, -0.000711641, -0.0160025, 0.011219,  0.0203672}, 
	 {-0.00333081, -0.000712013, -0.0160028, 0.0112204, 0.020366}, 
	 {-0.00333083, -0.000712033, -0.0160027, 0.0112204, 0.020366}, 
	 {-0.00333103, -0.000712117, -0.0160012, 0.01122,   0.0203647}, 
	 {-0.00333115, -0.000711863, -0.0160008, 0.0112198, 0.0203641}, 
	 {-0.00333115, -0.000711862, -0.0160008, 0.0112198, 0.0203641}}
 *
 * in order for this to work, the Vector &x must have 5 elements --
 * that is the user should load a PGA file with 5 principal geodesics defined in it
 *
 *
 */
double DEBUG_evaluate(const Vector& x) {

	if (x.size() < 5) { 
		//x needs to be a 5 vector) 
		// not enough data to optimize the quadratic so this will fail
		return 0;
	}

	double A[] = {55, 62, 50, 41, 34, 
								62, 92, 53, 42, 38, 
								50, 53, 55, 37, 39, 
								41, 42, 37, 37, 23, 
								34, 38, 39, 23, 31};


	double prod[5] = {0, 0, 0, 0, 0};
	double result = 0;
	int i;
	for (i=0; i<5; i++) {
		for (int j=0; j<5; j++) {
			prod[i] += x(j) * A[(5 * j) + i];  // x A
		}
		result += 0.5 * (prod[i] * x(i)) - (x(i) * 0.125); // 1/2 x A x' - b' x
	}

	return result;
}

#endif	/* USE_KNOWN_QUADRATIC */

#ifdef BINARY

void M3DSimilarityPGAProblem::initializePGA()
{
    pga->initializeMean(refObject);
}

#endif


