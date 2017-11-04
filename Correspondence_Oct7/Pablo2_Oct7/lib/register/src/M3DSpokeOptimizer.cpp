#include <math.h>
#include "OptimizerBase.h"
#include "M3DSpokeOptimizer.h"
#include "M3DSpokeProblem.h"
#include "M3DPGAPrimitiveStats.h"
#include "GeodesicSym.h"
#include "Tuning.h"
#include "BrentLinearMethod.h"

//#define COMPUTE_PROFILES

//#define DEBUG_SAVE_INTERMEDIATE_MODELS
#ifdef DEBUG_SAVE_INTERMEDIATE_MODELS
//#include "P3DControl.h"
#endif


/*  Because of the skipping of atoms when the optimization gets close to the
    best fit, it is probably best to refine in random order.

    Note that the optimization order is consistent between runs when the debug
    version of the program is built, because srand() is called with a constant
    argument in pablo.cpp.
*/
#ifndef BINARY
#define RANDOMIZE		/* Comment out for serial optimization by atom number */
#endif

const double MIN_PENALTY_WEIGHT = -1;
const double MIN_PENALTY_WEIGHT_EXP = 0.0;
const double MAX_PENALTY_WEIGHT_EXP = 10.0;
const double PENALTY_WEIGHT_SCALE_FACTOR = 10.0;
const double CONSTRAINT_PENALTY_SCALE_FACTOR = 0.02;////1, 10, 100?

#ifndef BINARY
const int SKIPPED_PASSES = 1;
// Defined elsewhere in M3DDeformationOptimizer
//const double DEFAULT_THRESHOLD = 1.0e-3;
//const int DEFAULT_NUM_CONJGRAD_ITERATIONS = 2;
//const int DEFAULT_SCHEDULE_MULTIPLIER = 2;
#endif

using namespace std;


extern int globalVerbosity;


M3DSpokeOptimizer::M3DSpokeOptimizer(int appID)
{
	type = Spoke;
	method = NULL;	// Note: method is a local variable in all cases below

    referenceObject = NULL;
    candidateObject = NULL;

	lastPenalty = 0.0;

	// Specify if the optimizer is used by Pablo (zero) or VSkelTool (nonzero).
	// This value will be passed on to the M3DSpokeProblem instance.
	applicationID = appID;

    counter = NULL;
    skip = NULL;
	lastResult = NULL;
	lastResultSize = 0;
#ifdef BINARY
	if (tuningWt(BpSpokeRandomOrder) == 1)
		perm = NULL;
#else
#ifdef RANDOMIZE
    perm = NULL;
#endif
#endif
    nextAtom = 0;
	lastFigureId = -1;

#ifndef BINARY
    numConjGradientIterations = DEFAULT_NUM_CONJGRAD_ITERATIONS;
    scheduleMultiplier = DEFAULT_SCHEDULE_MULTIPLIER;
    threshold = DEFAULT_THRESHOLD;
#endif

    totalEvaluationCost = 0;
}

#ifndef BINARY

void M3DSpokeOptimizer::setPenaltyWeightsInMatch(const int * penaltyNames,
	const double * penaltyWeights, int numPenalties)
{
	double wt;

	if (match == NULL) {
		cout << "Unable to transfer penalties before Match exists\n";
		return;
	}
	if (numPenalties <= 0)
		return;

	for (int i = 0; i < numPenalties; i++) {
		switch (penaltyNames[i]) {
			case GpSpokeHistogramMatch:
				wt = penaltyWeights[i];
				break;
			case SpokeMahalanobisMatch:
				wt = penaltyWeights[i];
				break;
			case GpSpokeGeometryWeight:
				wt = PENALTY_WEIGHT_SCALE_FACTOR*penaltyWeights[i];
				break;
			case GpSpokeNeighborWeight:
				wt = penaltyWeights[i];
				break;
			case GpSpokeInterWeight:
				wt = CONSTRAINT_PENALTY_SCALE_FACTOR*penaltyWeights[i];
				break;

			case GpSpokeProfileMatch:
				// Never set; see Match.cpp and gpTuning.cpp
			default:
				cout << "Unknown penalty weight: " << tuningShortName(penaltyNames[i]) << endl;
				return;
		}
		match->setWeight(type, penaltyNames[i], wt);
	}
}

#endif

M3DSpokeOptimizer::~M3DSpokeOptimizer()
{
    delete_internal();
}

void  M3DSpokeOptimizer::delete_internal()
{
	int i;
    if(referenceObject != NULL)
        delete referenceObject;
    referenceObject = NULL;

    if(candidateObject != NULL)
        delete candidateObject;
    candidateObject = NULL;

    if(lastBestX)
    {
        for(i = 0; i < lastBestXSize; i++) 
        {
            if (lastBestX[i] != NULL)
                delete lastBestX[i];
        }
        delete [] lastBestX;
        lastBestX = NULL;
    }

    if(counter != NULL)
    {
        delete [] counter;
        counter = NULL;
    }
    if(skip != NULL)
    {
        delete [] skip;
        skip = NULL;
    }
    if(lastResult != NULL)
    {
        delete [] lastResult;
        lastResult = NULL;
		lastResultSize = 0;
    }

#ifdef BINARY
	if (tuningWt(BpSpokeRandomOrder) == 1) {
		if(perm != NULL) {
			delete [] perm;
			perm = NULL;
		}
	}
#else
#ifdef RANDOMIZE
    if(perm != NULL)
    {
        delete [] perm;
        perm = NULL;
    }
#endif
#endif
}

bool M3DSpokeOptimizer::initialize(Match * _match, M3DObject * _referenceObject)

{
    Image3D * image;

    if(_match == NULL || _referenceObject == NULL)
        return false;

    delete_internal();

    referenceObject = _referenceObject->assign();

    candidateObject = referenceObject->assign();
    match = _match;

    match->resetReferenceObject(referenceObject);

    image = match->getTargetImage();
    if(image == NULL)
        return false;

    lastFigureId = -1;
    nextAtom = 0;

	return true;
}

void M3DSpokeOptimizer::reset(Match * _match, M3DObject * _referenceObject)
{
    if(_match == NULL || _referenceObject == NULL)
        return;

    delete_internal();

    referenceObject = _referenceObject->assign();
    candidateObject = _referenceObject->assign();

    match = _match;

    match->resetReferenceObject(referenceObject);

    lastFigureId = -1;
    nextAtom = 0;
}

// This is needed because of the definition of performIterations in
// OptimizerBase.h.
bool M3DSpokeOptimizer::performIterations(int figureId)
{
	assert(false);
	cout << "Error: wrong version of M3DSpokeOptimizer::performIterations() called\n";
	return false;
}

/* This function optimizes a single atom per call from the selected figure, so the
    display may be updated after each atom transformation.  When called repeatedly,
    successive atoms will be optimized.  The order of atoms will be fixed, unless
    RANDOMIZE is defined at compile time in regular Pablo, or for Binary Pablo, if
	tuningWt(BpSpokeRandomOrder) == 1.  Atoms may be skipped when they appear
    to not need optimization (see optimizeFigure() below).  Once all atoms are
    considered refined, the function returns true.  This function is most efficient
    when one figure is fully optimized, before moving on to another figure.

	The verbosity argument controls the amount of output printed, except with
	Binary Pablo.  When it is true, a dot will be printed for each atom optimized.
	When false, there will be no output.
*/
bool M3DSpokeOptimizer::performIterations(int figureId, bool preview,
	bool verbosity)
{
    int numFigs;

    if(match == NULL || referenceObject == NULL || candidateObject == NULL )
        return false;

    numFigs = candidateObject->getFigureCount();

#ifdef BINARY
	if (tuningWt(BpSpokeRSradPenalty) != 0.0 )
		match->setUpAtomPredictor(figureId);

	if (tuningWt(BpSpokeRandomOrder) == 1)
	{
		int i;
		M3DFigure * currFigure;

		// Set up to create a random permutation
		int numPrims = 0;
		if (perm == NULL) {
			currFigure = candidateObject->getFigurePtr(figureId);
			if(currFigure == NULL)
				return false;
			numPrims = currFigure->getPrimitiveCount();
			perm = new short[numPrims];
			for (i = 0; i < numPrims; i++)
				perm[i] = -1;
		}
	}
#else	/* BINARY */
#ifdef RANDOMIZE
	int i;
    M3DFigure * currFigure;

    // Set up to create a random permutation
	int numPrims = 0;
	if (perm == NULL) {
        currFigure = candidateObject->getFigurePtr(figureId);
        if(currFigure == NULL)
            return false;
        numPrims = currFigure->getPrimitiveCount();
		perm = new short[numPrims];
		for (i = 0; i < numPrims; i++)
			perm[i] = -1;
	}
#endif	/* RANDOMIZE */
#endif	/* BINARY */

#ifdef BINARY
    return optimizeFigure(figureId, preview, (globalVerbosity > -1));
#else
    return optimizeFigure(figureId, preview, verbosity);
#endif
}

inline int M3DSpokeOptimizer::selectRandomAtom(int numPrims)
{
	int skipCounter;
    int i;
    int primId;

    // Create a random permutation
	if (nextAtom == 0) {

		int index, temp;
		int lastAtom;

		lastAtom = perm[numPrims - 1];

		for (i = 0; i < numPrims; i++)
			perm[i] = i;
		for (i = 0; i < numPrims; i++) {
			index = i + (rand() % (numPrims - i));
			temp = perm[i];
			perm[i] = perm[index];
			perm[index] = temp;
		}

		if (lastAtom >= 0) {
			// Make sure the same atom is not repeated between changes of permutation
			if (perm[0] == lastAtom) {
				temp = perm[0];
				perm[0] = perm[numPrims - 1];
				perm[numPrims - 1] = temp;
			}
		}
	}

	// Stop if in previous pass all atoms were either skipped or scheduled to skip
	skipCounter = 0;
	if (nextAtom == 0) {
		for (i = 0; i < numPrims; i++) {
			if (counter[i] < 0)
				skipCounter++;
		}
	}
	if (skipCounter == numPrims) {
		allSkipped++;
	}

	// Select next atom to optimize
	primId = perm[nextAtom];
	nextAtom++;
	if (nextAtom == numPrims)
		nextAtom = 0;

	return primId;
}

inline int M3DSpokeOptimizer::selectAtom(int numPrims)
{
	int skipCounter;
    int i;
    int primId;

	// Stop if in previous pass all atoms were either skipped or scheduled to skip
	skipCounter = 0;
	if (nextAtom == 0) {
		for (i = 0; i < numPrims; i++) {
			if (counter[i] < 0)
				skipCounter++;
		}
	}
	if (skipCounter == numPrims) {
		allSkipped++;
	}

	// Select next atom to optimize
	primId = nextAtom;
	nextAtom++;
	if (nextAtom == numPrims)
		nextAtom = 0;

	return primId;
}

void M3DSpokeOptimizer::readyNextRun()
{
	for (int i = 0; i < lastResultSize; i++) {
		counter[i] = 0;
		skip[i] = 0;
        lastResult[i] = 0.0;
	}

	allSkipped = 0;
}

/*  This function is only called from the loop in performIterations().
    The way it terminates is as follows.  Initially, every atom will be
    optimized at least one cycle because the initial value returned by
    problem.evaluate() is very large.  Thereafter, the amount of change
    relative to the last value is compared to threshold.  Once it is
    below that value, the atom is skipped once and thereafter increasingly
    more often, as long as it continues to be below threshold.  The number
    of atoms skipped in each full pass through all atoms of the figure
    is continuously counted, and once it shows that all atoms have been
    skipped by the same, increasing amount, the allSkipped counter is
    incremented (it is initially 0).  This encourages all atoms to refine
    at about the same rate.  Once all atoms of the figure have been skipped
    SKIPPED_PASSES times, the function returns true, indicating that the
    figure's optimization is finished.

	The amount of output is controlled by the verbosity argument.  If it is
    true, a single dot is printed each time an atom is optimized.  If false,
    nothing will be printed.

    This function will hang if no atoms of the figure are selected.
*/

bool M3DSpokeOptimizer::optimizeFigure(int figId, bool preview, bool verbosity)
{
    int i;
    double cum_change;
    int primId;
	int numPrims;

    M3DFigure * referenceFigure = referenceObject->getFigurePtr(figId);
    M3DFigure * candidateFigure = candidateObject->getFigurePtr(figId);

    if (referenceFigure == NULL || candidateFigure == NULL )
            return false;

    numPrims = candidateFigure->getPrimitiveCount();

    if (! lastBestX)
    {
        // Initialize everything
        lastBestX = new Vector *[numPrims];
		lastBestXSize = numPrims;

        for (primId = 0; primId < numPrims; primId++)
        {
            Vector *init	= new Vector(referenceFigure->getNumberOfSpokes());
            lastBestX[primId] = init;
        }
    }

	if (lastFigureId != figId) {
		if (counter != NULL) {
			delete [] counter;
			delete [] skip;
			delete [] lastResult;
		}

        counter = new short[numPrims];
        skip = new short[numPrims];
        lastResult = new double[numPrims];
		lastResultSize = numPrims;

		for (i = 0; i < numPrims; i++) {
			counter[i] = 0;
			skip[i] = 0;
            lastResult[i] = 0.0;
		}

		allSkipped = 0;
		lastFigureId = figId;
	}

	primId = -1;
	// This block requires at least one atom to be selected
	while (primId < 0) {
		if (preview)	// Don't change next random number to be generated
			primId = selectAtom(numPrims);
		else {
#ifdef BINARY
			if (tuningWt(BpSpokeRandomOrder) == 1.0)
				primId = selectRandomAtom(numPrims);
			else
				primId = selectAtom(numPrims);
#else	/* BINARY */
#ifdef RANDOMIZE
			primId = selectRandomAtom(numPrims);
#else
			primId = selectAtom(numPrims);
#endif
#endif	/* BINARY */
		}

		candidatePrimitive = candidateFigure->getPrimitivePtr(primId);
		if (candidatePrimitive == NULL || ! candidatePrimitive->isSelected()) {
			// Unselected atoms are skipped
			counter[primId] = -1;
			primId = -1;
		}
	}

#ifdef BINARY
	if (allSkipped >= tuningWt(BpSpokeSkippedPasses))
#else
	if (allSkipped >= SKIPPED_PASSES)
#endif
		return true;	// Optimization completed

	if (skip[primId] > 0) {
		// Skip optimization of this atom
		skip[primId]--;
		return false;
	}

    referencePrimitive = referenceFigure->getPrimitivePtr(primId);

    if (referencePrimitive == NULL)
        return false;

	// Optimize a single atom
    if (verbosity)
		cout << '.';

	cum_change = optimizeAtom(figId, primId, candidateFigure, preview);

	double d_abs = fabs(cum_change - lastResult[primId]);
#ifdef BINARY
	if (d_abs < tuningWt(BpSpokeSkipThreshold)) {
		counter[primId]--;
		skip[primId] = -((int) tuningWt(BpSpokeSkipMultiplier))*counter[primId];
#else
	if (d_abs < threshold) {
		counter[primId]--;
		skip[primId] = -scheduleMultiplier * counter[primId];
#endif
	}
    else
        counter[primId] = 0;

	lastResult[primId] = cum_change;


	return false;
}

double M3DSpokeOptimizer::getLastBestVal(bool preview)
{
	if (preview) {
		if (lastFigureId < 0)
			return 0.0;
		M3DFigure * candidateFigure = candidateObject->getFigurePtr(lastFigureId);
		int numPrims = candidateFigure->getPrimitiveCount();
		int numPrimsSelected = candidateFigure->numberSelected();

		lastBestVal = 0.0;
		for (int i = 0; i < numPrims; i++) {
			lastBestVal += lastResult[i];
		}
		lastBestVal /= numPrimsSelected;
	}
	return lastBestVal;
}

double M3DSpokeOptimizer::optimizeAtom(int figId, int atomId,
	M3DFigure * candidateFigure, bool preview)
{
	int spoke;
    double initVal;
    double cum_change	= 0.0;

	M3DTubePrimitive* tCandidatePrimitive	= dynamic_cast<M3DTubePrimitive*>(candidatePrimitive);
	M3DTubePrimitive* tReferencePrimitive	= dynamic_cast<M3DTubePrimitive*>(referencePrimitive);

	const int numberOfSpokes	= tReferencePrimitive->getNumberOfSpokes();

	// Flush the stream
	cout << flush;
	if(globalVerbosity >= 2 ) {
		cout << "M3DSpokeOptimizer::optimizeAtom(atomId = " << atomId << " ) called.\n";
	}

	// Optimize a single atom
	const Vector& start = *(lastBestX[atomId]);
	int vecsize = start.size();

	Vector epsilon(vecsize);
    epsilon.setAll(tuningWt(SpokeEpsilon)); 
	M3DSpokeProblem problem(match, referenceObject, candidateObject,
		figId, atomId, -1, applicationID);

	if (preview) {
		double image_match = problem.evaluateImageMatch(start);
		lastBestVal = image_match;
		lastPenalty = 0.0;
		return image_match;
	}

	for( spoke = 0; spoke != numberOfSpokes; ++spoke ) {
		Vector spokeStart(1);
		spokeStart(0)	= start(spoke);
		problem.spokeId	= spoke;
		/*
		 * As we need to record the optimized dr in the spoke when optimizing other
		 * spokes, make sure that when we get back to optimizing this spoke, we
		 * are not doubly applying the optimum.
		 */
		tCandidatePrimitive->setRN(spoke, tReferencePrimitive->getRN(spoke) );

#ifdef COMPUTE_PROFILES
		// Gather match profiles.
		for(double dr = -0.1; dr < +0.1; dr += 0.001 ) {
			spokeStart(0)		= dr;
			const double val	= problem.evaluate(spokeStart);
			cout << "PROFILE: " << atomId << "\t" << spoke << "\t" << dr << "\t" << val << endl;
		}
#else
		NumericalFunction function(problem, epsilon);
		Vector gradient;
		initVal	= function.computeOneJet( spokeStart, gradient);
		BrentLinearMethod opt(function, spokeStart, gradient);
		opt.setBrentLinearSearchBoundFactor( tuningWt(BrentLinearSearchBoundFactor) );
		M3DSpokeProblemBoundsFunction boundsFunction;
	    opt.setXInitialLimit(boundsFunction.evaluate(gradient));
		opt.setBoundsFunction(&boundsFunction);
		opt.run();

		totalEvaluationCost += problem.getEvaluationCost();

		if (opt.value < initVal) {
			(*lastBestX[atomId])(spoke) = opt.x(0);
			cum_change += fabs(initVal - opt.value);
			if(globalVerbosity >= 2) {
				cout << "OPTIMAL: " << atomId << "\t" << spoke << "\t" << opt.x(0) << "\t" << opt.value << endl;
			}
			lastBestVal = opt.value;
#ifndef BINARY
			lastPenalty = problem.getGeometricPenalty();
#else
			lastPenalty = 0;
#endif
		}
		/*
		 * Restoring this primitive at its optimum is necessary because
		 * problem treats candidate primitive as the target
		 * primitive and is thus a temporary in there.
		 * If, we do not restore it now, then target Object
		 * will not have this atom at its optimal place
		 * when a new atom or spoke is selected.
		 */
		tCandidatePrimitive->setRN(spoke, tReferencePrimitive->getRN(spoke) );
		problem.applyTransform(*candidatePrimitive, (*lastBestX[atomId])(spoke) );
#endif
	}

	/*
	 * As Match loses track of the terms in the optimal objective function,
	 * we just recompute match values at the optimum so that it would
	 * report the optimum values.
	 */
	problem.spokeId	= -1;
	problem.evaluate(start);

#ifdef DEBUG_SAVE_INTERMEDIATE_MODELS
	char tempName[50];
    static int counterAtomIterations=0;
	static P3DControl * gControl	= NULL;
	if( gControl == NULL ) {
		gControl= new P3DControl(0);
	}

	M3DObject * temp = candidateObject->assign();
	gControl->newModel(temp);
	sprintf(tempName,"./spoke_iter%03d_atom%03d.m3d",counterAtomIterations++, atomId);
	gControl->write(tempName);
#endif

	// Flush the stream
	cout << flush;

	return cum_change;
}

