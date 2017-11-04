#include <fstream>
#include <iomanip>
#include <math.h>
#include "OptimizerBase.h"
#include "SurfacePatchEnsemble.h"
#include "M3DMainFigureOptimizer.h"
#include "Tuning.h"
#ifdef OPTIMIZATION_VISUALIZER
#include "LogManager.h"
#endif

#include <sstream>

//#define DEBUG

// Set define below to generate objective function plots
//#define MAIN_FIGURE_OUTPUT_TEST


using namespace std;

#ifdef BINARY
SimilarityTransform3D GLOBALMOMTransform;	// an ugly global to simplify where it's set and used
#endif

M3DMainFigureOptimizer::M3DMainFigureOptimizer(Match * m, int count) : epsilon(0)
{
	type = MainFigure;

    refObject = NULL;
	addToResidueObject	= NULL;
	match = m;

    problem = NULL;
    wrapperProblem = NULL;

	numParameters = count;
    counter = 0;

	lastBestX = new Vector *[1];
	lastBestXSize = 1;
	*lastBestX = new Vector;

#ifdef OPTIMIZATION_VISUALIZER
	if (match) {
		const Match::matchResult * results = match->getFigureStageResults();
		if (results) {
			globalLogManager.attributeNames[0] = "Objective Function";
			for (int i = 0; i < MAX_NUM_FIG_MATCH_RESULTS; i++) {
				globalLogManager.attributeNames[1 + i] = results[i].name;
			} //for loop
			globalLogManager.numAttributes = MAX_NUM_FIG_MATCH_RESULTS+1;
		} // results
	} // match
#endif /* OPTIMIZATION_VISUALIZER */
}

M3DMainFigureOptimizer::~M3DMainFigureOptimizer()
{
    if (refObject != NULL)
        delete refObject;
    if (problem != NULL)
        delete problem;
    if (wrapperProblem != NULL)
        delete wrapperProblem;
    if (boundsFunction != NULL)
        delete boundsFunction;
	// method belongs to OptimizerBase
	// match belongs to P3DControl
#ifdef BINARY
	extern int globalVerbosity;
	if (globalVerbosity > -1) {
		const int oldPrecision	= cout.precision(18);
		cout << "   bestx=";
		for (int i = 0; i < (lastBestX[0])->size(); i++)
			cout << ' ' << setw(20) << (*(lastBestX[0]))(i);
		cout << '\n';
		// FIXME: HACK: only works with M3DSimilarityElongationProblem.cpp
		// and should go in the destructor of the respective base classes anyway.
		const Vector& x = *lastBestX[0];
		if( x.size() == 8 || x.size() == 7 ) {
			Vector3D trans;
			Quat finalRot;
			Vector3D axis;
			double angle;
			double scale;
			double tubePhi;
			trans.set(x(0), x(1), x(2));
			trans	*= tuningWt(FigureTranslationFactor);

			axis.set(x(3), x(4), x(5));
			angle	= tuningWt(FigureRotationFactor) * axis.normalize();
			finalRot.setAxisAngle(axis, angle);
			scale	= exp(tuningWt(FigureScaleFactor) * x(7));
			tubePhi	= x(6) * tuningWt(FigureTubePhiFactor);
			cout << "bestSim: "
				 << setw(20) << trans.getX() << ' '
				 << setw(20) << trans.getY() << ' '
				 << setw(20) << trans.getZ() << ' '
				 << setw(20) << finalRot.getX() << ' '
				 << setw(20) << finalRot.getY() << ' '
				 << setw(20) << finalRot.getZ() << ' '
				 << setw(20) << finalRot.getW() << ' '
				 << setw(20) << scale << ' '
				 << setw(20) << tubePhi << endl;
		}
		cout.flush();
		cout.precision(oldPrecision);
	}
#endif
}

#ifndef BINARY

void M3DMainFigureOptimizer::setPenaltyWeightsInMatch(const int * penaltyNames,
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
			case GpFigureHistogramMatch:
				wt = penaltyWeights[i];
				break;
			case GpFigureMultiPatchMatch:
				wt = penaltyWeights[i];
				break;
			case FigureDQFMatch:
				wt = penaltyWeights[i];
				break;
			case FigureMahalanobisMatch:
				wt = M3DMainFigureProblem::PENALTY_SCALE_FACTOR*penaltyWeights[i];	// AGG: Move these constants to this class?
				break;
			case GpFigureGeometryWeight:
				wt = M3DMainFigureProblem::PENALTY_SCALE_FACTOR*penaltyWeights[i];
				break;
			case GpFigureInterObjectWeight:
				wt = M3DMainFigureProblem::CONSTRAINT_PENALTY_SCALE_FACTOR*penaltyWeights[i];
				break;

			case GpFigureProfileMatch:
				// Never set; see Match.cpp and gpTuning.cpp
			default:
				cout << "Unknown penalty weight: " << tuningShortName(penaltyNames[i]) << endl;
				return;
		}
		match->setWeight(type, penaltyNames[i], wt);
	}
}

#endif

bool M3DMainFigureOptimizer::initialize(M3DObject * referenceObject,
	int treeIndex)
{
    Image3D * image;
	double eps;

#ifdef OPTIMIZATION_VISUALIZER
	globalLogManager.setOptimizer(this);
#endif

#ifdef DEBUG
	cout << "M3DMainFigureOptimizer::initialize()\n";
#endif
    tree_index = treeIndex;

    if (match == NULL || referenceObject == NULL || problem == NULL) {
		cout << "Error: Match, referenceObject, or problem is NULL" << endl;
        return false;
	}

    refObject = referenceObject->assign();
    lastBestObject = refObject->assign();

	match->resetReferenceObject(refObject);

#ifdef BINARY
	if (tuningWt(BpFigureRSradPenalty) != 0.0) {
		match->setUpAtomPredictor(treeIndex);
	}
#endif 

    image = match->getTargetImage();
    if (image == NULL) {
		cout << "Error: Unable to initialize figure optimizer; no image" << endl;
        return false;
	}
	
    problem->initialize(match, refObject, treeIndex);
	problem->setAddToResidueObject( addToResidueObject );
	numParameters = problem->parameterCount();

    Vector start(numParameters);

	//
	// If a set of starting coefficients were specified, then use them.
	//

	const char* coefficientList = tuningStr(FigureStartingCoefficients);

	if( coefficientList == NULL ) {
    	start.setAll(0.0);
	}
	else {
		cout << "Using starting coefficients <";
		istringstream tokens(coefficientList);
		double* dstart	= new double[numParameters];
		int count	= 0;
		for( int i = 0; i != numParameters; ++i ) {
			tokens >> dstart[i];
			if( !tokens ) {
				break;
			}
			start(i) = dstart[i];
			cout << dstart[i] << ",";
			++count;
		}
		cout << "count:" << count << "(needed:" << numParameters << ")>" << endl << flush;
		if( count != numParameters ) {
			cout << "ERROR: Number of parameters do not match" << endl << flush;
			exit(-1);
		}

		//STOUGH hack: I want to do the starting values for multiple single figures in 
		//one run.  Simple solution.  Call tuneVals.setString after I'm done reading 
		//tokens for this figure. I'm not sure why, but this makes a repeat of the last 
		//coeff.  It's irrelevant though, since this function always knows how many params 
		//it's looking for.
		ostringstream remains;

		double temp;
	
		while(tokens)
		{
			tokens >> temp;
			remains << temp << " ";
		}

		(void) tuneVals.setString(FigureStartingCoefficients, remains.str().c_str());

		delete[] dstart;
	}

	epsilon.setSize(numParameters);
	eps = tuningWt(FigureEpsilon);
#ifdef BINARY
	epsilon.setAll(eps);
#else
	if (eps < 0.0)
		eps = getNumericalDerivativeEpsilon();
	epsilon.setAll(eps);
#endif


    if (wrapperProblem != NULL)
        delete wrapperProblem;

    wrapperProblem = new NumericalFunction(*problem, epsilon);

	setOptimizerPosition(start);


    lastBestVal = problem->evaluate(start);
    *(lastBestX[0]) = start;
    lastPenalty = problem->getLastPenalty();
    (dynamic_cast<ConjugateGradientMethod *>(method))->getFunction().setEvaluationCost(0);
	if( lastBestObject ) {
		delete lastBestObject;
	}
	lastBestObject = createTargetObject(*(lastBestX[0]));

	return true;
}

void M3DMainFigureOptimizer::reinitialize()					// AGG: Might match have changed?  If so, need to pass a new one in and reset penalties into it.
{
#ifdef DEBUG
	cout << "M3DMainFigureOptimizer::reinitialize()\n";
#endif
	if(wrapperProblem != NULL)
		delete wrapperProblem;

    epsilon *= 0.5;
    wrapperProblem = new NumericalFunction(*problem, epsilon);

	setOptimizerPosition(*(lastBestX[0]));
}

bool M3DMainFigureOptimizer::performIterations(int nIterations)
{
    int i;

#ifdef DEBUG
    cout << "M3DMainFigureOptimizer::performIterations() called for "
		<< nIterations << " iterations" << endl;
#endif

    if (method == NULL)
        return false;

    i = 0;
    while (i < nIterations && ! isFinished()) {
        ((ConjugateGradientMethod *) method)->performIteration();
        i++;
    }

    if (method->getBestEverValue() < lastBestVal) {
        lastBestVal = ((ConjugateGradientMethod *) method)->getBestEverValue();
        *(lastBestX[0]) = ((ConjugateGradientMethod *) method)->getBestEverX();
        lastPenalty = problem->getLastPenalty();
#ifdef DEBUG
		cout << "lastPenalty = " << lastPenalty << endl;
#endif

        if (lastBestObject != NULL)
            delete lastBestObject;
		lastBestObject = createTargetObject(*(lastBestX[0]));
#ifdef BINARY
		extern int globalVerbosity;
		if (globalVerbosity > -1) {
			const int oldPrecision	= cout.precision(18);
			cout << "   bestx=";
			for (int i = 0; i < (lastBestX[0])->size(); i++)
				cout << ' ' << setw(20) << (*(lastBestX[0]))(i);
			cout << '\n';
			cout.flush();
			cout.precision(oldPrecision);
		}
#endif
    }

#ifdef BINARY
    if (isFinished()) {
        if(counter < 2)
            reinitialize();
        counter++;
    }
#endif
	return true;
}

M3DObject * M3DMainFigureOptimizer::createTargetObject(const Vector & v)
{
	int figureId;

	if (problem == NULL)
		return NULL;

	return problem->createTargetObject(v, figureId, false);
}

void M3DMainFigureOptimizer::setOptimizerPosition(const Vector & start)
{
    if (method != NULL)
        delete method;

	method = new ConjugateGradientMethod(*wrapperProblem, start);
    ((ConjugateGradientMethod *) method)->setBrentLinearSearchBoundFactor(tuningWt(BrentLinearSearchBoundFactor));
    ((ConjugateGradientMethod *) method)->setBoundsFunction(boundsFunction);
    ((ConjugateGradientMethod *) method)->setParameterTolerance(epsilon);
#ifndef BINARY
    lastBestVal = problem->evaluate(start);
    *(lastBestX[0]) = start;
    lastPenalty = problem->getLastPenalty();
    ((ConjugateGradientMethod *) method)->getFunction().setEvaluationCost(0);
#endif
}


