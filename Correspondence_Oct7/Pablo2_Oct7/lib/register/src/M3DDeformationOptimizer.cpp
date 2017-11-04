#include <math.h>
#include "OptimizerBase.h"
#include "M3DDeformationOptimizer.h"
#include "M3DDeformationProblem.h"
#include "M3DPGAPrimitiveStats.h"
#include "M3DAtomPredictorQuad.h"
#include "M3DQuadPrimitive.h"
#include "Tuning.h"
#include "P3DControl.h"


//#define DEBUG

//#define DEBUG_2
//#define VERBOSE2

#define VERBOSE1
#define VERBOSE2
//#define DISPLAY_PENALTIES	/* See M3DDeformationProblem.cpp */


/*  Because of the skipping of atoms when the optimization gets close to the
    best fit, it is probably best to refine in random order.

    Note that the optimization order is consistent between runs when the debug
    version of the program is built, because srand() is called with a constant
    argument in pablo.cpp.
    */
#ifndef BINARY
#define RANDOMIZE		/* Comment out for serial optimization by atom number */
#endif

//#define DEFORMATION_OUTPUT_TEST


const double MIN_PENALTY_WEIGHT = -1;
const double MIN_PENALTY_WEIGHT_EXP = 0.0;
const double MAX_PENALTY_WEIGHT_EXP = 10.0;
const double PENALTY_WEIGHT_SCALE_FACTOR = 10.0;
const double PENALTY_SCALE_FACTOR =200;
const double CONSTRAINT_PENALTY_SCALE_FACTOR = 0.02;////1, 10, 100?

#ifndef BINARY
const int SKIPPED_PASSES = 1;
const double DEFAULT_THRESHOLD = 1.0e-3;
const int DEFAULT_NUM_CONJGRAD_ITERATIONS = 2;
const int DEFAULT_SCHEDULE_MULTIPLIER = 1;//xiaoxiao: originally 2
#endif

#ifdef DEFORMATION_OUTPUT_TEST
#include <fstream>
#endif


using namespace std;


extern int globalVerbosity;


M3DDeformationOptimizer::M3DDeformationOptimizer(int appID)
{
#ifdef DEBUG
	cout << "M3DDeformationOptimizer::M3DDeformationOptimizer()\n";
#endif
	type = Deformation;
	method = NULL;	// Note: method is a local variable in all cases below

	primPga = NULL;	// Atom PGA
	meanFigPgaObject = NULL;
	figResidualAtomDifObj = NULL;
	figureStageObj = NULL;
	startAtomObject = NULL;

	referenceObject = NULL;
	candidateObject = NULL;

	lastBestSigma = NULL;
	lastPenalty = 0.0;

	// Specify if the optimizer is used by Pablo (zero) or VSkelTool (nonzero).
	// This value will be passed on to the M3DDeformationProblem instance.
	applicationID = appID;

	counter = NULL;
	skip = NULL;
	lastResult = NULL;
	lastResultSize = 0;
#ifdef BINARY
	if (tuningWt(BpAtomRandomOrder) == 1)
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

void M3DDeformationOptimizer::setPenaltyWeightsInMatch(const int * penaltyNames,
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
			case GpAtomHistogramMatch:
				wt = penaltyWeights[i];
				break;
			case GpAtomMultiPatchMatch:
				wt = penaltyWeights[i];
				break;
			case AtomMahalanobisMatch:
				wt = PENALTY_SCALE_FACTOR*penaltyWeights[i];
				break;
			case GpAtomGeometryWeight:
				wt = PENALTY_WEIGHT_SCALE_FACTOR*penaltyWeights[i];
				break;
			case GpAtomNeighborWeight:
				wt = penaltyWeights[i];
				break;
			case GpAtomInterWeight:
				wt = CONSTRAINT_PENALTY_SCALE_FACTOR*penaltyWeights[i];
				break;

			case GpAtomProfileMatch:
				// Never set; see Match.cpp and gpTuning.cpp
			default:
				cout << "Unknown penalty weight: " << tuningShortName(penaltyNames[i]) << endl;
				return;
		}
		match->setWeight(type, penaltyNames[i], wt);
	}
}

#endif

M3DDeformationOptimizer::~M3DDeformationOptimizer()
{
#ifdef DEBUG
	cout << "M3DDeformationOptimizer::~M3DDeformationOptimizer()" << endl;
#endif
	delete_internal();
}

void  M3DDeformationOptimizer::delete_internal()
{
	int numLastBestSigma, i;

	if(referenceObject != NULL)
		delete referenceObject;
	referenceObject = NULL;

	if(candidateObject != NULL)
		delete candidateObject;
	candidateObject = NULL;


	if (meanFigPgaObject != NULL)
		delete meanFigPgaObject;
	meanFigPgaObject = NULL;


	if (figResidualAtomDifObj !=NULL)
		delete figResidualAtomDifObj;
	figResidualAtomDifObj = NULL;

	if (figureStageObj != NULL)
		delete figureStageObj;
	figureStageObj = NULL;


	if (startAtomObject !=NULL)
		delete startAtomObject;
	startAtomObject = NULL;


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

	if(lastBestSigma)
	{
		numLastBestSigma = lastBestSigma->size();
		for(i = 0; i < numLastBestSigma; i++) 
		{
			if((*lastBestSigma)[i] != NULL)
				delete (*lastBestSigma)[i];
		}
		delete lastBestSigma;
		lastBestSigma = NULL;
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
	if (tuningWt(BpAtomRandomOrder) == 1) {
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


bool M3DDeformationOptimizer::addDiffToPrimitive(M3DFigure * figure, const int atomId,   M3DFigure * alignedToFigure, const M3DPrimitive * deltaPrim)
{
	M3DPrimitive * prim = NULL;

	int numPrims;
	M3DPrimitive * alignSets[5];
	M3DPrimitive * alignToSets[5];

	M3DPrimitive * deltaPrimAlign;

	prim = figure->getPrimitivePtr(atomId);

	//Align the atom by its neigbors
	figure->getAtomsNeighbors(atomId, numPrims, alignSets, ALL_FIRST_NEIGHBORS);
	alignedToFigure->getAtomsNeighbors(atomId, numPrims, alignToSets, ALL_FIRST_NEIGHBORS);		

	//align by average normals
	alignSets[numPrims] = figure->getPrimitivePtr(atomId);
	alignToSets[numPrims] = alignedToFigure->getPrimitivePtr(atomId);//include itself	
	numPrims++;

	deltaPrimAlign = deltaPrim->copyPtr();

/*	//apply the inverse transform to the deltaM
	SimilarityTransform3D xform =  M3DPrimitive::atomsAlignmentSimtrans(numPrims,alignToSets,alignSets);

	//Apply rotation to the delta atom. For the difference between atom and it's neighbor, rotation is enough
    //	deltaPrimAlign->rotateBy(xform.getRotation());
	
	// //	deltaPrimAlign->scaleBy(scale);	
	//Apply translations too, if the stats are the difference between two figures. 
	//	cout<<"scale="<< scale<<" \n";
*/
	if (! M3DPrimitive::composePrimitives(deltaPrimAlign,prim) )	{ 
		cout<<"error at composePrimitives in M3DDeformationOptimizer\n"<<endl;
		return false;
	}
	delete deltaPrimAlign;

	return true;
}


M3DPrimitive *  M3DDeformationOptimizer ::getFigResidualAtomDifference(const M3DPrimitive * figStageAtom,
		const M3DPrimitive * predictedAtom, int figId, int atomId)
{
	M3DPrimitive * primFigOffset = figStageAtom->copyPtr();

	//primOffset = figStageAtom - predicatedAtom
	M3DPrimitive::subtractPrimitives( primFigOffset,predictedAtom);


	if (primPga ==NULL)
		return false;

	PGData * pgData = primPga->getPGDataPtr(figId,atomId);
	int numPGs = pgData->numPGs; //4,5
	int lenPG = pgData->lenPG;  //9

	Vector coeff = Vector (numPGs);
	int i,j;

	Vector vecFigOffset;
	M3DPrimitive::projectAtom2Vec(primFigOffset, vecFigOffset);

	for (i = 0; i < numPGs; i++)
		for(j = 0; j < lenPG; j++)			
			coeff(i) += (vecFigOffset(j) * ((pgData->pgVec[i])[j]));

	Vector vecProjAtom = Vector(9,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0);

	for(j = 0; j < lenPG; j++) 
		for (i = 0; i < numPGs; i++)				
			vecProjAtom(j) += coeff(i)* ((pgData->pgVec[i])[j]);			

	M3DPrimitive * projPrim;
	if (figStageAtom->type() == M3D_END_PRIMITIVE)
		projPrim = new M3DQuadEndPrimitive();	
	else
		projPrim = new M3DQuadPrimitive();


	projPrim->vec2Atom(vecProjAtom);     //  pojection coefficient i = cross product with  v(i)
	//save the project prim for alignments later

	M3DPrimitive::subtractPrimitives( primFigOffset,projPrim); 


	delete projPrim;
	return  primFigOffset;	

}



bool M3DDeformationOptimizer::initializeAtomPGAStage(M3DPGAPrimitiveStats * _primPga,int figId,M3DObject * _referObj)
{
	if (_referObj!=NULL)
	{
		referenceObject = _referObj->assign();		
		figureStageObj = _referObj->assign();//save for alignments in atom PGA stage
		startAtomObject = _referObj->assign();
		candidateObject = referenceObject->assign();
	}

	if (_primPga != NULL){
		meanFigPgaObject = _primPga->getAlignReferenceObj()->assign();	
		primPga = _primPga;
		return true;
	}

	return false;
}


bool M3DDeformationOptimizer::initialize(Match * _match, M3DObject * _referenceObject)

{
//	int numPrimitives,i;
//	M3DPrimitive * primPtr;
	Image3D * image;


#ifdef DEBUG
	cout << "M3DDeformationOptimizer::initialize()\n";
#endif
	if(_match == NULL || _referenceObject == NULL)
		return false;

	delete_internal();

	referenceObject = _referenceObject->assign();

	figureStageObj = _referenceObject->assign();//save for alignments in atom PGA stage
	startAtomObject = _referenceObject->assign();

	candidateObject = _referenceObject->assign();
	match = _match;

     match->resetReferenceObject(referenceObject);

	image = match->getTargetImage();
	if(image == NULL)
		return false;


//why set to zero?
//	numPrimitives = candidateObject->getPrimitiveCount();

/*	for (i = 0; i < numPrimitives; i++)
	{
		primPtr = candidateObject->getPrimitivePtr(i);
		primPtr->setX(0.0, 0.0, 0.0);
		primPtr->setQ(Quat(0.0, 0.0, 0.0, 0.0));
		primPtr->setR(0.0);
		primPtr->setTheta(0.0);

		if (primPtr->type() == M3D_END_PRIMITIVE)
			(dynamic_cast<M3DEndPrimitive *>( primPtr)) -> setElongation(0.0);
	}
*/
	lastFigureId = -1;
	nextAtom = 0;

#ifdef DEFORMATION_OUTPUT_TEST
	const int NUM_TICKS = 100;
	double START_TICK = -0.25;
	double END_TICK = 0.25;
	double TICK_SIZE = (END_TICK - START_TICK) / ((double) NUM_TICKS);

	ofstream outFile;
	int row, col;
	double x, val;
	int j;

	M3DFigure * targetFigure = candidateObject->getFigurePtr(0);
	if(targetFigure == NULL)
		return;
	numPrimitives = targetFigure->getPrimitiveCount();

	for (i = 0; i < numPrimitives; i++)
	{
		primPtr = targetFigure->getPrimitivePtr(i);
		if(primPtr == NULL || !(primPtr->isSelected()))
			continue;

		M3DDeformationProblem problem (match, referenceObject, candidateObject, 0, i, applicationID);

		int vecsize = 8;
		if (primPtr->type() == M3D_END_PRIMITIVE)
			vecsize = 9;
		Vector zv(vecsize);
		for (int k = 0; k < vecsize; k++)
			zv(k) = 0.0;

		outFile.open("xtrans.dat");
		for(j = 0; j < NUM_TICKS; j++)
		{
			x = START_TICK + ((double) j) * TICK_SIZE;
			//Vector v(8, x, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
			Vector v(zv);
			v(0) = x;
			val = problem.evaluate(v);

			outFile << x << "\t\t" << val << endl;
		}
		outFile.close();

		outFile.open("ytrans.dat");
		for(j = 0; j < NUM_TICKS; j++)
		{
			x = START_TICK + ((double) j) * TICK_SIZE;
			//Vector v(8, 0.0, x, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
			Vector v(zv);
			v(1) = x;
			val = problem.evaluate(v);

			outFile << x << "\t\t" << val << endl;
		}
		outFile.close();

		outFile.open("ztrans.dat");
		for(j = 0; j < NUM_TICKS; j++)
		{
			x = START_TICK + ((double) j) * TICK_SIZE;
			//Vector v(8, 0.0, 0.0, x, 0.0, 0.0, 0.0, 0.0, 0.0);
			Vector v(zv);
			v(2) = x;
			val = problem.evaluate(v);

			outFile << x << "\t\t" << val << endl;
		}
		outFile.close();

		outFile.open("xrot.dat");
		for(j = 0; j < NUM_TICKS; j++)
		{
			x = START_TICK + ((double) j) * TICK_SIZE;
			//Vector v(8, 0.0, 0.0, 0.0, x, 0.0, 0.0, 0.0, 0.0);
			Vector v(zv);
			v(3) = x;
			val = problem.evaluate(v);

			outFile << x << "\t\t" << val << endl;
		}
		outFile.close();

		outFile.open("yrot.dat");
		for(j = 0; j < NUM_TICKS; j++)
		{
			x = START_TICK + ((double) j) * TICK_SIZE;
			//Vector v(8, 0.0, 0.0, 0.0, 0.0, x, 0.0, 0.0, 0.0);
			Vector v(zv);
			v(4) = x;
			val = problem.evaluate(v);

			outFile << x << "\t\t" << val << endl;
		}
		outFile.close();

		outFile.open("zrot.dat");
		for(j = 0; j < NUM_TICKS; j++)
		{
			x = START_TICK + ((double) j) * TICK_SIZE;
			//Vector v(8, 0.0, 0.0, 0.0, 0.0, 0.0, x, 0.0, 0.0);
			Vector v(zv);
			v(5) = x;
			val = problem.evaluate(v);

			outFile << x << "\t\t" << val << endl;
		}
		outFile.close();

		outFile.open("scale.dat");
		for(j = 0; j < NUM_TICKS; j++)
		{
			x = START_TICK + ((double) j) * TICK_SIZE;
			//Vector v(8, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, x, 0.0);
			Vector v(zv);
			v(6) = x;
			val = problem.evaluate(v);

			outFile << x << "\t\t" << val << endl;
		}
		outFile.close();

		outFile.open("theta.dat");
		for(j = 0; j < NUM_TICKS; j++)
		{
			x = START_TICK + ((double) j) * TICK_SIZE;
			//Vector v(8, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, x);
			Vector v(zv);
			v(7) = x;
			val = problem.evaluate(v);

			outFile << x << "\t\t" << val << endl;
		}
		outFile.close();

		if (vecsize == 9)
		{
			outFile.open("eta.dat");
			for(j = 0; j < NUM_TICKS; j++)
			{
				x = START_TICK + ((double) j) * TICK_SIZE;
				//Vector v(8, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, x);
				Vector v(zv);
				v(8) = x;
				val = problem.evaluate(v);

				outFile << x << "\t\t" << val << endl;
			}
			outFile.close();
		}
	}
#endif
	return true;
}

void M3DDeformationOptimizer::reset(Match * _match, M3DObject * _referenceObject)
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
bool M3DDeformationOptimizer::performIterations(int figureId)
{
	cout << "Error: wrong version of M3DDeformationOptimizer::performIterations() called\n";
	return false;
}

/* This function optimizes a single atom per call from the selected figure, so the
   display may be updated after each atom transformation.  When called repeatedly,
   successive atoms will be optimized.  The order of atoms will be fixed, unless
   RANDOMIZE is defined at compile time in regular Pablo, or for Binary Pablo, if
   tuningWt(BpAtomRandomOrder) == 1.  Atoms may be skipped when they appear
   to not need optimization (see optimizeFigure() below).  Once all atoms are
   considered refined, the function returns true.  This function is most efficient
   when one figure is fully optimized, before moving on to another figure.

   The verbosity argument controls the amount of output printed, except with
   Binary Pablo.  When it is true, a dot will be printed for each atom optimized.
   When false, there will be no output.
   */
bool M3DDeformationOptimizer::performIterations(int figureId, bool preview,
		bool verbosity)
{
	int numFigs;

	if(match == NULL || referenceObject == NULL || candidateObject == NULL)
		return false;

	numFigs = candidateObject->getFigureCount();

#ifdef BINARY
	if (tuningWt(BpAtomRSradPenalty) != 0.0 || tuningWt(BpAtomDeltaNeighborMatch) != 0.0)
		match->setUpAtomPredictor(figureId);

	if (tuningWt(BpAtomRandomOrder) == 1)
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

inline int M3DDeformationOptimizer::selectRandomAtom(int numPrims)
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

#ifdef DEBUG
		cout << "Allocated perm[]: ";
		for (i = 0; i < numPrims; i++)
			cout << perm[i] << ' ';
		cout << endl;
#endif	/* VERBOSE */
	}

	// Stop if in previous pass all atoms were either skipped or scheduled to skip
	skipCounter = 0;
	if (nextAtom == 0) {
		for (i = 0; i < numPrims; i++) {
			if (counter[i] < 0)
				skipCounter++;
		}
	}
#ifdef VERBOSE2
	if(nextAtom == 0)
		cout << "skipCounter = " << skipCounter << '\n';
#endif
	if (skipCounter == numPrims) {
		allSkipped++;
#ifdef VERBOSE
		cout << "All skipped = " << allSkipped << '\n';
#endif
	}

	// Select next atom to optimize
	primId = perm[nextAtom];
	nextAtom++;
	if (nextAtom == numPrims)
		nextAtom = 0;

#ifndef DISPLAY_PENALTIES
#ifdef VERBOSE
	if (nextAtom == 1)
		cout << "==================\n";
	else
		cout << "------------------\n";
#endif	/* VERBOSE */
#endif	/* DISPLAY_PENALTIES */

#ifdef VERBOSE
	cout << "counter: ";
	for (i = 0; i < numPrims; i++)
		cout << counter[i] << ' ';
	cout << '\n';
	cout << "skip: ";
	for (i = 0; i < numPrims; i++)
		cout << skip[i] << ' ';
	cout << '\n';
#endif

	return primId;
}

inline int M3DDeformationOptimizer::selectAtom(int numPrims)
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
#ifdef VERBOSE
	if(nextAtom == 0)
		cout << "skipCounter = " << skipCounter << '\n';
#endif
	if (skipCounter == numPrims) {
		allSkipped++;
#ifdef VERBOSE
		cout << "All skipped = " << allSkipped << '\n';
#endif
	}

	// Select next atom to optimize
	primId = nextAtom;
	nextAtom++;
	if (nextAtom == numPrims)
		nextAtom = 0;

#ifndef DISPLAY_PENALTIES
#ifdef VERBOSE
	if (nextAtom == 1)
		cout << "==================\n";
	else
		cout << "------------------\n";
#endif	/* VERBOSE */
#endif	/* DISPLAY_PENALTIES */

#ifdef VERBOSE
	cout << "counter: ";
	for (i = 0; i < numPrims; i++)
		cout << counter[i] << ' ';
	cout << '\n';
	cout << "skip: ";
	for (i = 0; i < numPrims; i++)
		cout << skip[i] << ' ';
	cout << '\n';
#endif

	return primId;
}

void M3DDeformationOptimizer::readyNextRun()
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

bool M3DDeformationOptimizer::optimizeFigure(int figId, bool preview, bool verbosity,
		bool allAtomMoveTogather)

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
#ifdef VERBOSE
		cout << "initializing internal optimization structures\n";
#endif

		lastBestX = new Vector *[numPrims];
		lastBestXSize = numPrims;
		
		for (primId = 0; primId < numPrims; primId++)
		{		
			// The size of x is decided by numPGs
			if (primPga == NULL) {
				M3DPrimitive *  currPrimitive = candidateFigure->getPrimitivePtr(primId);
				Vector * init;
     		
				if ((tuningWt(AtomEtaFactor) == 0.0)|| currPrimitive->type() == M3D_STANDARD_PRIMITIVE)
					init = new Vector(8, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
				else // end primitive				
					init = new Vector(9, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,0.0);
				lastBestX[primId] = init;
			}
			else {
				int vecsize = (primPga->getPGDataPtr(figId,primId))->numPGs;
				lastBestX[primId] = new Vector(vecsize);
				lastBestX[primId]->setAll(0.0);		
			}
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
			if (tuningWt(BpAtomRandomOrder) == 1.0)
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
		if (candidatePrimitive == NULL){ // xiaoxiao:       
			//	|| ! candidatePrimitive->isSelected()) 
			// Unselected atoms are skipped
			counter[primId] = -1;
			primId = -1;
		}

	}

#ifdef BINARY
	if (allSkipped >= tuningWt(BpAtomSkippedPasses))
#else
		if (allSkipped >= SKIPPED_PASSES) //if all the atoms were skipped last round 
#endif
			return true;	// Optimization completed

	if (skip[primId] > 0) {
		// Skip optimization of this atom
		skip[primId]--;

		cout << "Skipping atom " << primId << '\n';

		return false;
	}

	referencePrimitive = referenceFigure->getPrimitivePtr(primId);

	if (referencePrimitive == NULL)
		return false;

	// Optimize a single atom

	cum_change = optimizeAtom(figId, primId, candidateFigure, preview,
			allAtomMoveTogather);

	// Report the accumulated change over the whole figure
#ifdef DISPLAY_PENALTIES
	cout << "After optimizing atom " << primId << ", cumulative change for figure "
		<< figId << " = " << cum_change << endl;
#endif

	double d_abs = fabs(cum_change - lastResult[primId]);
#ifdef VERBOSE2
	cout << "d_abs = " << d_abs;
#endif
#ifdef BINARY
	if (d_abs < tuningWt(BpAtomSkipThreshold)) {
		counter[primId]--;
		skip[primId] = -((int) tuningWt(BpAtomSkipMultiplier))*counter[primId];
#else
		if (d_abs < threshold) {
			counter[primId]--;
			skip[primId] = -scheduleMultiplier * counter[primId];
#endif
#ifdef VERBOSE
#ifdef BINARY
			cout << " < " << tuningWt(BpAtomSkipThreshold);
#else
			cout << " < " << threshold;
#endif
#endif
		}
		else
			counter[primId] = 0;
#ifdef VERBOSE2
		cout << '\n';
#endif

		lastResult[primId] = cum_change;


		return false;
	}

	double M3DDeformationOptimizer::getLastBestVal(bool preview)
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


	double M3DDeformationOptimizer::optimizeAtom(int figId, int atomId,	M3DFigure * candidateFigure,
			bool preview, bool allAtomMoveTogather)
	{
		int k,i,vecsize;
		double initVal, val;
		double cum_change;
		Vector * start;

		if (primPga!=NULL){// Optimize a single atom:if there is no atomPGA, starting from the last best answer from the last iteration;
		// if there is atomPGA, starting from zero coefficients every new iteration for that atom; 
			lastBestX[atomId]->setAll(0.0);	
		}		
		
		start = lastBestX[atomId];
		vecsize = start->size();	 	

		M3DDeformationProblem * problem =  new M3DDeformationProblem(match, referenceObject, candidateObject,figId,atomId,primPga, applicationID);		

		Vector epsilon(vecsize);
		epsilon.setAll(tuningWt(AtomEpsilon)); 	

		NumericalFunction function(* problem, epsilon);
		ConjugateGradientMethod m(function, *start);
		m.setBrentLinearSearchBoundFactor( tuningWt(BrentLinearSearchBoundFactor) );


		M3DDeformationProblemBoundsFunction boundsFunction(problem);
		m.setBoundsFunction(&boundsFunction);	

		if (primPga != NULL){
			initVal = problem->evaluate(*start);
			lastBestVal = initVal;
		}else{
			initVal = m.getBestEverValue();	
		}

		

#ifdef VERBOSE1
			cout << "\n AtomId = " << atomId << "  initVal = " << initVal << endl;
#endif

#ifdef BINARY
		for(k = 0; (k < tuningWt(BpAtomConjugateGradientIterations)) && ! m.isFinished(); k++)
#else
			for(k = 0; (k < numConjGradientIterations) && ! m.isFinished(); k++)
#endif
				m.performIteration();

		totalEvaluationCost += problem->getEvaluationCost();
		//cout << "OptimizeAtom - evaluationCost: " << problem.getEvaluationCost() << endl;

		Vector * endX = lastBestX[atomId];

		val = m.getBestEverValue();


		if (val < initVal)
		{
			Vector result = m.getBestEverX();
			for (k = 0; k < vecsize; k++)
			{
				(*lastBestX[atomId])(k) = result(k);
			}

#ifdef VERBOSE2		
			cout << " bestx=";
			for (i = 0; i < vecsize; i++)
				cout << ' ' << result(i);
			cout << '\n';			
			cout << "   val = " << val << '\n';
#endif


#ifdef VERBOSE1	
			if (primPga!=NULL){			
				double mpv = 0;  //mahamPenaltyVal		
				for (i = 0; i < vecsize; i++)
					mpv += result(i) * result(i);

				cout << " mahal Dist: "<< mpv<<'\n';		
			}
#endif
			cum_change = fabs(initVal - val);
			lastBestVal = val;
			
		}
		else
			cum_change = 0.0;	

       // the candidate Object is updated(modified) inside the createTargetPrimitive
	     problem->createTargetPrimitive( *endX );
	

	//	if (primPga!= NULL) {		
		//while having the statistics, atom by atom iteration, need to update the reference after each iteration
	//		*referencePrimitive	= *candidatePrimitive;
	//	}



#ifndef BINARY					
		lastPenalty = problem->getGeometricPenalty();// return geomPenalty, which is the best result's atom maham penalty
#else
		lastPenalty = 0;
#endif

		delete problem;

#ifdef DEBUG_2
		//save the prediction to see whether it's right
		char tempName[50];
		static int counterAtomIterations=0;
		P3DControl * gControl= new P3DControl(0);	
		M3DObject * saveTemp = candidateObject->assign();	
		gControl->newModel(saveTemp);
		sprintf(tempName,"./candidate%d%d.m3d",counterAtomIterations++,atomId);
		gControl->write(tempName);
       
		delete saveTemp;
#endif	

		return cum_change;
	}






