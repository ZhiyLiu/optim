#ifndef M3D_DEFORMATION_OPTIMIZER_H
#define M3D_DEFORMATION_OPTIMIZER_H

#include <vector>

#include "Match.h"
#include "M3DDeformationProblem.h"

#include "optima.h"
#include "ConjugateGradientMethod.h"

#include "OptimizerBase.h"


using namespace std;

#ifndef BINARY
extern const double DEFAULT_THRESHOLD;
extern const int DEFAULT_NUM_CONJGRAD_ITERATIONS;
extern const int DEFAULT_SCHEDULE_MULTIPLIER;
#endif


class M3DDeformationOptimizer : public OptimizerBase
{
public:
	M3DDeformationOptimizer(int appID = 0);
	virtual ~M3DDeformationOptimizer();

	M3DObject * getReferenceObject() { return referenceObject; }
	M3DObject * getCandidateObject() { return candidateObject; }

#ifndef BINARY
	void setPenaltyWeightsInMatch(const int * penaltyNames,
	const double * penaltyWeights, int numPenalties);

    void setConjGradientIterations(int i) { numConjGradientIterations = i; }
    void setScheduleMultiplier(int mult) { scheduleMultiplier = mult; }
    void setThreshold(double t) { threshold = t; }

    int getConjGradientIterations() { return numConjGradientIterations; }
    int getScheduleMultiplier() { return scheduleMultiplier; }
    double getThreshold() { return threshold; }

	// After the figure stage, apply the mean difference to the result
#endif
    
	/**
	
	T0 calculate figure residual atom difference  by the figure stage result
	before the first iteration in atom stage, and once for all the iterations.
	"predicatedAtom + residual = figStageAtom;"	
	*/
	M3DPrimitive * getFigResidualAtomDifference(const M3DPrimitive * figStageAtom,
					const M3DPrimitive * predictedAtom, int figId, int atomId);
	
					
	/**
    To add  some delta atom ( atom difference)  to an atom with alignment by its neighbors.
    Need both figures for alignment.
	The resulting atom is kept in "figure" (updated). "alignedtoFigure" 
	is the figure where you calculating the delta atom
	*/
	static bool addDiffToPrimitive(M3DFigure * figure, const int atomId,  M3DFigure * alignedToFigure, const M3DPrimitive * deltaPrim);
	
	/**
	Is still under debugging.The whole strategy may moved in to M3DDeformationProblem: doApplyVectorPGA()
	To calculate the predictedAtom to be ready for deform.
	predictedAtom = interpolatedAtom + mean Atom difference from the atom PGA
	*/
//    M3DPrimitive * getPredictedAtomWithAtomStatsMean(M3DFigure * figure,int figId, int atomId, M3DPGAPrimitiveStats * _primPga= NULL);
    

	//bool applyAllVectors();
   
	double getTotalEvaluationCost() {return totalEvaluationCost;}

    double getLastBestVal() { return getLastBestVal(false); }	// See OptimizerBase.h
    double getLastBestVal(bool preview);

	// MeanFigPgaObject is the mean object of all traning models to do alignmend in atom PGA stage
    bool initializeAtomPGAStage(M3DPGAPrimitiveStats * _primPga, int figId, M3DObject * _referObj = NULL);
	
	bool initialize(Match * _match, M3DObject * _referenceObject);
    

    void reset(Match * _match, M3DObject * _referenceObject);
	void readyNextRun();

	bool performIterations(int figureId);
	bool performIterations(int figureId, bool preview, bool verbosity = false);

	virtual M3DObject * createTargetObject(const Vector & x) { return NULL; }

	virtual Function * getProblem() { return NULL; }
	virtual DifferentiableFunction * getDifferentiableProblem() { return NULL; }
	virtual void setOptimizerPosition(const Vector & x) { throw "not implemented";  }
	M3DObject * getMeanFigPgaObject(){return meanFigPgaObject;}
	M3DObject * getFigureStageResult(){return figureStageObj;}
   
protected:
	M3DObject * meanFigPgaObject;	// PGA mean obj , save for Atom allignment  
	
	M3DObject * startAtomObject;//for debug, delete later	
	
	M3DObject * figResidualAtomDifObj;//save the fig residual for each atom by the first iteration of atom stage, and once for all
    
	M3DObject * figureStageObj;//save the projcted fig offset for each atom for alignment for adding up he fig residual 
	
	M3DPGAPrimitiveStats * primPga; //atomPGA
    
	Match * match; 
	
    vector<Vector *> * lastBestSigma;

	M3DObject * referenceObject;
	M3DObject * candidateObject;

    // HACK: current offset from the original refObject
    // The Euler angles for the rotation are saved in the x,y,z of the quaternions
    M3DObject * deltaObject;

	int applicationID;
	// The default value of applicationID is 0, which means the optimizer is used by
	// Pablo.  If applicationID = 1, the optimizer is used by VSkelTool, which requires
	// a different method of computing atom penalties.

    double penaltyWeight;
    double constraintsPenaltyWeight;
	double neighborPenaltyWeight;

    M3DPrimitive * referencePrimitive;
    M3DPrimitive * candidatePrimitive;
    M3DPrimitive * deltaPrimitive;
	bool optimizeFigure(int figId, bool preview, bool verbosity = false,
		bool allAtomMoveTogather = false);


    int selectAtom(int numPrims);
    int selectRandomAtom(int numPrims);

    // deletes internal structure methods (and deletes also structures it points to)
    void delete_internal();

    double dx, dy, dz;

    short * perm;
    short * counter;
    double * lastResult;
	int lastResultSize;
    short * skip;
	short allSkipped;
    short nextAtom;
	short lastFigureId;

#ifndef BINARY
    int numConjGradientIterations;
    int scheduleMultiplier;
    double threshold;
#endif

    double totalEvaluationCost;

private:

//#ifdef BINARY
	double optimizeAtom(int figId, int atomId, M3DFigure * candidateFigure, bool preview,
		bool allAtomMoveTogather = false);
//#else
//	double optimizeAtom(int figId, int atomId, M3DFigure * candidateFigure, bool preview);
//#endif

};


#endif

