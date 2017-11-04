#ifndef M3D_SPOKE_OPTIMIZER_H
#define M3D_SPOKE_OPTIMIZER_H

#include <vector>

#include "Match.h"
#include "M3DSpokeProblem.h"

#include "optima.h"
#include "ConjugateGradientMethod.h"

#include "OptimizerBase.h"


using namespace std;

#ifndef BINARY
extern const double DEFAULT_THRESHOLD;
extern const int DEFAULT_NUM_CONJGRAD_ITERATIONS;
extern const int DEFAULT_SCHEDULE_MULTIPLIER;
#endif

class M3DPGASpokeStats;


class M3DSpokeOptimizer : public OptimizerBase
{
public:
	M3DSpokeOptimizer(int appID = 0);
	virtual ~M3DSpokeOptimizer();

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
#endif
   
	double getTotalEvaluationCost() {return totalEvaluationCost;}

    double getLastBestVal() { return getLastBestVal(false); }	// See OptimizerBase.h
    double getLastBestVal(bool preview);

	// FIXME: This needs to be thought about
	// It's just a placeholder right now.
    bool initializeSpokePGAStage(M3DPGASpokeStats * _primPga, M3DObject * _meanFigPgaObject, int figId) { return false; }
	bool initialize(Match * _match, M3DObject * _referenceObject);
    

    void reset(Match * _match, M3DObject * _referenceObject);
	void readyNextRun();

	bool performIterations(int figureId);
	bool performIterations(int figureId, bool preview, bool verbosity = false);

	virtual M3DObject * createTargetObject(const Vector & x) { return NULL; }

	virtual Function * getProblem() { return NULL; }
	virtual DifferentiableFunction * getDifferentiableProblem() { return NULL; }
	virtual void setOptimizerPosition(const Vector & x) { throw "not implemented";  }

protected:
    Match * match;

	M3DObject * referenceObject;
	M3DObject * candidateObject;

	int applicationID;
	// The default value of applicationID is 0, which means the optimizer is used by
	// Pablo.  If applicationID = 1, the optimizer is used by VSkelTool, which requires
	// a different method of computing atom penalties.

    double penaltyWeight;
    double constraintsPenaltyWeight;
	double neighborPenaltyWeight;

    M3DPrimitive * referencePrimitive;
    M3DPrimitive * candidatePrimitive;
	bool optimizeFigure(int figId, bool preview, bool verbosity = false);


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
	double optimizeAtom(int figId, int atomId, M3DFigure * candidateFigure, bool preview);
};


#endif

