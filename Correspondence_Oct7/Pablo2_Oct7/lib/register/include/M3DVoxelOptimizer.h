#ifndef M3D_BOUNDARY_DISPLACEMENT_OPTIMIZER_H
#define M3D_BOUNDARY_DISPLACEMENT_OPTIMIZER_H

#include "M3DObject.h"
#include "Match.h"


class NeighborList
{
private:
	int degree;
	int list[MAX_NEIGHBOR_COUNT];
public:
	int& operator[](int id) { return list[id % MAX_NEIGHBOR_COUNT]; }
	void set_degree(int k) { if (k>0 && k<=MAX_NEIGHBOR_COUNT) degree = k;}
	int get_degree() const { return degree; } 
};


class M3DVoxelOptimizer : public OptimizerBase
{
public:
	M3DVoxelOptimizer();
	~M3DVoxelOptimizer();

    M3DObject * getTargetObject() { return targetObject; }

    int getFigureId() { return figureId; }

	double getLastMatchValue() { return lastMatchValue; };
	double getLastObjectiveFunctionValue() { return lastObjectiveFunctionValue; }

    void setPenaltyWeight(double w);
    void setConstraintsPenaltyWeight(double w);

    bool initialize(Match * matchPtr, M3DObject * object, int figId, 
		int surfaceLevel = MATCH_POINT_CLOUD_SUBDIVISIONS,
		bool useWindowing = true);

	virtual bool isFinished() { return finished; }

    bool performIterations(int nIterations);

	M3DObject * getReferenceObject();

	virtual M3DObject * createTargetObject(const Vector & x) { return NULL; }

	virtual Function * getProblem() { return NULL; }
	virtual DifferentiableFunction * getDifferentiableProblem() { return NULL; }
	virtual void setOptimizerPosition(const Vector & x) { throw "not implemented";  }

private:
	M3DObject * targetObject;
	DMask * dmask;
	Image3D * targetImage;
	int figureId;
	int surfaceLevel;
	Match * match;
	SubdivBoundary * boundary;
	double * dVals;

	/////////////////////////
	// Internal variables to store boundary information (Bpoint list and neighbor list)
	int numPts;
	Bpoint * Bpnt;
	NeighborList * neighbors;
	/////////////////////////

	bool finished;

    double cumChange;
	double lastMatchValue;
	double lastObjectiveFunctionValue;

    double penaltyWeight;
    double constraintsPenaltyWeight;

	//parameters
	double sigma1, sigma2;
	bool windowIntensities;

	void initializeBoundary(M3DObject * object);
	void clearWorkSpace();
	void deleteTargetObject();
	double computeLastMatchValue();
	double computeLastPenaltyValue();
};


#endif

