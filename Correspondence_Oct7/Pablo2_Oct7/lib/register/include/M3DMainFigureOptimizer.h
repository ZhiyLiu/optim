#ifndef M3D_MAIN_FIGURE_OPTIMIZER_H
#define M3D_MAIN_FIGURE_OPTIMIZER_H

#include "M3DMainFigureProblem.h"
#include "optima.h"
#include "ConjugateGradientMethod.h"


class M3DMainFigureOptimizer : public OptimizerBase
{
public:

    M3DMainFigureOptimizer(Match * m, int count = 8);
    ~M3DMainFigureOptimizer();

    M3DObject * getReferenceObject() { return refObject; }
    int getTreeIndex() { return tree_index; }

#ifndef BINARY
	void setPenaltyWeightsInMatch(const int * penaltyNames,
		const double * penaltyWeights, int numPenalties);
#endif
	// Only needed for similarity-elongation optimization of binary images
	void setParameterCount(int count = 8) { numParameters = count; }

    virtual bool initialize(M3DObject * referenceObject, int treeIndex);

	bool setAddToResidueObject( const M3DObject* object ) {
		if(addToResidueObject != NULL) {
			delete addToResidueObject;
		}
		addToResidueObject	= object;
		if(problem != NULL) {
			problem->setAddToResidueObject(object);
		}
		return true;
	}

	int parameterCount() const { return numParameters; }

    bool performIterations(int nIterations);

	virtual M3DObject * createTargetObject(const Vector & x);

	virtual Function * getProblem() { return problem; }
	virtual DifferentiableFunction * getDifferentiableProblem() {
		return wrapperProblem;
	}

	virtual void setOptimizerPosition(const Vector & x);

protected:

	virtual const double getNumericalDerivativeEpsilon() = 0;

    virtual void reinitialize();

    M3DObject * refObject;
	const M3DObject * addToResidueObject;
    int tree_index;

    int numParameters;
    Vector epsilon;

    NumericalFunction * wrapperProblem;
    M3DMainFigureProblem * problem;
	Match * match;

    BoundsFunction * boundsFunction;

    int counter;
};


#endif

