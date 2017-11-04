#ifndef M3D_REGISTRATION_OPTIMIZER_H
#define M3D_REGISTRATION_OPTIMIZER_H

#include "M3DObject.h"
#include "M3DRegistrationProblem.h"

#include "optima.h"
#include "ConjugateGradientMethod.h"
#include "OptimizerBase.h"


class M3DRegistrationOptimizer : public OptimizerBase
{
public:
    M3DRegistrationOptimizer(bool adaptive = false);
    ~M3DRegistrationOptimizer();

    virtual M3DObject * getReferenceObject()
    {
        if(problem == NULL || problem->getMatch() == NULL)
            return NULL;

        return problem->getMatch()->getReferenceObject();
    }
    SimilarityTransform3D & getLastBestTransformation() { return lastBestTransformation; }

    void setPenaltyWeight(double w);

    void initialize(Match * match);

    bool performIterations(int nIterations);

	virtual M3DObject * createTargetObject(const Vector & x);

	virtual Function * getProblem() { return problem; }
	virtual DifferentiableFunction * getDifferentiableProblem() { return wrapperProblem; }

	virtual void setOptimizerPosition(const Vector & x);

	bool isAdaptive() { return adapt; }

protected:

	virtual void reinitialize();

	int numParameters;

	NumericalFunction * wrapperProblem;
    M3DRegistrationProblem * problem;

	Vector epsilon;
	int counter;

    double penaltyWeight;
    BoundsFunction * boundsFunction;

    SimilarityTransform3D lastBestTransformation;

	bool adapt;
};


#endif

