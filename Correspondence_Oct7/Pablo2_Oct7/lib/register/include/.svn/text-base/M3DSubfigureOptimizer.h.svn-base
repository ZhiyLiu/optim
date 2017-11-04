#ifndef M3D_SUBFIGURE_OPTIMIZER_H
#define M3D_SUBFIGURE_OPTIMIZER_H

#define CONJ_GRADIENT

#include "M3DSubfigureProblem.h"
#include "optima.h"
#ifdef CONJ_GRADIENT
#include "ConjugateGradientMethod.h"
#else
#include "EvolutionaryStrategy.h"
#endif


class M3DSubfigureOptimizer : public OptimizerBase
{
public:

    M3DSubfigureOptimizer();
    ~M3DSubfigureOptimizer();

    M3DObject * getReferenceObject()
    {
        if(problem == NULL)
            return NULL;

        return problem->getReferenceObject();
    }

    void setPenaltyWeight(double w);

    void initialize(Match * match, M3DObject * referenceObject,
                    int figureId);

    bool performIterations(int nIterations);

	virtual M3DObject * createTargetObject(const Vector & x) { return NULL; }

	virtual Function * getProblem() { return problem; }
	virtual DifferentiableFunction * getDifferentiableProblem() { return wrapperProblem; }

	virtual void setOptimizerPosition(const Vector & x) { throw "not implemented"; }

private:

    NumericalFunction * wrapperProblem;
    M3DSubfigureProblem * problem;
    GaussianSS * ss;

    double penaltyWeight;
};


#endif

