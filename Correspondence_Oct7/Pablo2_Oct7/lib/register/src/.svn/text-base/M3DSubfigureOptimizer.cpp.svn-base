#include <math.h>
#include <fstream>
#include "OptimizerBase.h"
#include "M3DSubfigureOptimizer.h"
#include "Tuning.h"

const double MIN_PENALTY_WEIGHT = -1.0;
const double MIN_PENALTY_WEIGHT_EXP = 0.0;
const double MAX_PENALTY_WEIGHT_EXP = 10.0;

using namespace std;

M3DSubfigureOptimizer::M3DSubfigureOptimizer()
{
	type = Subfigure;

    wrapperProblem = NULL;
    problem = NULL;
    ss = NULL;
    lastPenalty = 0.0;

    // Hardcoded to match initial slider value
    setPenaltyWeight(0.5);

	lastBestX = new Vector *[1];
	lastBestXSize = 1;
	*lastBestX = new Vector;
}

M3DSubfigureOptimizer::~M3DSubfigureOptimizer()
{
    if(wrapperProblem != NULL)
        delete wrapperProblem;
    if(problem != NULL)
        delete problem;
    if(ss != NULL)
        delete ss;
	// method belongs to OptimizerBase
}

void M3DSubfigureOptimizer::setPenaltyWeight(double w)
{
    double penaltyExp = MIN_PENALTY_WEIGHT_EXP + (MAX_PENALTY_WEIGHT_EXP - MIN_PENALTY_WEIGHT_EXP) * w;
    penaltyWeight = MIN_PENALTY_WEIGHT + pow(10.0, penaltyExp);
//    printf("Penalty Weight: %f\n", penaltyWeight);

    if(method != NULL)
#ifdef CONJ_GRADIENT
        ((M3DSubfigureProblem &) ((ConjugateGradientMethod *) method)
			->getFunction()).setPenaltyWeight(penaltyWeight);
#else
        ((M3DSubfigureProblem &) ((EvolutionaryStrategy *) method)
			->getFunction()).setPenaltyWeight(penaltyWeight);
#endif
}

void M3DSubfigureOptimizer::initialize(Match * match, M3DObject * referenceObject,
                                       int figureId)
{
    if(method != NULL)
        return;

    if(match == NULL || match->getReferenceObject() == NULL)
        return;

#ifdef CONJ_GRADIENT
    Vector start(5, 0.0, 0.0, 0.0, 0.0, 0.0);
    Vector metric(5, 0.01, 0.01, 0.01, 0.02, 0.02);
#else
    Vector start(5, 0.0, 0.0, 0.0, 0.0, 0.0);
    Vector metric(5, 0.01, 0.01, 0.01, 0.02, 0.02);
#endif

    if(problem != NULL)
        delete problem;

    problem = new M3DSubfigureProblem(match, referenceObject, figureId);
    if(problem == NULL)
        return;

#ifdef CONJ_GRADIENT
    if(wrapperProblem != NULL)
        delete wrapperProblem;

    wrapperProblem = new NumericalFunction(*problem, metric);
#else
    if(ss != NULL)
        delete ss;

    ss = new GaussianSS(start, metric);
    if(ss == NULL)
        return;
#endif

    problem->setPenaltyWeight(penaltyWeight);

    if(method != NULL)
        delete method;

#ifdef CONJ_GRADIENT
    method = new ConjugateGradientMethod(*wrapperProblem, start);
#else
    method = new EvolutionaryStrategy(*problem, *ss, 4, 6, SELECTION_MuPlusLambda);
    method->setSigmaFactor(1.0);
#endif

    lastBestVal = problem->evaluate(start);
    *(lastBestX[0]) = start;
    lastPenalty = problem->computePenalty(start);

#ifdef CONJ_GRADIENT
	((ConjugateGradientMethod *) method)->getFunction().setEvaluationCost(0);
	((ConjugateGradientMethod *) method)->setBrentLinearSearchBoundFactor(tuningWt(BrentLinearSearchBoundFactor));
#else
	((EvolutionaryStrategy *) method)->getFunction().setEvaluationCost(0);
#endif
}

bool M3DSubfigureOptimizer::performIterations(int nIterations)
{
    if(method == NULL)
        return false;

    int i = 0;
	while(i < nIterations)
    {
#ifdef CONJ_GRADIENT
    if (((ConjugateGradientMethod *) method)->isFinished())
		break;
#endif

#ifdef CONJ_GRADIENT
		((ConjugateGradientMethod *) method)->performIteration();
#else
		((EvolutionaryStrategy *) method)->performIteration();
#endif
        i++;
    }

    if(method->getBestEverValue() < lastBestVal)
    {
#ifdef CONJ_GRADIENT
        lastBestVal = ((ConjugateGradientMethod *) method)->getBestEverValue();
        *(lastBestX[0]) = ((ConjugateGradientMethod *) method)->getBestEverX();
#else
        lastBestVal = ((EvolutionaryStrategy *) method)->getBestEverValue();
        *(lastBestX[0]) = ((EvolutionaryStrategy *) method)->getBestEverX();
#endif
        lastPenalty = problem->computePenalty(*(lastBestX[0]));
    }

//    (lastBestX[0])->print();
	return true;
}

