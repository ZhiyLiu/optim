#include <fstream>
#include <math.h>
#include "OptimizerBase.h"
#include "M3DRegistrationOptimizer.h"
#include "Tuning.h"

//#define DEBUG


const double MIN_PENALTY_WEIGHT = -1.0;
const double MIN_PENALTY_WEIGHT_EXP = 0.0;
const double MAX_PENALTY_WEIGHT_EXP = 10.0;
const double PENALTY_WEIGHT_SCALE_FACTOR = 80.0;//0.5;

using namespace std;

M3DRegistrationOptimizer::M3DRegistrationOptimizer(bool adaptive)
{
	type = Registration;
	adapt = adaptive;

    problem = NULL;
	wrapperProblem = NULL;
	boundsFunction = NULL;
	counter = 0;
	numParameters = 7;	//default when no PGs

    // Hardcoded to match initial slider value
    setPenaltyWeight(0.5);

    lastPenalty = 0.0;
    lastBestTransformation.setToIdentity();

	lastBestX = new Vector *[1];
	lastBestXSize = 1;
	*lastBestX = new Vector;
}

M3DRegistrationOptimizer::~M3DRegistrationOptimizer()
{
    if(problem != NULL)
        delete problem;
	if(wrapperProblem != NULL)
		delete wrapperProblem;
	if(boundsFunction != NULL)
		delete boundsFunction;
	// method belongs to OptimizerBase
}

void M3DRegistrationOptimizer::setPenaltyWeight(double w)
{
//    double penaltyExp = MIN_PENALTY_WEIGHT_EXP + 
//        (MAX_PENALTY_WEIGHT_EXP - MIN_PENALTY_WEIGHT_EXP) * w;
//    penaltyWeight = MIN_PENALTY_WEIGHT + pow(10.0, penaltyExp);
//    printf("Penalty Weight: %f\n", penaltyWeight);
      penaltyWeight = w * PENALTY_WEIGHT_SCALE_FACTOR;

    if(problem != NULL)
        problem->setPenaltyWeight(penaltyWeight);
}

void M3DRegistrationOptimizer::initialize(Match * match)
{
    Image3D * image;

#ifdef DEBUG
	cout << "M3DRegistrationOptimizer::initialize()\n";
#endif
    if(match == NULL || match->getReferenceObject() == NULL)
        return;

	if (isAdaptive()) {
		if(problem == NULL)
			return;

		problem->initialize(match, penaltyWeight, 0.0);
	}
	else {
		if(problem != NULL)
			delete problem;

		problem = new M3DRegistrationProblem(match, 0.0);
		if(problem == NULL)
			return;

		problem->setPenaltyWeight(penaltyWeight);
	}

    image = match->getTargetImage();
    if(image == NULL)
        return;

	Vector start(numParameters);
	epsilon.setSize(numParameters);

	if (isAdaptive()) {
		start.setAll(0.0);
		epsilon.setAll(1.0 / (64.0 * M3DRegistrationProblem::ENSEMBLE_PGA_SCALE_FACTOR));
	}
	else
		epsilon = Vector(7, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.003);

	if(wrapperProblem != NULL)
		delete wrapperProblem;

	wrapperProblem = new NumericalFunction(*problem, epsilon);

	setOptimizerPosition(start);
}

void M3DRegistrationOptimizer::reinitialize()
{
	if(method != NULL)
		delete method;
	if (wrapperProblem != NULL)
		delete wrapperProblem;

	epsilon *= 0.5;
	wrapperProblem = new NumericalFunction(*problem, epsilon);

	setOptimizerPosition(*(lastBestX[0]));
}

bool M3DRegistrationOptimizer::performIterations(int nIterations)
{
#ifdef DEBUG
	cout << "M3DRegistrationOptimizer::performIterations()\n";
#endif
    if(method == NULL)
        return false;

    int i = 0;
    while(i < nIterations && ! isFinished())
    {
        ((ConjugateGradientMethod *) method)->performIteration();
        i++;
    }

    if(method->getBestEverValue() < lastBestVal)
    {
        lastBestVal = method->getBestEverValue();
        *(lastBestX[0]) = ((ConjugateGradientMethod *) method)->getBestEverX();
		if (isAdaptive()) {
			lastPenalty = problem->getPenalty();

			if (lastBestObject != NULL)
				delete lastBestObject;
			lastBestObject = createTargetObject(*(lastBestX[0]));

#ifdef DEBUG
			double sumX = 0.0;
			for(int i = 0; i < (lastBestX[0])->size(); i++)
				std::cout << (*(lastBestX[0]))(i)*M3DRegistrationProblem::ENSEMBLE_PGA_SCALE_FACTOR << ' ';
			std::cout<<endl;
#endif

		}
		else {
			lastPenalty = problem->computePenalty(*(lastBestX[0]));
			problem->computeTransformation(lastBestTransformation, *(lastBestX[0]));
		}
    }

	if (isAdaptive()) {
#ifdef DEBUG
	    if(method->isFinished())
	    {
			cout << "method->isFinished == true. counter == " << counter << endl;
		if(counter < 4)
			counter++;
	    }
#endif
    }

	return true;
}

M3DObject * M3DRegistrationOptimizer::createTargetObject(const Vector & v)
{
	if (problem == NULL)
		return NULL;

	return problem->createTargetObject(v);
}

void M3DRegistrationOptimizer::setOptimizerPosition(const Vector & start)
{
	if (method != NULL)
        delete method;

    method = new ConjugateGradientMethod(*wrapperProblem, start);

	((ConjugateGradientMethod *) method)->setBrentLinearSearchBoundFactor(tuningWt(BrentLinearSearchBoundFactor));
	((ConjugateGradientMethod *) method)->setBoundsFunction(boundsFunction);
	((ConjugateGradientMethod *) method)->setParameterTolerance(epsilon);

    lastBestVal = problem->evaluate(start);
    *(lastBestX[0]) = start;

	if (isAdaptive()) {
		lastPenalty = problem->getPenalty();

		if (lastBestObject != NULL)
		delete lastBestObject;
		lastBestObject = problem->createTargetObject(*(lastBestX[0]));
	}
	else {
	    lastPenalty = problem->computePenalty(start);
	    problem->computeTransformation(lastBestTransformation, *(lastBestX[0]));

	//    cout << "Start Value: " << lastBestVal << endl;

	    ((ConjugateGradientMethod *) method)->getFunction().setEvaluationCost(0);
	}
}

