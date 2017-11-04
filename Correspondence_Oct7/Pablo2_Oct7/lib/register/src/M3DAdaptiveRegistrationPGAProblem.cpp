#include <math.h>
#include <iostream>
#include "M3DAdaptiveRegistrationPGAProblem.h"

//#define DISPLAY_PENALTIES	/* Uncomment to print details of penalty calculations */
//#define DISPLAY_SHIFTS	/* Uncomment to print scale and translation changes */

//#define NO_MATCH_PENALTY	/* Uncomment to eliminate use of the match penalty */

// Comment out for d/r geometry penalty, uncomment for PGA Mahalanobis penalty
#define USE_MAHALANOBIS_PENALTY


using namespace std;


M3DAdaptiveRegistrationPGAProblem::M3DAdaptiveRegistrationPGAProblem(M3DPGAStats * pgaPtr, int _order)
	: M3DRegistrationProblem()
{
	pga = pgaPtr;
	order = _order;

	M3DObject * tempObj = pga->getMeanObjPtr()->assign();
	pga->applyMeanResidue(tempObj, order);
	delete tempObj;
}

int M3DAdaptiveRegistrationPGAProblem::parameterCount() const
{
	return pga->getPGDataPtr(order)->numPGs;
}

double M3DAdaptiveRegistrationPGAProblem::evaluate(const Vector &x)
{
    double matchVal;
    double penaltyVal;
	double totalPenalty;
    M3DObject * targetObject;
    int i;

    evaluationCost += 1;

    if (match == NULL)
        return ENSEMBLE_MAX_RETURN;

    targetObject = createTargetObject(x);

    // TOM: HACK -- turn off match if constraint penalty is all the way up
    bool constraintOnly = true;
    if (constraintsPenaltyWeight == ENSEMBLE_CONSTRAINT_PENALTY_SCALE_FACTOR) {
        matchVal = 0.0;
        penaltyVal = 0.0;
        constraintOnly = false;
    }
    else {
#ifdef NO_MATCH_PENALTY
        matchVal = 0.0;
#else
        matchVal = match->computeModelMatch(targetObject);
#endif

#ifdef USE_MAHALANOBIS_PENALTY
        if (penaltyWeight > 0.0) {
            penaltyVal = 0.0;
            for(i = 0; i < x.size(); i++)
                penaltyVal += x(i) * x(i);
            penaltyVal *= penaltyWeight;

        }
#else
        if (penaltyWeight > 0.0)
            penaltyVal = match->computeModelPenalty(targetObject) * penaltyWeight;
#endif
        else
            penaltyVal = 0.0;
    }

    lastPenalty = penaltyVal;
    totalPenalty = (lastPenalty - matchVal)*ENSEMBLE_OBJECTIVE_SCALE_FACTOR;
#ifdef DISPLAY_PENALTIES
    std::cout << "   " << penaltyVal << " model - " << matchVal << " match = " << totalPenalty
        << std::endl;
#endif

    if (targetObject != NULL)
        delete targetObject;

    return totalPenalty;
}

M3DObject * M3DAdaptiveRegistrationPGAProblem::createTargetObject(const Vector & x)
{
	// JJ: It does not matter what object is assigned to 'targetObject'.
	// PG is applied to mean object implicitly. (M3DPGAStats class)	
	M3DObject * targetObject = pga->getMeanObjPtr()->assign();

    if(targetObject == NULL || x.size() != parameterCount())
        return NULL;

	// The principal geodesic weights 
	std::vector<double> vals;
	for(int i = 0; i < x.size(); i++)
		vals.push_back(x(i) * ENSEMBLE_PGA_SCALE_FACTOR);

	// Apply PGA deformation
	pga->doPGADeform(targetObject, vals, order);

    // Apply initial transformation
	SimilarityTransform3D initialTransform = match->getInitialTransform();
	targetObject->select();

    Quat q = initialTransform.getRotation();
#ifdef UNFLIPPED
    Vector3D center(0.0, 1.0, 0.0);
    targetObject->rotateBy(q, center);
#else
    Vector3D center(0.0, 0.0, 0.0);
    targetObject->rotateBy(q, center);
#endif
    targetObject->scaleBy(initialTransform.getScale(), center);
    targetObject->translateBy(initialTransform.getTranslation());

    return targetObject;
}

M3DAdaptiveRegistrationPGAProblemBoundsFunction::M3DAdaptiveRegistrationPGAProblemBoundsFunction()
{
	bounds	= Vector(1, 3.0);
	factors	= Vector(1, M3DAdaptiveRegistrationPGAProblem::ENSEMBLE_PGA_SCALE_FACTOR );
}

