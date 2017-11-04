#include <math.h>
#include <iostream>
#include "M3DRegistrationPGAProblem.h"

//#define DISPLAY_PENALTIES	/* Uncomment to print details of penalty calculations */
//#define DISPLAY_SHIFTS	/* Uncomment to print scale and translation changes */

//#define NO_MATCH_PENALTY	/* Uncomment to eliminate use of the match penalty */

//#define CHECK_SIMILARITY_TRANSFORM	/*JJ : tp print out the similarity transform */

// Comment out for d/r geometry penalty, uncomment for PGA Mahalanobis penalty
#define USE_MAHALANOBIS_PENALTY


using namespace std;


M3DRegistrationPGAProblem::M3DRegistrationPGAProblem() :
    M3DRegistrationProblem()
{
}

int M3DRegistrationPGAProblem::parameterCount() const
{
	return pga->getPGDataPtr(order)->numPGs + 7;
}

double M3DRegistrationPGAProblem::evaluate(const Vector &x)
{
    double matchVal;
    double penaltyVal;
    double constraintPenalty;
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
        matchVal = match->computeModelMatch(targetObject);								// AGG: This code needs to be moved into Match
#endif

#ifdef USE_MAHALANOBIS_PENALTY
        if (penaltyWeight > 0.0) {
            penaltyVal = 0.0;

            for(i = 7; i < x.size(); i++)
                penaltyVal += x(i) * x(i);
            penaltyVal *= penaltyWeight;

			//double mPen = match->computeModelPenalty(targetObject);
			//penaltyVal += mPen * penaltyWeight;
        }
#else
        if (penaltyWeight > 0.0)
            penaltyVal = match->computeModelPenalty(targetObject) * penaltyWeight;
#endif
        else
            penaltyVal = 0.0;
    }

    /*if(constraintsPenaltyWeight > 0.0 && targetObject != NULL)
    {
        constraintPenalty = 
            match->computeConstraintsPenalty(*targetObject, figureId, constraintOnly) * 
            constraintsPenaltyWeight;
    }
    else*/
        constraintPenalty = 0.0;

    lastPenalty = constraintPenalty + penaltyVal;
    totalPenalty = (lastPenalty - matchVal)*ENSEMBLE_OBJECTIVE_SCALE_FACTOR;
#ifdef DISPLAY_PENALTIES
    std::cout << "   " << penaltyVal << "(penW :" <<penaltyWeight<<") model + " << constraintPenalty
        <<  " constraint - " << matchVal << " match = " << totalPenalty
        << std::endl;
#endif

    if (targetObject != NULL)
        delete targetObject;

    return totalPenalty;
}


M3DObject * M3DRegistrationPGAProblem::createTargetObject(const Vector & x)
{
    //M3DObject * targetObject = match->getReferenceObject()->assign();
	//Get start object from pga.
	//M3DObject * targetObject = pga->getOriginalMean()->assign();
	M3DObject * targetObject = pga->getMeanObjPtr()->assign();

    if(targetObject == NULL || x.size() != parameterCount())
        return NULL;

	// The principal geodesic weights 
	std::vector<double> vals;
	for(int i = 7; i < x.size(); i++)
		vals.push_back(x(i) * ENSEMBLE_PGA_SCALE_FACTOR);

	// Apply PGA deformation
	////pga->pgaDeform(targetObject, vals);
	if (! pga->doPGADeform(targetObject, vals, order)) {
		std::cout << "Error in PGA deformation at M3DRegistrationPGAProblem" << endl;
		return NULL;
	}

    // Apply initial transformation
	SimilarityTransform3D initialTransform = *targetObject->getTransformation();
	targetObject->select();

#ifdef UNFLIPPED
	Vector3D center(0.0, 1.0, 0.0);
#else
	Vector3D center(0.0, 0.0, 0.0);
#endif
    Quat q = initialTransform.getRotation();
    targetObject->rotateBy(q, center);
    targetObject->scaleBy(initialTransform.getScale(), center);
    targetObject->translateBy(initialTransform.getTranslation());

	// Apply hand placement.
	SimilarityTransform3D handTransform = pga->getHandPlacementTransform();
	q = handTransform.getRotation();
	targetObject->rotateBy(q, center);
    targetObject->scaleBy(handTransform.getScale(), center);
    targetObject->translateBy(handTransform.getTranslation());

	//JJ : checking initialTransform (stored in m3d file)
#ifdef CHECK_SIMILARITY_TRANSFORM
	cout<<"InitialTransform in M3DRegistrationPGAProblem::createTargetObject()"<<endl;

	initialTransform.getTranslation().print();
	initialTransform.getRotation().print();
	cout<<"scale "<<initialTransform.getScale()<<endl;

	cout<<"handTransformation in  M3DRegistrationPGAProblem::createTargetObject()"<<endl;

	handTransform.getTranslation().print();
	handTransform.getRotation().print();
	cout<<"scale "<<handTransform.getScale()<<endl;
#endif

    // Apply similarity transform
    Vector3D trans;
    double angle;
    Vector3D axis;
    Quat finalRot;
    double scale;

    trans.set(x(0), x(1), x(2));
    trans *= ENSEMBLE_TRANSLATION_FACTOR;

    axis.set(x(3), x(4), x(5));
    angle = ENSEMBLE_ROTATION_FACTOR * axis.normalize();
    finalRot.setAxisAngle(axis, angle);

	// Hack: turn off scale for rigid transformation
	//scale = 1.0;
	scale = exp(ENSEMBLE_SCALE_FACTOR * x(6));

    targetObject->rotateBy(finalRot, center);
    targetObject->scaleBy(scale, center);
    targetObject->translateBy(trans);

	// Set best similarity transform to the current one
	pga->setBestTransform(SimilarityTransform3D(trans, finalRot, scale));

    return targetObject;
}


M3DRegistrationPGAProblemBoundsFunction::M3DRegistrationPGAProblemBoundsFunction()
{
    bounds	= Vector(8, 1.0, 1.0, 1.0, R_PI, R_PI, R_PI, 1.0, 3.0);
	factors	= Vector(8, M3DRegistrationPGAProblem::ENSEMBLE_TRANSLATION_FACTOR,
        M3DRegistrationPGAProblem::ENSEMBLE_TRANSLATION_FACTOR,
        M3DRegistrationPGAProblem::ENSEMBLE_TRANSLATION_FACTOR,
        M3DRegistrationPGAProblem::ENSEMBLE_ROTATION_FACTOR,
        M3DRegistrationPGAProblem::ENSEMBLE_ROTATION_FACTOR,
        M3DRegistrationPGAProblem::ENSEMBLE_ROTATION_FACTOR,
        M3DRegistrationPGAProblem::ENSEMBLE_SCALE_FACTOR,
        M3DRegistrationPGAProblem::ENSEMBLE_PGA_SCALE_FACTOR);
}

