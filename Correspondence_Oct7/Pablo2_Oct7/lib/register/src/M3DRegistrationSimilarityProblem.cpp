#include <math.h>
#include "M3DRegistrationSimilarityProblem.h"

//#define DISPLAY_PENALTIES	/* Uncomment to print details of penalty calculations */
//#define DISPLAY_SHIFTS	/* Uncomment to print scale and translation changes */

//#define NO_MATCH_PENALTY	/* Uncomment to eliminate use of the match penalty */


M3DRegistrationSimilarityProblem::M3DRegistrationSimilarityProblem() :
    M3DRegistrationProblem()
{
}

double M3DRegistrationSimilarityProblem::evaluate(const Vector &x)
{
    double matchVal;
///    double penaltyVal;
    double totalPenalty;
    SimilarityTransform3D transform;

    evaluationCost += 1;

    if(match == NULL)
        return ENSEMBLE_MAX_RETURN;

///	if(x.size() != PARAMETER_SIZE) 
///		return MAX_RETURN;

    computeTransformation(transform, x);
#ifdef NO_MATCH_PENALTY
        matchVal = 0.0;
#else
        matchVal = match->computeModelMatch(targetObject);								// AGG: This code needs to be moved into Match
#endif
//    matchVal = match->computeModelMatch(transform);
	lastPenalty = 0.0;
/// penaltyVal = 0.0;

    if(penaltyWeight > 0.0)
        lastPenalty = match->computeModelPenalty(transform) * penaltyWeight * ENSEMBLE_OBJECTIVE_SCALE_FACTOR;

    totalPenalty = lastPenalty - ENSEMBLE_OBJECTIVE_SCALE_FACTOR*matchVal;
///    totalPenalty = OBJECTIVE_SCALE_FACTOR*(penaltyVal - matchVal);
#ifdef DISPLAY_PENALTIES
    std::cout << "   (" << lastPenalty << " model - " << matchVal << " match)*"
		<< OBJECTIVE_SCALE_FACTOR << " = " << totalPenalty << std::endl;
#endif
    return totalPenalty;
}

/*double M3DRegistrationSimilarityProblem::computePenalty(const Vector &x)
{
    SimilarityTransform3D transform;

    if(match == NULL)
        return MAX_RETURN;

    computeTransformation(transform, x);
    if(penaltyWeight > 0.0)
        return match->computeModelPenalty(transform) * penaltyWeight *
               OBJECTIVE_SCALE_FACTOR;
    else
        return 0.0;
}*/

void M3DRegistrationSimilarityProblem::computeTransformation(
    SimilarityTransform3D & transform, const Vector &x)
{
    transform.setToIdentity();

    if(match == NULL)
        return;
    //if(x.size() != PARAMETER_SIZE)
    //    return;

    M3DObject * object = match->getReferenceObject();
    if(object == NULL)
        return;

    // Setup rotation quaternion
    Quat xRot, yRot, zRot, finalRot;
    xRot.setAxisAngle(Vector3D(1.0, 0.0, 0.0), ENSEMBLE_ROTATION_FACTOR * x(3));
    yRot.setAxisAngle(Vector3D(0.0, 1.0, 0.0), ENSEMBLE_ROTATION_FACTOR * x(4));
    zRot.setAxisAngle(Vector3D(0.0, 0.0, 1.0), ENSEMBLE_ROTATION_FACTOR * x(5));

    // Make a composite rotation and rotate the object
    finalRot = zRot * yRot * xRot;
    finalRot.normalize();

    // Scale the object
    double scale;
	/*
    if(x(6) >= 0.0)
        scale = 1.0 + 8.0 * x(6);
    else
        scale = 1.0 / (1.0 - 8.0 * x(6));
	*/
//#define EUCLIDEAN_ONLY
#ifdef EUCLIDEAN_ONLY
	scale = 1.0;
#else
	scale = exp(ENSEMBLE_SCALE_FACTOR * x(6));
#endif

    // Do translation
    Vector3D trans(x(0), x(1), x(2));
	trans*=ENSEMBLE_TRANSLATION_FACTOR;

    transform.scaleBy(scale);
    transform.rotateBy(finalRot);
    transform.translateBy(trans);
}

M3DObject * M3DRegistrationSimilarityProblem::createTargetObject(const Vector &x)
{
	SimilarityTransform3D newTransform;
	computeTransformation(newTransform, x);

	M3DObject * targetObject = match->getReferenceObject()->assign();

	SimilarityTransform3D transform = *targetObject->getTransformation();
	transform.multiplyBy(newTransform);

	targetObject->select();
#ifdef UNFLIPPED
	targetObject->applySimilarity(transform, Vector3D(0.0, 1.0, 0.0));
#else
	targetObject->applySimilarity(transform, Vector3D(0.0, 0.0, 0.0));
#endif

	return targetObject;
}


M3DRegistrationSimilarityProblemBoundsFunction::M3DRegistrationSimilarityProblemBoundsFunction()
{
    bounds	= Vector(7, 1.0, 1.0, 1.0, R_PI, R_PI, R_PI, 1.0);
    //bounds	= Vector(8, 0.25, 0.25, 0.25, 0.5*R_PI, 0.5*R_PI, 0.5*R_PI, 1.0, 1.0);
	factors	= Vector(8, M3DRegistrationSimilarityProblem::ENSEMBLE_TRANSLATION_FACTOR,
        M3DRegistrationSimilarityProblem::ENSEMBLE_TRANSLATION_FACTOR,
		M3DRegistrationSimilarityProblem::ENSEMBLE_TRANSLATION_FACTOR,
		M3DRegistrationSimilarityProblem::ENSEMBLE_ROTATION_FACTOR,
		M3DRegistrationSimilarityProblem::ENSEMBLE_ROTATION_FACTOR,
		M3DRegistrationSimilarityProblem::ENSEMBLE_ROTATION_FACTOR,
		M3DRegistrationSimilarityProblem::ENSEMBLE_SCALE_FACTOR);
}

