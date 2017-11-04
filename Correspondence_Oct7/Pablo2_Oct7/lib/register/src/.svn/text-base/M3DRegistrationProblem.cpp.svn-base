#include <math.h>
#include "M3DRegistrationProblem.h"

//#define DISPLAY_PENALTIES


const double M3DRegistrationProblem::ENSEMBLE_OBJECTIVE_SCALE_FACTOR = 1.0;
const double M3DRegistrationProblem::ENSEMBLE_CONSTRAINT_PENALTY_SCALE_FACTOR = 0.02;
const double M3DRegistrationProblem::ENSEMBLE_MAX_RETURN = 1.7e308;

const double M3DRegistrationProblem::ENSEMBLE_PGA_SCALE_FACTOR = 200.0;
const double M3DRegistrationProblem::ENSEMBLE_TRANSLATION_FACTOR = 1.0;
const double M3DRegistrationProblem::ENSEMBLE_ROTATION_FACTOR = 8.0;
const double M3DRegistrationProblem::ENSEMBLE_SCALE_FACTOR = 12.0;

const double OBJECTIVE_SCALE_FACTOR = 1.0;



using namespace std;


M3DRegistrationProblem::M3DRegistrationProblem()
{
	match = NULL;
	penaltyWeight = 0.0;
	constraintsPenaltyWeight = 0.0;	
	lastPenalty = 0.0;
}

void M3DRegistrationProblem::initialize(Match * _match, double _penaltyWeight,
	double _constraintPenaltyWeight)
{
    match = _match;
    penaltyWeight = _penaltyWeight;
	constraintsPenaltyWeight = _constraintPenaltyWeight;
#ifdef DISPLAY_PENALTIES
    std::cout << "Registration penalties:  \n"<<std::endl;
#endif
}

M3DRegistrationProblem::M3DRegistrationProblem(Match * _match,
                                               double _penaltyWeight)
{
    match = _match;
    penaltyWeight = _penaltyWeight;
#ifdef DISPLAY_PENALTIES
    std::cout << "Registration (simularity) penalties:  \n";
#endif
}

double M3DRegistrationProblem::evaluate(const Vector &x)
{
    double matchVal;
    double penaltyVal;
    double totalPenalty;
    SimilarityTransform3D transform;

    evaluationCost += 1;

    if(match == NULL)
        return MAX_RETURN;

    if(x.size() != PARAMETER_SIZE) 
        return MAX_RETURN;
/*
    cout << x(0) << ", " << x(1) << ", " << x(2) << ", "
         << x(3) << ", " << x(4) << ", " << x(5) << ", "
         << x(6) << endl;
*/
    computeTransformation(transform, x);
#ifndef BINARY
    matchVal = match->computeModelMatch(transform);
#else
	matchVal = 0.0;
#endif
    penaltyVal = 0.0;

    if(penaltyWeight > 0.0)
        penaltyVal = match->computeModelPenalty(transform) * penaltyWeight;

    totalPenalty = OBJECTIVE_SCALE_FACTOR*(penaltyVal - matchVal);
#ifdef DISPLAY_PENALTIES
    std::cout << "   (" << penaltyVal << " model - " << matchVal << " match)*"
		<< OBJECTIVE_SCALE_FACTOR << " = " << totalPenalty << std::endl;
#endif
    return totalPenalty;
}

double M3DRegistrationProblem::computePenalty(const Vector &x)
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
}

void M3DRegistrationProblem::computeTransformation(
    SimilarityTransform3D & transform, const Vector &x)
{
    transform.setToIdentity();

    if(match == NULL)
        return;
    if(x.size() != PARAMETER_SIZE)
        return;

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


