#include "M3DSpokeProblem.h"
#include "M3DRegistrationProblem.h"
#include "Tuning.h"
#include <iomanip>

// Enable only if you want an insane amount of intermediate debugging models generated
// (like tens of thousands)
//#define DEBUG_SAVE_INTERMEDIATE_MODELS

#ifdef DEBUG_SAVE_INTERMEDIATE_MODELS
#include "P3DControl.h"
#endif

using namespace std;
/* BINARY PABLO

  --- OBJECTIVE FUNCTION SCALE FACTORS ---

  The function scale factors were changed from const's to tuning().weight's so
  that they may be set at run-time.

  OLD COMMENTS:
  These scale values got tweaked slightly to work with the new objective function.

  There is a matlab script called CompareParams that takes output from the
  #define OUTPUT_TEST in the Optimizer and plots the parameters next to each
  other.  The goal is to get all the parameters to have more-or-less the same
  curvature at the bottom of the trough.

  To widen a parameter, reduce the constant value
  To narrow a parameter, increase the constant value

  I left the old values in, but commented out for comparison.  I lowered the scale
  (radius) by a factor of 10!, so that might indicate that radius is not scaled
  correctly regular Pablo.

  Also note that the full objective function is evaluated, so it takes into account
  things like slider settings.  That means that the scales that work for one slider
  might not be correct for another setting.  In Binary Pablo, it's clear that the
  regularizing penalty is much better behaved near the optimum than the distance
  match, for instance.
*/

#ifdef BINARY
extern int globalVerbosity;
#endif


M3DSpokeProblem::M3DSpokeProblem(int appID)
{
    referenceObject = NULL;
    match = NULL;
    targetObject = NULL;

    figureId = -1;
	atomId = -1;

    evaluationCost = 0;

	applicationID = appID;
}

M3DSpokeProblem::M3DSpokeProblem(Match * _match,
        M3DObject * _referenceObject, M3DObject * _targetObject,
	 int _figureId, int _atomId, int _spokeId, int appID)
{
    referenceObject = _referenceObject;
    match = _match;
    targetObject = _targetObject;

    figureId	= _figureId;
	atomId		= _atomId;
	spokeId		= _spokeId;

    evaluationCost	= 0;
	applicationID	= appID;
}

M3DSpokeProblem::~M3DSpokeProblem()
{
}


// Used only to determine the image match value before optimizing (preview operation)
double M3DSpokeProblem::evaluateImageMatch(const Vector & x)
{
    evaluationCost++;

    if (match == NULL || targetObject == NULL || referenceObject == NULL)
        return MAX_RETURN;

	//
	// Restore primitive only if we have a spoke id
	// that is not -1.
	//
	if( spokeId != -1 ) {
		M3DTubePrimitive* currPrimPtr	=
			dynamic_cast<M3DTubePrimitive*>(targetObject->getFigurePtr(figureId)->getPrimitivePtr(atomId));
		currPrimPtr->setRN( spokeId,
			dynamic_cast<const M3DTubePrimitive*>(referenceObject->getFigurePtr(figureId)->getPrimitivePtr(atomId))->getRN(spokeId)
		);
	}

#ifdef BINARY
	double geomPenalty;
	return match->computeSpokeMatchAndPenalties(targetObject, figureId, -1, x, geomPenalty);
#else
    return match->computeMainFigureMatch(targetObject, figureId);
#endif
}

double M3DSpokeProblem::evaluate(const Vector & x)
{
    evaluationCost++;

    if(match == NULL || targetObject == NULL || referenceObject == NULL)
        return MAX_RETURN;

	//
	// Restore primitive and use the x vector only if we have a spoke id that is not -1,
	// else silently ignore.
	//
	if( spokeId != -1 ) {
		M3DTubePrimitive* currPrimPtr	=
			dynamic_cast<M3DTubePrimitive*>(targetObject->getFigurePtr(figureId)->getPrimitivePtr(atomId));
		currPrimPtr->setRN( spokeId,
			dynamic_cast<const M3DTubePrimitive*>(referenceObject->getFigurePtr(figureId)->getPrimitivePtr(atomId))->getRN(spokeId)
		);
		applyTransform(*currPrimPtr, x(0));
	}

#ifdef DEBUG_SAVE_INTERMEDIATE_MODELS
	char tempName[50];
    static int counterAtomIterations=0;
	static P3DControl * gControl	= NULL;
	if( gControl == NULL ) {
		gControl	= new P3DControl(0);
	}

	M3DObject * temp = targetObject->assign();
	gControl->newModel(temp);
	sprintf(tempName,"./cg_atom%03d_spoke%02d_iter%05d.m3d", atomId, spokeId, counterAtomIterations++);
	printf("SAVE: %s\n", tempName );
	gControl->write(tempName);
#endif

	return match->computeSpokeMatchAndPenalties(targetObject, figureId, atomId, x,
		geomPenalty);
}

void M3DSpokeProblem::applyTransform(M3DPrimitive & prim, const double x)
{
	M3DTubePrimitive& tubePrim	= dynamic_cast<M3DTubePrimitive&>(prim);
	assert( spokeId >= 0 && spokeId < tubePrim.getNumberOfSpokes() );
	const double newR =
		(tubePrim.getRN(spokeId) - fabs(tubePrim.getD())) * exp(tuningWt(SpokeScaleFactor) * x)
		+ fabs(tubePrim.getD());
	if( globalVerbosity >= 2 ) {
		cout << setprecision(4)
			<< "spoke( " << setw(2) << spokeId << "): "
			<< setw(6) << tubePrim.getRN(spokeId)
			<< " - "
			<< setw(6) << newR
			<< ", x = " << x << endl;
	}
	tubePrim.setRN(spokeId, newR );
}


M3DSpokeProblemBoundsFunction::M3DSpokeProblemBoundsFunction()
{
	bounds	= Vector(1, log(tuningWt(SpokeScaleBound)));
	factors	= Vector(1, tuningWt(SpokeScaleFactor) );
}

