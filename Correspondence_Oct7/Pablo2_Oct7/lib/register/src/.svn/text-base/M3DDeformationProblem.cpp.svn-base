#include "M3DDeformationProblem.h"
#include "M3DRegistrationProblem.h"
#include "Tuning.h"
#include "M3DAtomPredictorQuad.h"
#include "P3DControl.h"
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

M3DDeformationProblem::M3DDeformationProblem(int appID)
{
	referenceObject = NULL;
	match = NULL;
	targetObject = NULL;
	primPga = NULL;

	figureId = -1;
	atomId = -1;

	evaluationCost = 0;

	applicationID = appID;
}

M3DDeformationProblem::M3DDeformationProblem(Match * _match,
		M3DObject * _referenceObject, M3DObject * _targetObject,
		int _figureId, int _atomId,
		M3DPGAPrimitiveStats * _primPga,int appID)
{ 
	//here most of the variables are using references from the optimizer class , not copy it over, so no need to delete them~
	referenceObject = _referenceObject;
	match = _match;
	targetObject = _targetObject;
	primPga = _primPga;

	figureId = _figureId;
	atomId = _atomId;

	evaluationCost = 0;
	applicationID = appID;
	figurePredictor = NULL;

}

M3DDeformationProblem::~M3DDeformationProblem()
{	
	if (figurePredictor !=NULL)
		delete figurePredictor ;	
}

void M3DDeformationProblem::setPGA(M3DPGAPrimitiveStats * primPgaStats)
{
	if (primPgaStats != NULL)
		primPga = primPgaStats;
}

// Used only to determine the image match value before optimizing (preview operation)
double M3DDeformationProblem::evaluateImageMatch(const Vector & x, M3DObject * _referenceObject)
{
		evaluationCost++;

	if (match == NULL || targetObject == NULL || referenceObject == NULL)
		return MAX_RETURN;


	if (_referenceObject !=NULL)
		referenceObject = _referenceObject->assign();


	M3DObject * candidateObject  = referenceObject->assign();
	// Apply initial transformation
	Vector3D zero(0.0, 0.0, 0.0);

	SimilarityTransform3D * xform = candidateObject->getTransformation();
	if (xform!= NULL)
	candidateObject->getFigurePtr(figureId)->applySimilarity( *xform, zero );
	

	// if an end-primitive, elongation is copied by the appropriate assignment operator
#ifdef BINARY
	double geomPenalty;
	return match->computeAtomMatchAndPenalties(candidateObject, figureId, atomId, x, geomPenalty);
#else
	 return match->computeAtomMatchAndPenalties(candidateObject, figureId,atomId, x,geomPenalty);
#endif

}

double M3DDeformationProblem::evaluate(const Vector & x)
{
	evaluationCost++;

	if(match == NULL || targetObject == NULL || referenceObject == NULL)
		return MAX_RETURN;

	// Replace a target atom with the corresponding reference atom - SIMPLIFY! //!Was wrong before!

    createTargetPrimitive(x);

	// Apply initial transformation
	Vector3D zero(0.0, 0.0, 0.0);

	SimilarityTransform3D * xform = targetObject->getTransformation();

	M3DObject * candidateObject = targetObject->assign(); //used to evalute the image match

	if (xform!= NULL && primPga!= NULL )
		candidateObject->getFigurePtr(figureId)->applySimilarity( *xform, zero );


#ifdef DEBUG_2
	    cout << " currentX=";
			for (int i = 0; i < x.size(); i++)
				cout << ' ' << x(i);
		cout << '\n';			
	
		//save the prediction to see whether it's right
		char tempName[50];
		static int counterAtomIterations=0;
		P3DControl * gControl= new P3DControl(0);	
		M3DObject * saveTemp = candidateObject->assign();	
		gControl->newModel(saveTemp);
		sprintf(tempName,"./iter%d_atom%d.m3d",counterAtomIterations++,atomId);
		gControl->write(tempName);	
		delete saveTemp;
	
#endif	

	double matchVal = match->computeAtomMatchAndPenalties(candidateObject, figureId,atomId, x,geomPenalty);

	delete candidateObject;

	return matchVal;
}



void M3DDeformationProblem ::createTargetPrimitive(const Vector & x)
{// apply the final X to get the best result in this iteration

	M3DFigure * referenceFigure = referenceObject->getFigurePtr(figureId);
	M3DFigure * targetFigure = targetObject->getFigurePtr(figureId);	
	
	M3DPrimitive  *  referencePrimitive = referenceFigure->getPrimitivePtr(atomId);	
    M3DPrimitive  *	 targetPrimitive  = targetFigure->getPrimitivePtr(atomId);
	
    (*targetPrimitive) = (*referencePrimitive);
	
	if (primPga != NULL) {	    
		
		doApplyVectorPGA(targetFigure,x); //need to pass in the figure to do the neighbor alignment
	}
	else{
		applyVector(*targetPrimitive, x);
	}

}




Components M3DDeformationProblem::evaluateTerms(const Vector & x)
{
    evaluationCost++;

	if(match == NULL || targetObject == NULL || referenceObject == NULL)
	 return Components(MAX_RETURN);

	// Replace a target atom with the corresponding reference atom - SIMPLIFY! //!Was wrong before!

   createTargetPrimitive(x);

   	M3DObject * candidateObject = targetObject->assign(); //used to evalute the image ma

	// Apply initial transformation
	Vector3D zero(0.0, 0.0, 0.0);

	SimilarityTransform3D * xform = targetObject->getTransformation();

	if (xform!= NULL && primPga!= NULL )
		candidateObject->getFigurePtr(figureId)->applySimilarity( *xform, zero );
//#ifdef DEBUG_2
					
		//save the prediction to see whether it's right
		//char tempName[50];
	
		//P3DControl * gControl= new P3DControl(0);	
		//M3DObject * saveTemp = candidateObject->assign();	
		//gControl->newModel(saveTemp);
		//sprintf(tempName,"debugbpointobjrelsampling.m3d");
		//gControl->write(tempName);	
		//delete saveTemp;
	
//#endif	

	Components comp =  match->computeAtomMatchAndPenaltiesAsComponents(candidateObject, figureId, atomId, x, geomPenalty);


	return comp;
}



void M3DDeformationProblem::doApplyVectorPGA(M3DFigure * currFigPtr, const Vector & x){

	int i,count = 0;
	for ( i = 0; i < x.size(); i++){
		if (x(i) == 0){
			count++;
		}
	}	

	if(count == x.size())//if x = (0,0,0,0)
		return;	

	M3DPrimitive * currPrimPtr = currFigPtr->getPrimitivePtr(atomId);
	M3DPrimitive * readyToDeformAtom = (primPga->getDeltaMeanPrimitivePtr(figureId,atomId))->copyPtr(); 

	applyVectorPGA(*readyToDeformAtom, x);

	M3DFigure * alignReferenceFigurePtr = primPga->getAlignReferenceObj()->getFigurePtr(figureId);
	M3DDeformationOptimizer::addDiffToPrimitive(currFigPtr,atomId, alignReferenceFigurePtr,readyToDeformAtom);

    //
	// Correct any consistency issues.
	//
	currFigPtr->fixGlobalConsistency();

	delete readyToDeformAtom;

}


void M3DDeformationProblem::applyVectorPGA(M3DPrimitive & prim, const Vector & x)
{
	int vecSize = x.size();
	int isZeroVector = 0;
	int j;
	for ( j= 0; j < vecSize; j++)
		if (x(j) ==0) 
			isZeroVector++;

	if (isZeroVector == vecSize)
		return;//all coefficients are zero, return;

	M3DPrimitive * newPrim = prim.copyPtr();

	Vector vals = x * ((int) tuningWt(AtomPGAFactor));

	primPga->doPGADeform(newPrim, vals, figureId, atomId, &prim);	//Deform newPrim and replace prim

	delete newPrim;
}

void M3DDeformationProblem::applyVector(M3DPrimitive & prim, const Vector & x)
{
	Vector3D b = prim.getB();
	Vector3D n = prim.getN();
	Vector3D bPerp = prim.getBPerp();

	Vector3D trans;
	Quat finalRot;
	Vector3D axis;
	double angle;
	double scale, eta;
	
	const bool axesAligned	= (0 != (int(tuningWt(AtomRestrictMovement)) & 0x0008));

	trans =  (x(0) * b);
	trans += (x(1) * n);
	trans += (x(2) * bPerp);
	trans *= tuningWt(AtomTranslationFactor);
	if( int(tuningWt(AtomRestrictMovement)) & 0x0001) {
		static bool reported	= false;
		if( !reported ) {
			cout << "Disabling atom movement along tangent/bisector.\n";
			reported = true;
		}
		if( axesAligned) {
			trans.setX(0.0);
		}
		else {
			// restrict translation along bisector to 0.
			trans	-= x(0) * b;
		}
	}
	if( int(tuningWt(AtomRestrictMovement)) & 0x0002) {
		static bool reported	= false;
		if( !reported ) {
			cout << "Disabling atom movement along normal.\n";
			reported = true;
		}
		if( axesAligned) {
			trans.setY(0.0);
		}
		else {
			// restrict translation along normal to 0.
			trans	-= x(1) * n;
		}
	}
	if( int(tuningWt(AtomRestrictMovement)) & 0x0004) {
		static bool reported	= false;
		if( !reported ) {
			cout << "Disabling atom movement along binormal.\n";
			reported = true;
		}
		if( axesAligned) {
			trans.setZ(0.0);
		}
		else {
			// restrict translation along binormal to 0.
			trans	-= x(2) * bPerp;
		}
	}
	prim.translateBy(trans);

	axis.set(x(3), x(4), x(5));
	angle = tuningWt(AtomRotationFactor) * axis.normalize();
	finalRot.setAxisAngle(axis, angle);
	prim.rotateBy(finalRot);

	// Scale the primitive
	scale = exp(tuningWt(AtomScaleFactor) * x(6));
	if( int(tuningWt(AtomRestrictMovement)) & 0x0010) {
		static bool reported	= false;
		if( !reported ) {
			cout << "Coupling r and theta changes.\n";
			reported = true;
		}
		// scale will cause a corresponding change in theta and theta is ignored.
		prim.setR( prim.getR() * scale);
		prim.setTheta( acos(max(min( cos(prim.getTheta())/scale, 1.0),-1.0)) );
	} else {
		prim.scaleBy(scale);
		// Change the primitive's theta
		prim.setTheta(prim.getTheta() + tuningWt(AtomThetaFactor) * x(7));
	}

	if ((tuningWt(AtomEtaFactor) != 0.0) && (prim.type() == M3D_END_PRIMITIVE)) {
		// Change the primitive's eta
		eta = (dynamic_cast<M3DEndPrimitive *>(&prim))->getElongation() * exp(tuningWt(AtomEtaFactor) * x(8));
		(dynamic_cast<M3DEndPrimitive *>(&prim))->setElongation(eta);
	}

	//
	// Correct any consistency issues.
	//
	M3DFigure* currFigPtr = targetObject->getFigurePtr(figureId);
	currFigPtr->fixGlobalConsistency();
}



double M3DDeformationProblemBoundsFunction::bound( int i )
{
	assert( i >= 0 );
	if( problem->primPga != NULL ) {
		return 3.0;
	}
	else {
		Vector bounds(9,
				tuningWt(AtomTranslationBound), tuningWt(AtomTranslationBound), tuningWt(AtomTranslationBound),
				tuningWt(AtomRotationBound)*R_PI/180.0, tuningWt(AtomRotationBound)*R_PI/180.0, tuningWt(AtomRotationBound)*R_PI/180.0,
				log(tuningWt(AtomScaleBound)),
				tuningWt(AtomThetaBound) * R_PI/180.0,
				log(tuningWt(AtomEtaBound)));
		return bounds( (i < bounds.size()) ? i : (bounds.size() - 1) );
	}
}

double M3DDeformationProblemBoundsFunction::factor( int i )
{
	assert( i >= 0 );
	if( problem->primPga != NULL ) {
		return tuningWt(AtomPGAFactor);
	}
	else {
		Vector factors(9, 
				tuningWt(AtomTranslationFactor), tuningWt(AtomTranslationFactor),
				tuningWt(AtomTranslationFactor), tuningWt(AtomRotationFactor),
				tuningWt(AtomRotationFactor), tuningWt(AtomRotationFactor),
				tuningWt(AtomScaleFactor), tuningWt(AtomThetaFactor),
				tuningWt(AtomEtaFactor));
		return factors( (i < factors.size()) ? i : (factors.size() - 1) );
	}
}

