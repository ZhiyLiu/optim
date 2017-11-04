// dibyendu

#include "M3DSRepProblem.h"
#include "M3DRegistrationProblem.h"
#include "Tuning.h"
#include "M3DAtomPredictorQuad.h"
#include "P3DControl.h"
using namespace std;


#ifdef BINARY
extern int globalVerbosity;
#endif

M3DSRepProblem::M3DSRepProblem(int appID)
{
	referenceObject = NULL;
	match = NULL;
	targetObject = NULL;

	figureId = -1 ;
	atomId = -1 ;
	spokeId = -1 ;

	evaluationCost = 0;

	applicationID = appID;
}

M3DSRepProblem::M3DSRepProblem(Match * _match,
		M3DObject * _referenceObject, M3DObject * _targetObject,
		int _figureId, int _atomId, int _spokeId, int appID)
{ 
	//here most of the variables are using references from the optimizer class , not copy it over, so no need to delete them~
	referenceObject = _referenceObject;
	match = _match;
	targetObject = _targetObject;

	figureId = _figureId ;
	atomId = _atomId ;
	spokeId = _spokeId ;

	evaluationCost = 0;
	applicationID = appID;
	figurePredictor = NULL;

}

M3DSRepProblem::M3DSRepProblem(Match * _match,
		M3DObject * _referenceObject, M3DObject * _targetObject,
		int _figureId, int _atomId, int appID)
{ 
	//here most of the variables are using references from the optimizer class , not copy it over, so no need to delete them~
	referenceObject = _referenceObject;
	match = _match;
	targetObject = _targetObject;

	figureId = _figureId ;
	atomId = _atomId ;	
	spokeId = -1 ;

	evaluationCost = 0;
	applicationID = appID;
	figurePredictor = NULL;

}


M3DSRepProblem::~M3DSRepProblem()
{	
	if (figurePredictor !=NULL)
		delete figurePredictor ;	
}


// Used only to determine the image match value before optimizing (preview operation)
double M3DSRepProblem::evaluateImageMatch(const Vector & x, M3DObject * _referenceObject)
{
		evaluationCost++;

	if (match == NULL || targetObject == NULL || referenceObject == NULL)
		return MAX_RETURN;


	if (_referenceObject !=NULL)
		referenceObject = _referenceObject->assign();


	M3DObject * candidateObject  = referenceObject->assign();
	

	// if an end-primitive, elongation is copied by the appropriate assignment operator
#ifdef BINARY
	double geomPenalty;
	return match->computeSRepMatchAndPenalties(candidateObject, figureId, atomId, spokeId, x, geomPenalty);
#else
	 return match->computeAtomMatchAndPenalties(candidateObject, figureId,atomId, x,geomPenalty);
#endif

}

double M3DSRepProblem::evaluate(const Vector & x)
{
	evaluationCost++;

	if(match == NULL || targetObject == NULL || referenceObject == NULL)
		return MAX_RETURN;

	// Replace a target atom with the corresponding reference atom - SIMPLIFY! //!Was wrong before!

    createTargetPrimitive(x);

	M3DObject * candidateObject = targetObject->assign(); //used to evalute the image match

	double matchVal = match->computeSRepMatchAndPenalties(candidateObject, figureId, atomId, spokeId, x, geomPenalty);

	delete candidateObject;

	return matchVal;
}



Components M3DSRepProblem::evaluateTerms(const Vector & x)
{
    evaluationCost++;

	if(match == NULL || targetObject == NULL || referenceObject == NULL)
	 return Components(MAX_RETURN);

	// Replace a target atom with the corresponding reference atom - SIMPLIFY! //!Was wrong before!

   createTargetPrimitive(x);

   	M3DObject * candidateObject = targetObject->assign(); //used to evalute the image match

	Components comp =  match->computeSRepMatchAndPenaltiesAsComponents(candidateObject, figureId, atomId, spokeId, x, geomPenalty);

	return comp;
}

void M3DSRepProblem ::createTargetPrimitive(const Vector & x)
{
	// apply the transformation X to get the result in this iteration

	// this function only changes one particular spoke (spokeId) of the atom

	// for applying an all-spoke transformation to the atom, use updateTargetPrimitive()

	M3DFigure * referenceFigure = referenceObject->getFigurePtr(figureId);
	M3DFigure * targetFigure = targetObject->getFigurePtr(figureId);	
	
	M3DPrimitive  *  referencePrimitive = referenceFigure->getPrimitivePtr(atomId);	
    M3DPrimitive  *	 targetPrimitive  = targetFigure->getPrimitivePtr(atomId);
	
    (*targetPrimitive) = (*referencePrimitive);

	applyVector(*targetPrimitive, x);

}


void M3DSRepProblem::applyVector(M3DPrimitive & prim, const Vector & x)
{
	/*

	At the s-rep stage, all spokes are optimized separately, so there are 3 parameters (length, orientation[2]) per spoke
	
	del_Spoke = ( log(del_r), del_theta, del_phi )  // check function Vector3D::rotateBy( delTheta, delPhi )

	*/	

	double scale	= exp(tuningWt(SRepScaleFactor) * x(0));
	double delTheta	= tuningWt(SRepThetaFactor) * x(1) ;
	double delPhi	= tuningWt(SRepThetaFactor) * x(2) ;

	// scale and rotate the spoke 

	prim.scaleSpokeBy( spokeId, scale ) ;
	prim.rotateSpokeBy( spokeId, delTheta, delPhi ) ;


	// Correct any consistency issues.	
	// M3DFigure* currFigPtr = targetObject->getFigurePtr(figureId);
	// currFigPtr->fixGlobalConsistency();
}

//#define DEBUG_DIBYENDU_2
#undef DEBUG_DIBYENDU2

#ifndef DEBUG_DIBYENDU_2

void M3DSRepProblem ::updateTargetPrimitive(const Vector & x)
{
	// apply the final X to get the best result in this iteration

	// this transformation is for all the 2 or 3 spokes and not just one spoke as in createPrimitive()

	M3DFigure * referenceFigure = referenceObject->getFigurePtr(figureId);
	M3DFigure * targetFigure = targetObject->getFigurePtr(figureId);	
	
	M3DPrimitive  *  referencePrimitive = referenceFigure->getPrimitivePtr(atomId);	
    M3DPrimitive  *	 targetPrimitive  = targetFigure->getPrimitivePtr(atomId);
	
    (*targetPrimitive) = (*referencePrimitive);

	applyVectorAtom(*targetPrimitive, x);

}


void M3DSRepProblem::applyVectorAtom(M3DPrimitive & prim, const Vector & x)
{
	/*

	At the s-rep stage, all spokes are optimized separately, so there are 3 parameters (length, orientation[2]) per spoke
	
	del_Spoke = ( log(del_r), del_theta, del_phi )  // check function Vector3D::rotateBy( delTheta, delPhi )

	The vector to be applied is,
	x = { del_Position, del_Spoke_0, del_Spoke_1, del_Spoke_End } = 3 + 3 + 3 + 3(end) = 9 or 12 entries

	del_Position = ( del_P_X, del_P_Y, del_P_Z )
	
	*/

	Vector3D trans ;

	double scale0, scale1, scaleEnd ;
	double delTheta0, delTheta1, delThetaEnd ;
	double delPhi0, delPhi1, delPhiEnd ;


	trans =  ( x(0) * Vector3D(1.0, 0.0, 0.0) ) ; 
	trans += ( x(1) * Vector3D(0.0, 1.0, 0.0) ) ;
	trans += ( x(2) * Vector3D(0.0, 0.0, 1.0) ) ;

	trans *= tuningWt(AtomTranslationFactor);

	// If atom movement is restricted, then the atom hub does not move
	if( int(tuningWt(SRepAtomRestrictMovement)) ) {

		static bool reported	= false;

		if( !reported ) {
			cout << "Disabling atom hub movement in s-rep stage" << endl ;
			reported = 1 ;
		}

		trans.setX( 0.0 ) ; trans.setY( 0.0 ) ; trans.setZ( 0.0 ) ;
	}

	prim.translateBy(trans);	
	
	scale0		= exp(tuningWt(SRepScaleFactor) * x(3));
	delTheta0	= tuningWt(SRepThetaFactor) * x(4) ;
	delPhi0		= tuningWt(SRepThetaFactor) * x(5) ;

	scale1		= exp(tuningWt(SRepScaleFactor) * x(6));
	delTheta1	= tuningWt(SRepThetaFactor) * x(7) ;
	delPhi1		= tuningWt(SRepThetaFactor) * x(8) ;

	// scale the radii

	prim.scaleSpokeBy( 0, scale0 ) ;
	prim.rotateSpokeBy( 0, delTheta0, delPhi0 ) ;

	prim.scaleSpokeBy( 1, scale1 ) ;
	prim.rotateSpokeBy( 1, delTheta1, delPhi1 ) ;

	// deal with the end spoke

	if( prim.type() == M3D_END_PRIMITIVE ) {		

		scaleEnd	= exp(tuningWt(SRepScaleFactor) * x(9));
		delThetaEnd = tuningWt(SRepThetaFactor) * x(10) ;
		delPhiEnd	= tuningWt(SRepThetaFactor) * x(11) ;		

		prim.scaleSpokeBy( 2, scaleEnd ) ;
		prim.rotateSpokeBy( 2, delThetaEnd, delPhiEnd ) ;

	}

	//
	// Correct any consistency issues.
	//
	// M3DFigure* currFigPtr = targetObject->getFigurePtr(figureId);
	// currFigPtr->fixGlobalConsistency();
}

#endif 

#ifdef DEBUG_DIBYENDU_2

void M3DSRepProblem ::updateTargetPrimitive(const Vector & x)
{
	// apply the final X to get the best result in this iteration

	// this transformation is for all the 2 or 3 spokes and not just one spoke as in createPrimitive()

	M3DFigure * referenceFigure = referenceObject->getFigurePtr(figureId);
	M3DFigure * targetFigure = targetObject->getFigurePtr(figureId);	
	
	M3DPrimitive  *  referencePrimitive = referenceFigure->getPrimitivePtr(atomId);	
    M3DPrimitive  *	 targetPrimitive  = targetFigure->getPrimitivePtr(atomId);
	
    (*targetPrimitive) = (*referencePrimitive);

	applyVectorAtom(*targetPrimitive, x);

	cout << "Atom ID = " << atomId << endl ;

	cout << endl << "Reference primitive was: " << endl << endl ;
	referencePrimitive->print() ;

	cout << endl << "Target primitive is: " << endl << endl ;
	targetPrimitive->print() ;

}

void M3DSRepProblem::applyVectorAtom(M3DPrimitive & prim, const Vector & x)
{
	/*

	At the s-rep stage, all spokes are optimized separately, so there are 3 parameters (length, orientation[2]) per spoke
	
	del_Spoke = ( log(del_r), del_theta, del_phi )  // check function Vector3D::rotateBy( delTheta, delPhi )

	The vector to be applied is,
	x = { del_Position, del_Spoke_0, del_Spoke_1, del_Spoke_End } = 3 + 3 + 3 + 3(end) = 9 or 12 entries

	del_Position = ( del_P_X, del_P_Y, del_P_Z )
	
	*/

	Vector3D trans ;

	double scale0, scale1, scaleEnd ;
	double delTheta0, delTheta1, delThetaEnd ;
	double delPhi0, delPhi1, delPhiEnd ;


	trans =  ( x(0) * Vector3D(1.0, 0.0, 0.0) ) ; 
	trans += ( x(1) * Vector3D(0.0, 1.0, 0.0) ) ;
	trans += ( x(2) * Vector3D(0.0, 0.0, 1.0) ) ;

	trans *= tuningWt(AtomTranslationFactor);

	// If atom movement is restricted, then the atom hub does not move
	if( int(tuningWt(SRepAtomRestrictMovement)) ) {

		static bool reported	= false;

		if( !reported ) {
			cout << "Disabling atom hub movement in s-rep stage" << endl ;
			reported = 1 ;
		}

		trans.setX( 0.0 ) ; trans.setY( 0.0 ) ; trans.setZ( 0.0 ) ;
	}

	prim.translateBy(trans);	
	
	scale0		= exp(tuningWt(AtomScaleFactor) * x(3));
	delTheta0	= tuningWt(AtomThetaFactor) * x(4) ;
	delPhi0		= tuningWt(AtomThetaFactor) * x(5) ;

	scale1		= exp(tuningWt(AtomScaleFactor) * x(6));
	delTheta1	= tuningWt(AtomThetaFactor) * x(7) ;
	delPhi1		= tuningWt(AtomThetaFactor) * x(8) ;

	// scale the radii

	cout << "Scaling spoke 0 by " << scale0 << endl ;
	prim.scaleSpokeBy( 0, scale0 ) ;

	cout << "Rotating spoke 0 by ( " << delTheta0 * 180 / 3.1416 << ", " << delPhi0  * 180 / 3.1416 << " ) degrees" << endl ;
	prim.rotateSpokeBy( 0, delTheta0, delPhi0 ) ;

	cout << "Scaling spoke 1 by " << scale1 << endl ;
	prim.scaleSpokeBy( 1, scale1 ) ;

	cout << "Rotating spoke 1 by ( " << delTheta1 * 180 / 3.1416 << ", " << delPhi1  * 180 / 3.1416 << " ) degrees" << endl ;
	prim.rotateSpokeBy( 1, delTheta1, delPhi1 ) ;

	// deal with the end spoke

	if( prim.type() == M3D_END_PRIMITIVE ) {		

		scaleEnd	= exp(tuningWt(AtomScaleFactor) * x(9));
		delThetaEnd = tuningWt(AtomThetaFactor) * x(10) ;
		delPhiEnd	= tuningWt(AtomThetaFactor) * x(11) ;		

		cout << "Scaling spoke End by " << scaleEnd << endl ;
		prim.scaleSpokeBy( 2, scaleEnd ) ;

		cout << "Rotating spoke End by ( " << delThetaEnd * 180 / 3.1416 << ", " << delPhiEnd  * 180 / 3.1416 << " ) degrees" << endl ;
		prim.rotateSpokeBy( 2, delThetaEnd, delPhiEnd ) ;
	}

	//
	// Correct any consistency issues.
	//
	//M3DFigure* currFigPtr = targetObject->getFigurePtr(figureId);
	//currFigPtr->fixGlobalConsistency();
}

#endif 

double M3DSRepProblemBoundsFunction::bound( int i )
{
	assert( i >= 0 );
	
	Vector bounds(3,	log(tuningWt(SRepScaleBound)), 
						tuningWt(SRepThetaBound) * R_PI/180.0, 
						tuningWt(SRepThetaBound) * R_PI/180.0) ;

	return bounds( (i < bounds.size()) ? i : (bounds.size() - 1) );
	
}

double M3DSRepProblemBoundsFunction::factor( int i )
{
	assert( i >= 0 );

	Vector factors(3,	tuningWt(SRepScaleFactor), 
						tuningWt(SRepThetaFactor),
						tuningWt(SRepThetaFactor) );

	return factors( (i < factors.size()) ? i : (factors.size() - 1) );

}

