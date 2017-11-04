#include "M3DPNSTransform.h"

#include "math.h"
#include <iomanip>
#include <vector>
#include <fstream>

/*
	Author: Dibyendusekhar Goswami, July 2011

	Description: A class to hold a PNS transformation and transform from spherical space to Euclidean space (PNS residuals)

	Theory: A PNSTransformation takes data from an n-dimensional spherical space to an n-dimensional Euclidean space.
			A description of this can be found in CPNS_sreps.PDF or Notes_on_CPNS.pdf
			This is known as h (hS, or hSpoke) in these documents

*/

M3DPNSTransform::M3DPNSTransform() {

	nDims = 0 ;
	sphereAxis.clear() ;
	sphereDist = NULL ;
}

M3DPNSTransform::M3DPNSTransform(int _nDims) {

	if( _nDims < 0 ) {
		cout << "Cannot initialize negative dimensional PNS transformation" << endl ;
		return ;
	}

	if( _nDims == 0 )
		M3DPNSTransform() ;

	// Initialize the transformation

	nDims = _nDims ;

	sphereAxis.reserve( nDims ) ;
	for( int i = 0 ; i < nDims ; i ++ ) 
		sphereAxis.push_back( NULL ) ;

	sphereDist = new double [nDims - 1] ;

}

//M3DPNSTransform::M3DPNSTransform(int _nDims, std::vector <double *> &_sphereAxis, double *_sphereDist) {
//
//	if( _nDims < 0 ) {
//		cout << "Cannot initialize negative dimensional PNS transformation" << endl ;
//		return ;
//	}
//
//	if( _nDims == 0 )
//		M3DPNSTransform() ;
//
//	// Initialize the transformation
//
//	nDims = _nDims ;
//
//	if( _sphereAxis.size() != nDims ) {
//		cout << "Cannot initialize PNS transformation with incorrect sphereAxis data" << endl ;
//		return ;
//	}
//
//	sphereAxis.reserve( nDims ) ;
//	sphereDist = new double [nDims - 1] ;
//
//	for( int nd = 0 ; nd < nDims ; nd ++ ) {
//
//		sphereAxis[nd] = _sphereAxis[nd] ;
//
//		if( nd < nDims - 1 )
//			sphereDist[nd] = _sphereDist[nd] ;
//
//	}
//
//}

void M3DPNSTransform::copy(M3DPNSTransform *thatTransform) {

	// cout << "M3DPNSTransform::copy started" << endl ;

	nDims = thatTransform->getNDimensions() ;

	// ----  First clear out all existing variables  ------

	if( ! sphereAxis.empty() ) {
		for( int i = 0 ; i < sphereAxis.size() ; i++ ) {
			if( sphereAxis[i] != NULL ) {
				delete [] sphereAxis[i] ;	
				sphereAxis[i] = NULL ;
			}
		}
		sphereAxis.clear() ;
	}	

	// cout << "SphereAxes cleared" << endl ;

	if( sphereDist != NULL ) {
		delete [] sphereDist ;
		sphereDist = NULL ;
	}

	// cout << "SphereDist cleared" << endl ;

	// ----------------------------------------------------

	sphereAxis.resize( nDims, NULL ) ;

	double * sphereAxisThis = NULL ;	// donot delete this pointer - its passed to this M3DPNSTransform object

	for( int i = 0 ; i < nDims ; i ++ ) {

		// cout << "sphereAxis[" << i << "] copying has started" << endl ;

		int lenAxis = nDims + 1 - i ;
		
		if( i == nDims-1 )
			lenAxis = 1 ;

		// cout << "length of axis = " << lenAxis << endl ;

		sphereAxisThis = NULL ;

		sphereAxisThis = new double [ lenAxis ] ;

		// cout << "space successfully alloted for sphereAxisThis" << endl ;

		for( int j = 0 ; j < lenAxis ; j++ ) 
			sphereAxisThis[j] = thatTransform->getSphereAxisVal( i, j ) ;

		sphereAxis[i] = sphereAxisThis ;		

		// cout << "sphereAxis[" << i << "] has been set" << endl ;
	}

	// cout << "SphereAxes set" << endl ;
	
	sphereDist = new double [nDims-1] ;

	for( int i = 0 ; i < nDims-1 ; i++ )
		sphereDist[i] = thatTransform->getSphereDist(i) ;

	// cout << "SphereDist set" << endl ;

	// cout << "M3DPNSTransform::copy completed successfully" << endl ;
}

//M3DPNSTransform * M3DPNSTransform::assign() {
//
//	// cout << "M3DPNSTransform::assign started" << endl ;
//
//	nDims = thatTransform->getNDimensions() ;
//
//	// ----  First clear out all existing variables  ------
//
//	if( ! sphereAxis.empty() ) {
//		for( int i = 0 ; i < sphereAxis.size() ; i++ ) {
//			if( sphereAxis[i] != NULL ) {
//				delete [] sphereAxis[i] ;	
//				sphereAxis[i] = NULL ;
//			}
//		}
//		sphereAxis.clear() ;
//	}	
//
//	// cout << "SphereAxes cleared" << endl ;
//
//	if( sphereDist != NULL ) {
//		delete [] sphereDist ;
//		sphereDist = NULL ;
//	}
//
//	// cout << "SphereDist cleared" << endl ;
//
//	// ----------------------------------------------------
//
//	sphereAxis.resize( nDims, NULL ) ;
//
//	double * sphereAxisThis = NULL ;	// donot delete this pointer - its passed to this M3DPNSTransform object
//
//	for( int i = 0 ; i < nDims ; i ++ ) {
//
//		cout << "sphereAxis[" << i << "] copying has started" << endl ;
//
//		int lenAxis = nDims + 1 - i ;
//		
//		if( i == nDims-1 )
//			lenAxis = 1 ;
//
//		// cout << "length of axis = " << lenAxis << endl ;
//
//		sphereAxisThis = NULL ;
//
//		sphereAxisThis = new double [ lenAxis ] ;
//
//		// cout << "space successfully alloted for sphereAxisThis" << endl ;
//
//		for( int j = 0 ; j < lenAxis ; j++ ) 
//			sphereAxisThis[j] = thatTransform->getSphereAxisVal( i, j ) ;
//
//		sphereAxis[i] = sphereAxisThis ;		
//
//		// cout << "sphereAxis[" << i << "] has been set" << endl ;
//	}
//
//	// cout << "SphereAxes set" << endl ;
//	
//	sphereDist = new double [nDims-1] ;
//
//	for( int i = 0 ; i < nDims-1 ; i++ )
//		sphereDist[i] = thatTransform->getSphereDist(i) ;
//
//	// cout << "SphereDist set" << endl ;
//
//	// cout << "M3DPNSTransform::copy completed successfully" << endl ;
//}


M3DPNSTransform::~M3DPNSTransform() {

	// delete each vector inside sphereAxis

	for( int i = 0 ; i < sphereAxis.size() ; i ++ ) {

		if( sphereAxis[i] != NULL ) {			
			delete [] sphereAxis[i] ;
			sphereAxis[i] = NULL ;
		}
	}

	sphereAxis.clear() ;

	if( sphereDist != NULL ) {
		delete [] sphereDist ;
		sphereDist = NULL ;
	}

}



void M3DPNSTransform::setNDimensions(int _nDims) {

	// this function necessarily performs the same function as the constructor(int _nDims)

	if( _nDims < 0 ) {
		cout << "Cannot initialize negative dimensional PNS transformation" << endl ;
		return ;
	}

	if( _nDims == 0 )
		M3DPNSTransform() ;

	// Initialize the transformation

	nDims = _nDims ;

	sphereAxis.clear() ;

	sphereAxis.reserve( nDims ) ;
	for( int i = 0 ; i < nDims ; i ++ ) 
		sphereAxis.push_back( NULL ) ;

	if( sphereDist != NULL )
		delete [] sphereDist ;

	sphereDist = new double [nDims - 1] ;

}

bool M3DPNSTransform::setSphereDist(double *_sphereDist, int lenSphereDist) {

	if( lenSphereDist != ( this->nDims-1 ) ) {	// wrong dimensions
		cout << "Wrong dimensional input to M3DPNSTransform::setSphereDist()" << endl ;
		return (0) ;				
	}

	if( this->sphereDist == NULL ) {
		cout << "Error in M3DPNSTransform::setSphereDist() !. sphereDist has to be set in a constructor or copy() function of PNSTransform" << endl ;
		return(0) ;
	}

	for( int i = 0 ; i < this->nDims - 1 ; i++ ) 
		sphereDist[i] = _sphereDist[i] ;

	return(1) ;

}
bool M3DPNSTransform::setSphereAxis(int nSphere, double *sphereAxisN, int lenSphereAxisN){

	// ---------------   input verification  ------------------------

	// verify sphere number
	if( nSphere > this->nDims ) {
		cout << "Error in M3DPNSTransform::setSphereAxis() ! " << endl ;
		cout << "Sphere number > nDims" << endl ;
		return(0) ;
	}

	// verify input length
	if( nSphere < nDims-1 ) {
		if( lenSphereAxisN != (this->nDims+1-nSphere) ) {
			cout << "Error in M3DPNSTransform::setSphereAxis() ! Wrong input length " << endl ;			
			return(0) ;
		}
	}
	else {	// the last sphere is of length-1
		if( lenSphereAxisN != 1 )  {
			cout << "Error in M3DPNSTransform::setSphereAxis() ! Wrong input length " << endl ;			
			return(0) ;
		}
	}

	// ---------------   set values			------------------------

	double * sphereAxisThis = new double[ lenSphereAxisN ] ;

	for( int i = 0 ; i < lenSphereAxisN ; i ++ )
		sphereAxisThis[i] = sphereAxisN[i] ;

	this->sphereAxis[nSphere] = sphereAxisThis ;

	return(1) ;

}

int M3DPNSTransform::getNDimensions() {
	return( nDims ) ;
}

double M3DPNSTransform::getSphereDist(int n) {

	if( sphereDist == NULL || n > nDims-1 || n < 0 ) {
		cout << "Error in M3DPNSTransform::getSphereDist(int n) !" << endl ;
		return(0.0) ;
	}

	return( this->sphereDist[n] ) ;
}

double M3DPNSTransform::getSphereAxisVal(int nSphere, int index) {

	if( sphereAxis.empty() || sphereAxis[nSphere] == NULL ) {
		cout << "sphereAxis not initialized yet! " << endl ;
		return( 0.0 ) ;
	}

	if( nSphere == nDims && index > 0 ) {
		cout << "last sphereAxis is of length 1, not 2" << endl ;
		return(0.0) ;
	}

	return( *(sphereAxis[nSphere] + index) ) ;

}

double * M3DPNSTransform::convertS2E( double * S ) {

	if( S == NULL ) {
		cout << "Error ! NULL input to M3DPNSTransform::convertS2E()" << endl ;
		return( NULL ) ;
	}

	if( this->nDims <= 0 ) {
		cout << "Error ! M3DPNSTransform::convertS2E() called before initializing PNSTransfom" << endl ;
		return( NULL ) ;
	}

	// If sphere is k-dimensional, then length of input vector S = k+1
	//									length of output vector E = k

	double * E = new double[ nDims ] ;

	for( int i = 0 ; i < nDims ; i ++ )
		E[i] = 0.0 ;

	// variables used during the conversion process

	vector <double> currentSphere( nDims+1, 0.0 ) ;

	vector <double> nestedSphere( nDims+2, 0.0 ) ;

	vector <double> currentAxis ;

	vector <double> currentRes( nDims, 0.0 ) ;

	double currentDist = 0.0 ;	

	double innerProd = 0.0 ;

	double scaleFactor = 1.0 ;	// the scaling factor to make residuals commensurate

	double ** rotMatrix ;		// dont forget to delete this at every iteration

	// Initialize currentSphere, nestedSphere, and currentRes


	for( int i = 0 ; i < nDims+1 ; i ++ )
		currentSphere[i] = S[i] ;

	// ------------------ COMPUTE FOR 1 : N-1 SPHERES  -------------------------------

	for( int nd = 0 ; nd < nDims-1 ; nd ++ ) {

		// ---------------------  load currentAxis  ----------------------------
	
		// getSphereAxis(): clear the contents of currentAxis and load into it the current Sphere

		if( ! getSphereAxis( currentAxis, nd ) ) {
			cout << "Error in M3DPNSTransform::convertS2E()! Cannot load sphereAxis for sphere no." << nd << endl ;
			currentSphere.clear() ;
			nestedSphere.clear() ;
			currentAxis.clear() ;
			currentRes.clear() ;
			delete [] E ;
			E = NULL ;
			return( NULL ) ;
		}

		// ---------------------  load currentDist  ----------------------------

		currentDist = getSphereDist( nd ) ;
		
		// --------------- calculate inner product -------------------------
		innerProd = 0.0 ;

		for( int i = 0 ; i < nDims+1 - nd ; i++ )
			innerProd += currentSphere[i] * currentAxis[i] ;

		// ----------  set the residue value (scaled) ----------------------

		currentRes[nd] =  acos( innerProd ) - currentDist  ;

		double resThis = currentRes[nd] ;

		currentRes[nd] *= scaleFactor ; 

		// ----------  update scale factor (prod(sin(ri))  -----------------

		scaleFactor *= sin(getSphereDist(nd)) ;

		// ----------  update currentSphere  -------------------------------

		// Get the rotation matrix to the north pole

		rotMatrix = getRotationMatrixToNorthPole( currentAxis ) ;


		if( rotMatrix == NULL ) {
			cout << "Error in M3DPNSTransform::convertS2E()! Rotation matrix is NULL for sphere no." << nd << endl ;
			currentSphere.clear() ;
			nestedSphere.clear() ;
			currentAxis.clear() ;
			currentRes.clear() ;
			delete [] E ;
			E = NULL ;
			return(NULL) ;
		}

		// Erase the last element of nestedSphere, thus reducing its dimension by one 
		// because this nestedSphere will be one-dimension smaller than the last one

		nestedSphere.erase( nestedSphere.end()-1 ) ;

		for( int i = 0 ; i < nDims+1-nd ; i++ ) {
			nestedSphere[i] = 0.0 ;
			for( int j = 0 ; j < nDims+1-nd ; j++ )
				nestedSphere[i] += rotMatrix[i][j] * currentSphere[j] ;
		}

		// Erase the last element of currentSphere, thus reducing its dimension by one 
		// because next currentSphere will be one-dimension smaller than this one

		currentSphere.erase( currentSphere.end()-1 ) ;

		double denominator = sqrt( 1 - nestedSphere[ nDims-nd ] * nestedSphere[ nDims-nd ] ) ;

		if( denominator == 0 ) {
			cout << "Error in M3DPNSTransform::convertS2E()! Denominator is zero for sphere no." << nd << endl ;
			currentSphere.clear() ;
			nestedSphere.clear() ;
			currentAxis.clear() ;
			currentRes.clear() ;
			delete [] E ;
			E = NULL ;
			for( int i = 0 ; i < nDims+1 - nd ; i++ )
				delete [] rotMatrix[i] ;

			delete [] rotMatrix ;
			return(NULL) ;
		}

		for( int i = 0 ; i < nDims-nd ; i++ )
			currentSphere[i] = nestedSphere[i] / denominator ; 

		// delete the rotationMatrix 
		for( int i = 0 ; i < nDims+1 - nd ; i++ )
			delete [] rotMatrix[i] ;
		delete [] rotMatrix ;

		rotMatrix = NULL ;

	}

	// ------------------ COMPUTE FOR THE LAST SPHERE  -------------------------------

	// Define variables

	double pi = 3.141592653589793 ;

	double s1ToRadian = atan2( currentSphere[1], currentSphere[0] ) ;

	double dividend = s1ToRadian - getSphereAxisVal(nDims-1, 0) + pi ;

	// find out the modulus for non-integer division	
	// currentRes[ last ] = mod( numerator, 2 * pi )

	while( dividend < 0 )
		dividend += 2*pi ;

	while( dividend >= 2*pi ) 
		dividend -= 2*pi ;

	dividend -= pi ;

	currentRes[ nDims-1 ] = scaleFactor * dividend ;

	// ------------------  Reverse the residuals and return E ---------------------------

	for( int nd = 0 ; nd < nDims ; nd ++ ) 
		E[nd] = currentRes[nDims-1-nd] ;
	

	// ---------------  Deallocate space used --------------------------------

	currentSphere.clear() ;
	nestedSphere.clear() ;
	currentAxis.clear() ;
	currentRes.clear() ;	

	return( E ) ;

}

double * M3DPNSTransform::convertE2S( double * E ) {

	if( E == NULL ) {
		cout << "Error ! NULL input to M3DPNSTransform::convertE2S()" << endl ;
		return(NULL) ;
	}

	if( this->nDims <= 0 ) {
		cout << "Error ! M3DPNSTransform::convertE2S() called before initializing PNSTransfom" << endl ;
		return(NULL) ;
	}

	// If sphere is k-dimensional, then length of input vector E = k
	//									length of output vector S = k+1

	double * S = new double[ nDims+1 ] ;

	for( int i = 0 ; i < nDims+1 ; i ++ )
		S[i] = 0.0 ;

	// variables used during the conversion process

	vector <double> currentAxis ;

	vector <double> prevSphere ;
	vector <double> currentSphere ;

	vector <double> residuals( nDims, 0.0 ) ;

	double ** rotMatrix ;		// dont forget to delete this at every iteration

	// -----------------  Unscale to recover residuals  -------------------------

	double scaleFactor = 1.0 ;	// the scaling factor to make residuals commensurate

	// residual[0] => for the highest-dimensional sphere, sphere[0]

	for( int nd = 0 ; nd < nDims ; nd ++ ) {

		residuals[nd] = E[ nDims-nd-1 ] / scaleFactor ;

		if( nd < nDims-1 )
			scaleFactor *= sin(getSphereDist(nd)) ; 		
	}

	// ---------------  Map from S^1 (sphere[nDims-1]) to S^2 (sphere[nDims-2]) --------------------------------------

	// get value of last sphereAxis (a single number, not a 2-vector)

	double lastSphereAxis = getSphereAxisVal( nDims-1, 0 ) ;

	// get value of second last sphereAxis
	
	if( ! getSphereAxis( currentAxis, nDims-2 ) ) {
		cout << "Error in M3DPNSTransform::convertE2S()! " ;
		cout << "Cannot load sphereAxis for second lowest-dimensional sphere " << endl ;	
		currentAxis.clear() ;
		residuals.clear() ;
		delete [] S ;
		S = NULL ;
		return( NULL ) ;
	}

	// get the rotation matrix for sphere no.[ nDims-2 ]

	rotMatrix = getRotationMatrixToNorthPole( currentAxis ) ;

	if( rotMatrix == NULL ) {
		cout << "Error in M3DPNSTransform::convertE2S()! Rotation matrix is NULL for sphere no." << nDims-2 << endl ;
		currentAxis.clear() ;
		residuals.clear() ;
		delete [] S ;
		S = NULL ;
		return( NULL ) ;
	}

	vector < vector <double> > rotMatTranspose (currentAxis.size(), vector <double> (currentAxis.size(), 0.0)) ;

	for( int i = 0 ; i < currentAxis.size() ; i ++) {
		for( int j = 0 ; j < currentAxis.size() ; j++ ) 
			rotMatTranspose[i][j] = rotMatrix[j][i] ;
	}

	// deallocate rotation matrix

	for( int i = 0 ; i < currentAxis.size() ; i++ )
		delete [] rotMatrix[i] ;
	delete [] rotMatrix ;
	rotMatrix = NULL ;

	/*
		sphere[k] = R'(v[k]) * [	sin(r[k] + res[k]) * sphere[k-1] ;
									cos(r[k] + res[k])					]
	*/

	// ------------ Construct S^2 sphere[nDims-2] ----------------------------------------------------------------------

	// Initialize prevSphere to S^2

	prevSphere.resize( currentAxis.size(), 0.0 ) ;	

	prevSphere[0] = sin( getSphereDist(nDims-2) + residuals[nDims-2] ) * cos( lastSphereAxis + residuals[nDims-1] ) ;
	prevSphere[1] = sin( getSphereDist(nDims-2) + residuals[nDims-2] ) * sin( lastSphereAxis + residuals[nDims-1] ) ;
	prevSphere[2] = cos( getSphereDist(nDims-2) + residuals[nDims-2] ) ;

	currentSphere.resize( currentAxis.size(), 0.0 ) ;

	for( int i = 0 ; i < currentAxis.size() ; i++ ) {
		for( int j = 0 ; j < prevSphere.size() ; j++ )
			currentSphere[i] += rotMatTranspose[i][j] * prevSphere[j] ;
	}

	// ------------ Construct S^k+1, sphere[nDims-(k+1)] from S^k, sphere[nDims-k]  ------------------------------------

	for( int nd = 3 ; nd <= nDims ; nd ++ ) {

		prevSphere.clear() ;

		prevSphere = currentSphere ;		

		// ---------  get sphereAxis[ nDims - nd ]  ---------------------

		// getSphereAxis(): clear the contents of currentAxis and load into it the current Sphere
		if( ! getSphereAxis( currentAxis, nDims-nd ) ) {
			cout << "Error in M3DPNSTransform::convertE2S()! " ;
			cout << "Cannot load sphereAxis for sphere no." << nDims-nd << endl ;	
			currentAxis.clear() ;
			residuals.clear() ;
			currentSphere.clear() ;
			prevSphere.clear() ;
			rotMatTranspose.clear() ;
			delete [] S ;
			S = NULL ;
			return( NULL ) ;
		}

		// ------------  get rotation matrix transpose  ----------------------

		rotMatrix = getRotationMatrixToNorthPole( currentAxis ) ;

		if( rotMatrix == NULL ) {
			cout << "Error in M3DPNSTransform::convertE2S()! Rotation matrix is NULL for sphere no." << nDims-nd << endl ;
			currentAxis.clear() ;
			residuals.clear() ;
			currentSphere.clear() ;
			prevSphere.clear() ;
			rotMatTranspose.clear() ;
			delete [] S ;
			S = NULL ;
			return( NULL ) ;
		}
		
		rotMatTranspose.clear() ;

		rotMatTranspose.resize( nd+1, vector <double> ( nd+1, 0.0 )) ;

		for( int i = 0 ; i < nd+1 ; i ++) {
			for( int j = 0 ; j < nd+1 ; j++ ) 
				rotMatTranspose[i][j] = rotMatrix[j][i] ;
		}

		// ----------------  get current sphere from previous sphere  ----------------------------------

		currentSphere.clear() ;
		currentSphere.resize( nd+1, 0.0 ) ;

		vector <double> tempSphere( nd+1, 0.0 ) ;
		for( int i = 0 ; i < nd ; i ++ )
			tempSphere[i] = sin( getSphereDist(nDims-nd) + residuals[nDims-nd] ) * prevSphere[i] ;
		tempSphere[nd] = cos( getSphereDist(nDims-nd) + residuals[nDims-nd] ) ;

		for( int i = 0  ; i < nd+1 ; i ++ ) {
			for( int j = 0 ; j < nd+1 ; j++ )
				currentSphere[i] += rotMatTranspose[i][j] * tempSphere[j] ;
		}

		tempSphere.clear() ;

		// ---------------  deallocate rotation matrix  -----------------------------

		for( int i = 0 ; i < nd+1 ; i++ )
			delete [] rotMatrix[i] ;
		delete [] rotMatrix ;
		rotMatrix = NULL ;

	}

	for( int i = 0 ; i < nDims+1 ; i ++ ) {
		S[i] = currentSphere[i] ;
	}

	// ------------------------------  debugging   -------------------------------------

#undef DEBUG_DIBYENDU_S2E

#ifdef DEBUG_DIBYENDU_S2E

	char debugFilename[100] ;
	sprintf( debugFilename, "../../bin/debug_E2S.txt" ) ;
	std::ofstream f ;
	f.open( (const char*) debugFilename,  ios::trunc);		

	if( f ) {
		f << setiosflags(ios::fixed) ;
		cout << debugFilename << " has been opened for writing debug info for convertE2S" << endl ;
		// writing residuals
		//f << "residuals = " ;
		for( int i = 0 ; i < residuals.size() ; i++)
			f << setw(12) << setprecision(8) << residuals[i] << endl ;		
	}
	f.close() ;

#endif 

	// ---------------  Deallocate space used --------------------------------

	currentAxis.clear() ;
	residuals.clear() ;
	currentSphere.clear() ;
	prevSphere.clear() ;
	rotMatTranspose.clear() ;

	return( S ) ;

}
bool M3DPNSTransform::getSphereAxis(std::vector<double> &sphereAxisN, int N) {

	if( N > nDims-1 ) {
		cout << "Error in M3DPNSTransform::getSphereAxis()! N exceeds dimensions" << endl ;
		return(0) ;
	}

	if( sphereAxis[N] == NULL ) {
		cout << "Error in M3DPNSTransform::getSphereAxis()! Queried sphereAxis is NULL" << endl ;
		return(0) ;
	}

	if( ! sphereAxisN.empty() )
		sphereAxisN.clear() ;

	sphereAxisN.reserve( nDims+1 - N ) ;

	if( N < nDims-1 ) {
		for( int i = 0 ; i < nDims+1-N ; i++ )
			sphereAxisN.push_back( *(sphereAxis[N] + i) ) ;
	}
	else
		sphereAxisN.push_back( *(sphereAxis[N]) ) ;

	return(1) ;

}

double * * M3DPNSTransform::getRotationMatrixToNorthPole( std::vector<double> &v ) {

	int lenVector = v.size() ;

	if( lenVector < 2 ) {
		cout << "Error! Less than 2 - dimensional input to M3DPNSTransform::getRotationMatrixToNorthPole()" << endl ;
		return( NULL ) ;
	}

	// -----------------  Initialize output rotation matrix  ----------------------------

	// rotMatrix[ lenVector ][ lenVector ]

	double* * rotMatrix = new double* [ lenVector ] ;

	for( int i = 0 ; i < lenVector ; i++ )
		rotMatrix[i] = new double [ lenVector ] ;

	for( int i = 0 ; i < lenVector ; i++ ) {
		for( int j = 0 ; j < lenVector ; j++ )
			rotMatrix[i][j] = 0.0 ;
	}

	// ------------------  Form intermediate variables  ---------------------------------

	// Normalize input vector 

	double normV = 0.0 ;

	for( int i = 0 ; i < lenVector ; i++ ) 
		normV += v[i] * v[i] ;

	normV = sqrt( normV ) ;

	for( int i = 0 ; i < lenVector ; i++ )
		v[i] /= normV ;

	// Take inner product with north pole (also called vector a in the following ) ev = (0, 0, ... 1), dot product = v[last]
	
	double innerProd = v[lenVector-1] ;

	// if innerProd is close to 1, return identity matrix

	if( fabs(innerProd - 1) < 1e-10 ) {
		for( int i = 0 ; i < lenVector ; i++ )
			rotMatrix[i][i] = 1.0 ;
		return( rotMatrix ) ;
	}
	else if( fabs(innerProd + 1) < 1e-10 ) {
		for( int i = 0 ; i < lenVector ; i++ )
			rotMatrix[i][i] = -1.0 ;
		return( rotMatrix ) ;
	}

	// angle betn the vectors

	double alpha = acos( innerProd ) ;

	double sinAlpha = sin( alpha ) ;
	double cosAlpha = innerProd ;

	// define: vector c = v - a * (innerProd) = [ v[1:n-1] ; 0.0 ]

	vector <double> c ; 

	c.assign( v.begin(), v.end() ) ;

	c[lenVector-1] = 0.0 ;

	// normalize the vector c

	double normC = sqrt( 1 - v[lenVector-1] * v[lenVector-1] ) ;

	for( int i = 0 ; i < lenVector ; i ++ )
		c[i] /= normC ;	

	// define: matrix A = a * c' - c * a'

	vector < vector<double> > A( lenVector, vector <double> ( lenVector, 0.0 ) ) ;

	/*	
		A = a * c' - c * a' 
			a * c' = [	zeros( n-1, n ) ; 
						c'					]
			c * a' = [ zeros( n, n-1 )   c  ]
		A = [	zeros( n-1, n-1 ) -c( 1 : n-1 ) ; 
				c'( 1 : n-1 )					]

		We later need sinAlpha * A, so premultiply by sinAlpha
		
	*/ 

	for( int i = 0 ; i < lenVector-1 ; i ++ ) {

		// fill up the last row
		A[ lenVector-1 ][i] += sinAlpha * c[i] ;

		// fill up the last column
		A[i][ lenVector-1 ] += - sinAlpha * c[i] ;
	}

	// fill up the diagonal elements (adding the identity matrix)		

	for( int i = 0 ; i < lenVector ; i++ )
		A[i][i] += 1.0 ;

	/* 
		rotMatrix	= Identity(lenVector) + sinAlpha * A + (cosAlpha - 1)*( v * v' + c * c' )
					= A + mat, where mat = (cosAlpha - 1)(a * a' + c * c')

		rotMatrix[i][j] = A[i][j] + (cosAlpha-1) * c[i]*c[j] ;
		rotMatrix[last][last] += cosAlpha - 1 ; // for a * a'
	*/ 

	for( int i = 0 ; i < lenVector ; i ++ ) {
		for( int j = 0 ; j < lenVector ; j++ ) 
			rotMatrix[i][j] = A[i][j] + (cosAlpha-1.0) * c[i]*c[j] ;
	}

	rotMatrix[ lenVector-1 ][ lenVector-1 ] += (cosAlpha - 1.0) ;	// for adding (cosAlpha-1) (a * a')

	// deallocate used memory

	if( ! c.empty() )
		c.clear() ;

	if( ! A.empty() )
		A.clear() ;

	return( rotMatrix ) ;
}
