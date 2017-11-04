#include "M3DCPNSStats.h"

#include <fstream>
#include <iomanip>

#include <iostream>

/*
	Author: Dibyendusekhar Goswami, July 2011

	Description: A class to hold CPNS statistics for a CPNS mean model

	Theory: The CPNS mean model contains eigenvectors and eigenvalues, and also transformations, known as PNSTransforms 
			that help to perform CPNS on an s-rep.
			A PNSTransformation takes data from an n-dimensional spherical space to an n-dimensional Euclidean space.
			A detailed description of this can be found in CPNS_sreps.PDF or Notes_on_CPNS.pdf
			
			CPNS can be used to observe the deformation of the mean model along principal eigenmodes
			and also to calculate Mahalanobis distance from the mean

*/

M3DCPNSStats::M3DCPNSStats() {

	nSpokes = -1 ;

	nEigenmodes = -1 ;

	eigenVectorLength = -1 ;

	nAtomRows = -1 ;
	nAtomCols = -1 ;

	PNSShape = NULL ;

	PNSSpoke.clear() ;

	scaleSpoke = NULL ;

	eigenVector.clear() ;

	eigenValue = NULL ;
}

M3DCPNSStats::M3DCPNSStats( int _nSpokes, int _nEigenmodes, int _eigenVectorLength, int _nAtomRows, int _nAtomCols ) {

	nSpokes = _nSpokes ;

	nEigenmodes = _nEigenmodes ;

	eigenVectorLength = _eigenVectorLength ;

	nAtomRows = _nAtomRows ;
	nAtomCols = _nAtomCols ;


	PNSShape = NULL ;	

	meanShape = Vector3D( 0.0, 0.0, 0.0 ) ;

	scaleShape = 0.0 ;	

	PNSSpoke.reserve( nSpokes ) ;
	for( int i = 0 ; i < nSpokes ; i ++ )
		PNSSpoke.push_back( NULL ) ;

	scaleSpoke = new double[ nSpokes ] ;
	for( int i = 0 ; i < nSpokes ; i++ )
		scaleSpoke[i] = 0.0 ;

	eigenValue = new double[ nEigenmodes ] ;
	for( int i = 0 ; i < nEigenmodes ; i++ )
		eigenValue[i] = 0.0 ;	

	eigenVector.reserve( nEigenmodes ) ;
	double* thisEigenVector = NULL ;
	for( int i = 0 ; i < nEigenmodes ; i++ ) {
		thisEigenVector = new double[ eigenVectorLength ] ;
		eigenVector.push_back( thisEigenVector ) ;
		thisEigenVector = NULL ;
	}
}

M3DCPNSStats::~M3DCPNSStats() {

	if( PNSShape != NULL ) {
		delete PNSShape ;
		PNSShape = NULL ;
	}

	if( ! PNSSpoke.empty() ) {
		for( int i = 0 ; i < PNSSpoke.size() ; i ++ ) {
			if( PNSSpoke[i] != NULL ) {
				delete PNSSpoke[i] ;
				PNSSpoke[i] = NULL ;
			}
		}
		PNSSpoke.clear() ;
	}

	if( scaleSpoke != NULL ) {
		delete [] scaleSpoke ;
		scaleSpoke = NULL ;
	}

	if( ! eigenVector.empty() ) {
		for( int i = 0 ; i < eigenVector.size() ; i ++ ) {
			if( eigenVector[i] != NULL ) {
				delete [] eigenVector[i] ;
				eigenVector[i] = NULL ;
			}
		}
		eigenVector.clear() ;
	}

	if( eigenValue != NULL ) {
		delete [] eigenValue ;
		eigenValue = NULL ;
	}

}

bool M3DCPNSStats::setScaleShape( double _scaleShape ) {

	if( _scaleShape <= 0 ) {
		cout << "Error in M3DCPNSStats::setScaleShape() ! Cannot set scaleShape <= 0." << endl ;
		return(0) ;
	}
	else
		this->scaleShape = _scaleShape ;

	return(1) ;

}
bool M3DCPNSStats::setMeanShape( double* _meanShape ) {

	if( _meanShape == NULL ) {
		cout << "Error in M3DCPNSStats::setMeanShape() ! Null double pointer passed" << endl ;
		return(0) ;
	}
	else
		this->meanShape = Vector3D( _meanShape[0], _meanShape[1], _meanShape[2] ) ;

	return(1) ;

}

bool M3DCPNSStats::setScaleSpoke( double* _scaleSpoke ) {

	if( _scaleSpoke == NULL ) {
		cout << "Error in M3DCPNSStats::setScaleSpoke() ! Null double pointer passed" << endl ;
		return(0) ;
	}
	
	if( scaleSpoke == NULL )
		scaleSpoke = new double[nSpokes] ;

	for( int i = 0 ; i < nSpokes ; i ++ )
		scaleSpoke[i] = _scaleSpoke[i] ;

	return(1) ;
}

// donot delete the pointer that is passed - its owned by cpnsStats

bool M3DCPNSStats::setPNSShape(M3DPNSTransform *_PNSShape) {

	// cout << "M3DCPNSStats::setPNSShape started" << endl ;

	if( this->PNSShape != NULL ) {		
		delete this->PNSShape ;
		this->PNSShape = NULL ;
	}

	this->PNSShape = _PNSShape ;

	// Activate the two following lines of you donot want the input argument to be owned by this object
	// Has caused unknown problems in the past with Visual Studio

	//this->PNSShape = new M3DPNSTransform() ;
	//this->PNSShape->copy( _PNSShape ) ;


	//cout << "M3DCPNSStats::setPNSShape completed successfully" << endl ;

	return(1) ;
}

bool M3DCPNSStats::setPNSSpoke(std::vector<M3DPNSTransform*> & _PNSSpoke) {

	// -------------------------  input verification  ------------------------------------

	if( _PNSSpoke.size() != nSpokes ) {
		cout << "Error ! Wrong size PNSSpoke[] vector passed to M3DCPNSStats::setPNSSpoke() !" << endl ;
		return(0) ;
	}

	for( int i = 0 ; i < nSpokes ; i++ ) {
		if( _PNSSpoke[i] == NULL ) {
			cout << "Error ! Null pointer (PNSSpoke) passed to M3DCPNSStats::setPNSSpoke()!" << endl ;
			return(0) ;
		}
	}

	// -----------------------------------------------------------------------------------

	M3DPNSTransform * pnsSpokeThis = NULL ;

	for( int ns = 0 ; ns < nSpokes ; ns ++ ) {
		 pnsSpokeThis = new M3DPNSTransform(2) ;

		 pnsSpokeThis->copy( _PNSSpoke[ns] ) ;

		 this->PNSSpoke[ns] = pnsSpokeThis ;

		 pnsSpokeThis = NULL ;
	}

	return(1) ;

}

bool M3DCPNSStats::setEigenModes(double *_eigenValue, std::vector<double*> & _eigenVector){

	// ---------------------- Input verification  ----------------------------------

	if( _eigenValue == NULL || _eigenVector.size() != nEigenmodes ) {
		cout << "Error in M3DCPNSStats::setEigenModes()!. Wrong inputs" << endl ;
		cout << "No. of eigenVectors is: " << _eigenVector.size() << endl ;
		return(0) ;
	}
	for( int i = 0 ; i < nEigenmodes ; i++ ) {
		if( _eigenVector[i] == NULL ) {
			cout << "Error in M3DCPNSStats::setEigenModes()!. eigenVector no." << i << " is NULL"  << endl ;
			return(0) ;
		}
	}

	// setting eigenValues

	if( eigenValue == NULL )
		eigenValue = new double[nEigenmodes] ;

	for( int i = 0 ; i < nEigenmodes ; i ++ )
		eigenValue[i] = _eigenValue[i] ;

	// setting eigenVectors

	if( eigenVector.empty() ) {
		eigenVector.reserve( nEigenmodes ) ;
		for( int i = 0 ; i < nEigenmodes ; i++ )
			eigenVector.push_back( NULL ) ;
	}

	double * thisEigenVector = NULL ;

	for( int i = 0 ; i < nEigenmodes ; i ++ ) {

		thisEigenVector = new double[eigenVectorLength] ;

		for( int j = 0 ; j < eigenVectorLength ; j++ )
			thisEigenVector[j] = *( _eigenVector[i] + j ) ;

		eigenVector[i] = thisEigenVector ;

		thisEigenVector = NULL ;

	}

	return(1) ;

}

int M3DCPNSStats::getNSpokes() {
	return( nSpokes ) ;
}

int M3DCPNSStats::getNEigenModes() {
	return( nEigenmodes ) ;
}

int M3DCPNSStats::getEigenVectorLength() {
	return( eigenVectorLength ) ;
}
M3DPNSTransform * M3DCPNSStats::getPNSShape() {
	return( this->PNSShape ) ;
}

Vector3D M3DCPNSStats::getMeanShape() {
	return( Vector3D( meanShape.getX(), meanShape.getY(), meanShape.getZ() ) ) ;
}

M3DPNSTransform * M3DCPNSStats::getPNSSpoke(int N) {

	if( N > nSpokes-1 ) {
		cout << "Error in M3DCPNSStats::getPNSSpoke()! input exceeds number of spokes" << endl ;
		return( NULL ) ;
	}

	return( this->PNSSpoke[N] ) ;
}

double M3DCPNSStats::getScaleSpoke( int N ) {

	if( N > nSpokes-1 ) {
		cout << "Error in M3DCPNSStats::getScaleSpoke()! input exceeds number of spokes" << endl ;
		return( -1 ) ;
	}

	return( scaleSpoke[N] ) ;

}

double M3DCPNSStats::getEigenValue( int N ) {

	if( N > nEigenmodes ) {
		cout << "Error in M3DCPNSStats::getEigenValue()! Input greater than nEigenmodes" << endl ;
		return(-1) ;
	}

	return( eigenValue[N] ) ;
}
bool M3DCPNSStats::writeToFile(const char *fileName) {	

	std::ofstream f;

	f.open(fileName,  ios::trunc);	

	if( ! f ) {
		string s("Unable to open " + string(fileName) + " for writing CPNS Stats!" ) ;
		cout << s << endl ;
		return(0) ;
	}
	else
		cout << fileName << " has been opened for writing CPNS Stats!" << endl ;	

	f << setiosflags( ios::fixed ) ;

	f << "CPNSStats {" << endl ;

		f << "scaleShape = " << setw(12) << setprecision(8) << scaleShape << " ; " << endl ; 

		f << "meanShape = { 3 2 "	<< setw(12) << setprecision(8) << meanShape.getX() << " "
									<< setw(12) << setprecision(8) << meanShape.getY() << " "
									<< setw(12) << setprecision(8) << meanShape.getZ() << " } ;" << endl ;

		// -----------------------  PNSShape ------------------------------------------

		f << "PNSShape {" << endl ;

			int nDimsThis = PNSShape->getNDimensions() ;

			f << "nDims = " << nDimsThis << " ; " << endl ;

			for( int ns = 0 ; ns < nDimsThis ; ns ++ ) {				

				f << "sphereAxis[" << ns << "] = { " ; 
				
				if( ns < nDimsThis - 1 )
					f << nDimsThis + 1 - ns << " 2 " ;
				else
					f << "1 2 " ;

				if( ns < nDimsThis-1 ) {

					for( int j = 0 ; j < nDimsThis+1 - ns ; j++ )
						f << setw(12) << setprecision(8) << PNSShape->getSphereAxisVal( ns, j ) << " " ;
				}
				else
					f << setw(12) << setprecision(8) << PNSShape->getSphereAxisVal( ns, 0 ) ;				

				f << "} ;" << endl ;	// end of sphereAxis[ns]
			}

			f << "sphereDist = { " << nDimsThis-1 << " 2 " ;

			for( int i = 0 ; i < nDimsThis-1 ; i++ )
				f << setw(12) << setprecision(8) << PNSShape->getSphereDist(i) << " " ;

			f << " } ; " << endl ; // end of sphereDist

		f << "}" << endl ;	// end of PNSShape

		// -----------------------------------------------------------------------------

		f << "nSpokes = " << nSpokes << " ; " << endl ;

		f << "scaleSpoke = { " << nSpokes << " 2 " ;

		for( int i = 0 ; i < nSpokes ; i ++ ) 
			f << setw(12) << setprecision(8) << scaleSpoke[i] << " " ;

		f << " } ; " << endl ;		// end of scaleSpoke

		// -----------------------  PNSSpokes ------------------------------------------

		for( int ns = 0 ; ns < nSpokes ; ns ++ ) {

			f << "PNSSpoke[" << ns << "] { " << endl ;

				nDimsThis = PNSSpoke[ns]->getNDimensions() ;

				f << "nDims = " << nDimsThis << " ; " << endl ;

				for( int d = 0 ; d <  nDimsThis ; d ++ ) {

					f << "sphereAxis[" << d << "] = { " ; 
					
					if( d < nDimsThis - 1 )
						f << nDimsThis + 1 - d << " 2 " ;
					else
						f << "1 2 " ;

					if( d < nDimsThis-1 ) {

						for( int j = 0 ; j < nDimsThis+1 - d ; j++ )
							f << setw(12) << setprecision(8) << PNSSpoke[ns]->getSphereAxisVal( d, j ) << " " ;
					}
					else
						f << setw(12) << setprecision(8) << PNSSpoke[ns]->getSphereAxisVal( d, 0 ) ;				

					f << "} ;" << endl ;	// end of sphereAxis[d]

				}

				f << "sphereDist = { " << nDimsThis-1 << " 2 " ;

				for( int i = 0 ; i < nDimsThis-1 ; i++ )
					f << setw(12) << setprecision(8) << PNSSpoke[ns]->getSphereDist(i) << " " ;

				f << " } ; " << endl ; // end of sphereDist

			f << " } ; " << endl ;	// end of PNSSpoke[ns]

		}

		// ---------------------------  Eigenmodes  ------------------------------------

		f << "nEigenModes = " << nEigenmodes << " ; " << endl ;		

		f << "eigenValue = { " << nEigenmodes << " 2 " ;

		for( int i = 0 ; i < nEigenmodes ; i ++ ) 
			f << setw(12) << setprecision(8) << eigenValue[i] << " " ;

		f << " } ; " << endl ;		// end of eigenValue

		f << "eigenVectorLength = " << eigenVectorLength << " ; " << endl ;		

		for( int v = 0 ; v <  nEigenmodes ; v ++ ) {

			f << "eigenVector[" << v << "] = { " << eigenVectorLength << " 2 " ; 		

			for( int j = 0 ; j < eigenVectorLength ; j++ )
				f << setw(12) << setprecision(8) << *( eigenVector[v] + j ) << " " ;

			f << "} ;" << endl ;	// end of eigenVector[v]
		}
		
		// -----------------------------------------------------------------------------

	f << "}" ;		// end of CPNSStats

	f.close() ;

	return(1) ;

}

// private functions

double * M3DCPNSStats::convertSRepToEComp( M3DObject * obj, Vector3D & meanOfPDM ) {


	// --------------------------- input verification  --------------------------------

	if( obj == NULL ) {
		cout << "Error in M3DCPNSStats::convertSRepToEComp()! Null object input." << endl ;
		return( NULL ) ;
	}

	if( obj->getFigureCount() > 1 ) {
		cout << "Error in M3DCPNSStats::convertSRepToEComp()! Multi-figure objects not allowed" << endl ;
		return( NULL ) ;
	}

	if( obj->getFigurePtr(0)->getSpokeCount() != this->nSpokes ) {
		cout << "Error in M3DCPNSStats::convertSRepToEComp()! No. of spokes donot match" << endl ;
		cout << "This object does not belong to the class of objects with which this CPNS Stats were created !" << endl ;
		return( NULL ) ;
	}

	// ------------------------  this CPNS stats verification --------------------------

	if( this->nSpokes < 0 || this->nEigenmodes < 0 || this->eigenVectorLength < 0 ) {
		cout << "Error in M3DCPNSStats::convertSRepToEComp()! CPNSStats not properly initialized yet" << endl ;		
		return( NULL ) ;
	}

	// ------------------------------------ CONVERSION PROCESS ----------------------------------------

	// ---- Vectorize the object (figure) into a vector of [ positions ; spoke-dirs ; radii ] ----

	M3DQuadFigure * qFigThis = dynamic_cast <M3DQuadFigure *>(obj->getFigurePtr(0)) ;

	if( qFigThis == NULL ) {
		cout << "Error in M3DCPNSStats::convertSRepToEComp()! Does not work for non-quad figures" << endl ;		
		return( NULL ) ;
	}
	
	double * sRepVector = qFigThis->vectorize() ;	

	int nAtoms = qFigThis->getPrimitiveCount() ;

	// separate containers for positions, spoke directions and spoke radii

	double * positions = new double[ 3*nAtoms ] ;
	for( int i = 0 ; i < 3*nAtoms ; i++ )
		positions[i] = 0.0 ;

	vector <double> spokeDirs( 3 * nSpokes, 0.0 ) ;
	vector <double> spokeRadii( nSpokes, 0.0 ) ;

	for( int i = 0 ; i < nAtoms ; i++ ) {
		positions[ 3*i + 0 ] = sRepVector[ 3*i + 0 ] ;
		positions[ 3*i + 1 ] = sRepVector[ 3*i + 1 ] ;
		positions[ 3*i + 2 ] = sRepVector[ 3*i + 2 ] ;
	}

	for( int i = 0 ; i < nSpokes ; i ++ ) {

		spokeRadii[ i ] = sRepVector[ 3*nAtoms + i ] ;
		
		spokeDirs[ 3*i + 0 ] = sRepVector[ 3*nAtoms + nSpokes + 3*i + 0 ] ;
		spokeDirs[ 3*i + 1 ] = sRepVector[ 3*nAtoms + nSpokes + 3*i + 1 ] ;
		spokeDirs[ 3*i + 2 ] = sRepVector[ 3*nAtoms + nSpokes + 3*i + 2 ] ;
	}

	if( sRepVector != NULL ) 
		delete [] sRepVector ;
	sRepVector = NULL ;

	// container for converted positions, spoke directions, spoke radii

	double * ePositions = NULL ;

	vector <double> eSpokeDirs ;

	// ----------------  initialize outputVector  --------------------

	double * EComp = new double[ 3*nAtoms + 3*nSpokes ] ;

	// -------------------------- convert positions  ---------------------------

	// Find out the mean of the PDM 

	double meanX = 0.0, meanY =0.0 , meanZ =0.0 ;

	for( int i = 0 ; i < nAtoms ; i++ ) {
		meanX += positions[ 3*i + 0 ] ;
		meanY += positions[ 3*i + 1 ] ;
		meanZ += positions[ 3*i + 2 ] ;
	}
	
	meanX /= nAtoms ;
	meanY /= nAtoms ;
	meanZ /= nAtoms ;

	meanOfPDM.setX(meanX) ;
	meanOfPDM.setY(meanY) ;
	meanOfPDM.setZ(meanZ) ;

	// Subtract the mean of the PDM

	for( int i = 0 ; i < nAtoms ; i++ ) {		
			positions[ 3*i + 0 ] -= meanOfPDM.getX() ;
			positions[ 3*i + 1 ] -= meanOfPDM.getY() ;
			positions[ 3*i + 2 ] -= meanOfPDM.getZ() ;
	}

	// calculate scale of this shape PDM 

	double scaleShapeThis = 0.0 ;

	for( int i = 0 ; i < 3*nAtoms ; i ++ )
		scaleShapeThis += positions[i] * positions[i] ;

	scaleShapeThis = sqrt( scaleShapeThis ) ;

	// scale the shape PDM: positions = positions / scaleShapeThis ;

	for( int i = 0 ; i < 3*nAtoms ; i++ )
		positions[i] /= scaleShapeThis ;

	ePositions = this->PNSShape->convertS2E( positions ) ;

	if( positions != NULL )
		delete [] positions ;
	positions = NULL ;

	// error check for position conversion 
	if( ePositions == NULL ) {
		cout << "Error in M3DCPNSStats::convertSRepToEComp()! Error in converting positions (PDM)" << endl ;
		spokeDirs.clear() ;
		spokeRadii.clear() ;
		if( EComp != NULL )
			delete [] EComp ;
		return( NULL ) ;
	}

	// ------------  output vector: positions  ----------------------

	// outputPositions = [ scaleShape * ePositions ; scaleShape * log( scaleShapeThis / scaleShape ) ]

	for( int i = 0 ; i < 3*nAtoms-1 ; i++ )
		EComp[i] = scaleShape * ePositions[i] ;

	EComp[ 3*nAtoms-1 ] = scaleShape * log( scaleShapeThis / scaleShape ) ; 


	// ------------  output vector: spoke-radii  ----------------------

	// output spoke radii[i] = scalepoke[i] * log( scaleSpokeThis[i] / scaleSpoke[i]

	for( int i = 0 ; i < nSpokes ; i++ ) 
		EComp[ 3*nAtoms + i ] = scaleSpoke[i] * log( spokeRadii[i] / scaleSpoke[i] ) ;	
	


	// ----------------  convert spoke directions ------------------------------

	double * spokeThis = new double[3] ;

	double * eSpokeThis = NULL ;

	eSpokeDirs.resize( 2*nSpokes, 0.0 ) ;

	for( int ns = 0 ; ns < nSpokes ; ns ++ ) {

		for( int j = 0 ; j < 3 ; j++ )
			spokeThis[ j ] = spokeDirs[ 3*ns + j ] ;

		eSpokeThis = this->PNSSpoke[ns]->convertS2E( spokeThis ) ;

		// error check for spoke dir conversions
		if( eSpokeThis == NULL ) {
			cout << "Error in M3DCPNSStats::convertSRepToEComp()! Error in converting spoke direction[" << ns << "]" <<  endl ;
			spokeDirs.clear() ;
			spokeRadii.clear() ;
			if( EComp != NULL )
				delete [] EComp ;
			if( spokeThis != NULL )
				delete [] spokeThis ;
			return( NULL ) ;
		}		

		eSpokeDirs[ 2*ns + 0 ] = eSpokeThis[0] ;
		eSpokeDirs[ 2*ns + 1 ] = eSpokeThis[1] ;

		delete [] eSpokeThis ;
		eSpokeThis = NULL ;

	}
	
	// ------------  output vector: spoke-dirs  ----------------------

	// output spoke dirs = [ scalepoke[i] * eSpokeDirs[ 2*i : 2*i+1 ] ] 

	for( int i = 0 ; i < nSpokes ; i++ ) {
		EComp[ 3*nAtoms + nSpokes + 2*i + 0 ] = scaleSpoke[i] * eSpokeDirs[ 2*i + 0 ] ;
		EComp[ 3*nAtoms + nSpokes + 2*i + 1 ] = scaleSpoke[i] * eSpokeDirs[ 2*i + 1 ] ;
	}



	// -----------------  deallocate space used ----------------------

	if( ePositions!= NULL )
		delete [] ePositions ;

	if( spokeThis != NULL )
		delete [] spokeThis ;

	if( eSpokeThis != NULL )
		delete [] eSpokeThis ;

	spokeDirs.clear() ;
	spokeRadii.clear() ;

	return( EComp ) ;
	
}

M3DFigure * M3DCPNSStats::convertECompToSRepFigure(double *eComp, Vector3D meanOfPDM ) {

	// --------------------------- input verification  --------------------------------

	if( eComp == NULL ) {
		cout << "Error in M3DCPNSStats::convertECompToSRepFigure()! Null object input." << endl ;
		return( NULL ) ;
	}

	// ------------------------  this CPNS stats verification --------------------------

	if( this->nSpokes < 0 || this->nEigenmodes < 0 || this->eigenVectorLength < 0 ) {
		cout << "Error in M3DCPNSStats::convertECompToSRepFigure()! CPNSStats not properly initialized yet" << endl ;		
		return( NULL ) ;
	}

	// ------------------------------------ CONVERSION PROCESS ----------------------------------------	

	int nAtoms = this->nAtomCols * this->nAtomRows ;

	// -------------------------- recover positions  ---------------------------

	// ePosition[i] = eComp[i] / gamma_bar 

	double * ePositions = new double[ 3*nAtoms-1 ] ;

	double * positions = NULL ;
	
	for( int i = 0 ; i < 3*nAtoms-1 ; i ++ )
		ePositions[i] = eComp[i] / scaleShape ;

	double scaleShapeThis = scaleShape * exp( eComp[ 3*nAtoms-1 ] / scaleShape ) ;

	positions = this->PNSShape->convertE2S( ePositions ) ;

	delete [] ePositions ;
	ePositions = NULL ;

	if( positions == NULL ) {
		cout << "Error in M3DCPNSStats::convertECompToSRepFigure()! PNSShape->convertE2S() did not work properly " << endl ;				
		return( NULL ) ;
	}

	// scale the positions and add mean

	for( int i = 0 ; i < 3*nAtoms ; i++ )
		positions[i] *= scaleShapeThis ;

	for( int i = 0 ; i < nAtoms ; i++ ) {
		positions[ 3*i + 0 ] += meanOfPDM.getX() ;
		positions[ 3*i + 1 ] += meanOfPDM.getY() ;
		positions[ 3*i + 2 ] += meanOfPDM.getZ() ;
	}


	// -------------------------- recover radii  ---------------------------

	vector <double> spokeRadii( nSpokes, 0.0 ) ;

	for( int ns = 0 ; ns < nSpokes ; ns ++ )
		spokeRadii[ns] = scaleSpoke[ns] * exp( eComp[ 3*nAtoms + ns ] / scaleSpoke[ns] ) ;

	// -------------------------- recover spoke directions  ---------------------------

	vector <double> spokeDirs( 3*nSpokes, 0.0 ) ;

	double * eSpokeThis = new double[2] ;

	double * spokeThis = NULL ;

	for( int ns = 0 ; ns < nSpokes ; ns ++ ) {

		for( int j = 0 ; j < 2 ; j++ )
			eSpokeThis[ j ] = eComp[ 3*nAtoms + nSpokes + 2*ns + j ] ;

		spokeThis = this->PNSSpoke[ns]->convertE2S( eSpokeThis ) ;

		// error check for spoke dir conversions
		if( spokeThis == NULL ) {
			cout << "Error in M3DCPNSStats::convertECompToSRepFigure()! Error in converting spoke direction[" << ns << "]" <<  endl ;
			spokeDirs.clear() ;
			spokeRadii.clear() ;
			delete [] positions ;			
			delete [] eSpokeThis ;
			return( NULL ) ;
		}		

		for( int j = 0 ; j < 3 ; j++ )
			spokeDirs[ 3*ns + j ] = spokeThis[j] ;	

		delete [] spokeThis ;
		spokeThis = NULL ;
	}

	// --------------------  Pack everything into a CPNS data vector  --------------------------------

	double * cpnsDataVector = new double[ 3*nAtoms + 4*nSpokes ] ;

	for( int i = 0 ; i < 3*nAtoms ; i++ )
		cpnsDataVector[i] = positions[i] ;

	for( int i = 0 ; i < nSpokes ; i++ ) {
		cpnsDataVector[ 3*nAtoms + i ] = spokeRadii[i] ;
		for( int j = 0 ; j < 3 ; j ++ )
			cpnsDataVector[ 3*nAtoms + nSpokes + 3*i + j ] = spokeDirs[ 3*i + j ] ;
	}

	// construct a new quad figure

	M3DFigure * quadFig = new M3DQuadFigure( cpnsDataVector, nAtomRows, nAtomCols ) ;

	// Set properties of the new quad figure

	float color[3] = { 0.0, 1.0, 1.0 } ;
	int numLandmarks = 0 ;
	int smoothness = 50;
    unsigned short stackedImageMask = 52685; 

    quadFig->setName("default");
    quadFig->setColor(color);
	quadFig->setStackedImageMask( stackedImageMask ) ;

	// -----------------  deallocate space used ----------------------

	if( cpnsDataVector != NULL )
	delete [] cpnsDataVector ;

	if( eSpokeThis != NULL )
		delete [] eSpokeThis ;

	if( positions != NULL )
		delete [] positions ;

	spokeDirs.clear() ;
	spokeRadii.clear() ;

	return( quadFig ) ;
}


M3DObject * M3DCPNSStats::eigenmodeDeformMean( double * scores, M3DObject * origObject  ) {

	// input verification

	if( scores == NULL || origObject == NULL ) {
		cout << "Error in M3DCPNSStats::eigenmodeDeformMean() ! Null input." << endl ;
		return( NULL ) ;
	}

	// if any of the scores are outside [-3, +3], then bring them within that bracket

	for( int i = 0 ; i < nEigenmodes ; i++ ) {
		if( scores[i] < -3 )
			scores[i] = -3 ;
		if( scores[i] >   3 )
			scores[i] =   3 ;
	}

	for( int i = 0 ; i < nEigenmodes ; i++ )
		scores[i] *= sqrt( eigenValue[i] ) ;

	int nAtoms = nAtomRows * nAtomCols ;

	//double * EComp = new double [ 3*nAtoms + 3*nSpokes ] ;


	double * EComp = convertSRepToEComp( origObject, this->meanShape );

	//for( int i = 0 ; i < 3*(nAtoms + nSpokes) ; i++ )
	//	EComp[i] = 0.0 ;

	for( int ne = 0 ; ne < nEigenmodes ; ne++ ) {
		for( int i = 0 ; i < 3*(nAtoms + nSpokes) ; i++ )
			EComp[i] += scores[ne] * ( *(eigenVector[ne] + i) ) ;

	}	

	// copy the original object 

	M3DObject * resultObj = origObject->assign() ;

	M3DFigure * newFig = convertECompToSRepFigure( EComp, this->meanShape ) ;

	resultObj->setFigurePtr( 0, newFig ) ;

	delete [] EComp ;

	return( resultObj ) ;
}