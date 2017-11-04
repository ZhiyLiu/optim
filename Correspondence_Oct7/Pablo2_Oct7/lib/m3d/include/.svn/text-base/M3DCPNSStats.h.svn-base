#ifndef _M3D_CPNS_STATS_H
#define _M3D_CPNS_STATS_H

#include <iostream>
using namespace std ;

#include <vector>
#include "Vector3D.h"

#include "M3DPNSTransform.h"
# include "M3DObject.h"

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

class M3DCPNSStats {

public:

	M3DCPNSStats() ;

	M3DCPNSStats( int _nSpokes, int _nEigenmodes, int _eigenVectorLength, int _nAtomRows, int _nAtomCols ) ;

	~M3DCPNSStats() ;


	bool read() ;			// to be written - currently read through the M3DObjectFile class
	bool write() ;			// to be written

	bool writeToFile( const char* fileName ) ;	// to be written to a text file
												// mainly for verification and debugging
	// get functions 

	// set functions ------------------------------------------------

	bool setScaleShape( double _scaleShape ) ;

	bool setMeanShape( double * meanShape ) ;

	bool setScaleSpoke( double * _scaleSpoke ) ;

	bool setPNSShape( M3DPNSTransform * _PNSShape ) ;

	// set PNSSpoke (entire)
	bool setPNSSpoke( vector<M3DPNSTransform *> & _PNSSpoke ) ;

	// set PNSSpoke (particular spoke) - still to be written
	bool setPNSSpoke( int nSpokeThis, M3DPNSTransform * PNSSpokeThis ) ;

	bool setEigenModes( double* _eigenValue, vector <double *> & _eigenVector ) ;	

	// get functions ------------------------------------------------

	int getNEigenModes() ;

	int getNSpokes() ;

	int getEigenVectorLength() ;

	M3DPNSTransform * getPNSShape() ;

	Vector3D getMeanShape() ;

	// get one particular spoke's PNS transform
	M3DPNSTransform * getPNSSpoke( int N ) ;

	// get a particular spoke's scaling factor (r_bar[N]): returns -1 for invalid spoke number
	double getScaleSpoke( int N ) ;

	// get a particular eigenvalue
	double getEigenValue( int N ) ;

	// --------------------------------------------------------------

	// For eigen mode deformations (for the mean model only)
	M3DObject * eigenmodeDeformMean( double * scores, M3DObject * origObject ) ;

	// For eigen mode deformations (general)
	M3DObject * eigenmodeDeform( double * scores, M3DObject * inModel ) ;

	// For mahalanobis distance between this object and the CPNS mean
	double getMahanobisDistance( M3DObject * inModel ) ;

	// private functions

	// convert an S-rep to EComp (composite Euclidean space for CPNS analysis)
	double * convertSRepToEComp( M3DObject * obj, Vector3D & meanOfPDM ) ;

	// convert from EComp to an S-rep
	M3DFigure * convertECompToSRepFigure( double * eComp, Vector3D meanOfPDM ) ;

protected:

	M3DPNSTransform *	PNSShape ;	// tranformation for the shape PDM (hub positions)

	Vector3D	meanShape ;			// mean of the entire PDM (x,y,z)

	double		scaleShape ;		// scale parameter for shape PDM (gamma_bar)

	int			nAtomRows ;			// no. of rows in the atom quad figure

	int			nAtomCols ;			// no. of cols in the atom quad figure

	int			nSpokes ;

	vector <M3DPNSTransform *>	PNSSpoke ;	// PNSSpoke[i] = PNS transformation for spoke-i

	double *	scaleSpoke ;		// an array of scaling factors for each spoke: scaleSpoke[i] = r_bar(i)

	int			nEigenmodes ;

	int			eigenVectorLength ;

	vector <double *>	eigenVector ;	// an array of eigenVectors: eigenVector[i] = pointer to i-th eigenVector

	double *			eigenValue ;		// eigenValue[i] = i-th eigenvalue

} ;

# endif // _M3D_CPNS_STATS_H