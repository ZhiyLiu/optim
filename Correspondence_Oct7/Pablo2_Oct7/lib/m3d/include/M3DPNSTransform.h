#ifndef _M3D_PNS_TRANSFORM_H
#define _M3D_PNS_TRANSFORM_H

#include <iostream>
using namespace std ;

#include <vector>

/*
	Author: Dibyendusekhar Goswami, July 2011

	Description: A class to hold a PNS transformation.

	Theory: A PNSTransformation takes data from an n-dimensional spherical space to an n-dimensional Euclidean space.
			A description of this can be found in CPNS_sreps.PDF or Notes_on_CPNS.pdf
			This is known as h (hS, or hSpoke) in these documents

*/

class M3DPNSTransform {

public:

	M3DPNSTransform() ;

	M3DPNSTransform( int _nDims ) ;

	M3DPNSTransform( int _nDims, std::vector <double*> & _sphereAxis, double * _sphereDist ) ;

	~ M3DPNSTransform() ;

	// copy from another PNSTransform

	void copy( M3DPNSTransform * that ) ;

	// assign this to another PNSTransform - make a copy

	M3DPNSTransform * assign() ;

	// -------------------  get functions --------------------------------

	// query the dimension of the sphere
	int getNDimensions() ;

	// query sphereDist (entire)
	double * getSphereDist() ;

	// query sphereDist (for a particular sphere)
	double getSphereDist( int n ) ;

	// query sphereAxis (one value )
	double getSphereAxisVal( int nSphere, int index ) ;

	//// query sphereAxis( entire )
	//bool getSphereAxis( vector <double *> & o_sphereAxis ) ;

	//// query sphereAxis( particular sphere )
	//bool getSphereAxis( int nSphere, double * o_shereAxis, int& o_nDims ) ;
	

	// -------------------  set functions --------------------------------

	// return value for the set methods with bool output { 1: success, 0: failure }

	// set the dimension of the largest sphere - sphereAxis and sphereDist will be set to NULL
	void setNDimensions( int _nDims ) ;

	// set sphereDist (entire) - checks whether the input _nDims is the same as nDims stored inside	
	bool setSphereDist( double * _sphereDist, int lenSphereDist ) ;

	// set sphereAxis (entire)
	bool setSphereAxis( std::vector <double*> & _sphereAxis ) ;

	// set sphereAxis (particular) - for sphere number nSphere (starts from 0, 0 is highest dimensional sphere)
	bool setSphereAxis( int nSphere, double * sphereAxisN, int lenSphereAxisN ) ;

	// -------------------  space conversion functions --------------------------------

	// convert from a d-dimensional spherical space to d-dimensional Euclidean space (PNS residuals)

	double * convertS2E( double * S ) ;

	// convert from a d-dimensional Euclidean space (PNS residuals) to d-dimensional spherical space
	
	double * convertE2S( double * E ) ;

	// --------------------------------------------------------------------------------

	// print to output screen

	void printToScreen() ;

protected:

	int nDims ;							// dimensionality of the largest sphere - d

	std::vector <double *> sphereAxis ;	// the axes of the successive spheres: v1, v2, ... vd

	double * sphereDist ;				// the dist of each sphere from its parent: r1, r2, ... rd-1

	// --------------------   private functions  --------------------------------

	// called from within convertS2E()

	// get a particular sphere axis
	bool getSphereAxis( std::vector <double> & sphereAxisN, int N ) ;

	// get the rotation matrix for rotating a vector v to the north pole
	double * * getRotationMatrixToNorthPole( std::vector <double> & v ) ;



} ;

#endif // _M3D_PNS_TRANSFORM_H