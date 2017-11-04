/*---------------------------------------------------------------------------
Class definitions for curvature at a point of the surface in R^3

Author:			Xiaojie Zhao, D Goswami
Version:		1.0
Last Modified:	03/17/2011
---------------------------------------------------------------------------*/

#ifndef CURVATURE_H
#define CURVATURE_H

#include <math.h>

#include "Vector3D.h"
#include "Vector2D.h"

class Curvature {

	public:		

		// principal directions: p1, p2

		Vector3D p1 ;		
		Vector3D p2 ;
	    
		// principal curvatures: kappa1, kappa2

		double kappa1 ;
		double kappa2 ;
	    
		// gradient of the principal curvatures

		// Vector3D gradKappa1 ;
		// Vector3D gradKappa2 ;
	    
		Curvature( ) ;

		Curvature( Vector3D pDir1, Vector3D pDir2, double k1, double k2 ) ;
	    
		~ Curvature( ) ;

} ;

#endif // CURVATURE_H