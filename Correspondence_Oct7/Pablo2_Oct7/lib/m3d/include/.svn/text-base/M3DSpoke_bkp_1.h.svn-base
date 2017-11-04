#ifndef _M3D_SPOKE
#define _M3D_SPOKE

#include "Vector3D.h"

class M3DSpoke {

	protected:

		Vector3D X ;	// the position of the tail of the spoke 

		Vector3D U ;	// the unit vector, that is direction of the spoke

		double  r ;		// the length of the spoke

	public:

		M3DSpoke() {
		}

		M3DSpoke( Vector3D _X, Vector3D _U, double _r ) ;

		M3DSpoke( double x, double y, double z, double uX, double uY, double uZ, double _r ) ;

		virtual Vector3D getX() { return X; }

		virtual Vector3D getU() { return U; } 

		virtual double getR() { return r; } 

		virtual Vector3D getS() { return (r * U); } // Returns the scaled spoke direction
 
		virtual Vector3D getB() { return ( X + (r * U) ); } // Returns implied boundary point
};


#endif // _M3D_SPOKE