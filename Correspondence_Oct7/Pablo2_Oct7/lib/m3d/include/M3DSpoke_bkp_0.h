#ifndef _M3D_SPOKE
#define _M3D_SPOKE

#include "Vector3D.h"

class M3DSpoke {

	private:

		Vector3D X ;	// the position of the tail of the spoke 

		Vector3D U ;	// the unit vector, that is direction of the spoke

		double  r ;		// the length of the spoke

	public:

		M3DSpoke() {
		}

		M3DSpoke( Vector3D _X, Vector3D _U, double _r ) ;

		M3DSpoke( double x, double y, double z, double uX, double uY, double uZ, double _r ) ;
} ;


#endif // _M3D_SPOKE