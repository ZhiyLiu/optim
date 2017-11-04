#include "M3DSpoke.h"

M3DSpoke::M3DSpoke( Vector3D _X, Vector3D _U, double _r ) {

	X = _X ;
	U = _U ;
	r = _r ;			
}

M3DSpoke::M3DSpoke( double x, double y, double z, double uX, double uY, double uZ, double _r ) {

	X = Vector3D( x, y, z ) ;
	U = Vector3D( uX, uY, uZ ) ;
	r = _r ;

}


// Liyun Jul 29, 2014
void M3DSpoke::setX(Vector3D _X){
    this->X = _X;
}
void M3DSpoke::setU(Vector3D _U){
    this->U = _U;
}
void M3DSpoke::setR(double _r){
    this->r = _r;
}
