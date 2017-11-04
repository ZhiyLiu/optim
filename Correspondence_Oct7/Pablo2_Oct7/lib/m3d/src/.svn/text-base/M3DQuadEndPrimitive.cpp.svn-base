#include "M3DPrimitive.h"
#include <iostream>

#include <typeinfo>

using namespace std;

// Constructor: Makes custom primitive.  
M3DQuadEndPrimitive::
M3DQuadEndPrimitive(const Vector3D& X, 
                    double R0, double R1, double REnd,
                    const Vector3D& U0, const Vector3D& U1,
                    const Vector3D& UEnd) :
M3DEndPrimitive(X, REnd),
M3DQuadPrimitive(X, R0, R1, REnd, U0, U1, UEnd),
M3DPrimitive(X)
{}

// Constructor: Makes custom primitive.  
//   x0..x2: Location 
//   Rr0..RrEnd: Radii of spokes
//   U0x..UEndz: Directions of the spokes 
M3DQuadEndPrimitive::
M3DQuadEndPrimitive(double x0, double x1, double x2, 
                    double R0, double R1, double REnd,
                    double U0x, double U0y, double U0z,
                    double U1x, double U1y, double U1z,
                    double UEndx, double UEndy, double UEndz) 
                    : M3DEndPrimitive(x0, x1, x2, REnd), 
                    M3DQuadPrimitive(x0, x1, x2, R0, R1, REnd, U0x, U0y, U0z,
                    U1x, U1y, U1z, UEndx, UEndy, UEndz),
                    M3DPrimitive(x0, x1, x2)
{}

// CHECK SEMANTICS OF ELONGATION
M3DQuadEndPrimitive::
M3DQuadEndPrimitive(double x0, double x1, double x2, double R,
                    const Quat &Q, double _theta, double _elongation)
                    : M3DEndPrimitive( x0, x1, x2, 1.0 ),
                    M3DQuadPrimitive( x0, x1, x2, R, Q, _theta),
                    M3DPrimitive(x0, x1, x2) 
{
    setElongation(_elongation);
}

M3DQuadEndPrimitive::M3DQuadEndPrimitive(const Vector3D &X, double R, const Quat &Q,
                                         double _theta, double _elongation)
                                         : M3DEndPrimitive( X, 1.0 ),
                                         M3DQuadPrimitive( X, R, Q, _theta),
                                         M3DPrimitive(X) 
{
    setElongation(_elongation);
}


/**
* A virtual copy constructor. Since C++ does not support
* a virtual copy constructor, all derived (non-abstract) classes
* need to implement this method. Replace M3DPrimitive by the
* correct class name.
*/
M3DPrimitive * M3DQuadEndPrimitive::copyPtr() const
{
    return new M3DQuadEndPrimitive(*this);
}

// The assignment operator series to take care of
// c++/compiler bugs.
M3DPrimitive & M3DQuadEndPrimitive::operator = (const M3DPrimitive & prim)
{
    assert( typeid(*this) == typeid(prim) );
    const M3DQuadPrimitive& qprim = dynamic_cast<const M3DQuadPrimitive&>(prim);
    const M3DEndPrimitive& eprim = dynamic_cast<const M3DEndPrimitive&>(prim);
    M3DPrimitive::copy(prim);
    M3DQuadPrimitive::copy(qprim);
    M3DEndPrimitive::copy(eprim);

#ifdef DEBUG
    cout << "M3DQuadEndPrimitive::operator =() called\n";
#endif

    return *this;
}

// Print method.  This doesn't call M3DQuadPrimitive::Print
// because that wouldn't format it right.
void M3DQuadEndPrimitive::print(std::ostream & out, char * prefix, bool marked) const
{
    M3DPrimitive::print(out, prefix, marked);

    std::string p;
    if (prefix) p = prefix;

    out << p << "r0, u0 = " << r0 << ", " << u0 << '\n';
    out << p << "r1, u1 = " << r1 << ", " << u1 << '\n';
    out << p << "rEnd, uEnd = " << rEnd << ", " << uEnd << '\n';
}

void M3DQuadEndPrimitive::setR(double r)
{
    M3DQuadPrimitive::setR(r);
    M3DEndPrimitive::setR(r);
}

void M3DQuadEndPrimitive::scaleBy(double mag)
{
    M3DQuadPrimitive::scaleBy(mag);
    rEnd *= mag;
}

void M3DQuadEndPrimitive::scaleBy(double mag0, double mag1, double mag2)
{
	M3DQuadPrimitive::scaleSpokeBy( 0, mag0 ) ;
	M3DQuadPrimitive::scaleSpokeBy( 1, mag1 ) ;

    rEnd *= mag2 ;
}

