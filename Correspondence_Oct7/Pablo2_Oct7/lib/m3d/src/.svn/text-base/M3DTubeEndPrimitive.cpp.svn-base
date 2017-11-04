#include "M3DPrimitive.h"
#include <iostream>

#include <typeinfo>

using namespace std;

/**
 * A virtual copy constructor. Since C++ does not support
 * a virtual copy constructor, all derived (non-abstract) classes
 * need to implement this method. Replace M3DPrimitive by the
 * correct class name.
 */
M3DPrimitive * M3DTubeEndPrimitive::copyPtr() const
{
	return new M3DTubeEndPrimitive(*this);
}

// The assignment operator series to take care of
// c++/compiler bugs.
M3DPrimitive & M3DTubeEndPrimitive::operator = (const M3DPrimitive & prim)
{
    assert( typeid(*this) == typeid(prim) );
    const M3DTubePrimitive& tprim = dynamic_cast<const M3DTubePrimitive&>(prim);
    const M3DEndPrimitive& eprim = dynamic_cast<const M3DEndPrimitive&>(prim);
    M3DPrimitive::copy(prim);
    M3DTubePrimitive::copy(tprim);
    M3DEndPrimitive::copy(eprim);
	
#ifdef DEBUG
    cout << "M3DTubeEndPrimitive::operator =() called\n";
#endif

    return *this;
}

void M3DTubeEndPrimitive::setR(double r)
{
    M3DTubePrimitive::setR(r);
    M3DEndPrimitive::setR(r);
}
