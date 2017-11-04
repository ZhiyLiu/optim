#include "M3DEndPrimitive.h"
#include <assert.h>

#include <typeinfo>

M3DPrimitive & M3DEndPrimitive::operator = (const M3DPrimitive & unknown_prim)
{
    assert( typeid(*this) == typeid(unknown_prim) );
    const M3DEndPrimitive& prim = 
        dynamic_cast<const M3DEndPrimitive&>(unknown_prim);

#ifdef DEBUG
    cout << "M3DEndPrimitive::operator=()" << endl;
#endif

    M3DPrimitive::copy(static_cast<const M3DPrimitive&>(prim));
    copy(prim);

    return (*this);
}

void M3DEndPrimitive::copy( const M3DEndPrimitive& prim )
{
    rEnd = prim.rEnd;
}

void M3DEndPrimitive::print(std::ostream & out, char * prefix, bool marked) const
{
    char * p;
    char str = '\0';

    if (prefix == NULL)
        p = &str;
    else
        p = prefix;

    M3DPrimitive::print(out, prefix, marked);
    out << p << "rEnd = " << rEnd << '\n';
}

void M3DEndPrimitive::setR(double r)
{
    rEnd *= r / getR();
}

