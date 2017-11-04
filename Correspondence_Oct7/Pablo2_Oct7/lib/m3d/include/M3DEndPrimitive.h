#ifndef M3DENDPRIMITIVE_H
#define M3DENDPRIMITIVE_H

#include "M3DPrimitive.h"

const char M3D_END_PRIMITIVE_STR[] = "EndPrimitive";

class M3DEndPrimitive : public virtual M3DPrimitive
{
public:

    M3DEndPrimitive() : M3DPrimitive()
    {
        rEnd = 1.0;
    }

    M3DEndPrimitive(double x0, double x1, double x2, double _rEnd) :
        M3DPrimitive(x0, x1, x2)
    {
        rEnd = _rEnd;
    }

    M3DEndPrimitive(const Vector3D &X, double _rEnd) :
        M3DPrimitive(X)
    {
        rEnd = _rEnd;
    }

    M3DEndPrimitive(const M3DEndPrimitive & prim)
        : M3DPrimitive(prim)
    {
        rEnd = prim.rEnd;
    }

 
    virtual ~M3DEndPrimitive() {}
    /**
     * This function returns whether this object is
     * an end primitive or a standard primitive
     * It is advised to use this function instead of
     * dynamic_cast for code simplicity, elegance and
     * better performance
     *
     * @return M3D_STANDARD_PRIMITIVE
     */
    virtual M3DPrimitiveType type() const { return M3D_END_PRIMITIVE; }

    virtual M3DPrimitive * copyPtr() const = 0;

    virtual void print(std::ostream & out, char * prefix, bool marked) const;

    virtual void setR(double r);

    double getREnd() const { return rEnd; }
    void setREnd(double _rEnd) { rEnd = _rEnd; }

    double getElongation() const { return rEnd / getR(); }
    void setElongation(double e) { rEnd = getR() * e; }

    virtual M3DPrimitive & operator = (const M3DPrimitive & prim);
    virtual M3DPrimitive & operator = (const M3DEndPrimitive& figure ) {
        return operator=((const M3DPrimitive&) figure);
    }

    virtual const char* name() const {
        return M3D_END_PRIMITIVE_STR;
    }

protected:
    double rEnd;

    // This member is deprecated, phase it out.
    double elongation;

    /**
     * Copy information relevant only to end primitives.
     */
    void copy(const M3DEndPrimitive& prim);
};

#include "M3DTubeEndPrimitive.h"
#include "M3DQuadEndPrimitive.h"

#endif

