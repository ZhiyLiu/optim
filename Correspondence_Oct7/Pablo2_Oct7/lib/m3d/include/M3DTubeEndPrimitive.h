#ifndef M3DTUBEENDPRIMITIVE_H
#define M3DTUBEENDPRIMITIVE_H

class M3DTubeEndPrimitive : public virtual M3DEndPrimitive, public virtual M3DTubePrimitive
{
  public:
	M3DTubeEndPrimitive() {
	}

    M3DTubeEndPrimitive(double x0, double x1, double x2,
		double rEnd, double d,
		double* radii, const Vector3D& U0, const Vector3D& U2, int numberOfSpokes) :
		M3DEndPrimitive(x0,x1,x2,rEnd),
		M3DTubePrimitive(x0,x1,x2, d, radii, U0, U2, numberOfSpokes),
		M3DPrimitive(x0,x1,x2) {
	}

	// Constructor: Makes custom primitive.
    M3DTubeEndPrimitive(const Vector3D &X,
		double rEnd, double d,
		double* radii, const Vector3D& U0, const Vector3D& U2, int numberOfSpokes) :
		M3DEndPrimitive(X,rEnd),
		M3DTubePrimitive(X, d, radii, U0, U2, numberOfSpokes),
		M3DPrimitive(X) {
	}

	/**
	 * Converts this atom from a lie group representation
	 * to a symmetric space representation which is used by
	 * the statistics code. The object is allocated
	 * by the function. So it should be deleted by the caller.
	 *
	 * @param none
	 * @return pointer to a M3DSymPrimitive
	 * @see M3DSymPrimitive
	 * @see M3DSymPrimitive::symAtomToAtom
	 */
	//virtual M3DSymPrimitive* atomToSymAtom() const;

	/**
	 * A virtual copy constructor. Since C++ does not support
	 * a virtual copy constructor, all derived (non-abstract) classes
	 * need to implement this method. Replace M3DPrimitive by the
	 * correct class name.
	 */
    virtual M3DPrimitive * copyPtr() const;

    virtual void scaleBy(double mag) {
		rEnd *= mag;
		M3DTubePrimitive::scaleBy(mag);
	}

    /**
     * Returns the extended B vector. Note, that
     * if theta is greater than 90 degrees, then this
     * vector will point in the reverse direction of B.
	 * Note the reverse direction difference. You cannot
	 * just multiply B by 90 degrees to get the extended B vector.
     */
    virtual Vector3D getExtendedB() const
    {
        Vector3D b  = getB();
        if( d < 0 ) {
            b   = -b;
        }
        return b * rEnd;
    }

    virtual void setR(double r);

	// The assignment operator series to take care of
	// c++/compiler bugs.
    virtual M3DPrimitive & operator = (const M3DPrimitive & prim);
	M3DPrimitive& operator = ( const M3DTubeEndPrimitive& prim ) {
		return operator = ((const M3DPrimitive&) prim);
	}

};

#endif

