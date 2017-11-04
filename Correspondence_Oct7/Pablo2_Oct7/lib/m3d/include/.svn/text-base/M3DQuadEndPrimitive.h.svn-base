#ifndef M3DQUADENDPRIMITIVE_H
#define M3DQUADENDPRIMITIVE_H

class M3DQuadEndPrimitive : public virtual M3DEndPrimitive, public virtual M3DQuadPrimitive
{
public:
    M3DQuadEndPrimitive() 
        : M3DQuadPrimitive(), M3DEndPrimitive(), M3DPrimitive() {}

    // Constructor: Makes custom primitive.  
    M3DQuadEndPrimitive(const Vector3D& X, 
                        double R0, double R1, double REnd,
                        const Vector3D& U0, const Vector3D& U1,
                        const Vector3D& UEnd);

    // Constructor: Makes custom primitive.  
    //   x0..x2: Location 
    //   Rr0..RrEnd: Radii of spokes
    //   U0x..UEndz: Directions of the spokes 
    M3DQuadEndPrimitive(double x0, double x1, double x2, 
                        double R0, double R1, double REnd,
                        double U0x, double U0y, double U0z,
                        double U1x, double U1y, double U1z,
                        double UEndx, double UEndy, double UEndz);

    M3DQuadEndPrimitive(double x0, double x1, double x2, double R,
                        const Quat &Q, double _theta, double _elongation = 1.0);

    M3DQuadEndPrimitive(const Vector3D &X, double R, const Quat &Q,
                        double _theta, double _elongation = 1.0);

    // Assignment
    // Non-virtual operator is needed to prevent implicit operator= 
    virtual M3DPrimitive & operator=(const M3DPrimitive & prim);
    M3DPrimitive& operator=( const M3DQuadEndPrimitive& prim ) {
        return operator=(static_cast<const M3DPrimitive&>(prim));
    }

    virtual void scaleBy(double mag);

	// dibyendu

	virtual void scaleBy(double mag0, double mag1, double mag3);



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

    virtual void print(std::ostream & out = std::cout, char * prefix = NULL,
                       bool marked = false) const;


    virtual void setR(double R);

    // This is retained for compatibility with old code.
    // Returns the "extended B vector."  Note: Formerly, the
    // extended B vector defined the crest.  Now, YEnd serves 
    // that function.  In general, YEnd may not bisect Y0 and Y1,
    // but B and ExtendedB always will.
    virtual Vector3D getExtendedB() const { return getB() * rEnd; }

    virtual Vector3D getYEnd() const { return uEnd * rEnd; }

};

#endif

