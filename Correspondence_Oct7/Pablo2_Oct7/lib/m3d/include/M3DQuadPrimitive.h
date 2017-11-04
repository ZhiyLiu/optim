#ifndef M3DQUADPRIMITIVE_H
#define M3DQUADPRIMITIVE_H

#include "M3DPrimitive.h"

#include <iostream>
using namespace std;

class M3DQuadPrimitive : public virtual M3DPrimitive
{
public:
    // Constructor: Default makes "stock" primitive.
    M3DQuadPrimitive();

    // Constructor: Makes custom primitive.  
    M3DQuadPrimitive(const Vector3D& X, 
                     double R0, double R1, double REnd,
                     const Vector3D& U0, const Vector3D& U1, const Vector3D& UEnd);

    // Constructor: Makes custom primitive.  
    //   x0..x2: Location 
    //   Rr0..RrEnd: Radii of spokes
    //   U0x..UEndz: Directions of the spokes 
    M3DQuadPrimitive(double x0, double x1, double x2, 
                     double R0, double R1, double REnd,
                     double U0x, double U0y, double U0z,
                     double U1x, double U1y, double U1z,
                     double UEndx, double UEndy, double UEndz, double deltaU0, double deltaU1, double deltaU2, double deltaV0, double deltaV1, double deltaV2);  //Liyun add deltaU,deltaV

    // Constructor: Makes custom primitive.
    M3DQuadPrimitive(double x0, double x1, double x2,
                     double R, const Quat &Q, double _theta);

    // Constructor: Makes custom primitive.
    M3DQuadPrimitive(const Vector3D &X, double R, const Quat &Q, double _theta);

    // Constructor: copy
    M3DQuadPrimitive(const M3DQuadPrimitive & prim);

    // Non-virtual operator is needed to prevent implicit operator= 
    virtual M3DPrimitive & operator=(const M3DPrimitive & prim);
    M3DPrimitive& operator=( const M3DQuadPrimitive& prim ) {
        return operator=(static_cast<const M3DPrimitive&>(prim));
    }

    // After calling this->setR(someR), a call to getR should yield someR,
    // but the proportions of the lengths of the spokes should stay the same
    virtual void setR(double r);

	// This version of setR allows to set the length of every spoke independently
	virtual void setR( int spokeId, double r ) ;

    virtual void setQ(const Quat & Q);
    virtual void setTheta(double _theta);

	virtual void setU0( Vector3D _u0 ) ;
	virtual void setU1( Vector3D _u1 ) ;

	virtual void setUEnd( Vector3D _uEnd ) ;

    double getR0() const { check_r0(); return r0; }
    double getR1() const { return r1; }

    virtual void setR0( double _r0 ) ;//add by liyun Oct 23, 2013
    virtual void setR1( double _r1 ) ;//add by liyun Oct 23, 2013

    virtual Vector3D getU0() const { return u0; }
    virtual Vector3D getU1() const { return u1; }
    // Synonyms for the above two functions; remove?
    virtual Vector3D getNormalizedY0() const { return u0; }
    virtual Vector3D getNormalizedY1() const { return u1; }

    virtual Vector3D getUEnd() const { return uEnd; }

    virtual Vector3D getY0() const { return u0 * r0; }
    virtual Vector3D getY1() const { return u1 * r1; }
    virtual Vector3D getYEnd() const { return uEnd * getR(); }

    // Calculate a common radius for all spokes
    virtual double getR() const;

	// Calculate the radius for a particular spoke
	// virtual double getR( int spokeId ) ;

    // Return a quaternion describing the orientation
    virtual Quat getQ() const;

    // Return the angle between a representative spoke
    // and the bisector
    virtual double getTheta() const;

    // Return B, the normalized bisector of U0 and U1.  B is 
    // calculated directly from U0 and U1.  If this object 
    // is an m-rep, then B should equal UEnd.
    virtual Vector3D getB() const;

    // This is retained for compatibility with old code.
    // Returns the "extended B vector."  Note: Formerly, the
    // extended B vector defined the crest.  Now, UEnd serves 
    // that function.  In general, YEnd might not bisect Y0 and Y1,
    // but B and ExtendedB always will.
    // FIXME: This function is wrong. It's sometimes in the
    // opposite direction of getB().
    virtual Vector3D getExtendedB() const { return getB() * getR(); }

    // N is a unit vector that, in a slab, is normal to the 
    // medial sheet and on the same side as U1.
    virtual Vector3D getN() const;

    virtual void rotateBy(Quat dQ);
    virtual void scaleBy(double mag);

	// dibyendu
	// Function to dilate an atom - same effect as scaling all the spokes by a factor

	virtual void dilate( const double dilationFactor ) ;

	// dibyendu
	virtual void scaleSpokeBy( int spokeId, double mag ) ;

	// dibyendu
	virtual void rotateSpokeBy( int spokeId, double delTheta, double delPhi ) ;

    static M3DQuadPrimitive* readPrimitive(Registry& registry, 
                                           const char * regStr, ...);

    virtual void writePrimitive(Registry& registry,
                                const char * regStr, ...) const;

    /**
     * A virtual copy constructor. Since C++ does not support
     * a virtual copy constructor, all derived (non-abstract) classes
     * need to implement this method. Replace M3DPrimitive by the
     * correct class name.
     */
    virtual M3DPrimitive * copyPtr() const {
        return new M3DQuadPrimitive(*this);
    }

    virtual void print(std::ostream & out = std::cout, char * prefix = NULL,
                       bool marked = false) const;

    /**
     * Interpolates between two atoms. Interpolates between
     * m1 and m2 by t in [0,1] and modifies the current object
     * to match the result. Derived classes need to provide
     * an implementation for this method.
     *
     * @param	m1	One of the primitives for interpolation
     * @param	m2	The other primitive for the interpolation
     * @param	t	The weight given to m1/m2 for interpolation
     * @return	true/false depending upon whether the operation
     * succeeded or failed.
     */
    virtual bool atomInterp( double t, const M3DPrimitive* m1,
                             const M3DPrimitive* m2);

    /**
     * Returns the distance^2 between two atoms. Derived classes
     * need to provide an implementation for this method.
     * @param	m	The second atom 
     * @param	radius	Some radius normalization(?) used
     *					for geodesic distance computation.
     * @return	double, geodesic distance^2.
     */
    virtual double atomSquareDistance( const M3DPrimitive* m, const double* radius = NULL ) const;

    /**
     * Sets the current atom to the average of the atoms
     * passed in. Note, that if the type on which this
     * function is called is an M3DEndPrimitive, then the
     * average will take the elongation part of any end
     * primitives passed into the array into account,
     * else the elongation part is simply ignored.
     *
     * @param	mNum	number of atoms to take average of.
     * @param	mList	array of pointers to the atoms
     * @param	weights	weights, if a weighted average is
     *					desired, NULL
     * @return	true/false depending upon success (currently
     *			returns true all the time)
     */
    virtual bool atomAverage( const int mNum, M3DPrimitive** mList,
                              const double* weights = NULL );

    /**
     * Converts the current primitive to a primitive in
     * symmetric space representation.
     * Note, this operation carries over only the basic
     * elements of geometry and not stuff such as hinge et. al
     *
     * @param	none
     * @return	pointer to QuadSymPrimitive object
     *			*NOTE* the callee is responsible for clearing memory.
     */
    virtual SymPrimitive* convert2Sym() const {
        return new QuadSymPrimitive(this);
    }
	
    virtual bool vec2Atom(const Vector & vecPrim){
        QuadSymPrimitive symPrim ;
        symPrim.vec2Atom(vecPrim);
        return( symPrim.convert2Lie(this));
    }

    // DEBUGGING
    Vector3D dgetExtendedB() const;
    Vector3D dgetB() const;
    Vector3D dgetY0() const;
    Vector3D dgetY1() const;
    Vector3D dgetN() const;
    Vector3D dgetBPerp() const;
    Vector3D dgetNormalizedY0() const;
    Vector3D dgetNormalizedY1() const;


    //---------------------Liyun
    virtual double getDeltaU0() const { return deltaU0; }
    virtual double getDeltaU1() const { return deltaU1; }
    virtual double getDeltaU2() const { return deltaU2; }
    virtual double getDeltaV0() const { return deltaV0; }
    virtual double getDeltaV1() const { return deltaV1; }
    virtual double getDeltaV2() const { return deltaV2; }

    virtual void setDeltaU0( double _du0 );
    virtual void setDeltaU1( double _du1 );
    virtual void setDeltaU2( double _du2 );
    virtual void setDeltaV0( double _dv0 );
    virtual void setDeltaV1( double _dv1 );
    virtual void setDeltaV2( double _dv2 );


    //------------------------------


protected:

    //---------------------Liyun
    double deltaU0;
    double deltaU1;
    double deltaU2;
    double deltaV0;
    double deltaV1;
    double deltaV2;
    //------------------------------

    // New s-rep members.  These can be calculated from the old
    // representation, and calculation can go the other way if opposing
    // spokes are consistently the same length.
    double r0;
    double r1;
    //double rEnd; //Add by Liyun Tu, Mar 6, 2014. rEnd has already been defined in M3DEndPrimitive.h

    Vector3D u0;
    Vector3D u1;
    Vector3D uEnd;

    void copy(const M3DQuadPrimitive & prim);

    void setQTheta(Quat q, double theta);

    // This function is for debugging purposes.  There was a problem with 
    // radius values becoming negative.  That problem has been solved by
    // appropriately choosing parameters, but I am leaving the testing code
    // in for now.  Eventually it should be removed. -- MSF
    void check_r0() const
    {
#ifdef _DEBUG

		// dibyendu: commenting out the check_r0 ... screen output

        // cout << "check_r0 at line " << __LINE__ << endl;
        if (r0 < 0) {
            std::cout << "r0 = " << r0 << std::endl;
        }
#endif
    }

    /**
     * This class is a placeholder for the symmetric space representation
     * for quads.  It is meant to be used internally only by quads.  Note
     * that the "protected" declaration has no effect on this class; you
     * can access it from anywhere as QuadPrimitive::QuadSymPrimitive.
     * 
     * Foskey: This will likely eventually be removed as the main
     * representation will be s-reps.
     */
    class QuadSymPrimitive : public SymPrimitive {
    public:
        /**
         * Default constructor.
         */
        QuadSymPrimitive() : x(0.0, 0.0, 0.0), r(1.0), n0(0.0,-1.0,0.0),
                             n1(0.0, 1.0, 0.0), elongation(1.0), ignorePrimType(false) {
        }
        /**
         * Constructor for creating a Symmetric primitive from a quad primitive.
         * It handles both ordinary quad primitives as well as end primitives.
         *
         * @param	prim	A quad primitive.
         */
        QuadSymPrimitive( const M3DQuadPrimitive* prim );

        /**
         * Computes the Exponential map of the current primitive.
         * @param	none
         * @result	saved in the object on which the function was called.
         */
        virtual void Exp();
        /**
         * Computes the Log map of the current primitive.
         * @param	none
         * @result	saved in the object on which the function was called.
         */
        virtual void Log();
        /**
         * Converts from symmetric representation back to
         * lie group representation.
         *
         * @param	prim	The result is saved in prim.
         * @return	true/false depending upon whether it was
         *			successful or not.
         */
        virtual bool convert2Lie( M3DPrimitive* prim ) const;
        virtual bool projectAtom2Vec(Vector & vecPrim);
        virtual bool vec2Atom(const Vector & vecPrim);
        /**
         * Computes the difference between this and the other primitive
         * and stores the result in the current primitive. This works
         * in manifold space and not in tangent space. No checks are performed
         * to ensure the correct state.
         */
        virtual SymPrimitive& operator-=(const SymPrimitive& b);

        /**
         * Composethe two primitives
         * and stores the result in the current primitive. This works
         * in manifold space and not in tangent space. No checks are performed
         * to ensure the correct state.
         */
        virtual SymPrimitive& operator+=(const SymPrimitive& b);

        Vector2D getTn0() {	return tn0;}
        Vector2D getTn1() { return tn1;}

        /** Primitive medial site position */
        Vector3D x;

        /**  The spoke radius */
        double r;

        /** The two spoke directions */
        Vector3D n0, n1;

        /** The two spoke directions in tangent space.
         * These members have to be present because
         * Vector3D and Vector2D cannot change length
         */
        Vector2D tn0, tn1;

        /** The elongation of this atom's B-vector */
        double elongation;

        const bool ignorePrimType;
        // Forcibly setting to false for now ...
        // no way to change it
        // Why do we need this kludge anyway?
        // - rrs
    };
};

#endif

