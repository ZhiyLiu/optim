#ifndef M3DTUBEPRIMITIVE_H
#define M3DTUBEPRIMITIVE_H
#include <assert.h>
#include <Vector2D.h>
#include <vector>

using std::vector;

#define EPSILON		1.0e-3

class M3DTubePrimitive : public virtual M3DPrimitive
{
public:
    // Constructor: Default makes "stock" primitive.
    M3DTubePrimitive() : baseAtom(false), d(0.0), sr(DEFAULT_NUMBER_OF_SPOKES, 0.0 ), U2(1,0,0), U0(0,1,0) {}

    // Constructor: Makes custom primitive.
    M3DTubePrimitive(double x0, double x1, double x2, double _d,
                     double* radii, const Vector3D& _U0, const Vector3D& _U2, 
                     int numberOfSpokes) :
        M3DPrimitive(x0,x1,x2), baseAtom(false),
        d(_d), sr(radii, radii + numberOfSpokes), U0(_U0), U2(_U2) {}

    // Constructor: Makes custom primitive.
    M3DTubePrimitive(const Vector3D& X, double _d,
                     double* radii, const Vector3D& _U0, 
                     const Vector3D& _U2, int numberOfSpokes) :
        M3DPrimitive(X), baseAtom(false),
        d(_d), sr(radii, radii + numberOfSpokes), U0(_U0), U2(_U2) {}

    // copy constructor.
    M3DTubePrimitive(const M3DTubePrimitive& prim) : 
        M3DPrimitive(prim), baseAtom(prim.baseAtom), sr(prim.sr),
        U2(prim.U2), U0(prim.U0), d(prim.d) {}

    // The assignment operator series to take care of
    // c++/compiler bugs.
    virtual M3DPrimitive & operator=(const M3DPrimitive & prim);
    M3DPrimitive& operator=( const M3DTubePrimitive& prim ) {
        return operator=((const M3DPrimitive&) prim);
    }

    static M3DTubePrimitive* readPrimitive(Registry& registry, int numberOfSpokes, const char * regStr, ...);
    virtual void writePrimitive(Registry& registry, const char * regStr, ...) const;

    virtual void writePrimitiveOldFormat(Registry& registry,
                                         const char * regStr, ...) const {};

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
    virtual M3DPrimitive * copyPtr() const {
        return new M3DTubePrimitive(*this);
    }

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
     * @param	none
     * @return	pointer to TubeSymPrimitive object
     *			*NOTE* the callee is responsible for clearing memory.
     */
    virtual SymPrimitive* convert2Sym() const {
        return new TubeSymPrimitive(this);
    }


    virtual bool vec2Atom(const Vector & vecPrim){
        TubeSymPrimitive symPrim ;
        symPrim.vec2Atom(vecPrim);
        return( symPrim.convert2Lie(this));
    }


    /**
     * Returns whether this is a base atom or not.
     * @return	true if a base atom, false otherwise
     */
    bool isBaseAtom() const {
        return baseAtom;
    }

    /**
     * Sets the base atom flag for this atom.
     * Note that only atom *should* be a base atom in a tube at any time.
     * @param	value	the new base atom flag value.
     * @return	none
     */
    void setBaseAtom( const bool value ) {
        baseAtom	= value;
    }

    /**
     * Returns the number of spokes per atom.
     * @return	number of spokes per atom.
     */
    int getNumberOfSpokes() const {
        return sr.size();
    }

    // These getters and setters are deprecated and should be phased out over time.
    virtual void setR(double r);
    virtual void setQ(const Quat & Q);
    virtual void setTheta(double _theta);

	// dibyendu
	// These setters are for the M3DQuadPrimitives only
	virtual void setU0( Vector3D _u0 ) { /* do nothing */ }
	virtual void setU1( Vector3D _u1 ) { /* do nothing */ }

    // Calculate a common radius for all spokes
    virtual double getR() const;
    // Return a quaternion describing the orientation
    virtual Quat getQ() const;
    // Return the angle between a representative spoke
    // and the bisector
    virtual double getTheta() const;

    virtual Vector3D getY0() const { return U0 * sr[0]; }
    virtual Vector3D getY1() const { return getYPhi(R_PI); }
    virtual Vector3D getU0() const { return U0; }
    virtual Vector3D getU1() const { return getNormalizedYPhi(R_PI); }
    virtual Vector3D getNormalizedY0() const { return U0; }
    virtual Vector3D getNormalizedY1() const { return getNormalizedYPhi(R_PI); }
    virtual Vector3D getB() const { return U2; }
    virtual Vector3D getN() const { return U0; }

    /**
     * Returns the radius of the nth spoke.
     * @param	n	The spoke number [0,NUM_SPOKES)
     * @return	The length of the corresponding spoke.
     */
    double getRN( const int n ) const {
        //assert( n >= 0 && n < sr.size() );

        return sr[n];
    }

    /**
     * Returns the angle of the nth spoke.
     * @param	phi	 The angle of the spoke [0, NUM_SPOKES)
     * @return	The angle (theta) of the corresponding spoke.
     */
    double getThetaN( const int n ) const {
        //assert( n >= 0 && n < sr.size() );

        return acos( d / sr[n] );
    }

    /**
     * Sets the radius of the nth spokes. In doing so it may change the angle of the nth spoke.
     * @param	n		The spoke number [0,NUM_SPOKES)
     * @param	radius	The new radius of the nth spoke
     * @return	true (false if unsuccessful, reserved for future)
     */
    bool setRN( const int n, const double radius ) {
        //assert( n >= 0 && n < sr.size() );
        sr[n]	= radius;
        return true;
    }

    /**
     * Sets the angle of the nth spoke. In doing so it may change the radius of the nth spoke.
     * @param	n		The spoke number [0,NUM_SPOKES)
     * @param	radius	The new radius of the nth spoke
     * @return	true if successful, false otherwise. The failure happens around angles of pi/2.
     */
    bool setThetaN( const int n, const double newTheta ) {
        //assert( n >= 0 && n < sr.size() );

        const double cnewTheta	= cos(newTheta);
        if( fabs(cnewTheta) < EPSILON ) {
            return false;
        }
        else {
            sr[n]	= d / cnewTheta;
            return true;
        }
    }


    /**
     * Returns the radius of the spoke at the angle phi.
     * @param	phi	The angle of the spoke [0, 2pi)
     *				The angle can be between (-2pi,+2pi) and is then remapped.
     * @return	The length of the corresponding spoke (linear interpolation)
     */
    double getRPhi( double phi ) const {
        if( phi < 0.0 ) {
            const int temp	= 1+int(-phi/(2*R_PI));
            phi += temp*2*R_PI;
        }
        const double n	= sr.size() * phi/(2.0*R_PI);
        const int n0	= int(n);
        const int n1	= (n0+1) % sr.size();
        const double w1	= n - n0;
        const double w0	= 1.0 - w1;
        const double srphi	= exp(log(sr[n0]) * w0 + log(sr[n1]) * w1);
        return srphi;
    }

    /**
     * Returns the angle of the spoke at the angle phi.
     * @param	phi	The angle of the spoke [0, 2pi)
     *				The angle can be between (-2pi,+2pi) and is then remapped.
     * @return	The angle (theta) of the corresponding spoke (linear interpolation)
     */
    double getThetaPhi( const double phi ) const {
        //return acos( getR() / getRPhi(phi) * cos(getTheta()) );
        return acos( d / getRPhi(phi) );
    }

    /**
     * Computes the nth unit spoke.
     * @param	n		The spoke number [0,NUM_SPOKES)
     * @return	The corresponding unit spoke.
     */
    Vector3D getNormalizedYN( const int n ) const {
        //assert( n >= 0 && n < sr.size() );
        const double theta	= getThetaN(n);
        const double phi	= 2 * R_PI * n / sr.size();
        // The spoke in the local frame.
    	Vector3D yn(cos(theta), - sin(theta) * cos(phi), sin(theta) * sin(phi) );
        // perform a rotation of yn by our orientation q
    	q().rotateVector(yn);
        return yn;
    }

    /**
     * Computes the nth spoke.
     * @param	n		The spoke number [0,NUM_SPOKES)
     * @return	The corresponding spoke.
     */
    Vector3D getYN( const int n ) const {
        //assert( n >= 0 && n < sr.size() );
        return getNormalizedYN(n) * getRN(n);
    }

    /**
     * Computes the unit spoke corresponding to the position phi [0,2pi] on the
     * circumference.
     * @param	phi	an angle from [0,2pi) that corresponds to v [0,1]
     *				on the circumference of the tube.
     *				The angle can be between (-2pi,+2pi) and is then remapped.
     * @return	The corresponding unit spoke.
     */
    Vector3D getNormalizedYPhi( const double phi ) const {
        const double theta	= getThetaPhi( phi );
        // The spoke in the local frame.
    	Vector3D yphi(cos(theta), - sin(theta) * cos(phi), sin(theta) * sin(phi) );
        // perform a rotation of yphi by our orientation q
    	q().rotateVector(yphi);
        return yphi;
    }

    /**
     * Computes the spoke corresponding to the position phi [0,2pi] on the
     * circumference.
     * @param   phi an angle from [0,2pi) that corresponds to v [0,1]
     *              on the circumference of the tube.
     *              The angle can be between (-2pi,+2pi) and is then remapped.
     * @return  The corresponding spoke.
     */
    Vector3D getYPhi( const double phi ) const {
        return getNormalizedYPhi(phi) * getRPhi(phi);
    }

    // Primitive transformation functions.
    // Scaling and rotation are about medial site location.
    virtual void rotateBy(Quat dQ);
    virtual void scaleBy(double mag);
	

	// dibyendu
	// overloaded scaleBy function for scaleing a particular spoke 	
	virtual void scaleSpokeBy(int spokeId, double mag) {
		// do nothing : function still to be written
	}

	// dibyendu
	// overloaded rotateSpokeBy function for rotating a particular spoke 	
	virtual void rotateSpokeBy(int spokeId, double delTheta, double delPhi) {
		// do nothing : function still to be written
	}


    /**
     * The number of spokes per tube atom.
     */
    // VC++ 6.0 does not let us define the default value here
    //static const int DEFAULT_NUMBER_OF_SPOKES	= 8;
    static const int DEFAULT_NUMBER_OF_SPOKES;

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
        return b * getR();
    }

private:
    friend class M3DSpokeProblem;
    friend class M3DSpokeOptimizer;
    /**
     * Provides a direct access to set spoke radius. This function is
     * only supposed to be used by M3DSpokeProblem.
     * @param	n	Spoke number.
     * @param	rn	Value of radius for the nth spoke.
     * @return	none
     */
    void setR( const unsigned int n, const double rn )
    {
        //assert( n >= 0 && n < sr.size() );
        // cos-1(1/1.02) = 11.36 degrees (approx 10)
        if( rn > 1.02 * d ) {
            sr[n]	= rn;
        }
        else {
            // If radius were to become this small,
            // then the spoke would come very close to the tangent and things
            // would become bad.
            // So just clamp ...
            sr[n]	= 1.02 * d;
        }
    }
    /**
     * Provides a direct access for spoke radius. This function is
     * only supposed to be used by M3DSpokeProblem.
     * @param	n	Spoke number.
     * @return	value of radius for nth spoke.
     */
    double getR( const unsigned int n ) const
    {
        //assert( n >= 0 && n < sr.size() );
        return sr[n];
    }

    double getD() const
    {
        return d;
    }

    void setD( const double d )
    {
        this->d	= d;
    }


protected:
    /**
     * Computes the quaternion corresponding to the orientation of this hub atom.
     * @return	A Quat object representing the orientation of this hub atom.
     */
    Quat q() const
    {
        return getQ();
    }

    /**
     * Whether this atom is a base atom or not. The base atom's rotational orientation
     * is used to enforce the orientation for the rest of the tube.
     */
    bool baseAtom;

    /**
     * List of radii of the spokes (how many number of spokes)
     */
    vector<double> sr;

    /**
     * Central (bisector) unit vector.
     */
    Vector3D U2;

    /**
     * The unit vector corresponding to the 'zero' spoke.
     */
    Vector3D U0;

    /**
     * distance from hub to plane 
     */
    double d;

    /**
     * Copy information relevant only to tube primitives.
     */
    void copy(const M3DTubePrimitive & prim);

    /**
     * This class is a placeholder for the symmetric space representation for tubes.
     * It is meant to be used internally only by tubes.
     */
    class TubeSymPrimitive : public SymPrimitive {
    public:
        /**
         * Default constructor.
         */
        TubeSymPrimitive() : x(0.0, 0.0, 0.0), d(1.0), sr(DEFAULT_NUMBER_OF_SPOKES, 0.0), U2(1.0,0.0,0.0), rEnd(0), ignorePrimType(false)
        {}

        /**
         * Constructor for creating a Symmetric primitive from a tube primitive.
         * It handles both ordinary tube primitives as well as end primitives.
         *
         * @param	prim	A tube primitive.
         */
        TubeSymPrimitive( const M3DTubePrimitive* prim );

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




        /** Primitive medial site position */
        Vector3D x;

        /**  The distance of the cone base plane from the hub */
        double d;

        /** Individual spoke radius multiplies */
        vector<double> sr;

        /** The tangent direction of the atom */
        Vector3D U2;

        /** The tangent direction in tangent space.
         * This member has to be present because
         * Vector3D and Vector2D cannot change length
         */
        Vector2D tU2;

        /** The length of this atom's B-vector */
        double rEnd;

        const bool ignorePrimType;
        // Forcibly setting to false for now ...
        // no way to change it
        // Why do we need this kludge anyway?
        // - rrs
    };

    friend class M3DPrimitive;
    friend class TubeSymPrimitive;
};


#undef EPSILON

#endif

