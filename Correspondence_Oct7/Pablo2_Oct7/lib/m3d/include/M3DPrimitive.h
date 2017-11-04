#ifndef M3DPRIMITIVE_H
#define M3DPRIMITIVE_H

#include <math.h>
#include "Quat.h"
#include "Registry.h"
#include "VectorND.h"
#include "matrix.h"
#include "SimilarityTransform3D.h"

// dibyendu
#include "M3DSpoke.h"

#include <vector>

// Disable warning C4250 on Windows (inheritance via dominance)
#ifdef WIN32
#pragma warning ( disable : 4250 )
#endif

extern const int INVALID_PRIMITIVE_ID;

enum M3DPrimitiveType
{
    M3D_STANDARD_PRIMITIVE,
    M3D_END_PRIMITIVE
};

const char M3D_STANDARD_PRIMITIVE_STR[] = "StandardPrimitive";

class M3DPrimitive;
class Image3D;


/**
 * class SymPrimitive
 * @author Rohit Saboo
 *
 * This class is an abstract base class for all representing the
 * different primitive types in their symmetric space representation.
 * As of now there are two primitive types: quad and tube.
 * @see M3DQuadPrimitive::SymPrimitive
 * @see M3DTubePrimitive::SymPrimitive
 *
 */
class SymPrimitive {
public:
    /**
     * Computes the Exponential map of the current primitive.
     * @param	none
     * @result	saved in the object on which the function was called.
     */
    virtual void Exp()	= 0;
    /**
     * Computes the Log map of the current primitive.
     * @param	none
     * @result	saved in the object on which the function was called.
     */
    virtual void Log()	= 0;
    /**
     * Converts from symmetric representation back to
     * lie group representation.
     *
     * @param	prim	The result is saved in prim.
     * @return	true/false depending upon whether it was
     *			successful or not.
     */
    virtual bool convert2Lie( M3DPrimitive* prim ) const	= 0;
    /**
     * Computes the difference between two primitives
     */
    virtual SymPrimitive& operator-=(const SymPrimitive& b) = 0;
    virtual SymPrimitive& operator+=(const SymPrimitive& b) = 0;
     
    virtual bool projectAtom2Vec(Vector & vecPrim) = 0;
    virtual bool vec2Atom( const Vector & vecPrim) = 0;
};


/**
 * class M3DPrimitive
 * @author Rohit Saboo
 * 
 * This class is an abstract base class for all the different primitive types.
 * As of now, the two primitive types are quad and tube
 * @see M3DQuadPrimitive
 * @see M3DTubePrimitive
 *
 * Also, the M3DEndPrimitive class (@see M3DEndPrimitive) derives from this
 * class However, it too is an abstract class.  The end primitives for
 * tubes and quads derive from both M3DEndPrimitive and M3DPrimitive using
 * virtual inheritance, so that two base copies of M3DPrimitive do not get
 * inherited.
 */
class M3DPrimitive
{
public:

    M3DPrimitive();

    M3DPrimitive(double x0, double x1, double x2);

    M3DPrimitive(const Vector3D& X);

    M3DPrimitive(const M3DPrimitive & prim);

    virtual ~M3DPrimitive() {}

    // Assignment operator
    virtual M3DPrimitive & operator = (const M3DPrimitive & prim);

    /**
     * This function returns whether this object is
     * an end primitive or a standard primitive
     * It is advised to use this function instead of
     * dynamic_cast for elegance and better performance
     *
     * @return M3D_STANDARD_PRIMITIVE
     */

    virtual M3DPrimitiveType type() const { return M3D_STANDARD_PRIMITIVE; }

    virtual const char* name() {
        return M3D_STANDARD_PRIMITIVE_STR;
    }
    /**
     * A virtual copy constructor. Since C++ does not support
     * a virtual copy constructor, all derived (non-abstract) classes
     * need to implement this method. Replace M3DPrimitive by the
     * correct class name.
     */
    virtual M3DPrimitive * copyPtr() const = 0; /* {
        return new M3D<correct class name>Primitive(*this);
    } */

#ifndef PRODUCTION_VERSION
    static void setWorld(Image3D * image);
#endif
    virtual void print(std::ostream & out = std::cout, char * prefix = NULL,
                       bool marked = false) const;

    // Primitive set and get functions
    void setX(const Vector3D & X)               {x = X;}
    void setX(double X[3])                      {x.set(X);}
    void setX(double x0, double x1, double x2)  {x.set(x0, x1, x2);}

    Vector3D getX() const {return x;}


    // The following functions are conceptually specific to
    // M3DQuadPrimitives and should only appear at that level, but
    // they are called by outside code that was written before tubes.
    // Quads also have getR{0,1,End} methods, but they are not
    // included at this level because they are not required for
    // compilation.  The setters corresponding to these getters are
    // also not required at this level for compilation.
    virtual Vector3D getU0() const = 0;
    virtual Vector3D getU1() const = 0;
    // Synonyms for the above two functions; remove?
    virtual Vector3D getNormalizedY0() const = 0;
    virtual Vector3D getNormalizedY1() const = 0;

    virtual Vector3D getY0() const = 0;
    virtual Vector3D getY1() const = 0;
    // Calculate a common radius for all spokes
    virtual double getR() const = 0;

    // The following methods are remnants of the old Lie Group 
    // representation.  They should ideally be phased out.

    // After calling this->setR(someR), a call to getR should yield someR,
    // but the proportions of the lengths of the spokes should stay the same
    virtual void setR(double r) = 0;

    virtual void setQ(const Quat & Q) = 0;
    virtual void setTheta(double _theta) = 0;

	// Dibyendu
	// Set functions: to be primarily used for s-reps because there is no concept of bisector, angle etc. for s-reps

	virtual void setU0( Vector3D _U0 ) = 0;
	virtual void setU1( Vector3D _U1 ) = 0;

    // Return a quaternion describing the orientation
    virtual Quat getQ() const = 0;

    // Return the angle between a representative spoke
    // and the bisector
    virtual double getTheta() const = 0;

	// Return the angles between the two spokes and the B vector

    // Return B, the normalized bisector of U0 and U1.  B is 
    // calculated directly from U0 and U1.  If this object 
    // is an m-rep, then B should equal UEnd.
    virtual Vector3D getB() const = 0;

    // This is retained for compatibility with old code.
    // Returns the "extended B vector."  Note: Formerly, the
    // extended B vector defined the crest.  Now, UEnd serves 
    // that function.  In general, YEnd may not bisect Y0 and Y1,
    // but B and ExtendedB always will.
    // FIXME: This function is wrong. It's sometimes in the
    // opposite direction of getB().
    virtual Vector3D getExtendedB() const = 0;

    // BPerp == B cross N.  In an m-rep slab, BPerp and B are both
    // tangent to the medial sheet and normal to each other.  
    Vector3D getBPerp() const;

    // N is a unit vector that, in a slab, is normal to the 
    // medial sheet and on the same side as U1.
    virtual Vector3D getN() const = 0;

    // Selection functions
    void select()           {selected |= 0x1;}
    void deselect()         {selected &= ~(0x1);}
    void toggleSelected()   {selected ^= 0x1;}
    bool isSelected() const {return (selected & 0x1) ? true : false;}
    short getSelectionFlags() const { return selected; }

    // optimizer can treat 0x2 atoms differently
    void selectForRegularity()           {selected &= ~(0x2);}
    void deselectForRegularity()         {selected |= 0x2;}
    void toggleSelectedForRegularity()   {selected ^= 0x2;}
    bool isSelectedForRegularity() const 
    {return (selected & 0x2) ? false : true;}

    // Is the atom used as a hinge atom in a link?
    void toggleHinge(bool yesNo) { hinge = yesNo; }
    bool isHinge() const { return hinge; }

    // Primitive transformation functions.
    // Scaling and rotation are about medial site location.
    void translateBy(const Vector3D & dv) {x += dv;}
    virtual void rotateBy(Quat dQ) = 0;
    virtual void scaleBy(double mag) = 0;

	// Dibyendu
	// Function to scale a particular spoke of an s-rep
	virtual void scaleSpokeBy( int spokeId, double mag ) = 0 ;

	// Dibyendu
	// Function to rotate a particular spoke by an amount (delTheta, delPhi)
	virtual void rotateSpokeBy( int spokeId, double delTheta, double delPhi ) = 0 ;

    // Compiler trouble: This function was causing
    // problems with gcc version 3.2 and lower.
    // Ideally it should be a virtual function,
    // but given that the world is less than ideal,
    // we need to do this.  (We need to do what???)
    // 
    virtual void writePrimitive(Registry& registry, 
                                const char * regStr, ...) const = 0;
    void writePrimitiveOld(Registry& registry, const char * regStr, ...) const;

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
                             const M3DPrimitive* m2)	= 0;

    /**
     * Returns the distance^2 between two atoms. Derived classes
     * need to provide an implementation for this method.
     * @param	m	The second atom 
     * @param	radius	Some radius normalization(?) used
     *					for geodesic distance computation.
     * @return	double, geodesic distance^2.
     */
    virtual double atomSquareDistance( const M3DPrimitive* m, const double* radius = NULL ) const	= 0;

    /**
     * Returns the distance^2 between two atoms.
     * @param	m	The second atom 
     * @param	radius	Some radius normalization(?) used
     *					for geodesic distance computation.
     * @return	double, geodesic distance^2.
     */
    double atomDistance( const M3DPrimitive* m, const double* radius = NULL ) const {
        const double sqrDist	= atomSquareDistance( m, radius );
        return sqrt(sqrDist);
    }

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
                              const double* weights = NULL ) = 0;

    /**
     * Converts the current primitive to a primitive in
     * symmetric space representation.
     * Note, this operation carries over only the basic
     * elements of geometry and not stuff such as hinge et. al
     * @param	none
     * @return	A pointer to a SymPrimitive object.
     *			*NOTE* the callee is responsible for deallocating it.
     */
    virtual SymPrimitive* convert2Sym() const = 0;
    // align two sets of atoms by a similarty (or a subset of it) transformation
    static bool atomsAlignment(
        int nPrims, M3DPrimitive** mSet1, M3DPrimitive** mSet2,
        Vector3D &tran, Quat &q, double &s, bool alignAll= false);

    static bool atomsAlignmentByNormal(int nPrims, M3DPrimitive** mSet1,
                                       M3DPrimitive** mSet2, Vector3D &tran,
                                       Quat &q, double &s, bool alignAll);

    static SimilarityTransform3D atomsAlignmentSimtrans(
        int nPrims, M3DPrimitive** mSet1, M3DPrimitive** mSet2);

    static bool projectAtom2Vec(M3DPrimitive * prim, Vector & vecPrim){
        SymPrimitive * symPrim = prim->convert2Sym();
        return(symPrim->projectAtom2Vec(vecPrim));
    }
    // to project an atom into vector in its atom PGA space

    virtual bool vec2Atom(const Vector & vecPrim) = 0;
    // to turn a vector in its atom PGA space into an atom

    static bool composePrimitives( const M3DPrimitive * m1,  M3DPrimitive * m2) ;
    //m1 + m2 (preserving the order) and saves the result in m2.

    static bool subtractPrimitives(M3DPrimitive * m1, const M3DPrimitive * m2) ;
    // m1 - m2 (preserving the order) and saves the result in m1.

protected:
    // to calculate the sum of square distances between two sets of atoms
    static double atomSetSquareDistance(int nPrims, M3DPrimitive **mSet1, 
                                        M3DPrimitive **mSet2, double *radius = NULL);
    // apply f as a similarity transformation to mSetBase and store results in mSet
    static bool applyVector(VectorND & f, int nPrims, M3DPrimitive ** mSetBase,
                            M3DPrimitive ** mSet);

    virtual void calc_elongation() {} // For end primitives this isn't a no-op.

    // Data members.

    // Same for both old and new representations.
    Vector3D x;     // Primitive medial site position
    bool hinge;     // Is atom a hinge atom in a subfigure?
    short selected; // 0x1:atom is selected for modification

    // 0x2: atom is NOT selected for regularity calc
    // (flags &&'ed together?)
    /*
    // Old "Lie group representation".  
    // WARNING: These are not valid and are not set to anything meaningful.  
    // Tube code won't work until it's rewritten to use the new stuff. 
    double r;
    mutable Quat q;
    mutable double theta;
    */

#ifndef PRODUCTION_VERSION
    // Used so print() can report world coordinates
    // Why is this not in the production version?
    static Image3D * worldImage;
#endif

    /**
     * Functions such as these should be implemented
     * here and by base classes, to do the copy specific
     * to the class, this eases writing copy code where
     * multiple inheritence is involved.
     */
    void copy(const M3DPrimitive & prim);
};

inline std::ostream& operator<<(std::ostream& os, const M3DPrimitive& p)
{ p.print(os); return os; }

#include "M3DTubePrimitive.h"
#include "M3DQuadPrimitive.h"
#include "M3DEndPrimitive.h"

#endif

