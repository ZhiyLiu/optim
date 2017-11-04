#ifndef M3DFIGURE_H
#define M3DFIGURE_H

#include <string.h>

#include "M3DPrimitive.h"
#include "M3DEndPrimitive.h"
#include "M3DFigureStats.h"
#include "InterfiguralConstraints.h"
#include "SimilarityTransform3D.h"
#include "Geodesic.h"	//Could get buried in the primitives by modifying fig->fig gedesic dist
#include "GeodesicSym.h"	// We could probably bury this all the way down in primitive
						// if we moved the atom->neighbor prediction code down there too.
#include "SubdivBoundary.h"
//#include "../../paul_code/include/Registry.h"
#include "Registry.h"
//#include "M3DDisplayGlobals.h"
extern int globalVerbosity;
#include <vector>

extern const float DEFAULT_FIGURE_COLOR[3];

#define MAX_ATOMS_IN_FIGURE 400

// How to handle atom-neighborhoods for computing regularity in different fig types.
// DO NOT CHANGE THE ORDER OF THESE: they are tied to bpTune values
enum PrimNeighborhoodDefn {
	// QuadFigure Types
	ALL_FIRST_NEIGHBORS,
	PIN_CORNERS,
	PIN_CORNERS_AND_EDGES,
	EDGES_HAVE_2_NEIGHBORS,
	CORNERS_HAVE_3_NEIGHBORS,
};

// DO NOT CHANGE THE ORDER OF THESE: they are tied to bpTune values
enum DistanceType {
	EUCLIDEAN_DIST,
	AVE_EUCLIDEAN_DIST,
	GEODESIC_DIST,
	AVE_GEODESIC_DIST
};

/*  Class M3DFigure.
    This is a base class for figures, from which right now M3DQuadFigure and
	M3DTubeFigure are derived.

    For every instance of the same figure, only a copy of some features of the
    figure are specified.  These are the figure's name, visibility, polarity,
    color, type and tolerance.

    Assignment of a figure is the same as copying it, except that in assignment,
    the figure name and other invariant properties are not reallocated; copies
    of a figure created by the assignment operator refer to the same invariants.
    Copies created by the copy constructor contain a new name string and other
    invariants.  This means that figures should not be passed to or returned
    from functions, unless it is intended that duplication of the name and other
	invariants occur.  Instead, a reference or pointer to the figure should be
    passed.

    Note that because this is an abstract class, it cannot be instantiated.
    The pure virtual functions newFigure(), assign(), and clone() provide a
    way to produce a second object having the same derived type.  Member
    newFigure() can be used as a default constructor, producing an object of
    the same derived type as *this.  Members assign() and clone() do the same
    thing, but also copy the contents of *this into the new object.  The
    difference between them is that assign() acts as if the operator=()
    function were used for the copying, while clone() acts as if the copy
    constructor were used.  For example, when assign() is used, the resultant
    object shares the same name string as *this.  The string is actually
    duplicated if clone() is used.

    The point of this complexity is revealed by an example; when optimizing
    a figure, there is no need to copy the name, visibility, etc., since each
    refinement of the figure represents the a modification of the same figure.
    When a figure is copied and pasted, however, the result is two figures,
    and their invariants may later be changed, so these properties must be
    represented by different objects.

*/

#if defined(LM_METHOD_OBJ) && defined(SWIG)
// These functions have been declared but not defined, so make sure swig does not try to wrap them.
%ignore M3DFigure::clearLandmark(int);
%ignore M3DFigure::setLandmarkName(int, char const*);
%ignore M3DFigure::markLandmark(int);
#endif

class M3DFigure
{
public:
    M3DFigure();
    M3DFigure(const M3DFigure & fig);
    virtual ~M3DFigure();

    virtual M3DFigure & operator =(const M3DFigure & fig) {
		return copy(fig);
	}

	virtual M3DFigure * newFigure() const = 0; // Default constructor for derived classes
	virtual M3DFigure * assign() const = 0;	// Partial copy constructor for derived classes
	virtual M3DFigure * clone() const = 0;	// Copy constructor for derived classes

    virtual void print(bool dump = true, std::ostream & out = std::cout,
		int markedPrimIndex = -1) const;

    virtual int getPrimitiveCount() const = 0 ;

	// dibyendu
	virtual int getSpokeCount() const = 0 ;

    virtual bool verifyInBounds() const = 0;

    virtual M3DPrimitive * getPrimitivePtr(int primIndex) const = 0;
    virtual void setPrimitivePtr(int primIndex, M3DPrimitive * newPrimPtr) = 0;

	SubdivBoundary * getBoundaryPtr() const { return boundary; }
	void setBoundary(SubdivBoundary * refBoundary) {
		if (refBoundary == NULL)
		return;

		if (boundary != NULL)
			delete boundary;

		boundary = new SubdivBoundary(refBoundary);
	}
	// Image stats
	void setFigureStatsPtr(M3DFigureStats * newFigStatPtr);
	M3DFigureStats * getFigureStatsPtr();

	// Stacked image functions
	void setStackedImageMask(unsigned short mask) { stackedImageBit = mask; }
	void setStackedImageNumber(int num);	// num should be in the range [0, 15]
	unsigned short getStackedImageMask() { return stackedImageBit; }

    const char * getName() const { return properties->name; }
    char * copyName() const;
    void setName(const char * newName);

    const float * getColor() { return &(properties->color[0]); }
    void setColor(const float c[3]) {
        properties->color[0] = c[0];
        properties->color[1] = c[1];
        properties->color[2] = c[2];
    }

#ifdef LM_METHOD_OBJ
	// Landmark functions - new method
	void addLandmark(int atomIndex, const char *name, double atomT = 0);	// method objLM
	void addLandmark(double atomU, double atomV, const char *name, double atomT = 0);	// method objLM

	// For making changes after addition of a landmark
	void setLandmarkAtomU(int lmIndex, double u) { landmarkAtomUs[lmIndex] = u; }
	void setLandmarkAtomV(int lmIndex, double v) { landmarkAtomVs[lmIndex] = v; }
	void setLandmarkAtomT(int lmIndex, double t) { landmarkAtomTs[lmIndex] = t; }

	int getLandmarkCount() const { return landmarkCount; }
	int    getLandmarkAtomIndex(int lmIndex) { return landmarkAtomIndices[lmIndex];}
	double getLandmarkAtomT (int lmIndex) { return landmarkAtomTs[lmIndex];}

// For Quads, atomU is the column. Tubes: atomU is the index. In either case it is really an int
	double getLandmarkAtomU (int lmIndex) { return landmarkAtomUs[lmIndex];}

// For Quads, atomV is the row (really an int).  Tubes: atomV is [0..1] how far around the circumference
	double getLandmarkAtomV (int lmIndex) { return landmarkAtomVs[lmIndex];}

	char * getLandmarkName (int lmIndex) { return landmarkNames[lmIndex];}
//	int	countLandmarkIndicesByName(const char *name);
//	int	findLandmarkAtomIndexByName(const char *name, int instance = -1);
	int	findLandmarkAtomIndexByName(const char * name);
	const char *findLandmarkAtomNameByIndex(const int atomIndex);
	int	findLandmarkIndexByName(const char * name);  // lmIndex is index into public landmark* arrays

	// Landmark functions - original method
	// These are needed to compile Binary Pablo, but do not appear to ever be used.
	Vector3D getLandmark(int index);
	void addLandmark(Vector3D & landmark);
 	void clearLandmarks();
	bool clearLandmark(int index);
    bool setLandmarkName(int index, const char * newName);
    bool markLandmark(int index);
    void unmarkLandmark() { markedLandmark = -1; }
#endif

#ifndef LM_METHOD_OBJ
	// Landmark functions - original method
	int getLandmarkCount() const { return landmarks.size()/3; }
	void addLandmark(Vector3D & landmark);
	Vector3D getLandmark(int index);
	void clearLandmarks();
	bool clearLandmark(int index);
	const char * getLandmarkName(int index) const;
    bool setLandmarkName(int index, const char * newName);
    bool markLandmark(int index);
    void unmarkLandmark() { markedLandmark = -1; }
#endif
	
    int getMarkedLandmark() { return markedLandmark; }
	void translateLandmarks(const Vector3D & vTrans);
	void rotateLandmarks(const Quat &q, const Vector3D &vCenter);
	void scaleLandmarks(double scalefact, const Vector3D &vCenter);

    // Selection functions
    virtual void select() = 0;
    virtual void deselect() = 0;
    virtual void toggleSelected() = 0;
    virtual bool isSelected() = 0;
    virtual bool isAnySelected() = 0;
	virtual int numberSelected() = 0;

    bool isModified() { return properties->modified; }
    void setModified(bool flag) { properties->modified = flag; }

    bool isPositiveSpace() { return properties->positiveSpace; }
    void setPositiveSpace(bool flag) { properties->positiveSpace = flag; }

    void setTolerance(int tol) { properties->tolerance = tol; }

	static void setDefaultSurfaceTolerance(int tol) { defaultFigureTolerance = tol; }
	static int getDefaultSurfaceTolerance() { return defaultFigureTolerance; }

    bool hasPositivePolarity() { return properties->positivePolarity; }
    void setPositivePolarity(bool flag) { properties->positivePolarity = flag; }

    virtual Vector3D getCOG(bool all = false) const = 0;

    virtual void scaleWidth(double scalefact) = 0;
    virtual void scaleBy(double scalefact) = 0;
    virtual void scaleBy(double scalefact, const Vector3D &vCenter) = 0;

    virtual void rotateBy(const Quat &q) = 0;
    virtual void rotateBy(const Quat &q, const Vector3D &vCenter) = 0;

    virtual void translateBy(const Vector3D &vTrans, bool selectAll = false) = 0;

	void applySimilarity(const SimilarityTransform3D & transform,
		const Vector3D & center);

	void applySimilarityAboutCOG(const SimilarityTransform3D * transform) {
		applySimilarity(*transform, getCOG());
	}

    void setVisibility(bool toggle) { properties->visibility = toggle; }
    bool getVisibility() const { return properties->visibility; } 

    int getTolerance() const { return properties->tolerance; }

    virtual void deleteAll() = 0;

	// Constraints imposed by this figure on other figures
	InterfiguralConstraints & constraints() { return ifconstraints; }

	// Constraints imposed by all governors of this figure
	InterfiguralConstraints & inverseConstraints() { return inverse_constraints; }

	// REGULARITY CALCULATION

	double dist2FromFigure(M3DFigure * figure, DistanceType DT = EUCLIDEAN_DIST,
		bool verbose = false);
	// Calculates the Euclidian or Geodesic distance^2 between this figure and 
	// the arg figure as the average of the atom-by-atom squared distance sum

	virtual void getAtomsNeighbors(int index, int &numNeighbors, M3DPrimitive* neighbors[], PrimNeighborhoodDefn PND) = 0;

	M3DFigure* fromAveOfNeighbors(PrimNeighborhoodDefn PND = PIN_CORNERS_AND_EDGES);
	// Return a figure such that each atom is the average of its neighbors
	//
	// Requires in QuadFigure, or whatever derrived figure type a method to getAtomsNeighbors()

	M3DFigure* fractionalStepToFigure(M3DFigure* targets[], int weights[], int numTargets, DistanceType DT);
	//
	// Persistent object support:
	// derived classes should have a static const char* public member called friendlyFigureName
	// which contains err a friendly name of the type of figure
	//
	// For reading and writing figures:
    virtual void writeFigure(const char * figureStr, Registry& registry)	= 0;
	// readFigure needs to be static and cannot be virtual,
	// derived classes should implement a method to read the object
    static M3DFigure* readFigure(int figureIndex, Registry& registry);


	//
	// PGA helpers (FIXME: This should be in primitive class)
	//
	virtual int NUM_PGA_PARAMS() const	= 0;

	//
	// This function should be called whenever a figure is converted back from symmetric space
	// representation.
	//
	virtual void fixGlobalConsistency() {};


	/**
	 * This function returns the atom corresponding to the medial coord(s) m.
	 * The atom is interpolated from the neighbors. Note this function is
	 * different from just interpolating using M3DPrimitive's atomInterp.
	 * This function interpolates position along some curve, as opposed to 
	 * a simple linear interpolation. We happen to use Hermite splines for
	 * this purpose.
	 * 
	 * @param	m	medial co-ord(s)
	 * @param	prim	Saves the result in the reference variable prim.
	 */
	virtual void atomAtCoordinates( const double* m, M3DPrimitive& prim ) const = 0;

	/**
	 * Resamples the entire model.
	 */
	virtual void resampleForRegularSpacing();

	/**
	 * Subdivide the current figure.
	 */
	virtual void subdivide() = 0;

	/**
	 * This function returns the curviness of the medial sheet.
	 * @param	none
	 * @return	the curviness of the medial sheet.
	 */
	virtual double curviness() const = 0;

	/**
	 * Returns the number of spokes per primitive.
	 * @param	none
	 * @return	The number of spokes per primitive.
	 */
	virtual unsigned int getNumberOfSpokes() const = 0;

protected:

	// Invariant properties of the figure
	struct FigProps_t {

		char * name;

		float color[3];
		bool visibility;

		bool modified;

		bool positiveSpace;
		bool positivePolarity;

		int tolerance;	// Tolerance in angle between spoke vector and boundary surface

		int counter;
	} * properties;

#ifdef LM_METHOD_OBJ

public:
	// Landmarks (method objLM).
	// LATER: combine {atomIndex, name and atomT} into a struct
	int landmarkCount;
	int landmarkAtomIndices[MAX_ATOMS_IN_FIGURE];		// Atoms of landmarks
	double landmarkAtomTs[MAX_ATOMS_IN_FIGURE];			// t coordinate within atom;
	double landmarkAtomUs[MAX_ATOMS_IN_FIGURE];			// u coordinate within atom;
	double landmarkAtomVs[MAX_ATOMS_IN_FIGURE];			// v coordinate within atom;
										// as in (u,v,t) coordinate system
#endif

protected:

	// Landmarks (method 1)
	std::vector<double> landmarks;		// Coordinates of landmarks
#ifdef LM_METHOD_OBJ
	char *landmarkNames[MAX_ATOMS_IN_FIGURE];			// Landmark names
#else
	std::vector<char *> landmarkNames;	// Names of landmarks
#endif
	int markedLandmark;					// Number of designated landmark (for rendering)

	InterfiguralConstraints ifconstraints;			// Governed figures
	InterfiguralConstraints inverse_constraints;	// Governing figures

	static int defaultFigureTolerance;

	M3DFigureStats * figureStats;				//Statistics information for image match.
	unsigned short stackedImageBit;

	SubdivBoundary* boundary;

    // Copy function used in child classes
    M3DFigure & copy(const M3DFigure & prim);
};

inline std::ostream& operator<<(std::ostream& os, const M3DFigure& f)
{ f.print(true, os); return os; }

// Now include all the other types of figures
#include "M3DQuadFigure.h"
#include "M3DTubeFigure.h"

#endif

