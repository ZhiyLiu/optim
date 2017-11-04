#ifndef M3DOBJECT_H
#define M3DOBJECT_H

// MSVC6: disable compiler warning about debugging symbol names too long
#pragma warning (disable:4786)


#include "M3DFigure.h"
#include "M3DFigureTreeNode.h"


/*  Class M3DObject.

	An object is a a collection of figures and a collection of figure trees,
	specifying connectivity between the figures.

    Assignment of an object is the same as copying it, except that in assignment,
    the object's name is allocated only once; copies of a object created by the
    assignment operator refer to the same string.  Copies created by the copy
    constructor contain a new name string.  This means that objects should not be
    passed to or returned from functions, unless it is intended that object name
	duplication occur.  Instead, a reference or pointer to the object should be
	passed.

    These same rules apply recursively to the figures contained in an object, but
    in addition to the figural names, their colors, tolerances, and other invariant
	properties (stored in class M3DFigure) are similarly handled.

    The newObject(), assign() and clone() functions are provided to simplify
    the creation of objects.  Thus, newObject() should be used to allocate a new
    object, assign() should be used to create a second instance of an existing
    object without duplication of the invariant fields, and clone() should be used
    to copy an existing object to produce a new object.

    For additional details about this, see M3DFigure.h and classes derived from it.

*/



class M3DPGAStats ;
class M3DPGAPrimitiveStats ;

// dibyendu - cpns stats

class M3DCPNSStats ;

class SimilarityTransform3D;
class WorldSystem;
class Image3D;


class M3DObject
{
	friend class M3DObjectFile;

public:
	// Type of PGA statistics to be applied
	enum modelType { NotAdaptive, Adaptive };

    M3DObject();
    M3DObject(const M3DObject & obj);
    M3DObject & operator = (const M3DObject & obj);

    ~M3DObject();

	M3DObject * newObject() const { return new M3DObject; }
	M3DObject * assign() const {
		M3DObject * object = new M3DObject;
		*object = *this;
		return object;
	}
	M3DObject * clone() const { return new M3DObject(*this); }

    void print(int markedFigIndex = -1, int markedPrimIndex = -1,
		bool dump = true, std::ostream & out = std::cout) const;
    void printTree(std::ostream & out = std::cout) const;

    const char * getName() const { return name; }
    char * copyName() const;
    void setName(const char * newName);

    int getFigureCount() const { return figures.size(); }
    int getPrimitiveCount() const;
    int getLandmarkCount() const;
    int getMarkedLandmark(int & figureId);

    M3DFigure * getFigurePtr(int figIndex) const;	// Locate by index
    M3DFigure * getFigurePtr(char * figName) const;	// Locate by name
    int getFigureIndex(M3DFigure * figurePtr) const;
    int getFigureAtomIndexes(int & primIndex) const;
	void setFigurePtr(int figIndex, M3DFigure* fig_in);

    M3DPrimitive * getPrimitivePtr(int primIndex) const;
    M3DPrimitive * getPrimitivePtr(int figIndex, int primIndex) const;

	// Traditionally, pablo's images used a left-handed world coordinate
	// system.  However, they were converted from the right-handed system
	// of PLUNC by flipping the image data along the Y-axis.  Models
	// built on the flipped images (most Comp. Sci. images) have a flipped
	// value of false, while those built on unflipped, PLUNC images, have a
	// value of true.  The print() function reports CompSci models as
	// left-handed, and models for PLUNC images, such as those used by
	// the ConStruct program, as right-handed.  See classes Image3D and
	// M3DObjectfile for more information.
	bool orientation() const { return flipped; }
	void orient(bool yFlip) { flipped = yFlip; }

    void addFigure(M3DFigure * figPtr);
    M3DFigure * removeFigure(int index);
    void deleteFigure(int index);
    bool verifyInBounds() const;
	void renumber(int * newFigureNums);

	bool testForPotentialLoop(int parentId, int childId);
	void addTree(M3DFigureTreeNode * tree, M3DObject * oldObject,
		int * correspondence);

    void replaceFigure(int index, M3DFigure * figurePtr);

    int getFigureTreeCount() const {
        return figureTrees.size();
    }
    M3DFigureTreeNode * getFigureTreeRoot(int treeId) const {
        if(figureTrees.size() > treeId)
            return figureTrees[treeId];

        return NULL;
    }
    M3DFigureTreeNode * getFigureTreeNode(int figureId);
    void addFigureTree(M3DFigureTreeNode * treeNode);
	bool attachFigureTreeNode(M3DFigureTreeNode * treeNode, int figureId);
	M3DFigureTreeNode * detachFigureTreeNode(int figureId);
	bool isRootNode(int figureId);

	// This class is merely a container for this SimilarityTransform3D object;
	// the class initializes it, but does not directly use it.
    SimilarityTransform3D * getTransformation() const { return transformation; }
    void setTransformation(SimilarityTransform3D * trans) {
		transformation = trans;
	}

	void setPGAStats(M3DPGAStats * pgaStats, bool keep = false);
	M3DPGAStats * getPGAStats() const { return pga_stats; }

	// dibyendu - cpns
	void setCPNSStats( M3DCPNSStats * cpnsStats, bool keep = false ) ;
	M3DCPNSStats * getCPNSStats() const { return cpns_stats ; }

	void setAtomPGAStats(M3DPGAPrimitiveStats * pgaStats, bool keep = false);
	M3DPGAPrimitiveStats * getAtomPGAStats() const { return atom_pga_stats; }
	//
	// Morphological operations:
	// The units for the dilationFactor are in model co-ordinates
	// i.e. 0-1 space
	//
	void dilate( const double dilationFactor );
	void erode( const double erosionFactor );

	// Subdivide the medial sheet(s) of this object
	void subdivide();

	// resample the medial sheet(s) of this object so that the atoms are regularly spaced.
	void resampleForRegularSpacing();

    // Selection functions
    void select();
    void deselect();
    void toggleSelected();
    bool isSelected();
    bool isAnySelected();

    bool isModified();
    void setModified(bool flag);

	// dibyendu
	bool isObjectDilated() { return isDilated ; } 
	double getDilationFactorInModelUnits() { return dilationFactorInModelUnits ; } 
	void setDilationInfo( double _dilationFactorInModelUnits ) ;

    Vector3D getCOG(bool all = false) const;

    void scaleWidth(double scaleFactor);
    void scaleBy(double scaleFactor);
    void scaleBy(double scaleFactor, const Vector3D & center);

    void translateBy(const Vector3D & vTrans, bool selectAll = false);

    void rotateBy(const Quat &q);
    void rotateBy(const Quat &q, const Vector3D & center);

    void applySimilarity(const SimilarityTransform3D & transform,
                         const Vector3D & center);
    void applySimilarityAboutCOG(const SimilarityTransform3D * transform);

    void deleteAll();

    void invertConstraints();

	void regularize(double stepsize, int iterations, bool verbose = false);

	// Boundary information
	////SubdivBoundary * const getSubdivBoundary() const { return boundary; }
	/*SubdivBoundary * const getSubdivBoundary(int figId) const 
	{ 
		if (boundaries.size() > figId && figId >= 0)
			return boundaries[figId];
		else
			return NULL;
	}
	void setDisplacements(Displacements * disp);*/

	M3DObject * loadedObject() const { return that; }
	void restore();		// Reverts this to that

	void setWorld(WorldSystem * world, bool original = false) {
		if (original)
			that->setWorld(world);
		wrld = world;
	}
	bool applyWorld(Image3D * image);
	bool unApplyWorld(Image3D * image);
	WorldSystem * getWorld() const { return wrld; }

	// Type of PGA statistics to be applied
	modelType isAdaptive()  { return type; }
	void setModelType(modelType t) { type = t; }

#ifdef BINARY

	/*  REGULARITY CALCULATION

		Calculates the distance between 2 mrep objects of the same atom
		count and type.
	*/
	double dist2FromObject(M3DObject * object, int figureID,
		DistanceType DT = AVE_EUCLIDEAN_DIST, bool verbose = false);

	/*  Create an M3DFigure such that each atom is at the average position
		predicted by its neighbors.

		Such a model is a complete step to the average of its neighbors,
		should such be desired.

		Also, it turns out that calculating the distance to such a figure
		is identical to calculating the geoedesic distance over this figure
		of each atom with the atom as predicted by its neighbors (below).

		For a quad figure, we need to know what neighborhood type we're 
		interested in, i.e., how to include the corner and edge atoms.

		Should this have a Geodesic vs. Euclidean average?
	*/
	M3DObject* fromAveOfNeighbors(PrimNeighborhoodDefn PND, int figureID);

	/*  This is a frequently used function to query optimizer performance
		whenever you are minimizing distance to neighbor prediction; it is
		nothing but a composition of fromAveOfNeighbors and dist2FromObject.
		Use DistType to get geodesic or Euclidean distances and averages.

		Note that THIS FUNCTION SCALES THE DISTANCE BY 100,000 when returning it!
	*/
	double dist2FromAveOfNeighbors(PrimNeighborhoodDefn PND, int figureID,
		DistanceType DT = AVE_EUCLIDEAN_DIST, bool verbose = false);

	/*  This takes a fractional step from here to the weighted average of some
		other model(s).

		Given an array of models <A B ... R> and an array of weights
		<x y ... z>, this returns the model where each atom i is
		(xAi + yBi + ... + zRi ) / (x + y + z) for all i

		DT is defined for GEODESIC_DIST or EUCLIDEAN_DIST only

		Useful for moving a model 1/2 way to the average of its neighbors,
		or 1/3 way to the average of neighbors and 1/3 of the way to a reference
		object, etc.
	*/
	M3DObject* fractionalStepToObject(M3DObject * targets[], int weights[],
		int numTargets, DistanceType DT, int figureID);

#endif

protected:


    char * name;
	int * nameCount;
	WorldSystem * wrld;
	bool flipped;

	bool isDilated ;
	double dilationFactorInModelUnits ;

	M3DPGAStats * pga_stats;
	M3DPGAPrimitiveStats * atom_pga_stats;

	// dibyendu - cpns
	M3DCPNSStats * cpns_stats ;

	M3DObject * that;	// A duplicate of this, set by class M3DObjectFile
	int * counter;

	modelType type;	// Adaptive or non-adaptive

    std::vector<M3DFigure *> figures;

    // Forest of trees of figures 
    std::vector<M3DFigureTreeNode *> figureTrees;

    // Similarity transformation -- used to initialize registration
    SimilarityTransform3D * transformation;

	// Boundaries for figures;
	////SubdivBoundary * boundary;
	//std::vector<SubdivBoundary *> boundaries;

	void copyTreeFigures(M3DFigureTreeNode * parent, M3DObject * oldObject,
		M3DFigureTreeNode * newTree, int * correspondence);
	void print_tree_nodes(std::ostream & out = std::cout) const;

    void addFigureWithoutTree(M3DFigure * figPtr) {
	    figures.insert(figures.end(), figPtr);
	}

	void renumberSubTree(M3DFigureTreeNode * node, int & newFigNum);

	void markHingeAtoms(M3DFigureTreeNode * node);
	int verifyConnectivity(M3DFigureTreeNode * node);
};


#endif

