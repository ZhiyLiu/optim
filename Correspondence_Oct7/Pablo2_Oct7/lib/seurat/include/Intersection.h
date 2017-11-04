#ifndef INTERSECTION_H
#define INTERSECTION_H


/*
    Class Intersection.  This class constructs a bounding region about
	a set of hinge atoms.  Points inside the region must be tested, by
	calling testPoint(), to see of they are inside the blending region.
	Points outside can be assumed to not require blending.

	The blending region consists of a series of truncated cones with
	rounded ends, two per pair of link atoms.  The centers of the
	end caps are offset from the link atoms so that one truncated cone
	lies on each side of the interatom line at the figural intersection.

    The default radius used at the centers of the end caps is the radius
	of the corresponding atom.  The setMultiplier() function allows the
	caller to specify a factor to be used to expand, or potentially
	contract, this radius, and the radius of the bounding regions between
	the atoms.  The value of the multiplier thus affects the results of
	function testPoint().

  Limitation:  This class works poorly in the case attachment of the child
  figure to the end of the parent figure at an angle appoaching 90 degrees,
  so that the B vector of the child's link atoms is essentially perpendicular
  to the B vector of the parent's end atoms.

*/

#define Point	Vector3D

class M3DObject;

class Intersection 
{

public:

    Intersection(M3DFigureTreeNode & subNode, M3DObject * subObject,
		void * parentPList, void * childPList);
	void setMultiplier(double m) { factorSq = m*m; }

    Intersection(const Intersection & isc);
    Intersection & operator= (const Intersection & isc);
    virtual ~Intersection();

	bool testPoint(const Point & p);
	int parentId() { return mainFigureId; }
	int childId() { return subFigureId; }

	void render(double currentScale);	// For debugging
	void print();	// For debugging

protected:

	int nlinks;
	Point * vertex;		// End points of inter-atom line segments
	double * radius;	// Radii of atoms
	double factorSq;
	int mainFigureId;
	int subFigureId;

	bool test(const Point & p, int m, int n);
	bool testPoint(const Point & p, int n);
	void destroy();

};

#endif

