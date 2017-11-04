#ifndef M3DQUADFIGURE_H
#define M3DQUADFIGURE_H

#include <vector>
#include "M3DFigure.h"
#include "SubdivBoundary.h"
#include "Geodesic.h"	// We could probably bury this all the way down in primitive
						// if we moved the atom->neighbor prediction code down there too.


/*  Class M3DQuadFigure.

    This is a derived class for figures represented as a rectangular mesh
	of atoms.

    Assignment of a figure is the same as copying it, except that in assignment,
    the figure name, color, tolerance, and other invariant properties (stored in
    class M3DFigure) are allocated only once; copies of a figure created by the
    assignment operator thus refer to the same name string and other invariants.

    Copies created by the copy constructor contain a new name string and new
	allocations of the other invariant properties.

    The newFigure(), assign() and clone() functions are provided to simplify
    the creation of figures.  So, newFigure() should be used to allocate a
    new quad-figure, assign() should be used to create a second instance of an
    existing figure, and clone() should be used to copy an existing figure to
    produce a new figure.

    For additional details about this, see M3DFigure.h.

*/

using namespace std;


class M3DQuadFigure : public M3DFigure
{
public:
    M3DQuadFigure();
    M3DQuadFigure(int _numRows, int _numColumns);
    M3DQuadFigure(const M3DQuadFigure & fig);
	// dibyendu - construct a quad figure from CPNS data
	M3DQuadFigure( double * CPNSDataVector, int _numRows, int _numCols ) ;
    virtual ~M3DQuadFigure();
    M3DFigure & operator =(const M3DFigure & figure);
    M3DQuadFigure & operator =(const M3DQuadFigure & figure) {
		// This function is needed to avoid bitwise copying,
		// when assigning one M3DQuadFigure to another. 
		*this = (M3DFigure &) figure;
		return *this;
	}

	M3DFigure * newFigure() const { return (M3DFigure *) new M3DQuadFigure; }
	M3DFigure * assign() const {
		M3DFigure * figure = new M3DQuadFigure;
		*figure = *this;
		return figure;
	}
	M3DFigure * clone() const { return (M3DFigure *) new M3DQuadFigure(*this); }

    virtual void print(bool dump = true, std::ostream & out = std::cout,
		int markedPrimIndex = -1) const;

    void makeRectangularStockFigure(int rows, int cols);
	void makeLandmarkStockFigure(int rows, int cols);
    void makeEllipticalStockFigure(int rows, int cols);

    int getPrimitiveCount() const;
	// dibyendu - get the total no. of spokes in a figure
	int getSpokeCount() const;
    int getRowCount() const { return numRows; }
    int getColumnCount() const { return numColumns; }

	// dibyendu - vectorize this figure to a vector that is used by the CPNS program
	// [ positions ; radii ; spoke_directions ]
	double * vectorize() ;


    M3DPrimitive * getPrimitivePtr(int primIndex) const;
    M3DPrimitive * getPrimitivePtr(int rowIndex, int colIndex) const;

	int getPrimitiveID(int rowIndex, int colIndex) const;
	int indexOfAtom(int rowIndex, int colIndex) {
		return rowIndex*numColumns + colIndex;
	}
	void figuralCoordinates(int index, int & rowIndex, int & colIndex) const;

    void setPrimitivePtr(int primIndex, M3DPrimitive * newPrimPtr);
    void setPrimitivePtr(int rowIndex, int colIndex, M3DPrimitive * newPrimPtr);

    void addRow(std::vector<M3DPrimitive *> row, int rowIndex);
    void addColumn(std::vector<M3DPrimitive *> column, int columnIndex);

    void select();
    void deselect();
    void toggleSelected();
    bool isSelected();
    bool isAnySelected();
	int numberSelected();
	int rowSelected();
	int columnSelected();

	void reverseRows();
	void reverseColumns();
	void transpose(int ** map = NULL);

    bool verifyInBounds() const;

    Vector3D getCOG(bool all = false) const;

    void scaleWidth(double scalefact);
    void scaleBy(double scalefact);
    void scaleBy(double scalefact, const Vector3D &vCenter);

    void rotateBy(const Quat &q);
    void rotateBy(const Quat &q, const Vector3D &vCenter);

    void translateBy(const Vector3D &vTrans, bool selectAll = false);

    void deleteAll();


	// REGULARITY CALCULATIONS
	void getAtomsNeighbors(int index, int & numNeighbors, M3DPrimitive * neighbors[],
		PrimNeighborhoodDefn PND);
    virtual void writeFigure(const char * figureStr, Registry& registry);
	static M3DFigure* readFigure(int figureIndex, Registry& registry );

	//
	// PGA helpers (FIXME: This should be in primitive class)
	//
	virtual int NUM_PGA_PARAMS() const {
		return 9;
	}

	/**
	 * Subdivide the current figure. In quads it measures the distance between rows and that between
	 * columns, if avg difference is greater than 1.5 times in one direction, then only that
	 * direction is subdivded, else the entire figure is subdivided.
	 *
	 * In slab figures it is implemented using geodesic interpolation of primitives
	 */
	virtual void subdivide();

	/**
	 * This function returns the curviness of the medial sheet.
	 * Does nothing for a quad figure right now
	 * TODO: Implement this function for quads.
	 *
	 * @param	none
	 * @return	the curviness of the medial sheet.
	 */
	virtual double curviness() const;


	/**
	 * Returns an interpolated atom for a slab figure at the medial coord
	 * given by (u,v) = (m[0], m[1] ).
	 */
	void atomAtCoordinates( const double* m, M3DPrimitive& prim ) const;

	/**
	 * Returns the number of spokes per primitive.
	 * @param	none
	 * @return	The number of spokes per primitive.
	 */
	virtual unsigned int getNumberOfSpokes() const {
		return 2;
	}

public:
	static const char* friendlyFigureName;

protected:

    int numRows;
    int numColumns;

    std::vector<M3DPrimitive *> primitives;
};

#endif

