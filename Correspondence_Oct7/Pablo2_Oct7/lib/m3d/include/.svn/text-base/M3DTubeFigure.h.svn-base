#ifndef M3DTUBEFIGURE_H
#define M3DTUBEFIGURE_H

#include <vector>
#include "M3DFigure.h"
#include "Geodesic.h"	// We could probably bury this all the way down in primitive
						// if we moved the atom->neighbor prediction code down there too.

/*  Class M3DTubeFigure.

    This is a derived class for figures represented as a curve of atoms.

    Assignment of a figure is the same as copying it, except that in assignment,
    the figure name, color, tolerance, and other invariant properties (stored in
    class M3DFigure) are allocated only once; copies of a figure created by the
    assignment operator thus refer to the same name string and other invariants.

    Copies created by the copy constructor contain a new name string and new
	allocations of the other invariant properties.

    The newFigure(), assign() and clone() functions are provided to simplify
    the creation of figures.  So, newFigure() should be used to allocate a
    new tube-figure, assign() should be used to create a second instance of an
    existing figure, and clone() should be used to copy an existing figure to
    produce a new figure.

    For additional details about this, see M3DFigure.h.

*/

class M3DTubeFigure : public M3DFigure
{
public:
    M3DTubeFigure();
    M3DTubeFigure(int _numColumns, int _numberOfSpokes	= M3DTubePrimitive::DEFAULT_NUMBER_OF_SPOKES);
    M3DTubeFigure(const M3DTubeFigure & fig);
    virtual ~M3DTubeFigure();
    virtual M3DFigure & operator = (const M3DFigure & figure);
	M3DFigure& operator =( const M3DTubeFigure& figure ) {
		return operator=((const M3DFigure&) figure);
	}

	M3DFigure * newFigure() const { return (M3DFigure *) new M3DTubeFigure; }
	M3DFigure * assign() const {
		M3DTubeFigure * figure = new M3DTubeFigure;
		*figure = *this;
		return (M3DFigure *) figure;
	}
	M3DFigure * clone() const { return (M3DFigure *) new M3DTubeFigure(*this); }

    virtual void print(bool dump = true, std::ostream & out = std::cout, int markedPrimIndex = -1) const;

    void makeRectangularStockFigure(int cols);
	void makeLandmarkStockFigure(int cols);
    void makeEllipticalStockFigure(int cols);

	/**
	 * This function fixes the directions of the tangents
	 * so that all the bisector vectors and tangent vectors
	 * point along increasing arc length.
	 *
	 * This function is expected to stick around only
	 * temporarily to convert from old tube formats to new
	 * ones.
	 *
	 * @return	true if any changes were made, else false.
	 */
	bool orientTubeTangents();

    int getPrimitiveCount() const;

	//dibyendu
	int getSpokeCount() const;

    int getColumnCount() const { return numColumns; }

    M3DPrimitive * getPrimitivePtr(int primIndex) const;
	// primitive index and column index are the same for a tube
	void figuralCoordinates(int index, int & colIndex) const { colIndex	= index; };

    void setPrimitivePtr(int primIndex, M3DPrimitive * newPrimPtr);

    void addColumn(M3DPrimitive* column, int columnIndex);

    void select();
    void deselect();
    void toggleSelected();
    bool isSelected();
    bool isAnySelected();
	int numberSelected();
	int columnSelected();

	void reverseRows();

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
	void getAtomsNeighbors(int index, int &numNeighbors, M3DPrimitive* neighbors[], PrimNeighborhoodDefn PND);

    virtual void writeFigure(const char * figureStr, Registry& registry);
	static M3DFigure* readFigure(int figureIndex, Registry& registry );
	virtual int NUM_PGA_PARAMS() const {
		return 8;
	}


	//
	// Call this function each time the tube needs a rotational re-alignment.
	// The re-alignment is performed by taking the base atom position
	// into account.
	//
	virtual void fixGlobalConsistency();

	/**
	 * Returns an interpolated atom for this tube at the medial coord
	 * given by *m.
	 */
	void atomAtCoordinates( const double* m, M3DPrimitive& prim ) const;

	/**
	 * Subdivide the current figure.
	 */
	virtual void subdivide();


	/**
	 * Resamples all the atoms so that they are evenly
	 * spaced.
	 */
	virtual void resampleForRegularSpacing();

	/**
	 * Returns the rSrad interpolated primitive at co-ordinates.
	 * 
	 * @param	u	co-ordinate(s) of interpolated primitive.
	 * 				Should be a double array of appropriate size.
	 * @param	prim	Where the interpolated primitive will be saved.
	 * @return	true if the operation succeeded, false otherwise.
	 */
	virtual bool getInterpolatedAtom( M3DPrimitive* prim, const double* u ) const;

	/**
	 * This function returns the curviness of the medial sheet.
	 *
	 * @param	none
	 * @return	the curviness of the medial sheet.
	 */
	virtual double curviness() const;

	/**
	 * Returns the index of the base atom in this tube.
	 * If no atom is marked as the base atom, then
	 * the center atom is chosen by default.
	 * @param	none
	 * @return index of base atom.
	 */
	int getBaseAtomIndex() const;

	/**
	 * Rotates the base atom by the specified angle
	 * in the plane defined by the plane defined by 
	 * N and Bperp. It then sympathetically rotates
	 * the rest of the tube primitives by the same
	 * amount, making sure that alignment is preserved.
	 *
	 * @param	phi	Amount to rotate the base atom by
	 * @return	none
	 */
	void rotateAlongAxisBy( double phi );

	/**
	 * Returns the number of spokes per primitive.
	 * @param	none
	 * @return	The number of spokes per primitive.
	 */
	virtual unsigned int getNumberOfSpokes() const {
		return numberOfSpokes;
	}

	static double curveLength(const Vector3D& p0, const Vector3D& p1,
		const Vector3D& d0, const Vector3D& d1, const double t);

public:
	static const char* friendlyFigureName;

protected:
    int numColumns;
	int numberOfSpokes;

    std::vector<M3DPrimitive *> primitives;
};

#endif

