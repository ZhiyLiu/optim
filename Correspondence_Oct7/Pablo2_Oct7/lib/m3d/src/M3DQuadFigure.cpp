//#include <typeinfo.h>
#include <math.h>
#include <stdio.h>
#include "M3DQuadFigure.h"
#include <assert.h>
//#define DEBUG

#include <typeinfo>

const char* M3DQuadFigure::friendlyFigureName	= "QuadFigure";

using namespace std;

M3DQuadFigure::M3DQuadFigure() : M3DFigure()
{
#ifdef DEBUG
	cout << "M3DQuadFigure::M3DQuadFigure()" << endl;
#endif
    numRows = 0;
    numColumns = 0;

	boundary = new SubdivBoundary;
}

M3DQuadFigure::M3DQuadFigure(int _numRows, int _numColumns)
	: M3DFigure()
{
    int size;

#ifdef DEBUG
	cout << "M3DQuadFigure::M3DQuadFigure(int, int)" << endl;
#endif

    if(_numRows <= 0 || _numColumns <= 0)
    {
        numRows = 0;
        numColumns = 0;
        return;
    }
    else
    {
        numRows = _numRows;
        numColumns = _numColumns;
    }

    size = numRows * numColumns;
    primitives.resize(size);
    for(int i = 0; i < size; i++)
		primitives[i] = NULL;

	boundary = new SubdivBoundary;
}

M3DQuadFigure::M3DQuadFigure(const M3DQuadFigure & figure) : M3DFigure(figure)
{
    M3DPrimitive * primPtr,
                 * newPrim;
    int size;

#ifdef DEBUG
	cout << "M3DQuadFigure::M3DQuadFigure(const M3DQuadFigure &)" << endl;
#endif
    numRows = figure.numRows;
    numColumns = figure.numColumns;

    size = numRows * numColumns;
    for(int i = 0; i < size; i++)
    {
        primPtr = figure.primitives[i];
        if(primPtr != NULL)
            newPrim = primPtr->copyPtr();
        else
            newPrim = NULL;

        primitives.insert(primitives.end(), newPrim);
    }

	if (figure.getBoundaryPtr() != NULL)
		boundary = new SubdivBoundary(figure.getBoundaryPtr());
	else
		boundary = new SubdivBoundary;

}

// Dibyendu

/* 
Function to construct a figure from a vector of parameters as used in the CPNS program
The following is the format: [ positions ; radii ; spoke_directions ]
Standard primitive: (Px, Py, Pz, r1, r2, U1x, U1y, U1z, U2x, U2y, U2z)
End primitive: (Px, Py, Pz, r1, r2, r3, U1x, U1y, U1z, U2x, U2y, U2z, U3x, U3y, U3z)
*/ 

M3DQuadFigure::M3DQuadFigure(double * CPNSDataVector, int _numRows, int _numCols) {

	if( CPNSDataVector == NULL || _numRows <= 0 || _numCols <= 0 ) {
		cout << "Error in initializing M3DQuadFigure from CPNS data vector! Invalid input" << endl ;
		return ;
	}

	// set the number of rows and columns for this figure

	numRows = _numRows ;

	numColumns = _numCols ;

    int nTotalAtoms = numRows * numColumns;

	int nEndPrimitives = 2 * numRows + 2 * ( numColumns - 2 ) ;

	int nStdPrimitives = numRows * numColumns - nEndPrimitives ;

	// set the size of the primitives vector

	primitives.clear() ;

    primitives.resize(nTotalAtoms, NULL);
	

	// ------------------------  PARSE THE INPUT VECTOR  ---------------------------------

	// position vectors: n_total_positions 	= 3 * ( n_end_atoms + n_std_atoms )

	double* p = new double [ 3 * nTotalAtoms ] ;

	for( int ii = 0 ; ii < 3 * nTotalAtoms ; ii ++ )
		p[ii] = CPNSDataVector[ii] ;

	// n_total_radii		= 3 * n_end_atoms + 2 * n_std_atoms

	double* r = new double [ 3 * nEndPrimitives + 2 * nStdPrimitives  ] ;

	for( int ii = 0 ; ii < 3 * nEndPrimitives + 2 * nStdPrimitives ; ii ++ ) {		

		r[ii] = CPNSDataVector[ 3*nTotalAtoms + ii ] ;

	}	

	// n_total_unit_spokes	= 9 * n_end_atoms + 6 * n_std_atoms

	double* u = new double [ 9 * nEndPrimitives + 6 * nStdPrimitives ] ;

	for( int ii = 0 ; ii < 9 * nEndPrimitives + 6 * nStdPrimitives ; ii ++ ) {

		int uOffset = 3 * nTotalAtoms + 3 * nEndPrimitives + 2 * nStdPrimitives ;

		u[ii] = CPNSDataVector[ uOffset + ii ] ;
	}

	// Parse the radii and unit vectors as you traverse the column

	M3DPrimitive*			prim = NULL ;

	int pIndex = 0 ;
	int rIndex = 0 ;
	int uIndex = 0 ;	

	for( int i = 0 ; i < numRows ; i ++ ) {

		for( int j = 0 ; j < numColumns ; j ++ ) {

			if( i == 0 || i == numRows-1 || j == 0 || j == numColumns-1 ) {

				// end primitive

				prim = new M3DQuadEndPrimitive(	p[pIndex+0], p[pIndex+1], p[pIndex+2], 
												r[rIndex+0], r[rIndex+1], r[rIndex+2],
												u[uIndex+0], u[uIndex+1], u[uIndex+2], 
												u[uIndex+3], u[uIndex+4], u[uIndex+5],
                                                u[uIndex+6], u[uIndex+7], u[uIndex+8],0.0,0.0,0.0,0.0,0.0,0.0) ;//----Liyun keep old method,deltaU,deltaV set to 0

				prim->select() ;

				pIndex += 3 ;
				rIndex += 3 ;
				uIndex += 9 ;				

				setPrimitivePtr( i, j, prim ) ;				

			}
			else {

				// std primitive

				prim = new M3DQuadPrimitive(	p[pIndex+0], p[pIndex+1], p[pIndex+2], 
                                                r[rIndex+0], r[rIndex+1], r[rIndex+2],
												u[uIndex+0], u[uIndex+1], u[uIndex+2], 
												u[uIndex+3], u[uIndex+4], u[uIndex+5],
                                                u[uIndex+3], u[uIndex+4], u[uIndex+5],0.0,0.0,0.0,0.0,0.0,0.0) ; //----Liyun keep old method,deltaU,deltaV set to 0
				prim->select() ;

				pIndex += 3 ;
				rIndex += 2 ;
				uIndex += 6 ;

				setPrimitivePtr( i, j, prim ) ;
			}			
		}
	}

	boundary = new SubdivBoundary ;

}

M3DFigure & M3DQuadFigure::operator = (const M3DFigure & unknown_figure)
{
    M3DPrimitive * primPtr,
                 * newPrim;
    int size;

	assert( typeid(*this) == typeid(unknown_figure));

	const M3DQuadFigure& figure	= dynamic_cast<const M3DQuadFigure&>(unknown_figure);
#ifdef DEBUG
	cout << "M3DQuadFigure::operator=()" << endl;
#endif
	if (&figure == this)
		return *this;

    deleteAll();

	M3DFigure::copy(figure);

	M3DQuadFigure & fig = (M3DQuadFigure &) figure;
    numRows = fig.numRows;
    numColumns = fig.numColumns;

    size = numRows * numColumns;
    for(int i = 0; i < size; i++)
    {
        primPtr = fig.primitives[i];
        if(primPtr != NULL)
            newPrim = dynamic_cast<M3DQuadPrimitive*>(primPtr->copyPtr());
        else
            newPrim = NULL;

        primitives.insert(primitives.end(), newPrim);
    }

	if (fig.getBoundaryPtr() != NULL)
		this->boundary = new SubdivBoundary(fig.getBoundaryPtr());

    return (*this);
}

M3DQuadFigure::~M3DQuadFigure()
{
#ifdef DEBUG
	cout << "M3DQuadFigure::~M3DQuadFigure()" << endl;
#endif
    deleteAll();
}

void M3DQuadFigure::print(bool dump, ostream & out, int markedPrimIndex) const
{
    int i, j;
	int index;

    M3DFigure::print(dump, out, markedPrimIndex);

	out << "Shape: " << numRows << " rows by " << numColumns << " columns\n";
	index = 0;
    for(i = 0; i < numRows; i++)
    {
        for(j = 0; j < numColumns; j++)
        {
            out << "    Primitive(" << i << ", " << j << ')';
			out << " at 0x" << hex << primitives[i*numColumns + j] << dec;
			if (dump) {
			    out << ":\n";
				primitives[i*numColumns + j]->print(out, "\t",
					(index == markedPrimIndex) ? true : false);
			}
			else
				out << '\n';
			index++;
        }
    }

	int n = ifconstraints.size();
	if (n > 0) {
		out << "Constrained figures: ";
		for (i = 0; i < n; i++)
			out << ifconstraints.figure(i) << '(' << ifconstraints.distance(i) << ") ";
		out << '\n';
	}

	n = inverse_constraints.size();
	if (n > 0) {
		out << "Governing figures: ";
		for (i = 0; i < n; i++)
			out << inverse_constraints.figure(i) << '(' << inverse_constraints.distance(i) << ") ";
		out << '\n';
	}
	out << flush;
}

/*
void M3DQuadFigure::subdivideColumns()
{
	// Just subdivide in one of the directions.
	int nColsx2	= numColumns * 2 - 1;
	int i, r;
	std::vector<M3DPrimitive *> newPrimitives(nColsx2 * numRows);
	for(r = 0; r != numRows; ++r ) {
		for(i = 0; i < numColumns; ++i ) {
			newPrimitives[r*nColsx2 + 2*i]	= primitives[r*numColumns + i];
		}
		for(i = 0; i < numColumns - 1; ++i ) {
			const double u[2] = {r, i + 0.5};
			const int index	= r*nColsx2 + 2*i + 1;
			if( r == 0 || r == numRows - 1 ) {
				newPrimitives[index]	= new M3DQuadEndPrimitive(*dynamic_cast<M3DQuadEndPrimitive*>(newPrimitives[index-1]));
			}
			else {
				newPrimitives[index]	= new M3DQuadPrimitive(*dynamic_cast<M3DQuadPrimitive*>(newPrimitives[index-1]));
			}
	//		getInterpolatedAtom( newPrimitives[index], u );
			atomAtCoordinates(u,*(newPrimitives[index]));
		}
	}
	primitives	= newPrimitives;
	numColumns	= nColsx2;

	fixGlobalConsistency();
}
*/

void M3DQuadFigure::subdivide()
{
	// Just subdivide in one of the directions.
	int nRowsx2	= numRows * 2 - 1;
	int c, r;
	std::vector<M3DPrimitive *> newPrimitives(nRowsx2 * numColumns);
	for(c = 0; c != numColumns; ++c ) {
		for(r = 0; r < numRows; ++r ) {
			newPrimitives[2*r*numColumns + c]	= primitives[r*numColumns + c];
		}
		for(r = 0; r < numRows - 1; ++r ) {
			const double u[2] = {r + 0.5, c};
			const int index	= (2*r+1)*numColumns + c;
			if( c == 0 || c == numColumns - 1 ) {
				newPrimitives[index]	= new M3DQuadEndPrimitive(*dynamic_cast<M3DQuadEndPrimitive*>(newPrimitives[index-numColumns]));
			}
			else {
				newPrimitives[index]	= new M3DQuadPrimitive(*dynamic_cast<M3DQuadPrimitive*>(newPrimitives[index-numColumns]));
			}
	//		getInterpolatedAtom( newPrimitives[index], u );
			atomAtCoordinates(u,*(newPrimitives[index]));
		}
	}
	primitives	= newPrimitives;
	numRows		= nRowsx2;

	fixGlobalConsistency();
}


void M3DQuadFigure::atomAtCoordinates( const double* m, M3DPrimitive& prim ) const
{
	if( m[0] < 0.0 || m[0] > getRowCount() - 1 ||
		m[1] < 0.0 || m[1] > getColumnCount() - 1 ) {
		cout << "m = [" << m[0] << "," << m[1] << "], "
			 << getRowCount() << "," << getColumnCount() << endl;
		assert(false);
		return;
	}
	const int u		= int(m[0]);
	const int v		= int(m[1]);
	double t;
	const M3DQuadPrimitive *prim0, *prim1;
	if( m[0] - u > 1e-9 && m[1] - v > 1e-9 ) {
		cout << "M3DQuadFigure::atomAtCoordinates(...): Only one of u or v can be non-integer at present." << endl;
		assert(false);
	}
	if( m[0] - u > 1e-9 ) {
		t	= m[0] - u;
		prim0	= dynamic_cast<const M3DQuadPrimitive*>(getPrimitivePtr(u, v));
		prim1	= dynamic_cast<const M3DQuadPrimitive*>(getPrimitivePtr(u+1, v));
	}
	else {
		t	= m[1] - v;
		prim0	= dynamic_cast<const M3DQuadPrimitive*>(getPrimitivePtr(u, v));
		prim1	= dynamic_cast<const M3DQuadPrimitive*>(getPrimitivePtr(u, v+1));
	}

	if(!prim.atomInterp( t, prim0, prim1 )) {
		assert(false);
	}
}


void M3DQuadFigure::makeRectangularStockFigure(int rows, int cols)
{
    M3DPrimitive * primPtr;
    Quat q;
    int size;

    if(rows <= 0 || cols <= 0)
        return;

    size = rows * cols;

    q.setAxisAngle(Vector3D(1.0, 0.0, 0.0), R_HALF_PI);

    deleteAll();

    for(int i = 0; i < rows; i++)
    {
        for(int j = 0; j < cols; j++)
        {
            if(i == 0 || i == rows - 1 || j == 0 || j == cols - 1)
                primPtr = new M3DQuadEndPrimitive((double) j, (double) i, 0.0, 0.5, q, R_HALF_PI, 1.0);
            else
                primPtr = new M3DQuadPrimitive((double) j, (double) i, 0.0, 0.5, q, R_HALF_PI);
            primitives.insert(primitives.end(), primPtr);
        }
    }

    numRows = rows;
    numColumns = cols;

    select();

    if(numRows >= numColumns)
        scaleBy(1.0 / (double) (numRows + 1));
    else
        scaleBy(1.0 / (double) (numColumns + 1));

    q.setAxisAngle(Vector3D(0.0, 0.0, 1.0), R_HALF_PI);
    rotateBy(q);

    translateBy(Vector3D(0.5, 0.5, 0.5) - getCOG());

    properties->modified = true;
}

void M3DQuadFigure::makeLandmarkStockFigure(int rows, int cols)
{	
	makeEllipticalStockFigure(rows, cols);

	// set all atoms' positions to the position of the first atom so that all
	// will be either selected or deselected together. This makes them act
	// as a single point in space, suitable for landmarks.  -gst 20041103
	Vector3D posFirstPrim = getPrimitivePtr(0, 0)->getX();
	for(int i = 0; i < rows; i++)
	{
        for(int j = 0; j < cols; j++)
        {
            getPrimitivePtr(i, j)->setX(posFirstPrim);
        }
	}
}

#ifdef BINARY

void M3DQuadFigure::makeEllipticalStockFigure(int rows, int cols)
{
    M3DPrimitive * primPtr;
    Quat q,q2;
    int size,
        i,
        j;


    if(rows <= 0 || cols <= 0)
        return;

    size = rows * cols;

	q.setAxisAngle(Vector3D(1.0, 0.0, 0.0), R_HALF_PI);
	double angle	= 0.0;

    deleteAll();

	// FIXME: ROHIT: Changed the code to make a nice looking slab as the stock figure.
    for(i = 0; i < rows; i++)
    {
        for(j = 0; j < cols; j++)
        {
            if(i == 0 || i == rows - 1 || j == 0 || j == cols - 1) {
				if( j == 0 ) {
					angle	= (i==0) ? -1.5*R_HALF_PI : ( (i==rows-1) ? +1.5*R_HALF_PI : 2.0*R_HALF_PI);
				}
				else if( i == rows-1 ) {
					angle	= (j == cols-1) ? 0.5*R_HALF_PI : 1.0*R_HALF_PI;
				}
				else if( j == cols-1 ) {
					angle	= ( i == 0 ) ? -0.5*R_HALF_PI : 0.0;
				}
				else if( i == 0 ) {
					angle	= -R_HALF_PI;
				}
				q2.setAxisAngle(Vector3D(0.0,0.0,1.0), angle);
                primPtr = new M3DQuadEndPrimitive((double) j, (double) i, 0.0, 0.5, q2*q, R_HALF_PI * 2.0/3.0, 1.0);
			}
            else {
                primPtr = new M3DQuadPrimitive((double) j, (double) i, 0.0, 0.5, q, R_HALF_PI);
			}
            primitives.insert(primitives.end(), primPtr);
        }
    }

    numRows = rows;
    numColumns = cols;

    select();

    if(numRows >= numColumns)
        scaleBy(1.0 / (double) (numRows + 1));
    else
        scaleBy(1.0 / (double) (numColumns + 1));

    q.setAxisAngle(Vector3D(0.0, 0.0, 1.0), R_HALF_PI);
    rotateBy(q);
    translateBy(Vector3D(0.5, 0.5, 0.5) - getCOG());

    properties->modified = true;
#ifndef DEBUG
    cout << "Added stock figure with " << rows << " rows and " << cols << " columns" << endl;
#endif
}

#else

void M3DQuadFigure::makeEllipticalStockFigure(int rows, int cols)
{
    M3DPrimitive * primPtr;
    Quat q;
    int size,
        i,
        j;


    if(rows <= 0 || cols <= 0)
        return;

    size = rows * cols;

    q.setAxisAngle(Vector3D(1.0, 0.0, 0.0), R_HALF_PI);

    deleteAll();

    for(i = 0; i < rows; i++)
    {
        for(j = 0; j < cols; j++)
        {
            if(i == 0 || i == rows - 1 || j == 0 || j == cols - 1)
                primPtr = new M3DQuadEndPrimitive((double) j, (double) i, 0.0, 0.5, q, R_HALF_PI * 2.0/3.0, 1.0);
            else
                primPtr = new M3DQuadPrimitive((double) j, (double) i, 0.0, 0.5, q, R_HALF_PI);
            primitives.insert(primitives.end(), primPtr);
        }
    }

    numRows = rows;
    numColumns = cols;

    select();

    if(numRows >= numColumns)
        scaleBy(1.0 / (double) (numRows + 1));
    else
        scaleBy(1.0 / (double) (numColumns + 1));

    q.setAxisAngle(Vector3D(0.0, 0.0, 1.0), R_HALF_PI);
    rotateBy(q);

    translateBy(Vector3D(0.5, 0.5, 0.5) - getCOG());

    Vector3D cog = getCOG();

    for(i = 0; i < rows; i++)
    {
        for(j = 0; j < cols; j++)
        {
            primPtr = getPrimitivePtr(i, j);
            if(primPtr == NULL)
                continue;

            Vector3D pos = primPtr->getX();
            Vector3D diff = pos - cog;
            double length = diff.norm();
            if(length >= R_SMALL_TOLERANCE)
            {
                diff.normalize();
                Vector3D b = primPtr->getB();
                b.normalize();
                double cosAng = diff * b;
                Vector3D crossProd = b.cross(diff);

                if(cosAng < 0.0)
                    q.setAxisAngle(crossProd, acos(cosAng));
                else
                    q.setAxisAngle(crossProd, acos(cosAng));

                primPtr->rotateBy(q);
            }

#ifdef DEBUG
            printf("Primitive (%d, %d):\n", i, j);
            printf("    B: "); b.print();
            printf("    diff: "); diff.print();
            printf("    Dot: %f", cosAng);
            printf("    Angle: %f\n", acos(cosAng));
            printf("\n");
#endif
        }
    }

    properties->modified = true;
#ifndef DEBUG
    cout << "Added stock figure with " << rows << " rows and " << cols << " columns" << endl;
#endif
}

#endif	/* BINARY */

int M3DQuadFigure::getPrimitiveCount() const
{
    return primitives.size();
}

// dibyendu

int M3DQuadFigure::getSpokeCount() const
{

	/*

	With the assumption of a rectangular grid

		n_end_atoms = 2 * n_rows  + 2 * (n_cols-2)	

		n_std_atoms = n_rows * n_cols - n_end_atoms	

	*/	

	int nEndPrims = 2 * numRows + 2 * ( numColumns - 2 ) ;

	int nStdPrims = numRows * numColumns - nEndPrims ;

	return( 3 * nEndPrims + 2 * nStdPrims ) ;


}

M3DPrimitive * M3DQuadFigure::getPrimitivePtr(int primIndex) const
{
    if(primIndex < 0 || primIndex >= primitives.size())
        return NULL;

    return (primitives[primIndex]);
}

void M3DQuadFigure::figuralCoordinates(int index, int & rowIndex, int & colIndex) const
{
    rowIndex = index/numColumns;
	colIndex = index - rowIndex*numColumns;
}

M3DPrimitive * M3DQuadFigure::getPrimitivePtr(int rowIndex, int colIndex) const
{
    int index;

    index = rowIndex * numColumns + colIndex;
    if(index < 0 || index >= primitives.size())
        return NULL;

    return (primitives[index]);
}

int M3DQuadFigure::getPrimitiveID(int rowIndex, int colIndex) const
{
    int index;

    index = rowIndex * numColumns + colIndex;
    if(index < 0 || index >= primitives.size())
        return -1;

    return index;
}

void M3DQuadFigure::setPrimitivePtr(int primIndex, M3DPrimitive * newPrimPtr)
{
    if(primIndex < 0 || primIndex > primitives.size())
        return;

    if(primitives[primIndex] != NULL)
        delete primitives[primIndex];

    primitives[primIndex] = newPrimPtr;

    properties->modified = true;
}

void M3DQuadFigure::setPrimitivePtr(int rowIndex, int colIndex,
                                    M3DPrimitive * newPrimPtr)
{
    int primIndex;

    primIndex = rowIndex * numColumns + colIndex;
    if(primIndex < 0 || primIndex >= primitives.size())
        return;

    if(primitives[primIndex] != NULL)
        delete primitives[primIndex];

    primitives[primIndex] = newPrimPtr;

    properties->modified = true;
}

void M3DQuadFigure::addRow(vector<M3DPrimitive *> row, int rowIndex)
{
    vector<M3DPrimitive *>::iterator it,
                                     index;

    if(rowIndex < 0 || rowIndex > numRows || row.size() != numColumns)
        return;

    index = primitives.begin() + rowIndex * numColumns;
    for(it = row.begin(); it != row.end(); it++)
    {
        primitives.insert(index, *it);
        index++;
    }
#ifdef BINARY
	numRows++;
#endif

    properties->modified = true;
}

void M3DQuadFigure::addColumn(vector<M3DPrimitive *> column, int columnIndex)
{
    vector<M3DPrimitive *>::iterator it,
                                     index;

    if(columnIndex < 0 || columnIndex > numColumns || column.size() != numRows)
        return;

    index = primitives.begin() + columnIndex;
    for(it = column.begin(); it != column.end(); it++)
    {
        primitives.insert(index, *it);
#ifdef BINARY
        index += numColumns + 1;	// "+1" accounts for the atom just inserted
#else
        index += numColumns;
#endif
    }
#ifdef BINARY
	numColumns++;
#endif

    properties->modified = true;
}

void M3DQuadFigure::select()
{
    vector<M3DPrimitive *>::iterator it;
    M3DPrimitive * primPtr;

    for(it = primitives.begin(); it != primitives.end(); it++)
    {
        primPtr = *it;
        if(primPtr != NULL)
            primPtr->select();
    }
}

void M3DQuadFigure::deselect()
{
    vector<M3DPrimitive *>::iterator it;
    M3DPrimitive * primPtr;

    for(it = primitives.begin(); it != primitives.end(); it++)
    {
        primPtr = *it;
        if(primPtr != NULL)
            primPtr->deselect();
    }
}

void M3DQuadFigure::toggleSelected()
{
    vector<M3DPrimitive *>::iterator it;
    M3DPrimitive * primPtr;

    for(it = primitives.begin(); it != primitives.end(); it++)
    {
        primPtr = *it;
        if(primPtr != NULL)
            primPtr->toggleSelected();
    }
}

bool M3DQuadFigure::isSelected()
{
    vector<M3DPrimitive *>::iterator it;
    M3DPrimitive * primPtr;

    for(it = primitives.begin(); it != primitives.end(); it++)
    {
        primPtr = *it;
        if(primPtr != NULL && (! primPtr->isSelected()))
            return false;
    }

    return true;
}

bool M3DQuadFigure::isAnySelected()
{
    vector<M3DPrimitive *>::iterator it;
    M3DPrimitive * primPtr;

    for(it = primitives.begin(); it != primitives.end(); it++)
    {
        primPtr = *it;
        if(primPtr != NULL && primPtr->isSelected())
            return true;
    }

    return false;
}

// Return the number of selected atoms
int M3DQuadFigure::numberSelected()
{
    vector<M3DPrimitive *>::iterator it;
    M3DPrimitive * primPtr;
	int count;

	count = 0;
    for(it = primitives.begin(); it != primitives.end(); it++)
    {
        primPtr = *it;
        if(primPtr != NULL && primPtr->isSelected())
            count++;
    }

    return count;
}

// Return index of any fully selected row; -1 if anything else is found
int M3DQuadFigure::rowSelected()
{
    M3DPrimitive * primPtr;
	int row, col;
	bool full;

	full = false;
    for (row = 0; row < numRows; row++)
    {
        primPtr = getPrimitivePtr(row, 0);
        if (primPtr != NULL && primPtr->isSelected()) {
			full = true;
			for (col = 1; col < numColumns; col++) {
				primPtr = getPrimitivePtr(row, col);
				if (primPtr != NULL && ! primPtr->isSelected()) {
					full = false;
					break;
				}
			}
			if (full)
				break;
		}
    }
	if (full)
		return row;
	else
		return -1;
}

// Return index of any fully selected column; -1 if anything else is found
int M3DQuadFigure::columnSelected()
{
    M3DPrimitive * primPtr;
	int row, col;
	bool full;

	full = false;
    for (col = 0; col < numColumns; col++)
    {
        primPtr = getPrimitivePtr(0, col);
        if (primPtr != NULL && primPtr->isSelected()) {
			full = true;
			for (row = 1; row < numRows; row++) {
				primPtr = getPrimitivePtr(row, col);
				if (primPtr != NULL && ! primPtr->isSelected()) {
					full = false;
					break;
				}
			}
			if (full)
				break;
		}

    }
	if (full)
		return col;
	else
		return -1;
}

// Reverse the indexing along rows
void M3DQuadFigure::reverseRows()	// Will not work with NULL atoms
{
    M3DPrimitive * atom;
    M3DPrimitive * p1;
    M3DPrimitive * p2;
	int row, col;

	if (numColumns <= 1)
		return;

    for (row = 0; row < numRows; row++)
    {
		for (col = 0; col < numColumns/2; col++) {
			p1 = getPrimitivePtr(row, col);

			if (p1->type() == M3D_END_PRIMITIVE) {
				atom = new M3DQuadEndPrimitive();
			}
			else {
				atom = new M3DQuadPrimitive();
			}

			*atom = *p1;
			p2 = getPrimitivePtr(row, numColumns - col - 1);
			*p1 = *p2;
			*p2 = *atom;

			delete atom;
		}
	}
}

// Reverse the indexing along columns
void M3DQuadFigure::reverseColumns()	// Will not work with NULL atoms
{
    M3DPrimitive *atom;
    M3DPrimitive * p1;
    M3DPrimitive * p2;
	int row, col;

	if (numRows <= 1)
		return;

    for (col = 0; col < numColumns; col++)
    {
		for (row = 0; row < numRows/2; row++) {
			p1 = getPrimitivePtr(row, col);

			if (p1->type() == M3D_END_PRIMITIVE) {
				atom = new M3DQuadEndPrimitive();
			}
			else {
				atom = new M3DQuadPrimitive();
			}

			*atom = *p1;
			p2 = getPrimitivePtr(numRows - row - 1, col);
			*p1 = *p2;
			*p2 = *atom;

			delete atom;
		}
	}
}


// Exchange row indexing and column indexing about the major diagonal
void M3DQuadFigure::transpose(int ** map)
{
	int * x_map;
    M3DQuadPrimitive atom;
	int row, col, i, n;
    std::vector<M3DPrimitive *> temp;

	n = numRows*numColumns;
	x_map = new int[n];
	i = 0;

	for (row = 0; row < numRows; row++)
		for (col = 0; col < numColumns; col++)
			x_map[i++] = row + numRows*col;

	// Brute-force transpose
	temp = primitives;
	for (i = 0; i < n; i++)
		primitives[x_map[i]] = temp[i];
	i = numRows;
	numRows = numColumns;
	numColumns = i;

	if (map != NULL)
		*map = x_map;
	else
		delete [] x_map;
}

// Verifies that the spoke ends and extended B-vector end lie inside the unit cube
bool M3DQuadFigure::verifyInBounds() const
{
    M3DPrimitive * atom;
	int atomId;

    for (atomId = 0; atomId != primitives.size(); atomId++)
    {
        atom = primitives[atomId];
        if (atom != NULL) {
			Vector3D involute = atom->getX() + atom->getExtendedB();
			// Could use < and > in place of <= and >= here, but this is a bit safer
			if (involute.getX() <= 0.0 || involute.getX() >= 1.0
			    || involute.getY() <= 0.0 || involute.getY() >= 1.0
			    || involute.getZ() <= 0.0 || involute.getZ() >= 1.0)
					return false;
			involute = atom->getX() + atom->getY0();
			if (involute.getX() <= 0.0 || involute.getX() >= 1.0
			    || involute.getY() <= 0.0 || involute.getY() >= 1.0
			    || involute.getZ() <= 0.0 || involute.getZ() >= 1.0)
					return false;
			involute = atom->getX() + atom->getY1();
			if (involute.getX() <= 0.0 || involute.getX() >= 1.0
			    || involute.getY() <= 0.0 || involute.getY() >= 1.0
			    || involute.getZ() <= 0.0 || involute.getZ() >= 1.0)
					return false;
		}
    }
    return true;
}

Vector3D M3DQuadFigure::getCOG(bool all) const
{
    Vector3D center(0.0, 0.0, 0.0);
    M3DPrimitive * primPtr;
    int numPrimitives,
        count,
        i;

    count = 0;
    numPrimitives = primitives.size();

	if (all) {
		for(i = 0; i < numPrimitives; i++)
		{
			primPtr = primitives[i];

			if(primPtr != NULL)
			{
				center += primPtr->getX();
				count++;
			}
		}
	}
	else {
		for(i = 0; i < numPrimitives; i++)
		{
			primPtr = primitives[i];

			if(primPtr != NULL && primPtr->isSelected())
			{
				center += primPtr->getX();
				count++;
			}
		}
    }

    if(count != 0)
        center /= (double) count;

    return center;
}

void M3DQuadFigure::scaleWidth(double scalefact)
{
    M3DPrimitive * primPtr;
    int numPrimitives,
        i;

    numPrimitives = primitives.size();

    for(i = 0; i < numPrimitives; i++)
    {
        primPtr = primitives[i];

        if(primPtr != NULL && primPtr->isSelected())
        {
            primPtr->scaleBy(scalefact);
            properties->modified = true;
        }
    }
}

void M3DQuadFigure::scaleBy(double scalefact)
{
    scaleBy(scalefact, getCOG());
}

void M3DQuadFigure::scaleBy(double scalefact, const Vector3D &vCenter)
{
    M3DPrimitive * primPtr;
    Vector3D x;
    int numPrimitives, i, n;

    numPrimitives = primitives.size();
	n = 0;
    for(i = 0; i < numPrimitives; i++)
    {
        primPtr = primitives[i];

        if(primPtr != NULL && primPtr->isSelected())
        {
            primPtr->scaleBy(scalefact);

            x = primPtr->getX() - vCenter;
            x *= scalefact;
            primPtr->setX(vCenter + x);

            properties->modified = true;
			n++;
        }
    }
	if (n == numPrimitives)		// Scaling entire figure?
		scaleLandmarks(scalefact, vCenter);
}

void M3DQuadFigure::rotateBy(const Quat &q)
{
    rotateBy(q, getCOG());
}

void M3DQuadFigure::rotateBy(const Quat &q, const Vector3D &vCenter)
{
    M3DPrimitive * primPtr;
    Vector3D x;
    int numPrimitives, i, n;

    numPrimitives = primitives.size();
	n = 0;
    for(i = 0; i < numPrimitives; i++)
    {
        primPtr = primitives[i];
        if(primPtr != NULL && primPtr->isSelected())
        {
            primPtr->rotateBy(q);

            x = primPtr->getX() - vCenter;
            q.rotateVector(x);
            primPtr->setX(vCenter + x);

            properties->modified = true;
			n++;
        }
    }
	if (n == numPrimitives)		// Rotating entire figure?
		rotateLandmarks(q, vCenter);
}

void M3DQuadFigure::translateBy(const Vector3D &vTrans, bool selectAll)
{
    M3DPrimitive * primPtr;
    int numPrimitives, i, n;

    numPrimitives = primitives.size();
	n = 0;
    for(i = 0; i < numPrimitives; i++)
    {
        primPtr = primitives[i];
        if(primPtr != NULL && (selectAll || primPtr->isSelected()))
        {
            primPtr->translateBy(vTrans);
            properties->modified = true;
			n++;
        }
    }
	if (n == numPrimitives)		// Moving entire figure?
		translateLandmarks(vTrans);
}

void M3DQuadFigure::deleteAll()
{
    int numPrimitives,
        i;

    numPrimitives = primitives.size();

    for(i = 0; i < numPrimitives; i++)
    {
        if(primitives[i] != NULL) {
            delete primitives[i];
			primitives[i] = NULL;
		}
    }

    primitives.clear();

    numRows = 0;
    numColumns = 0;

	if (boundary != NULL)
		delete boundary;
	boundary = NULL;

    properties->modified = true;
}

//#ifdef BINARY

// REGULARITY CALCULATION

#define WEST i,j-1
#define NORTH i+1,j
#define EAST i,j+1
#define SOUTH i-1,j
#define NE i+1,j+1
#define SE i-1,j+1
#define NW i+1,j-1
#define SW i-1,j-1

// This is an accessor for the regularity calculations done in M3DFigure.  In a quad-mesh, there
// are several different ways of defining how many neighbors a given atom has.  These are enum'd in
// M3DQuadFigure's header.
void M3DQuadFigure::getAtomsNeighbors(int index, int & numNeighbors, M3DPrimitive * neighbors[],
	PrimNeighborhoodDefn PND)
{
	int nRows = getRowCount();
	int nCols = getColumnCount();

	int i = index/nCols;  // Which row you're on
	int j = index%nCols;  // Which column

	// Internal
	if (i > 0 && i < nRows - 1 && j > 0 && j < nCols - 1) {
		numNeighbors = 4;
		neighbors[0] = getPrimitivePtr(EAST);
		neighbors[1] = getPrimitivePtr(NORTH);
		neighbors[2] = getPrimitivePtr(WEST);
		neighbors[3] = getPrimitivePtr(SOUTH);
	}

	// Corners
	else if (i == 0 && j == 0) {
		if (PND == ALL_FIRST_NEIGHBORS || PND == EDGES_HAVE_2_NEIGHBORS) {
			numNeighbors = 2;
			neighbors[0] = getPrimitivePtr(EAST);
			neighbors[1] = getPrimitivePtr(NORTH);
		}
		else if (PND == CORNERS_HAVE_3_NEIGHBORS){
			numNeighbors = 3;
			neighbors[0] = getPrimitivePtr(EAST);
			neighbors[1] = getPrimitivePtr(NORTH);
			neighbors[2] = getPrimitivePtr(NE);		
		}
		else
			numNeighbors = 0;
	}
	else if (i == 0 && j == nCols - 1) {
		if (PND == ALL_FIRST_NEIGHBORS || PND == EDGES_HAVE_2_NEIGHBORS) {
			numNeighbors = 2;
			neighbors[0] = getPrimitivePtr(NORTH);
			neighbors[1] = getPrimitivePtr(WEST);
		}	
		else if (PND == CORNERS_HAVE_3_NEIGHBORS){
			numNeighbors = 3;
			neighbors[0] = getPrimitivePtr(NORTH);
			neighbors[1] = getPrimitivePtr(WEST);
			neighbors[2] = getPrimitivePtr(NW);
		}
		else
			numNeighbors = 0;
	}
	else if (i == nRows - 1 && j == 0) {
		if (PND == ALL_FIRST_NEIGHBORS || PND == EDGES_HAVE_2_NEIGHBORS) {
			numNeighbors = 2;
			neighbors[0] = getPrimitivePtr(SOUTH);
			neighbors[1] = getPrimitivePtr(EAST);
		}
		else if (PND == CORNERS_HAVE_3_NEIGHBORS){
			numNeighbors = 3;
			neighbors[0] = getPrimitivePtr(SOUTH);
			neighbors[1] = getPrimitivePtr(EAST);
			neighbors[2] = getPrimitivePtr(SE);
		}
		else
			numNeighbors = 0;
	}
	else if (i == nRows - 1 && j == nCols - 1) {
		if (PND == ALL_FIRST_NEIGHBORS || PND == EDGES_HAVE_2_NEIGHBORS) {
			numNeighbors = 2;
			neighbors[0] = getPrimitivePtr(WEST);
			neighbors[1] = getPrimitivePtr(SOUTH);
		}
		else if (PND == CORNERS_HAVE_3_NEIGHBORS){
			numNeighbors = 3;
			neighbors[0] = getPrimitivePtr(WEST);
			neighbors[1] = getPrimitivePtr(SOUTH);
			neighbors[2] = getPrimitivePtr(SW);
		}
		else
			numNeighbors = 0;
	}

	// Edges
	else if (i == 0) {
		numNeighbors = 2;
		neighbors[0] = getPrimitivePtr(EAST);
		neighbors[1] = getPrimitivePtr(WEST);

		if (PND == ALL_FIRST_NEIGHBORS || PND == CORNERS_HAVE_3_NEIGHBORS || PND == PIN_CORNERS) {
			neighbors[2] = getPrimitivePtr(NORTH);
			numNeighbors = 3;
		}
	}
	else if (j == 0) {
		numNeighbors = 2;
		neighbors[0] = getPrimitivePtr(NORTH);
		neighbors[1] = getPrimitivePtr(SOUTH);

		if (PND == ALL_FIRST_NEIGHBORS || PND == CORNERS_HAVE_3_NEIGHBORS || PND == PIN_CORNERS) {
			neighbors[2] = getPrimitivePtr(EAST);
			numNeighbors = 3;
		}
	}
	else if (i == nRows - 1) {
		numNeighbors = 2;
		neighbors[0] = getPrimitivePtr(EAST);
		neighbors[1] = getPrimitivePtr(WEST);

		if (PND == ALL_FIRST_NEIGHBORS || PND == CORNERS_HAVE_3_NEIGHBORS || PND == PIN_CORNERS) {
			neighbors[2] = getPrimitivePtr(SOUTH);
			numNeighbors = 3;
		}
	}
	else if (j == nCols - 1) {
		numNeighbors = 2;
		neighbors[0] = getPrimitivePtr(NORTH);
		neighbors[1] = getPrimitivePtr(SOUTH);

		if (PND == ALL_FIRST_NEIGHBORS || PND == CORNERS_HAVE_3_NEIGHBORS || PND == PIN_CORNERS) {
			neighbors[2] = getPrimitivePtr(WEST);
			numNeighbors = 3;
		}
	}
}

//#endif


// Persistent object support
void M3DQuadFigure::writeFigure( const char * figureStr, Registry& registry )
{
    int numRows,
        numColumns;
    int i, j;
    const char * figureName;

    char * format;
    char * subordinateStr;

    const float * color;

    if(figureStr == NULL)
        return;

    format = new char[strlen(figureStr) + 128];
    subordinateStr = new char[strlen(figureStr) + 128];

    strcpy(format, figureStr);
    strcat(format, ".%s");

    numRows = getRowCount();
    numColumns = getColumnCount();

    figureName = getName();
    if(figureName != NULL)
        registry.setStringValue(format, (char *) figureName, "name");
    registry.setStringValue(format, friendlyFigureName, "type");

    registry.setIntValue(format, numRows, "numRows");
    registry.setIntValue(format, numColumns, "numColumns");

    registry.setBooleanValue(format, isPositiveSpace(), "positiveSpace");
    registry.setBooleanValue(format, hasPositivePolarity(), "positivePolarity");

    registry.setIntValue(format, getTolerance(), "smoothness");

    color = getColor();
    registry.setDoubleValue(format, color[0], "color.red");
    registry.setDoubleValue(format, color[1], "color.green");
    registry.setDoubleValue(format, color[2], "color.blue");

#ifdef LM_METHOD_OBJ
	const int numLandmarks = getLandmarkCount();
    registry.setIntValue(format, numLandmarks, "numLandmarks");
	if (numLandmarks > 0) {
        sprintf(subordinateStr, format, "landmark[%d].%s");
		for (i = 0; i < numLandmarks; i++) {
			const char * name = getLandmarkName(i);
			if (name != NULL && name[0] != '\0')
				registry.setStringValue(subordinateStr, (char *) name, i, "name");
            int atomIndex = getLandmarkAtomIndex(i);
			registry.setIntValue(subordinateStr, atomIndex, i, "atomIndex");
            double atomT = getLandmarkAtomT(i);
			registry.setDoubleValue(subordinateStr, atomT, i, "atomT");
            double atomU = getLandmarkAtomU(i);
			registry.setDoubleValue(subordinateStr, atomU, i, "atomU");
            double atomV = getLandmarkAtomV(i);
			registry.setDoubleValue(subordinateStr, atomV, i, "atomV");
		}
	}
#endif

    sprintf(subordinateStr, format, "primitive[%d][%d]");
    for(i = 0; i < numRows; i++)
    {
        for(j = 0; j < numColumns; j++)
        {
            getPrimitivePtr(i,j)->writePrimitive(registry, subordinateStr, i, j);

            //cout<<"-------------------------------------------------"<<registry.getDoubleValue("model.figure[1].primitive[%d][%d].r[0]",i,j);
        }
    }

    InterfiguralConstraints & ifc = constraints();
    int n = ifc.size();
    if (n > 0) {
        registry.setIntValue(format, n, "numConstraints");
        sprintf(subordinateStr, format, "constraint[%d]");
        for (i = 0; i < n; i++) {
            sprintf(format, subordinateStr, i);
            strcat(format, ".%s");
            registry.setIntValue(format, ifc.figure(i), "figure");
            registry.setDoubleValue(format, ifc.distance(i), "distance");
        }
    }

	// FIXME
	// Write boundary information
//	writeBoundary(figureStr, getBoundaryPtr());

	// Write the stackedImageMask, if the model was optimized against a 
	// stacked image.
	sprintf(format, figureStr);
	strcat(format, ".%s");

	if (getStackedImageMask() != 0)
		registry.setIntValue(format, (int) getStackedImageMask(),
			"stackedImageMask");

	// Write figure statistics
	//if (getFigureStatsPtr() != NULL)
	//	writeFigureStats(figureStr, getFigureStatsPtr());

    delete [] format;
    delete [] subordinateStr;


    //test if write into registry correctly?
    //cout<<"-------------------------------------------------"<<registry.getDoubleValue(format, 0.0, "r[0]")<<endl;






}

M3DFigure * M3DQuadFigure::readFigure(int figureIndex, Registry& registry )
{
    M3DQuadFigure * figure;
    const char * name;
    int numRows,
        numColumns;
    int i, j;
    float color[3];
	int numLandmarks;
    bool positiveSpace;
    bool positivePolarity;
    int tolerance;

    name = registry.getStringValue("model.figure[%d].name", NULL, figureIndex);

    color[0] = registry.getDoubleValue("model.figure[%d].color.red",
		DEFAULT_FIGURE_COLOR[0], figureIndex);
    color[1] = registry.getDoubleValue("model.figure[%d].color.green",
		DEFAULT_FIGURE_COLOR[1], figureIndex);
    color[2] = registry.getDoubleValue("model.figure[%d].color.blue",
		DEFAULT_FIGURE_COLOR[2], figureIndex);

    numRows = registry.getIntValue("model.figure[%d].numRows", 0, figureIndex);
    numColumns = registry.getIntValue("model.figure[%d].numColumns", 0, figureIndex);

    positiveSpace = registry.getBooleanValue("model.figure[%d].positiveSpace",
		true, figureIndex);
    positivePolarity = registry.getBooleanValue("model.figure[%d].positivePolarity",
		true, figureIndex);

    tolerance = registry.getIntValue("model.figure[%d].smoothness",
		M3DFigure::getDefaultSurfaceTolerance(), figureIndex);

    figure = new M3DQuadFigure(numRows, numColumns);
    if(figure == NULL)
        return NULL;

    figure->setName(name);
    figure->setColor(color);
    figure->setPositiveSpace(positiveSpace);
    figure->setPositivePolarity(positivePolarity);
    figure->setTolerance(tolerance);

    for(i = 0; i < numRows; i++)
    {
        for(j = 0; j < numColumns; j++)
        {
            figure->setPrimitivePtr(i, j, M3DQuadPrimitive::readPrimitive(registry, "model.figure[%d].primitive[%d][%d]", figureIndex, i, j));
        }
    }

	// read landmarks after atoms, since landmarks refer to atoms
#ifdef LM_METHOD_OBJ
    numLandmarks = registry.getIntValue("model.figure[%d].numLandmarks", 0, figureIndex);
	if (numLandmarks > 0) {
		for (i = 0; i < numLandmarks; i++) {
            int atomIndex = 
				registry.getIntValue("model.figure[%d].landmark[%d].atomIndex", -1, figureIndex, i);
            double atomT = 
				registry.getDoubleValue("model.figure[%d].landmark[%d].atomT", 0, figureIndex, i);
            double atomU = 
				registry.getDoubleValue("model.figure[%d].landmark[%d].atomU", -1, figureIndex, i);
            double atomV = 
				registry.getDoubleValue("model.figure[%d].landmark[%d].atomV", -1, figureIndex, i);
			const char * atomName = registry.getStringValue("model.figure[%d].landmark[%d].name",
				"", figureIndex, i);
			char *permName = new char[strlen(atomName) + 1];
			strcpy(permName, atomName);
			if ((atomU != -1)  && (atomV != -1)) {
			  figure->addLandmark(atomU, atomV, permName, atomT);
            }
			else { 
			  figure->addLandmark(atomIndex, permName, atomT);
            }

			// warn if landmark placed on b-spoke of non-end primitive, but accept landmark anyway.
			if (globalVerbosity > -1)
			{
				if (atomT == 0 && atomIndex >= 0)	// b-spoke and valid atomIndex
				{
					M3DPrimitive *atom = figure->getPrimitivePtr(atomIndex);

					// valid atom but not an end atom
					if (atom && atom->type() != M3D_END_PRIMITIVE)
						cout << "WARNING: invalid landmark found on the bisector spoke of an internal atom: "
							 << "lm on atom #" << atomIndex << "=" << atomName << ", figure=" << name
							 << endl;
				}
			}
		}
	}
#endif

    int n = registry.getIntValue("model.figure[%d].numConstraints", 0, figureIndex);
    if (n > 0) {
        InterfiguralConstraints & ifc = figure->constraints();
        for (i = 0; i < n; i++) {
            j = registry.getIntValue("model.figure[%d].constraint[%d].figure",
                -1, figureIndex, i);
            double d = registry.getDoubleValue("model.figure[%d].constraint[%d].distance",
                0.0, figureIndex, i);
            ifc.addFigure(j, d);
        }
    }

	// Pick up stacked image info here, if it is available.
	// The figure may contain a tag called stackedImageMask, an integer
	// between 0 and 15, that indicates which bit(s) of the image this
	// figure's binary image should be taken from.  This is used in Match.
	int mask = registry.getIntValue(
	    "model.figure[%d].stackedImageMask", 0, figureIndex);
	figure->setStackedImageMask(mask);

	// FIXME
//	figure->setFigureStatsPtr(readFigureStats("model.figure[%d]",
//	    figureIndex));

    return figure;
}






const double pnorm = 10.0;

/**
 * A helper function that returns the total curvature at t for the section of the medial curve
 * between prim0 (t=0) and prim1 (t=1).
 */
inline double totalCurvatureSquared( const M3DQuadPrimitive* prim0, const M3DQuadPrimitive* prim1, const double t)
{
	const Vector3D p0	= prim0->getX();
	const Vector3D p1	= prim1->getX();
	//project the  deltaP along u/v direction to the tangent plane  indicated by the atom normal
	//
	Vector3D deltaP = p1-p0;
	Vector3D d0	= ( deltaP -  (deltaP * prim0->getN()) * prim0->getN() );
	d0.normalize();
    d0 = d0 * deltaP.norm();
	
	Vector3D d1	= ( deltaP -  (deltaP * prim1->getN()) * prim1->getN() );
    d1.normalize();
 	d1 = d1 * deltaP.norm();

   // 1st, 2nd and 3rd order derivatives of the interpolated curve.
//	const Vector3D g	= (t*t * ( 2.0*t-3) + 1.0) * p0
//						+ (t*t * (-2.0*t+3)      ) * p1
//						+ (t * (t * (t-2.0) + 1.0)) * d0
//						+ (t * (t * (t - 1.0) )  ) * d1;
//	printf("pts = [ pts; %f, %f, %f];\n", g.getX(), g.getY(), g.getZ());
	const Vector3D dg	= (-6*t+6*t*t) * p0 + (6*t-6*t*t) * p1
						+ (1-4*t+3*t*t) * d0 + (-2*t+3*t*t)	* d1;
	const Vector3D d2g	= (-6+12*t) * p0 + (6-12*t)*p1
						+ (-4+6*t)*d0 + (-2+6*t)*d1;
	const Vector3D d3g	= 12 * p0 - 12 * p1 + 6 * d0 + 6 * d1;

	const double normdg		= dg.norm();
	const double normdgxd2g	= dg.cross(d2g).norm();
	const double kappa	=  normdgxd2g / (normdg*normdg*normdg);
//	cout <<"kappa = "<<kappa<<"\n"<<endl;
//	const double tau	= (fabs(normdgxd2g) <= EPSILON) ?
//		0.0 : dg * (d2g.cross(d3g)) / (normdgxd2g * normdgxd2g) ;
	const double tau	= 0.0;

	return pow(fabs(kappa),pnorm) + pow(fabs(tau),pnorm);
}

/**
 * Reeturns the total curvature for the section of the medial curve
 * between prim0 (t=0) and prim1 (t=1).
 */

inline double  segmentCurviness(const M3DQuadPrimitive * prim0, const M3DQuadPrimitive * prim1)
{
	// Evaluate the total using Simpson's 1/3rd rule.
	
	const int n	= 20;				// number of points of evaluation = 1 + 2*n
	const double h	= 1.0/(2.0*n);			// spacing between points.
	const double x0	= totalCurvatureSquared(prim0,prim1,0.0);
	double segmentCurviness	= x0;
	for( int iStep = 1; iStep != n; ++iStep ) {
		const double x2i	= totalCurvatureSquared(prim0,prim1, 2*iStep*h);
		const double x2ip1	= totalCurvatureSquared(prim0,prim1, (2*iStep+1)*h);
		segmentCurviness += 4 * x2ip1 + 2 * x2i;
	}
	const double x2n	= totalCurvatureSquared(prim0,prim1,1.0);
	segmentCurviness	+= x2n;
	segmentCurviness = segmentCurviness * h/3.0;
//	cout <<"segment Curviness: "<<segmentCurviness<<"\n"<<endl;
	return segmentCurviness;

}



/**
 * This function returns the curviness of the medial sheet.
 *
 * @param	none
 * @return	the curviness of the medial sheet.
 */
double M3DQuadFigure::curviness() const
{	const M3DQuadPrimitive* prim0;
	const M3DQuadPrimitive* prim1;
	double totalCurviness	= 0.0;
	double curvinessU = 0.0;
	double curvinessV = 0.0;
	int i, j ;	
	//double meanRadius;
	double curviness;
	double approxCurveLength;
	
//along the u directions:

	for ( i = 0; i < numRows; i++){

		curviness = 0.0;
		approxCurveLength = 0.0;

		prim0 = dynamic_cast<const M3DQuadPrimitive*>(getPrimitivePtr(i,0));
		//meanRadius = log(prim0->getR());

		for( j = 1; j< numColumns; j++) {

			prim1	= dynamic_cast<const M3DQuadPrimitive*>(getPrimitivePtr(i,j));
			approxCurveLength	+= (prim1->getX() - prim0->getX()).norm();

			//meanRadius	+= log(prim1->getR());
			curviness += segmentCurviness(prim0, prim1);
			
			prim0	= prim1;
		}

		//meanRadius	/= getPrimitiveCount();
		//meanRadius	= exp(meanRadius);

		//	cout << "t: " << totalCurviness << ", mr: " << meanRadius << endl;
		//curviness = meanRadius * pow(curviness, 1.0/pnorm);
		curviness = approxCurveLength * pow(curviness, 1.0/pnorm);
		curvinessU += pow(curviness,pnorm);

	}

//cout <<"curvinessU ="<<curvinessU<<"\n"<<endl;
//along the v directions:
	
	for( j = 0; j< numColumns ; j++) {

		curviness = 0.0;
		approxCurveLength = 0.0;

		prim0 = dynamic_cast<const M3DQuadPrimitive*>(getPrimitivePtr(0,j));
		//meanRadius = log(prim0->getR());

		for( i = 1; i< numRows; i++) {

			prim1	= dynamic_cast<const M3DQuadPrimitive*>(getPrimitivePtr(i,j));
			approxCurveLength	+= (prim1->getX() - prim0->getX()).norm();

			//meanRadius	+= log(prim1->getR());
			curviness += segmentCurviness(prim0, prim1);
		
			prim0	= prim1;
		}

		//meanRadius	/= getPrimitiveCount();
		//meanRadius	= exp(meanRadius);

		//	cout << "t: " << totalCurviness << ", mr: " << meanRadius << endl;
		//curviness = meanRadius * pow(curviness, 1.0/pnorm);
		curviness = approxCurveLength * pow(curviness, 1.0/pnorm);

		curvinessV += pow(curviness,pnorm);

	}
	
//	cout <<"curvinessV ="<<curvinessV<<"\n"<<endl;
	totalCurviness = pow(curvinessU + curvinessV,1.0/pnorm) / getPrimitiveCount();
	return totalCurviness * totalCurviness ;
}

// Dibyendu

/* 
Function to convert a figure into a vector of parameters as required by the CPNS program
The following is the format: [positions ; radii ; spoke directions]
Standard primitive: (Px, Py, Pz, r1, r2, U1x, U1y, U1z, U2x, U2y, U2z)
End primitive: (Px, Py, Pz, r1, r2, r3, U1x, U1y, U1z, U2x, U2y, U2z, U3x, U3y, U3z)
*/ 

double*  M3DQuadFigure::vectorize() {

	/* 

	If no. of rows = n, no. of cols = m

	Data per object:

		No. of end atoms = 2 * n  + 2 * (m-2)
		No. of std atoms = m * n - ( 2 * n  + 2 * (m-2) )

		No. of entities per std atom = 11 
		No. of entities per end atom = 15 

		Total entries = 15 * nEndAtoms + 11 * nStdAtoms + 2 ;

		2 values are for nRows, nCols

	*/ 


	int nRows = getRowCount() ;
	int nCols = getColumnCount() ;

	int nEndPrims = 2 * nRows + 2 * ( nCols - 2 ) ;
	int nStdPrims = nRows * nCols - nEndPrims ;

	int nValues = 15 * nEndPrims + 11 * nStdPrims + 2 ;

	double *  values = new double [ nValues ] ;

	for( int ii = 0 ; ii < nValues ; ii ++ )
		values[ii] = 0 ;

	for( int ii = 0 ; ii < nRows ; ii ++ ) {
		for( int jj = 0 ; jj < nCols ; jj ++ ) {

			M3DPrimitive* p = getPrimitivePtr(ii, jj) ;

			// calculating the number of end atoms before

			int nEndAtomsBefore = 0 ;

			if( ii > 0 ) {

				nEndAtomsBefore += nCols + 2 * ( ii - 1 ) ;

				if( jj > 0 )  {

					if( ii < nRows-1 )
						nEndAtomsBefore += 1 ;
					else
						nEndAtomsBefore += jj ;
				}
			}
			else {
				if( jj > 0 ) 
					nEndAtomsBefore += jj ;
			}

			// calculate the number of standard atoms before

			int nStdAtomsBefore = ( ii * nCols + jj ) - nEndAtomsBefore ;

			// position start index

			int pStartIndex =	3 * ( ii * nCols + jj ) ;

			// radii start index

			int rStartIndex = 		3 * nRows * nCols		// total positions
								+	3 * nEndAtomsBefore	+ 2 * nStdAtomsBefore ;

			// unit direction start index

			int uStartIndex =		3 * nRows * nCols				// total positions
								+	3 * nEndPrims + 2 * nStdPrims	// total radii
								+	9 * nEndAtomsBefore + 6 * nStdAtomsBefore ;

			// pointers for the std and end atoms

			M3DQuadPrimitive* stdP = NULL ;

			M3DQuadEndPrimitive* endP = NULL ;

			switch( p->type() ) {

				case M3D_STANDARD_PRIMITIVE:

					stdP = dynamic_cast <M3DQuadPrimitive*> ( p ) ;

					values[ pStartIndex + 0 ] = stdP->getX().getX() ;
					values[ pStartIndex + 1 ] = stdP->getX().getY() ;
					values[ pStartIndex + 2 ] = stdP->getX().getZ() ;
					
					values[ rStartIndex + 0 ] = stdP->getR0() ;
					values[ rStartIndex + 1 ] = stdP->getR1() ;

					values[ uStartIndex + 0 ] = stdP->getU0().getX() ;
					values[ uStartIndex + 1 ] = stdP->getU0().getY() ;
					values[ uStartIndex + 2 ] = stdP->getU0().getZ() ;
					values[ uStartIndex + 3 ] = stdP->getU1().getX() ;
					values[ uStartIndex + 4 ] = stdP->getU1().getY() ;
					values[ uStartIndex + 5 ] = stdP->getU1().getZ() ;

					break ;

				case M3D_END_PRIMITIVE:

					endP = dynamic_cast <M3DQuadEndPrimitive*> ( p ) ;

					values[ pStartIndex + 0 ] = endP->getX().getX() ;
					values[ pStartIndex + 1 ] = endP->getX().getY() ;
					values[ pStartIndex + 2 ] = endP->getX().getZ() ;
					
					values[ rStartIndex + 0 ] = endP->getR0() ;
					values[ rStartIndex + 1 ] = endP->getR1() ;
					values[ rStartIndex + 2 ] = endP->getREnd() ;

					values[ uStartIndex + 0 ] = endP->getU0().getX() ;
					values[ uStartIndex + 1 ] = endP->getU0().getY() ;
					values[ uStartIndex + 2 ] = endP->getU0().getZ() ;
					values[ uStartIndex + 3 ] = endP->getU1().getX() ;
					values[ uStartIndex + 4 ] = endP->getU1().getY() ;
					values[ uStartIndex + 5 ] = endP->getU1().getZ() ;
					values[ uStartIndex + 6 ] = endP->getUEnd().getX() ;
					values[ uStartIndex + 7 ] = endP->getUEnd().getY() ;
					values[ uStartIndex + 8 ] = endP->getUEnd().getZ() ;

					break ;

			}
		}
	}

	return( values ) ;

}







