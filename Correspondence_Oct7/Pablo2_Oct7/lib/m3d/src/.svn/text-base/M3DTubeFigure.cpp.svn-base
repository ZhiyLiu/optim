//#include <typeinfo.h>
#include <math.h>
#include <stdio.h>
#include "M3DTubeFigure.h"
#include "M3DEndPrimitive.h"
#include "M3DAtomPredictor.h"
#include <assert.h>

#include <typeinfo>

#ifndef EPSILON
#define EPSILON 1e-8
#endif

#define max(a,b) ((a) >= (b) ? (a) : (b))

//#define DEBUG

const char* M3DTubeFigure::friendlyFigureName	= "TubeFigure";

using namespace std;

M3DTubeFigure::M3DTubeFigure() : M3DFigure()
{
#ifdef DEBUG
	cout << "M3DTubeFigure::M3DTubeFigure()" << endl;
#endif
    numColumns		= 0;
	numberOfSpokes	= M3DTubePrimitive::DEFAULT_NUMBER_OF_SPOKES;

	boundary		= new SubdivBoundary();
}

M3DTubeFigure::M3DTubeFigure(int _numColumns, int _numberOfSpokes)
	: M3DFigure()
{
#ifdef DEBUG
	cout << "M3DTubeFigure::M3DTubeFigure(int)" << endl;
#endif
	if(_numberOfSpokes < 3 )
	{
		cout << "M3DTubeFigure(int,int): numberOfSpokes(" << _numberOfSpokes << ") should be at least 3." << endl;
		_numberOfSpokes	= M3DTubePrimitive::DEFAULT_NUMBER_OF_SPOKES;
		numColumns		= 0;
		return;
	}
	else if(_numColumns <= 0)
    {
		numColumns		= 0;
		numberOfSpokes	= M3DTubePrimitive::DEFAULT_NUMBER_OF_SPOKES;
        return;
    }
    else
    {
        numColumns		= _numColumns;
		numberOfSpokes	= _numberOfSpokes;
    }

    primitives.resize(numColumns);
    for(int i = 0; i < numColumns; i++)
		primitives[i] = NULL;

	boundary = new SubdivBoundary();
}

M3DTubeFigure::M3DTubeFigure(const M3DTubeFigure & figure) : M3DFigure(figure)
{
    M3DPrimitive * primPtr,
                 * newPrim;
#ifdef DEBUG
	cout << "M3DTubeFigure::M3DTubeFigure(const M3DTubeFigure &)" << endl;
#endif
    numColumns		= figure.numColumns;
	numberOfSpokes	= figure.numberOfSpokes;

    for(int i = 0; i < numColumns; i++)
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
		boundary = new SubdivBoundary();

}

M3DFigure & M3DTubeFigure::operator = (const M3DFigure & unknown_figure)
{
    M3DPrimitive * primPtr,
                 * newPrim;

	assert( typeid(*this) == typeid(unknown_figure) );
	const M3DTubeFigure& figure	= dynamic_cast<const M3DTubeFigure&>(unknown_figure);

#ifdef DEBUG
	cout << "M3DTubeFigure::operator=()" << endl;
#endif
    deleteAll();

	M3DFigure::operator =(static_cast<const M3DFigure&>(figure));

    numColumns		= figure.numColumns;
	numberOfSpokes	= figure.numberOfSpokes;

    for(int i = 0; i < numColumns; i++)
    {
        primPtr = figure.primitives[i];
        if(primPtr != NULL)
            newPrim = primPtr->copyPtr();
        else
            newPrim = NULL;

        primitives.insert(primitives.end(), newPrim);
    }

	if (figure.getBoundaryPtr() != NULL)
		this->boundary = new SubdivBoundary(figure.getBoundaryPtr());

    return (*this);
}

M3DTubeFigure::~M3DTubeFigure()
{
#ifdef DEBUG
	cout << "M3DTubeFigure::~M3DTubeFigure()" << endl;
#endif
    deleteAll();
}

void M3DTubeFigure::print(bool dump, ostream & out, int markedPrimIndex) const
{
    M3DFigure::print(dump, out, markedPrimIndex);

	out << "Shape: " << numColumns << " columns\n";
	for(int j = 0; j < numColumns; j++)
	{
		out << "    Primitive(" << j << ')';
		out << " at 0x" << hex << primitives[j] << dec;
		if (dump) {
			out << ":\n";
			primitives[j]->print(out, "\t");
		}
		else
			out << '\n';
	}

	int n = ifconstraints.size();
	if (n > 0) {
		out << "Constrained figures: ";
		for (int i = 0; i < n; i++)
			out << ifconstraints.figure(i) << '(' << ifconstraints.distance(i) << ") ";
		out << '\n';
	}

	n = inverse_constraints.size();
	if (n > 0) {
		out << "Governing figures: ";
		for (int i = 0; i < n; i++)
			out << inverse_constraints.figure(i) << '(' << inverse_constraints.distance(i) << ") ";
		out << '\n';
	}
	out << flush;
}

void M3DTubeFigure::makeRectangularStockFigure(int cols)
{
    M3DPrimitive * primPtr;
    Quat q;

    if(cols <= 0)
        return;

    deleteAll();

	double* radii;
	radii = new double[M3DTubePrimitive::DEFAULT_NUMBER_OF_SPOKES];
	for( int i = 0; i != M3DTubePrimitive::DEFAULT_NUMBER_OF_SPOKES; ++i ) {
		radii[i] = 0.5;
	}

	for(int j = 0; j < cols; j++)
	{
		if(j == 0 || j == cols - 1)
			primPtr = new M3DTubeEndPrimitive((double) j, 0.0, 0.0, 0.5, 0, radii,
				Vector3D(0, 1.0, 0.0), Vector3D(1.0, 0.0, 0.0), M3DTubePrimitive::DEFAULT_NUMBER_OF_SPOKES);
		else
			primPtr = new M3DTubePrimitive((double) j, 0.0, 0.0, 0, radii,
				Vector3D(0, 1.0, 0.0), Vector3D(1.0, 0.0, 0.0), M3DTubePrimitive::DEFAULT_NUMBER_OF_SPOKES);
		primitives.insert(primitives.end(), primPtr);
	}

    numColumns = cols;

    select();

    scaleBy(1.0 / (double) (numColumns + 1));

    q.setAxisAngle(Vector3D(0.0, 0.0, 1.0), R_HALF_PI);
    rotateBy(q);

    translateBy(Vector3D(0.5, 0.5, 0.5) - getCOG());

	//
	// Correct orientation of tube tangent vectors.
	//
	orientTubeTangents();

    properties->modified = true;

	delete[] radii;
}

void M3DTubeFigure::makeLandmarkStockFigure(int cols)
{
	makeEllipticalStockFigure(cols);

	// set all atoms' positions to the position of the first atom so that all
	// will be either selected or deselected together. This makes them act
	// as a single point in space, suitable for landmarks.  -gst 20041103
	Vector3D posFirstPrim = getPrimitivePtr(0)->getX();
	for(int j = 0; j < cols; j++)
	{
		getPrimitivePtr(j)->setX(posFirstPrim);
	}
}

void M3DTubeFigure::makeEllipticalStockFigure(int cols)
{
    M3DPrimitive * primPtr;
    Quat q;
	int j;

    if(cols <= 0)
        return;

    deleteAll();

	double* radii;
	radii = new double[M3DTubePrimitive::DEFAULT_NUMBER_OF_SPOKES];
	for( int i = 0; i != M3DTubePrimitive::DEFAULT_NUMBER_OF_SPOKES; ++i ) {
		radii[i] = 0.5;
	}

	for(j = 0; j < cols/2; j++)
	{
		if(j == 0 || j == cols - 1)
			primPtr = new M3DTubeEndPrimitive((double) j, 0.0, 0.0, 0.5, 0.5*cos(4.0/3.0*R_HALF_PI), radii,
				Vector3D(0, 1.0, 0.0), Vector3D(1.0, 0.0, 0.0), M3DTubePrimitive::DEFAULT_NUMBER_OF_SPOKES);
		else
			primPtr = new M3DTubePrimitive((double) j, 0.0, 0.0, 0.5*cos(3.5/3.0*R_HALF_PI), radii,
				Vector3D(0, 1.0, 0.0), Vector3D(1.0, 0.0, 0.0), M3DTubePrimitive::DEFAULT_NUMBER_OF_SPOKES);
		primitives.insert(primitives.end(), primPtr);
	}
//	q.setAxisAngle(Vector3D(1.0,0.0,0.0),R_PI);
	for(; j < cols; j++)
	{
		if(j == 0 || j == cols - 1)
			primPtr = new M3DTubeEndPrimitive((double) j, 0.0, 0.0, 0.5, 0.5*cos(2.0/3.0*R_HALF_PI), radii,
				Vector3D(0, 1.0, 0.0), Vector3D(1.0, 0.0, 0.0), M3DTubePrimitive::DEFAULT_NUMBER_OF_SPOKES);
		else
			primPtr = new M3DTubePrimitive((double) j, 0.0, 0.0, 0.5*cos(2.5/3.0*R_HALF_PI), radii,
				Vector3D(0, 1.0, 0.0), Vector3D(1.0, 0.0, 0.0), M3DTubePrimitive::DEFAULT_NUMBER_OF_SPOKES);
		primitives.insert(primitives.end(), primPtr);
	}

    numColumns = cols;

    select();

    scaleBy(1.0 / (double) (numColumns + 1));

    q.setAxisAngle(Vector3D(0.0, 0.0, 1.0), R_HALF_PI);
    rotateBy(q);

    translateBy(Vector3D(0.5, 0.5, 0.5) - getCOG());

	//
	// Correct orientation of tube tangent vectors.
	//
	orientTubeTangents();
    
	properties->modified = true;
    cout << "Added stock figure with " <<  cols << " columns" << endl;
	delete[] radii;
}



bool M3DTubeFigure::orientTubeTangents() {
	bool changed	= false;
	Vector3D dx;
	for( int i = 0; i < primitives.size(); ++i ) {
		if( i < primitives.size()-1 ) {
			dx	= getPrimitivePtr(i+1)->getX() - getPrimitivePtr(i)->getX();
		}
		// Very approximate.
		if( dx * primitives[i]->getB() <= 0 ) {
			changed	= true;
			Quat q;
			q.setAxisAngle( getPrimitivePtr(i)->getN(), R_PI );
			primitives[i]->rotateBy( q );
			// B vector pointing along negative arc length?
			primitives[i]->setTheta( R_PI - primitives[i]->getTheta() );
		}
	}
	return changed;
}


//
// Just a helper function for aligning an atom with the previous atom.
// Takes a quaternion representing the orientation of the "previous atom"
// and just the bisector/tangent vector for the current atom.
// Returns the appropriate quaternion for the current atom.
//
// See the code in @TubePrimitive/alignTube.m in matlab for explanation.
//
Quat phiAlignAtom( const Quat& prevQ, const Vector3D& T )
{
	Vector3D prevT, dT;
	Vector3D normal, prevNormal;
	Quat R, P;
	prevT.set(1,0,0);
	prevQ.rotateVector(prevT);
	dT	= T - prevT;
	if( dT.norm() <= EPSILON ) {
		return prevQ;
	}
	else {
		prevNormal	= dT - (dT * prevT) * prevT;
		prevNormal.normalize();
		normal		= dT - (dT * T) * T;
		normal.normalize();
		R.buildFromFrame( prevT, prevNormal );
		P.buildFromFrame( T, normal );
		return P * ( R.conj() * prevQ );
	}
}

//
// This function should be called whenever a figure is converted back from symmetric space
// representation.
// It uses the helper function phiAlignAtom defined above.
//
void M3DTubeFigure::fixGlobalConsistency()
{
	int ibaseAtom	= getBaseAtomIndex();
	int starti;
	int endi;
	int stepi;
	Vector3D T;
	M3DPrimitive *prevPrim, *prim;

	for( int side = -1; side <= +1; side += 2 ) {
		if( side < 0 ) {
			starti	= ibaseAtom - 1;
			endi	= -1;
			stepi	= -1;
		}
		else {
			starti	= ibaseAtom+1;
			endi	= getPrimitiveCount();
			stepi	= +1;
		}
		prevPrim	= getPrimitivePtr(ibaseAtom);
		for( int i = starti; i != endi; i += stepi ) {
			prim	= getPrimitivePtr(i);
			T.set(1,0,0);
			prim->getQ().rotateVector(T);
			prim->setQ( phiAlignAtom( prevPrim->getQ(), T ) );
			prevPrim	= prim;
		}
	}
}


int M3DTubeFigure::getPrimitiveCount() const
{
    return primitives.size();
}

// dibyendu

// Function not defined for tubular m-reps yet
// Its simple to write but I dont use tube m-reps; so havent written it yet

int M3DTubeFigure::getSpokeCount() const
{

	cout << "Function to return spoke count for tubular m-reps not yet defined" << endl ;

	return( 1 ) ;

}

M3DPrimitive * M3DTubeFigure::getPrimitivePtr(int primIndex) const
{
    if(primIndex < 0 || primIndex >= primitives.size())
        return NULL;
    return (primitives[primIndex]);
}

void M3DTubeFigure::setPrimitivePtr(int primIndex, M3DPrimitive * newPrimPtr)
{
    if(primIndex < 0 || primIndex > primitives.size())
        return;

    if(primitives[primIndex] != NULL)
        delete primitives[primIndex];

    primitives[primIndex] = newPrimPtr;

    properties->modified = true;
}

void M3DTubeFigure::addColumn(M3DPrimitive * column, int columnIndex)
{
    if(columnIndex < 0 || columnIndex > numColumns)
        return;

    primitives.insert(primitives.begin() + columnIndex, column);
    properties->modified = true;
}

void M3DTubeFigure::select()
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

void M3DTubeFigure::deselect()
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

void M3DTubeFigure::toggleSelected()
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

bool M3DTubeFigure::isSelected()
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

bool M3DTubeFigure::isAnySelected()
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
int M3DTubeFigure::numberSelected()
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

// Return index of any selected column; -1 if anything else is found
int M3DTubeFigure::columnSelected()
{
    M3DPrimitive * primPtr;
	int col;
	bool full;

	full = false;
    for (col = 0; col < numColumns; col++)
    {
        primPtr = getPrimitivePtr(col);
        if (primPtr != NULL && primPtr->isSelected()) {
			return col;
		}
	}
	return -1;
}

// Reverse the indexing along "the" row
void M3DTubeFigure::reverseRows()
{
	M3DPrimitive* temp;
	for( int j = 0, i = primitives.size() - 1; j < primitives.size()/2; ++j, --i ) {
		temp	= primitives[j];
		primitives[j]	= primitives[i];
		primitives[i]	= temp;
	}
}


// Verifies that the sail ends and extended B-vector end lie inside the unit cube
bool M3DTubeFigure::verifyInBounds() const
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

Vector3D M3DTubeFigure::getCOG(bool all) const
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

void M3DTubeFigure::scaleWidth(double scalefact)
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

void M3DTubeFigure::scaleBy(double scalefact)
{
    scaleBy(scalefact, getCOG());
}

void M3DTubeFigure::scaleBy(double scalefact, const Vector3D &vCenter)
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

void M3DTubeFigure::rotateBy(const Quat &q)
{
    rotateBy(q, getCOG());
}

void M3DTubeFigure::rotateBy(const Quat &q, const Vector3D &vCenter)
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

void M3DTubeFigure::translateBy(const Vector3D &vTrans, bool selectAll)
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

void M3DTubeFigure::deleteAll()
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

    numColumns = 0;

	if (boundary != NULL)
		delete boundary;
	boundary = NULL;

    properties->modified = true;
}

// REGULARITY CALCULATION
// This is an accessor for the regularity calculations done in M3DFigure.  In a TubeMesh, there
// are less number of ways of defining how many neighbors a given atom has compared to QuadMesh.
// These are enum'd in M3DFigure's header.
void M3DTubeFigure::getAtomsNeighbors(int j, int &numNeighbors, M3DPrimitive* neighbors[], PrimNeighborhoodDefn PND)
{

	if (PND == CORNERS_HAVE_3_NEIGHBORS)// this is only for quad figure
		PND = ALL_FIRST_NEIGHBORS;

	// FIXME: I don't know what the Pnds mean here - rohit
	numNeighbors	= 0;
	if( j == 0 ) {
		// west corner
		if( PND == ALL_FIRST_NEIGHBORS || PND == EDGES_HAVE_2_NEIGHBORS || PND == PIN_CORNERS ) {
			neighbors[0]	= getPrimitivePtr(j+1);
			numNeighbors	= 1;
		}
	}
	else if( j == numColumns-1 ) {
		// east corner
		if( PND == ALL_FIRST_NEIGHBORS || PND == EDGES_HAVE_2_NEIGHBORS || PND == PIN_CORNERS ) {
			neighbors[0]	= getPrimitivePtr(j-1);
			numNeighbors	= 1;
		}
	}
	else {
		// internal
		numNeighbors	= 2;
		neighbors[0]	= getPrimitivePtr(j+1);
		neighbors[1]	= getPrimitivePtr(j-1);
	}
}


// Persistent object support
void M3DTubeFigure::writeFigure( const char * figureStr, Registry& registry )
{
    int numColumns;
    int i;
	int numLandmarks;
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

    numColumns = getColumnCount();

    figureName = getName();
    if(figureName != NULL)
        registry.setStringValue(format, (char *) figureName, "name");
    registry.setStringValue(format, friendlyFigureName, "type");

    registry.setIntValue(format, numColumns, "numColumns");
	registry.setIntValue(format, 1, "numRows");
	registry.setIntValue(format, numberOfSpokes, "numberOfSpokes");

    registry.setBooleanValue(format, isPositiveSpace(), "positiveSpace");
    registry.setBooleanValue(format, hasPositivePolarity(), "positivePolarity");

    registry.setIntValue(format, getTolerance(), "smoothness");

    color = getColor();
    registry.setDoubleValue(format, color[0], "color.red");
    registry.setDoubleValue(format, color[1], "color.green");
    registry.setDoubleValue(format, color[2], "color.blue");

	numLandmarks = getLandmarkCount();
    registry.setIntValue(format, numLandmarks, "numLandmarks");
#ifdef BINARY
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

    sprintf(subordinateStr, format, "primitive[%d]");
	for(i = 0; i < numColumns; i++)
	{
		getPrimitivePtr(i)->writePrimitive(registry, subordinateStr, i);
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

//	FIXME @see M3dObjectFile::write
//	// Write boundary information
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
}

M3DFigure * M3DTubeFigure::readFigure(int figureIndex, Registry& registry)
{
    M3DTubeFigure * figure;
    const char * name;
    int  numColumns;
	int numberOfSpokes;
    int i;
    float color[3];
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

    numColumns = registry.getIntValue("model.figure[%d].numColumns", 0, figureIndex);
	numberOfSpokes	= registry.getIntValue("model.figure[%d].numberOfSpokes", M3DTubePrimitive::DEFAULT_NUMBER_OF_SPOKES, figureIndex);

    positiveSpace = registry.getBooleanValue("model.figure[%d].positiveSpace",
		true, figureIndex);
    positivePolarity = registry.getBooleanValue("model.figure[%d].positivePolarity",
		true, figureIndex);

    tolerance = registry.getIntValue("model.figure[%d].smoothness",
		M3DFigure::getDefaultSurfaceTolerance(), figureIndex);

    figure = new M3DTubeFigure(numColumns, numberOfSpokes);
    if(figure == NULL)
        return NULL;

    figure->setName(name);
    figure->setColor(color);
    figure->setPositiveSpace(positiveSpace);
    figure->setPositivePolarity(positivePolarity);
    figure->setTolerance(tolerance);

	for(i = 0; i < numColumns; i++)
	{
		figure->setPrimitivePtr(i, M3DTubePrimitive::readPrimitive(registry, numberOfSpokes, "model.figure[%d].primitive[%d]", figureIndex, i));
	}

#ifdef LM_METHOD_OBJ
	// read landmarks after atoms, since landmarks refer to atoms
    const int numLandmarks = registry.getIntValue("model.figure[%d].numLandmarks", 0, figureIndex);
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

            // for tubes we need to be explicit about U,V for a landmark
			// we can not specify a landmark based on atom index alone
			if ((atomU != -1)  || (atomV != -1)) {
			  figure->addLandmark(atomU, atomV, permName, atomT);
            }
			//figure->addLandmark(atomIndex, permName, atomT);

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
            int j = registry.getIntValue("model.figure[%d].constraint[%d].figure",
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

	// FIXME @see read in object file
//	figure->setFigureStatsPtr(readFigureStats("model.figure[%d]",
//	    figureIndex));

    return figure;
}


/**
 * Returns the index of the base atom in this tube.
 * If no atom is marked as the base atom, then
 * the center atom is chosen by default.
 * @param	none
 * @return index of base atom.
 */
int M3DTubeFigure::getBaseAtomIndex() const
{
	int ibaseAtom = 0;
	for(;ibaseAtom != getPrimitiveCount(); ibaseAtom++ ) {
		if((dynamic_cast<const M3DTubePrimitive*>(getPrimitivePtr(ibaseAtom)))->isBaseAtom()) {
			break;
		}
	}
	if( ibaseAtom >= getPrimitiveCount() ) {
		ibaseAtom	= getPrimitiveCount()/2;
	}
	return ibaseAtom;
}

void M3DTubeFigure::atomAtCoordinates( const double* m, M3DPrimitive& prim ) const
{
	assert( *m >= 0.0 && *m < getPrimitiveCount() - 1 );
	const int u		= int(*m);
	const double t	= *m - u;
//	cout << "Interpolating " << u << " and " << (u+1) << ", t = " << t << endl;
	const M3DTubePrimitive* prim0	= dynamic_cast<const M3DTubePrimitive*>(getPrimitivePtr(u));
	const M3DTubePrimitive* prim1	= dynamic_cast<const M3DTubePrimitive*>(getPrimitivePtr(u+1));

	const Vector3D p0	= prim0->getX();
	const Vector3D p1	= prim1->getX();
	const Vector3D d0	= prim0->getB()*(p0-p1).norm();
	const Vector3D d1	= prim1->getB()*(p0-p1).norm();

	// position of the interpolated atom (Hermite spline).
	const Vector3D x	= (t*t * ( 2.0*t-3) + 1.0) * p0
						+ (t*t * (-2.0*t+3)      ) * p1
						+ (t * (t * (t-2.0) + 1.0)) * d0
						+ (t * (t * (t - 1.0) )  ) * d1;
	// interpolate the rest using ordinary geodesic interpolation
	if(!prim.atomInterp( t, prim0, prim1 )) {
		assert(false);
	}
	prim.setX(x);
}

void M3DTubeFigure::subdivide()
{
	int nColsx2	= numColumns * 2 - 1;
	int i;
	std::vector<M3DPrimitive *> newPrimitives(nColsx2);
	for(i = 0; i != numColumns; ++i ) {
		newPrimitives[2*i]	= primitives[i];
	}
	for(i = 0; i != numColumns - 1; ++i ) {
		const double u	= i + 0.5;
		newPrimitives[2*i+1]	= new M3DTubePrimitive(*dynamic_cast<M3DTubePrimitive*>(newPrimitives[2*i]));
		dynamic_cast<M3DTubePrimitive*>(newPrimitives[2*i+1])->setBaseAtom(false);
		getInterpolatedAtom( newPrimitives[2*i+1], &u );
//		atomAtCoordinates(&u,*(newPrimitives[2*i+1]));
	}
	primitives	= newPrimitives;
	numColumns	= nColsx2;

#ifdef LM_METHOD_OBJ
	// Change atomIndex and u for the landmarks.
	for( i = 0; i != landmarkCount; ++i ) {
		landmarkAtomUs[i]		*= 2;
		landmarkAtomIndices[i]	*= 2;
	}
#endif
	fixGlobalConsistency();
}


inline double hermiteLength(const Vector3D& p0, const Vector3D& p1, const Vector3D& d0, const
Vector3D& d1, const double t)
{
	const Vector3D dg	= (-6*t+6*t*t) * p0 + (6*t-6*t*t) * p1
						+ (1-4*t+3*t*t) * d0 + (-2*t+3*t*t)	* d1;
	return dg.norm();
}

double M3DTubeFigure::curveLength(const Vector3D& p0, const Vector3D& p1, const Vector3D& d0, const Vector3D& d1, const double t) 
{
		//
		// Evaluate the total using Simpson's 1/3rd rule.
		//
		const int n	= 20;				// number of points of evaluation = 1 + 2*n
		const double h	= t/(2.0*n);			// spacing between points.
		const double x0	= hermiteLength(p0, p1, d0, d1, 0.0);
		double segmentLength	= x0;
		for( int i = 1; i != n; ++i ) {
			const double x2i	= hermiteLength(p0, p1, d0, d1, 2*i*h);
			const double x2ip1	= hermiteLength(p0, p1, d0, d1, (2*i+1)*h);
			segmentLength += 4 * x2ip1 + 2 * x2i;
		}
		const double x2n	= hermiteLength(p0, p1, d0, d1, t);
		segmentLength	+= x2n;
		segmentLength = segmentLength * h/3.0;

		return segmentLength;


}


void M3DTubeFigure::resampleForRegularSpacing()
{
    int iprim;

	// Get total length of the curve and then
	// do the resampling ...
	const M3DTubePrimitive* prim0	= dynamic_cast<const M3DTubePrimitive*>(getPrimitivePtr(0));
	const M3DTubePrimitive* prim1;
	double curveLength	= 0.0;
	double* cumCurveLengths	= new double[getPrimitiveCount()];
	cumCurveLengths[0]	= 0.0;
	for (iprim = 1; iprim != getPrimitiveCount(); ++iprim) {
		prim1	= dynamic_cast<const M3DTubePrimitive*>(getPrimitivePtr(iprim));

		const Vector3D p0	= prim0->getX();
		const Vector3D p1	= prim1->getX();
		const Vector3D d0	= prim0->getB()*(p0-p1).norm();
		const Vector3D d1	= prim1->getB()*(p0-p1).norm();

		double segmentLength = M3DTubeFigure::curveLength(d0, d1, p0, p1, 1.0);
		curveLength	+= segmentLength;
		cumCurveLengths[iprim]	= curveLength;

		prim0	= prim1;
	}

	cout << "Total curve length = " << curveLength << endl;
	// Reposition all other intermediate atoms...
    std::vector<M3DPrimitive *> primitives(getPrimitiveCount());
	primitives[0]	= getPrimitivePtr(0);
	primitives[getPrimitiveCount()-1]	= getPrimitivePtr(getPrimitiveCount()-1);
	for (iprim = 1; iprim < getPrimitiveCount() - 1; ++iprim) {
		int iu	= getPrimitiveCount() - 1;
		const double s = curveLength/(getPrimitiveCount()-1) * iprim;
		while( cumCurveLengths[iu] > s )
			iu--;
		const double u = (s - cumCurveLengths[iu])/(cumCurveLengths[iu+1] - cumCurveLengths[iu]) + iu;
		cout << "atom " << iprim << " -> " << u << endl;
		M3DTubePrimitive* prim	= new M3DTubePrimitive();
		primitives[iprim]	= prim;
		getInterpolatedAtom(prim, &u);
		prim->setBaseAtom(dynamic_cast<M3DTubePrimitive*>(getPrimitivePtr(iprim))->isBaseAtom());
	}
	for (iprim = 1; iprim < getPrimitiveCount() - 1; ++iprim) {
		delete getPrimitivePtr(iprim);
	}
	this->primitives	= primitives;
	
	delete[] cumCurveLengths;
}

// Interpolates rSrad between 0 and 1.
inline double rSrad(const double u, const double rSrad0, const double rSrad1)
{
	return 1.0 - pow(1.0 - rSrad0, u) * pow(1.0 - rSrad1, 1.0-u);
}

inline double drSrad(const double u, const double rSrad0, const double rSrad1)
{
	return (log(1.0 - rSrad0) - log(1.0 - rSrad1)) *
		pow(1.0 - rSrad0, u) * pow(1.0 - rSrad1, 1.0 - u);
}


Vector3D spokeAtPosition(double u, const Vector3D& S0,
	const Vector3D& p0, const Vector3D& d0,
	const Vector3D& p1, const Vector3D& d1,
	const double rSrad0, const double rSrad1)
{
	// position of the interpolated atom (Hermite spline) (at position t=0).
	const Vector3D x	= p0;
	const Vector3D dx   = d0;
	// FIXME: Hermite splines are only C1 continuous!!
	// So is this correct to do?
	const Vector3D d2x	= -6 * p0 + 6 * p1
						+ (-4)*d0 + (-2)*d1;
	//const Vector3D d3x	= 12 * p0 - 12 * p1 + 6 * d0 + 6 * d1;

	if( rSrad0 >= 1.0 || rSrad1 >= 1.0 ) {
		cout << __FILE__ << ":" << __LINE__ << ": spokeAtPosition(), rSrad out of bound (" << rSrad0 << "," << rSrad1 << ")\n";
	}
	const double phi	= rSrad(0, rSrad0, rSrad1);
	const double dphi	= drSrad(0, rSrad0, rSrad1);
	// See mathematica notebook for details of the expression.
	const double S0S0  = S0 * S0;
	const double dxS0  = dx * S0;
	const double a      = (phi-1)/S0S0;
	const Vector3D dS   = ((dxS0)*a) * S0 - phi*dx;
	const Vector3D d2S  = ((((dx *dS) + (d2x*S0))*a)
		- (dxS0*(2*(dS*S0))*a/S0S0)) * S0
		+ (dxS0*a) * dS
		+ (dxS0 / S0S0 * dphi) * S0
		- phi * d2x - dphi * dx;
	// Taylor series expansion of Su (just to see it in
	// the debugger)
	Vector3D Su = S0;
	Su += u * dS;
	Su += u*u/2 * d2S;
	return Su;
}

extern
double rSrad( const double phi, const double dt, 
	const M3DTubePrimitive* prev, const M3DTubePrimitive* prim, const M3DTubePrimitive* next);

bool
M3DTubeFigure::getInterpolatedAtom( M3DPrimitive* _prim, const double* _u ) const
{
	int iphi;

	M3DTubePrimitive* primu	= dynamic_cast<M3DTubePrimitive*>(_prim);
	double u	= *_u;
	if( u < 0 && u > getPrimitiveCount() - 1 ) {
		return false;
	}
	const int uprim0	= int(u);
	u	= u - uprim0;
	const M3DTubePrimitive* prim0 = dynamic_cast<const M3DTubePrimitive*>(getPrimitivePtr(uprim0));
	if( u <= 1e-5 ) {
		// very close to first primitive? Discard M3DEndPrimitive type qualifiers and go ahead.
		*primu	= M3DTubePrimitive(*prim0);
		return true;
	}
	const M3DTubePrimitive* prim1 = dynamic_cast<const M3DTubePrimitive*>(getPrimitivePtr(uprim0+1));
	if( u >= 1-1e-5) {
		// very close to second primitive? Discard M3DEndPrimitive type qualifiers and go ahead.
		*primu	= M3DTubePrimitive(*prim1);
		return true;
	}

	*primu	= M3DTubePrimitive(*prim0);
	// First find interpolated hub position.
	const Vector3D p0	= prim0->getX();
	const Vector3D p1	= prim1->getX();
	const Vector3D d0	= prim0->getB()*(p0-p1).norm();
	const Vector3D d1	= prim1->getB()*(p0-p1).norm();

	const Vector3D x	= (u*u * ( 2.0*u-3) + 1.0) * p0
						+ (u*u * (-2.0*u+3)      ) * p1
						+ (u * (u * (u-2.0) + 1.0)) * d0
						+ (u * (u * (u - 1.0) )  ) * d1;
	primu->setX(x);

	// Get the previous and next incremental primitives by
	// geodesic interpolation. We will use these for the estimates.
	// (Ideally we should recurse with the incremental primitives
	// initialized by the estimates by the geodesic.)
	const double dt = 0.02;
	double t;
	// Initialize with the nearest primitives as when the conversion takes place from sym space to
	// lie space, we want a phi orientation that is close to
	// what we wish.
	M3DTubePrimitive prev0(*prim0), prev1(*prim1), next0(*prim0), next1(*prim1);
	if( uprim0 > 0 ) {
		// prim0 is not the 1st primitive
		t	= uprim0 - dt;
		atomAtCoordinates( &t, prev0 );
	}
	// This should be true.
	//if( uprim0 < getPrimitiveCount()-1 ) {
		t	= uprim0 + dt;
		atomAtCoordinates( &t, next0 );
	//}
	// This should be true.
	//if( uprim0+1 > 0 ) {
		t	= uprim0+1 - dt;
		atomAtCoordinates( &t, prev1 );
	//}
	if( uprim0+1 < getPrimitiveCount()-1 ) {
		// prim1 is not the last primitive
		t	= uprim0+1 + dt;
		atomAtCoordinates( &t, next1 );
	}
	Vector3D* spokes	= new Vector3D[getNumberOfSpokes()];
	// Now interpolate the spokes
	for (iphi = 0; iphi != getNumberOfSpokes(); ++iphi) {
		const double phi	= 2.0*R_PI*iphi/getNumberOfSpokes();
		// compute rSrads
		const double rSrad0 = rSrad(phi, dt,
			(uprim0 > 0) ? &prev0 : NULL,
			prim0, &next0);
		const double rSrad1 = rSrad(phi, dt,
			&prev1, prim1,
			(uprim0+1 < getPrimitiveCount() - 1) ? &next1 : NULL);
		// average the two results to create the spoke.
		// (This is to make the construction more robust)
		Vector3D S[2];
		// Use a continous weighting, weighting the result that is closer to a primitive more.
		// Why this? no good reason, just need something with a smooth ramp in the central region.
		double weights[] = { max(0.6-u,0.0)*max(0.6-u,0.0), max(u-0.4,0.0)*max(u-0.4,0.0) };
		double tmp	= weights[0] + weights[1];
		weights[0]	/= tmp;
		weights[1]	/= tmp;
		// Interpolate spokes in both directions.
		S[0] = spokeAtPosition(u, prim0->getYPhi(phi),
			p0, d0, p1, d1, rSrad0, rSrad1 );
		S[1] = spokeAtPosition(1-u, prim1->getYPhi(phi),
			p1, -d1, p0, -d0, rSrad1, rSrad0 );
		const double r = exp( weights[0]*log(S[0].normalize()) + weights[1]*log(S[1].normalize()) );
		spokes[iphi] = r * ShapeSpace::S2::mean( S, 2, weights );
	}

	/*
	// Set the bisector vector as the tangent at u
	// and phi align the atom with the previous one.
	Vector3D dx	= (-6*u+6*u*u) * p0 + (6*u-6*u*u) * p1
					+ (1-4*u+3*u*u) * d0 + (-2*u+3*u*u)	* d1;
	dx.normalize();
	primu->setQ( phiAlignAtom(prim0->getQ(), dx) );
	*/
	// Set the bisector vector to the average of the two bisector vectors. This gives more
	// stability ... and just happens to be mathematically a fully bogus thing ...
	// anyway this function is completely mathematically bogus, and the representation is bogus
	// too ... oh well, whatever ...
	const Vector3D bisectors[]	= { prim0->getB(), prim1->getB() };
	const Vector3D dx = ShapeSpace::S2::mean( bisectors, 2, NULL );
	primu->setQ( phiAlignAtom(prim0->getQ(), dx) );

	// Finally set the cone angle and r as the "average" of all
	// the spokes and adjust their deviations accordingly.
	// Get the mean projection length and the angle the spokes make
	// with the bisector.
	double hca	= 0.0;
	for (iphi = 0; iphi != getNumberOfSpokes(); ++iphi) {
		const double coneAngle	= acos( spokes[iphi] * primu->getB() / spokes[iphi].norm());
//		cout << "coneAngle[" << iphi << "] = " << coneAngle*180/R_PI;
		hca += ShapeSpace::RP1::Log( coneAngle );
	}
	hca	= ShapeSpace::RP1::Exp(hca/getNumberOfSpokes());
	double r	= 0.0;
	primu->setTheta(hca);
	// The radius of the cone is the radius obtained by projecting all the spokes onto this cone.
	// However, if the cone is close to a plane, then we can have numerical issues, in which case,
	// we just take the mean of the radii. (0.1 ~ 5 degrees off 90)
	// FIXME: Just use the first one, the negative entries come out all too often and we don't want
	// a sudden change in the method employed in between. Oh yeah, by doing this, we no longer know
	// whether spokes will cross or not.
	//if( fabs(cos(hca)) < 0.1 ) {
	for (iphi = 0; iphi != getNumberOfSpokes(); ++iphi) {
		r += log( spokes[iphi].norm() );
	}
	//}
	//else {
	//	for( int iphi = 0; iphi != getNumberOfSpokes(); ++iphi ) {
	//		const double radius	= spokes[iphi] * primu->getB() / cos(hca);
	//		if( radius < 0.0 ) {
	//			// Ouch!
	//			cout << __FILE__ << ":" << __LINE__ << " getInterpolatedAtom(...), radius of base cone has negative entry = " << radius << endl;
	//		}
	//		r += log( radius );
	//	}
	//}
	
	r	= exp( r / getNumberOfSpokes() );
	primu->setR(r);

	// Now set each of the individual spokes ...
	for (iphi = 0; iphi != getNumberOfSpokes(); ++iphi) {
		primu->setRN( iphi, spokes[iphi].norm() );
	}

	delete[] spokes;
	return true;
}


//
// This is the norm used for averaging the curviness values ...
// A higher value will strongly penalize local aberrations.
//
const double pnorm = 10.0;

/**
 * A helper function that returns the total curvature at t for the section of the medial curve
 * between prim0 (t=0) and prim1 (t=1).
 */
inline double totalCurvatureSquared( const M3DTubePrimitive* prim0, const M3DTubePrimitive* prim1, const double t)
{
	const Vector3D p0	= prim0->getX();
	const Vector3D p1	= prim1->getX();
	const Vector3D d0	= prim0->getB()*(p0-p1).norm();
	const Vector3D d1	= prim1->getB()*(p0-p1).norm();

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
//	const double tau	= (fabs(normdgxd2g) <= EPSILON) ?
//		0.0 : dg * (d2g.cross(d3g)) / (normdgxd2g * normdgxd2g) ;
	const double tau	= 0.0;

	return pow(fabs(kappa),pnorm) + pow(fabs(tau),pnorm);
}


/**
 * This function returns the curviness of the medial sheet.
 *
 * @param	none
 * @return	the curviness of the medial sheet.
 */
double M3DTubeFigure::curviness() const
{
	const M3DTubePrimitive* prim0	= dynamic_cast<const M3DTubePrimitive*>(getPrimitivePtr(0));
	const M3DTubePrimitive* prim1;
	double totalCurviness	= 0.0;
	double approxCurveLength	= 0.0;
	//double meanRadius	= log(prim0->getR());
	for( int iprim = 1; iprim != getPrimitiveCount(); ++iprim) {
		prim1	= dynamic_cast<const M3DTubePrimitive*>(getPrimitivePtr(iprim));
		approxCurveLength	+= (prim1->getX() - prim0->getX()).norm();
		//meanRadius	+= log(prim1->getR());
		//
		// Evaluate the total using Simpson's 1/3rd rule.
		//
		const int n	= 20;				// number of points of evaluation = 1 + 2*n
		const double h	= 1.0/(2.0*n);			// spacing between points.
		const double x0	= totalCurvatureSquared(prim0,prim1,0.0);
		double segmentCurviness	= x0;
		for( int i = 1; i != n; ++i ) {
			const double x2i	= totalCurvatureSquared(prim0,prim1, 2*i*h);
			const double x2ip1	= totalCurvatureSquared(prim0,prim1, (2*i+1)*h);
			segmentCurviness += 4 * x2ip1 + 2 * x2i;
		}
		const double x2n	= totalCurvatureSquared(prim0,prim1,1.0);
		segmentCurviness	+= x2n;
		segmentCurviness = segmentCurviness * h/3.0;
		totalCurviness	+= segmentCurviness;

		prim0	= prim1;
	}
	//meanRadius	/= getPrimitiveCount();
	//meanRadius	= exp(meanRadius);
	//
	// Scale by the mean radius / inter atom spacing.
	//
//	cout << "t: " << totalCurviness << ", mr: " << meanRadius << endl;
	//totalCurviness	= meanRadius * pow(totalCurviness, 1.0/pnorm);
	totalCurviness	= approxCurveLength * pow(totalCurviness, 1.0/pnorm);
	totalCurviness	/= getPrimitiveCount();
	totalCurviness	/= getPrimitiveCount();
	return totalCurviness * totalCurviness;
}



void M3DTubeFigure::rotateAlongAxisBy( double phi )
{
	// Get base atom from the tube, and do the phi alignment.
	M3DPrimitive* prim	= getPrimitivePtr( getBaseAtomIndex() );
	Quat q;
	const Vector3D newN	= prim->getN() * cos(phi) + prim->getBPerp() * sin(phi);
	q.buildFromFrame( prim->getB(), newN );
	/* FIXME: Debug stuff
	cerr << "Old B: ";
	prim->getB().print(cerr);
	cerr << "Old q: ";
	prim->getQ().print(cerr);
	cerr << "q phi value :" << phi << endl;
	cerr << "new q: ";
	q.print(cerr);
	*/
	prim->setQ(q);
	//prim->getB().print(cerr);
	//cerr << flush;
	fixGlobalConsistency();
}

