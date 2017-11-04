#include <math.h>
#include <stdarg.h>
#include <assert.h>
#include <limits>

#include "M3DPrimitive.h"
#include "VectorND.h"
#include "M3DPGAStats.h"

#include <typeinfo>

#ifndef EPSILON
#define EPSILON 1e-10
#endif

inline
int min(const int& a, const int& b) {
	return (a < b) ? a : b;
}

using namespace std;

const int M3DTubePrimitive::DEFAULT_NUMBER_OF_SPOKES	= 8;

M3DPrimitive & M3DTubePrimitive::operator = (const M3DPrimitive & unknown_prim)
{
    assert( typeid(*this) == typeid(unknown_prim) );
    const M3DTubePrimitive& prim = 
        dynamic_cast<const M3DTubePrimitive&>(unknown_prim);

#ifdef DEBUG
    cout << "M3DTubePrimitive::operator=()" << endl;
#endif

    M3DPrimitive::copy(static_cast<const M3DPrimitive&>(prim));
    copy(prim);

    return (*this);
}

void M3DTubePrimitive::copy(const M3DTubePrimitive & prim)
{
#ifdef DEBUG
    cout << "M3DTubePrimitive::copy() called\n";
#endif
    baseAtom    = prim.baseAtom;
    sr          = prim.sr;
    U2          = prim.U2;
    U0          = prim.U0;
    d           = prim.d;
}

M3DTubePrimitive * M3DTubePrimitive::readPrimitive(Registry& registry, int numberOfSpokes, const char * regStr, ...)
{
    char newStr[1024];
    const char * typeStr;
	int i;
	double* sr, rEnd, d;
	sr = new double[numberOfSpokes];
    double x, y, z;
	double ux[3], uy[3], uz[3];
    bool baseAtom;
	short selected;
	Vector3D U0, U2;

    va_list val;
    va_start(val, regStr);
    vsprintf(newStr, regStr, val);
    va_end(val);

    strcat(newStr, ".%s");

    M3DTubePrimitive * primPtr;

    typeStr = registry.getStringValue(newStr, NULL, "type");
    x =  registry.getDoubleValue(newStr, 0.0, "x");
    y =  registry.getDoubleValue(newStr, 0.0, "y");
    z =  registry.getDoubleValue(newStr, 0.0, "z");
    d =  registry.getDoubleValue(newStr, std::numeric_limits<double>::infinity(), "d");
    selected = registry.getIntValue(newStr, 0, "selected");
	baseAtom = registry.getBooleanValue(newStr, false, "baseAtom");

	if( d == std::numeric_limits<double>::infinity() ) {
		// d was not found in the model file ...
		// old format, compute values for sr, d, ux, uy and uz.
		Quat q;
		double theta, r, elongation;
		r =  registry.getDoubleValue(newStr, 1.0, "r");
		q.setX(registry.getDoubleValue(newStr, 0.0, "qx"));
		q.setY(registry.getDoubleValue(newStr, 0.0, "qy"));
		q.setZ(registry.getDoubleValue(newStr, 0.0, "qz"));
		q.setW(registry.getDoubleValue(newStr, 1.0, "qw"));
		theta = registry.getDoubleValue(newStr, 90.0, "theta") * R_DEGREES_TO_RADIANS;
	    if(typeStr != NULL && strcmp(typeStr, M3D_END_PRIMITIVE_STR) == 0) {
    	    elongation = registry.getDoubleValue(newStr, 1.0, "elongation");
			// compute rEnd
			rEnd = r * elongation;
		}
		
		// compute d
		d = r * cos(theta);
		// compute U0
		U0.set(0, 1, 0);
		q.rotateVector(U0);
		// compute U2
		U2.set(1, 0, 0);
		q.rotateVector(U2);

		strcat(newStr, "[%d]");
		for( i = 0; i < numberOfSpokes; ++i ) {
			sr[i] = registry.getDoubleValue(newStr, 0.0, "dr", i);
			// convert from logarithmic deviation to actual radius value.
			sr[i] = r * exp(sr[i]);
		}
	}
	else {
		strcat(newStr, "[%d]");
		sr[0] = 0.5;
		// We resort to the first spoke length, when additional data
		// is not found.
		for( i = 0; i < numberOfSpokes; ++i ) {
			sr[i] = registry.getDoubleValue(newStr, sr[0], "r", i);
		}
		for( i = 0; i < 3; ++i ) {
			ux[i] = registry.getDoubleValue(newStr, 1.0, "ux", i);
			uy[i] = registry.getDoubleValue(newStr, 1.0, "uy", i);
			uz[i] = registry.getDoubleValue(newStr, 1.0, "uz", i);
		}
    	if(typeStr != NULL && strcmp(typeStr, M3D_END_PRIMITIVE_STR) == 0) {
			rEnd = registry.getDoubleValue(newStr, sr[0], "r", numberOfSpokes);
		}
		U0.set(ux[0], uy[0], uz[0]);
		U2.set(ux[2], uy[2], uz[2]);
	}

    if(typeStr != NULL && strcmp(typeStr, M3D_END_PRIMITIVE_STR) == 0) {
        primPtr = new M3DTubeEndPrimitive(x, y, z, rEnd, d, sr,
			U0, U2, numberOfSpokes);
    }
    else {
        primPtr = new M3DTubePrimitive(x, y, z, d, sr,
			U0, U2, numberOfSpokes);
	}
	primPtr->setBaseAtom(baseAtom);

    if(primPtr == NULL)
        return NULL;

    if(selected & 0x1)
        primPtr->select();
    if(selected & 0x2)
        primPtr->deselectForRegularity();

	delete[] sr;
    return primPtr;
}

void M3DTubePrimitive::writePrimitive(Registry& registry, const char * regStr, ...) const
{
    char newStr[1024];

    va_list val;
    va_start(val, regStr);
    vsprintf(newStr, regStr, val);
    va_end(val);

    strcat(newStr, ".%s");

    registry.setBooleanValue(newStr, isBaseAtom(), "baseAtom");
    registry.setIntValue(newStr, getSelectionFlags(), "selected");

    registry.setDoubleValue(newStr, getX().getX(), "x");
    registry.setDoubleValue(newStr, getX().getY(), "y");
    registry.setDoubleValue(newStr, getX().getZ(), "z");
    registry.setDoubleValue(newStr, d, "d");

	const Vector3D U1 = U2.cross(U0);
    registry.setDoubleValue(newStr, U0.getX(), "ux[0]");
    registry.setDoubleValue(newStr, U0.getY(), "uy[0]");
    registry.setDoubleValue(newStr, U0.getZ(), "uz[0]");
    registry.setDoubleValue(newStr, U1.getX(), "ux[1]");
    registry.setDoubleValue(newStr, U1.getY(), "uy[1]");
    registry.setDoubleValue(newStr, U1.getZ(), "uz[1]");
    registry.setDoubleValue(newStr, U2.getX(), "ux[2]");
    registry.setDoubleValue(newStr, U2.getY(), "uy[2]");
    registry.setDoubleValue(newStr, U2.getZ(), "uz[2]");

    if(type() == M3D_STANDARD_PRIMITIVE) {
        registry.setStringValue(newStr, (char*)M3D_STANDARD_PRIMITIVE_STR, "type");
		strcat(newStr, "[%d]");
	}
    else if(type() == M3D_END_PRIMITIVE) {
        registry.setStringValue(newStr, (char*)M3D_END_PRIMITIVE_STR, "type");
		strcat(newStr, "[%d]");
		registry.setDoubleValue(newStr,
			dynamic_cast<const M3DEndPrimitive*>(this)->getREnd(), "r", sr.size());
	}

	for( int i = 0; i < sr.size(); ++i ) {
	    registry.setDoubleValue(newStr, sr[i], "r", i);
	}
}


void M3DTubePrimitive::setR(double r)
{
	double mean = 0.0;
    int i = 0;  // Needed out here for old compilers
	// Compute the 'radius' in the cone base plane and compute the geometric mean.
	for( i = 0; i < sr.size(); ++i ) {
		mean += log(sqrt( sr[i]*sr[i] - d*d ));
	}
	mean /= sr.size();
	mean = exp(mean);

	const double ratio = mean/r;
	// Now multiply all the in-plane radii with this ratio.
    for( i = 0; i < sr.size(); ++i ) {
		const double inPlaneR = sqrt(sr[i]*sr[i] - d*d) * ratio;
		sr[i] = sqrt(d*d + inPlaneR*inPlaneR);
	}
}

void M3DTubePrimitive::setQ(const Quat & Q)
{
	U2 = Vector3D(1,0,0);
	U0 = Vector3D(0,1,0);
	Q.rotateVector(U2);
	Q.rotateVector(U0);
}

void M3DTubePrimitive::setTheta(double theta)
{
	double mean = 0.0;
	// Compute the angle of each spoke.
	int i;
	for( i = 0; i < sr.size(); ++i ) {
		mean += ShapeSpace::RP1::Log(acos(d/sr[i]));
	}
	mean /= sr.size();
	theta = ShapeSpace::RP1::Log(theta);
	const double delta = theta - mean;
	for( i = 0; i < sr.size(); ++i ) {
		sr[i] = d / cos(ShapeSpace::RP1::Exp(ShapeSpace::RP1::Log(acos(d/sr[i])) + delta));
	}
}

// Calculate a common radius for all spokes
double M3DTubePrimitive::getR() const
{
	double mean = 0.0;
	// Compute the 'radius' in the cone base plane and compute the geometric mean.
	for( int i = 0; i < sr.size(); ++i ) {
		mean += log(sqrt( sr[i]*sr[i] - d*d ));
	}
	mean /= sr.size();
	mean = exp(mean);
	// Now convert this mean in-plane 'radius' to the one for the inclined cone spokes.
	mean = sqrt(mean*mean + d*d);
	return mean;
}

// Return a quaternion describing the orientation
Quat M3DTubePrimitive::getQ() const
{
	return Quat( U2, U0, U2.cross(U0) );
}

// Return the angle between a representative spoke
// and the bisector
double M3DTubePrimitive::getTheta() const
{
	double mean = 0.0;
	// Compute the angle of each spoke.
	for( int i = 0; i < sr.size(); ++i ) {
		mean += ShapeSpace::RP1::Log(acos(d/sr[i]));
	}
	mean /= sr.size();
	mean = ShapeSpace::RP1::Exp(mean);
	return mean;
}


// Primitive transformation functions.
// Scaling and rotation are about medial site location.
void M3DTubePrimitive::rotateBy(Quat dQ)
{
	dQ.rotateVector(U2);
	dQ.rotateVector(U0);
}

void M3DTubePrimitive::scaleBy(double mag)
{
	for(int i = 0; i < sr.size(); ++i ) {
		sr[i] *= mag;
	}
	d *= mag;
}


bool M3DTubePrimitive::atomInterp( double t, const M3DPrimitive* m1, const M3DPrimitive* m2) {
	if( m1->type() == M3D_END_PRIMITIVE && m2->type() == M3D_END_PRIMITIVE && type() != M3D_END_PRIMITIVE) {
		if( globalVerbosity >= 0 ) {
			cout << "Warning: Interpolating two end primitives and saving result in an ordinary primitive." << endl;
		}
	}
	TubeSymPrimitive sM1(dynamic_cast<const M3DTubePrimitive*>(m1));
	TubeSymPrimitive sM2(dynamic_cast<const M3DTubePrimitive*>(m2));
	if( sM1.sr.size() != sM2.sr.size() ) {
		if( globalVerbosity >= 0 ) {
			cout << "M3DTubePrimitive::atomInterp: Interpolating two primitives with different number of spokes is not supported." << endl;
		}
	}
	TubeSymPrimitive deltaSM, interpolatedSM;
	deltaSM.sr.resize(min(sM1.sr.size(),sM2.sr.size()));
	interpolatedSM.sr.resize(min(sM1.sr.size(),sM2.sr.size()));
	Vector3D rotatedSpoke;
	Quat q;
	int i;

	deltaSM.x = (sM2.x-sM1.x)*t;
	for(i = 0; i < deltaSM.sr.size(); ++i ) {
		deltaSM.sr[i] = (log(sM2.sr[i] - fabs(sM2.d)) - log(sM1.sr[i] - fabs(sM1.d)))*t;
	}
	deltaSM.rEnd = (log(sM2.rEnd) - log(sM1.rEnd))*t;
	deltaSM.d = (sM2.d - sM1.d)*t;

	interpolatedSM.x = sM1.x+deltaSM.x;
	interpolatedSM.rEnd = exp( log(sM1.rEnd) + deltaSM.rEnd);
	interpolatedSM.d = sM1.d + deltaSM.d;
	for(i = 0; i < deltaSM.sr.size(); ++i ) {
		interpolatedSM.sr[i] = fabs(interpolatedSM.d) + exp( log(sM1.sr[i] - fabs(sM1.d)) + deltaSM.sr[i]);
	}

	q=ShapeSpace::S2::rotationToOrigin(sM1.U2);
	q=q.conj();
	rotatedSpoke = sM2.U2;
	q.conj().rotateVector(rotatedSpoke);
	Vector2D logDatat =ShapeSpace::S2::Log(rotatedSpoke);
	logDatat *= t;
	interpolatedSM.U2 = ShapeSpace::S2::Exp(logDatat);
	q.rotateVector(interpolatedSM.U2);

	interpolatedSM.convert2Lie(this);

	if(m1->isSelected())
		select();
	else
		deselect();
	if(m1->isSelectedForRegularity())
		selectForRegularity();
	else
		deselectForRegularity();
	toggleHinge(m1->isHinge());
	return true;
}


double M3DTubePrimitive::atomSquareDistance( const M3DPrimitive* m, const double* radius ) const {
	double effectiveR = (radius == NULL)
		? effectiveR = sqrt(m->getR()*getR())
		: *radius;

	TubeSymPrimitive sM1(this);
	TubeSymPrimitive sM2(dynamic_cast<const M3DTubePrimitive*>(m));
	TubeSymPrimitive deltaSM;
	int i;
	if( sM1.sr.size() != sM2.sr.size() ) {
		if( globalVerbosity >= 0 ) {
			cout << "M3DTubePrimitive::atomSquareDistance: Taking differences between two primitives with different number of spokes is not supported." << endl;
		}
	}
	deltaSM.sr.resize(min(sM1.sr.size(),sM2.sr.size()));

	deltaSM.x = (sM2.x-sM1.x);
	for(i = 0; i < deltaSM.sr.size(); ++i ) {
		// We divide by the sqrt of the size so that the distance is not influenced by the amount of
		// discretization we choose for the cone.
		deltaSM.sr[i] = (log(sM2.sr[i] - fabs(sM2.d)) - log(sM1.sr[i] - fabs(sM1.d)))*effectiveR /
			sqrt(double(deltaSM.sr.size()));
	}
	deltaSM.rEnd = (log(sM2.rEnd) - log(sM1.rEnd))*effectiveR;
	deltaSM.d = (sM2.d - sM1.d);

	Quat q=ShapeSpace::S2::rotationToOrigin(sM1.U2);
	q.normalize();
//	q=q.conj();
//	(q.conj()).rotateVector(sM2.U2);
	q.rotateVector(sM2.U2);
	Vector2D logDatat = ShapeSpace::S2::Log(sM2.U2)*effectiveR;
//	logDatat = sphereLogMap(sM2.U2)*sqrt(m->getR()*getR());


	double sumSquare = 0;
	sumSquare  = deltaSM.x.norm()*deltaSM.x.norm();
	sumSquare += deltaSM.d*deltaSM.d;
	sumSquare += deltaSM.rEnd*deltaSM.rEnd;
	sumSquare += logDatat.norm()*logDatat.norm();
	for( i = 0; i != deltaSM.sr.size(); ++i ) {
		sumSquare	+= deltaSM.sr[i] * deltaSM.sr[i];
	}

	return sumSquare;
}


bool M3DTubePrimitive::atomAverage( const int mNum, M3DPrimitive** mList,
	const double* _weights)
{
	assert( mNum > 0 && mList != NULL );

	double* weights	= new double[mNum];
	int i,j;

	if( _weights == NULL ) {
		for( i = 0; i != mNum; ++i ) {
			weights[i]	= 1.0/mNum;
		}
	}
	else {
		// Normalize the weights.
		double totalWeight	= 0.0;
		for( i = 0; i != mNum; ++i ) {
			totalWeight	+= _weights[i];
		}
		for( i = 0; i != mNum; ++i ) {
			weights[i]	= _weights[i]/totalWeight;
		}
	}

	Vector3D* spokes	= new Vector3D[mNum];
	TubeSymPrimitive mean;
	mean.x.set(0.0,0.0,0.0);
	mean.rEnd	= 0.0;
	mean.d		= 0.0;
	mean.sr.resize( dynamic_cast<const M3DTubePrimitive*>(mList[0])->getNumberOfSpokes(), 0.0 );
	for( i = 0; i != mNum; ++i ) {
		const M3DTubePrimitive* tPrim	= dynamic_cast<const M3DTubePrimitive*>(mList[i]);
		TubeSymPrimitive sym	= TubeSymPrimitive(tPrim);
		if( sym.sr.size() != mean.sr.size() ) {
			if( globalVerbosity >= 0 ) {
				cout << "M3DTubePrimitive::atomAverage: Computing a mean of primitives with different number of spokes is not supported." << endl;
			}
		}
		sym.Log();
		mean.x		+= sym.x * weights[i];
		mean.rEnd	+= sym.rEnd * weights[i];
		mean.d		+= sym.d * weights[i];
		for( j = 0; j != min(mean.sr.size(), sym.sr.size()); ++j ) {
			mean.sr[j]	+= sym.sr[j] * weights[i];
		}

		spokes[i]	= tPrim->U2;
		spokes[i].normalize();
	}
	mean.Exp();

	// Compute average of spoke U2.
	mean.U2	= ShapeSpace::S2::mean(spokes, mNum, weights);

	mean.convert2Lie(this);

	// Clean up.
	delete[] weights;
	delete[] spokes;

	return true;
}



/*
 * The logic for this code is copied from GeodesicSym::atomToSymAtom.
 */
M3DTubePrimitive::TubeSymPrimitive::TubeSymPrimitive( const M3DTubePrimitive* prim ) :
	ignorePrimType(false) {
	x	= prim->x;
	d	= prim->d;
	U2	= prim->U2;
	U2.normalize();
	sr	= prim->sr;

	if(!ignorePrimType && prim->type()==M3D_END_PRIMITIVE) {
		rEnd = (dynamic_cast<const M3DTubeEndPrimitive *>(prim))->getREnd();
	}
	else {
		// The particular value here shouldn't matter as long as it is > 0.
		rEnd = 1;
	}
}


/*
 * The logic for this code is copied from GeodesicSym::symAtomToAtom.
 */
bool M3DTubePrimitive::TubeSymPrimitive::convert2Lie( M3DPrimitive* _m ) const {
	M3DTubePrimitive* m	= dynamic_cast<M3DTubePrimitive*>(_m);

	// Use the normal from the previous spoke directions
	// so that the rotation about the tangent is minimized.

	Vector3D n	= m->U0;

	// We really don't expect the old normal to become
	// parallel to the new tangent direction, but still ...
	// Get the new bi-normal direction.
	Vector3D b = U2.cross(n);
	if(b.norm() <= EPSILON) {
		std::cout << __FILE__ << ":" << __LINE__
			<< ": warning: new tangent and normal are parallel" << std::endl
			<< "will try some standard directions for bi-normal." << std::endl;
		// TODO: Use some other 'sensible' direction for computing the bi-normal
		b	= U2.cross( Vector3D(1.0,0.0,0.0) );
		if( b.norm() <= EPSILON ) {
			b	= U2.cross( Vector3D(0.0,1.0,0.0));
		}
	}
	b.normalize();
	// Now compute the new normal.
	n	= b.cross(U2);
	n.normalize();

	m->U2	= U2;
	m->U0	= n;

	m->setX(x);
	m->sr	= sr;
	m->d	= d;

	if(!ignorePrimType && m->type() == M3D_END_PRIMITIVE) {
		(dynamic_cast<M3DTubeEndPrimitive *>(m))->setREnd( rEnd );
	}
	
	return true;
}


/*
 * The logic for this code is based on GeodesicSym::symAtomExp()
 */
void M3DTubePrimitive::TubeSymPrimitive::Exp() {
	// x stays the same.
	// d stays the same.
	rEnd		= exp(rEnd);
	U2			= ShapeSpace::S2::Exp(tU2);
	for(int i = 0; i < sr.size(); ++i ) {
		sr[i]	= exp(sr[i]);
	}
}

/*
 * The logic for this code is based on GeodesicSym::symAtomLogMap()
 */
void M3DTubePrimitive::TubeSymPrimitive::Log() {
	// x stays the same.
	// d stays the same.
	rEnd		= log(rEnd);
	tU2			= ShapeSpace::S2::Log(U2);
	for(int i = 0; i < sr.size(); ++i ) {
		sr[i]	= log(sr[i]);
	}
}


/*
 * The logic for this code is based on GeodesicSym::symAtomDiff()
 */
SymPrimitive& M3DTubePrimitive::TubeSymPrimitive::operator-=(const SymPrimitive& _b) {
	const TubeSymPrimitive& b	= dynamic_cast<const TubeSymPrimitive&>(_b);

	if( sr.size() != b.sr.size() ) {
		if( globalVerbosity >= 0 ) {
			cout << "M3DTubePrimitive::TubeSymPrimitive::operator -=: Cannot compute differences of two primitives with different number of spokes." << endl;
		}
	}

	x		= x - b.x;
	d		= d - b.d;
	rEnd	= rEnd/b.rEnd;
	for( int i = 0; i != min(sr.size(),b.sr.size()); ++i ) {
		sr[i]	= sr[i]/b.sr[i];
	}

	Vector3D U2(b.U2);

	Quat q=ShapeSpace::S2::rotationToOrigin(U2);
	q.normalize();
	q.rotateVector(this->U2);
	this->U2.normalize();

	return *this;
}


SymPrimitive& M3DTubePrimitive::TubeSymPrimitive::operator+=(const SymPrimitive& _b) {
	const TubeSymPrimitive& b	= dynamic_cast<const TubeSymPrimitive&>(_b);

	if( sr.size() != b.sr.size() ) {
		if( globalVerbosity >= 0 ) {
			cout << "M3DTubePrimitive::TubeSymPrimitive::operator +=: Cannot compose two primitives with different number of spokes." << endl;
		}
	}

	x		= x + b.x;
	d		= d + b.d;
	rEnd	= rEnd*b.rEnd;
	for( int i = 0; i != min(sr.size(),b.sr.size()); ++i ) {
		sr[i]	= sr[i]*b.sr[i];
	}

	Vector3D U2(b.U2);

	Quat q=ShapeSpace::S2::rotationFromOrigin(U2);
	q.normalize();
	q.rotateVector(this->U2);
	this->U2.normalize();

	return *this;
}


bool M3DTubePrimitive::TubeSymPrimitive::projectAtom2Vec( Vector & v )
{
	int i = 0;
	Log();

	v		= Vector(3 + sr.size() + 1 + 1 + 2);
	v(i++)	= x.getX();
	v(i++)	= x.getY();
	v(i++)	= x.getZ();
	while(i-3 < sr.size()) {
		v(i++)	= sr[i-3];
	}
	v(i++)	= rEnd;
	v(i++)	= d;
	v(i++)	= tU2.getX();
	v(i)	= tU2.getY();
	return true;
}

bool M3DTubePrimitive::TubeSymPrimitive::vec2Atom( const Vector & v )
{
	int i;

	x	= Vector3D(v(0), v(1), v(2));
	const int nSpokes	= v.size() - 7;
	for(i = 3; i-3 < nSpokes; ++i ) {
		sr[i-3]	= v(i);
	}
	rEnd	= v(i++);
	d		= v(i++);
	tU2		= Vector2D( v(i++), v(i) );

	Exp();
	return true;
}
