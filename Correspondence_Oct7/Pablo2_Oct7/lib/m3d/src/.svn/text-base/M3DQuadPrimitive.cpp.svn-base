#include <math.h>
#include <stdarg.h>
#include <assert.h>
#include "M3DPrimitive.h"
#include "Vector2D.h"
#include "M3DPGAStats.h"

#include <typeinfo>

const double epsilon = 1.0 / 1024.0;
const double sqeps = epsilon * epsilon;

using namespace std;

M3DQuadPrimitive::
M3DQuadPrimitive(const Vector3D& X, 
                 double R0, double R1, double REnd,
                 const Vector3D& U0, const Vector3D& U1, 
                 const Vector3D& UEnd)
    : M3DPrimitive(X), r0(R0), r1(R1), u0(U0), u1(U1), uEnd(UEnd)
{
    //check_r0();
}

M3DQuadPrimitive::
M3DQuadPrimitive(double x0, double x1, double x2, 
                 double R0, double R1, double REnd,
                 double U0x, double U0y, double U0z,
                 double U1x, double U1y, double U1z,
                 double UEndx, double UEndy, double UEndz)
                 : M3DPrimitive(x0, x1, x2), r0(R0), r1(R1), 
                   u0(U0x, U0y, U0z), u1(U1x, U1y, U1z), uEnd(UEndx, UEndy, UEndz)
{ 
    //check_r0();
}

M3DQuadPrimitive::
M3DQuadPrimitive() : M3DPrimitive(), r0(1.0), r1(1.0),
                     u0(1.0, 0.0, 0.0), u1(-1.0, 0.0, 0.0), uEnd(0.0, 1.0, 0.0)
{ 
    //check_r0();
}

M3DQuadPrimitive::
M3DQuadPrimitive(double x0, double x1, double x2,
                 double R, const Quat &Q, double _theta)
    : M3DPrimitive(x0, x1, x2), r0(1.0), r1(1.0),
      u0(1.0, 0.0, 0.0), u1(-1.0, 0.0, 0.0), uEnd(0.0, 1.0, 0.0)
{
    setR(R);
    setQTheta(Q, _theta);
    //check_r0();
}

M3DQuadPrimitive::
M3DQuadPrimitive(const Vector3D &X, double R, const Quat &Q,
                 double _theta)
    : M3DPrimitive(X), r0(1.0), r1(1.0), u0(1.0, 0.0, 0.0),
      u1(-1.0, 0.0, 0.0), uEnd(0.0, 1.0, 0.0)
{
    setR(R);
    setQTheta(Q, _theta);
    //check_r0();
}

// Copy constructor
M3DQuadPrimitive::M3DQuadPrimitive(const M3DQuadPrimitive & prim)
: M3DPrimitive(prim), r0(prim.r0), r1(prim.r1), 
  u0(prim.u0), u1(prim.u1), uEnd(prim.uEnd)
{
    //check_r0();
}

M3DPrimitive & M3DQuadPrimitive::operator=(const M3DPrimitive & unknown_prim)
{
    assert( typeid(*this) == typeid(unknown_prim) );
    const M3DQuadPrimitive& prim = 
        dynamic_cast<const M3DQuadPrimitive&>(unknown_prim);

#ifdef DEBUG
    cout << "M3DQuadPrimitive::operator=()" << endl;
#endif

    M3DPrimitive::copy(static_cast<const M3DPrimitive&>(prim));
    copy(prim);

    //check_r0();

    return (*this);
}

void M3DQuadPrimitive::copy(const M3DQuadPrimitive & prim)
{
    r0 = prim.r0;
    r1 = prim.r1;
    u0 = prim.u0;
    u1 = prim.u1;
    uEnd = prim.uEnd;

    //check_r0();
}

M3DQuadPrimitive* readSrepPrimitive(Registry& registry, const char * newStr);
M3DQuadPrimitive* readLiePrimitive(Registry& registry, const char * newStr);

void M3DQuadPrimitive::setR(double r)
{
    double ratio = r / getR();
    r0 *= ratio;
    r1 *= ratio;

    //check_r0();
}

void M3DQuadPrimitive::setR(int spokeId, double r)
{
	// validate spokeId 

	if( ( spokeId > 2 ) || 
		( ( this->type() == M3D_STANDARD_PRIMITIVE ) && ( spokeId > 1 ) ) )
		return ;

	switch( spokeId ) {

		case 0 :	r0 = r ;
					break ;
		case 1 :	r1 = r ;
					break ;
		case 2 :	(dynamic_cast<M3DEndPrimitive*>(this))->setREnd( r ) ;					
					//(dynamic_cast<M3DEndPrimitive*>(& *this))->setREnd( r ) ;					
					break ;
	}
}


void M3DQuadPrimitive::setQ(const Quat & Q) 
{ 
    setQTheta(Q, getTheta());

#ifdef _DEBUG
    Quat qnew = getQ();
    Quat qdiff = qnew - Q;
    Quat qsum = qnew + Q;
    if (qdiff.length() > .000001 && qsum.length() > .000001) {
        std::cout << "setQ didn't work right." << std::endl;
    }
#endif
}

void M3DQuadPrimitive::setTheta(double _theta) 
{ 
    setQTheta(getQ(), _theta); 
#ifdef _DEBUG
    double th = getTheta();
    if (fabs(th - _theta) > .000001) {
        //std::cout << "setTheta didn't work right. DEL Theta = " << fabs(th - _theta) << std::endl;
    }
#endif
}

void M3DQuadPrimitive::setQTheta(Quat q, double theta)
{
    u0.set(cos(theta), -sin(theta), 0);
    q.rotateVector(u0);

    u1.set(cos(theta), sin(theta), 0);
    q.rotateVector(u1);

    uEnd.set(1, 0, 0);
    q.rotateVector(uEnd);
}

void M3DQuadPrimitive::setU0( Vector3D _u0 ) 
{
	_u0.normalize() ;

	u0 = _u0 ;	
}

void M3DQuadPrimitive::setU1( Vector3D _u1 ) 
{
	_u1.normalize() ;

	u1 = _u1 ;	
}

void M3DQuadPrimitive::setUEnd( Vector3D _uEnd )
{
	_uEnd.normalize() ;

	uEnd = _uEnd ;
}

double M3DQuadPrimitive::getR() const
{
    return sqrt(r0 * r1);
}

//double M3DQuadPrimitive::getR( int spokeId ) 
//{
//	switch( spokeId ) {
//
//		case 0 :	return( this->r0 ) ;
//					break ;
//
//		case 1 :	return( this->r1 ) ;
//					break ;
//		
//		case 2 :	if( (dynamic_cast<M3DQuadEndPrimitive*>(this)) != NULL )
//						return( (dynamic_cast<M3DQuadEndPrimitive*>(this))->rEnd ) ;
//					else
//						return( this->getR() ) ;
//					break ;
//
//		default :	return( this->getR() ) ;
//					break ;
//
//	}
//
//}

Quat M3DQuadPrimitive::getQ() const
{
    Vector3D b = this->getB();
    Vector3D n = this->getN();
    Vector3D bPerp = b.cross(n);
    return Quat(b, n, bPerp);
}

double M3DQuadPrimitive::getTheta() const
{
    return acos(this->getB() * u0);
}

// Returns B, a normalized bisector of u0 and u1.  It is calculated directly from
// u0 and u1.  If this object is an m-rep, then B should equal uEnd.
Vector3D M3DQuadPrimitive::getB() const
{
    Vector3D b = u0 + u1;
    if ( b * b < sqeps ) {
        // printf("Two spokes back to back in M3DQuadPrimitive::getB()\n");
        if (fabs(uEnd * u0) < sqeps) {
            b = uEnd;
        } else {

            // Find an arbitrary vector not too close to u1
            double x = fabs(u1.getX());
            double y = fabs(u1.getY());
            double z = fabs(u1.getZ());
            if (x < y && x < z) b.set(1.0, 0.0, 0.0);
            else if (y < z)     b.set(0.0, 1.0, 0.0);
            else                b.set(0.0, 0.0, 1.0);

            // Make sure b is normal to n, which equals u1 in this case
            b = b - (b * u1) * u1;  
        }
    }
    b.normalize();
    return b;
}

// N is a unit vector that, in a slab, is normal to the 
// medial sheet and in the direction of u1.
Vector3D M3DQuadPrimitive::getN() const
{
    // Choose n to be the difference of the two spokes
    Vector3D n = u1 - u0;
    if (n * n < sqeps) {
        printf("Two spokes identical in M3DQuadPrimitive::getN!\n");

        // Handle the case where u0 = u1;
        // choose a direction normal to the pair
        double x = fabs(u1.getX());
        double y = fabs(u1.getY());
        double z = fabs(u1.getZ());
        if (x < y && x < z) n.set(1.0, 0.0, 0.0);
        else if (y < z)     n.set(0.0, 1.0, 0.0);
        else                n.set(0.0, 0.0, 1.0);
        n = n - (n * u1) * u1; // Make sure n is normal to u1.
    }
    n.normalize();
    return n;
}

void M3DQuadPrimitive::rotateBy(Quat dQ)               
{
    dQ.rotateVector(u0);
    dQ.rotateVector(u1);
    dQ.rotateVector(uEnd);
}

void M3DQuadPrimitive::scaleBy(double mag)
{
    r0 *= mag;
    r1 *= mag;

    //check_r0();
}

void M3DQuadPrimitive::dilate(const double dilationFactor) {			
	
	this->setR( 0, this->getR0() + dilationFactor ) ;
	this->setR( 1, this->getR1() + dilationFactor ) ;

	if( this->type() == M3D_END_PRIMITIVE ) {

		double rEndEarlier = (dynamic_cast <M3DQuadEndPrimitive *> (this))->getREnd() ;
		(dynamic_cast <M3DQuadEndPrimitive *> (this))->setREnd( rEndEarlier + dilationFactor ) ;
	}
}

#undef DEBUG_DIBYENDU_1

#ifndef DEBUG_DIBYENDU_1

// dibyendu
// overloaded scaleBy function to facilitate scaling of a particular spoke in an s-rep
void M3DQuadPrimitive::scaleSpokeBy(int spokeId, double mag) {
	
	// validate spokeId 

	if( ( spokeId > 2 ) || 
		( ( this->type() == M3D_STANDARD_PRIMITIVE ) && ( spokeId > 1 ) ) )
		return ;

	switch( spokeId ) {

		case 0 :	r0 *= mag ;
					break ;
		case 1 :	r1 *= mag ;
					break ;
		case 2 :	double _rEnd = (dynamic_cast<M3DEndPrimitive*>(this))->getREnd() ;
					(dynamic_cast<M3DEndPrimitive*>(this))->setREnd( _rEnd * mag ) ;										
					break ;					
	}
}

// dibyendu
// overloaded rotateSpokeBy function to facilitate rotation of a particular spoke in an s-rep
void M3DQuadPrimitive::rotateSpokeBy( int spokeId, double delTheta, double delPhi ) {
	
	// validate spokeId 

	if( ( spokeId > 2 ) || 
		( ( this->type() == M3D_STANDARD_PRIMITIVE ) && ( spokeId > 1 ) ) )
		return ;

	switch( spokeId ) {

		case 0 :	u0.rotateBy( delTheta, delPhi ) ;
					u0.normalize() ;
					break ;
		case 1 :	u1.rotateBy( delTheta, delPhi ) ;
					u1.normalize() ;
					break ;
		case 2 :	uEnd.rotateBy( delTheta, delPhi ) ;								
					uEnd.normalize() ;
					break ;
	}
}

#endif // DEBUG_DIBYENDU_1

#ifdef DEBUG_DIBYENDU_1

void M3DQuadPrimitive::scaleSpokeBy(int spokeId, double mag) {
	
	// validate spokeId 

	if( ( spokeId > 2 ) || 
		( ( this->type() == M3D_STANDARD_PRIMITIVE ) && ( spokeId > 1 ) ) )
		return ;

	switch( spokeId ) {

		case 0 :	r0 *= mag ;	cout << "Scaling spoke 0 by " << mag << endl ;
					break ;
		case 1 :	r1 *= mag ;	cout << "Scaling spoke 1 by " << mag << endl ;
					break ;
		case 2 :	double _rEnd = (dynamic_cast<M3DEndPrimitive*>(this))->getREnd() ;
					(dynamic_cast<M3DEndPrimitive*>(this))->setREnd( _rEnd * mag ) ;
					cout << "Scaling spoke End by " << mag << endl ;
					break ;					
	}
}

// dibyendu
// overloaded rotateSpokeBy function to facilitate rotation of a particular spoke in an s-rep
void M3DQuadPrimitive::rotateSpokeBy( int spokeId, double delTheta, double delPhi ) {
	
	// validate spokeId 

	if( ( spokeId > 2 ) || 
		( ( this->type() == M3D_STANDARD_PRIMITIVE ) && ( spokeId > 1 ) ) )
		return ;

	switch( spokeId ) {

		case 0 :	u0.rotateBy( delTheta, delPhi ) ;
					u0.normalize() ;
					cout << "Rotating spoke 0 by ( " << delTheta * 180 / 3.1416 << ", " << delPhi << " ) degrees" << endl ;
					break ;
		case 1 :	u1.rotateBy( delTheta, delPhi ) ;
					u1.normalize() ;
					cout << "Rotating spoke 1 by ( " << delTheta * 180 / 3.1416 << ", " << delPhi << " ) degrees" << endl ;
					break ;
		case 2 :	uEnd.rotateBy( delTheta, delPhi ) ;								
					uEnd.normalize() ;
					cout << "Rotating spoke End by ( " << delTheta * 180 / 3.1416 << ", " << delPhi << " ) degrees" << endl ;
					break ;					
	}
}

#endif	// DEBUG_DIBYENDU

void M3DQuadPrimitive::print(ostream & out, char * prefix, bool marked) const
{
    M3DPrimitive::print(out, prefix, marked);

    std::string p;
    if (prefix) p = prefix;

    out << p << "r0, u0 = " << r0 << ", " << u0 << '\n';
    out << p << "r1, u1 = " << r1 << ", " << u1 << '\n';
    out << p << "uMid = " << uEnd << '\n';
}

M3DQuadPrimitive * M3DQuadPrimitive::readPrimitive(Registry& registry, const char * regStr, ...)
{
    char newStr[1024];
    va_list val;
    va_start(val, regStr);
    vsprintf(newStr, regStr, val);
    va_end(val);

    strcat(newStr, ".%s");

    if (registry.hasKey(newStr, "qx")) {
        return readLiePrimitive(registry, newStr);
    } else {
        return readSrepPrimitive(registry, newStr);
    }
}

M3DQuadPrimitive* readSrepPrimitive(Registry& registry, const char * newStr)
{
    const char* typeStr = registry.getStringValue(newStr, NULL, "type");
    double x =  registry.getDoubleValue(newStr, 0.0, "x");
    double y =  registry.getDoubleValue(newStr, 0.0, "y");
    double z =  registry.getDoubleValue(newStr, 0.0, "z");
    double r0 =  registry.getDoubleValue(newStr, 1.0, "r[0]");
    double r1 =  registry.getDoubleValue(newStr, 1.0, "r[1]");
    double rEnd = registry.getDoubleValue(newStr, 1.0, "r[2]");
    double u0x = registry.getDoubleValue(newStr, 1.0, "ux[0]"); 
    double u0y = registry.getDoubleValue(newStr, 1.0, "uy[0]"); 
    double u0z = registry.getDoubleValue(newStr, 1.0, "uz[0]"); 
    double u1x = registry.getDoubleValue(newStr, 1.0, "ux[1]"); 
    double u1y = registry.getDoubleValue(newStr, 1.0, "uy[1]"); 
    double u1z = registry.getDoubleValue(newStr, 1.0, "uz[1]"); 
    double uEndx = registry.getDoubleValue(newStr, 1.0, "ux[2]"); 
    double uEndy = registry.getDoubleValue(newStr, 1.0, "uy[2]"); 
    double uEndz = registry.getDoubleValue(newStr, 1.0, "uz[2]"); 

    short selected = registry.getIntValue(newStr, 0, "selected");

    M3DQuadPrimitive * primPtr;
    if(typeStr != NULL && 
       strcmp(typeStr, M3D_END_PRIMITIVE_STR) == 0) 
    {
        primPtr = new M3DQuadEndPrimitive(
            x, y, z, r0, r1, rEnd, 
            u0x, u0y, u0z, u1x, u1y, u1z, 
            uEndx, uEndy, uEndz);
    } else {
        primPtr = new M3DQuadPrimitive(
            x, y, z, r0, r1, rEnd, 
            u0x, u0y, u0z, u1x, u1y, u1z, 
            uEndx, uEndy, uEndz);
    }

    if(selected & 0x1) primPtr->select();
    if(selected & 0x2) primPtr->deselectForRegularity();

    return primPtr;
}

M3DQuadPrimitive* readLiePrimitive(Registry& registry, const char * newStr)
{

    const char* typeStr = registry.getStringValue(newStr, NULL, "type");
    double x =  registry.getDoubleValue(newStr, 0.0, "x");
    double y =  registry.getDoubleValue(newStr, 0.0, "y");
    double z =  registry.getDoubleValue(newStr, 0.0, "z");
    double r =  registry.getDoubleValue(newStr, 1.0, "r");
    Quat q;
    q.setX(registry.getDoubleValue(newStr, 0.0, "qx"));
    q.setY(registry.getDoubleValue(newStr, 0.0, "qy"));
    q.setZ(registry.getDoubleValue(newStr, 0.0, "qz"));
    q.setW(registry.getDoubleValue(newStr, 1.0, "qw"));
    double theta = registry.getDoubleValue(newStr, 0.0, "theta") * R_DEGREES_TO_RADIANS;
    short selected = registry.getIntValue(newStr, 0, "selected");

    M3DQuadPrimitive * primPtr;
    if(typeStr != NULL && strcmp(typeStr, M3D_END_PRIMITIVE_STR) == 0) {
        double elongation = registry.getDoubleValue(newStr, 1.0, "elongation");
        primPtr = new M3DQuadEndPrimitive(x, y, z, r, q, theta, elongation);
        // std::cout << *primPtr << std::endl;
    } else {
        primPtr = new M3DQuadPrimitive(x, y, z, r, q, theta);
        // std::cout << *primPtr << std::endl;
    }

    if(selected & 0x1) primPtr->select();
    if(selected & 0x2) primPtr->deselectForRegularity();

    return primPtr;
}

void M3DQuadPrimitive::writePrimitive(Registry& registry, const char * regStr, ...) const
{
    char newStr[1024];

    va_list val;
    va_start(val, regStr);
    vsprintf(newStr, regStr, val);
    va_end(val);

    strcat(newStr, ".%s");

    if(type() == M3D_STANDARD_PRIMITIVE) {
        registry.setStringValue(newStr, (char*)M3D_STANDARD_PRIMITIVE_STR, "type");
    } else if(type() == M3D_END_PRIMITIVE) {
        registry.setStringValue(newStr, (char*)M3D_END_PRIMITIVE_STR, "type");
        registry.setDoubleValue(
            newStr, dynamic_cast<const M3DEndPrimitive*>(this)->getREnd(), "r[2]");
    }

    registry.setDoubleValue(newStr, getX().getX(), "x");
    registry.setDoubleValue(newStr, getX().getY(), "y");
    registry.setDoubleValue(newStr, getX().getZ(), "z");
    registry.setDoubleValue(newStr, getR0(), "r[0]");
    registry.setDoubleValue(newStr, getR1(), "r[1]");
    registry.setDoubleValue(newStr, getU0().getX(), "ux[0]");
    registry.setDoubleValue(newStr, getU0().getY(), "uy[0]");
    registry.setDoubleValue(newStr, getU0().getZ(), "uz[0]");
    registry.setDoubleValue(newStr, getU1().getX(), "ux[1]");
    registry.setDoubleValue(newStr, getU1().getY(), "uy[1]");
    registry.setDoubleValue(newStr, getU1().getZ(), "uz[1]");
    registry.setDoubleValue(newStr, getUEnd().getX(), "ux[2]");
    registry.setDoubleValue(newStr, getUEnd().getY(), "uy[2]");
    registry.setDoubleValue(newStr, getUEnd().getZ(), "uz[2]");
    registry.setIntValue(newStr, getSelectionFlags(), "selected");
}


bool M3DQuadPrimitive::atomInterp( double t, const M3DPrimitive* m1, const M3DPrimitive* m2) 
{
	if( m1->type() == M3D_END_PRIMITIVE && m2->type() == M3D_END_PRIMITIVE && type() != M3D_END_PRIMITIVE) {
		if( globalVerbosity >= 0 ) {
			cout << "Warning: Interpolating two end primitives and saving result in an ordinary primitive." << endl;
		}
	}
	QuadSymPrimitive sM1(dynamic_cast<const M3DQuadPrimitive*>(m1));
	QuadSymPrimitive sM2(dynamic_cast<const M3DQuadPrimitive*>(m2));
	QuadSymPrimitive deltaSM, interpolatedSM;
	Vector3D rotatedSpoke;

	deltaSM.x=(sM2.x-sM1.x)*t;
	deltaSM.r=(log(sM2.r)-log(sM1.r))*t;
	deltaSM.elongation=(log(sM2.elongation)-log(sM1.elongation))*t;

	interpolatedSM.x=sM1.x+deltaSM.x;
	interpolatedSM.r=exp(log(sM1.r)+deltaSM.r);
	interpolatedSM.elongation=exp(log(sM1.elongation)+deltaSM.elongation);

	Quat q=ShapeSpace::S2::rotationToOrigin(sM1.n0);
	q=q.conj();
	rotatedSpoke = sM2.n0;
	q.conj().rotateVector(rotatedSpoke);
	Vector2D logDataN0=ShapeSpace::S2::Log(rotatedSpoke);
	logDataN0*=t;
	interpolatedSM.n0=ShapeSpace::S2::Exp(logDataN0);
	q.rotateVector(interpolatedSM.n0);

	q=ShapeSpace::S2::rotationToOrigin(sM1.n1);
	q=q.conj();
	rotatedSpoke = sM2.n1;
	q.conj().rotateVector(rotatedSpoke);
	Vector2D logDataN1=ShapeSpace::S2::Log(rotatedSpoke);
	logDataN1*=t;
	interpolatedSM.n1=ShapeSpace::S2::Exp(logDataN1);
	q.rotateVector(interpolatedSM.n1);
	
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

// dibyendu
// Modified this function to handle s-reps 

double M3DQuadPrimitive::atomSquareDistance( const M3DPrimitive* m, const double* radius ) const 
{
	double effectiveR = (radius == NULL)
		? sqrt(m->getR()*getR())
		: *radius;		

	M3DQuadPrimitive * quadM = dynamic_cast <M3DQuadPrimitive *> (const_cast <M3DPrimitive*> (m)) ;
		
	// I - differences in position	

	Vector3D delX = this->getX() - quadM->getX() ;

	// II - differences in spoke radii

	double delR0 = ( log(this->getR0()) - log(quadM->getR0()) ) * effectiveR ;
	double delR1 = ( log(this->getR1()) - log(quadM->getR1()) ) * effectiveR ;

	double delREnd = 0.0 ;

	if( ( this->type() == M3D_END_PRIMITIVE ) && ( quadM->type() == M3D_END_PRIMITIVE ) ) {

		double rEnd0 = ( dynamic_cast <M3DEndPrimitive *> (const_cast <M3DQuadPrimitive *> (this)) )->getREnd() ;
		double rEnd1 = ( dynamic_cast <M3DEndPrimitive *> (quadM) )->getREnd() ;

		delREnd = ( log( rEnd0 ) - log( rEnd1 ) ) * effectiveR ;
	}

	// III - differences in spoke directions

	Quat q ;

	Vector3D u0_0 = this->getU0() ;
	Vector3D u0_1 = quadM->getU0() ;

	Vector3D u1_0 = this->getU1() ;
	Vector3D u1_1 = quadM->getU1() ;

	Vector3D uEnd_0 = this->getUEnd() ;
	Vector3D uEnd_1 = quadM->getUEnd() ;
	
	q = ShapeSpace::S2::rotationToOrigin( u0_0 ) ;
	q.normalize() ;
	q.rotateVector( u0_1 ) ;

	Vector2D logDataU0 = ShapeSpace::S2::Log( u0_1 ) * effectiveR ;

	q = ShapeSpace::S2::rotationToOrigin( u1_0 ) ;
	q.normalize() ;
	q.rotateVector( u1_1 ) ;

	Vector2D logDataU1 = ShapeSpace::S2::Log( u1_1 ) * effectiveR ;

	Vector2D logDataUEnd = Vector2D(0.0, 0.0) ;	

	if( ( this->type() == M3D_END_PRIMITIVE ) && ( quadM->type() == M3D_END_PRIMITIVE ) ) {

		q = ShapeSpace::S2::rotationToOrigin( uEnd_0 ) ;
		q.normalize() ;
		q.rotateVector( uEnd_1 ) ;

		logDataUEnd = ShapeSpace::S2::Log( uEnd_1 ) * effectiveR ;
	}

	// Sum of square of each difference = total difference

	double sumSquare = 0;

	sumSquare += delX.norm() * delX.norm() ;
	sumSquare += delR0 * delR0 + delR1 * delR1 + delREnd * delREnd ;
	sumSquare += logDataU0.norm() * logDataU0.norm() 
			  +  logDataU1.norm() * logDataU1.norm()
			  +  logDataUEnd.norm() * logDataUEnd.norm() ;

	//if( thatM != NULL )
	//	delete thatM ;
	
	return sumSquare;
}


bool M3DQuadPrimitive::atomAverage( const int mNum, M3DPrimitive** mList,
	const double* _weights)
{
	assert( mNum > 0 && mList != NULL );

	double* weights	= new double[mNum];
	int i;

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

	QuadSymPrimitive mean;
	mean.elongation	= 0.0;
	mean.r			= 0.0;
	for( i = 0; i != mNum; ++i ) {
		QuadSymPrimitive sym	= QuadSymPrimitive(dynamic_cast<const M3DQuadPrimitive*>(mList[i]));
		sym.Log();
		mean.x          += sym.x * weights[i];
		mean.r			+= sym.r * weights[i];
		mean.elongation	+= sym.elongation * weights[i];
	}
	mean.Exp();

	Vector3D* spokes	= new Vector3D[mNum];
	// Compute average of spoke n0.
	for( i = 0; i != mNum; ++i ) {
        spokes[i] = dynamic_cast<const M3DQuadPrimitive*>(mList[i])->getU0();
	}
	mean.n0	= ShapeSpace::S2::mean(spokes, mNum, weights);

	// Compute average of spoke n1.
	for( i = 0; i != mNum; ++i ) {
		spokes[i] = dynamic_cast<const M3DQuadPrimitive*>(mList[i])->getU1();
	}
	mean.n1	= ShapeSpace::S2::mean(spokes, mNum, weights);

	mean.convert2Lie(this);

	// Clean up.
	delete[] weights;
	delete[] spokes;

	return true;
}


/*
 * The logic for this code is copied from GeodesicSym::atomToSymAtom.
 */
M3DQuadPrimitive::QuadSymPrimitive::QuadSymPrimitive( const M3DQuadPrimitive* prim ) :
ignorePrimType(false) 
{
    x	= prim->getX();
    r	= prim->getR();
    n0	= prim->getNormalizedY0();
    n1	= prim->getNormalizedY1();

    if(!ignorePrimType && prim->type()==M3D_END_PRIMITIVE) {
        elongation = (dynamic_cast<const M3DQuadEndPrimitive *>(prim))->
            getElongation();
    }
    else {
        elongation = 1.0;
    }
}

/*
 * The logic for this code is copied from GeodesicSym::symAtomToAtom.
 */
bool M3DQuadPrimitive::QuadSymPrimitive::convert2Lie( M3DPrimitive* _m ) const {
	Quat q;
	double theta;

	M3DQuadPrimitive* m	= dynamic_cast<M3DQuadPrimitive*>(_m);

	m->setX(x);
	m->setR(r);

    Vector3D b, n, bPerp;

    // Choose n to be the difference of the two spokes
	// (n is in the same direction as Y1, not Y0?)

    if (n1 == n0) {
        printf("Two spokes identical in QuadSymPrimitive::convert2Lie!\n");

        // Handle the case where n0 = n1
        if ( fabs(n0.getX()) < fabs(n0.getY()) )
            if ( fabs(n0.getX()) < fabs(n0.getZ()) )
                n.set(1.0, 0.0, 0.0);
            else
                n.set(0.0, 0.0, 1.0);
        else
            if ( fabs(n0.getY()) > fabs(n0.getZ()) )
                n.set(0.0, 0.0, 1.0);
            else
                n.set(0.0, 1.0, 0.0);
		n = n - (n * n1) * n1;
    }
    else
        n = n1 - n0;
    n.normalize();

    b = n0 + n1;
    // Handle the case where n0 = -n1
    if (b.norm() == 0) {
		//printf("Two spokes back to back in QuadSymPrimitive::convert2Lie!\n");

		// If n0 and n1 are back to back,
		// pick a set of local frames
        if ( fabs(n0.getX()) < fabs(n0.getY()) )
            if ( fabs(n0.getX()) < fabs(n0.getZ()) )
                b.set(1.0, 0.0, 0.0);
            else
                b.set(0.0, 0.0, 1.0);
        else
            if ( fabs(n0.getY()) > fabs(n0.getZ()) )
                b.set(0.0, 0.0, 1.0);
            else
                b.set(0.0, 1.0, 0.0);
    }
    b = b - (b * n) * n;
    b.normalize();

    bPerp = b.cross(n);

    theta = acos(b * n0);
    q = Quat(b, n, bPerp);

	m->setTheta(theta);
	m->setQ(q);
	
	if(!ignorePrimType && m->type() == M3D_END_PRIMITIVE) {
		(dynamic_cast<M3DQuadEndPrimitive *>(m))->setElongation( elongation );
	}
	
	return true;
}


/*
 * The logic for this code is based on GeodesicSym::symAtomExp()
 */
void M3DQuadPrimitive::QuadSymPrimitive::Exp() {
	// x stays the same.
	r			= exp(r);
	elongation	= exp(elongation);
	n0			= ShapeSpace::S2::Exp(tn0);
	n1			= ShapeSpace::S2::Exp(tn1);
}

/*
 * The logic for this code is based on GeodesicSym::symAtomLogMap()
 */
void M3DQuadPrimitive::QuadSymPrimitive::Log() {
	// x stays the same.
	r			= log(r);
	elongation	= log(elongation);
	tn0			= ShapeSpace::S2::Log(n0);
	tn1			= ShapeSpace::S2::Log(n1);
}

/*
 * The logic for this code is based on GeodesicSym::symAtomDiff()
 */
SymPrimitive& M3DQuadPrimitive::QuadSymPrimitive::operator-=(const SymPrimitive& _b) {
	const QuadSymPrimitive& b	= dynamic_cast<const QuadSymPrimitive&>(_b);

	x	= x - b.x;
	r	= r/b.r;
	elongation	= elongation/b.elongation;

	Vector3D n0(b.n0), n1(b.n1);

	Quat q=ShapeSpace::S2::rotationToOrigin(n0);
	q.normalize();
	q.rotateVector(this->n0);


	q=ShapeSpace::S2::rotationToOrigin(n1);
	q.normalize();
	q.rotateVector(this->n1);


	this->n0.normalize();
	this->n1.normalize();

	return *this;


}


/*
 * 
 */
SymPrimitive& M3DQuadPrimitive::QuadSymPrimitive::operator+=(const SymPrimitive& _b) {
	const QuadSymPrimitive& b	= dynamic_cast<const QuadSymPrimitive&>(_b);

	x	= x + b.x;
	r	= r * b.r;
	elongation	= elongation * b.elongation;//???xiaoxiao?

	Vector3D n0(b.n0), n1(b.n1);

	Quat q=ShapeSpace::S2::rotationFromOrigin(n0);
	q.normalize();
	q.rotateVector(this->n0);


	q=ShapeSpace::S2::rotationFromOrigin(n1);
	q.normalize();
	q.rotateVector(this->n1);


	this->n0.normalize();
	this->n1.normalize();

	return *this;
}


bool  M3DQuadPrimitive::QuadSymPrimitive::projectAtom2Vec( Vector & vecPrim){
	//prjected into atom pga space
	
    Log();
    
	vecPrim = Vector(9,x.getX(),x.getY(),x.getZ(),r,tn1.getX(),tn1.getY(),tn0.getX(),tn0.getY(),elongation);
   // tn0 = U(-1);   tn1 = U(+1)
	return true;

}

bool  M3DQuadPrimitive::QuadSymPrimitive::vec2Atom( const Vector & vecPrim){
    
	x = Vector3D(vecPrim(0),vecPrim(1),vecPrim(2));
	r = vecPrim(3);
	tn1 = Vector2D(vecPrim(4),vecPrim(5));
	tn0 = Vector2D(vecPrim(6),vecPrim(7));
	elongation = vecPrim(8);
	
	Exp();
	
	return true;
}

// DEBUGGING

Vector3D M3DQuadPrimitive::dgetExtendedB() const
{
    // Vector to represent the x-axis
    Vector3D b(1, 0, 0);

    // Perform a rotation of x-axis by our orientation q
    getQ().rotateVector(b);

    // Scale b vector by the radius
    b *= getR();

//	if( getTheta() > R_HALF_PI ) {
//		b	= -b;
//	}

    return b;
}

Vector3D M3DQuadPrimitive::dgetB() const
{
    // Vector to represent the x-axis
    Vector3D b(1, 0, 0);

    // Perform a rotation of x-axis by our orientation q
    getQ().rotateVector(b);

    return b;
}

Vector3D M3DQuadPrimitive::dgetY0() const
{
    double r = getR();
    double theta = getTheta();
    // Vector to represent the x-axis rotated clockwise by theta
    Vector3D y0(r * cos(theta), -r * sin(theta), 0);

    // perform a rotation of y0 by our orientation q
    getQ().rotateVector(y0);

    return y0;
}

Vector3D M3DQuadPrimitive::dgetY1() const
{
    double r = getR();
    double theta = getTheta();
    // Vector to represent the x-axis rotated counterclockwise by theta
    Vector3D y1(r * cos(theta), r * sin(theta), 0);

    // perform a rotation of y1 by our orientation q
    getQ().rotateVector(y1);

    return y1;
}

Vector3D M3DQuadPrimitive::dgetN() const
{
    // Vector to represent the y-axis
    Vector3D n(0, 1, 0);

    // perform a rotation of y-axis by our orientation q
    getQ().rotateVector(n);

    return n;
}

Vector3D M3DQuadPrimitive::dgetBPerp() const
{
    // Vector to represent the z-axis
    Vector3D bPerp(0, 0, 1);

    // perform a rotation of z-axis by our orientation q
    getQ().rotateVector(bPerp);

    return bPerp;
}

Vector3D M3DQuadPrimitive::dgetNormalizedY0() const
{
    // Vector to represent the x-axis rotated clockwise by theta
    double theta = getTheta();
    Vector3D y0(cos(theta), -sin(theta), 0);

    // perform a rotation of y0 by our orientation q
    getQ().rotateVector(y0);

    return y0;
}

Vector3D M3DQuadPrimitive::dgetNormalizedY1() const
{
    // Vector to represent the x-axis rotated counterclockwise by theta
    double theta = getTheta();
    Vector3D y1(cos(theta), sin(theta), 0);

    // perform a rotation of y1 by our orientation q
    getQ().rotateVector(y1);

    return y1;
}
