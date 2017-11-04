#ifndef XFERLIST_H
#define XFERLIST_H

/**
 * We need to get *rid* of both Xferlist and Diatomgrid and use M3DFigure directly
 * in seurat, however, it's too much of work to try and rectify other's mistakes,
 * so I'll simply include a pointer to M3DFigure in both these classes and
 * henceforth everyone should use functions/data from that pointer for all new things.
 *
 * - Rohit Saboo
 */

/********************************************************************************/
/*                                                                              */
/*      File    :  Xferlist.H													*/
/*                                                                              */
/*      Description:  function to load from M3DQuadfigure into Xferlist to		*/
/*			transfer data to Diatomgrids in ThallCode							*/
/*																				*/
/*      Project :  Seurat                                                       */
/*                                                                              */
/*      Author  :  A. Thall                                                     */
/*                                                                              */
/*      Date    :  28. May 2000													*/
/*                                                                              */
/*      Modifications:                                                          */
/*			3. June -- added NPC_Bpoint and Bpoint to declarations.  This file	*/
/*				will be holder for all includes outside of ThallCode namespace.	*/
/*			July 2000 -- added (u, v, t) to Bpoint struct (Joshi)				*/
/********************************************************************************/

// forward declare this for a pointer
class M3DQuadFigure;
class M3DTubeFigure;
class M3DFigure;

typedef struct {
	double X_p[3];
	double X_q[4];	// assume q = (<x, y, z>, w) is our quaternion indexed [0, 1, 2, 3]
	double X_r;
	double X_theta;	// assume pi/2 is parallel sides (NOT my convention, where phi=0 is parallel)
	double X_eta;	// elongation factor for end atoms.  Assumed = 1 there and elsewhere
	double X_rho;	// nobody else uses rho, but I make it 0.25 by default

	// dibyendu - on 06/17/2011, added the following variables to enable the use of s-reps
	// I donot know about tubular m-reps and donot know how this change would affect tube m-reps
	// The old variables still remain to maintain compatibility
	// This structure should probably be removed and the M3DPrimitive structures should be used instead

	double X_r0 ;
	double X_r1 ;
	double X_rEnd ;

	double X_u0[3] ;
	double X_u1[3] ;
	double X_uEnd[3] ;
        

} XferAtom;

typedef struct {
	XferAtom *atomlist;
	int numrows, numcols;	// assume in row major order---need to reverse this for my code
	int tolerance;
	const M3DFigure* figure;
} Xferlist;

Xferlist *convertM3DtoXfer(const M3DQuadFigure *thisfigure);
Xferlist *convertM3DtoXfer(const M3DTubeFigure *thisfigure);
Xferlist *convertM3DtoXfer(const M3DFigure *figure);

class Bpoint {
public:
	double pnt[3];
	double norm[3];
	double rad;
	double u;
	double v;
	double t;

	void print(std::ostream & os = std::cout) const {
		os << "pnt = " << pnt[0] << ", " << pnt[1] << ", " << pnt[2] << "   ";
		os << "norm = " << norm[0] << ", " << norm[1] << ", " << norm[2] << '\n';
		os << "\trad = " << rad << "   uvt = " << u << ", " << v << ", " << t << '\n';
	}
	Bpoint& operator=(const Bpoint& bpnt) {
		pnt[0] = bpnt.pnt[0];
		pnt[1] = bpnt.pnt[1];
		pnt[2] = bpnt.pnt[2];
		norm[0] = bpnt.norm[0];
		norm[1] = bpnt.norm[1];
		norm[2] = bpnt.norm[2];
		rad = bpnt.rad;
		u = bpnt.u;
		v = bpnt.v;
		t = bpnt.t;

		return (*this);
	}

	Bpoint& setvals(double px, double py, double pz, double nx, double ny, double nz, double rval, \
		double uval, double vval, double tval) 
	{
		pnt[0] = px; pnt[1] = py; pnt[2] = pz; norm[0] = nx; norm[1] = ny; norm[2] = nz;
	    rad = rval; u = uval; v = vval; t = tval;
		return (*this);
	}

	bool operator==(const Bpoint& bpnt) {
		return ((u == bpnt.u) && (v == bpnt.v) && (t == bpnt.t));
	}

	bool operator!=(const Bpoint& bpnt) const{
		return ((u != bpnt.u) || (v != bpnt.v) || (t != bpnt.t));
	}
};

//		a variation class from the Bpoint is defined to accomodate another set of UVTR from
//		the sub-figure for the vertices on the intersection curve
class Bpoint2 
{
public:
	double pnt[3];
	double norm[3];
	double rad, rad_2;
	double u, u_2;
	double v, v_2;
	double t, t_2;

	void print(std::ostream & os = std::cout) const {
		os << "pnt = " << pnt[0] << ", " << pnt[1] << ", " << pnt[2] << "   ";
		os << "norm = " << norm[0] << ", " << norm[1] << ", " << norm[2] << '\n';
		os << "\trad = " << rad << "   uvt = " << u << ", " << v << ", " << t << '\n';
		os << "\trad2 = " << rad_2 << "   uvt_2 = " << u_2 << ", " << v_2 << ", " << t_2 << '\n';
	}

	Bpoint2& operator=(const Bpoint2& bpnt) {
		pnt[0] = bpnt.pnt[0];
		pnt[1] = bpnt.pnt[1];
		pnt[2] = bpnt.pnt[2];
		norm[0] = bpnt.norm[0];
		norm[1] = bpnt.norm[1];
		norm[2] = bpnt.norm[2];
		rad = bpnt.rad;
		u = bpnt.u;
		v = bpnt.v;
		t = bpnt.t;

		rad_2 = bpnt.rad_2;
		u_2 = bpnt.u_2;
		v_2 = bpnt.v_2;
		t_2 = bpnt.t_2;

		return (*this);
	}

	//		to assign a Bpoint object to a Bpoint2 object, only the values of
	//		pnt/norm/rad are tranferred so far!
	Bpoint2& operator=(const Bpoint& bpnt) {
		pnt[0] = bpnt.pnt[0];
		pnt[1] = bpnt.pnt[1];
		pnt[2] = bpnt.pnt[2];
		norm[0] = bpnt.norm[0];
		norm[1] = bpnt.norm[1];
		norm[2] = bpnt.norm[2];
		rad = bpnt.rad;
		return (*this);
	}

	Bpoint2& setvals(double px, double py, double pz, double nx, double ny, double nz, double rval, \
		double uval, double vval, double tval) 
	{
		pnt[0] = px; pnt[1] = py; pnt[2] = pz; norm[0] = nx; norm[1] = ny; norm[2] = nz;
	    rad = rval; u = uval; v = vval; t = tval;
		return (*this);
	}

	Bpoint2& setvals(double px, double py, double pz, double nx, double ny, double nz, double rval, \
		double uval, double vval, double tval, double rval2, double uval2, double vval2, double tval2) 
	{
		pnt[0] = px; pnt[1] = py; pnt[2] = pz; norm[0] = nx; norm[1] = ny; norm[2] = nz;
	    rad = rval; u = uval; v = vval; t = tval;
	    rad_2 = rval2; u_2 = uval2; v_2 = vval2; t_2 = tval2;
		return (*this);
	}

	Bpoint2& setvals(const Bpoint& bpnt) 
	{
		pnt[0] = bpnt.pnt[0]; pnt[1] = bpnt.pnt[1]; pnt[2] = bpnt.pnt[2]; 
		norm[0] = bpnt.norm[0]; norm[1] = bpnt.norm[1]; norm[2] = bpnt.norm[2];
	    rad = bpnt.rad; u = bpnt.u; v = bpnt.v; t = bpnt.t;
		return (*this);
	}

	Bpoint2& setvals(const Bpoint2& bpnt) 
	{
		pnt[0] = bpnt.pnt[0]; pnt[1] = bpnt.pnt[1]; pnt[2] = bpnt.pnt[2]; 
		norm[0] = bpnt.norm[0]; norm[1] = bpnt.norm[1]; norm[2] = bpnt.norm[2];
	    rad = bpnt.rad; u = bpnt.u; v = bpnt.v; t = bpnt.t;
	    rad_2 = bpnt.rad_2; u_2 = bpnt.u_2; v_2 = bpnt.v_2; t_2 = bpnt.t_2;
		return (*this);
	}
};

class MyPoint2D
{
public:
	MyFloat x, y;

public:
	//////////////////////////////////////////////////////////////
	// Constructors
	MyPoint2D()
	{
		x=0; y=0;
	}
	MyPoint2D(const MyFloat v)
	{
		x=v; y=v;
	}
	MyPoint2D(const MyFloat xx, const MyFloat yy)
	{
		x=xx; y=yy;
	}
	MyPoint2D(const MyPoint2D &A)
	{
		x=A.x; y=A.y;
	}
	MyPoint2D(MyFloat *v)
	{
		x=v[0];	y=v[1];
	}


	//////////////////////////////////////////////////////////////
	// Setters and Getters
	void Set(MyFloat xx, MyFloat yy)
	{
		x=xx; y=yy;
	}
	MyFloat GetX()
	{
		return x;
	}
	MyFloat GetY()
	{
		return y;
	}

	//////////////////////////////////////////////////////////////
	// Overloaded Operators
	void operator=(const MyPoint2D &A)
	{
		x=A.x; y=A.y;
	}
	void operator=(const MyFloat v)
	{
		x=v; y=v;
	}
	MyPoint2D operator+(const MyPoint2D &A)
	{
		return MyPoint2D(x+A.x, y+A.y);
	}
	void operator+=(const MyPoint2D &A)
	{
		x+=A.x; y+=A.y;
	}
	MyPoint2D operator-(const MyPoint2D &A)
	{
		return MyPoint2D(x-A.x, y-A.y);
	}
	void operator-=(const MyPoint2D A)
	{
		x-=A.x; y-=A.y;
	}
	friend MyPoint2D operator*(const MyFloat a, const MyPoint2D &A)
	{
		return MyPoint2D(A.x*a, A.y*a);
	}
	friend MyPoint2D operator*(const MyPoint2D &A, const MyFloat a)
	{
		return MyPoint2D(A.x*a, A.y*a);
	}
	// dot product
	friend MyFloat operator*(const MyPoint2D &A, const MyPoint2D &B)
	{
		//return (x*A.x+y*A.y+z*A.z);
		return B.x*A.x+B.y*A.y;
	}
	void operator*=(const MyFloat a)
	{
		x*=a; y*=a;;
	}
	MyPoint2D operator/(const MyFloat a)
	{
		return MyPoint2D(x/a, y/a);
	}
	void operator/=(const MyFloat a)
	{
		x/=a; y/=a;;
	}
	MyPoint2D operator-(void)
	{
		return MyPoint2D(-x, -y);
	}
	//friend MyPoint2D operator^(const MyPoint2D &A, const MyFloat a)
	//{
	//	return MyPoint2D((MyFloat)pow(A.x, a), (MyFloat)pow(A.y, a));
	//}
	//friend MyPoint2D operator^(const MyPoint2D &A, const MyPoint2D &B)
	//{
	//	return MyPoint2D(A.y*B.z-A.z*B.y, A.z*B.x-A.x*B.z, A.x*B.y-A.y*B.x);
	//}
	MyPoint2D operator|(const MyPoint2D &A)
	{
		return MyPoint2D(x*A.x, y*A.y);
	}

	//////////////////////////////////////////////////////////////
	// Other functions
	//MyFloat Mag(void)
	//{
	//	return (MyFloat)sqrt(x*x+y*y);
	//}
	MyFloat MagSquare(void)
	{
		return x*x+y*y;
	}
	MyFloat Sum(void)
	{
		return x+y;
	}
	//MyFloat Distance(const MyPoint2D &A)
	//{
	//	return (MyFloat)sqrt((A.x-x)*(A.x-x)+(A.y-y)*(A.y-y));
	//}
	MyFloat DistanceSquare(const MyPoint2D &A)
	{
		return ((A.x-x)*(A.x-x)+(A.y-y)*(A.y-y));
	}

	MyFloat CircularDistance(const MyPoint2D &A)
	{
		MyFloat vDelta, tDelta;
		vDelta=x-A.x>=0?x-A.x:-(x-A.x);
		tDelta=y-A.y>=0?y-A.y:-(y-A.y);
		return vDelta+tDelta;
		//return ((A.x-x)*(A.x-x)+(A.y-y)*(A.y-y));
	}

	/*
	// axb = a2b3-a3b2 a3b1-a1b3 a1b2-a2b1
	static MyPoint2D Normal(MyPoint2D &v0, MyPoint2D &v1, MyPoint2D &v2)
	{
		MyPoint2D a(v1.x-v0.x, v1.y-v0.y, v1.z-v0.z);
		MyPoint2D b(v1.x-v2.x, v1.y-v2.y, v1.z-v2.z);
		return Normalize(MyPoint2D(a.y*b.z-a.z*b.y, a.z*b.x-a.x*b.z, a.x*b.y-a.y*b.x));
	}
	static void Normal(MyFloat *n, MyFloat* v0, MyFloat* v1, MyFloat* v2)
	{
		MyPoint2D a(v1[0]-v0[0], v1[1]-v0[1], v1[2]-v0[2]);
		MyPoint2D b(v1[0]-v2[0], v1[1]-v2[1], v1[2]-v2[2]);
		MyPoint2D temp = Normalize(MyPoint2D(a.y*b.z-a.z*b.y, a.z*b.x-a.x*b.z, a.x*b.y-a.y*b.x));
		n[0]=temp.x;
		n[1]=temp.y;
		n[2]=temp.z;
	}

	void Normalize(void)
	{
		MyFloat mag;//=Mag();		Uninitialized variable!
		if(mag>0)
		{
	//		x/=mag;
	//		y/=mag;
		}
		else
			fprintf(stderr, "\tMyPoint2D::Normalize, the magnitude is 0\n");
	}
	static MyPoint2D Normalize(MyPoint2D &A)
	{
		MyFloat mag;//=A.Mag();		Uninitialized variable!
		if(mag>0)
			return MyPoint2D(A.x/mag, A.y/mag);
		else
			return MyPoint2D(0, 0);
	}
	*/

	//////////////////////////////////////////////////////////////
	// Debugging tools
	void print()
	{
		std::cout << "The x/y are " << x << '/' << y << "%f\n";
	}
};

#endif

