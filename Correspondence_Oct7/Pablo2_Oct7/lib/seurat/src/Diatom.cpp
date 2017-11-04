/********************************************************************************/
/*																				*/
/*  	File	:  Diatom.cpp													*/
/*																				*/
/*	Description:  class functions for medial atom base class					*/
/*	    Defined functions are:													*/
/*		void Diatom::rot_q(Quat qrot, Diatom& newDA);							*/
/*		void Diatom::rot_q(DbMatrix3& matrot, Diatom& newDA);					*/
/*		void Diatom::rot_v1(double theta, Diatom& newDA);						*/
/*		void Diatom::rot_v2(double theta, Diatom& newDA);						*/
/*		void Diatom::rot_Diatom(Quat q_rot, Diatom& newDA);						*/
/*		void Diatom::rot_Diatom(DbMatrix3& matrot, Diatom& newDA);				*/
/*		void xlate_p(DbVector3 xlation, Diatom& newDA);							*/
/*		void xlate_p(double x, double y, double z, Diatom& newDA);				*/
/*		void zoom_r(double delta_r, Diatom& newDA);								*/
/*		Diatom& selfrot_q(Quat qrot);											*/
/*		Diatom& selfrot_q(DbMatrix3& matrot); 									*/
/*		Diatom& selfrot_v1(double theta);										*/
/*		Diatom& selfrot_v2(double theta);										*/
/*		Diatom& selfrot_Diatom(Quat q_rot); 									*/
/*		Diatom& selfrot_Diatom(DbMatrix3& matrot);								*/
/*		Diatom& selfxlate_p(DbVector3 xlation);									*/
/*		Diatom& selfxlate_p(double x, double y, double z);						*/
/*		Diatom& selfzoom_r(double delta_r);										*/
/*		DbVector3 bvector()														*/
/*		DbVector3 b_val()														*/
/*		DbVector3 bperp_val()													*/
/*		DbVector3 n_val()														*/
/*		DbVector3 edgeboundary_t()												*/
/*		void GLrender()															*/
/*																				*/
/*	Project :  Seurat															*/
/*	Author  :  A. Thall															*/
/*	Date	:  4. February 1998													*/
/*																				*/
/*	Modifications:  															*/
/*		22. Sept 98 -- added bvector() code										*/
/*		29. May 99 -- modifed bvector() to return x-axis of q even				*/
/*			when phi = 0.														*/
/*		19. May 2000 -- added alpha variable describing multiplier				*/
/*			for r in bvector direction to get implied boundary at				*/
/*			edge, corner, and end-Diatoms.										*/
/*		30. Oct 00 -- added from3DtoDframe_basis() to give coordinates of a		*/
/*			vector in R^3 in [b, bperp, n] coordinates.  Useful for projection	*/
/*			onto tangent space of medial surface at Diatom, etc.				*/
/*		14. Nov -- added GLrender() routine, to draw Diatom as  tiny sphere,	*/
/*			v1 and v2 vectors, and frame										*/
/*		27. Nov -- fixed bugs in b_val(), bperp_val(), n_val()					*/
/*		4. Dec 00 -- changed name of alpha variable to 'eta'					*/
/*		16. Dec 00 -- added edgeboundary_t() function, returning (x, y, z) loc	*/
/*			of the boundary for a given t-value parameterizing the arc between	*/
/*			(p + v1) and (p + v2).												*/
/*		8. Feb 01 -- added from3DtoV{1,2}_bperp_basis(), to give coordinates	*/
/*			of vector in R^3 as projected onto v1_bperp plane, in [v1, bperp]	*/
/*			coordinates (normalized basis)										*/
/*		17. May 02 -- added Push/PopAttrib() to GLrender()						*/
/********************************************************************************/

/********************************************************************************/
/* Class Diatom:  a medial atom base class, giving direction and displacement 	*/
/*	to paired BASOC regions at distance r from the Diatom, and having their		*/
/*	surface normals parallel to the Diatom-BASOC vector.						*/
/********************************************************************************/

/********************************************************************************/
/* NOTE on transformation of diatoms:											*/
/*    a Diatom can be translated and rotated normally, but all scaling must		*/
/*    be done uniformly---non-uniform stretchings change the medial structure,  */
/*    and can only be done by manipulation of the medial surface, not by		*/
/*    changing individual Diatoms.												*/
/********************************************************************************/

#include <iostream>
#include <stdio.h>
#include <math.h>

#include <stdlib.h>

#ifdef WIN32
#include <windows.h>
#endif
#ifdef __APPLE__
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#else
#include <GL/gl.h>
#include <GL/glu.h>
#endif

#define D_DIATOM
#define D_XFERLIST
#include "Shapedepend.h"

using namespace ThallCode;

typedef GLfloat Fastvec[3];
const Quat NULLROT(0.0, 0.0, 0.0, 1.0);


inline void loadfastvec(Fastvec fvec, DbVector3 *vec)
{
    fvec[0] = (float) vec->X();
    fvec[1] = (float) vec->Y();
    fvec[2] = (float) vec->Z();
}

const double EPSveclength = 1e-8;
/********************************************************************************/
/* constructor for Diatom with given location, radius, and perp. directions		*/
/*    note:  doesn't assume v1 and v2 are normalized, but can be.				*/
/********************************************************************************/
/*
Diatom::Diatom(const DbVector3& location, double radius,
			   const DbVector3& v1dir, const DbVector3& v2dir)
{
	p = location;
	r = radius;

	DbVector3 v1vec = v1dir.normalize();
	DbVector3 v2vec = v2dir.normalize();
	DbVector3 bperpvec = v1vec.cross(v2vec);

	if (fabs(bperpvec.length()) < EPSveclength) {
		// parallel sides, so just rotate zero quaternion to put z axis along v1dir,
		//   and let other axes be free in plane
		DbVector3 zaxis(0.0, 0.0, 1.0);
		DbVector3 rotaxis1 = zaxis.cross(v1vec);
		double rottheta1 = acos(zaxis.dot(v1vec));

		if (rottheta1 < EPSveclength)
			// zaxis is parallel to v1vec, so set null quaternion
			q = NULLROT;
		else if (rottheta1 - M_PI/2.0 < EPSveclength) {
			// zaxis is reversed from v1vec, so set neg. quaternion
			DbVector3 one (1.0, 0.0, 0.0);
			//q.rot(one, M_PI);
			q.setAxisAngle(Vector3D(one.X(), one.Y(), one.Z()), M_PI);
		}
		else
			//q.rot(rotaxis1, rottheta1);
			q.setAxisAngle(Vector3D(rotaxis1.X(), rotaxis1.Y(), rotaxis1.Z()), rottheta1);

		phi = 0;
	}
	else {
		DbVector3 bvecvec = v1vec + v2vec;
		bvecvec.selfnormalize();
		DbVector3 normvec = bvecvec.cross(bperpvec);
		normvec.selfnormalize();
		bperpvec.selfnormalize();

		// These vectors need to be normalized!!!!
		//q.q_from_framevecs(bvecvec, bperpvec, normvec);
		q.buildFromFrame(Vector3D(bvecvec.X(), bvecvec.Y(), bvecvec.Z()), 
			Vector3D(bperpvec.X(), bperpvec.Y(), bperpvec.Z()));

		phi = acos(v1vec.dot(normvec));
	}
	rho = DEFAULT_RHO;
	eta = DEFAULT_ETA;
	v_vals_computed = false;
	dtype = SOLO;

}
*/

/********************************************************************************/
/* v1_val() -- return DbVector3 giving direction and distance to the BASOC		*/
/*	computing it if necessary													*/
/********************************************************************************/
DbVector3& Diatom::v1_val()
{
    // if values valid, no computation is necessary, else
    if (!v_vals_computed) {

        // values invalid, so compute values for both v1 and v2
        //v1.set(sin(phi), 0, cos(phi));
		//v2.set(sin(phi), 0, -cos(phi));
		v2.set(cos(M_PI/2.0-phi), sin(M_PI/2.0-phi),  0);
		v1.set(cos(M_PI/2.0-phi), -sin(M_PI/2.0-phi), 0);

		//q.rotate_vec(v1);
			Vector3D v1Tmp(v1.X(), v1.Y(), v1.Z()); 
			q.rotateVector(v1Tmp);
			v1.set(v1Tmp.getX(), v1Tmp.getY(), v1Tmp.getZ());
		//q.rotate_vec(v2);
			Vector3D v2Tmp(v2.X(), v2.Y(), v2.Z()); 
			q.rotateVector(v2Tmp);
			v2.set(v2Tmp.getX(), v2Tmp.getY(), v2Tmp.getZ());

			//std::cout << "r: " << r << "; r1: " << r1 << "; r2: " << r2 << std::endl;
		v1 *= r1;
		v2 *= r2;
                v1 = u1 * r1;
                v2 = u2 * r2;              
                
                
        v_vals_computed = true;
    }
    return v1;
}

/********************************************************************************/
/* v2_val() -- return DbVector3 giving direction and distance to the BASOC		*/
/*	computing it if necessary													*/
/********************************************************************************/
DbVector3& Diatom::v2_val()
{
    // if values valid, no computation is necessary, else
    if (!v_vals_computed) {
        // values invalid, so compute values for both v1 and v2
        //v1.set(sin(phi), 0, cos(phi));
		//v2.set(sin(phi), 0, -cos(phi));
		v2.set(cos(M_PI/2.0-phi), sin(M_PI/2.0-phi),  0);
		v1.set(cos(M_PI/2.0-phi), -sin(M_PI/2.0-phi), 0);

		//q.rotate_vec(v1);
			Vector3D v1Tmp(v1.X(), v1.Y(), v1.Z()); 
			q.rotateVector(v1Tmp);
			v1.set(v1Tmp.getX(), v1Tmp.getY(), v1Tmp.getZ());
		//q.rotate_vec(v2);
			Vector3D v2Tmp(v2.X(), v2.Y(), v2.Z()); 
			q.rotateVector(v2Tmp);
			v2.set(v2Tmp.getX(), v2Tmp.getY(), v2Tmp.getZ());

		v1 *= r1;
		v2 *= r2;
                v1 = u1 * r1;
                v2 = u2 * r2;
        v_vals_computed = true;
    }
    return v2;
}

/********************************************************************************/
/* v_val() -- return DbVector3 giving direction and distance to the BASOC		*/
/*	computing it if necessary													*/
/* NOTE: v_val() returns v1 or v2 scaled by r to boundary locations				*/
/********************************************************************************/
DbVector3& Diatom::v_val(int v_select)
{
    // if values valid, no computation is necessary, else
    if (!v_vals_computed) {
        // values invalid, so compute values for both v1 and v2
        //v1.set(sin(phi), 0, cos(phi));
		//v2.set(sin(phi), 0, -cos(phi));
		v2.set(cos(M_PI/2-phi), sin(M_PI/2.0-phi),  0);
		v1.set(cos(M_PI/2-phi), -sin(M_PI/2.0-phi), 0);

		//q.rotate_vec(v1);
			Vector3D v1Tmp(v1.X(), v1.Y(), v1.Z()); 
			q.rotateVector(v1Tmp);
			v1.set(v1Tmp.getX(), v1Tmp.getY(), v1Tmp.getZ());
		//q.rotate_vec(v2);
			Vector3D v2Tmp(v2.X(), v2.Y(), v2.Z()); 
			q.rotateVector(v2Tmp);
			v2.set(v2Tmp.getX(), v2Tmp.getY(), v2Tmp.getZ());

		v1 *= r1;
		v2 *= r2;
                v1 = u1 * r1;
                v2 = u2 * r2;
        v_vals_computed = true;
    }
    if (v_select == 1)
        return v1;
    else
        return v2;
}

// return the bisector vector between v1 and v2, the scale-gradient direction
//   of the medial sheet in its tangent plane
//    DbVector3 bvector();

/********************************************************************************/
/* bvector() -- return the bisector vector between v1 and v2, the scale-gradient*/
/*	direction of the medial sheet in its tangent plane. 						*/
/*	return NULLVEC if phi is 0.													*/
/* MODIFICATION:  now return quaternion X-axis direction even for phi=0.		*/
/* NOTE: bvector is scaled to position of involute along bvec direction			*/
/* MODIFICATION:  scales by eta (= 1 for internal atoms, varies for boundary)	*/
/********************************************************************************/
DbVector3 Diatom::bvector()
{
    DbVector3 b(1.0, 0.0, 0.0);

    // b is the x axis of the quaternion frame, scaled by r
    b *= eta*r;

    //q.rotate_vec(b);
			Vector3D bTmp(b.X(), b.Y(), b.Z()); 
			q.rotateVector(bTmp);
			b.set(bTmp.getX(), bTmp.getY(), bTmp.getZ());


    return b;
}

/********************************************************************************/
/* b_val(), n_val(), bperp_val() -- return coordinate frame vector for Diatom	*/
/********************************************************************************/
DbVector3 Diatom::b_val()
{
	DbVector3 b(1.0, 0.0, 0.0);

	// b is the x axis of the quaternion frame
    //q.rotate_vec(b);
			Vector3D bTmp(b.X(), b.Y(), b.Z()); 
			q.rotateVector(bTmp);
			b.set(bTmp.getX(), bTmp.getY(), bTmp.getZ());

	return b;
}

DbVector3 Diatom::n_val()
{
	DbVector3 n(0.0, 1.0, 0.0);

	// n is the z axis of the quaternion frame
	//q.rotate_vec(n);
			Vector3D nTmp(n.X(), n.Y(), n.Z()); 
			q.rotateVector(nTmp);
			n.set(nTmp.getX(), nTmp.getY(), nTmp.getZ());

	return n;
}

DbVector3 Diatom::bperp_val()
{
	DbVector3 bperp(0.0, 0.0, 1.0);

	// bperp is the y axis of the quaternion frame
	//q.rotate_vec(bperp);
			Vector3D bpTmp(bperp.X(), bperp.Y(), bperp.Z()); 
			q.rotateVector(bpTmp);
			bperp.set(bpTmp.getX(), bpTmp.getY(), bpTmp.getZ());

	return bperp;
}


/****************************************************************************/
/* Functions to create Diatoms from boundary configurations				    */
/****************************************************************************/

/********************************************************************************/
/* sets diatom to have given location, radius, and perp. directions				*/
/*    note:  doesn't assume v1 and v2 are normalized, but can be.				*/
/********************************************************************************/
/*
void Diatom::diatom_from_involutes(const DbVector3& location, double radius,
								   const DbVector3& v1dir, const DbVector3& v2dir)
{
	p = location;
	r = radius;

	DbVector3 v1vec = v1dir.normalize();
	DbVector3 v2vec = v2dir.normalize();
	DbVector3 bperpvec = v1vec.cross(v2vec);

	if (fabs(bperpvec.length()) < EPSveclength) {
		// parallel sides, so just rotate zero quaternion to put z axis along v1dir,
		//   and let other axes be free in plane
		DbVector3 zaxis(0.0, 0.0, 1.0);
		DbVector3 rotaxis1 = zaxis.cross(v1vec);
		double rottheta1 = acos(zaxis.dot(v1vec));

		if (rottheta1 < EPSveclength)
			// zaxis is parallel to v1vec, so set null quaternion
			q = NULLROT;
		else if (rottheta1 - M_PI/2.0 < EPSveclength) {
			// zaxis is reversed from v1vec, so set neg. quaternion
			DbVector3 one (1.0, 0.0, 0.0);

			//q.rot(one, M_PI);
			q.setAxisAngle(Vector3D(one.X(), one.Y(), one.Z()), M_PI);
		}
		else
			//q.rot(rotaxis1, rottheta1);
			q.setAxisAngle(Vector3D(rotaxis1.X(), rotaxis1.Y(), rotaxis1.Z()), rottheta1);

		phi = 0;
	}
	else {
		DbVector3 bvecvec = v1vec + v2vec;
		bvecvec.selfnormalize();
		DbVector3 normvec = bvecvec.cross(bperpvec);
		normvec.selfnormalize();
		bperpvec.selfnormalize();

		// These vectors need to be normalized!!!!
		//q.q_from_framevecs(bvecvec, bperpvec, normvec);
		q.buildFromFrame(Vector3D(bvecvec.X(), bvecvec.Y(), bvecvec.Z()),
			Vector3D(bperpvec.X(), bperpvec.Y(), bperpvec.Z()));

		phi = acos(v1vec.dot(normvec));
	}
	rho = DEFAULT_RHO;
	eta = DEFAULT_ETA;
	v_vals_computed = false;
	dtype = SOLO;
}

// Given two boundary points and their normals, compute an "average" medial
//   atom by intersecting the normals with the equidistant plane and finding
//   the average point on this plane.  Can include a measure of normal deviation
//   at A and B, if needed---NOT DONE YET.
// PRECONDITION:  normA and normB are normalized, and are assumed outward pointing
void Diatom::averageDiatom(const DbVector3& ptA, const DbVector3& normA,
						   const DbVector3& ptB, const DbVector3& normB)
{
	DbVector3 midpoint, midnormal, normAintercept, normBintercept;
	double middist;

	midpoint = ptA + ptB;
	midpoint *= 0.5;

	midnormal = midpoint - ptA;
	middist = midnormal.length();
	midnormal.selfnormalize();

	// use similar triangles to find intercepts
	// 
	double adjacentlength = normA.dot(-midnormal);  // cos(angle) = length ratio,
	// nvec projected onto midnormal
	normAintercept = ptA - normA*(middist/adjacentlength);

	adjacentlength = normB.dot(midnormal);
	normBintercept = ptB - normB*(middist/adjacentlength);

	DbVector3 newmedialpoint(normAintercept + normBintercept);
	newmedialpoint *= 0.5;
    DbVector3 newv1(ptA - newmedialpoint);
	DbVector3 newv2(ptB - newmedialpoint);
	double newradius = newv1.length();

	diatom_from_involutes(newmedialpoint, newradius, newv1, newv2);
	GLrender(RAK_atom | RAK_v1v2 | RAK_frame);
}

// as above, but instead, go out distance rval inward from pt{A,B} along their
//    normals, and compute the intersection of the segment connecting these two
//    interior points with the bisecting plane.  This avoids problems with 
//    normal-to-medial-plane intersections that go to dist(infinity) as the
//    boundary normals near coplanarity with the bisecting plane.
// NO, that doesn't work either.  Project midpoint of their inverse-involutes
//    onto the plane.
void Diatom::averageDiatom2(const DbVector3& ptA, const DbVector3& normA, double rA,
							const DbVector3& ptB, const DbVector3& normB, double rB)
{
	DbVector3 midpoint, midnormal, normAintercept, normBintercept;
	double middist;
	DbVector3 newmedial;

	midpoint = ptA + ptB;
	midpoint *= 0.5;

	midnormal = ptA - midpoint;
	middist = midnormal.length();
	midnormal.selfnormalize();

	normAintercept = ptA - rA*normA;
	normBintercept = ptB - rB*normB;

	DbVector3 newmidinvolute = normAintercept + normBintercept;
	newmidinvolute *= 0.5;

	double projNMI = midnormal.dot(newmidinvolute - midpoint);

	newmedial = newmidinvolute - projNMI*midnormal;

    DbVector3 newv1(ptA - newmedial);
	DbVector3 newv2(ptB - newmedial);
	double newradius = newv1.length();

	diatom_from_involutes(newmedial, newradius, newv1, newv2);
}
*/

/*
/********************************************************************************
/* change orientation by rotating about q_rot to New_q = qrot*q			
/******************************************************************************** /
void Diatom::rot_q(Quat qrot, Diatom& newDA)
{
	newDA.p = myself.p;
	newDA.q = qrot*myself.q;
	newDA.phi = myself.phi;
	newDA.r = myself.r;
	newDA.rho = myself.rho;
	v_vals_computed = false;
}

/********************************************************************************
/* as above, but use change-of-basis matrix matrot instead of quaternion	*
/******************************************************************************** /
void Diatom::rot_q(DbMatrix3& matrot, Diatom& newDA)
{
	Quat newq;

	newq.makequat(matrot);
	newDA.p = myself.p;
	newDA.q = newq*myself.q;
	newDA.phi = myself.phi;
	newDA.r = myself.r;
	newDA.rho = myself.rho;
	v_vals_computed = false;
}


/********************************************************************************
/* translate and rotate Diatom appropriately for a rotation of qrot		*
/*   about the world coordinate origin, returning self				*
/******************************************************************************** /
Diatom& Diatom::selfrot_Diatom(Quat q_rot)
{
    q = q_rot*q;
    v_vals_computed = false;

    return myself;
}

/********************************************************************************
/* translate p by the give xlation vector, returning self			*
/******************************************************************************** /
Diatom& Diatom::selfxlate_p(DbVector3 xlation)
{
    //p += xlation;
    return myself;
}

/********************************************************************************
/* as above, but use triple <x, y, z>						*
/******************************************************************************** /
Diatom& Diatom::selfxlate_p(double x, double y, double z)
{
    double newx = p.X() + x;
    double newy = p.Y() + y;
    double newz = p.Z() + z;
    p.set(newx, newy, newz);

	return myself;
}

/********************************************************************************
/* scale r by the multipler delta_r, returning self				*
/******************************************************************************** /
Diatom& Diatom::selfzoom_r(double delta_r)
{
    if (delta_r == 0.0)
		fprintf(stderr, "WARNING:  setting Diatom radius to 0\n");
    r *= delta_r;

    v_vals_computed = false;

	return myself;

}

/********************************************************************************
/* rotate right-handed about the v1 vector					*
/******************************************************************************** /
void Diatom::rot_v1(double theta, Diatom& newDA)
{
	fprintf(stderr, "ERROR: rot_v1 not implemented\n");
}

/********************************************************************************
/* rotate right-handed about the v2 vector					*
/******************************************************************************** /
void Diatom::rot_v2(double theta, Diatom& newDA)
{
	fprintf(stderr, "ERROR: rot_v2 not implemented\n");
}

/********************************************************************************
/* translate and rotate Diatom appropriately for a rotation of qrot		*
/*   about the world coordinate origin.						*
/******************************************************************************** /
void Diatom::rot_Diatom(Quat q_rot, Diatom& newDA)
{
	fprintf(stderr, "ERROR: rot_Diatom not implemented\n");
}

/********************************************************************************
/* as above, but use change-of-basis matrix rather than quaternion		*
/******************************************************************************** /
void Diatom::rot_Diatom(DbMatrix3& matrot, Diatom& newDA)
{
	fprintf(stderr, "ERROR: rot_Diatom not implemented\n");
}

/********************************************************************************
/* translate p by the given xlation vector					*
/******************************************************************************** /
void Diatom::xlate_p(DbVector3 xlation, Diatom& newDA)
{
	fprintf(stderr, "ERROR: xlate_p not implemented\n");
}

/********************************************************************************
/* translate p by the <x, y, z> triple						*
/******************************************************************************** /
void Diatom::xlate_p(double x, double y, double z, Diatom& newDA)
{
	fprintf(stderr, "ERROR: xlate_p not implemented\n");
}

/********************************************************************************
/* scale r by the multipler delta_r						*
/******************************************************************************** /
void Diatom::zoom_r(double delta_r, Diatom& newDA)
{
	fprintf(stderr, "ERROR: zoom_r not implemented\n");
}

/********************************************************************************
/* change orientation by rotating about q_rot to New_q = qrot*q, returning self	*
/******************************************************************************** /
Diatom& Diatom::selfrot_q(Quat qrot)
{
	fprintf(stderr, "ERROR: selfrot_q not implemented\n");
    return myself;
}

/********************************************************************************
/* as above, but use change-of-basis matrix matrot instead of quaternion	*
/******************************************************************************** /
Diatom& Diatom::selfrot_q(DbMatrix3& matrot)
{
	fprintf(stderr, "ERROR: selfrot_q not implemented\n");
    return myself;
}

/********************************************************************************
/* rotate right-handed about the v1 direction					*
/******************************************************************************** /
Diatom& Diatom::selfrot_v1(double theta)
{
	fprintf(stderr, "ERROR: selfrot_v1 not implemented\n");
    return myself;
}

/********************************************************************************
/* rotate right-handed about the v2 direction					*
/******************************************************************************** /
Diatom& Diatom::selfrot_v2(double theta)
{
	fprintf(stderr, "ERROR: selfrot_v2 not implemented\n");
    return myself;
}

/********************************************************************************
/* as above, but use change of basis matrix matrot				*
/******************************************************************************** /
Diatom& Diatom::selfrot_Diatom(DbMatrix3& matrot)
{
	fprintf(stderr, "ERROR: selfrot_Diatom(DbMatrix3& ) not implemented\n");
    return myself;
}

*/

/********************************************************************************/
/* from3DtoDframe_basis() -- return coordinates of a vector (x, y, z)			*/
/*		in [b_val, bperp_val, n_val] coordinates								*/
/********************************************************************************/
/*
DbVector3 Diatom::from3DtoDframe_basis(const DbVector3& vecR3)
{
	DbVector3 changedbasis;
	changedbasis.X() = vecR3.dot(b_val());
	changedbasis.Y() = vecR3.dot(bperp_val());
	changedbasis.Z() = vecR3.dot(n_val());

	return changedbasis;
}
*/

/********************************************************************************/
/* from3DtoV1_bperp_basis() -- return coordinates of a vector (x, y, z)			*/
/*    projected to plane [v1, bperp_val()] in normalized [v1, bperp] coordinates*/
/********************************************************************************/
/*
DbVector2 Diatom::from3DtoV1_bperp_basis(const DbVector3& vecR3)
{
	DbVector2 changedbasis;
	changedbasis.X() = vecR3.dot(v1_val().normalize());
	changedbasis.Y() = vecR3.dot(bperp_val());

	return changedbasis;
}
*/

/********************************************************************************/
/* from3DtoV2_bperp_basis() -- return coordinates of a vector (x, y, z)			*/
/*    projected to plane [v2, bperp_val()] in normalized [v2, bperp] coordinates*/
/********************************************************************************/
/*
DbVector2 Diatom::from3DtoV2_bperp_basis(const DbVector3& vecR3)
{
	DbVector2 changedbasis;
	changedbasis.X() = vecR3.dot(v2_val().normalize());
	changedbasis.Y() = vecR3.dot(bperp_val());

	return changedbasis;
}
*/

/********************************************************************************/
/* etafun() -- compute a multiplier for the etavec for a given theta and		*/
/*		a given rotation rotangle from v1 to v2 where angle(v1,v2) = 2*theta.	*/
/********************************************************************************/
static inline double etafun(double eta, double theta, double rotangle)
{
	// try a simple sinusoid---window a cosine between -pi and pi, shift and scale
	// So, rotangle goes from 0 to 2*theta, so subtract theta and multiply by pi/theta

	double sinusend = ((rotangle - theta)*M_PI/theta);
	double multiplier = (cos(sinusend) + 1.0) / 2.0;
	return 1.0 + multiplier*(eta - 1.0);
}
/********************************************************************************/
/* edgeboundary_t() function, returning (x, y, z) loc of the boundary for		*/
/*    a given t-value parameterizing the arc between (p + v1) and (p + v2).		*/
/********************************************************************************/
void Diatom::edgeboundary_t(double t, DbVector3& pnt, DbVector3& normal)
{
	Quat q_v1_b_v2;

	if (t > 1.0 || t < -1.0) {
		std::cerr << "Fatal Error:  in Diatom::edgeboundary_t() -- t value outside of [-1, 1].\n";
		exit(0);
	}

	DbVector3 bperp = myself.bperp_val();	// unit length perpendicular
	DbVector3 this_v1 = myself.v1_val();

	// amount of rotation to produce gridsubdivedge new samples between
	//   v1 or v2 and bvec
	double theta = myself.theta_val();

	// t goes from 1 (at v1) to 0 (at b*eta) to -1 (at v2)
	//   --the rotation from b+theta to b-theta is linearly parameterized by t.
	double rotangle = (1.0 - t) * theta;

	// rotate the v1 vector to the correct angle for given t
	//q_v1_b_v2.rot(bperp, rotangle);
	q_v1_b_v2.setAxisAngle(Vector3D(bperp.X(), bperp.Y(), bperp.Z()), rotangle);
	//q_v1_b_v2.rotate_vec(this_v1);
	Vector3D this_vTmp(this_v1.X(), this_v1.Y(), this_v1.Z()); 
	q_v1_b_v2.rotateVector(this_vTmp);
	this_v1.set(this_vTmp.getX(), this_vTmp.getY(), this_vTmp.getZ());

	double etamult = etafun(eta, theta, rotangle);

	normal = this_v1*etamult;
	pnt = p + normal;
	normal.selfnormalize(); 
}

/********************************************************************************/
/* GLrender() -- OpenGL output function, to draw Diatom as tiny sphere (1/8 radius)	*/
/*    with v1 & v2 vectors and optional Diatom Frame vectors					*/
/********************************************************************************/
void Diatom::GLrender(const unsigned int draw_flag)
{
	DbVector3 this_p = myself.p_val();
	double this_r = myself.r_val();
	Fastvec pt1, pt2, pt3, pt4;

	DbVector3 this_v1;
	DbVector3 this_v2;
	DbVector3 N, B, BPERP;

	// Save current color
	glPushAttrib(GL_CURRENT_BIT);
	// First, draw atom
	if (draw_flag & RAK_atom) {

		glEnable(GL_LIGHTING);
		glPushMatrix();
		glTranslated(this_p.X(), this_p.Y(), this_p.Z());
		GLUquadric *qobj = gluNewQuadric();
		glColor3d(1.0, 1.0, 0.0);
		gluSphere(qobj, this_r/12.0, 8, 8);
		gluDeleteQuadric(qobj);
		glPopMatrix();
		glDisable(GL_LIGHTING);
	}
	if (draw_flag & RAK_v1v2) {
		this_v1 = myself.v1_val();
		this_v2 = myself.v2_val();

		this_v1 = this_p + this_v1;
		this_v2 = this_p + this_v2;

		loadfastvec(pt1, &this_p);
		loadfastvec(pt2, &this_v1);
		loadfastvec(pt3, &this_v2);

		glColor3d(0.0, 1.0, 1.0);
		glBegin(GL_LINES);
		glVertex3fv(pt1);
		glVertex3fv(pt2);
		glVertex3fv(pt1);
		glVertex3fv(pt3);
		glEnd();
	}
	if (draw_flag & RAK_frame) {
		N = myself.n_val()*this_r;
		B = myself.b_val()*this_r;
		BPERP = myself.bperp_val()*this_r;

		N = this_p + N;
		B = this_p + B;
		BPERP = this_p + BPERP;

		loadfastvec(pt1, &this_p);
		loadfastvec(pt2, &B);
		loadfastvec(pt3, &BPERP);
		loadfastvec(pt4, &N);

		glColor3d(1.0, 0.0, 0.0);
		glBegin(GL_LINES);
		glVertex3fv(pt1);
		glVertex3fv(pt2);
		glColor3d(0.0, 1.0, 0.0);
		glVertex3fv(pt1);
		glVertex3fv(pt3);
		glColor3d(0.0, 0.0, 1.0);
		glVertex3fv(pt1);
		glVertex3fv(pt4);
		glEnd();
	}

	glPopAttrib();

}

/********************************************************************************/
/* Diagnostic
/********************************************************************************/
void Diatom::printvals(char *message)
{
    if (message != NULL)
        std::cerr << message << '\n';
    p.printvals((char *)"p = ");
    //q.printvals((char *)"q = ");
	q.print();

    std::cerr << "phi = " << phi << ", r = " << r << ", rho = " << rho
		<< "\n, eta = " << eta << '\n';
    std::cerr << "v_vals_computed = " << v_vals_computed << '\n';
    if (v_vals_computed) {
        v1.printvals((char *)"v1 value");
		v2.printvals((char *)"v2 value");
    }
    switch (dtype) {
	case EDGE_M:
		std::cerr << "Diatom type is EDGE_M\n";
		break;
	case CORNER_M:
		std::cerr << "Diatom type is CORNER_M\n";
		break;
	case INTERNAL_M:
		std::cerr << "Diatom type is INTERNAL_M\n";
		break;
	case END_L:
		std::cerr << "Diatom type is END_L\n";
		break;
	case INTERNAL_L:
		std::cerr << "Diatom type is INTERNAL_L\n";
		break;
	case SOLO:
		std::cerr << "Diatom type is SOLO\n";
		break;
	default:
		std::cerr << "unknown Diatom type\n";
    }
}


