/********************************************************************************/
/*																				*/
/*  	File	:  Diatom.H														*/
/*																				*/
/*	Description:  class declaration for medial atom base class					*/
/*		Diatoms have as components:												*/
/*			a vector position p = (x,y,z)										*/
/*			a Quaternion orientation q = a + b*i + c*j + d*k					*/
/*			an angle phi giving the radian angle from the z axis				*/
/*			    in the positive x direction in the x-z plane					*/
/*			a radius r giving the distance from Diatom to BASOC					*/
/*			ratios rho1 & rho2 giving the conversion between r and				*/
/*			     sigma, the BASOC's tolerance std. dev.							*/
/*																				*/
/*	Project :  Seurat															*/
/*	Author  :  A. Thall															*/
/*	Date	:  4. February 1998													*/
/*																				*/
/*	Modifications:  															*/
/*		16. April -- moved typedef Diatomgrid here from Quadfig 				*/
/*		22. Sept 98 -- added EDGE_M, CORNER_M, INTERNAL_M, END_L,				*/
/*			INTERNAL_L, SOLO DiatomTYPE flag to diatom's field.					*/
/*			---probably better in Fig* class, but can change later.				*/
/*		22. Sept 98 -- added bvector() function									*/
/*		29. May 99 -- modified bvector() to give value even for phi=0,			*/
/*			based on quaternion frame x-axis direction.							*/
/*		19. May 2000 -- modified Diatom to have eta variable for b-vector		*/
/*			multiplier for edge and end-atoms.	For other atoms, eta = 1.0 	*/
/*		21. May -- added b_val(), bperp_val(), n_val() to get unit vectors		*/
/*			in cardinal directions of Diatom.									*/
/*		24. May -- added DEFAULT_RHO and DEFAULT_eta							*/
/*		30. May -- added Parentset field to Diatomgrid							*/
/*		 1. June -- class encapsulated Diatomgrid as an open struct (the same)	*/
/*			but with constructor (pointers = NULL) and destructors.				*/
/*		 2. June -- added set_theta() and theta_val(), == M_PI/2 - phi.			*/
/*		 4. June -- made Diatomgrid a separate class with own header file		*/
/*		 1. Aug -- made copy constructor & operator= take (const Diatom&)		*/
/*		30. Oct 00 -- added from3DtoDframe_basis() to give coordinates of a		*/
/*			vector in R^3 in [b, bperp, n] coordinates.  Useful for projection	*/
/*			onto tangent space of medial surface at Diatom, etc.				*/
/*		14. Nov -- added GLrender() routine, to draw Diatom as  tiny sphere,	*/
/*			v1 and v2 vectors, and frame										*/
/*		27. Nov -- fixed bugs in b_val(), bperp_val(), n_val()					*/
/*		 4. Dec 00 -- changed name of alpha variable to 'eta'					*/
/*		16. Dec 00 -- added edgeboundary_t() function, returning (x, y, z) loc	*/
/*			of the boundary for a given t-value parameterizing the arc between	*/
/*			(p + v1) and (p + v2).												*/
/*		8. Feb 01 -- added from3DtoV{1,2}_bperp_basis(), to give coordinates	*/
/*			of vector in R^3 as projected onto v1_bperp plane, in [v1, bperp]	*/
/*			coordinates (normalized basis)										*/
/*		23. Sept 01 -- wrote new constructor Diatom(pnt, rad, v1, v2)			*/
/*			note:  doesn't assume v1 and v2 are normalized, but can be.			*/
/*																				*/
/*		25. June 2011 - Dibyendu												*/
/*			Note: added members to incoporate boundary calculation for s-reps	*/
/*																				*/
/*																				*/
/********************************************************************************/

// Flags for Diatom rendering
const unsigned int RAK_atom = 1;
const unsigned int RAK_v1v2 = 2;
const unsigned int RAK_frame = 4;

const DbVector3 NULLVEC	= DbVector3(0.0, 0.0, 0.0);
const double DEFAULT_RHO = 0.25;
const double DEFAULT_ETA = 1.0;

// Define a struct to hold orientation, scaling, and position of EVERYTHING
typedef struct {
	DbVector3 Objpos;
      Quat Objorient;
   	double Objscale;
} spatialstate;

typedef enum { EDGE_M, CORNER_M, INTERNAL_M, END_L, INTERNAL_L, SOLO } DiatomTYPE;

/********************************************************************************/
/* Class Diatom:  a medial atom base class, giving direction and displacement 	*/
/*	to paired BASOC regions at distance r from the Diatom, and having their */
/*	surface normals parallel to the Diatom-BASOC vector.			*/
/********************************************************************************/

/********************************************************************************/
/* NOTE on transformation of diatoms:						*/
/*    a Diatom can be translated and rotated normally, but all scaling must	*/
/*    be done uniformly---non-uniform stretchings change the medial structure,  */
/*    and can only be done by manipulation of the medial surface, not by	*/
/*    changing individual Diatoms.						*/
/********************************************************************************/

class Diatom 
{
    DbVector3 p; // location of Diatom
    Quat q;    	// orientation of Diatom
    double phi;	// the radian angle from the z-axis in the positive x direction
		//   under orientation q
    double r, r1, r2;	// the radius of the perpendicular distance to the BASOC from p
    double rho;	// ratio giving conversion factor between r and sigma--the std.dev
                //   of the BASOC tolerance region
	double eta; // a multiplier for end and edge-Diatoms, so that distance to
				  //   medially-implied boundary is eta*r in b-vector direction
    DbVector3 v1, v2;	// (non-unit) vectors in the perpendicular directions from p 
                        //   to BASOC
    // Jared 2012/10/31: fixing problems with spoke ends not matching boundary mesh
    DbVector3 u1, u2, uEnd;
    double rEnd;
    
    bool v_vals_computed;	// flag indicating above vectors computed
                                //   from q, phi and r
    DiatomTYPE dtype;  // is the diatom an EDGE_M, CORNER_M, or INTERNAL_M atom
  		       //   in a mesh, is it an END_L or INTERNAL_L in a linear
  		       //   figure, or is it a SOLO diatom?
	//double r1, r2;

public:
    // default constructor -- static storage only
    Diatom() { ; }

    // copy constructor
    Diatom(const Diatom& DA);

    // destructor
    ~Diatom() { ; }

    // copy operator
    Diatom& operator=(const Diatom& DA);

    // constructor for Diatom with given location, radius, and perp. directions
    //Diatom(const DbVector3& location, double radius, const DbVector3& v1dir, const DbVector3& v2dir);

    /****************************************************************************/
    /* data access functions							*/
    /****************************************************************************/

    DbVector3& p_val() 	{ return p; }
    Quat& q_val() 	{ return q; }
    double& phi_val() 	{ return phi; }
    double& r_val() 	{ return r; }
    double& rho_val() 	{ return rho; }
	double& eta_val() { return eta; }
	double& r1_val() {return r1;}
	double& r2_val() {return r2;}
        double& rEnd_val() {return rEnd;}
        DbVector3& u1_val() { return u1; }
        DbVector3& u2_val() { return u2; }
        DbVector3& uEnd_val() { return uEnd; }
        

	// obviously, doesn't return a reference
	double theta_val() { return M_PI/2.0 - phi; }

    void set_p(DbVector3 vec) { p = vec; }
    void set_q(Quat q1) { q = q1; v_vals_computed = false; }
    void set_phi(double val) { phi = val; v_vals_computed = false; }
	void set_theta(double val) { phi = M_PI/2.0 - val; v_vals_computed = false; }
    void set_r(double val) { r = val; v_vals_computed = false; }
	void set_r1(double val) { r1 = val; v_vals_computed = false; }
	void set_r2(double val) { r2 = val; v_vals_computed = false; }
        void set_rEnd(double val) { rEnd = val; v_vals_computed = false; }
    void set_rho(double val) { rho = val; }
	void set_eta(double val) { eta = val; }
        void set_u1(DbVector3 vec1) { u1 = vec1; v_vals_computed = false; }
        void set_u2(DbVector3 vec2) { u2 = vec2; v_vals_computed = false; }
        void set_uEnd(DbVector3 vec) { uEnd = vec; v_vals_computed = false; }

    void Diatomtype(DiatomTYPE whatkind) { dtype = whatkind; }
    DiatomTYPE Diatomtype() { return dtype; }

	// v1 and v2 returned are NOT unit length, but give direction and distance
	//   to medial involutes on boundary
    DbVector3& v1_val();
    DbVector3& v2_val();
    DbVector3& v_val(int v_select);

    // return the bisector vector between v1 and v2, the scale-gradient direction
    //   of the medial sheet in its tangent plane; return NULLVEC if phi is 0.
    // ---MODIFICATION:  since just using x-axis of quaternion frame, return this
    //   axis even if phi=0
	// NOTE: bvector NOT unit length, but gives direction and distance
	//   to medial involutes on boundary, including eta multiplier
    DbVector3 bvector();

	// return cardinal frames of Diatom
	DbVector3 b_val();
	DbVector3 bperp_val();
	DbVector3 n_val();

	// return grad(r) = - b_val * sin(phi);
	DbVector3 grad_r() { return b_val()*(-sin(phi_val())); }

	// return coordinates of a vector (x, y, z) in [b_val, bperp_val, n_val] coordinates
	//DbVector3 from3DtoDframe_basis(const DbVector3& vecR3);

	// return coordinates of a vector (x, y, z) in [v{1,2}_val().normalize(), bperp_val]
	//    coordinates (projection onto that plane)
	//DbVector2 from3DtoV1_bperp_basis(const DbVector3& vecR3);
	//DbVector2 from3DtoV2_bperp_basis(const DbVector3& vecR3);

    void set_vflag() { v_vals_computed = true; }
    void clear_vflag() { v_vals_computed = false; }

	/****************************************************************************/
    /* Functions to create Diatoms from boundary configurations				    */
    /****************************************************************************/

	/********************************************************************************/
	/* sets diatom to have given location, radius, and perp. directions				*/
	/*    note:  doesn't assume v1 and v2 are normalized, but can be.				*/
	/*    Same as constructor above.												*/
	/********************************************************************************/
	//void diatom_from_involutes(const DbVector3& location, double radius,
	//								   const DbVector3& v1dir, const DbVector3& v2dir);

	// Given two boundary points and their normals, compute an "average" medial
	//   atom by intersecting the normals with the equidistant plane and finding
	//   the average point on this plane.  Can include a measure of normal deviation
	//   at A and B, if needed---NOT DONE YET.
    // PRECONDITION:  normA and normB are normalized, and are assumed outward pointing
	//void averageDiatom(const DbVector3& ptA, const DbVector3& normA,
	//	                               const DbVector3& ptB, const DbVector3& normB);

	// as above, but instead, go out distance rval inward from pt{A,B} along their
	//    normals, and compute the intersection of the segment connecting these two
	//    interior points with the bisecting plane.  This avoids problems with 
	//    normal-to-medial-plane intersections that go to dist(infinity) as the
	//    boundary normals near coplanarity with the bisecting plane.
	//void averageDiatom2(const DbVector3& ptA, const DbVector3& normA, double rA,
	//	                const DbVector3& ptB, const DbVector3& normB, double rB);

    /****************************************************************************/
    /* Functions to translate, rotate, and scale Diatom, returning new Diatom   */
    /****************************************************************************/

    // change orientation by rotating about q_rot to New_q = qrot*q
//    void rot_q(Quat qrot, Diatom& newDA);
    // as above, but use change-of-basis matrix matrot instead of Quaternion
//    void rot_q(DbMatrix3& matrot, Diatom& newDA); 

    // rotate the Diatom theta radians in a right-handed direction about the
    //   v1 or v2 axis, returning the new Diatom
//    void rot_v1(double theta, Diatom& newDA);
//    void rot_v2(double theta, Diatom& newDA);

    // translate and rotate Diatom appropriately for a rotation of qrot or matrot
    //   about the world coordinate origin.
//    void rot_Diatom(Quat q_rot, Diatom& newDA); 
//    void rot_Diatom(DbMatrix3& matrot, Diatom& newDA);

    // translate p by the give xlation vector or <x, y, z> triple
//	void xlate_p(DbVector3 xlation, Diatom& newDA);
//	void xlate_p(double x, double y, double z, Diatom& newDA);	

    // scale r by the delta_r multiplier
//    void zoom_r(double delta_r, Diatom& newDA);

    /********************************************************************/
    /* Functions to translate, rotate, and scale Diatom, modifying self */
    /********************************************************************/


    // change orientation by rotating about q_rot to New_q = qrot*q
//    Diatom& selfrot_q(Quat qrot);
    // as above, but use change-of-basis matrix matrot instead of Quaternion
//    Diatom& selfrot_q(DbMatrix3& matrot); 

    // rotate the Diatom theta radians in a right-handed direction about the
    //   v1 or v2 axis
//    Diatom& selfrot_v1(double theta);
//    Diatom& selfrot_v2(double theta);

    // translate and rotate Diatom appropriately for a rotation of qrot or matrot
    //   about the world coordinate origin.
//    Diatom& selfrot_Diatom(Quat q_rot);			//////
//    Diatom& selfrot_Diatom(DbMatrix3& matrot);	//////

    // translate p by the give xlation vector or <x, y, z> triple
//    Diatom& selfxlate_p(DbVector3 xlation);		//////
//    Diatom& selfxlate_p(double x, double y, double z);	

    // scale r by the delta_r multiplier
//    Diatom& selfzoom_r(double delta_r);			//////

    /********************************************************************/
    /* output of values to stderr, along with message header		*/
    /********************************************************************/
    void printvals(char *message = NULL);

	// GLrender() -- OpenGL output function, to draw Diatom as tiny sphere (1/8 radius)
	//    with v1 & v2 vectors and Diatom Frame vectors
	void GLrender(const unsigned int draw_flag = RAK_atom);

	// edgeboundary_t() function, returning (x, y, z) loc of the boundary for 
	//   a given t-value parameterizing the arc between	(p + v1) and (p + v2).
	void edgeboundary_t(double t, DbVector3& pnt, DbVector3& normal);
};

/********************************************************************************/
/* copy constructor																*/
/********************************************************************************/
inline Diatom::Diatom(const Diatom& DA)
{
    p = DA.p;
    q = DA.q;
    phi = DA.phi;
    r = DA.r;
	r1 = DA.r1;
	r2 = DA.r2;
    rho = DA.rho;
	eta = DA.eta;
    v1 = DA.v1;
    v2 = DA.v2;
    v_vals_computed = DA.v_vals_computed;
    dtype = DA.dtype;
    u1 = DA.u1;
    u2 = DA.u2;
    uEnd = DA.uEnd;
    rEnd = DA.rEnd;
}

/********************************************************************************/
/* copy operator								*/
/********************************************************************************/
inline Diatom& Diatom::operator=(const Diatom& DA)
{
    p = DA.p;
    q = DA.q;
    phi = DA.phi;
    r = DA.r;
	r1 = DA.r1;
	r2 = DA.r2;
    rho = DA.rho;
	eta = DA.eta;
    v1 = DA.v1;
    v2 = DA.v2;
    v_vals_computed = DA.v_vals_computed;
    dtype = DA.dtype;
    u1 = DA.u1;
    u2 = DA.u2;
    uEnd = DA.uEnd;
    rEnd = DA.rEnd;
    return myself;
}

