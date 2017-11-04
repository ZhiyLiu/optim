/********************************************************************************/
/*																				*/
/*  	File	:  Zerofinder.cpp												*/
/*																				*/
/*	Description:  functions for class to find zeros of T(u, v) in the range		*/
/*					[0, 1]x[0, 1] for T as defined below.  This used to find	*/
/*					nearpoint (according to Phong-derived normal-vector) to		*/
/*					a triangle.													*/
/*				   See Fig. 4.3 of Thall dissertation for (u, v) parameterization*/
/*						description.											*/
/*																				*/
/*	Project :  Seurat															*/
/*																				*/
/*	Dependencies:  Bpoints and DbVectors only, so need Xferlist.h and LinAlg.h	*/
/*																				*/
/*	Author  :  A. Thall															*/
/*																				*/
/*	Date	:  6. August 2001													*/
/*	Modifications:																*/
/*				15. Sept01 -- modified to use Conjgrad code						*/
/*				12. Dec01 -- fixed bug, was returning TRUE when minima found	*/
/*					that was not = 0, true minima but not a ROOT as advertised  */
/********************************************************************************/

#define D_XFERLIST
#define D_LINALG
#define D_CONJGRAD
#define D_ZEROFINDER
#include "Shapedepend.h"
#include <limits.h>
#include <float.h>

using namespace ThallCode;

inline DbVector3 Zerofinder::pbar_fun(double u, double v)
{
	return A + (B - A)*u + (C - A)*v + (A - C)*u*v;
}

inline DbVector3 Zerofinder::pbar_fun_u(double u, double v)
{
	return (B - A) + (A - C)*v;
}

inline DbVector3 Zerofinder::pbar_fun_v(double u, double v)
{
	return (C - A) + (A - C)*u;
}

inline DbVector3 Zerofinder::nvec_fun(double u, double v)
{
	return nA + (nB - nA)*u + (nC - nA)*v + (nA - nC)*u*v;
}

inline DbVector3 Zerofinder::nvec_fun_u(double u, double v)
{
	return (nB - nA) + (nA - nC)*v;
}

inline DbVector3 Zerofinder::nvec_fun_v(double u, double v)
{
	return (nC - nA) + (nA - nC)*u;
}

/****************************************************************************/
/* T_fun() -- compute T(u,v), setting pbar, nvec and their derivs as well	*/
/****************************************************************************/
double Zerofinder::T_fun(double u, double v)
{
	DbVector3 pbar, nvec;
	pbar = pbar_fun(u, v);
	nvec = nvec_fun(u, v);

	DbVector3 mvec = xbar - pbar;
	double firstterm = mvec.dot(mvec) * nvec.dot(nvec);
	double subterm = mvec.dot(nvec);

	return firstterm - subterm*subterm;
}

// same as above, but assumes pbar and nvec already set
double Zerofinder::T_fun(const DbVector3& pbar, const DbVector3& nvec)
{
	DbVector3 mvec = xbar - pbar;
	double firstterm = mvec.dot(mvec) * nvec.dot(nvec);
	double subterm = mvec.dot(nvec);

	return firstterm - subterm*subterm;
}

/****************************************************************************/
/* T_fun_{u,v}() -- assuming pbar and nvec, pbar_{u,v}, nvec_{u,v} already  */
/*    computed, compute partial derivatives	T_{u, v}(u,v)					*/
/****************************************************************************/
inline double Zerofinder::T_fun_u(double u, double v)
{
	DbVector3 pbar, nvec, pbar_u, nvec_u;
	pbar = pbar_fun(u, v);
	nvec = nvec_fun(u, v);
	pbar_u = pbar_fun_u(u, v);
	nvec_u = nvec_fun_u(u, v);


	DbVector3 mvec = xbar - pbar;
	double firstterm = - (mvec.dot(pbar_u))*(nvec.dot(nvec));
	double secondterm = (mvec.dot(mvec))*(nvec.dot(nvec_u));
	double thirdterm = - (mvec.dot(nvec))*(mvec.dot(nvec_u) - nvec.dot(pbar_u));

	return 2.0*(firstterm + secondterm + thirdterm);
}

inline double Zerofinder::T_fun_v(double u, double v)
{
	DbVector3 pbar, nvec, pbar_v, nvec_v;
	pbar = pbar_fun(u, v);
	nvec = nvec_fun(u, v);
	pbar_v = pbar_fun_v(u, v);
	nvec_v = nvec_fun_v(u, v);

	DbVector3 mvec = xbar - pbar;
	double firstterm = - (mvec.dot(pbar_v))*(nvec.dot(nvec));
	double secondterm = (mvec.dot(mvec))*(nvec.dot(nvec_v));
	double thirdterm = - (mvec.dot(nvec))*(mvec.dot(nvec_v) - nvec.dot(pbar_v));

	return 2.0*(firstterm + secondterm + thirdterm);
}

int fun_calledtimes = 0;
int GRADfun_calledtimes = 0;

void Zerofinder::ComputeTandDerivs(double u, double v)
{
	// Compute pbar and derivs, and nvec and derivs, FAST
	DbVector3 BminusA = B - A;
	DbVector3 CminusA = C - A;
	DbVector3 nBminusnA = nB - nA;
	DbVector3 nCminusnA = nC - nA;
	DbVector3 pbarfast = A + BminusA*u + CminusA*(v - u*v);
	DbVector3 pbar_ufast = BminusA - CminusA*v;
    DbVector3 pbar_vfast = CminusA*(1.0 - u);
	DbVector3 nvecfast = nA + nBminusnA*u + nCminusnA*(v - u*v);

	DbVector3 nvec_ufast = nBminusnA - nCminusnA*v;
	DbVector3 nvec_vfast = nCminusnA*(1.0 - u);

	DbVector3 mvecfast = xbar - pbarfast;
	double MdotM = mvecfast.dot(mvecfast);
	double MdotN = mvecfast.dot(nvecfast);
	double NdotN = nvecfast.dot(nvecfast);

	Tfast = MdotM*NdotN - MdotN*MdotN;

	double firstterm = - NdotN*(mvecfast.dot(pbar_ufast));
	double secondterm = MdotM*(nvecfast.dot(nvec_ufast));
	double thirdterm = - MdotN*(mvecfast.dot(nvec_ufast) - nvecfast.dot(pbar_ufast));

	T_ufast = 2.0*(firstterm + secondterm + thirdterm);

	firstterm = - NdotN*(mvecfast.dot(pbar_vfast));
	secondterm = MdotM*(nvecfast.dot(nvec_vfast));
	thirdterm = - MdotN*(mvecfast.dot(nvec_vfast) - nvecfast.dot(pbar_vfast)); 

	T_vfast = 2.0*(firstterm + secondterm + thirdterm);
}

/********************************************************************************/
/* findroot() -- given a triangular patch ABC with normals, assuming order is	*/
/*	 counterclockwise, and given a point x_bar in space, find u, v, _pbar, _nvec*/
/*   of the point on the patch with Phong-interpolated normal passing through	*/
/*	 the point in space.  Does so by finding a zero of the sine function between*/
/*   the normal and vector to x_bar over all (u, v) on the triangle.			*/
/********************************************************************************/
/********************************************************************************/
/* NOTE:  for now, just do brute force for a minimum sine value, by dicing the	*/
/*   (u, v) parameter space and checking at discrete grid values.				*/
/********************************************************************************/
const double SQRT_inv_numsamples = 0.1;			// will test 25 (u, v) values


// encapsulated versions of T_fun() and gradT_fun() = [T_fun_u(), T_fun_v()]
	//   for passing to conjugate gradient routines
Zerofinder *thisfinder;

double funfunfun(DbVector2& uv)
{
	return thisfinder->T_fun(uv.X(), uv.Y());
}

void GRADfunfunfun(DbVector2& uv, DbVector2& gradval)
{
	gradval.X() = thisfinder->T_fun_u(uv.X(), uv.Y());
	gradval.Y() = thisfinder->T_fun_v(uv.X(), uv.Y());
}


// NOTE:  these FAST functions should only be called in tandem,
//   first, funfunfunFAST(), then GRADfunfunfunFAST().
//   The first calls ComputeTandDerivs() to compute values
//   as efficiently as possible.  We avoid testing (u, v) for
//   new values by simply insisting that funfunfunFAST() ALWAYS
//   be called first.
double FAST_funfunfun(DbVector2& uv)
{
//	fun_calledtimes++;
	thisfinder->ComputeTandDerivs(uv.X(), uv.Y());
	return thisfinder->Tfast;
}

void FAST_GRADfunfunfun(DbVector2& uv, DbVector2& gradval)
{
//	GRADfun_calledtimes++;
	gradval.X() = thisfinder->T_ufast;
	gradval.Y() = thisfinder->T_vfast;
}

const double FTOL = 1e-10;
const double TESTTOL = 1e-5;

bool Zerofinder::drawpointPandtval(double u, double v,
						  const DbVector3& x_bar,
						  const DbVector3& a, const DbVector3& b, const DbVector3& c,
						  const DbVector3& n_A, const DbVector3& n_B, const DbVector3& n_C)
{
	DbVector3 p_bar;
	double t_val;

	xbar = x_bar;
	A = a;
	B = b;
	C = c;
	nA = n_A;
	nB = n_B;
	nC = n_C;

	t_val = T_fun(u, v);
	p_bar = pbar_fun(u, v);

/*
	glBegin(GL_POINTS);
	glColor3d(t_logval, t_logval, t_logval);
	glVertex3d(p_bar.X(), p_bar.Y(), p_bar.Z());
	glEnd();
 */
	if (t_val < 0.000001) {

//		fprintf(stderr, "tval = %.15f < 0.000001 at P(%f, %f) =\n", t_val, u, v);
//		p_bar.printvals();


		glBegin(GL_POINTS);
		glColor3d(1.0 - t_val, 1.0 - t_val, 1.0 - t_val);
		glVertex3d(p_bar.X(), p_bar.Y(), p_bar.Z());
		glEnd();

		DbVector3 n_bar = nvec_fun(u, v);
		n_bar.selfnormalize();
		n_bar *= 0.99*(p_bar - x_bar).length();

		DbVector3 ptA = p_bar + n_bar;
		DbVector3 ptB = p_bar - n_bar;

		DbVector3 ptC = ((x_bar - ptA).length() < (x_bar - ptB).length()) ? ptA : ptB;

		glColor3d(0.0, 1.0, 0.0);
		glBegin(GL_LINES);
		glVertex3d(p_bar.X(), p_bar.Y(), p_bar.Z());
		glVertex3d(ptC.X(), ptC.Y(), ptC.Z());
		glEnd();

		return true;
	}
	else
		return false;
}

// if thispoint is on the positive side of triangle (A-(A+nA)-B or A-(B+nB)-B)
//     and likewise for the other 2 edges BC, CA, then thispoint
//     has a possible Phong nearpoint on the triangle, else return FALSE
bool possible_phongneartriangle(const DbVector3& thispoint,
		const DbVector3& aA, const DbVector3& bB, const DbVector3& cC,
		const DbVector3& anA, const DbVector3& bnB, const DbVector3& cnC)
{
	bool retval;
	if (thispoint.outsideface(aA, aA + anA, bB) == -1 && thispoint.outsideface(aA, bB + bnB, bB) == -1)
		retval = false;
	else if (thispoint.outsideface(cC, cC + cnC, aA) == -1 && thispoint.outsideface(cC, aA + anA, aA) == -1)
		retval = false;
	else if (thispoint.outsideface(bB, bB + bnB, cC) == -1 && thispoint.outsideface(bB, cC + cnC, cC) == -1)
		retval = false;
	else
		retval = true;

	return retval;
}

bool Zerofinder::findroot(double &u, double &v, DbVector3& p_bar, DbVector3& n_vec,
	  const DbVector3& x_bar, 
	  const DbVector3& a, const DbVector3& b, const DbVector3& c,
	  const DbVector3& n_A, const DbVector3& n_B, const DbVector3& n_C)
{
	A = a;
	B = b;
	C = c;
	nA = n_A;
	nB = n_B;
	nC = n_C;

	xbar = x_bar;

	if ((A - B).cross(A - C) .dot(nA) < 0.0) {
		if (!possible_phongneartriangle(xbar, B, A, C, nB, nA, nC))
			return false;
	}
	else if (!possible_phongneartriangle(xbar, A, B, C, nA, nB, nC))
		return false;

	DbVector2 uv(0.0, 0.0);
	int num_iterations = 0;
	double minval = DBL_MAX;

	thisfinder = this;

	fun_calledtimes = 0;
	GRADfun_calledtimes = 0;
	Conjgrad minimizer(FAST_funfunfun, FAST_GRADfunfunfun);
    minimizer.frprmn(uv, FTOL, &num_iterations, &minval);

	if (num_iterations >= minimizer.maxiterations()) {

//		fprintf(stderr, "didn't find zero of function after %d iterations\n", num_iterations);
		return false;
	}
	else {
		u = uv.X();
		v = uv.Y();
		if (u > 1.0 || u < 0.0 || v > 1.0 || v < 0.0) {
//			fprintf(stderr, "converged on values outside [0,1]x[0,1] after %d iterations\n", num_iterations);
//			if (u > 1.1 || u < -0.1 || v > 1.1 || v < -0.1) {
//				fprintf(stderr, "converged on values (%f, %f) outside [0,1]x[0,1] after %d iterations\n", u, v, num_iterations);

				return false;
			}
	//		fprintf(stderr, "outside of bounds, but very near, (u, v) = (%.15f, %.15f)\n", u, v);
/*			fprintf(stderr, ".");
			if (u > 1.0) u = 1.0;
			if (u < 0.0) u = 0.0;
			if (v > 1.0) v = 1.0;
			if (v < 0.0) v = 0.0;
		}
 */
		p_bar = pbar_fun(u, v);
		n_vec = nvec_fun(u, v);

		// NOW, been a problem with false roots found, simply functional minima.
		//   Return false if n_vec.dot(xbar - p_bar) != 0;
		//   ---alternate to return false if minval > some_eps, but not sure what limit to use
		if (fabs(n_vec.normalize().dot((xbar - p_bar).normalize())) < 1.0 - TESTTOL) {
//		if (fabs(minval) > TESTTOL) {
//			fprintf(stderr, "bad root found, just a minval = %f\n", minval);
			return false;
		}
/*		fprintf(stderr, "found zero of function after %d iterations\n", num_iterations);
		fprintf(stderr, "   T(%f, %f) = %f\n", u, v, T_fun(u, v));
		fprintf(stderr, "   grad_T() = [%f, %f]\n", T_fun_u(u, v), T_fun_v(u, v));
		fprintf(stderr, "fun was called %d times and gradfun was called %d times\n", fun_calledtimes,
			    GRADfun_calledtimes);
 */
//		fun_calledtimes = 0;
//		GRADfun_calledtimes = 0;
/*
		glBegin(GL_POINTS);
		glColor3d(1.0, 0.0, 0.0);
		glVertex3d(xbar.X(), xbar.Y(), xbar.Z());
		glEnd();

		glBegin(GL_POINTS);
		glColor3d(1.0, 1.0, 1.0);
		glVertex3d(p_bar.X(), p_bar.Y(), p_bar.Z());
		glEnd();

		DbVector3 n_bar = n_vec.normalize();
		n_bar *= 0.99*(p_bar - xbar).length();

		DbVector3 ptA = p_bar + n_bar;
		DbVector3 ptB = p_bar - n_bar;

		DbVector3 ptC = ((xbar - ptA).length() < (xbar - ptB).length()) ? ptA : ptB;

		glColor3d(0.0, 1.0, 0.0);
		glBegin(GL_LINES);
		glVertex3d(p_bar.X(), p_bar.Y(), p_bar.Z());
		glVertex3d(ptC.X(), ptC.Y(), ptC.Z());
		glEnd();
*/
		return true;
	}
}

bool Zerofinder::findroot2(double &u, double &v, DbVector3& p_bar, DbVector3& n_vec,
	  const DbVector3& x_bar, 
	  const DbVector3& a, const DbVector3& b, const DbVector3& c,
	  const DbVector3& n_A, const DbVector3& n_B, const DbVector3& n_C)
{
	A = a;
	B = b;
	C = c;
	nA = n_A;
	nB = n_B;
	nC = n_C;
	xbar = x_bar;

	if ((A - B).cross(A - C) .dot(nA) < 0.0) {
		if (!possible_phongneartriangle(xbar, B, A, C, nB, nA, nC))
			return false;
	}
	else if (!possible_phongneartriangle(xbar, A, B, C, nA, nB, nC))
		return false;

	int ucnt, vcnt;
	double t_val;
	double min_tval = DBL_MAX;
	double minu, minv;
	double testu, testv;

	const int NUMITER1 = 20;
	const int NUMITER2 = 10;
	const int NUMITER3 = 10;
	const double TVALTOL = 1.0e-5;

	const double FRAC1 = 1.0/NUMITER1;	
	const double FRAC2 = 1.0/NUMITER2;
	const double FRAC3 = 1.0/NUMITER3;

	for (ucnt = 0; ucnt <= NUMITER1; ucnt++) {

		testu = ucnt*FRAC1;

		for (vcnt = 0; vcnt <= NUMITER1; vcnt++) {	

			testv = vcnt*FRAC1;
			t_val = T_fun(testu, testv);

			if (t_val < min_tval) {
				min_tval = t_val;
				minu = testu;
				minv = testv;
			}
		}
	}

	double startu, stopu, startv, stopv;

	// if need to refine
	if (min_tval > TVALTOL) {
//		fprintf(stderr, "   level 1 refining %.15f\n", min_tval);

		startu = MIN(0.0, minu - FRAC1);
		stopu = MAX(1.0, minu + FRAC1);
		startv = MIN(0.0, minv - FRAC1);
		stopv = MAX(1.0, minv + FRAC1);

		for (ucnt = 0; ucnt <= NUMITER2; ucnt++) {

			// compute test u parameter
			testu = startu + (stopu - startu)*(ucnt*FRAC2);

			for (vcnt = 0; vcnt <= NUMITER2; vcnt++) {

				// compute test v parameter
				testv = startv + (stopv - startv)*(vcnt*FRAC2);

				t_val = T_fun(testu, testv);
				if (t_val < min_tval) {
					min_tval = t_val;
					minu = testu;
					minv = testv;
				}
			}
		}
//		fprintf(stderr, "   level 1 ended up %.15f\n", min_tval);
	}



	if (min_tval > TVALTOL) {

		startu = MIN(0.0, minu - FRAC2);
		stopu = MAX(1.0, minu + FRAC2);
		startv = MIN(0.0, minv - FRAC2);
		stopv = MAX(1.0, minv + FRAC2);

		for (ucnt = 0; ucnt <= NUMITER3; ucnt++) {

			// compute test u parameter
			testu = startu + (stopu - startu)*(ucnt*FRAC3);

			for (vcnt = 0; vcnt <= NUMITER3; vcnt++) {

				// compute test v parameter
				testv = startv + (stopv - startv)*(vcnt*FRAC3);

				t_val = T_fun(testu, testv);
				if (t_val < min_tval) {
					min_tval = t_val;
					minu = testu;
					minv = testv;
				}
			}
		}
//		fprintf(stderr, "   level 2 ended up %.15f\n", min_tval);
	}


/*	glBegin(GL_POINTS);
	glColor3d(1.0, 0.0, 0.0);
	glVertex3d(xbar.X(), xbar.Y(), xbar.Z());
	glEnd();
*/
	if (min_tval > TVALTOL) {
//		fprintf(stderr, "did NOT find zero of function---min was %.15f\n", min_tval);
//		fprintf(stderr, "tfun gradient = [%.15f, %.15f]\n", T_fun_u(minu, minv), T_fun_v(minu, minv));
		return false;
	}
	else {
//		fprintf(stderr, "hey DID find zero of function---min was %.15f\n", min_tval);
		u = minu;
		v = minv;
		p_bar = pbar_fun(u, v);
		n_vec = nvec_fun(u, v);

/*		
		glBegin(GL_POINTS);
		glColor3d(1.0 - t_val, 1.0 - t_val, 1.0 - t_val);
		glVertex3d(p_bar.X(), p_bar.Y(), p_bar.Z());
		glEnd();

		DbVector3 n_bar = n_vec.normalize();
		n_bar *= 0.99*(p_bar - x_bar).length();

		DbVector3 ptA = p_bar + n_bar;
		DbVector3 ptB = p_bar - n_bar;

		DbVector3 ptC = ((x_bar - ptA).length() < (x_bar - ptB).length()) ? ptA : ptB;

		glColor3d(0.0, 1.0, 0.0);
		glBegin(GL_LINES);
		glVertex3d(p_bar.X(), p_bar.Y(), p_bar.Z());
		glVertex3d(ptC.X(), ptC.Y(), ptC.Z());
		glEnd();
*/
		return true;
	}
}

bool Zerofinder::findroot3(double &u, double &v, DbVector3& p_bar, DbVector3& n_vec,
	  const DbVector3& x_bar, 
	  const DbVector3& a, const DbVector3& b, const DbVector3& c,
	  const DbVector3& n_A, const DbVector3& n_B, const DbVector3& n_C)
{
	A = a;
	B = b;
	C = c;
	nA = n_A;
	nB = n_B;
	nC = n_C;
	xbar = x_bar;

	int ucnt, vcnt;
	double t_val;
	double min_tval = DBL_MAX;
	double minu, minv;
	double testu, testv;

	const int NUMITER1 = 100;
	const double TVALTOL = 1.0e-4;

	const double FRAC1 = 1.0/NUMITER1;	

	for (ucnt = 0; ucnt <= NUMITER1; ucnt++) {

		testu = ucnt*FRAC1;

		for (vcnt = 0; vcnt <= NUMITER1; vcnt++) {	

			testv = vcnt*FRAC1;
			t_val = T_fun(testu, testv);

			if (t_val < min_tval) {
				min_tval = t_val;
				minu = testu;
				minv = testv;
			}
		}
	}



/*	fprintf(stderr, "\n");
	fprintf(stderr, "   I couldn't think of one clever way to stop this guy,\n");
	fprintf(stderr, "   so I just trusted to mindless violence.\n");
*/
	if (min_tval > TVALTOL) {
		return false;
	}
	else {
//		fprintf(stderr, "hey DID find zero of function---min was %.15f\n", min_tval);
		u = minu;
		v = minv;
		p_bar = pbar_fun(u, v);
		n_vec = nvec_fun(u, v);
/*
		glBegin(GL_POINTS);
		glColor3d(1.0, 0.0, 0.0);
		glVertex3d(xbar.X(), xbar.Y(), xbar.Z());
		glEnd();

		glBegin(GL_POINTS);
		glColor3d(1.0 - t_val, 1.0 - t_val, 1.0 - t_val);
		glVertex3d(p_bar.X(), p_bar.Y(), p_bar.Z());
		glEnd();

		DbVector3 n_bar = n_vec.normalize();
		n_bar *= 0.99*(p_bar - x_bar).length();

		DbVector3 ptA = p_bar + n_bar;
		DbVector3 ptB = p_bar - n_bar;

		DbVector3 ptC = ((x_bar - ptA).length() < (x_bar - ptB).length()) ? ptA : ptB;

		glColor3d(1.0, 0.0, 1.0);
		glBegin(GL_LINES);
		glVertex3d(p_bar.X(), p_bar.Y(), p_bar.Z());
		glVertex3d(ptC.X(), ptC.Y(), ptC.Z());
		glEnd();
*/ 
		return true;
	}
}

// Bail out of recursive findroot4() if more than 15 levels deep
const int MAXRECUR = 15;

// This uses a Loop style triangle subdivision (without averaging) to subdivide
//    a possible triangle recursively until a good root is found
// Assumes vertices in counterclockwise (right-hand) order about outward-pointed
//    normal
bool Zerofinder::findroot4(Bpoint *bvals, const DbVector3& x_bar, 
						   DbVector3& a, DbVector3& b, DbVector3& c,
						   DbVector3& n_A, DbVector3& n_B, DbVector3& n_C,
						   DbVector3& uvt_A, DbVector3& uvt_B, DbVector3& uvt_C,
						   double r_A, double r_B, double r_C)
{
	rd += 1;
	if (rd > MAXRECUR) {
//		fprintf(stderr, "Bailing out of level 1: recursion too deep at depth %d\n", rd);
		rd -= 1;
		return false;
	}

	xbar = x_bar;

	if (!possible_phongneartriangle(xbar, a, b, c, n_A, n_B, n_C)){
		rd -= 1;
		return false;
	}

	const double TVALTOL = 1.0e-7;

	double testval;
	if ((testval = T_fun(a, n_A)) <  TVALTOL) {

		bvals->pnt[0] = a.X();
		bvals->pnt[1] = a.Y();
		bvals->pnt[2] = a.Z();
		bvals->norm[0] = n_A.X();
		bvals->norm[1] = n_A.Y();
		bvals->norm[2] = n_A.Z();
		bvals->u = uvt_A.X();
		bvals->v = uvt_A.Y();
		bvals->t = uvt_A.Z();
		bvals->rad = r_A;
		rd -= 1;

//		fprintf(stderr, "(u, v, t) = (%f, %f, %f)\n", bvals->u, bvals->v, bvals->t);
		return true;
	}

	if ((testval = T_fun(b, n_B)) <  TVALTOL) {

		bvals->pnt[0] = b.X();
		bvals->pnt[1] = b.Y();
		bvals->pnt[2] = a.Z();
		bvals->norm[0] = n_B.X();
		bvals->norm[1] = n_B.Y();
		bvals->norm[2] = n_B.Z();
		bvals->u = uvt_B.X();
		bvals->v = uvt_B.Y();
		bvals->t = uvt_B.Z();
		bvals->rad = r_B;
		rd -= 1;

//		fprintf(stderr, "(u, v, t) = (%f, %f, %f)\n", bvals->u, bvals->v, bvals->t);
		return true;
	}

	if ((testval = T_fun(c, n_C)) < TVALTOL) {

		bvals->pnt[0] = c.X();
		bvals->pnt[1] = c.Y();
		bvals->pnt[2] = c.Z();
		bvals->norm[0] = n_C.X();
		bvals->norm[1] = n_C.Y();
		bvals->norm[2] = n_C.Z();
		bvals->u = uvt_C.X();
		bvals->v = uvt_C.Y();
		bvals->t = uvt_C.Z();
		bvals->rad = r_C;
		rd -= 1;

//		fprintf(stderr, "(u, v, t) = (%f, %f, %f)\n", bvals->u, bvals->v, bvals->t);
		return true;
	}

	// Split into 4 subtriangles and run recursively on each
	DbVector3 d = (b + c)*0.5;
	DbVector3 n_D = (n_B + n_C)*0.5;
	DbVector3 uvt_D = (uvt_B + uvt_C)*0.5;
	double r_D = (r_B + r_C)*0.5;
	DbVector3 e = (c + a)*0.5;
	DbVector3 n_E = (n_C + n_A)*0.5;
	DbVector3 uvt_E = (uvt_C + uvt_A)*0.5;
	double r_E = (r_C + r_A)*0.5;
	DbVector3 f = (a + b)*0.5;
	DbVector3 n_F = (n_A + n_B)*0.5;
	DbVector3 uvt_F = (uvt_A + uvt_B)*0.5;
	double r_F = (r_A + r_B)*0.5;

	// AFE, FBD, EDC, FDE
	if (findroot4(bvals, xbar, a, f, e, n_A, n_F, n_E, uvt_A, uvt_F, uvt_E, r_A, r_F, r_E)) {
		rd -= 1;
		return true;
	}
	else if (findroot4(bvals, xbar, f, b, d, n_F, n_B, n_D, uvt_F, uvt_B, uvt_D, r_F, r_B, r_D)) {
		rd -= 1;
		return true;
	}
	else if (findroot4(bvals, xbar, e, d, c, n_E, n_D, n_C, uvt_E, uvt_D, uvt_C, r_E, r_D, r_C)) {
		rd -= 1;
		return true;
	}
	else if (findroot4(bvals, xbar, f, d, e, n_F, n_D, n_E, uvt_F, uvt_D, uvt_E, r_F, r_D, r_E)) {
		rd -= 1;
		return true;
	}
	else {
		rd -= 1;
		return false;
	}
}

/*
bool Zerofinder::findroot4(double &u, double &v, DbVector3& p_bar, DbVector3& n_vec,
	  const DbVector3& x_bar, 
	  const DbVector3& a, const DbVector3& b, const DbVector3& c,
	  const DbVector3& n_A, const DbVector3& n_B, const DbVector3& n_C,
	  DbVector2& uv_A, const DbVector2& uv_B, const DbVector2& uv_C)
{

	rd += 1;
	if (rd > 12) {
		fprintf(stderr, "Bailing out of level 1: recursion too deep at depth %d\n", rd);
		rd -= 1;
		return false;
	}

	DbVector3 locA = a;
	DbVector3 locB = b;
	DbVector3 locC = c;
	DbVector3 locnA = n_A;
	DbVector3 locnB = n_B;
	DbVector3 locnC = n_C;

	xbar = x_bar;

	if (!possible_phongneartriangle(xbar, locA, locB, locC, locnA, locnB, locnC)){
		rd -= 1;
		return false;
	}

	const double TVALTOL = 1.0e-5;

	double testval;
	if ((testval = T_fun(locA, locnA)) <  TVALTOL) {

		u = uv_A.X();
		v = uv_A.Y();
		p_bar = locA;
		n_vec = locnA;
		rd -= 1;
		return true;
	}

	if ((testval = T_fun(locB, locnB)) <  TVALTOL) {

		u = 1.0;
		v = 0.0;
		p_bar = locB;
		n_vec = locnB;
		rd -= 1;
		return true;
	}

	if ((testval = T_fun(locC, locnC)) < TVALTOL) {

		u = 0.0;
		v = 1.0;
		p_bar = locC;
		n_vec = locnC;
		rd -= 1;
		return true;
	}

	// Split into 4 subtriangles and run recursively on each
	DbVector3 locD = (locB + locC)*0.5;
	DbVector3 locnD = (locnB + locnC)*0.5;
	DbVector3 locE = (locC + locA)*0.5;
	DbVector3 locnE = (locnC +locnA)*0.5;
	DbVector3 locF = (locA + locB)*0.5;
	DbVector3 locnF = (locnA + locnB)*0.5;

	// AFE, FBD, EDC, FDE
	if (findroot4(u, v, p_bar, n_vec, xbar, locA, locF, locE, locnA, locnF, locnE)) {
		rd -= 1;
		return true;
	}
	else if (findroot4(u, v, p_bar, n_vec, xbar, locF, locB, locD, locnF, locnB, locnD)) {
		rd -= 1;
		return true;
	}
	else if (findroot4(u, v, p_bar, n_vec, xbar, locE, locD, locC, locnE,locnD, locnC)) {
		rd -= 1;
		return true;
	}
	else if (findroot4(u, v, p_bar, n_vec, xbar, locF, locD, locE, locnF, locnD, locnE)) {
		rd -= 1;
		return true;
	}
	else {
		rd -= 1;
		return false;
	}
}
*/
// v = 0, u = [0, 1]
inline DbVector3 Zerofinder::pbar_funV0(double u, double v)
{
	return A + (B - A)*u;
}

inline DbVector3 Zerofinder::nvec_funV0(double u, double v)
{
	return nA + (nB - nA)*u;
}

// v = 1, u = [0, 1]
inline DbVector3 Zerofinder::pbar_funV1(double u, double v)
{
	return C + (B - C)*u;
}

inline DbVector3 Zerofinder::nvec_funV1(double u, double v)
{
	return nC + (nB - nC)*u;
}

// v = [0, 1], u = 0
inline DbVector3 Zerofinder::pbar_funU0(double u, double v)
{
	return A + (C - A)*v;
}

inline DbVector3 Zerofinder::nvec_funU0(double u, double v)
{
	return nA + (nC - nA)*v;
}



// v = 0, u = [0, 1], direction A to B
void Zerofinder::ComputeTandDerivs_V0(double u)
{
	// Compute pbar and derivs, and nvec and derivs, FAST
	DbVector3 BminusA = B - A;
	DbVector3 nBminusnA = nB - nA;
	DbVector3 pbarfast = A + BminusA*u;
	DbVector3 pbar_ufast = BminusA;
	DbVector3 nvecfast = nA + nBminusnA*u;
	DbVector3 nvec_ufast = nBminusnA;

	DbVector3 mvecfast = xbar - pbarfast;
	double MdotM = mvecfast.dot(mvecfast);
	double MdotN = mvecfast.dot(nvecfast);
	double NdotN = nvecfast.dot(nvecfast);

	Tfast = MdotM*NdotN - MdotN*MdotN;

	double firstterm = - NdotN*(mvecfast.dot(pbar_ufast));
	double secondterm = MdotM*(nvecfast.dot(nvec_ufast));
	double thirdterm = - MdotN*(mvecfast.dot(nvec_ufast) - nvecfast.dot(pbar_ufast));

	T_ufast = 2.0*(firstterm + secondterm + thirdterm);
}

double Zerofinder::T_V0(double u)
{
	DbVector3 BminusA = B - A;
	DbVector3 nBminusnA = nB - nA;
	DbVector3 pbarfast = A + BminusA*u;
	DbVector3 pbar_ufast = BminusA;
	DbVector3 nvecfast = nA + nBminusnA*u;
	DbVector3 nvec_ufast = nBminusnA;

	DbVector3 mvecfast = xbar - pbarfast;
	double MdotM = mvecfast.dot(mvecfast);
	double MdotN = mvecfast.dot(nvecfast);
	double NdotN = nvecfast.dot(nvecfast);

	return MdotM*NdotN - MdotN*MdotN;
}

double Zerofinder::dT_V0(double u)
{
	// Compute pbar and derivs, and nvec and derivs, FAST
	DbVector3 BminusA = B - A;
	DbVector3 nBminusnA = nB - nA;
	DbVector3 pbarfast = A + BminusA*u;
	DbVector3 pbar_ufast = BminusA;
	DbVector3 nvecfast = nA + nBminusnA*u;
	DbVector3 nvec_ufast = nBminusnA;

	DbVector3 mvecfast = xbar - pbarfast;
	double MdotM = mvecfast.dot(mvecfast);
	double MdotN = mvecfast.dot(nvecfast);
	double NdotN = nvecfast.dot(nvecfast);

	double firstterm = - NdotN*(mvecfast.dot(pbar_ufast));
	double secondterm = MdotM*(nvecfast.dot(nvec_ufast));
	double thirdterm = - MdotN*(mvecfast.dot(nvec_ufast) - nvecfast.dot(pbar_ufast));

	return 2.0*(firstterm + secondterm + thirdterm);
}

// v = 1, u = [0, 1], direction C to B
void Zerofinder::ComputeTandDerivs_V1(double u)
{
	// Compute pbar and derivs, and nvec and derivs, FAST
	DbVector3 BminusC = B - C;
	DbVector3 nBminusnC = nB - nC;
	DbVector3 pbarfast = C + BminusC*u;
	DbVector3 pbar_ufast = BminusC;
	DbVector3 nvecfast = nC + nBminusnC*u;

	DbVector3 nvec_ufast = nBminusnC;

	DbVector3 mvecfast = xbar - pbarfast;
	double MdotM = mvecfast.dot(mvecfast);
	double MdotN = mvecfast.dot(nvecfast);
	double NdotN = nvecfast.dot(nvecfast);

	Tfast = MdotM*NdotN - MdotN*MdotN;

	double firstterm = - NdotN*(mvecfast.dot(pbar_ufast));
	double secondterm = MdotM*(nvecfast.dot(nvec_ufast));
	double thirdterm = - MdotN*(mvecfast.dot(nvec_ufast) - nvecfast.dot(pbar_ufast));

	T_ufast = 2.0*(firstterm + secondterm + thirdterm);
}

double Zerofinder::T_V1(double u)
{
	// Compute pbar and derivs, and nvec and derivs, FAST
	DbVector3 BminusC = B - C;
	DbVector3 nBminusnC = nB - nC;
	DbVector3 pbarfast = C + BminusC*u;
	DbVector3 pbar_ufast = BminusC;
	DbVector3 nvecfast = nC + nBminusnC*u;

	DbVector3 nvec_ufast = nBminusnC;

	DbVector3 mvecfast = xbar - pbarfast;
	double MdotM = mvecfast.dot(mvecfast);
	double MdotN = mvecfast.dot(nvecfast);
	double NdotN = nvecfast.dot(nvecfast);

	return MdotM*NdotN - MdotN*MdotN;
}

double Zerofinder::dT_V1(double u)
{
	// Compute pbar and derivs, and nvec and derivs, FAST
	DbVector3 BminusC = B - C;
	DbVector3 nBminusnC = nB - nC;
	DbVector3 pbarfast = C + BminusC*u;
	DbVector3 pbar_ufast = BminusC;
	DbVector3 nvecfast = nC + nBminusnC*u;

	DbVector3 nvec_ufast = nBminusnC;

	DbVector3 mvecfast = xbar - pbarfast;
	double MdotM = mvecfast.dot(mvecfast);
	double MdotN = mvecfast.dot(nvecfast);
	double NdotN = nvecfast.dot(nvecfast);

	double firstterm = - NdotN*(mvecfast.dot(pbar_ufast));
	double secondterm = MdotM*(nvecfast.dot(nvec_ufast));
	double thirdterm = - MdotN*(mvecfast.dot(nvec_ufast) - nvecfast.dot(pbar_ufast));

	return 2.0*(firstterm + secondterm + thirdterm);
}

// v = [0, 1], u = 0, direction A to C
void Zerofinder::ComputeTandDerivs_U0(double v)
{
	// Compute pbar and derivs, and nvec and derivs, FAST
	DbVector3 CminusA = C - A;
	DbVector3 nCminusnA = nC - nA;
	DbVector3 pbarfast = A + CminusA*v;
    DbVector3 pbar_vfast = CminusA;
	DbVector3 nvecfast = nA + nCminusnA*v;

	DbVector3 nvec_vfast = nCminusnA;

	DbVector3 mvecfast = xbar - pbarfast;
	double MdotM = mvecfast.dot(mvecfast);
	double MdotN = mvecfast.dot(nvecfast);
	double NdotN = nvecfast.dot(nvecfast);

	Tfast = MdotM*NdotN - MdotN*MdotN;

	double firstterm = - NdotN*(mvecfast.dot(pbar_vfast));
	double secondterm = MdotM*(nvecfast.dot(nvec_vfast));
	double thirdterm = - MdotN*(mvecfast.dot(nvec_vfast) - nvecfast.dot(pbar_vfast)); 

	T_vfast = 2.0*(firstterm + secondterm + thirdterm);
}

double Zerofinder::T_U0(double v)
{
	// Compute pbar and derivs, and nvec and derivs, FAST
	DbVector3 CminusA = C - A;
	DbVector3 nCminusnA = nC - nA;
	DbVector3 pbarfast = A + CminusA*v;
    DbVector3 pbar_vfast = CminusA;
	DbVector3 nvecfast = nA + nCminusnA*v;

	DbVector3 nvec_vfast = nCminusnA;

	DbVector3 mvecfast = xbar - pbarfast;
	double MdotM = mvecfast.dot(mvecfast);
	double MdotN = mvecfast.dot(nvecfast);
	double NdotN = nvecfast.dot(nvecfast);

	return MdotM*NdotN - MdotN*MdotN;
}

double Zerofinder::dT_U0(double v)
{
	// Compute pbar and derivs, and nvec and derivs, FAST
	DbVector3 CminusA = C - A;
	DbVector3 nCminusnA = nC - nA;
	DbVector3 pbarfast = A + CminusA*v;
    DbVector3 pbar_vfast = CminusA;
	DbVector3 nvecfast = nA + nCminusnA*v;

	DbVector3 nvec_vfast = nCminusnA;

	DbVector3 mvecfast = xbar - pbarfast;
	double MdotM = mvecfast.dot(mvecfast);
	double MdotN = mvecfast.dot(nvecfast);
	double NdotN = nvecfast.dot(nvecfast);

	double firstterm = - NdotN*(mvecfast.dot(pbar_vfast));
	double secondterm = MdotM*(nvecfast.dot(nvec_vfast));
	double thirdterm = - MdotN*(mvecfast.dot(nvec_vfast) - nvecfast.dot(pbar_vfast)); 

	return 2.0*(firstterm + secondterm + thirdterm);
}

typedef enum {V0edge, V1edge, U0edge} edgeTYPE;
edgeTYPE whichedge;

const double LINMINTOL = 1.0e-7;
/********************************************************************************/
/* findroot_boundary() -- find a root along the boundary of the given triangle,	*/
/*    according to above functions, assuming an error tolerance---to be used	*/
/*    when above routine fails to converge in [0,1]x[0,1].						*/
/********************************************************************************/
bool Zerofinder::findroot_boundary(double &u, double &v, DbVector3& p_bar, DbVector3& n_vec,
	  const DbVector3& x_bar, 
	  const DbVector3& a, const DbVector3& b, const DbVector3& c,
	  const DbVector3& n_A, const DbVector3& n_B, const DbVector3& n_C)
{
	A = a;
	B = b;
	C = c;
	nA = n_A;
	nB = n_B;
	nC = n_C;

	xbar = x_bar;


	// Search for minimum along side AB, AC, BC
    double alpha, beta, gamma;
	double funmin0, funmin1, funmin2;
	double umin0, umin1, umin2, vmin0, vmin1, vmin2;
	double funmin, umin, vmin;

	// Side AB is T_V0, u = [0, 1] side
	alpha = T_V0(0.0);
	beta = T_V0(0.5);
	gamma = T_V0(1.0);

	if (beta <= alpha && beta <= gamma) {
		double (Zerofinder::*f)(double) = &ThallCode::Zerofinder::T_V0;
		double (Zerofinder::*df)(double) = &ThallCode::Zerofinder::dT_V0;
		funmin0 = dbrent(0.0, 0.5, 1.0, LINMINTOL, f, df, &umin0);
		vmin0 = 0.0;
	}
	else {
		funmin0 = DBL_MAX;
	}

	// Side AC is T_U0, v = [0, 1] side
	alpha = T_U0(0.0);
	beta = T_U0(0.5);
	gamma = T_U0(1.0);

	if (beta <= alpha && beta <= gamma) {
		double (Zerofinder::*f)(double) = &ThallCode::Zerofinder::T_U0;
		double (Zerofinder::*df)(double) = &ThallCode::Zerofinder::dT_U0;
		funmin1 = dbrent(0.0, 0.5, 1.0, LINMINTOL, f, df, &vmin1);
		umin1 = 0.0;
	}
	else {
		funmin1 = DBL_MAX;
	}

	// Side CB is T_V1, u = [0, 1] side
	alpha = T_V1(0.0);
	beta = T_V1(0.5);
	gamma = T_V1(1.0);

	if (beta <= alpha && beta <= gamma) {
		double (Zerofinder::*f)(double) = &ThallCode::Zerofinder::T_V1;
		double (Zerofinder::*df)(double) = &ThallCode::Zerofinder::dT_V1;
		funmin2 = dbrent(0.0, 0.5, 1.0, LINMINTOL, f, df, &umin2);
		vmin2 = 1.0;
	}
	else {
		funmin2 = DBL_MAX;
	}

	if (funmin0 < DBL_MAX || funmin1 < DBL_MAX || funmin2 < DBL_MAX) {

		if (funmin0 < funmin1) {
			if (funmin0 < funmin2) {
				funmin = funmin0;
				umin = umin0;
				vmin = vmin0;
				p_bar = pbar_funV0(umin, vmin);
				n_vec = nvec_funV0(umin, vmin);
			}
			else {
				funmin = funmin2;
				umin = umin2;
				vmin = vmin2;
				p_bar = pbar_funV1(umin, vmin);
				n_vec = nvec_funV1(umin, vmin);
			}
		}
		else if (funmin1 < funmin2) {
			funmin = funmin1;
			umin = umin1;
			vmin = vmin1;
			p_bar = pbar_funU0(umin, vmin);
			n_vec = nvec_funU0(umin, vmin);
		}
		else {
			funmin = funmin2;
			umin = umin2;
			vmin = vmin2;
			p_bar = pbar_funV1(umin, vmin);
			n_vec = nvec_funV1(umin, vmin);
		}

		u = umin;
		v = vmin;

/*		if (funmin < 0.001) {

//			fprintf(stderr, "in rootfind_boundary, got funmin = %.15f at (u, v) = (%f, %f)\n",
//			        funmin, umin, vmin);

			glBegin(GL_POINTS);
			glColor3d(1.0 - funmin, 1.0 - funmin, 1.0 - funmin);
			glVertex3d(p_bar.X(), p_bar.Y(), p_bar.Z());
			glEnd();

			DbVector3 n_bar = n_vec.normalize();
			n_bar *= 0.99*(p_bar - x_bar).length();

			glColor3d(0.0, 1.0, 0.0);
			glBegin(GL_LINES);
			glVertex3d(p_bar.X(), p_bar.Y(), p_bar.Z());
			glVertex3d(p_bar.X() + n_bar.X(), p_bar.Y() + n_bar.Y(), p_bar.Z() + n_bar.Z());
			glEnd();

			glColor3d(0.0, 1.0, 1.0);
			glBegin(GL_LINES);
			glVertex3d(p_bar.X() + n_bar.X(), p_bar.Y() + n_bar.Y(), p_bar.Z() + n_bar.Z());
			glVertex3d(xbar.X(), xbar.Y(), xbar.Z());
			glEnd();
		}
*/
		// NOW, been a problem with false roots found, simply functional minima.
		//   Return false if n_vec.dot(xbar - p_bar) != 0;
		//   ---alternate to return false if minval > some_eps, but not sure what limit to use
		if (fabs(n_vec.normalize().dot((xbar - p_bar).normalize())) < 1.0 - TESTTOL) {
//		if (fabs(minval) > TESTTOL) {
//			fprintf(stderr, "bad root found, just a minval = %f\n", minval);
			return false;
		}
		return true;
	}
	else
		return false;
}

#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);

#define MOV3(a,b,c, d,e,f) (a)=(d);(b)=(e);(c)=(f);

const int DBRENT_ITMAX = 100;
const double ZEPS = 1.0e-10;

double Zerofinder::dbrent(double ax, double bx, double cx, double tol, 
						  double (Zerofinder::*f)(double), double (Zerofinder::*df)(double), double *xmin)
{
	int iter,ok1,ok2;
	double a,b,d1,d2,du,dv,dw,dx,e=0.0, d=0.0;
	double fu,fv,fw,fx,olde,tol1,tol2,u,u1,u2,v,w,x,xm;

	a=(ax < cx ? ax : cx);
	b=(ax > cx ? ax : cx);
	x=w=v=bx;


//	double (Zerofinder::*f)(double) = Zerofinder::T_V0;
//	double (Zerofinder::*df)(double) = Zerofinder::dT_V0;

	fw = fv = fx = (this->*f)(x);
	dw = dv = dx = (this->*df)(x);

	for (iter = 1; iter <= DBRENT_ITMAX; iter++) {
		xm = 0.5*(a + b);
		tol1 = tol*fabs(x) + ZEPS;
		tol2 = 2.0*tol1;
		if (fabs(x - xm) <= (tol2 - 0.5*(b - a))) {
			*xmin = x;
			return fx;
		}
		if (fabs(e) > tol1) {
			d1=2.0*(b-a);
			d2=d1;
			if (dw != dx) d1=(w-x)*dx/(dx-dw);
			if (dv != dx) d2=(v-x)*dx/(dx-dv);
			u1=x+d1;
			u2=x+d2;
			ok1 = (a-u1)*(u1-b) > 0.0 && dx*d1 <= 0.0;
			ok2 = (a-u2)*(u2-b) > 0.0 && dx*d2 <= 0.0;
			olde=e;
			e=d;
			if (ok1 || ok2) {
				if (ok1 && ok2)
					d=(fabs(d1) < fabs(d2) ? d1 : d2);
				else if (ok1)
					d=d1;
				else
					d=d2;
				if (fabs(d) <= fabs(0.5*olde)) {
					u=x+d;
					if (u-a < tol2 || b-u < tol2)
						d=SIGN(tol1,xm-x);
				}
				else {
					d=0.5*(e=(dx >= 0.0 ? a-x : b-x));
				}
			}
			else {
				d=0.5*(e=(dx >= 0.0 ? a-x : b-x));
			}
		} else {
			d=0.5*(e=(dx >= 0.0 ? a-x : b-x));
		}
		if (fabs(d) >= tol1) {
			u=x+d;
			fu = (this->*f)(u);
		} else {
			u=x+SIGN(tol1,d);
			fu = (this->*f)(u);
			if (fu > fx) {
				*xmin=x;
				return fx;
			}
		}
		du = (this->*df)(u);
		if (fu <= fx) {
			if (u >= x) a=x; else b=x;
			MOV3(v,fv,dv, w,fw,dw)
			MOV3(w,fw,dw, x,fx,dx)
			MOV3(x,fx,dx, u,fu,du)
		} else {
			if (u < x) a=u; else b=u;
			if (fu <= fw || w == x) {
				MOV3(v,fv,dv, w,fw,dw)
				MOV3(w,fw,dw, u,fu,du)
			} else if (fu < fv || v == x || v == w) {
				MOV3(v,fv,dv, u,fu,du)
			}
		}
	}
//	fprintf(stderr, "Too many iterations in routine dbrent");
	return 0.0;
}


