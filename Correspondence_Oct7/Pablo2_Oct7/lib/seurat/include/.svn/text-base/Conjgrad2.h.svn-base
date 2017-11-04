/********************************************************************************/
/*																				*/
/*  	File	:  Conjgrad.h													*/
/*																				*/
/*	Description:  class header file for encapsulated versions of conjugate		*/
/*		gradient routines derived from NumRecinC and previous ThallCode			*/
/*		implementation.															*/
/*																				*/
/*	Project :  Pablo															*/
/*																				*/
/*	Author  :  A. Thall															*/
/*																				*/
/*	Date	:  15. September 2001												*/
/*	Modifications:																*/
/********************************************************************************/

namespace CGconstants {

	const int DBRENT_ITMAX = 100;
	const double CGOLD = 0.3819660;
	const double ZEPS = 1.0e-10;

	// For mbrak()
	const double GOLD = 1.618034;
	const double GLIMIT = 100.0;
	const double TINY = 1.0e-20;

	// For dlinmin()
	const double DLINMIN_TOL = 2.0e-4;

	// for frprmn()
	const int FRPRMN_ITMAX = 15;
	const double EPS = 1.0e-10;
}

class Conjgrad
{
	double (*myfun)(DbVector2&);
	void (*mydfun)(DbVector2&, DbVector2&);

	DbVector2 pcom, xicom;

	double dbrent(double ax, double bx, double cx, double tol, double *xmin);
	void mnbrak(double *ax, double *bx, double *cx, double *fa, double *fb, double *fc);
	void dlinmin(DbVector2& p, DbVector2& xi, double *fret);
	double f1dim(double x);
	double df1dim(double x);

/*	double SIGN(double a, double b) {
	    return (b >= 0.0 ? fabs(a) : -fabs(a));
	}
	void SHFT(double& a, double& b, double& c, double d) {
		a = b; b = c; c = d;
	}
	void MOV3(double& a, double& b, double& c, double d, double e, double f) {
		a = d; b = e; c = f;
	}
 */
public:
	//Conjgrad();
	Conjgrad(double (*func)(DbVector2&), void (*dfunc)(DbVector2&, DbVector2&)) {
		myfun = func;
		mydfun = dfunc;
	};
	void init(double (*func)(DbVector2&), void (*dfunc)(DbVector2&, DbVector2&)) {
		myfun = func;
		mydfun = dfunc;
	};

	int maxiterations() { return CGconstants::FRPRMN_ITMAX; }

    void frprmn(DbVector2& p, double ftol, int *iter, double *fret);
};


