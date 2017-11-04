/******************************************************************
 * OPTIMA Library                                                 *
 ******************************************************************
 * Author:						Paul Yushkevich
 *
 * Date:							Apr 1, 1999
 *
 * Description					Multidimensional optimization algorithms
 *									See http://www.cs.unc.edu/~pauly/optima
 *									
 *	Sources:						"Numerical Recepies in C",
 *									Michaelewitz, "Genetic Algorithms + Data
 *									Structures = Evolutionary Programs"
 *
 * Dependencies:				PY Matrix library, CLAPACK
 ******************************************************************
 * BrentLinearMethod.cpp
 *	---------------------
 * This method from NRC-10.2 allows us to optimize along a vector in
 * n-space.
 ******************************************************************/


#include <math.h>
#include "BrentLinearMethod.h"
#include <stdio.h>

// Begin namespace
NAMESPACE_PAULY_START

#define GOLD 1.618034
#define GLIMIT 100.0
//#define GLIMIT .05
#define TINY 1.0e-20
#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);

#define ITMAX 100
//#define ITMAX 10
#define CGOLD 0.3819660
#define ZEPS 1.0e-10
#define TOL 2.0e-4
#define MOV3(a,b,c, d,e,f) (a)=(d);(b)=(e);(c)=(f);

#define MAX_TRIES_FOR_BRACKETING	10

//#define dprintf printf

/*
#include <stdarg.h>
template<class T>
inline static void dprintf(const char* format, const T& arg1, ...) {
	if(globalVerbosity >= 2) {
		va_list args;
		va_start(args, arg1);
		vprintf( format, args );
		va_end(args);
	}
}
*/
#ifndef dprintf
#define dprintf donothing
inline static void donothing(...) {}
#endif

inline bool ISNAN(const double a) { return a != a; }

inline double SIGN(double a,double b) { return (b<0) ? -a : a;}
inline double FMAX(double a,double b) { return (a<b) ? b : a;}


double BrentLinearMethod::defaultBrentLinearSearchBoundFactor	= GOLD;

/*******************************************************************
 Solve a problem of finding s that optimizes f(P+sn)
 ******************************************************************/
BrentLinearMethod::BrentLinearMethod(Function &problem,Vector inX,Vector inN)
	:  p(problem)
{
	x = inX;
	n = inN;
	value = 0;
    xInitialLimit	= GLIMIT / (GOLD * 10.0);
	boundsFunction	= NULL;
	brentLinearSearchBoundFactor	= defaultBrentLinearSearchBoundFactor;
}

double BrentLinearMethod::func(double t)
{
	Vector xt = x + t*n;
	return p.evaluate(xt);
}

bool BrentLinearMethod::run()
{
//	double xx = /*1.0*/GLIMIT/(GOLD*10.),
    double xmin, fc, fb, fa, bx, cx, ax = 0.0;
    bx = xInitialLimit / brentLinearSearchBoundFactor;
	dprintf("OPT: Calling mnbrak\n");
	if(!mnbrak(&ax, &bx, &cx, &fa, &fb, &fc) ) {
		// Could not find a bracket.
		n *= cx;
		x += n;
		value = fc;
		// Signal that minimum along direction could not
		// be found ( because we reached a bound. )
		return false;
	}
	dprintf("OPT: Calling brent\n");
	value = brent(ax, bx, cx, fb, TOL, &xmin);
	n *= xmin;
	x += n;
	return true;
}

// Given a function func, and given distinct initial points ax and bx, this routine searches in
// the downhill direction (defned by the function as evaluated at the initial points) and returns
// new points ax, bx, cx that bracket a minimum of the function. Also returned are the function
// values at the three points, fa, fb, and fc.
// Here GOLD is the default ratio by which successive intervals are magnifed;
// GLIMIT is the maximum magnifcation allowed for a parabolic-ft step.
bool BrentLinearMethod::mnbrak(double *ax, double *bx, double *cx,
									  double *fa, double *fb, double *fc)
{
	double ulim,u,r,q,fu,dum;
	u	 = -1.0;	// any value that doesn't make sense. The value is never used in actual
					// computation.
	const double orig_ax	= *ax;
	const double orig_bx	= *bx;
	int iters	= 0;
	*fa=func(*ax);
	*fb=func(*bx);
	if (*fb > *fa) {
		//Switch roles of a and b so that we can go
		// downhill in the direction from a to b.
		SHFT(dum,*ax,*bx,dum);
		SHFT(dum,*fb,*fa,dum);
	}

	// First guess for c.
	*cx=(*bx)+GOLD*(*bx-*ax);
	*fc=func(*cx);

	dprintf("MNBRAK: %6f %6f %6f %6f %6f %6f %6f\n", *ax, *bx, *cx, u, *fa, *fb, *fc );

	while (*fb > *fc) {
		if( iters++ > MAX_TRIES_FOR_BRACKETING ) {
			printf("MNBRAK: <<< bracketing failed >>>\n");
			return false;
		}
		// Keep returning here until we bracket.
		// Compute u by parabolic extrapolation from a; b; c. TINY is used to
		// prevent any possible division by zero.
		r=(*bx-*ax)*(*fb-*fc);
		q=(*bx-*cx)*(*fb-*fa);
		u=(*bx)-((*bx-*cx)*q-(*bx-*ax)*r)/(2.0*SIGN(FMAX(fabs(q-r),TINY),q-r));
		if( boundsFunction ) {
			ulim	= boundsFunction->limit(x, n, u);
		}
		else {
			ulim=(*bx)+GLIMIT*(*cx-*bx);
		}
		// We won't go farther than this limit. Test different possibilities:

		dprintf("MNBRAK: %6f %6f %6f %6f %6f %6f %6f\n", *ax, *bx, *cx, u, *fa, *fb, *fc );

		if ((*bx-u)*(u-*cx) > 0.0) {
			// Parabolic u is between b and c: try it.
			fu=func(u);
			dprintf("MNBRAK: fu %6f\n", fu);
			if (fu < *fc) {
				// Got a minimum between b and c.
				*ax=(*bx);
				*bx=u;
				*fa=(*fb);
				*fb=fu;
				dprintf("MNBRAK: %6f %6f %6f %6f %6f %6f %6f\n", *ax, *bx, *cx, u, *fa, *fb, *fc );
				return true;
			}
			else if (fu > *fb) {
				// Got a minimum between between a and u.
				*cx=u;
				*fc=fu;
				dprintf("MNBRAK: %6f %6f %6f %6f %6f %6f %6f\n", *ax, *bx, *cx, u, *fa, *fb, *fc );
				return true;
			}

			u=(*cx)+GOLD*(*cx-*bx);
			// Parabolic fit was no use. Use default mag-nification.
			fu=func(u);
		} else if ((*cx-u)*(u-ulim) > 0.0) {
			//Parabolic fit is between c and its allowed limit.
			fu=func(u);
			if (fu < *fc) {
				SHFT(*bx,*cx,u,*cx+GOLD*(*cx-*bx));
				SHFT(*fb,*fc,fu,func(u));
			}
		} else if ((u-ulim)*(ulim-*cx) >= 0.0) {
			//Limit parabolic u to maximum allowed value.
			printf("MNBRAK: <<< bracketing limited by supplied bound >>>\n");
			u=ulim;
			fu=func(u);
		} else {
			// Reject parabolic u, use default magnifica-tion.
			u=(*cx)+GOLD*(*cx-*bx);
			fu=func(u);
		}
		SHFT(*ax,*bx,*cx,u); // Eliminate oldest point and continue.
		SHFT(*fa,*fb,*fc,fu);
	}
	return true;
}


// Given a function f, and given a bracketing triplet of abscissas ax, bx, cx (such that bx is
// between ax and cx, and f(bx) is less than both f(ax) and f(cx)), this routine isolates
// the minimum to a fractional precision of about tol using Brent's method. The abscissa of
// the minimum is returned as xmin, and the minimum function value is returned as brent, the
// returned function value.
double BrentLinearMethod::brent(double ax, double bx, double cx, double fbx,
										double tol, double *xmin)
{
	int iter;
	double a,b,d,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm;
	u	 = -1.0;	// any value that doesn't make sense. The value is never used in actual
					// computation.
	double e=0.0;								
	// This will be the distance moved on the step before last.
	// a and b must be in ascending order,
	// but input abscissas need not be.

	// Initializations...
	d = 0.0;
	a=(ax < cx ? ax : cx);					
	b=(ax > cx ? ax : cx);
	x=w=v=bx;									
	fw=fv=fx=fbx;

	// Main program loop.
	for (iter=1;iter<=ITMAX;iter++) {	
		dprintf("BRENT: %6f %6f %6f %6f %6f %6f\n", a, b, x, w, v, u);
		xm=0.5*(a+b);
		tol2=2.0*(tol1=tol*fabs(x)+ZEPS);

		// Test for done here.
		if (fabs(x-xm) <= (tol2-0.5*(b-a))) {	
			*xmin=x;
			return fx;
		}
		if (fabs(e) > tol1) {				
			// Construct a trial parabolic fit
			r=(x-w)*(fx-fv);
			q=(x-v)*(fx-fw);
			p=(x-v)*q-(x-w)*r;
			q=2.0*(q-r);
			if (q > 0.0) p = -p;
			q=fabs(q);
			etemp=e;
			e=d;
			if (fabs(p) >= fabs(0.5*q*etemp) || p <= q*(a-x) || p >= q*(b-x))
				d=CGOLD*(e=(x >= xm ? a-x : b-x));

			// We check for bracketing by (ax,cx) here among others.
			// The above conditions determine the acceptability of the parabolic fit. Here we
			// take the golden section step into the larger of the two segments.
			else {
				// Take the parabolic step.
				d=p/q;							
				u=x+d;
				if (u-a < tol2 || b-u < tol2)
					d=SIGN(tol1,xm-x);
			}
		}
		else {
			d=CGOLD*(e=(x >= xm ? a-x : b-x));
		}
		u=(fabs(d) >= tol1 ? x+d : x+SIGN(tol1,d));
		fu=func(u);
		// This is the one function evaluation per iteration.
		if (fu <= fx) {						
			// Now decide what to do with our function evaluation.
			if (u >= x) a=x; else b=x;
			SHFT(v,w,x,u);						
			SHFT(fv,fw,fx,fu);
		}
		else {
			//Housekeeping follows:
			if (u < x) a=u; else b=u;
			if (fu <= fw || w == x) {
				v=w;
				w=u;
				fv=fw;
				fw=fu;
			} else if (fu <= fv || v == x || v == w) {
				v=u;
				fv=fu;
			}
		} //Done with housekeeping. Back for another iteration.
	}

	dassert(0);
	*xmin=x; //Never get here.
	return fx;
}

/*******************************************************************
 Solve a problem of finiding s that optimizes f(P+sn) with f'
 ******************************************************************/
BrentDerivativeMethod::BrentDerivativeMethod(DifferentiableFunction &problem,Vector inX,Vector inN) :
p(problem)
{
	x = inX;
	n = inN;
	value = 0;
}

double BrentDerivativeMethod::func(double t) {
	Vector xt = x + t*n;
	return p.evaluate(xt);
}

double BrentDerivativeMethod::dfunc(double t,double &fn) {
	Vector xt = x + t*n;
	return p.computeDirectionalJet(xt,n,fn);
}

void BrentDerivativeMethod::run() {
	double xx = 1.0,xmin,fx,fb,fa,bx,ax = 0.0;
	mnbrak(&ax,&xx,&bx,&fa,&fx,&fb);
	value = dbrent(ax,xx,bx,TOL,&xmin);
	n *= xmin;
	x += n;
}

// Given a function func, and given distinct initial points ax and bx, this routine searches in
// the downhill direction (defned by the function as evaluated at the initial points) and returns
// new points ax, bx, cx that bracket a minimum of the function. Also returned are the function
// values at the three points, fa, fb, and fc.
// Here GOLD is the default ratio by which successive intervals are magnifed;
// GLIMIT is the maximum magnifcation allowed for a parabolic-ft step.
void BrentDerivativeMethod::mnbrak(double *ax, double *bx, double *cx,
									  double *fa, double *fb, double *fc)
{
	double ulim,u,r,q,fu,dum;
	*fa=func(*ax);
	*fb=func(*bx);
	if (*fb > *fa) {
		//Switch roles of a and b so that we can go
		// downhill in the direction from a to b.
		SHFT(dum,*ax,*bx,dum);
		SHFT(dum,*fb,*fa,dum);
	}

	// First guess for c.
	*cx=(*bx)+GOLD*(*bx-*ax);
	*fc=func(*cx);

	while (*fb > *fc) {
		// Keep returning here until we bracket.
		// Compute u by parabolic extrapolation from a; b; c. TINY is used to
		// prevent any possible division by zero.
		r=(*bx-*ax)*(*fb-*fc);
		q=(*bx-*cx)*(*fb-*fa);
		// This expression is an estimate for x that minimizes -
		// y = a * (x-b)^2 + c
		u=(*bx)-((*bx-*cx)*q-(*bx-*ax)*r)/(2.0*SIGN(FMAX(fabs(q-r),TINY),q-r));
		ulim	= (*bx)+GLIMIT*(*cx-*bx);
		// We won't go farther than the above limits. Test various possibilities:

		if ((*bx-u)*(u-*cx) > 0.0) {
			// Parabolic u is between b and c: try it.
			fu=func(u);
			if (fu < *fc) {
				// Got a minimum between b and c.
				*ax=(*bx);
				*bx=u;
				*fa=(*fb);
				*fb=fu;
				return;
			}
			else if (fu > *fb) {
				// Got a minimum between between a and u.
				*cx=u;
				*fc=fu;
				return;
			}

			u=(*cx)+GOLD*(*cx-*bx);
			// Parabolic fit was no use. Use default magnification.
			fu=func(u);
		} else if ((*cx-u)*(u-ulim) > 0.0) {
			//Parabolic fit is between c and its allowed limit.
			fu=func(u);
			if (fu < *fc) {
				SHFT(*bx,*cx,u,*cx+GOLD*(*cx-*bx));
				SHFT(*fb,*fc,fu,func(u));
			}
		} else if ((u-ulim)*(ulim-*cx) >= 0.0) {
			// Limit parabolic u to maximum allowed value
			u=ulim;
			fu=func(u);
		} else {
			// Reject parabolic u, use default magnification.
			u=(*cx)+GOLD*(*cx-*bx);
			fu=func(u);
		}
		SHFT(*ax,*bx,*cx,u); // Eliminate oldest point and continue.
		SHFT(*fa,*fb,*fc,fu);
	}
}


// Given a function f, and given a bracketing triplet of abscissas ax, bx, cx (such that bx is
// between ax and cx, and f(bx) is less than both f(ax) and f(cx)), this routine isolates
// the minimum to a fractional precision of about tol using Brent's method. The abscissa of
// the minimum is returned as xmin, and the minimum function value is returned as brent, the
// returned function value.
double BrentDerivativeMethod::dbrent(double ax, double bx, double cx,
												 double tol, double *xmin)
{
	// Will be used as ags for whether proposed steps are acceptable or not.
	int iter,ok1,ok2;
	double a,b,d,d1,d2,du,dv,dw,dx,e=0.0;
	double fu,fv,fw,fx,olde,tol1,tol2,u,u1,u2,v,w,x,xm;

	d = 0.0;
	a=(ax < cx ? ax : cx);
	b=(ax > cx ? ax : cx);
	x=w=v=bx;
	fw=fv=fx=dfunc(x,dx);
	dw=dv=dx;
	// All our housekeeping chores are doubled by the necessity of moving
	// derivative values around as well as function values.
	for (iter=1;iter<=ITMAX;iter++) {
		xm=0.5*(a+b);
		tol1=tol*fabs(x)+ZEPS;
		tol2=2.0*tol1;
		if (fabs(x-xm) <= (tol2-0.5*(b-a))) {
			*xmin=x;
			return fx;
		}
		if (fabs(e) > tol1) {
			// Initialize these d's to an out-of-bracket value.
			d1=2.0*(b-a);
			d2=d1;
			// Secant method with one point.  Which of these two estimates of d shall we take?
			// We will insist that they be within the bracket, and on the side pointed to
			// by the derivative at x:
			if (dw != dx) d1=(w-x)*dx/(dx-dw);

			if (dv != dx) d2=(v-x)*dx/(dx-dv);

			u1=x+d1;
			u2=x+d2;
			ok1 = (a-u1)*(u1-b) > 0.0 && dx*d1 <= 0.0;
			ok2 = (a-u2)*(u2-b) > 0.0 && dx*d2 <= 0.0;
			// Movement on the step before last.
			olde=e;
			e=d;

			if (ok1 || ok2) {
				// Take only an acceptable d, and if both are acceptable, then take
				// the smallest one.
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
					// Bisect, not golden section.
					d=0.5*(e=(dx >= 0.0 ? a-x : b-x));
					// Decide which segment by the sign of the derivative.
				}
			}
			else {
				d=0.5*(e=(dx >= 0.0 ? a-x : b-x));
			}
		}
		else {
			d=0.5*(e=(dx >= 0.0 ? a-x : b-x));
		}

		// Optimization (paul)
		double _dfu,_fu;
		_fu = dfunc(u,_dfu);

		if (fabs(d) >= tol1) {
			u=x+d;
			fu=_fu;
		}
		else {
			u=x+SIGN(tol1,d);
			fu=_fu;
			if (fu > fx) {
				// If the minimum step in the downhill direction takes us uphill, then
				// we are done.
				*xmin=x;
				return fx;
			}
		}

		// Now all the housekeeping, sigh.
		du=_dfu;
		if (fu <= fx) {
			if (u >= x) a=x; else b=x;
			MOV3(v,fv,dv, w,fw,dw)
				MOV3(w,fw,dw, x,fx,dx)
				MOV3(x,fx,dx, u,fu,du)
		}
		else {
			if (u < x) a=u; else b=u;
			if (fu <= fw || w == x) {
				MOV3(v,fv,dv, w,fw,dw)
					MOV3(w,fw,dw, u,fu,du)
			}
			else if (fu < fv || v == x || v == w) {
				MOV3(v,fv,dv, u,fu,du)
			}
		}
	}

	// printf("Too many iterations in routine dbrent");
	dassert(0);
	return 0.0; // Never get here.
}

// End namespace
NAMESPACE_PAULY_END

