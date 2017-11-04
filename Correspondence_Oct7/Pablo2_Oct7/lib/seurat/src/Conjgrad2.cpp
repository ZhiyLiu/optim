/********************************************************************************/
/*																				*/
/*  	File	:  Conjgrad.cpp													*/
/*																				*/
/*	Description:  class encapsulated versions of conjugate gradient routines	*/
/*		derived from NumRecinC and previous ThallCode implementation.			*/
/*																				*/
/*	Project :  Pablo															*/
/*																				*/
/*	Author  :  A. Thall															*/
/*																				*/
/*	Date	:  15. September 2001												*/
/*	Modifications:																*/
/********************************************************************************/
#include <math.h>
#include <stdio.h>
#define D_CONJGRAD
#include "Shapedepend.h"

using namespace ThallCode;
using namespace CGconstants;

// Some annoying macros (maxarg1 & 2 are in namespace CGconstants)
//   tried to replace with inline functions, but got screwy results.
//   Should fix, but fuck it hell who has the time?
// (See the problem:  there's a SHFT I think that's call by reference
//  but got inlined as a call by value, so a new value didn't get evaluated
//  for a later term in the shift.  Reminds me of that Programming Languages
//  course I took 8 or 10 years ago.)  (What fucked up code, gotta love those
//  NumRec guys, Scorpius needs to pay them a visit.)

// used in macro calls (want to avoid collisions)
static double CGmaxarg1, CGmaxarg2;
#define FMAX(a,b) (CGmaxarg1=(a),CGmaxarg2=(b),(CGmaxarg1) > (CGmaxarg2) ?\
        (CGmaxarg1) : (CGmaxarg2))

#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);

#define MOV3(a,b,c, d,e,f) (a)=(d);(b)=(e);(c)=(f);


double Conjgrad::dbrent(double ax, double bx, double cx, double tol, double *xmin)
{
	int iter,ok1,ok2;
	double a,b,d,d1,d2,du,dv,dw,dx,e=0.0;
	double fu,fv,fw,fx,olde,tol1,tol2,u,u1,u2,v,w,x,xm;

	d = 0.0;
	a=(ax < cx ? ax : cx);
	b=(ax > cx ? ax : cx);
	x=w=v=bx;
	fw = fv = fx = f1dim(x);
	dw = dv = dx = df1dim(x);
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
				} else {
					d=0.5*(e=(dx >= 0.0 ? a-x : b-x));
				}
			} else {
				d=0.5*(e=(dx >= 0.0 ? a-x : b-x));
			}
		} else {
			d=0.5*(e=(dx >= 0.0 ? a-x : b-x));
		}
		if (fabs(d) >= tol1) {
			u=x+d;
			fu = f1dim(u);
		} else {
			u=x+SIGN(tol1,d);
			fu = f1dim(u);
			if (fu > fx) {
				*xmin=x;
				return fx;
			}
		}
		du = df1dim(u);
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
	fprintf(stderr, "Too many iterations in routine dbrent");
	return 0.0;
}

void Conjgrad::mnbrak(double *ax, double *bx, double *cx, double *fa, double *fb, double *fc)
{
	double ulim,u,r,q,fu;
	double dum = 0.0;

	*fa = f1dim(*ax);
	*fb = f1dim(*bx);
	if (*fb > *fa) {
		SHFT(dum,*ax,*bx,dum)
		SHFT(dum,*fb,*fa,dum)
	}
	*cx = (*bx)+GOLD*(*bx-*ax);
	*fc = f1dim(*cx);
	while (*fb > *fc) {
		r = (*bx-*ax)*(*fb-*fc);
		q = (*bx-*cx)*(*fb-*fa);
		u = (*bx)-((*bx-*cx)*q-(*bx-*ax)*r)/
			(2.0*SIGN(FMAX(fabs(q-r),TINY), q-r));
		ulim = (*bx) + GLIMIT*(*cx-*bx);
		if ((*bx-u)*(u-*cx) > 0.0) {
			fu = f1dim(u);
			if (fu < *fc) {
				*ax=(*bx);
				*bx=u;
				*fa=(*fb);
				*fb=fu;
				return;
			} else if (fu > *fb) {
				*cx=u;
				*fc=fu;
				return;
			}
			u=(*cx)+GOLD*(*cx-*bx);
			fu = f1dim(u);
		} else if ((*cx-u)*(u-ulim) > 0.0) {
			fu = f1dim(u);
			if (fu < *fc) {
				SHFT(*bx, *cx, u, *cx+GOLD*(*cx-*bx))
				SHFT(*fb, *fc, fu, f1dim(u))
			}
		} else if ((u-ulim)*(ulim-*cx) >= 0.0) {
			u = ulim;
			fu = f1dim(u);
		} else {
			u = (*cx) + GOLD*(*cx-*bx);
			fu = f1dim(u);
		}
		SHFT(*ax,*bx,*cx,u)
		SHFT(*fa,*fb,*fc,fu)
	}
}

double Conjgrad::f1dim(double x)
{
	DbVector2 xt;

	xt = pcom + xicom*x;

	return (*myfun)(xt);
}

double Conjgrad::df1dim(double x)
{
	DbVector2 xt, df;

	xt = pcom + xicom*x;

	(*mydfun)(xt,df);

	return df.dot(xicom);
}

void Conjgrad::dlinmin(DbVector2& p, DbVector2& xi, double *fret)
{
	double xx,xmin,fx,fb,fa,bx,ax;

	pcom = p;
	xicom = xi;

	ax=0.0;
	xx=1.0;
	mnbrak(&ax,&xx,&bx,&fa,&fx,&fb);
	*fret=dbrent(ax,xx,bx, DLINMIN_TOL,&xmin);

	xi *= xmin;
	p += xi;
}

void Conjgrad::frprmn(DbVector2& p, double ftol, int *iter, double *fret)
{
	int its;
	double gg,gam,fp,dgg;
	DbVector2 g, h, xi;

	fp = (*myfun)(p);
	(*mydfun)(p,xi);

	g = -xi;
	xi = h = g;

	for (its=1; its <= FRPRMN_ITMAX; its++) {
		*iter=its;
		dlinmin(p, xi, fret);
		if (2.0*fabs(*fret-fp) <= ftol*(fabs(*fret)+fabs(fp)+EPS)) {
			return;
		}
		fp = (*myfun)(p);
		(*mydfun)(p,xi);
		dgg=gg=0.0;

		gg = g.dot(g);
		dgg = xi.dot(xi + g);

		if (gg == 0.0) {
			return;
		}

		gam=dgg/gg;

		g = -xi;
		xi = h = g + h*gam;
	}
//	fprintf(stderr, "Too many iterations in frprmn");
}


