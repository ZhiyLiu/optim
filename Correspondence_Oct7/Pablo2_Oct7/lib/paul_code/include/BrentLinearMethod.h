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
 * BrentLinearMethod.h
 *	-------------------
 * This method from NRC-10.2 allows us to optimize along a vector in
 * n-space.
 ******************************************************************/
#ifndef _OPTIMA_LINMIN_
#define _OPTIMA_LINMIN_

// Include files
#include "optima.h"

// Begin namespace
NAMESPACE_PAULY_START


/*******************************************************************
 Solve a problem of finiding s that optimizes f(P+sn)

 Use this by creating a BrentLinearMethod and calling run.  Check
 the values of x,n and 'value' after the algorithm has finished.
 ******************************************************************/
 class BrentLinearMethod {
protected:
	Function &p;
    BoundsFunction * boundsFunction;

    double xInitialLimit;

	double brentLinearSearchBoundFactor;
	static double defaultBrentLinearSearchBoundFactor;

	double func(double);
	bool mnbrak(double *ax, double *bx, double *cx, double *fa, double *fb, double *fc);
	double brent(double ax, double bx, double cx, double fbx, double tol,double *xmin);

public:
	// Position of the solution
	Vector x;

	// Vector from old solution to new solution
	Vector n;

	// Optimal valur
	double value;

	// Constructor (does not require a solution space because there is a single
	// defined starting point.
	BrentLinearMethod(Function &p, Vector x, Vector n);

	// Runs the method until optimum is found.
	bool run();

    void setXInitialLimit(double _xInitialLimit)
    {
        if(_xInitialLimit < xInitialLimit)
        {
            //printf("Set xInitialLimit: %f -> %f\n", xInitialLimit, _xInitialLimit);
            xInitialLimit = _xInitialLimit;
        }
    }
    void setBoundsFunction(BoundsFunction * _boundsFunction)
    {
        boundsFunction = _boundsFunction;
    }

	void setBrentLinearSearchBoundFactor( const double _factor )
	{
		brentLinearSearchBoundFactor	= _factor;
	}

	static void setDefaultBrentLinearSearchBoundFactor( const double _factor )
	{
		defaultBrentLinearSearchBoundFactor	= _factor;
	}
};

/*******************************************************************
 Solve a problem of finiding s that optimizes f(P+sn) This uses derivatives
 to speed up the process.

 Use this by creating a BrentDerivativeMethod and calling run.  Check
 the values of x,n and 'value' after the algorithm has finished.
 ******************************************************************/
class BrentDerivativeMethod {
protected:
	DifferentiableFunction &p;

	double func(double t);
	double dfunc(double t,double &fn);
	void mnbrak(double *ax, double *bx, double *cx, double *fa, double *fb, double *fc);
	double dbrent(double ax, double bx, double cx, double tol,double *xmin);

public:
	// Position of the solution
	Vector x;

	// Vector from old solution to new solution
	Vector n;

	// Optimal valur
	double value;

	// Constructor (does not require a solution space because there is a single
	// defined starting point.
	BrentDerivativeMethod(DifferentiableFunction &p,Vector x,Vector n);

	// Runs the method until optimum is found.
	void run();
};



// End namespace
NAMESPACE_PAULY_END

#endif /* _OPTIMA_LINMIN_ */

