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
 * ConjugateGradientMethod.h
 *	-------------------------
 * This method from NRC-10.6 allows us to optimize nicely using first
 * partial derivatives in n-space.
 ******************************************************************/
#ifndef _OPTIMA_CONJGRAD_
#define _OPTIMA_CONJGRAD_

// Include files
#include "optima.h"

// Begin namespace
NAMESPACE_PAULY_START



/*******************************************************************
 Optimize f(x) using conjugate gradient descent 

 Usage:
	-	Create a ConjugateGradientMethod (with a function and a sample x)
	-	Set options with setXXX methods
	-   Call performIteration() until isFinished==true
	-	Check results with getBestEverXXX()
 ******************************************************************/
class ConjugateGradientMethod : public NumericalMethod {
protected:
	NumericalSolution bestEver,current;
	DifferentiableFunction &p;
    BoundsFunction * boundsFunction;

	// Input parameters
	double fTolerance; 
    Vector xTolerance;

	// Variables used between iterations:
	//   g is grad f at x
	//   h is the direction to minimize in 1D
	//   prevG is the gradient at the previous step
	Vector g, h, prevG; 

	// Count the number of iterations performed.
	// isFinished() should return false until we've done N iterations
	// for an N dimensional problem
	int itersCount;
	int nDims;

	// Initialize
	void initialize();

	// Are we finished
	bool finished;

	double brentLinearSearchBoundFactor;

public:
	// Initialize the method
	ConjugateGradientMethod(DifferentiableFunction &problem,const Vector &v);

	// Set the value for epsilon, the terminating condition
	void setValueTolerance(double fTol) {
		fTolerance = fTol;
	}

	void setParameterTolerance(Vector xTol) {
		xTolerance = xTol;
	}


   // Perform an iteration of the method
   void performIteration();

	// Check termination condition (termination determined by tolerance)
	bool isFinished() { return finished; };

	// Solution report methods
	double getBestEverValue() {return bestEver.getValue();}
	const Vector &getBestEverX() {return bestEver.x;}

	// Get the problem on which we are operating
	virtual Function &getFunction() {
		return p;
	}

    void setBoundsFunction(BoundsFunction * _boundsFunction)
    {
        boundsFunction = _boundsFunction;
    }

	void setBrentLinearSearchBoundFactor( const double _factor )
	{
		brentLinearSearchBoundFactor	= _factor;
	}
};

// End namespace
NAMESPACE_PAULY_END

#endif	/* _OPTIMA_CONJGRAD_ */

