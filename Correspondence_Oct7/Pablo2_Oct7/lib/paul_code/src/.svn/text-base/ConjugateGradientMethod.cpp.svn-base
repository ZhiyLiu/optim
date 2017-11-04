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
#include <math.h>
#include <stdlib.h>
#include <memory.h>
#include <iostream>

#include "support.h"
#include "ConjugateGradientMethod.h"
#include "BrentLinearMethod.h"
#include "LogManager.h"

extern int globalVerbosity;

// Begin namespace
NAMESPACE_PAULY_START

/*****************************************************************
 ConjugateGradientMethod
 ****************************************************************/
ConjugateGradientMethod::ConjugateGradientMethod(DifferentiableFunction &problem, const Vector &v)
	: p(problem), current(v), fTolerance(1.0e-5), bestEver(v), brentLinearSearchBoundFactor(0.0)
{
	type = ConjugateGradient;

	// Since no optional parameters affect initialization, we perform
	// initialization right here and now.
	finished = false;
    boundsFunction = NULL;
	initialize();
    xTolerance.setSize(v.size());
    xTolerance.setAll(0.0);

	// We will count the number of iterations performed, and use that 
	// to say we can't be finshed until we've taken N steps,
	// where N is the number of dimensions of the problem
	itersCount = 0;
	nDims = v.size();
}

void ConjugateGradientMethod::initialize()
{
	// Compute gradient and function

#ifdef OPTIMIZATION_VISUALIZER
	globalLogManager.setCategory(0); // calculate gradient;
#endif
	double val = p.computeOneJet(current.x,g);
	current.setValue(val);

	// first direction to search is gradient 
    h = g;

	// Set best ever
	bestEver = current;
}

void ConjugateGradientMethod::performIteration() {

	itersCount++;

	// Perform linear minimization problem (non-derivative Brent seems to
	// outperform derivative brent by an order of 2.
/*
	Before the bug fixes in March 2005, h was a unit-vector.
	h needs to keep its length for conjugate gradient search
	in the brent linear search, when a bounds function exists
	the scaling of this vector has no effect.
	However when the bounds function is undefined (M3DRegistrationOptimizer)
	the scaling of this vector can send the brent search into 
	invalid regions of the search space, which leads to mathematical
	errors.

	The solution is to run the brent search along the unit-vector
	but to keep the scaled vector around for calculating the next
	search direction
*/
	Vector searchDir(h);
	searchDir.normalize();

	BrentLinearMethod opt(p,current.x, searchDir);
	if(brentLinearSearchBoundFactor != 0.0 ) {
		opt.setBrentLinearSearchBoundFactor(brentLinearSearchBoundFactor);
	}
	//
	// Evaluating the searchDir gives us the
	// the maximum number of steps of size |searchDir|
	// we can take before hitting any bound.
	// This bounds evaluation assumes that x is zero!!!
	// and also that the bounds are symmetric about
	// the origin.
	//
    if(boundsFunction) {
        opt.setXInitialLimit(boundsFunction->evaluate(searchDir));
		opt.setBoundsFunction(boundsFunction);
	}
#ifdef OPTIMIZATION_VISUALIZER
	globalLogManager.setCategory(1); // search for 1D minimum
#endif
	const bool brentSucceeded	= opt.run();

	// the next solution is the result of the 1d optimizer
	NumericalSolution next(opt.x);
	next.setValue(opt.value);

	// Update the best ever solution
	if(bestEver.getValue() > next.getValue())
		bestEver = next;

	// Stop Condition: If improvement in objective function is less than "nothing"
	double fThreshold = fTolerance*(fabs(next.getValue()) + fabs(current.getValue()) + 1.0e-10);
	double fChange = 2.0 * fabs(next.getValue() - current.getValue());
	if ((itersCount >= nDims) && (fChange <= fThreshold))
	{
		// The solution is equally good at current point and next point...
		finished = true;

		if (globalVerbosity >= 1)
		{
			std::cout << "Conjugate Gradient Method finshed because value isn't changing: " 
				 << fChange << " <= " << fThreshold << endl;
		}

		return;
	}

    // Stop Condition: If we move less than epsilon in parameter space
    Vector dif = current.x - next.x;
//    printf("Checking stop in conjGrad: movement %f, tolerance %f \n", dif.twoNorm(), xTolerance.twoNorm());
	// only stop after searching each dimension once
	if (itersCount >= nDims) {
		bool noChange = true;
		double factor = 1.0; // set to zero to skip this stopping condition - 0.5 to be conservative
		for (int i=0; i < xTolerance.size(); i++) {
			if (fabs(dif(i)) > factor * xTolerance(i))
				noChange = false;
			else {
				// printf("Conjugate Gradient Method  no change at pos %d : %g  tolerance: %g\n", i, dif(i), xTolerance(i));
			}
		}
		if (noChange) {
			if (globalVerbosity >= 1)
			{
				// FIXME: This is true only if both components are not changing much.
				std::cout << "Conjugate Gradient Method finshed because no change in x overall: " << endl;
				for (int i=0; i < xTolerance.size(); i++) 
					std::cout << fabs(dif(i)) << "<" << factor * xTolerance(i) << " ";
				std::cout << endl;
			}

		finished = true;
		return;
		} // no change
	} //final dimension

	// Update the current solution
	current = next;

	// save the previous gradient dir, and calculate a new one
	prevG = g;
#ifdef OPTIMIZATION_VISUALIZER
	globalLogManager.setCategory(0); // calculate gradient;
#endif
	current.setValue(p.computeOneJet(current.x,g));
	p.evaluate(current.x);

	if ( !brentSucceeded ) {
		// Brent Linear method failed to find a minima along the search
		// direction most likely due to the minima lying out of bounds.
		// As cg iterations rely on the fact that we have achieved the
		// minimum along a search direction, we should restart with
		// the gradient being computed at this new direction.
		h	= g;
		// ... also reset itersCount
		// FIXME: What is correct - reset itersCount or don't count this iteration?
		itersCount	= 0;
	}
	else if ((itersCount % nDims)== 0) {
		// After going through n dimensions, I want to start over with a gradient
		h	= g;
	}
	else {
		//otherwise calculate the next direction to search

		// The next direction search is 
		// g + gamma * h
		// where g is the new gradient, h is the previous search direction
		// and gamma is defined as  dgg/gg (with dgg and gg) defined below
		double gg = prevG.dotProduct(prevG);
		double dgg = g.dotProduct(g) - g.dotProduct(prevG);		
		// for Fletcher-Reeves algorithm, set dgg to be g.dotProduct(g)

		// Stop Condition: If gradient is zero, bail
		if(gg == 0.0) {
			
			if (globalVerbosity >= 1)
				std::cout << "Conjugate Gradient Method finshed because gradient is zero" << endl;

			finished = true;
			return;
		}

		double gam = dgg/gg;

		// Update the vectors
		// h = g + gamma * h
		// I do this in two steps to avoid creating a new vector
		h *= gam ;
		h += g ;

		// for Fletcher-Reeves, h = -g + gamma h

	} // else (iteratons % nDims != 0) 

}


// End namespace
NAMESPACE_PAULY_END
