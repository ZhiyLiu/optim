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
 * SimplexMethod.h
 *	-------------------------
 * This method similar to NRC-10.4 allows us to optimize nicely 
 * using no derivatives in n-space.  This method is nice but dumb...
 ******************************************************************/
#ifndef _OPTIMA_SIMPLEX_
#define _OPTIMA_SIMPLEX_

// Include files
#include "optima.h"

// Begin namespace
NAMESPACE_PAULY_START

/*******************************************************************
 Optimize f(x) using the simplex method 

 Usage:
	-	Create a SimplexMethod (with a function and a sample x)
	-	Set options with setXXX methods
	-  Call performIteration() until happy 
	-	Check results with getBestEverXXX()
 ******************************************************************/
class SimplexMethod : public NumericalMethod {
protected:
	NumericalSolution **trials;
	NumericalSolution bestEver;

   Function &p;
	SolutionSpace &ss;

   double rc,ec,pcc,ncc;
   double prr;
   int minAge;

   int worst,best,secondWorst;
   int simp;

   bool initialized;

public:
	// Create a new method
	SimplexMethod(Function &problem,SolutionSpace &ss);
   ~SimplexMethod();

	// Reflection constants determine how much the simplex grows or shrinks by depending
	// on the operation that it performs.
   void setReflectionConstants(double reflection,
   									 double expansion,
                               double posContration,
                               double negContration);

   // Random replacement policy (after solutions reach a minimum age, they are expelled with
	// a given probability).
   void setPRandReplacement(double pRep,int minimumAge) {
   	prr = pRep;
      minAge = minimumAge;
   }

	// This method allows you to add a solution to the initial set (or add one along the
	// way).  Typically all solutions but one are pulled out from SolutionSpace with calls
	// to getFeasibleSolution.  One solution (at index 0) is the mean of the feasible 
	// solution space.
   void setSolution(Vector v,int i);

   // Run method for one iteration
   void performIteration();

	// Solution report methods
   double getBestEverValue() {return bestEver.getValue();}
   const Vector &getBestEverX() {return bestEver.x;}

	// Get the problem on which we are operating
	virtual Function &getFunction() {
		return p;
	}
};

// End namespace
NAMESPACE_PAULY_END

#endif	/* _OPTIMA_SIMPLEX_ */

