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
 * EvolutionaryStrategy.h
 *	-------------------------
 * This method from Michaelewitz allows us to optimize using 
 * genetic-type algorithms
 ******************************************************************/
#ifndef _OPTIMA_EVOLUTIONARY_STRATEGY_
#define _OPTIMA_EVOLUTIONARY_STRATEGY_

// Include files
#include "optima.h"

// Begin namespace
NAMESPACE_PAULY_START


// There are two different selection modes, m+l and m,l.
#define SELECTION_MuCommaLambda 0
#define SELECTION_MuPlusLambda 1
#define SELECTION_Invalid 2

/**
 * A solution for ES has a sigma vector as well as the x vector.
 */
class ESSolution : public NumericalSolution
{
public:
	Vector sigma;
	ESSolution(int dims);
	ESSolution(const ESSolution &copy);
   ESSolution(const Vector &x,const Vector &sigma);
	virtual ~ESSolution();
};

/*******************************************************************
 Optimize f(x) using the simplex method 

 Usage:
	-	Create an EvolutionStrategy (with a function and a sample x)
	-	Set options with setXXX methods
	-  Call performIteration() until happy 
	-	Check results with getBestEverXXX()
 ******************************************************************/
class EvolutionaryStrategy : public NumericalMethod
{
protected:
   // Function that we are optimizing
   Function &problem;

   // Solution space that defines the range of possible solutions and possible penalty
   // to be assessed
	SolutionSpace &ss;

   // A collection of solutions
   ESSolution **mlArray;

	// Parameters of the method
   int mu,lambda,mode;

   // Probability of random solution being introduces into the population
	double pRandSolution;

   // Best ever solution
   ESSolution bestEver;

   // Upper and lower bounds on x's and sigmas of the method
	Vector xUpperBound,xLowerBound;
   Vector sUpperBound,sLowerBound;
   bool xBounded,sBounded;

   // Factor that relates the sigma's of the method with the sigmas of the solution space
	double sigmaFactor;

   // Compute a random sigma vector based on the solutions space, and sigmaFactor 
   void computeRandomSigma(Vector &sigma);

   // This defines how much the sigmas change each iteration
   Vector deltaSigma;

	// Are we initialized?
   bool initialized;

	void setBestEver(ESSolution *);

	// virtual void crossover(ESSolution *s1,ESSolution *s2,ESSolution *child);
	virtual void mutate(ESSolution *s1);
	virtual void selection();
	virtual void makeRandomSolution(ESSolution *s);

   // I have implemented support for multi-parent crossover with custom
   // crossover functions.  This will be interesting
   int nCrossoverParents;
   void (*crossover)(Function &f,int nParents,int *parentIndex,ESSolution **parentArray,ESSolution *child);
   // void (mutate*)(ESSolution *s);

public:
	// Create a new evolutionary strategy
	EvolutionaryStrategy(Function &p,SolutionSpace &ss,int inMu,int inLambda,int mode);
	virtual ~EvolutionaryStrategy();

	// Set the sigma factor. The factor determines the proportion between the sigma
	// vector reported by SolutionSpace and upper bound on initial mutability of 
	// solutions in the population.  Default value of this factor is 0.1
	void setSigmaFactor(double sf) {
		sigmaFactor = sf;
	}

	// Set Nth solution while the algorithm is running
	void setNth(int n,const Vector &x);
   void setNth(int n,const Vector &x,const Vector &sigma);

	// Vector of deltaSigmas determines how much mutation the delta sigmas undergo
	void setDeltaSigma(const Vector &ds);

	// Set the probability that a random solution is introduced to the population
	void setPRandSolution(double p) {
		pRandSolution = p;
	}

	// Set the bounds on the solution space.  The method checks against this bounding box.
	// Bounds can also me manipulated via the function and the solution space.
   void setXBounds(const Vector &xLower,const Vector &xUpper) {
      xLowerBound = xLower;
      xUpperBound = xUpper;
      xBounded = true;
   }
   void setSigmaBounds(const Vector &sigmaLower,const Vector &sigmaUpper) {
      sLowerBound = sigmaLower;
      sUpperBound = sigmaUpper;
      sBounded = true;
   }

	// Peform an iteration of the algorithms
	void performIteration() {
		selection();
	}

	// These methods let us get the best values from the populations
   double getBestEverValue() {return getBestEver()->getValue();}
   const Vector &getBestEverX() {return getBestEver()->x;}

	ESSolution *getBestEver();
	ESSolution *getBest();
	ESSolution *getNth(int n);

	// Get the problem on which we are operating
	virtual Function &getFunction() {
		return problem;
	}

   // Force use of a custom crossover function
   // The first parameter is the number of parents that undergo crossover and the second
   // is the function that performs the crossover
   void useCustomCrossover(int nParents,void (*xo)(Function &,int,int *,ESSolution **,ESSolution *)) {
      assert(nParents <= mu);
      nCrossoverParents = nParents;
      crossover = xo;
   }
};

void discreteCrossover(Function &,int,int *,ESSolution **,ESSolution *);
void intermediateCrossover(Function &,int,int *,ESSolution **,ESSolution *);
void simplexCrossover(Function &,int,int *,ESSolution **,ESSolution *);
void noCrossover(Function &f,int nParents,int *index,ESSolution **array,ESSolution *child);

// End namespace
NAMESPACE_PAULY_END

#endif	/* _OPTIMA_EVOLUTIONARY_STRATEGY_ */

