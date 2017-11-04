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
 * problems.h
 *	----------
 * Sample problems, mostly from Michaelewitz
 ******************************************************************/
#ifndef _OPTIMA_PROBLEMS_
#define _OPTIMA_PROBLEMS_

// Include files
#include "optima.h"

// Begin namespace
NAMESPACE_PAULY_START

#define constraintPenalty 1000000.0

/**
 * A simple test problem - squared distance to center, aka DeJongf1
 */
class CenterDistancePbm : public DifferentiableFunction {
public:
   Vector center;
   CenterDistancePbm(Vector center);

	double evaluate(const Vector &x);
	double computeOneJet(const Vector &x,Vector &dxi);
};

/**
 * DeJongF4 Problem.  Similar to distance to center.
 */ 
class DeJongF4Pbm : public DifferentiableFunction {
public:
   Vector center;
   DeJongF4Pbm(Vector center);

	double evaluate(const Vector &x);
	double computeOneJet(const Vector &x,Vector &dxi);
};


/**
 * A noise envelope for other problems.  Adds white noise with given
 * intensity standard deviation to all points.  Makes it very hard to 
 * find the optimum using traditional methods.
 */ 
class NoiseProblem : public Function {
private:
	double sigma;
	Function &function;
public:
	NoiseProblem(Function &function,double sigma);
	double evaluate(const Vector &v);

	// This function does not incur any cost on its own so it
	// asks the client function for the cost accounting information
	virtual void setEvaluationCost(double c) {
		function.setEvaluationCost(c);
	}
	virtual double getEvaluationCost() {
		return function.getEvaluationCost();
	}
};

/**
 * SixHumpCamelBack Problem.  
 */ 
class SixHumpCamelBackPbm : public DifferentiableFunction {
public:
	double evaluate(const Vector &x);
	double computeOneJet(const Vector &x,Vector &dxi);
};

/**
 * Shubert Problem.  
 */ 
class ShubertPbm : public DifferentiableFunction {
public:
	double evaluate(const Vector &v);
	double computeOneJet(const Vector &x,Vector &dxi);
};

// End namespace
NAMESPACE_PAULY_END

#endif	/* _OPTIMA_PROBLEMS_ */

