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
 * optima.h
 *	--------
 * Header declares base classes in the library
 ******************************************************************/
#ifndef _OPTIMA_CLASSES_
#define _OPTIMA_CLASSES_

// Include files
#include "libs.h"
#include "matrix.h"
#include "support.h"
#include <iostream>
using namespace std;

// Begin namespace
NAMESPACE_PAULY_START


/**
 * Divides a function value into two components
 * - one that is smooth wrt to the argument
 * - and another that is not.
 * This aids in derivative computation.
 */
class Components
{
  public:
	double smooth;
	double rough;
	bool smoothSet;
	bool roughSet;

	Components() : smooth(0.0), rough(0.0), smoothSet(true), roughSet(false) {}

	// A helper constructor
	explicit Components(const double& _smooth)
		: smooth(_smooth), rough(0.0), smoothSet(true), roughSet(false) {}

	Components( const double& _smooth, const double& _rough)
		: smooth(_smooth), rough(_rough), smoothSet(true), roughSet(true) {}

	friend std::ostream& operator << (std::ostream& out, const Components& c)
	{
		out << "( smooth: ";
		if( c.smoothSet ) {
			out << c.smooth;
		}
		else {
			out << "-";
		}
		out << ", rough: ";
		if( c.roughSet ) {
			out << c.rough;
		}
		else {
			out << "-";
		}
		out << " )";
		return out;
	}

	enum FunctionClass_t {
		SMOOTH = 0,
		ROUGH = 1,
	};
};

/******************************************************************
 Problems
 ******************************************************************/

/**
 * A function is an optimization problem defined on Euclidean N-space.
 * Many problems can map onto the space of euclidean problems
 */
class Function {
protected:
	// Evaluation cost.  It's the functions responsibility to account for
	// evaluation cost, but we provide the variable for usability
	double evaluationCost;
public:
	Function() {
		evaluationCost = 0;
	}

	// Evaluate the function at a position
	virtual double evaluate(const Vector &v) = 0;

	// Evaluate the function at a position
	// and return terms that vary smoothly with v
	// and those that do not.
	virtual Components evaluateTerms(const Vector &v) {
		// Default behavior - only smooth terms exist.
		return Components(evaluate(v));
	}

	// Get and set evaluation cost.  For functions that are wrappers around
	// other functions these methods should be overridden
	virtual void setEvaluationCost(double c) {
		evaluationCost = c;
	}
	virtual double getEvaluationCost() {
		return evaluationCost;
	}

	virtual ~Function() {
	};
};

class BoundsFunction : public Function
{
  public:
	/**
	 * The bound to be used for the ith component of the vector.
	 * This bound is assumed to be symmetric about 0.
	 */
	virtual double bound(int i)
	{
		assert(i >= 0);
		return bounds( (i < bounds.size()) ? i : (bounds.size() - 1) );
	}
	/**
	 * The factor for the ith component of the vector.
	 */
	virtual double factor(int i)
	{
		assert(i >= 0);
		return factors( (i < factors.size()) ? i : (factors.size() - 1) );
	}
	/**
	 * Given an initial point x, a direction dir, it returns
	 * the maximum step one can take in sign(t)*dir direction
	 * w/o violating any bounds.
	 */
	virtual double limit(const Vector& x, const Vector& dir,
		const double& t );

	/**
	 * The evaluate function for this class works the following way
	 * Given a direction v, give the maximum step that one can take along
	 * v or -v in units of |v| so that we do not go out of bounds.
	 * The formulation implicitly assumes that the steps are taken in
	 * the direction of v from the origin.
	 */
	double evaluate(const Vector& v);
  protected:
	/**
	 * Usually the default implementations provided up suffice.
	 * Derived classes should set the values in these arrays
	 * when they are constructed.
	 */
	Vector bounds;
	Vector factors;
};

/**
 * A Gradient problem is an optimization problem defined on Euclidean N-space.
 * for which there must be a way to compute the gradient.
 */
class DifferentiableFunction : public Function {
public:
	// Compute function at x and the partial derivatives at x
	virtual double computeOneJet(const Vector &x,Vector &dxi)=0;

	// Computes directional derivative and function value.  Return function value
	// Defaults to <u,grad f>
	virtual double computeDirectionalJet(const Vector &x,const Vector &u,double &fu) {
		Vector grad;
		double fx = computeOneJet(x,grad);
		fu = grad.dotProduct(u);
		return fx;
	}

   virtual ~DifferentiableFunction() {};
};

/**
 * There will be cases when you want to run a derivative based method on a 
 * function that does not have derivatives.  I provide a method that you
 * can call to do the computation numerically.  The parameter is the 
 * epsilon value to be used in computing the derivatives.
 *
 * Use the first method if you want to specify different epsilons in all
 * directions and the second if you just want to specify on epsilon
 */
double approxOneJet(Function &f,const Vector &x,Vector &dxi,const Vector &epsili);
double approxOneJet(Function &f,const Vector &x,Vector &dxi,double epsilon=1.0e-10);

/**
 * And just to save people from having to implement a dummy computeOneJet everywhere,
 * we let you use a class NumericalFunction to wrap a non differentiable function and
 * make it differentiable
 */
class NumericalFunction : public DifferentiableFunction {
private:
	Function &function;
	double epsilon;
	Vector epsili;
public:
	NumericalFunction(Function &function,const Vector &epsili);
	NumericalFunction(Function &function,double epsilon = 1.0e-10);

	double computeOneJet(const Vector &v,Vector &dxi);
	double evaluate(const Vector &v) {
		return function.evaluate(v);
	}

	// Computes directional derivative and function value.  Return function value
	// Does just one computation for derivative instead of three otherwize required
	// to compute the gradient
	virtual double computeDirectionalJet(const Vector &x,const Vector &u,double &fu);

	// This function does not incur any cost on its own so it
	// asks the client function for the cost accounting information
	virtual void setEvaluationCost(double c) {
		function.setEvaluationCost(c);
	}
	virtual double getEvaluationCost() {
		return function.getEvaluationCost();
	}

   virtual ~NumericalFunction() {};
};


/******************************************************************
 Solutions
 ******************************************************************/

/**
 * A potential solution to any problem
 */
class Solution {
public:
   bool hasValue;
	double value;

   Solution();
   Solution(const Solution &copy);
   virtual ~Solution() {};

	double getValue() {
      return value;
   }

   void setValue(double val) {
      value=val;
      hasValue=true;
   }
};

/**
 * A solution to a numerical problem is a vector
 */
class NumericalSolution : public Solution {
protected:
public:
	Vector x;

	NumericalSolution(int dim);
	NumericalSolution(const Vector &inX);
	NumericalSolution(const NumericalSolution &ns);
};



/******************************************************************
 Solution Spaces
 ******************************************************************/

/**
 * A Solution space is a prior distribution on Rn that describes the feasibility
 * of finding a solution at any point.  A typical solution space used is the Gaussian
 * solution space.  Solution spaces can also describe constraints on the problem. 
 */
class SolutionSpace {
protected:
public:
	// Get the mean of the distribution, also the most feasible solution
	virtual const Vector &getMean() = 0;

	// Get the standard deviations of the feasible space.  These standard deviations
	// describe a metric on the space of potential solutions.
	virtual const Vector &getStandardDeviations() = 0;

	// Get a feasible solution
	virtual Vector getFeasibleSolution() = 0;

	// Check the feasibility of the solution (feasibility 0 means that the solution may not 
	// exist at a given position
	virtual double getFeasibility(const Vector &x) = 0;

	// Get dimentionality
	virtual int dimensions() = 0;

   virtual ~SolutionSpace() {};
};

/**
 * A GaussianSS is a space in which all points in space are feasible solutions, but 
 * some of the solutions are more feasible than others, with the mean being the most
 * feasible solution of them all
 */
class GaussianSS : public SolutionSpace {
protected:
	Vector mean,stdev;
public:
	// Constructor for the solution space
	GaussianSS(const Vector &mean,const Vector &stdev);

	// Get the mean of the distribution, also the most feasible solution
	const Vector &getMean() {
		return mean;
	}

	// Get the standard deviations of the feasible space.  These standard deviations
	// describe a metric on the space of potential solutions.
	const Vector &getStandardDeviations() {
		return stdev;
	}

	// Get a feasible solution
	Vector getFeasibleSolution();

	// Check the feasibility of the solution (feasibility 0 means that the solution may not 
	// exist at a given position
	double getFeasibility(const Vector &x);

	// Get dimentionality
	int dimensions() {
		return mean.size();
	}
};


/**
 * A UniformSS is a space in which all points in a box are feasible solutions,  
 * with all solutions as likely to be optimal as others.  Points outside the
 * box are infeasible
 */
class UniformSS : public SolutionSpace {
protected:
	Vector center,radius;
public:
	// Constructor for the solution space
	UniformSS(const Vector &center,const Vector &radius);

	// Get the mean of the distribution, also the most feasible solution
	const Vector &getMean() {
		return center;
	}

	// Get the standard deviations of the feasible space.  These standard deviations
	// describe a metric on the space of potential solutions.
	const Vector &getStandardDeviations() {
		return radius;
	}

	// Get a feasible solution
	Vector getFeasibleSolution();

	// Check the feasibility of the solution (feasibility 0 means that the solution may not 
	// exist at a given position
	double getFeasibility(const Vector &x);

	// Get dimentionality
	int dimensions() {
		return center.size();
	}
};


/******************************************************************
 Methods
 ******************************************************************/

/**
 * A way to solve a problem 
 */
class Method {
public:
	enum method_t { ConjugateGradient, Evolutionary, Simplex };

	virtual double getBestEverValue() = 0;
	virtual ~Method() {};
	method_t methodType() { return type; }

protected:
	method_t type;
};

/**
 * A numerical method is a method that deals with numerical problems
 */
class NumericalMethod : public Method {
public:
	virtual const Vector &getBestEverX() = 0;

	// Perform an iteration of the method
	virtual void performIteration() = 0;

	// Check if we are finished (some algorithms may have termination conditions
	// but most dont so we default by returning false
	virtual bool isFinished() {
		return false;
	}

	// A method always operates on a function, so we should be able to return
	// the current function
	virtual Function &getFunction() = 0;
};

/**
 * A list of defined methods as integers.  Use these constants in switch statements
 * to run different methods
 */
#define OPTIMA_CGM			1							// Conjugate Gradient Method
#define OPTIMA_SIMPLEX		2							// Simplex Method
#define OPTIMA_ES				3							// Evolutionary strategy
#define OPTIMA_BLM			4							// Brent linear method

#define OPTIMA_ConjugateGradientMethod		OPTIMA_CGM
#define OPTIMA_SimplexMethod					OPTIMA_SIMPLEX
#define OPTIMA_EvolutionaryStrategy			OPTIMA_ES
#define OPTIMA_BrentLinearMethod				OPTIMA_BLM

// End namespace
NAMESPACE_PAULY_END

//---------------------------------------------------------------------------

#endif	/* _OPTIMA_CLASSES_ */

