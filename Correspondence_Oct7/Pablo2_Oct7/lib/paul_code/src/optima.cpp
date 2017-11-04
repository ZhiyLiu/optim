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
 * optima.cpp
 *	----------
 * Header declares base classes in the library
 ******************************************************************/
#include "optima.h"
#include <iostream>
#include <iomanip>
#include <math.h>
using namespace std;

extern int globalVerbosity;

/**
 * How to chose this tolerance value
 * (later, make this a parameter)
 * When we step by epsilon, the mahalnobis distance changes
 * by 2 * epsilon * PGAFactor (comes from derivative of the
 * objective function)
 * This component should always be there in the derivative
 * of the objective function, however since mahalanobis
 * distance is not the only term in there and the change
 * here may be cancelled out by another term in the objective
 * function, this may be misleading. If Mahal dist
 * is weighted by something else, then use that weight.
 *
 * Similarly, if there is any term in the objective function
 * that is computed directly from the vector the optimizer is
 * working on, it should be similarly taken into account.
 *
 * LATER: Somehow, transfer just the components of the
 * objective function that are not like a shape mahalanobis
 * penalty.
 *
 * Note this is used for checking df/dx and not just df.
 */
#define OBJ_FUNCTION_EPSILON 0.01


// Begin namespace
NAMESPACE_PAULY_START

/**
 * Given an initial point x, a direction dir, it returns
 * the maximum step one can take in sign(t)*dir direction
 * w/o violating any bounds.
 */
double BoundsFunction::limit(const Vector& x, const Vector& _dir,
	const double& t )
{
	const Vector dir = (t < 0.0) ? (const Vector&)(-_dir) : (_dir);
	double minT	= 1e20;
	assert( x.size() == _dir.size() );
	for( int i = 0; i != x.size(); ++i )
	{
		const double lim	= bound(i);
		const double step	= factor(i)*dir(i);
		double val;
		if( step > 0.0 ) {
			val	= (lim - x(i)) / step;
		}
		else {
			val	= (-lim - x(i)) / step;
		}
		val	= (t < 0.0) ? -val : val;
		if( val < minT ) {
			minT	= val;
		}
	}
	return minT;
}

/**
 * The evaluate function for this class works the following way
 * Given a direction v, give the maximum step that one can take along
 * v or -v in units of |v| so that we do not go out of bounds.
 * The formulation implicitly assumes that the steps are taken in
 * the direction of v from the origin.
 */
double BoundsFunction::evaluate(const Vector& v)
{
	int i;
	double val;
	double minValue	= 1.0e20;
	for( i = 0; i < v.size(); i++ ) {
		val	= bound(i) / (factor(i) * fabs(v(i)));
		if(val < minValue) {
			minValue	= val;
		}
	}
	return minValue;
}

/******************************************************************
 Gaussian Solution Space
 ******************************************************************/
// Create a Gaussian solution space
GaussianSS::GaussianSS(const Vector &mean,const Vector &stdev) {
	this->mean = mean;
	this->stdev = stdev;
}

// Get a feasible solution
Vector GaussianSS::getFeasibleSolution() {
	Vector x(mean.size());

   for(int i=0;i<mean.size();i++) {
      x(i) = getGaussianRnd(mean(i),stdev(i));
   }

	return x;
}

// Check the feasibility of the solution (feasibility 0 means that the solution may not 
// exist at a given position
double GaussianSS::getFeasibility(const Vector &x) {
   double feasibility = 1.0;
	for(int i=0;i<mean.size();i++) {
		feasibility *= penaltyFunction(x(i),mean(i),stdev(i));
	}
   return feasibility;
}

/******************************************************************
 Uniform Solution Space
 ******************************************************************/
// Create a uniform solution space
UniformSS::UniformSS(const Vector &center,const Vector &radius) {
	this->center = center;
	this->radius = radius;
}

// Get a feasible solution
Vector UniformSS::getFeasibleSolution() {
	Vector x(center.size());

   for(int i=0;i<center.size();i++) {
      x(i) = center(i) - radius(i) + rand(2*radius(i));
   }

	return x;
}

// Check the feasibility of the solution (feasibility 0 means that the solution may not 
// exist at a given position
double UniformSS::getFeasibility(const Vector &x) {
	for(int i=0;i<center.size();i++) {
		if(x(i) < center(i) - radius(i) || x(i) > center(i) + radius(i))
			return 0;
	}

	return 1.0;
}

const int MAX_TRIES_FOR_NONZERO_DERIVATIVE	= 5;

double approxOneJetNoTerms(Function &f,const Vector &x,Vector &dxi,
	const Vector &epsili, const double fx)
{
	double fPlus, fMinus;

	// The gradient vector
	dxi.setSize(x.size());

	// Compute each partial derivative
	Vector x1 = x;
	for(int i=0;i<x.size();i++) {
		x1(i) += 0.5 * epsili(i);
		fPlus = f.evaluate(x1);
		x1(i) -= epsili(i);
		fMinus = f.evaluate(x1);
		x1(i) += 0.5 * epsili(i);

		dxi(i) = (fPlus - fMinus) / epsili(i);
	}

	return fx;
}


int numDirectionIgnoredComplaints	= 20;

/******************************************************************
 Jet approximators
 ******************************************************************/
double approxOneJet(Function &f,const Vector &x,Vector &dxi,const Vector &epsili)
{
	// The function value
	Components fx = f.evaluateTerms(x);
	if( !fx.roughSet ) {
		return approxOneJetNoTerms(f, x, dxi, epsili, fx.smooth);
	}
	Components fPlus, fMinus;
	int trials;

	Vector adaptiveEpsili	= epsili / 2;

	// The gradient vector
	dxi.setSize(x.size());

	// Compute each partial derivative
	Vector x1 = x;
	for(int i=0;i<x.size();i++) {
		for( trials = 0; trials < MAX_TRIES_FOR_NONZERO_DERIVATIVE; ++trials ) {
			adaptiveEpsili(i) = adaptiveEpsili(i) * 2;
			x1(i) += 0.5 * adaptiveEpsili(i);
			fPlus = f.evaluateTerms(x1);
			x1(i) -= adaptiveEpsili(i);
			fMinus = f.evaluateTerms(x1);
    	    x1(i) += 0.5 * adaptiveEpsili(i);

			if( fabs(fPlus.rough - fMinus.rough) < 1e-9 && fabs(fPlus.smooth - fMinus.smooth) < 1e-9) {
				// The factor for this direction is set to zero. So don't try any adaptive epsilon
				// business.
				if( numDirectionIgnoredComplaints-- > 0 ) {
					cout << "Warning: |df(x)_" << i << "| = 0. Ensure that the corresponding factor is set to zero." << endl;
					if( numDirectionIgnoredComplaints == 0 ) {
						cout << "(From here on, the above warning message will be suppressed.)" << endl;
					}
				}
				break;
			}
			//
			// If rough derivative > epsilon, then we have a good derivative,
			// else if rough derivative <= epsilon, and smooth derivative is also <= epsilon (and
			// smooth is set), then we are at optimum, break.
			//
			if( fabs(fPlus.rough - fMinus.rough) > OBJ_FUNCTION_EPSILON * adaptiveEpsili(i) ) {
				break;
			}
		}
		if( trials == MAX_TRIES_FOR_NONZERO_DERIVATIVE ) {
//			if( globalVerbosity >= 1 )
				cout << "WARNING: max trials reached, |df(x)_" << i << "| = " << fabs(fPlus.rough - fMinus.rough) << " < " << OBJ_FUNCTION_EPSILON * adaptiveEpsili(i) << "\n";
		}

		dxi(i) = (fPlus.rough + fPlus.smooth - fMinus.rough - fMinus.smooth) / adaptiveEpsili(i);
	}

	if(globalVerbosity >= 1 ) {
		cout << "df(x) = ";
		for(int i=0;i<x.size();i++) {
			cout << setw(8) << dxi(i) << ' ';
		}
		cout << endl;
	}

	return fx.rough + fx.smooth;
}

double approxOneJetNoTerms(Function &f,const Vector &x,Vector &dxi,
	double epsilon, double fx)
{
	double fPlus, fMinus;

	// The gradient vector
	dxi.setSize(x.size());	

	// Compute each partial derivative
	Vector x1 = x;
	for(int i=0;i<x.size();i++) {
		x1(i) += 0.5 * epsilon;
		fPlus = f.evaluate(x1);
		x1(i) -= epsilon;
		fMinus = f.evaluate(x1);
		x1(i) += 0.5 * epsilon;
		dxi(i) = (fPlus - fMinus) / epsilon;
	}

	return fx;
}


double approxOneJet(Function &f,const Vector &x,Vector &dxi,double epsilon)
{
	// The function value
	Components fx = f.evaluateTerms(x);
	if( !fx.roughSet ) {
		return approxOneJetNoTerms(f, x, dxi, epsilon, fx.smooth);
	}
	Components fPlus, fMinus;
	int trials;
	double adaptiveEpsilon;

	// The gradient vector
	dxi.setSize(x.size());	

	// Compute each partial derivative
	Vector x1 = x;
	for(int i=0;i<x.size();i++) {
		adaptiveEpsilon	= epsilon / 2.0;
		for( trials = 0; trials < MAX_TRIES_FOR_NONZERO_DERIVATIVE; ++trials ) {
			adaptiveEpsilon = adaptiveEpsilon * 2;
			x1(i) += 0.5 * adaptiveEpsilon;
			fPlus = f.evaluateTerms(x1);
			x1(i) -= adaptiveEpsilon;
	        fMinus = f.evaluateTerms(x1);
    	    x1(i) += 0.5 * adaptiveEpsilon;
			if( fabs(fPlus.rough - fMinus.rough) < 1e-9 && fabs(fPlus.smooth - fMinus.smooth) < 1e-9) {
				// The factor for this direction is set to zero. So don't try any adaptive epsilon
				// business.
				if( numDirectionIgnoredComplaints-- > 0 ) {
					cout << "Warning: |df(x)_" << i << "| = 0. Ensure that the corresponding factor is set to zero." << endl;
					if( numDirectionIgnoredComplaints == 0 ) {
						cout << "(From here on, the above warning message will be suppressed.)" << endl;
					}
				}
				break;
			}
			//
			// If rough derivative > epsilon, then we have a good derivative,
			// else if rough derivative <= epsilon, and smooth derivative is also <= epsilon (and
			// smooth is set), then we are at optimum, break.
			//
			if( fabs(fPlus.rough - fMinus.rough) > OBJ_FUNCTION_EPSILON * adaptiveEpsilon ) {
				break;
			}
		}
		if( trials == MAX_TRIES_FOR_NONZERO_DERIVATIVE ) {
//			if( globalVerbosity >= 1 )
				cout << "WARNING: max trials reached, |df(x)_" << i << "| = " << fabs(fPlus.rough - fMinus.rough) << " < " << OBJ_FUNCTION_EPSILON * adaptiveEpsilon << "\n";
		}
		dxi(i) = (fPlus.smooth + fPlus.rough - fMinus.smooth - fMinus.rough) / adaptiveEpsilon;
	}

	return fx.smooth + fx.rough;
}

NumericalFunction::NumericalFunction(Function &f,const Vector &epsili) : function(f) {
	this->epsili = epsili;
	epsilon = 0.0;
}

NumericalFunction::NumericalFunction(Function &f,double epsilon) : function(f) {
	assert(epsilon > 0.0);
	this->epsilon = epsilon;
}


double NumericalFunction::computeOneJet(const Vector &v,Vector &dxi) {
	if(epsilon > 0.0) {
		return approxOneJet(function,v,dxi,epsilon);
	}
	else {
		return approxOneJet(function,v,dxi,epsili);
	}
}

double NumericalFunction::computeDirectionalJet(const Vector &x,const Vector &u,double &fu) {
	double eps = epsilon;
	if(eps==0) {
		// Must compute a new epsilon, lets use epsili dot u
		eps = epsili.dotProduct(u);
	}

	// The function value
	double fx = evaluate(x);
	double fX = evaluate(x + eps*u);
	fu = (fX-fx) / eps;

	return fx;
}


/******************************************************************
 Solution stuff
 ******************************************************************/
Solution::Solution() {
	value = 0;
   hasValue = false;
}

Solution::Solution(const Solution &copy) {
	setValue(copy.value);
}

NumericalSolution::NumericalSolution(int dims) : Solution(),x(dims) {
	assert(dims > 0);
}

NumericalSolution::NumericalSolution(const NumericalSolution &ns) : Solution(ns),x(ns.x) {
}

NumericalSolution::NumericalSolution(const Vector &inX) : Solution(),x(inX) {
}

// End namespace
NAMESPACE_PAULY_END

