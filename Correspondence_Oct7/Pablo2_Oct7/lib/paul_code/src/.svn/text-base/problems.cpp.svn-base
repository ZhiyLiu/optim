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
 * problems.cpp
 *	------------
 * Sample problems, mostly from Michaelewitz
 ******************************************************************/
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "problems.h"

// Begin namespace
NAMESPACE_PAULY_START

/******************************************************************
 CenterDistancePbm
 ******************************************************************/
CenterDistancePbm::CenterDistancePbm(Vector center) {
   this->center = center;
}

double CenterDistancePbm::evaluate(const Vector &v) {
	double result = 0;
	for(int i=0;i<v.size();i++) {
      double xi = v(i)-center(i);
		result += xi*xi;
	}

	evaluationCost += 1;

	return result;
}

double CenterDistancePbm::computeOneJet(const Vector &x,Vector &dxi) {
	for(int i=0;i<x.size();i++) {
		dxi(i) = 2.0 * (x(i) - center(i));
	}

	evaluationCost += 1;

	return evaluate(x);
}


/******************************************************************
 DeJongF4Pbm
 ******************************************************************/
DeJongF4Pbm::DeJongF4Pbm(Vector center) {
   this->center = center;
}

double DeJongF4Pbm::evaluate(const Vector &v) {
	double result = 0;

	for(int i=0;i<v.size();i++) {
		double xi = v(i)-center(i);
		result += (i+1)*xi*xi;
	}

	evaluationCost += 1;

	return result;
}

double DeJongF4Pbm::computeOneJet(const Vector &x,Vector &dxi) {
	for(int i=0;i<x.size();i++) {
		dxi(i) = i * 2.0 * (x(i) - center(i));
	}

	evaluationCost += 1;

	return evaluate(x);
}


/******************************************************************
 SixHumpCamelBackPbm
 ******************************************************************/
NoiseProblem::NoiseProblem(Function &f,double inSigma) : function(f),sigma(inSigma) {
}

double NoiseProblem::evaluate(const Vector &v) {
	return function.evaluate(v) + rand(sigma);
}



/******************************************************************
 SixHumpCamelBackPbm
 ******************************************************************/
double SixHumpCamelBackPbm::evaluate(const Vector &v) {
	double x = v(0);
	double y = v(1);

	evaluationCost += 1;

	return (4 - 2.1*x*x + x*x*x*x/3)*x*x +x*y+(-4+4*y*y)*y*y;
}

double SixHumpCamelBackPbm::computeOneJet(const Vector &v,Vector &dxi) {
	double x = v(0);
	double y = v(1);

	dxi(0) = 8*x-8.4*x*x*x + 2*x*x*x*x*x + y;
	dxi(1) = x - 8*y + 16*y*y*y;

	evaluationCost += 1;

	return evaluate(v);
}


/******************************************************************
 ShubertPbm
 ******************************************************************/
double ShubertPbm::evaluate(const Vector &v)  {
	double x = v(0);
	double y = v(1);

	double p1=0,p2=0;

	for(int i=1;i<=5;i++) {
		p1 += i*cos((i+1)*x+i);
		p2 += i*cos((i+1)*y+i);
	}

	evaluationCost += 1;

	return p1*p2;
}

double ShubertPbm::computeOneJet(const Vector &v,Vector &dxi) {
	double x = v(0);
	double y = v(1);

	double p1=0,p2=0;
	double p1x=0,p2y=0;

	for(int i=1;i<=5;i++) {
		p1 += i*cos((i+1)*x+i);
		p2 += i*cos((i+1)*y+i);
		p1x += -i*(i+1)*sin((i+1)*x+i);
		p2y += -i*(i+1)*sin((i+1)*y+i);
	}

	dxi(0) = p1x*p2;
	dxi(1) = p1*p2y;

	evaluationCost += 2;

	return p1*p2;
}

// End namespace
NAMESPACE_PAULY_END

