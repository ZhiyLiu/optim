//---------------------------------------------------------------------------
#ifndef _OPTIMA_CLASSES_
#define _OPTIMA_CLASSES_

#include "matrix.h"
#include "support.h"


class Problem {
};

class Solution {
public:
   bool hasValue;
	double value;

   Solution();
   Solution(const Solution &copy);

	double getValue() const {
      return value;
   }

   void setValue(double val) {
      value=val;
      hasValue=true;
   }

};

class Method {
public:
	virtual double getBestEverValue() = 0;
};

class NumericalProblem : public Problem {
public:
	virtual double evaluate(const Vector &v) = 0;
   virtual void randomX(Vector &v) = 0;
   virtual int getDimensions() = 0;
};

class NumericalSolution : public Solution {
protected:
public:
	Vector x;

	NumericalSolution(int dim);
   NumericalSolution(const Vector &inX);
	NumericalSolution(const NumericalSolution &ns);
};

class NumericalMethod : public Method {
	virtual Vector getBestEverX() = 0;
};

//---------------------------------------------------------------------------

#endif	/* _OPTIMA_CLASSES_ */

