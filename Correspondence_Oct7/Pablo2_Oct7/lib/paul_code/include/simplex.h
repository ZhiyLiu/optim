//---------------------------------------------------------------------------
#ifndef _OPTIMA_SIMPLEX_
#define _OPTIMA_SIMPLEX_

#include "Classes.h"

class SimplexProblem : public virtual NumericalProblem {
public:
	virtual void reflect(Vector &victim,Vector &centroid,double factor);
   virtual void randomSimplexX(Vector &x);
};

class SimplexSolution : public NumericalSolution {
public:
	int age;
   SimplexSolution(int dims) : NumericalSolution(dims) {
   	age = 0;
   }
   SimplexSolution(Vector &v) : NumericalSolution(v) {
   	age = 0;
   }
	SimplexSolution(const SimplexSolution *s) : NumericalSolution(*s) {
   	age = 0;
   }

};

class Simplex : public NumericalMethod {
protected:
	SimplexSolution **trials;
	SimplexSolution *bestEver;
   SimplexProblem &p;

   double rc,ec,pcc,ncc;
   double prr;
   int minAge;

   int worst,best,secondWorst;
   int initialized;
   int simp;


public:
	// Setup methods
	Simplex(SimplexProblem &problem);
   void setReflectionConstants(double reflection,
   									 double expansion,
                               double posContration,
                               double negContration);
   void setSolution(Vector v,int i);

   // Optimization methods
   void performIteration();

	// Solution report methods
   double getBestEverValue() {return bestEver->getValue();}
   Vector getBestEverX() {return bestEver->x;}

   // Random replacement probability
   void setPRandReplacement(double pRep,int minimumAge) {
   	prr = pRep;
      minAge = minimumAge;
   }
};

//---------------------------------------------------------------------------

#endif	/* _OPTIMA_SIMPLEX_ */

