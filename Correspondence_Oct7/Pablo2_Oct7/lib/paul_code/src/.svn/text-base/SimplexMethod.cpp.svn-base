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
#include <stdlib.h>
#include <memory.h>
#include <math.h>
#include "SimplexMethod.h"
#include "support.h"

// Begin namespace
NAMESPACE_PAULY_START

// Reflect function used a lot in the method
void reflect(Vector &victim,const Vector &centroid,double factor) {
	victim = (centroid+(centroid-victim)*factor);
}

// Performs an iteration of the simplex method.  The input solutions
// must be evaluated already.  Returns the index of the lowest solution
// CODE from NRC
/*
int performSimplexIteration(Function &problem,NumericalSolution p[]) {
   int i,ihi,ilo,inhi,j;
   int ndim=p[0].x.size();
   int mpts=ndim+1;
   double sum,swap,ysave,ytry;
   Vector psum(ndim);

   // Compute the sum of the points
   for(i=0;i<mpts;i++) {
      psum += p[i].x;
   }

   // This section of the code computes the highest, 2nd highest and lowest
   // priced solutions
   ilo = 0;
   ihi = p[0].getValue() > p[1].getValue() ? (inhi=1,0) : (inhi=0,1);
   for(i=0;i<mpts;i++) {
      if(p[i].getValue() <= p[ilo].getValue())
         ilo = i;
      if(p[i].getValue() > p[ihi].getValue()) {
         inhi = ihi;
         ihi = i;
      }
      else if(p[i].getValue() > p[inhi].getValue() && i!=ihi)
         inhi = i;
   }

   // Extrapolate by the factor of -1
   ytry = amotry(problem,p,psum,ihi,-1.0);
   if(ytry <= y[ilo]) {
      // Try extrapolation
      ytry = amotry(problem,p,psum,ihi,2.0);
   }
   else if(ytry >= y[inhi]) {
      // Try simple contraction
      ysave = p[ihi].getValue();
      ytry = amotry(problem,p,psum,ihi,0.5);
      if(ytry >= ysave) {
         // Must perform complete contraction
         for(i=0;i<mpts;i++) {
            if(i != ilo) {
               p[i].x = psum = 0.5 * (p[i].x+p[ilo].x);
               p[i].setValue(problem.evaluate(psum));
            }
         }
      }
   }
}
*/

SimplexMethod::SimplexMethod(Function &problem,SolutionSpace &space) :
p(problem),ss(space),bestEver(space.dimensions())
{
	type = Simplex;

	setReflectionConstants(1.0,2.0,0.5,-0.5);
   simp = ss.dimensions()+1;

   trials = new NumericalSolution *[simp];
   for(int i=0;i<simp;i++) {
      if(i==0) 
         trials[i] = new NumericalSolution(ss.getMean());
      else 
         trials[i] = new NumericalSolution(ss.getFeasibleSolution());
   }

   // First iteration has not been performed
   initialized = false;

	// Random replacement disabled by default
	prr = 0.05;
   minAge = 5;
}

SimplexMethod::~SimplexMethod() {
   for(int i=0;i<simp;i++)
      delete trials[i];
   delete[] trials;
}

void SimplexMethod::setReflectionConstants(double r,double e,double nc,double pc) {
	rc = r;
   ec = e;
   ncc = nc;
   pcc = pc;
}

int compareNumericalSolutions(const void *s1,const void *s2) {
	// UGHhhh
	NumericalSolution *p1 = (NumericalSolution *)*((NumericalSolution **)(s1));
	NumericalSolution *p2 = (NumericalSolution *)*((NumericalSolution **)(s2));

	return (p1->value==p2->value) ? 0 :
	(p1->value < p2->value) ? -1 : 1;
}

void SimplexMethod::performIteration() {
   int i;

   if(!initialized) {
      // Rank trials by value - I use qsort because it's fast enough
		qsort(trials,simp,sizeof(NumericalSolution *),compareNumericalSolutions);
      best = 0;
      worst = simp-1;
      secondWorst = simp-2;

      // Compute the best-ever solution
      bestEver = *trials[best];

      // Now we are initialized
      initialized = true;
   }

	// Find centroid of all but the worstt point
   Vector centroid(simp-1);
   for(i=0;i<simp;i++) {
		if(i!=worst)
      	centroid = centroid + trials[i]->x;
   }
   centroid /= (simp-1);

   // Perform reflection
   NumericalSolution ref(trials[worst]->x);
   reflect(ref.x,centroid,rc);
   ref.setValue(p.evaluate(ref.x));

   // Compare reflection to best
   if(trials[best]->value > ref.value) {
   	// This is an improvement.  Consider an expansion
		NumericalSolution expan(trials[worst]->x);
      reflect(expan.x,centroid,ec);
      expan.setValue(p.evaluate(expan.x));
      *trials[worst] = (expan.value < ref.value) ? expan : ref;
   } 

   else if(trials[secondWorst]->value > ref.value) {
    	// Reflection better than last two vertices but not good enough for expansion
      *trials[worst] = ref;
   }

   else if(trials[worst]->value > ref.value) {
   	// Solution is better than the worst, but marginally.  Try a positive contraction
		NumericalSolution pc(trials[worst]->x);
      reflect(pc.x,centroid,pcc);
      pc.setValue(p.evaluate(pc.x));
      *trials[worst] = (ref.value < pc.value) ? ref : pc;
   }
   else {
		// Reflection is strictly worse than all other solutions.  Force negative contraction
		NumericalSolution nc(trials[worst]->x);
      reflect(nc.x,centroid,ncc);
      nc.setValue(p.evaluate(nc.x));
      *trials[worst] = nc;
   }

   // Find a place for the new solution
   int newPosition = 0;
   while(newPosition < simp-2) {
      if(trials[newPosition]->value > trials[worst]->value) 
         break;
      newPosition++;
   }

   // Shift down the solutions and insert the new one
   NumericalSolution *tmp = trials[worst];
   for(i=simp-1;i>newPosition;i--) {
      trials[i] = trials[i-1];
   }
   trials[newPosition] = tmp;

	// Remember best ever
	if(bestEver.value > trials[best]->value) {
	   bestEver = *trials[best];
   }
}

void SimplexMethod::setSolution(Vector v,int i) {
	trials[i]->x = v;
   trials[i]->setValue(p.evaluate(v));
}

// End namespace
NAMESPACE_PAULY_END

