//---------------------------------------------------------------------------
//#include <vcl.h>
//#pragma hdrstop

#include <stdlib.h>
#include <memory.h>
#include <math.h>
#include "simplex.h"
#include "support.h"

void SimplexProblem::reflect(Vector &victim,Vector &centroid,double factor) {
	victim = (centroid+(centroid-victim)*factor);
}

void SimplexProblem::randomSimplexX(Vector &v) {
	randomX(v);
}



Simplex::Simplex(SimplexProblem &problem) :
p(problem)
{
	setReflectionConstants(1.0,2.0,0.5,-0.5);
   simp = p.getDimensions()+1;
   trials = new SimplexSolution *[simp];
	memset(trials,0,sizeof (SimplexSolution *)*simp);
   initialized = 0;
   prr = 0.05;
   minAge = 5;
}

void Simplex::setReflectionConstants(double r,double e,double nc,double pc) {
	rc = r;
   ec = e;
   ncc = nc;
   pcc = pc;
}

int compareSimplexSolutions(const void *s1,const void *s2) {
	// UGHhhh
	SimplexSolution *p1 = (SimplexSolution *)*((SimplexSolution **)(s1));
	SimplexSolution *p2 = (SimplexSolution *)*((SimplexSolution **)(s2));

	return (p1->value==p2->value) ? 0 :
	(p1->value < p2->value) ? -1 : 1;
}

void Simplex::performIteration() {
   int i;

   if(!initialized) {
   	for(int i=0;i<simp;i++) {
      	if(trials[i]==NULL) {
         	trials[i] = new SimplexSolution(p.getDimensions());
            p.randomSimplexX(trials[i]->x);
            trials[i]->setValue(p.evaluate(trials[i]->x));
         }
      }
      initialized = true;

      // Rank trials by value - I use qsort because it's fast enough
		qsort(trials,simp,sizeof(SimplexSolution *),compareSimplexSolutions);
      best = 0;
      worst = simp-1;
      secondWorst = simp-2;

      /*
      SimplexSolution *origin = new SimplexSolution(p.getDimensions());
      origin->setValue(p.evaluate(origin->x));

      if(origin->getValue() < trials[best]->getValue()) {
         bestEver = origin;
      }
      else {*/
         bestEver = new SimplexSolution(trials[best]);
      /*   delete origin;
      }
      */
   }
   else {
	   // If required, introduce random solution in place of existing solution
      // based on age.
      // ...
   }

	// Find centroid of the vectors
   Vector centroid(simp-1);
   for(i=0;i<simp;i++) {
		if(i!=worst)
      	centroid = centroid + trials[i]->x;
   }
   centroid /= (simp-1);

   // Perform reflection
   SimplexSolution ref(simp-1);
   ref.x = trials[worst]->x;
   p.reflect(ref.x,centroid,rc);
   ref.setValue(p.evaluate(ref.x));

   // Compare reflection to best
   if(trials[best]->value > ref.value) {
   	// This is an improvement
		SimplexSolution expan(simp-1);
      expan.x = trials[worst]->x;
      p.reflect(expan.x,centroid,ec);
      expan.setValue(p.evaluate(expan.x));
      if(expan.value < ref.value) {
      	// Expansion is an improvement
         delete trials[worst];
         trials[worst] = new SimplexSolution(&expan);
      }
      else {
      	// Reflection beats expansion
         delete trials[worst];
         trials[worst] = new SimplexSolution(&ref);
      }
   } // if reflection better that best vertex

   else if(trials[secondWorst]->value > ref.value) {
    	// Reflection better than last two vertices
      delete trials[worst];
      trials[worst] = new SimplexSolution(&ref);
   }

   else if(trials[worst]->value > ref.value) {
   	// Positive contraction
		SimplexSolution pc(simp-1);
      pc.x = trials[worst]->x;
      p.reflect(pc.x,centroid,pcc);
      pc.setValue(p.evaluate(pc.x));
      delete trials[worst];
      trials[worst] = new SimplexSolution(&pc);
   }
   else {
		// Negative contraction
		SimplexSolution nc(simp-1);
      nc.x = trials[worst]->x;
      p.reflect(nc.x,centroid,ncc);
      nc.setValue(p.evaluate(nc.x));
      delete trials[worst];
      trials[worst] = new SimplexSolution(&nc);
   }

   // Update rankings
   SimplexSolution *swap = trials[worst];
   trials[worst] = trials[secondWorst];
   trials[secondWorst] = swap;

   for(i=0;i<simp-1;i++) {
    	trials[i]->age++;
		if(trials[i]->age > minAge) {
        	if(rand(1.0) < prr) {
            delete(trials[i]);
 	      	trials[i] = new SimplexSolution(p.getDimensions());
            p.randomSimplexX(trials[i]->x);
     	      trials[i]->setValue(p.evaluate(trials[i]->x));
         }
      }
   }

 	qsort(trials,simp-1,sizeof(SimplexSolution *),compareSimplexSolutions);

	// Remember best ever
	if(bestEver->value > trials[best]->value) {
   	delete bestEver;
	   bestEver = new SimplexSolution(trials[best]);
   }
}

void Simplex::setSolution(Vector v,int i) {
	delete trials[i];
   trials[i] = new SimplexSolution(v);
   trials[i]->setValue(p.evaluate(v));
}











//---------------------------------------------------------------------------

