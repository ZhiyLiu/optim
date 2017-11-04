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
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <memory.h>
#include <math.h>
#include "EvolutionaryStrategy.h"

// Begin namespace
NAMESPACE_PAULY_START

#undef dmax
double dmax(double a,double b) {
	return (a<b) ? b : a;
}
#undef dmin
double dmin(double a,double b) {
	return (a<b) ? a : b;
}

/******************************************************************
 * ESSolution methods
 ******************************************************************/
ESSolution::ESSolution(int dims) : NumericalSolution(dims),sigma(dims)
{
}

ESSolution::ESSolution(const ESSolution &copy) : NumericalSolution(copy), sigma(copy.sigma)
{
}

ESSolution::ESSolution(const Vector &x,const Vector &sigma) : 
NumericalSolution((int)x.size()) 
{
   this->x = x;
   this->sigma = sigma;
}

ESSolution::~ESSolution()
{
}

/******************************************************************
 * Some crossover functions
 ******************************************************************/
void discreteCrossover(Function &f,
                       int nParents,int *index,ESSolution **array,ESSolution *child) 
{
   int ndim = child->x.size();

   // In case of 2 parents we use random bits for speed
   if(nParents==2) {
      ESSolution *s1 = array[index[0]];
      ESSolution *s2 = array[index[1]];
      for(int i=0;i<ndim;i++) {
         child->x(i) = randBit() ? s1->x(i) : s2->x(i);
         child->sigma(i) = randBit() ? s1->sigma(i) : s2->sigma(i);
      }
   }

   // Choose a value from each parent
   else {
      for(int i=0;i<ndim;i++) {
         int parent = rand(nParents);
         child->x(i) = array[index[parent]]->x(i);
         parent = rand(nParents);
         child->sigma(i) = array[index[parent]]->sigma(i);
      }
   }
}

void intermediateCrossover(Function &f,
                           int nParents,int *index,ESSolution **array,ESSolution *child) 
{
   int ndim = child->x.size();

   // In case of 2 parents we use random bits for speed
   if(nParents==2) {
      ESSolution *s1 = array[index[0]];
      ESSolution *s2 = array[index[1]];
      for(int i=0;i<ndim;i++) {
         child->x(i) = 0.5 * (s1->x(i) + s2->x(i));
         child->sigma(i) = 0.5 * (s1->sigma(i) + s2->sigma(i));
      }
   }

   // Choose a value from each parent
   else {
      for(int i=0;i<ndim;i++) {
         child->x(i) = 0;
         child->sigma(i) = 0;
         for(int parent = 0;parent<nParents;parent++) {
            child->x(i) += array[index[parent]]->x(i);
            child->sigma(i) += array[index[parent]]->sigma(i);
         }
         child->x(i) /= nParents;
         child->sigma(i) /= nParents;
      }
   }
}

void noCrossover(Function &f,
                 int nParents,int *index,ESSolution **array,ESSolution *child) 
{
   child->x = array[index[0]]->x;
   child->sigma = array[index[0]]->sigma;
}

void simplexCrossover(Function &f,
                      int nParents,int *index,ESSolution **array,ESSolution *child) 
{
   int i;

   // Simplex must have n+1 points
   int dim = child->x.size();
   dassert(nParents==dim+1);

   // Find the worst solution
   int worst = 0;
   for(i=1;i<nParents;i++) {
      if(array[index[i]]->getValue() > array[index[worst]]->getValue()) {
         worst = i;
      }
   }

   // Find the centroid of other solutions
   Vector centroid(dim),cSigma(dim);
   for(i=0;i<nParents;i++) {
      if(worst!=i) {
         centroid += array[index[i]]->x;
         cSigma += array[index[i]]->sigma;
      }
   }
   centroid /= dim;
   cSigma /= dim;

   // OKay, now perform the reflection
   child->x = centroid*2-array[index[worst]]->x;
   child->sigma = centroid*2-array[index[worst]]->sigma;
}

/******************************************************************
 * EvolutionaryStrategy methods
 ******************************************************************/
EvolutionaryStrategy::EvolutionaryStrategy(Function &p,SolutionSpace &space,
														 int inMu,int inLambda, int inMode) :
problem(p),ss(space),bestEver(ss.dimensions())
{
	int i;

	type = Evolutionary;

	assert(inMu > 0 && inLambda > 0 && ss.dimensions()>0);
	assert(inMode	< SELECTION_Invalid);

	mu = inMu;
	lambda = inLambda;
	mode = inMode;
   sigmaFactor = 0.2;
   pRandSolution = 0.0;
   initialized = false;

	// Initialize the solution array
   mlArray = new ESSolution *[mu+lambda];
   for(i=0;i<mu+lambda;i++) {
      mlArray[i] = new ESSolution(ss.dimensions());
   }

   // Compute best ever solution
   bestEver = *mlArray[0];
   for(i=1;i<mu;i++) 
      setBestEver(mlArray[i]);

	// Preset the delta sigma vector
   deltaSigma.setSize(ss.dimensions());	
   deltaSigma.setAll(1.0);

   // Set the upper and lower bounds
   xBounded = sBounded = false;

   // Set crossover function
   useCustomCrossover(2,discreteCrossover);
}

EvolutionaryStrategy::~EvolutionaryStrategy()
{
	for(int i=0;i<mu+lambda;i++)
		delete mlArray[i];
   delete[] mlArray;
}

void EvolutionaryStrategy::computeRandomSigma(Vector &sigma) {
   Vector sigmaBase = ss.getStandardDeviations()*sigmaFactor;
   for(int i=0;i<ss.dimensions();i++) {
		sigma(i) = sigmaBase(i) * fastExp(getGaussianRnd(0,1.0));
	}
}

void EvolutionaryStrategy::makeRandomSolution(ESSolution *s) {
	s->x = ss.getFeasibleSolution();
   computeRandomSigma(s->sigma);
   s->setValue(problem.evaluate(s->x));
}

// Perform the mutation step
void EvolutionaryStrategy::mutate(ESSolution *s1)
{
   int i;

   // Mutate sigma
   for(i=0;i<s1->x.size();i++) {
      double newSigma;
      while(1) {
         // Compute the new sigma
         newSigma = s1->sigma(i) * fastExp(getGaussianRnd(0,deltaSigma(i)),false);

         // Check that we are within bounds
         if(!sBounded || (newSigma >= sLowerBound(i) && newSigma <= sUpperBound(i))) 
            break;
      }

		s1->sigma(i) = newSigma;
	}

   // Mutate x
	for(i=0;i<ss.dimensions();i++) {
      double newX;
		while(1) {
         // Compute the new x
         newX = getGaussianRnd(s1->x(i),s1->sigma(i));

         // Check that we are within bounds
         if(!xBounded || (newX >= xLowerBound(i) && newX <= xUpperBound(i))) 
            break;
      }

		s1->x(i) = newX;
	}
}

// Perform the crossover step
/*
void EvolutionaryStrategy::crossover(ESSolution *s1, ESSolution *s2, ESSolution *child)
{
	// Which crossover do we want to use.  Options are discrete and intermediate,
   if(true) {   // randBit()
		// Discrete crossover
		for(int i=0;i<ss.dimensions();i++) {
			child->x(i) = randBit() ? s1->x(i) : s2->x(i);
			child->sigma(i) = randBit() ? s1->sigma(i) : s2->sigma(i);
		}
	}
	else {
		// Intermediate crossover
		for(int i=0;i<ss.dimensions();i++) {
			child->x(i) = (s1->x(i) + s2->x(i))/2.0;
			child->sigma(i) = (s1->sigma(i) + s2->sigma(i))/2.0;
		}
	}
}
*/
int compareESSolutions(const void *s1,const void *s2) {
	// UGHhhh
	ESSolution *p1 = (ESSolution *)*((ESSolution **)(s1));
	ESSolution *p2 = (ESSolution *)*((ESSolution **)(s2));

	return (p1->getValue()==p2->getValue()) ? 0 :
	(p1->getValue() < p2->getValue()) ? -1 : 1;
}

void EvolutionaryStrategy::selection()
{
  int i, j;
	if(!initialized) {
      for(i=0;i<mu;i++) {
         if(!mlArray[i]->hasValue)		   // Set the first mu solutions to random values
            makeRandomSolution(mlArray[i]);
      }
		initialized = true;
	}

   // First generate the next lambda solutions.
   for(i=mu;i<mu+lambda;i++) {
      if(nCrossoverParents==2) {
         // Select two solutions for crossover
         int s[2];
         s[0] = rand(mu-1),s[1] = rand(mu-2);
         if(s[1] >= s[0]) s[1]++;

         // Crossover the two solutions and store into the child solution
         crossover(problem,nCrossoverParents,s,mlArray,mlArray[i]);

         // Mutate the solution
         mutate(mlArray[i]);

         // Set the value of the child solution
         mlArray[i]->setValue(problem.evaluate(mlArray[i]->x));

         // If the child is no different from its parents remove it from the solution space.
         // This is needed so that the solution space does not get overpopulated with the same
         // 'best' solution.  Replace with a random solution
         if(mlArray[i]->value == mlArray[0]->value ||
            mlArray[i]->value == mlArray[s[0]]->value ||
            mlArray[i]->value == mlArray[s[1]]->value)
         {
            makeRandomSolution(mlArray[i]);
         }
      }
      else {
         // Select n solutions for crossover
         bool *taken = new bool[mu];
         int *index = new int[nCrossoverParents];
         memset(taken,0,sizeof(bool)*mu);

         for(j=0;j<nCrossoverParents;j++) {
            index[j] = rand(mu);
            while(taken[index[j]]) index[j] = (index[j]+1)%mu;
            taken[index[j]] = true;            
         }

         // Crossover the two solutions and store into the child solution
         crossover(problem,nCrossoverParents,index,mlArray,mlArray[i]);

         delete taken;
         delete index;

         // Mutate the solution
         mutate(mlArray[i]);

         // Set the value of the child solution
         mlArray[i]->setValue(problem.evaluate(mlArray[i]->x));         
      }
   }

   // Now sort the solutions.
	if(mode == SELECTION_MuPlusLambda) {
      // In this mode we sort both the parent solutions and the child solutions
		qsort(mlArray,mu+lambda,sizeof(ESSolution *),compareESSolutions);
	} 
	else {
		// In this mode we discard the parent solutions and only sort the child solutions
		qsort(mlArray+mu,lambda,sizeof(Solution *),compareESSolutions);

      // Swap the first mu solutions with the second mu solutions
      ESSolution *tmp;
      for(i=0;i<mu;i++) {
         tmp = mlArray[mu+i];
         mlArray[mu+i] = mlArray[mu];
         mlArray[mu] = tmp;
      }
	}

   // Now we will compute the mean and the standard deviation of the surviving population   
   /*
	if(pRandSolution > 0.0) {
		for(int i=0;i<mu;i++) {
			if(rand(1.0) < pRandSolution) {
				makeRandomSolution(mlArray[i]);
			}
		}
	}
   */

   // Update the best ever solution if nessesary
   setBestEver(mlArray[0]);
}

void EvolutionaryStrategy::setDeltaSigma(const Vector &v)
{
	deltaSigma = v;
}

void EvolutionaryStrategy::setNth(int index,const Vector &x,const Vector &sigma)
{
	assert(index >=0 && index < mu);
	mlArray[index]->x = x;
   mlArray[index]->sigma = sigma;
   mlArray[index]->setValue(problem.evaluate(x));

	setBestEver(mlArray[index]);
}

void EvolutionaryStrategy::setNth(int index,const Vector &x) {
   assert(index >=0 && index < mu);
	mlArray[index]->x = x;
   computeRandomSigma(mlArray[index]->sigma);
   mlArray[index]->setValue(problem.evaluate(x));

	setBestEver(mlArray[index]);
}

ESSolution *EvolutionaryStrategy::getBest()
{
	return mlArray[0];
}

ESSolution *EvolutionaryStrategy::getNth(int n)
{
	return mlArray[n];
}

ESSolution * EvolutionaryStrategy::getBestEver()
{
   return &bestEver;
}

void EvolutionaryStrategy::setBestEver(ESSolution *s)
{
   if((!bestEver.hasValue) || (bestEver.value > s->value)) {
		bestEver = *s;
	}
}






// End namespace
NAMESPACE_PAULY_END

