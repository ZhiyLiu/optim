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
 * support.cpp
 *	---------
 * Support functions for Optima library 
 ******************************************************************/
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "libs.h"
#include <assert.h>
#include "support.h"

// Begin namespace
NAMESPACE_PAULY_START

#define PI 3.141592

// N(0,1) generation:
// Wow - Checked this against Excels NORMINV function - 
// we are very close with this crude function!

// Lookup table (goes to 0.5)
// Get Uniform normalized gaussian value
double gaussian(double x) {
	static double gaussianCoeff = 1.0/sqrt(2*PI);
	return gaussianCoeff*exp(-x*x/2);
}

// Compute a the exponent via lookup tables, in range -5 to 5.  The second parameter
// computeOutOfRange determines whether the values out of range are computed by calling
// exp or whether a constant is returned.  Default is to compute.
double fastExp(double x,bool computeOutOfRange) {
   static bool firstTime = true;
   static const double xMin = -5.0,xMax = 5.0,xStep = 0.001;
   static int tableSize = (int)((xMax-xMin) / xStep) + 1;
   static double *table;
   static double xToIndexMult = (tableSize-1) / (xMax-xMin);

   if(firstTime) {
      // Compute all the values in the table
      table = new double[tableSize];
      int tableIndex = 0;
      for(double xx = xMin;xx < xMax + xStep/2; xx+= xStep) {
         table[tableIndex] = exp(xx);
         tableIndex++;
      }
      firstTime = false;
   }

   if(x >= xMin && x<= xMax) {
      int index = (int)((x - xMin) * xToIndexMult);
      return table[index];
   }

   else if(computeOutOfRange) {
      return exp(x);
   }

   else if(x < xMin) {
      return 0;
   }

   else {
      return 1000000;
   }
}

// Sample the Gaussian distribution via lookup tables.
double getGaussianRnd(double mean,double sigma) {
   static bool firstTime = true;
   static const double xMax = 0.5,xStep = 0.0001;
   static int tableSize = (int)((xMax) / xStep) + 1;
   static double *table;
   static double xToIndexMult = (tableSize-1) / (xMax);

   if(firstTime) {
      // Compute all the values in the table
      table = new double[tableSize];
      double integral = 0.0;
      int tableIndex = 0;
      while(tableIndex < tableSize) {
         table[tableIndex] = integral;
         integral += xStep / gaussian(integral);
         tableIndex++;
      }

      firstTime = false;
   }

	double rnd = rand(2.0) - 1.0;
	double offset = (tableSize-1) * fabs(rnd);
	double result = table[(int)offset];
	if(rnd < 0)
		result = -result;
	return sigma*result+mean;
}


// Compute the 'smoothed box' penalty function.  This function
// is similar to the Gaussian, but uses a higher power
// for x in order to penalize values within one std. dev. less
// and ones outside of it more.
double penaltyFunction(double x,double mu,double sigma) {
   static bool firstTime = true;
   static const double xMax = 5.0,xStep = 0.001;
   static int tableSize = (int)((xMax) / xStep) + 1;
   static double *table;
   static double xToIndexMult = (tableSize-1) / (xMax);

   if(firstTime) {
      // Compute all the values in the table
      table = new double[tableSize];
      int tableIndex = 0;
      for(double xx = 0;xx < xMax + xStep/2; xx+= xStep) {
         table[tableIndex] = exp(- xx * xx * xx * xx / 2.0);
         tableIndex++;
      }
      firstTime = false;
   }

   double z = fabs((x-mu) / sigma);

   if(z<= xMax) {
      int index = (int)(z * xToIndexMult);
      return table[index];
   }

   return 0;
}

double penaltyFunction(double x,double mu,double sigma1,double sigma2) {
   double z = fabs((x-mu));
   if(z <= sigma1)
      return 1;
   else
      return penaltyFunction(z-sigma1,0,sigma2-sigma1);
}

// Get a random bit.  The current version assumes that rand_max is 0x7fff which seems to be pretty standard
bool randBit() {
   static short randomShort;
   static short bitIndex = 15;

   if(bitIndex == 15) {
      // Just in case someone compiles this on a screwed up operating system
      assert(2 == sizeof(short));

      // Get a random value
      randomShort = ::rand();

      // Bit index is 0
      bitIndex = 0;
   }

   bool rtn = ((0x0001 & randomShort) == 0x0001) ;
   randomShort >>= 1;
   bitIndex++;

   return rtn;
}



// End namespace
NAMESPACE_PAULY_END


