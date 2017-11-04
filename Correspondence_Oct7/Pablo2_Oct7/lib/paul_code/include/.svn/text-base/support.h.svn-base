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
 * support.h
 *	---------
 * Support functions for Optima library 
 ******************************************************************/
#ifndef _OPTIMA_SUPPORT_
#define _OPTIMA_SUPPORT_

#include "libs.h"
#include <time.h>

// Begin namespace
NAMESPACE_PAULY_START

// Compute a the exponent via lookup tables, in range -5 to 5.  The second parameter
// computeOutOfRange determines whether the values out of range are computed by calling
// exp or whether a constant is returned.  Default is to compute.
double fastExp(double x,bool computeOutOfRange=true);

// Sample the gaussian distribution
double getGaussianRnd(double mean,double sigma);

// Compute the 'smoothed box' penalty function.  This function
// is similar to the Gaussian, but uses a higher power
// for x in order to penalize values within one std. dev. less
// and ones outside of it more.
double penaltyFunction(double x,double mu,double sigma);

/* Same as above, but the region between mu and sigma1 has value 1 
 and smoothes down to zero after sigma1 with s.d. sigma2-sigma1


       -------
     /
   /
--/
    s2  s1  mu  
*/
double penaltyFunction(double x,double mu,double sigma1,double sigma2);

// Return a random bit
bool randBit();

// Get a uniformly distrubuted number
inline double rand(double max) {
   return max * ((double)::rand()) / RAND_MAX;
}

inline int rand(int max) {
   return (max * ::rand()) / RAND_MAX;
}


// An inline timer class for profiling programs
class iTimer {
private:
   clock_t tStart;

   long tTotal;
   int laps;

   bool running;
public:
   void start() {
      if(!running) {
         running = true;
         tStart = clock();
      }
   }

   void stop() {
      if(running) {
         tTotal += clock() - tStart;
         running = false;
         laps++;
      }
   }

   void nextLap() {
      stop();
      start();
   }

   void reset() {
      stop();
      tTotal = 0;
      laps = 0;
   }

   double ms() {
      return 1000.0 * tTotal / CLOCKS_PER_SEC;
   }

   int nLaps() {
      return laps;
   }

   iTimer() {
      running = false;
      reset();
   }
};



// End namespace
NAMESPACE_PAULY_END

#endif	/* _OPTIMA_SUPPORT_ */

