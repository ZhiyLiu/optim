#include "Hermite.h"


namespace Hermite {

inline double sqr(const double x) { return x * x; }

// evaluates h(t) where h(t) is hermite cubic
// with values p0, p1, and derivatives d0, d1
double cubic(double t, double p0, double p1, double d0, double d1) {
  double h00 = 2*(t*t*t) - 3 * (t*t) + 1;
  double h10 = (t*t*t) - 2*(t*t) + t;
  double h01 = -2*(t*t*t) + 3 * (t*t);
  double h11 = (t*t) * (t - 1);

  return (h00 * p0) + (h01 * p1) + (h10 * d0) + (h11 * d1);

}

// finds t: 0<=t<=1 s.t. the h(t) == val, where h(t) is hermite cubic
// with values p0, p1, and derivatives d0, d1
double invertCubic(double val, double p0, double p1,  double d0, double d1) {
  if (val == p0) {
    return 0;
  }
  else if (val == p1) {
    return 1;
  } 
  double a = 0;
  double b = 1;

  if (p0 > p1) {
    a = 1;
	b = 0;
  }

  double bestT = 0;
  double bestDif = (val - p0) * (val - p0);
  double temp;
  double diff;
  for (double t=0; t<=1; t+=0.005) {
    temp = cubic(t, p0, p1, d0, d1);
	diff = (temp - val) * (temp - val);
	if (diff < bestDif) {
	  bestDif = diff;
	  bestT = t;
    }
  }
}
	

};
