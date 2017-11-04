#ifndef UTILITY_H
#define UTILITY_H


#include <iostream>
#include <math.h>
#include <vector>
#include "M3DPrimitive.h"
#include "Image3D.h"


//#define DEBUG		/* Uncomment to get lots of printed output */


#ifndef M_PI
#define M_PI	3.1415926535897932384626433832795
#endif
const double RAD = M_PI/180.0;

#define Point	Vector3D


double inter_pt_distance(const Point & p1, const Point & p2);	// Distance between points

Point project_point_to_line(const Point & p, const Point & l0, const Point & l1, double & t);
inline Point project_point_to_line(const Point & p, const Point & l0, const Point & l1) {
	double t;
	return project_point_to_line(p, l0, l1, t);
}
double pt_to_line_distance(const Point & p, const Point & l0, const Point & l1);

Point project_point_to_plane(const Point & p, const Vector3D & n, double d);
Point project_point_to_plane(const Point & p, const Point & v0,
							 const Point & v1, const Point & v2);

Point intersect_line_with_plane(const Point & l0, const Point & l1,
								const Point & v0, const Point & v1, const Point & v2);

void print_image_info(Image3D & im, const char * name = NULL);

bool above_line(const Point & m0, const Point & m1, const Point & p);
bool inside_triangle(const Point & p, const Point & v0, const Point & v1,
					 const Point & v2);

#ifdef SWIG
// This function has not been defined
%ignore print_image_conversions();
#endif
void print_image_conversions();
void calc_intensity_range(Image3D & image, GreyValue & min, GreyValue & max);
void report_intensity_range(Image3D & image, char * label = NULL);	// Diagnostic output
void report_intensity_range(GreyValue * pixels, int size, char * label = NULL);	// Diagnostic output

void calc_intensity_range(GreyValue * pixels, int size, GreyValue & min, GreyValue & max);
void calc_intensity_range(float * pixels, int size, float & min, float & max);
void calc_intensity_range(unsigned char * pixels, int size, unsigned char & min,
    unsigned char & max);

void calc_rangeWithSwap(bool findMin, int & min, bool findMax, int & max,
	bool swapBytes, short * data, int size);


void verify_squareness(const Point & origin, const Point & corner,
					   const Point & end, char * label);

void beep();
#ifdef _WIN32
void reverseBackSlashes(char * str);
#endif

// Returns true if the file exists and is a regular file
bool fileExists(const char * filename);

// A buffered fread() that prevents some failures
int safeFread(void * buffer, unsigned long size, FILE * fp);


#endif	/* UTILITY_H */

