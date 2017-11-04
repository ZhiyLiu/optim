/********************************************************************************/
/*																				*/
/*  	File	:  Samplestats.cpp												*/
/*																				*/
/*	Description:  statistical analysis class---member functions for class		*/
/*																				*/
/*	Project :  Seurat															*/
/*																				*/
/*	Author  :  A. Thall															*/
/*																				*/
/*	Date	:  16. April 2002													*/
/********************************************************************************/

#include <iostream>
#include <math.h>
#define D_SAMPLESTATS
#include "Shapedepend.h"

using namespace ThallCode;
/*
 * getmean() -- compute mean of initialized samples...will perform on whatever subarray
 *    samples[0..samploc-1] has been assigned values
 */
double Samplestats::getmean()
{
	mean = 0;

	for (int samp = 0; samp < samploc; samp++)
		mean += samples[samp];
	mean /= samploc;
	meancomputed = true;
	varcomputed = false;
	stddevcomputed = false;

	return mean;
}

/*
 * getvariance() -- return samplevariance based on computed mean
 */
double Samplestats::getvariance()
{
	double diff;

	if (!meancomputed)
		getmean();

	samplevariance = 0;
	for (int samp = 0; samp < samploc; samp++) {
		diff = samples[samp] - mean;
		samplevariance += diff*diff;
	}
	samplevariance /= (samploc - 1);
	varcomputed = true;
	stddevcomputed = false;

	return samplevariance;
}

/*
 * getstddev() -- return samplestddev based on computed variance
 */
double Samplestats::getstddev()
{
	if (!varcomputed)
		getvariance();

	samplestddev = sqrt(samplevariance);
	stddevcomputed = true;

	return samplestddev;
}


