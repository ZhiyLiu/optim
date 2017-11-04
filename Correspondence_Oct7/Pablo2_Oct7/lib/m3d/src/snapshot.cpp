
#include <iostream>
#include "snapshot.h"
#include "Image3D.h"
#include "RAWImageFile.h"
#include "utility.h"


#define MAX_INTENS	((float) MAX_GREY_VALUE)


using namespace std;


void snapshot(const char * filename, Image3D & image, GreyValue min,
	GreyValue max)
{
	snapshot(filename, image.getVoxels(), image.getXDim(), image.getYDim(),
		image.getZDim(), min, max);
}


void snapshot(const char * filename, float * pixels, int x, int y, int z,
	float min, float max, float xSpacing, float ySpacing, float zSpacing)
{
	int size;
	register int k;
	float lo, hi;
	float range;

	if (filename == NULL || pixels == NULL) {
		cout << "Could not produce snapshot image" << endl;
		return;
	}

	Image3D dump(x, y, z);
	size = x*y*z;
	GreyValue * im = dump.getVoxels();
	if (xSpacing > 0.0f)
		dump.setSpacingAndOrigin(xSpacing, ySpacing, zSpacing);

	if (min == FLT_MAX || max == FLT_MIN) {
		// Compute the intensity range
		calc_intensity_range(pixels, size, lo, hi);
	}
	else {
		lo = min;
		hi = max;
	}
	range = hi - lo;
	if (range == 0.0) range = 1.0;

	for (k = 0; k < size; k++)
		im[k] = (GreyValue) (0.5f + MAX_INTENS*(pixels[k] - lo)/range);
	GreyValue mn, mx;
	mn = 0;
	mx = (GreyValue) MAX_INTENS;

	RAWImageFile file;
	file.write(filename, dump, mn, mx);

	cout << "Snapshot of image with intensity range [" << lo << ", " << hi
		<< "] written to\n  " << filename << " with range [0, "
		<< mx << ']' << endl;
}


void snapshot(const char * filename, GreyValue * pixels, int x, int y, int z,
	GreyValue min, GreyValue max, float xSpacing, float ySpacing, float zSpacing)
{
	int size;
	register int k;
	GreyValue lo, hi;

	if (filename == NULL || pixels == NULL) {
		cout << "Could not produce snapshot image" << endl;
		return;
	}

	Image3D dump(x, y, z);
	size = x*y*z;
	GreyValue * im = dump.getVoxels();
	if (xSpacing > 0.0f)
		dump.setSpacingAndOrigin(xSpacing, ySpacing, zSpacing);

	if (min == MAX_GREY_VALUE || max == MIN_GREY_VALUE) {
#ifdef AE2_BUILD
		lo = pixels[0]&IMAGE_BITS;
		hi = pixels[0]&IMAGE_BITS;
#else
		lo = pixels[0];
		hi = pixels[0];
#endif

		int start = size & 0x01;
		im[0] = pixels[0];

		// Compute the intensity range
		for(k = start; k < size; k += 2)
		{
			register GreyValue val0, val1;

#ifdef AE2_BUILD
			val0 = pixels[k]&IMAGE_BITS;
			val1 = pixels[k + 1]&IMAGE_BITS;
#else
			val0 = pixels[k];
			val1 = pixels[k + 1];
#endif
			im[k] = val0;
			im[k + 1] = val1;

			if(val0 < val1) {
				if(val0 < lo)
					lo = val0;
				if(val1 > hi)
					hi = val1;
			}
			else {
				if(val1 < lo)
					lo = val1;
				if(val0 > hi)
					hi = val0;
			}
		}
	}
	else {
		lo = min;
		hi = max;
		for (k = 0; k < size; k++)
#ifdef AE2_BUILD
			im[k] = pixels[k]&IMAGE_BITS;
#else
			im[k] = pixels[k];
#endif
	}

	RAWImageFile file;
	file.write(filename, dump, lo, hi);

	cout << "Snapshot of image with intensity range [" << lo << ", " << hi
		<< "] written to\n  " << filename << endl;
}


void snapshot(const char * filename, unsigned char * pixels, int x, int y, int z,
	unsigned char min, unsigned char max, float xSpacing, float ySpacing, float zSpacing)
{
	int size;
	register int k;
	GreyValue lo, hi;

	if (filename == NULL || pixels == NULL) {
		cout << "Could not produce snapshot image" << endl;
		return;
	}

	Image3D dump(x, y, z);
	size = x*y*z;
	GreyValue * im = dump.getVoxels();
	if (xSpacing > 0.0f)
		dump.setSpacingAndOrigin(xSpacing, ySpacing, zSpacing);

	int cmin = CHAR_MIN;
	int cmax = CHAR_MAX;
	if (((int) min) == cmax || ((int) max) == cmin)
	{
		lo = (int) pixels[0];
		hi = (int) pixels[0];

		int start = size & 0x01;
		im[0] = (GreyValue) pixels[0];

		// Compute the intensity range
		for(k = start; k < size; k += 2)
		{
			register GreyValue val0 = (GreyValue) pixels[k];
			register GreyValue val1 = (GreyValue) pixels[k + 1];
			im[k] = val0;
			im[k + 1] = val1;

			if(val0 < val1) {
				if(val0 < lo)
					lo = val0;
				if(val1 > hi)
					hi = val1;
			}
			else {
				if(val1 < lo)
					lo = val1;
				if(val0 > hi)
					hi = val0;
			}
		}
	}
	else {
		lo = min;
		hi = max;
		for (k = 0; k < size; k++)
			im[k] = (GreyValue) pixels[k];
	}

	RAWImageFile file;
	file.write(filename, dump, lo, hi);

	cout << "Snapshot of image with intensity range [" << lo << ", " << hi
		<< "] written to\n  " << filename << " with range [0, "
		<< MAX_INTENS << ']' << endl;
}

