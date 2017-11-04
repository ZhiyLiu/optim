#include <math.h>
#include "GaussianBlur3D.h"
#include "ControlParms.h"

using namespace std;

void GaussianBlur3D::blur(Image3D & image, double sigma, int kernelSize)
{
    double * kernel;
    GreyValue * data;
    GreyValue * newData;
    GreyValue * dataPtr;
    GreyValue * newDataPtr;
    int minKernelIndex, maxKernelIndex;
    int xDim, yDim, zDim;
    int imageSize;
    int i, j, xIndex, yIndex, zIndex;
    int stepSize;
    double val;
    double factor;

	if (sigma < 1.0e-6)
		return;

	if (globalVerbosity >= 1)
		cout << "Starting Gaussian blurring ... " << endl;
	else if (globalVerbosity == 0)
		cout << "Gaussian blurring image" << endl;

    // Make sure we have an odd-sized kernel
    if (kernelSize % 2 == 0)
        kernelSize++;

    maxKernelIndex = kernelSize/2;
    minKernelIndex = -maxKernelIndex;
    kernel = new double[kernelSize];

    data = image.getVoxels();
    xDim = image.getXDim();
    yDim = image.getYDim();
    zDim = image.getZDim();
    imageSize = xDim * yDim * zDim;

    newData = new GreyValue[imageSize];
    if (newData == NULL) {
        cerr << "Cannot allocate space for image" << endl;
        return;
    }

	double minspacing = image.getAbsXSpacing();
	if (minspacing > image.getAbsYSpacing())
		minspacing = image.getAbsYSpacing();
	if (minspacing > image.getAbsZSpacing())
		minspacing = image.getAbsZSpacing();

	double sum = 0.0;

	// AGG: the 3 large loops below probably can be combined for speed

    // Convolve in X; Change sigma
	double scaledsigma = sigma*(minspacing/image.getAbsXSpacing()); 
	factor = 1.0 / (sqrt(R_TWO_PI) * scaledsigma);
	for (i = minKernelIndex; i <= maxKernelIndex; i++) {
        kernel[i - minKernelIndex] =
			factor * exp((-0.5 / (scaledsigma * scaledsigma)) * (double)(i * i));
		sum += kernel[i - minKernelIndex];
		if (globalVerbosity >= 1)
			cout << kernel[i - minKernelIndex];
		if (i < maxKernelIndex) cout << ", ";
	}
	for (i = minKernelIndex; i <= maxKernelIndex; i++)
		kernel[i - minKernelIndex] /= sum;

    dataPtr = data;
    newDataPtr = newData;
    for (i = 0; i < imageSize; i++) {
        xIndex = i % xDim;

        val = 0;
        for (j = minKernelIndex; j <= maxKernelIndex; j++) {
            if (xIndex + j < 0)
                val += (double)(*(dataPtr - xIndex))*kernel[j - minKernelIndex];
            else if (xIndex + j >= xDim)
                val += (double)(*(dataPtr + xDim - xIndex - 1))*kernel[j - minKernelIndex];
            else
                val += (double)(*(dataPtr + j))*kernel[j - minKernelIndex];
        }

        (*newDataPtr) = (GreyValue) val;
        dataPtr++;
        newDataPtr++;
    }

    // Convolve in Y; Change sigma
	scaledsigma = sigma*(minspacing/image.getAbsYSpacing()); 
	factor = 1.0 / (sqrt(R_TWO_PI) * scaledsigma);
	if (globalVerbosity >= 1) cout << '\n';
	sum = 0.0;
	for (i = minKernelIndex; i <= maxKernelIndex; i++) {
        kernel[i - minKernelIndex] =
			factor * exp((-0.5/(scaledsigma * scaledsigma)) * (double)(i * i));
		sum += kernel[i - minKernelIndex];
		if (globalVerbosity >= 1)
			cout << kernel[i - minKernelIndex];
		if (i < maxKernelIndex) cout << ", ";
	}
	for (i = minKernelIndex; i <= maxKernelIndex; i++)
		kernel[i - minKernelIndex] /= sum;

    dataPtr = newData;
    newDataPtr = data;
    for (i = 0; i < imageSize; i++) {
        xIndex = i % xDim;
        yIndex = ((i - xIndex) / xDim) % yDim;

        val = 0;
        for (j = minKernelIndex; j <= maxKernelIndex; j++) {
            if (yIndex + j < 0)
                val += (double)(*(dataPtr - yIndex*xDim))*kernel[j - minKernelIndex];
            else if (yIndex + j >= yDim)
                val += (double)(*(dataPtr + (yDim - yIndex - 1)*xDim))*kernel[j - minKernelIndex];
            else
                val += (double)(*(dataPtr + j*xDim))*kernel[j - minKernelIndex];
        }

        (*newDataPtr) = (GreyValue) val;
        dataPtr++;
        newDataPtr++;
    }

    // Convolve in Z; Change sigma
	scaledsigma = sigma*(minspacing/image.getAbsZSpacing()); 
	factor = 1.0 / (sqrt(R_TWO_PI) * scaledsigma);
	if (globalVerbosity >= 1) cout << '\n';
	sum = 0.0;
	for (i = minKernelIndex; i <= maxKernelIndex; i++) {
        kernel[i - minKernelIndex] =
			factor * exp((-0.5 / (scaledsigma * scaledsigma)) * (double)(i * i));
		sum += kernel[i - minKernelIndex];
		if (globalVerbosity >= 1)
			cout << kernel[i - minKernelIndex];
		if (i < maxKernelIndex) cout << ", ";
	}
	for (i = minKernelIndex; i <= maxKernelIndex; i++)
		kernel[i - minKernelIndex] /= sum;

    dataPtr = data;
    newDataPtr = newData;
    stepSize = xDim * yDim;
    for (i = 0; i < imageSize; i++) {
        xIndex = i % xDim;
        yIndex = ((i - xIndex) / xDim) % yDim;
        zIndex = (i - xIndex - xDim * yIndex) / stepSize;

        val = 0;
        for (j = minKernelIndex; j <= maxKernelIndex; j++) {
            if (zIndex + j < 0)
                val += (double)(*(dataPtr - zIndex*stepSize))*kernel[j - minKernelIndex];
            else if (zIndex + j >= zDim)
                val += (double)(*(dataPtr + (zDim - zIndex - 1)*stepSize))*kernel[j - minKernelIndex];
            else
                val += (double)(*(dataPtr + j*stepSize))*kernel[j - minKernelIndex];
        }

        (*newDataPtr) = (GreyValue) val;
        dataPtr++;
        newDataPtr++;
    }

    delete [] kernel;

    (void) image.replaceVoxels(newData);
	if (globalVerbosity >= 1)
		cout << "\n  ... done." << endl;
}


