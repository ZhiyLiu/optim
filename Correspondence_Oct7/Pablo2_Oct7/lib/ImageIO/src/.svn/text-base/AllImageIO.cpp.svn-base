
// This file should include no files from the application using AllImageIO
#include <iostream>
#include "BasicException.h"
#include "AllImageIO.h"
#include "ImageStruct.h"
#include "ImageIO.h"


//#define DEBUG



using namespace std;


extern int globalVerbosity;			// Current verbosity level of Pablo


Image3D * AllImageIO::read(const char * filename, bool stacked, bool headerOnly)
{
	ImageIO mix;
	Image3D * image;

#ifdef DEBUG
	cout << "AllImageIO::read()" << endl;
#endif
	string fn = filename;
    last_format = NoImageFormat;

	ImageIO::ImageType extension = mix.guessImageFormat(fn);
    if (extension == ImageIO::raw3) {
	image = read_raw3(filename, headerOnly);                // AGG: How can this handle stacked images?
        if (image != NULL)
            last_format = Raw3;
        return image;
    }

	ImageStruct image_struct;
	image_struct.voxels = NULL;
	image_struct.len = 0;
	image_struct.modality = 0;	// UNKNOWN_MODALITY

	try {
		mix.loadThisImage(fn, image_struct, headerOnly);
	}
	catch (BasicException x) {
		cout << x.getMessage() << endl;
		return NULL;
	}

	if (image_struct.voxels == NULL && ! headerOnly)
		return NULL;

	image = convertInputImage(image_struct, stacked, globalVerbosity == 1);
    last_format = (image_t) extension;

	delete [] image_struct.voxels;
	image_struct.voxels = NULL;				// Should not be needed

	return image;
}

// Write an image to the specified file.
bool AllImageIO::write(const char * filename, const Image3D & image,
	int minIntensity, int maxIntensity, image_t ifmt)
{
	ImageIO mix;
	unsigned short * voxels;
	unsigned short min;
	unsigned short max;
	ImageIO::ImageType extension;
	ImageStruct image_struct;

#ifdef DEBUG
	cout << "AllImageIO::write()" << endl;
#endif
	string fn = filename;

	if (ifmt != NoImageFormat)
		extension = (ImageIO::ImageType) ifmt;
	else
		extension = mix.guessImageFormat(fn);
	if (extension == ImageIO::unknown || extension == ImageIO::raw3) {
		if (write_raw3(filename, &image, minIntensity, maxIntensity)) {
            last_format = Raw3;
			return true;
		}
		else {
            last_format = NoImageFormat;
			return false;
		}
	}

	min = minIntensity;   // AGG: minIntensity may be < 0, since the voxels are type u_short
	max = maxIntensity;                         // AGG: See P3DControl::saveImage()
	extractOutputImageInfo(image_struct, &min, &max, image);
	if (image_struct.voxels == NULL) {
		cout << "Image contains no voxels\n";
		return false;
	}

	voxels = new u_short[image_struct.len];
	remapImageForOutput(image_struct, min, max, voxels, globalVerbosity == 1);
	//	When writing, the original image_struct.voxels belongs to the caller
	image_struct.voxels = voxels;

	try {
		mix.saveImage(fn, image_struct);
	}
	catch (BasicException x) {
		delete [] image_struct.voxels;
        last_format = NoImageFormat;
		cout << x.getMessage() << endl;
		return false;
	}
    last_format = (image_t) extension;

	delete [] image_struct.voxels;
	image_struct.voxels = NULL;				// Should not be needed

	return true;
}


