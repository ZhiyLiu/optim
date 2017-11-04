#ifndef RAW_IMAGE_FILE_H
#define RAW_IMAGE_FILE_H

#include <stdio.h>
#include "Image3D.h"
#include "Vector3D.h"


extern const int UNKNOWN_GREY_VALUE;
extern const int CT_MIN_INTENS;
extern const int CT_MAX_INTENS;


/*  Global function to determine the byte order of the machine running
    Pablo.  Returns 0 on Solaris, 1 on Intel PCs, and -1 in all other
	cases.
*/
int cpuByteOrder();


/*  Structure used for storing optional header information.  This was
	absent from the original version of the raw3 format.  This could be
	a class, but until it contains more variables, there is little point.
*/
struct ExtendedHeader
{
	bool haveModality;	// true if modality != UNKNOWN_MODALITY
	modality_t modality;
	// Fields used only for DQF images
	int sWidth;			// Spectral width (items per voxel)
	// ROI origin (first voxel's image coordinates)
	int xOrigin, yOrigin, zOrigin;
	// Parameters for converting file intensities to DQF values
	double scale;
	int shift;
	// Dimensions of the image to which the ROI applies
	int separateXDim, separateYDim, separateZDim;
};


/*  Class RAWImageFile.

	See file pablo_format.txt for a description of the .raw3 image file
	format.

	This class simply reads and writes .raw3 files to and from Image3D
	objects.  There are two parameters affecting input operations, which
	are set by member functions setImageScaling() and setImageShifting().

	If setImageScaling() is called with true, then the intensities will be
	scaled from their input range to the maximum allowed by the GreyValue
	pixel type defined in Image3D.  If the value is false, no scaling is
	applied.  When the argument is true, the minimum intensity will be
	subtracted from all pixels, making them positive valued.  When it is
	false, any shift is determined by the setImageShifting() function.

	If setImageShifting() is called with true, then all intensities will be
	shifted by adding a constant.  The amount may be provided, but the
	standard amount, that is automatically used when no shift amount is
	minimum GreyValue of zero in the Image3D object returned by this class.
	provided, is the image's minimum intensity value.  This will result in a
	Intensity shifting is ignored when scaling is being done.

	Instead of using functions setImageScaling() and setImageShifting(),
	function setImageBitLength() may be used to specify the number of
	bits to which the image should be scaled and shifted.  If the range
	of voxels is already within nbits, then no remapping will be done
	when the image is read.  If a remapping occurs, the image will be
	stored to cover the full range of [0, (nbits << 1) - 1].

	Scaling, shifting, and setting the number of bits are only controlled
	through these functions and are not affected by construction or
	destruction.

    Images with particular modalities (presently CT and Shifted_CT images) are
    assigned the predefined intensity range of the modality before any
    remapping occurs, regardless of what intensity range is found in the
    image header.  This assignment will be omitted if setImageMapActual() is
    called with an argument of true before the I/O.

    Functions setImageShifting(), setImageBitLength(), and setImageMapActual()
    are only provided for special applications and should never be called from
    inside Pablo, except by the calls to them in class AllImageIO.
*/

class RAWImageFile
{
public:
    RAWImageFile() {}
    ~RAWImageFile() {}

    Image3D * read(const char * filename, bool headerOnly = false);
    bool write(const char * filename, const Image3D & image,
		GreyValue min = MIN_GREY_VALUE, GreyValue max = MAX_GREY_VALUE);
	static void setImageScaling(bool scaleInput);
	static void setImageShifting(bool shiftInput, int shiftAmount = 0);
	static void setImageBitLength(int nbits);
	static void setImageMapActual(bool yesNo);

private:

    bool skipComments(FILE * fp);
	void copyComments(const char * old_filename, FILE * fp);
	bool write(const char * filename, int xDim, int yDim, int zDim,
		double xSpacing, double ySpacing, double zSpacing,
		double xOrig, double yOrig, double zOrig, int minVal,
		int maxVal, const short * voxels, int byte_order, bool useGzip,
		ExtendedHeader * eh, const char * old_filename = NULL);

	static bool scale_input;
	static bool shift_input;
	static int shift;
	static int bits;
        static bool map_actual;
};

#endif

