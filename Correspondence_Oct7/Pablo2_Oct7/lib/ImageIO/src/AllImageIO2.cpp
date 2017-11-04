
// This file must include the files needed for performing *.raw3 I/O
#include "RAWImageFile.h"
#include "ImageStruct.h"
#include "Image3D.h"
#include "utility.h"
#include "AllImageIO.h"


//#define DEBUG


using namespace std;


bool AllImageIO::scale_input = true;
bool AllImageIO::shift_input = true;
int AllImageIO::shift = 0;
int AllImageIO::bits = 0;
bool AllImageIO::map_actual = false;
AllImageIO::image_t AllImageIO::last_format = NoImageFormat;


void AllImageIO::setImageScaling(bool scaleInput)
{
	scale_input = scaleInput;
    bits = 0;
	RAWImageFile::setImageScaling(scaleInput);
}

void AllImageIO::setImageShifting(bool shiftInput, int shiftAmount)
{
	shift_input = shiftInput;
	shift = shiftAmount;
    bits = 0;
	RAWImageFile::setImageShifting(shiftInput, shiftAmount);
}

void AllImageIO::setImageBitLength(int nbits)
{
    bits = nbits;
    shift_input = false;
    shift = 0;
    scale_input = false;
	RAWImageFile::setImageBitLength(nbits);
}

void AllImageIO::setImageMapActual(bool yesNo)
{
    map_actual = yesNo;
	RAWImageFile::setImageMapActual(yesNo);
}

Image3D * AllImageIO::read_raw3(const char * filename, bool headerOnly)
{

    RAWImageFile im_file;

	return im_file.read(filename, headerOnly);
}

// Write an image to the specified file.
bool AllImageIO::write_raw3(const char * filename, const Image3D * image,
	unsigned short min, unsigned short max)
{
    RAWImageFile im_file;

	return im_file.write(filename, *image, min, max);
}

Image3D * AllImageIO::convertInputImage(ImageStruct & image_struct, bool stacked,
	bool verbose)
{
	int min, max, diff;
	unsigned long size;
	double scale;
    bool scaling;
    bool shifting;

	size = image_struct.len;
	if (size == 0)
		size = image_struct.dims[0]*image_struct.dims[1]*image_struct.dims[2];

	if (image_struct.modality == CT && ! map_actual) {
		min = CT_MIN_INTENS;
		max = CT_MAX_INTENS;
#ifdef DEBUG
		cout << "CT image detected; using intensity range of [" << min << ", "
			<< max << "]\n";
#endif
	}
	else if (image_struct.modality == SHIFTED_CT && ! map_actual) {
		min = 0;
		max = CT_MAX_INTENS - CT_MIN_INTENS;
#ifdef DEBUG
		cout << "Shifted CT image detected; using intensity range of [" << min << ", "
			<< max << "]\n";
#endif
	}
	else if (image_struct.min < image_struct.max) {
		min = image_struct.min;
		max = image_struct.max;
#ifdef DEBUG
		cout << "Using file's intensity range of [" << min << ", " << max << "]\n";
#endif
	}
    else {
		if (image_struct.voxels == NULL) {
			min = 0;
			max = 0;
		}
        else if (image_struct.dataIsShort)
            calc_rangeWithSwap(true, min, true, max, false,
                (short *) image_struct.voxels, size);
        else {
            GreyValue n, x;
	    	calc_intensity_range(image_struct.voxels, (int) size, n, x);
            min = n;
            max = x;
        }
    }
	if (stacked)
		max = MAX_GREY_VALUE;
	// At this point, the data's intensity range will be [min, max], unless the image
	// was CT or SHIFTED_CT without map_actual having been set to true (when a fixed
	// range is used), or unless the image's header contained  incorrect values for
	// the range (not likely).  The max is probably wrong in the case of stacked
	// images, and I need to look into this.				    AGG: NOTE!

	if (verbose) {
		if (image_struct.modality != UNKNOWN_MODALITY)
			cout << "Modality is " << modalityString((modality_t) image_struct.modality) << '\n';
		cout << "Voxel intensity range = [" << min << ", " << max << "]\n";
	}

	Vector3D worldOriginPos(image_struct.origin[0], image_struct.origin[1],
		image_struct.origin[2]);
#ifdef DEBUG
	cout << "Original image origin: ";
	worldOriginPos.print();
#endif

	// Construct the image object
    Image3D * image = new Image3D(image_struct.dims[0], image_struct.dims[1],
		image_struct.dims[2]);
#ifdef DEBUG
	cout << "Dimensions: " << image_struct.dims[0] << " x " << image_struct.dims[1]
		<< " x " << image_struct.dims[2] << '\n';
#endif

    if (bits > 0) {
        int maxIntens = (1 << bits) - 1;
        // Compare current and future ranges
        if (max - min > maxIntens) {
            scaling = true;
            diff = min;
            shifting = (diff != 0);
            scale = maxIntens/(double) (max - min);
        }
        else {
            // Only need possibly to shift
            scaling = false;
            scale = 1.0;
            if (max > maxIntens) {
                diff = max - maxIntens;
                shifting = true;
            }
            else if (min < 0) {
                diff = min;
                shifting = true;
                if (diff == 0)
                    shifting = false;
            }
			else {
				shifting = false;
				diff = 0;
			}
        }
        if (verbose && image_struct.voxels != NULL)
            cout << "Voxels will be mapped to fall within [0, " << max << "], " << bits << " bits\n";
    }
    else {
	    if (image_struct.modality != CT && scale_input)
		    scale = (double) MAX_GREY_VALUE / (max - min);
	    else
		    scale = 1.0;

	    if (! scale_input && shift_input && shift != 0)
		    diff = shift;
	    else
		    diff = min;

	    if (verbose && image_struct.voxels != NULL) {
			if (image_struct.modality == SHIFTED_CT)
			    cout << "Voxels will be neither shifted nor scaled\n";
			else {
				if (scale_input)
					cout << "Voxels will be mapped to [0, " << MAX_GREY_VALUE << "]\n";
				else if (shift_input)
					cout << "Voxels will be mapped to [0, " << max - diff << "]\n";
				else
					cout << "Voxels will be neither shifted nor scaled\n";
			}
	    }

	    scaling = scale_input;
	    shifting = shift_input;
	    if (scale == 1.0)
		    scaling = false;
	    if (diff == 0)
		    shifting = false;
    }
#ifdef DEBUG
	cout << "stacked: " << stacked << ",  scale: " << scale << ",  scale_input: "
		<< scale_input << ",  shift: " << shift << ",  shift_input: "
		<< shift_input << endl;
#endif

	if (image_struct.voxels == NULL) {
		// Loading the header only
		image->clear(0, false);
	}
	else {
		double xlatepixval;
		GreyValue finalpixval;

		GreyValue * voxels = image->getVoxels();
        //xiaojie 
		//TOMODIFY: for downloading antialised distance image for calculating, we don't want to scale or shift it
		bool adjust = scaling || shifting ; // false;
		for (int index = 0; index < size; index++) {
			if (adjust) {
		//            ?????????? AGG:  how can this work when voxels contains shorts?
				//Xiaojie Comment
				//xlatepixval = image_struct.voxels[index] - diff;
				//finalpixval = (GreyValue) (scale * xlatepixval + 0.5);
				//Xiaojie
				xlatepixval = (short)image_struct.voxels[index] - diff;
                finalpixval = (GreyValue) (scale * xlatepixval + 0.5);
			}
			else
				finalpixval = image_struct.voxels[index];

			voxels[index] = finalpixval;
		}

		if (adjust) {
		    min = (int) ((min - diff)*scale + 0.5);
		    max = (int) ((max - diff)*scale + 0.5);
		}
		// Now the data's current range is [min, max]
	}

	image->setSpacingAndOrigin(image_struct.spacing[0], image_struct.spacing[1],
		image_struct.spacing[2], &worldOriginPos);
    image->setRange(min, max);      // Range of voxels in image
	image->setIntensityMapping(scale, diff);
#ifdef DEBUG
	cout << "setRange()called with min: " << min << ",  max: " << max << endl;
	cout << "setIntensityMapping() called with scale: " << scale << ",  diff: "
		<< diff << endl;
#endif
	image->setModality((modality_t) image_struct.modality);

	return image;
}

void AllImageIO::remapImageForOutput(ImageStruct & image_struct,
	unsigned short min, unsigned short max, unsigned short * outVoxels,
	bool verbose)
{
	unsigned short diff;
	double scale;
	int i;

#ifdef DEBUG
	cout << "Dimensions: " << image_struct.dims[0] << " x " << image_struct.dims[1]
		<< " x " << image_struct.dims[2] << '\n';
#endif

	// If the range is more than 15 bits, it must also be reduced,
	// since most formats store voxels as short integers.
	if (max > MAX_GREY_VALUE_WRITTEN) {
		if (max - min <= MAX_GREY_VALUE_WRITTEN) {
			diff = max - MAX_GREY_VALUE_WRITTEN;
			for (i = 0; i < image_struct.len; i++)
				outVoxels[i] = image_struct.voxels[i] - diff;
			min -= diff;
			max -= diff;
		}
		else {
			scale = MAX_GREY_VALUE_WRITTEN/((double) (max - min));
			for (i = 0; i < image_struct.len; i++)
				outVoxels[i] = (unsigned short)
					(0.5 + scale*(image_struct.voxels[i] - min));
			min = 0;
			max = (unsigned short) (0.5 + scale*(max - min));
		}
	}
	else {
			for (i = 0; i < image_struct.len; i++)
			outVoxels[i] = image_struct.voxels[i];
	}
	image_struct.min = min;
	image_struct.max = max;

	if (verbose)
		cout << "Voxel intensity range = [" << min << ", " << max << "]\n";
}

void AllImageIO::extractOutputImageInfo(ImageStruct & image_struct,
	unsigned short * min, unsigned short * max, const Image3D & image3d)
{
    image_struct.dims[0] = image3d.getXDim();
	image_struct.dims[1] = image3d.getYDim();
	image_struct.dims[2] = image3d.getZDim();
	image_struct.len = image_struct.dims[0]*image_struct.dims[1]*image_struct.dims[2];

	Vector3D origin = image3d.getWorldOrigin();
	image_struct.origin[0] = origin.getX();
	image_struct.origin[1] = origin.getY();
	image_struct.origin[2] = origin.getZ();

	image_struct.spacing[0] = image3d.getXSpacing();
	image_struct.spacing[1] = image3d.getYSpacing();
	image_struct.spacing[2] = image3d.getZSpacing();

	image_struct.voxels = image3d.getVoxels();
    image_struct.dataIsShort = false;

    image_struct.modality = image3d.modality();

	if (*min > *max)
		calc_intensity_range((GreyValue *)image_struct.voxels, image_struct.len,
			*((GreyValue *) min), *((GreyValue *) max));
#ifdef DEBUG
	cout << "Voxel intensity range = [" << *min << ", " << *max << "]\n";
#endif
	image_struct.min = *min;
	image_struct.max = *max;
}


