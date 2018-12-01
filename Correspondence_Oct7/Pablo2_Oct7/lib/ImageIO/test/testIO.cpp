#include <iostream>
#include "Image3D.h"
#include "utility.h"
#include "AllImageIO.h"
#include "ControlParms.h"
#include "ControlKeys.h"
#include "control_parms_defaults.h"

using namespace std;

class ControlParms * globalControl;
int globalVerbosity;

void report(Image3D * image)
{
	cout << "Dimensions: " << image->getXDim() << ", " << image->getYDim()
		<< ", " << image->getZDim() << endl;
	cout << "Spacing: " << image->getXSpacing() << ", " << image->getYSpacing()
		<< ", " << image->getZSpacing() << endl;
	cout << "Origin: " << image->getXWorldOrigin() << ", " << image->getYWorldOrigin()
		<< ", " << image->getZWorldOrigin() << endl;

	GreyValue min, max;
	calc_intensity_range(*image, min, max);
	cout << "Intensity range = [" << min << ", " << max << "]\n";

	int z = image->getZDim();
	z /= 2;
	cout << "slice " << z << " trace:\n";
	int n = 0;
	int l = 0;
	int i;
	GreyValue * voxels = image->getVoxels();
	for (int j = 0; j < image->getYDim(); j++) {
		i = j;
		if (l >= 12)	// Print only 12 lines
			break;
		cout << voxels[i + image->getXDim()*(j + z*image->getYDim())];
		n++;
		if (n % 20 == 0) {
			cout << '\n';
			l++;
		}
		else
			cout << ' ';
	}
	cout << endl;
}

int mainI()
{
//	char in_filename[] = "../../../../data/pat/plan_im";
	char in_filename[] = "../../../../../data/3106.fr01.pim";
//	char in_filename[] = "../../../../data/images/junk.gipl";
//	char in_filename[] = "../../../../data/637.4.slice4-28.gmap.blur.gipl";
//	char in_filename[] = "../../../../Meta/exampleData/angio.mhd";

	char out_filename[] = "./junk.pim";

	Image3D * image;
	AllImageIO im_file;

	globalControl = new ControlParms(NULL, true);
	// Set parameter defaults
	setGlobalControlDefaults();

	cout << "\nReading image:\n";
	image = im_file.read(in_filename);
	if (image == NULL) {
		cout << "Read of " << in_filename << " failed" << endl;
		return 1;
	}
	report(image);

	// Write a copy of the image
	cout << "\nWriting image:\n";
	if (im_file.write(out_filename, *image) == false) {
		cout << "Write failed\n";
		return 1;
	}

	cout << "\nReading image:\n";
	image = im_file.read(out_filename);
	if (image == NULL) {
		cout << "Read of " << out_filename << " failed" << endl;
		return 1;
	}
	report(image);

	return 0;
}


