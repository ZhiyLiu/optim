
#include <stdio.h>
#include <string.h>
#include <iostream>
#include "BasicException.h"
#include "ImageIO.h"
#include "ImageStruct.h"
#include "Analyze.h"
#include "analyze_io.h"
#include "mayo_analyze.h"


using namespace std;


bool Analyze::load(string filename, ImageStruct & image3D)
{
	struct dsr * hdr;
	char img_filename[MAXPATHLEN];
	char hdr_filename[MAXPATHLEN];
	char fprefix[MAXPATHLEN];
	bool swapBytes;

	constructAnalyzeNames(filename.c_str(), fprefix, hdr_filename, img_filename);

	hdr = read_analyze_header(hdr_filename, ANDT_SIGNED_SHORT, swapBytes);
	if (hdr == NULL)
		return false;

	image3D.dims[0] = hdr->dime.dim[1];
	image3D.dims[1] = hdr->dime.dim[2];
	image3D.dims[2] = hdr->dime.dim[3];

	image3D.spacing[0] = hdr->dime.pixdim[1];
	image3D.spacing[1] = hdr->dime.pixdim[2];
	image3D.spacing[2] = hdr->dime.pixdim[3];

	// Apparently, Analyze files do not contain the world-coordinate origin
	image3D.origin[0] = 0.0;
	image3D.origin[1] = 0.0;
	image3D.origin[2] = 0.0;

	image3D.min = hdr->dime.glmin;
	image3D.max = hdr->dime.glmax;

	if (! read_analyze_image(img_filename, hdr->dime.datatype, &image3D, swapBytes)) {
		printf("Error reading %s\n", img_filename);
		return false;
	}
	delete hdr;

	return true;
}

bool Analyze::save(string filename, ImageStruct & image3D)
{
	char img_filename[MAXPATHLEN];
	char hdr_filename[MAXPATHLEN];
	char dbase[MAXPATHLEN];

	constructAnalyzeNames(filename.c_str(), dbase, hdr_filename, img_filename);
	if (! write_analyze(hdr_filename, img_filename, dbase, &image3D)) {
		printf("Error writing %s\n", filename.c_str());
		return false;
	}

	return true;
}


// Given a filename return :   filename.img (for the data)
//                             filename.hdr (for the header)
void Analyze::constructAnalyzeNames(const char * filename, char * database,
	char * hdr_path, char * img_path)
{
	char * dot = strrchr((char *) filename, '.');
	if (0 == strcmp(dot + 1, "hdr")) {
		strcpy(hdr_path, filename);
		*dot = '\0';
		strcpy(img_path, filename);
		strcat(img_path, ".img");
	}
	else if (0 == strcmp(dot + 1, "img")) {
		strcpy(img_path, filename);
		*dot = '\0';
		strcpy(hdr_path, filename);
		strcat(hdr_path, ".hdr");
	}
	char * s1 = strrchr((char *) filename, '/');
#ifdef _MSC_VER
	char * s2 = strrchr((char *)filename, '\\');
	if (s2 > s1)
		s1 = s2;
#endif
	if (s1 != NULL)
		strcpy(database, s1 + 1);
	else
		strcpy(database, filename);
	*dot = '.';
}


