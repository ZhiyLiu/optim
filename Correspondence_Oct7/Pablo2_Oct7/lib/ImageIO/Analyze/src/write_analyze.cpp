/*

Some of the contents of this file are modifications of files from the
AFNI (Analyses of Functional Neuroimages) package copyrighted by the
	Medical College of Wisconsin
	8701 Watertown Plank Road
	Milwaukee, WI 53226
See http://afni.nimh.nih.gov/afni/.
*/

#include <stdio.h>
#include <string.h>
#include "ImageStruct.h"
#include "analyze_io.h"
#include "mayo_analyze.h"


/*-----------------------------------------------------------------
  Write an image in ANALYZE 7.5 format;
  fname.hdr and fname.img will be written, in native byte order.
                                                29 Nov 2001 - RWCox
-------------------------------------------------------------------*/

bool write_analyze(char * hdr_name, char * img_name, char * database,
	ImageStruct * im)
{
	struct dsr hdr ;  /* ANALYZE .hdr */
	void *ip ;
	FILE *fp ;
	int len;


	/*-- check inputs for sanity --*/

	if (hdr_name == NULL || hdr_name[0] == '\0' || img_name == NULL
	   || img_name[0] == '\0' || im == NULL)
		   return false;

	ip = im->voxels;
	len = im->dims[0]*im->dims[1]*im->dims[2];

	/*-- clear header --*/

	memset(&hdr, 0, sizeof(struct dsr));

	/*-- set random header stuff --*/

	hdr.hk.sizeof_hdr = sizeof(struct dsr);
	hdr.hk.regular = 'r';

	sprintf(hdr.hk.db_name, "%.17s", database);

	/*-- set data dimensions --*/

	hdr.dime.dim[0] = 4;
	hdr.dime.dim[1] = im->dims[0]; hdr.dime.pixdim[1] = im->spacing[0];
	hdr.dime.dim[2] = im->dims[1]; hdr.dime.pixdim[2] = im->spacing[1];
	hdr.dime.dim[3] = im->dims[2]; hdr.dime.pixdim[3] = im->spacing[2];
	hdr.dime.dim[4] = 1; hdr.dime.pixdim[4] = 0;
	hdr.dime.dim[5] = 0; hdr.dime.pixdim[5] = 0;
	hdr.dime.dim[6] = 0; hdr.dime.pixdim[6] = 0;
	hdr.dime.dim[7] = 0; hdr.dime.pixdim[7] = 0;

	hdr.dime.funused1 = 1.0;  // SPM scale factor
	strcpy(hdr.hk.data_type, "dsr");						// AGG added this

	/*-- set data type --*/

	hdr.dime.datatype = ANDT_SIGNED_SHORT;

	hdr.dime.bitpix = 8*sizeof(short);

	int min = im->min;
	int max = im->max;
	if (max > min) {
		hdr.dime.glmin = min;
		hdr.dime.glmax = max;
	}
	else {
		hdr.dime.glmin = 0;
		hdr.dime.glmax = 0;
	}

	/*-- write header --*/

	fp = fopen(hdr_name , "wb");
	if (fp == NULL) {
		fprintf(stderr, "** Can't open file %s for output!\n", hdr_name);
		return false;
	}
	fwrite(&hdr, sizeof(struct dsr), 1, fp);
	fclose(fp);

	if (ip == NULL)      /* 30 Sep 2002: skip .img if data not present */
		return false;

	/*-- write image --*/

	fp = fopen(img_name, "wb");
	if (fp == NULL) {
		fprintf(stderr, "** Can't open file %s for output!\n", img_name);
		return false;
	}
	fwrite(ip, len, sizeof(short), fp);
	fclose(fp);

	return true;
}


