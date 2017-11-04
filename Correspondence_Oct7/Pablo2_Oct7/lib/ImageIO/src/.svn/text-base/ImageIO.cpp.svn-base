
#include <iostream>
#include <string>
#include <fstream>
#include <stdio.h>
#include <unistd.h>
#include <fcntl.h>
#include "BasicException.h"
#include "ImageStruct.h"
#include "Analyze.h"
#include "giplio.h"
#include "metaImage.h"
#include "AllImageIO.h"
#include "ImageIO.h"
#include "gen.h"
#include "plan_im.h"
#include "plan_file_io.h"
#include "extbm.h"
#include "libplanio.h"
#include "utility.h"
#include "ControlParms.h"

#ifndef WIN32
#define O_BINARY 0x0
#endif


using namespace std;

const double fuzz = 1e-5;    // Floating (not double) point comparison fuzz

extern void get_date(PDATE *date);


//   ImageIO allow to load images of the following format :
//         -> Analyze (divided into two files, the header (*.hdr) with 
//				the information about the image, and the data (*.img)) 
//         -> Meta (*.mha or *.mhd + *.raw)     
//         -> Gipl (*.gipl)
//         -> PlanIM (no specific name, most of the time you have plan_im 
//				in the filename, can be *.pim)
//         -> Dicom (one file per slice, no specific name) 
//
//   Byte swapping functions have written for analyze and gipl images
// We assume (it has been tested ! ) that PlanIm, Meta and Dicom are already
// handling byte order.



/*****************************************************************************/


// Guess the format of an image depending 
//    1. on its extension
//           if ( *.hdr ) return analyze
//           if ( *.gipl ) return gipl
//           if ( *.[mha|mhd] ) return meta
//           if ( *.dcm ) return dicom
//           if ( *.pim ) return planIM
//   2. if none of the above match, check if there is "plan_im" in 
//      the filename, if so return planIM
//   3. if there is a * in the filename the image is considered as a Dicom


ImageIO::ImageType ImageIO::guessImageFormat(string & filename)
{
	string extension = filename;
	int pos = extension.find_last_of(".");
        if (pos != string::npos)
	    extension.erase(extension.begin(), extension.begin() + pos + 1 );

	/*if (extension.compare("img")==0 ||
		extension.compare("hdr")==0) { 
			return analyze;
	}*/
	if (extension.compare("hdr")==0) { 
		return analyze;
	}
	if (extension.compare("mha")==0 || extension.compare("mhd")==0)
		return meta;  
	if (extension.compare("gipl")==0) return gipl; 
	if (extension.compare("dcm")==0) return dicom;
//	if (extension.compare("img")==0) return dicom;
	if (extension.compare("pim")==0) return planIM;
	if (extension.compare("raw3")==0) return raw3;


	string::size_type posPlan_im = filename.find("plan_im",0);
	if (posPlan_im != string::npos) return planIM;

	string::size_type posStar = filename.find("*",0);
	if (posStar != string::npos) return dicom;

	string::size_type posGiplGz = filename.find("gipl.gz",0);
	if (posGiplGz != string::npos) return gipl;
	else cout << "This file isn't a gipl.gz:  [" << filename << "]" << endl;

	return unknown;
}


// Load an image guessing its format
// With PlanIM, take the default parameter by guessParameters(...))
// With Dicom, we browse the different iamges of the volume (nb_DICOM)
// and we select the one with the maximun nb of slice in Z.
// We assume that the "interesting" part is the one with the max nb of slice
// We throw an exception if the file format is unknown.

void ImageIO::loadThisImage(string & filename, ImageStruct & imStr,
	bool headerOnly, ImageType extension) throw (BasicException)
{

	if (extension == unknown){
		extension = guessImageFormat(filename);
	}

	switch(extension){
		case gipl :
			loadGipl(filename, imStr, headerOnly);
			break;
		case analyze :
			loadAnalyze(filename, imStr, headerOnly);
			break;
		case meta :
			loadMeta(filename, imStr, headerOnly);
			break;
		case planIM :
			loadPlanim(filename, imStr, headerOnly);
			break;
		case dicom :
			throw bException("Not implemented");
			break;
		case raw3:
			throw bException("Load Failed: Raw3 input not implemented in ImageIO");
			break;
		case unknown :
			throw bException("Load Failed: Unknown file format");
			break;
	}
}

// Save an image guessing its format
// Don't work with PlanIM and Dicom because of the user interface part
// We throw an exception in such cases.

void ImageIO::saveImage(string & filename, ImageStruct & imStr){

  ImageType extension;
  extension = guessImageFormat(filename);

  switch(extension){
	  case gipl :
		saveGipl(filename, imStr);
		break;
	  case analyze :
		saveAnalyze(filename, imStr);
		break;
	  case meta :
		saveMeta(filename, imStr);
		break;
	  case planIM :
		savePlanim(filename, imStr);
		break;
	  case dicom :
		throw bException("Not implemented");
		break;
	  case raw3:
		throw bException("Save Failed: Raw3 output not implemented in ImageIO");
		break;
	  case unknown :
		throw  bException("Save failed: Unknown file format");
		break;
  }
}



/**********************************************************************/
/**********************************************************************/
/*                                                                    */
/*                              ANALYZE                               */
/*                                                                    */
/**********************************************************************/
/**********************************************************************/

void ImageIO::loadAnalyze(string & filename, ImageStruct & imStr,
	bool /*headerOnly*/)
{
	Analyze a;

	if (! a.load(filename, imStr))
		throw bException("Unable to load analyze file");
}

void ImageIO::saveAnalyze(string & filename, ImageStruct & imStr)
{
	Analyze a;

	if (! a.save(filename, imStr))
		throw bException("Unable to save analyze file");
}



/**********************************************************************/
/**********************************************************************/
/*                                                                    */
/*                               META                                 */
/*                                                                    */
/**********************************************************************/
/**********************************************************************/


// LoadMeta load a meta image in 3 steps :
//    1. We copy the file in a MetaImage object
//    2. Fill the header from this object
//    3. Fill the data (imStr) from this object

void ImageIO::loadMeta(string & filename, ImageStruct & imStr,
                       bool /*headerOnly*/)
{
    int i;

    // Read the file in meta (class MetaImage)   
    MetaImage meta;
    if (! meta.Read(filename.c_str()))
        throw bException("Unable to open meta file");

    // Get header data
    imStr.dims[0] = meta.DimSize(0);
    imStr.dims[1] = meta.DimSize(1);
    imStr.dims[2] = meta.DimSize(2);

    int len = imStr.dims[0]*imStr.dims[1]*imStr.dims[2];

    const double * o = meta.Origin();
    imStr.origin[0] = (float) o[0];
    imStr.origin[1] = (float) o[1];
    imStr.origin[2] = (float) o[2];

    const float * es = meta.ElementSize();
    imStr.spacing[0] = es[0];
    imStr.spacing[1] = es[1];
    imStr.spacing[2] = es[2];

    // Then we fill imStr depending on the data type
    MET_ValueEnumType t = meta.ElementType();
    if (imStr.voxels != NULL) {
        delete [] imStr.voxels;
        imStr.voxels = NULL;
    }

    if (meta.ElementData() == NULL) {
        imStr.len = 0;
        return;
    }
    switch(t) {
        case MET_UCHAR:
            imStr.voxels = new u_short[len];
            for (i = 0; i < meta.Quantity(); i++)
                imStr.voxels[i] = (unsigned short) ((char *) meta.ElementData())[i];
            break;

        case MET_FLOAT:
            if (! meta.ConvertIntensityDataToElementData(MET_USHORT))
                throw bException(
                "Could not convert MET_FLOAT voxels to type MET_USHORT");
            // pixels = meta.ElementData();
            imStr.voxels = new u_short[len];
            memcpy(imStr.voxels, meta.ElementData(), sizeof(u_short)*meta.Quantity());
            break;

        case MET_SHORT:
            if (! meta.ConvertIntensityDataToElementData(MET_USHORT))
                throw bException(
                "Could not convert MET_SHORT voxels to type MET_USHORT");
            // Fall through

        case MET_USHORT:
            imStr.voxels = new u_short[len];
            memcpy(imStr.voxels, meta.ElementData(), sizeof(u_short)*meta.Quantity());
            break;

        case MET_CHAR:
            throw bException("MET_CHAR not supported");
            break;

        case MET_INT:
            throw bException("MET_INT not supported");
            break;

        case MET_UINT:
            throw bException("MET_UINT not supported");
            break;

        case MET_DOUBLE:
            throw bException("MET_DOUBLE not supported");
            break;

        default:
            throw bException("Meta Image Load Failed: Unknown file format");
            break;
    }
    imStr.len = len;

    imStr.min = (int) meta.ElementMin();
    imStr.max = (int) meta.ElementMax();
}

// To save an image as a Meta, we create an object of the class
// MetaImageClass and write it to a file.

void ImageIO::saveMeta(string & filename, ImageStruct & imStr)
{
	double Meta_offset[3];
//	float Meta_orientation[6];
        bool compress;

	// We create meta (class MetaImageClass) and add other info to it
	MetaImage meta(imStr.dims[0], imStr.dims[1], imStr.dims[2],
                imStr.spacing[0], imStr.spacing[1], imStr.spacing[2], 
		MET_USHORT, 1, imStr.voxels);

	compress = globalControl->readBool(CompressImages);
	meta.CompressedData(compress);

	Meta_offset[0] = imStr.origin[0];
	Meta_offset[1] = imStr.origin[1];
	Meta_offset[2] = imStr.origin[2];
	meta.Position(Meta_offset);
//  Meta_orientation[0] = (float)imStr.getOrientation(); 
//	meta.Orientation((const float*)&Meta_orientation);
	if (imStr.max > imStr.min) {
		meta.ElementMin(imStr.min);
		meta.ElementMax(imStr.max);
	}

	// Write the meta file
	if (! meta.Write(filename.c_str()))
		throw bException("Unable to save meta file");
}




/**********************************************************************/
/**********************************************************************/
/*                                                                    */
/*                               GIPL                                 */
/*                                                                    */
/**********************************************************************/
/**********************************************************************/


// LoadGipl load a gipl image in 3 steps :
//    1. We copy the file in a ipGiplImage object
//    2. Fill the header from this object
//    3. Fill the data (imStr) from this object


void ImageIO::loadGipl(string & filename, ImageStruct & imStr,
	bool /*headerOnly*/)
{
	// We write the info of the file in gipl
	// Byte order issue is taking care of in this function : ipReadGipl()

	ipGiplImage *gipl = ipReadGipl(filename.c_str());  
	if (gipl == NULL) {
		throw bException("Unable to open gipl file");
	}

	// Then we get gipl infos

	ipVectorint* dims = ipGetImageDims(gipl->image);
	ipDataType gtype = ipGetImageType(gipl->image);
	ipVectorfloat* pixdims = ipGetImagePixDims(gipl->image);

	if (dims->dims[0] != 3) {
		throw bException("Gipl image must be 3 dimensions");
	}	

	// And write it in our header

	imStr.dims[0] = dims->data[0];
	imStr.dims[1] = dims->data[1];
	imStr.dims[2] = dims->data[2];

	int len = imStr.dims[0]*imStr.dims[1]*imStr.dims[2];

	imStr.origin[0] = gipl->header->origin[0];
	imStr.origin[1] = gipl->header->origin[1];
	imStr.origin[2] = gipl->header->origin[2];

	imStr.spacing[0] = pixdims->data[0];
	imStr.spacing[1] = pixdims->data[1];
	imStr.spacing[2] = pixdims->data[2];

	// Then we fill imStr depending on the data type
	if (imStr.voxels != NULL) {
		delete [] imStr.voxels;
		imStr.voxels = NULL;
	}

	switch(gtype) {
		case IP_BYTE:
			imStr.voxels = new u_short[len];
			memcpy(imStr.voxels, gipl->image->field.data._byte, len);
			break;
		case IP_FLOAT:
			ipConvertImageData(gipl->image, 0, IP_SHORT);
			// Fall through
		case IP_SHORT:
			imStr.voxels = new u_short[len];
			memcpy(imStr.voxels, gipl->image->field.data._void,
				len*sizeof(u_short));
			break;
		default:
			break;
    }
	imStr.len = len;

	imStr.min = 0;
	imStr.max = 0;

	ipDeleteGipl(gipl);
	ipDeleteVectorint(dims);
	ipDeleteVectorfloat(pixdims);
}


// To save an image as a Gipl, we create an object of the class
// ipGiplImage and write it to a file.

void ImageIO::saveGipl(string & filename, ImageStruct & imStr)
{

	// Create gipl
	ipGiplImage *gipl = ipNewGipl();
	gipl->image = ipNewImage();
	gipl->header = ipNewGiplHeader();
	if (!gipl || !gipl->image)
		throw bException("Unable to allocate gipl image memory");

	ipSetImageFilename(gipl->image, filename.c_str());
	ipVectorint *dimsvec = ipNewVectorint(3);

	// Set the dimensions
	dimsvec->data[0] = imStr.dims[0];
	dimsvec->data[1] = imStr.dims[1];
	dimsvec->data[2] = imStr.dims[2];
	int len = imStr.len;

	ipDataType gtype;

	gtype = IP_SHORT;
	ipAllocateImageData(gipl->image, dimsvec, gtype);
	memcpy(gipl->image->field.data._void, imStr.voxels,
		len*sizeof(u_short));

	// Set the origin and voxel spacings
	gipl->header->origin[0] = imStr.origin[0];
	gipl->header->origin[1] = imStr.origin[1];
	gipl->header->origin[2] = imStr.origin[2];

	gipl->image->pixdims[0] = imStr.spacing[0];
	gipl->image->pixdims[1] = imStr.spacing[1];
	gipl->image->pixdims[2] = imStr.spacing[2];

	// Write the file
	if (-1 == ipWriteGipl(gipl, 0)) {
		ipDeleteGipl(gipl);
		ipDeleteVectorint(dimsvec);
		throw bException("Unable to save gipl file");
	}

	ipDeleteGipl(gipl);
	ipDeleteVectorint(dimsvec);

}




/**********************************************************************/
/**********************************************************************/
/*                                                                    */
/*                             PLAN_IM                                */
/*                                                                    */
/**********************************************************************/
/**********************************************************************/

bool fuzzyEq(double a, double b)
{
    double diff = a - b;
    if (diff < 0.0)
        diff = -diff;
    if (diff > fuzz)
        return false;
    return true;
}

void ImageIO::loadPlanim(string & filename, ImageStruct & imStr,
	bool headerOnly)
{
    plan_im_header hdr;
    int i;
    double pixelZSize;
    bool uniform;
    int slice_size;
    short * data;

    int fdes = open(filename.c_str(), O_RDONLY|O_BINARY, 0);
    if (fdes < 0)
	throw bException("Unable to load plan_im file");

    if (read_image_header(fdes, &hdr)) {
		(void) close(fdes);
		throw bException("Can't read header");
    }

    slice_size = hdr.x_dim * hdr.y_dim;
	if (headerOnly)
		data = NULL;
	else {
		int data_size = slice_size * hdr.slice_count;
		data = new short[data_size];
		for (i = 0; i < hdr.slice_count; i++) {
			if (read_scan_xy(fdes, data + i*slice_size, &hdr, i))
			{
				(void) close(fdes);
				throw bException("Can't read image data");
			}
		}
		(void) close(fdes);
	}

    imStr.dims[0] = hdr.x_dim;
    imStr.dims[1] = hdr.y_dim;
    imStr.dims[2] = hdr.slice_count;

    imStr.len = slice_size*hdr.slice_count;
    imStr.max = hdr.max;
    imStr.min = hdr.min;

    imStr.origin[0] = hdr.pixel_to_patient_TM[3][0];
    imStr.origin[1] = hdr.pixel_to_patient_TM[3][1];
    imStr.origin[2] = hdr.per_scan[0].z_position;

    uniform = true;
    if (hdr.slice_count < 2)
        pixelZSize = 1.0;
    else {
        int numgaps = hdr.slice_count - 1;
        pixelZSize = (hdr.per_scan[numgaps].z_position 
            - hdr.per_scan[0].z_position) / numgaps;
	if (fuzzyEq(pixelZSize, 0.0)) pixelZSize = 1.0;
	if (hdr.slice_count > 2)
	    for (i = 1; i < hdr.slice_count; i++) {
		double size = (double) hdr.per_scan[i].z_position
		    - (double) hdr.per_scan[i - 1].z_position;
		if (! fuzzyEq(size, pixelZSize))
		    uniform = false;    // Slice thickness varies
	    }
    }

    if (! uniform)
        cout << "Warning: file " << filename <<
            "\n    appears to contain slices of variable thickness or in non-ascending order."
            << endl;

    imStr.spacing[0] = hdr.pixel_to_patient_TM[0][0];
    imStr.spacing[1] = hdr.pixel_to_patient_TM[1][1];   // This is normally negative
    imStr.spacing[2] = (float) pixelZSize;
    imStr.voxels = (unsigned short *) data;
	imStr.dataIsShort = true;
	if (hdr.pixel_type == mri_number)
		imStr.modality = 3;    // MRI
	else if (hdr.pixel_type == ct_number) {
		imStr.modality = 1;    // CT
		if (imStr.min >= 0 && imStr.max > 1024 && imStr.max <= 4095)
		    imStr.modality = 2;    // SHIFTED_CT
	}
	else
		imStr.modality = 0;    // Unknown modality
}

void ImageIO::savePlanim(string & filename, ImageStruct & imStr)
{
    int fdes, i;
    plan_im_header hdr;
    float pos;
    int slicecnt;
	bool compress;

    hdr.max = imStr.max;
    hdr.min = imStr.min;

    hdr.x_size = fabs(imStr.spacing[0]);
    hdr.y_size = fabs(imStr.spacing[1]);
    hdr.pixel_size = hdr.x_size;

    hdr.resolution = imStr.dims[0];
    hdr.x_dim = imStr.dims[0];
    hdr.y_dim = imStr.dims[1];
    hdr.slice_count = imStr.dims[2];
	int slize_size = hdr.x_dim*hdr.y_dim*sizeof(short);

    hdr.unit_number[0] = '\0';
    hdr.patient_name[0] = '\0';
    strcpy(hdr.comment, "Generated by Pablo2.");
    get_date(&hdr.date);

    for (i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++)
            hdr.pixel_to_patient_TM[i][j] = 0.0f;
    hdr.pixel_to_patient_TM[0][0] = imStr.spacing[0];
	// This already will be negative for planIm input files
    hdr.pixel_to_patient_TM[1][1] = imStr.spacing[1];
//    hdr.pixel_to_patient_TM[3][0] = -(imStr.dims[0] >> 1)*imStr.spacing[0];
    hdr.pixel_to_patient_TM[3][0] = imStr.origin[0];
    hdr.pixel_to_patient_TM[3][1] = imStr.origin[1];
    hdr.pixel_to_patient_TM[3][3] = 1.0f;
    hdr.table_height = (int) (0.5 +
        ((hdr.y_dim - 1) - imStr.origin[1]/hdr.y_size));

    hdr.machine_id = bogus_scanner;
    hdr.patient_position = bogus_position;
    hdr.pixel_type = bogus_pixel;
    hdr.whats_first = bogus_entry;
    hdr.time_count = 1;

    switch (imStr.modality) {
        default:
        case 0: hdr.pixel_type = bogus_pixel;
                break;
        case 1:
        case 2: hdr.pixel_type = ct_number;
                break;
        case 3: hdr.pixel_type = mri_number;
                break;
    }

//    pos = -(imStr.dims[2] >> 1)*imStr.spacing[2];
    pos = imStr.origin[2];
    for (i = 0; i < imStr.dims[2]; i++) {
        hdr.per_scan[i].scan_number = i;
        hdr.per_scan[i].z_position = pos;
        pos += imStr.spacing[2];
        hdr.per_scan[i].gantry_angle = 0.0f;
        hdr.per_scan[i].table_angle = 0.0f;
        hdr.per_scan[i].nbytes = slize_size;
    }

	hdr.cut = axial;
	compress = globalControl->readBool(CompressImages);
	if (compress)
		set_image_compression(1);
	else
		set_image_compression(0);	// The libplanio default is 1

    fdes = open(filename.c_str(), O_RDWR|O_TRUNC|O_CREAT|O_BINARY, 0666);
    if (fdes < 0) {
		throw bException("Cannot open file for writing");
			return;
    }

    if (write_image_header(fdes, &hdr)) {
		throw bException("Can't write output image header");
		return;
    }
    slicecnt = imStr.dims[0]*imStr.dims[1];
    PIXELTYPE * voxels = (short *) imStr.voxels;
    for (i = 0; i < hdr.slice_count; i++) {
		if (write_scan_xy(fdes, voxels + i*slicecnt,
			   &hdr, i, false)) {
			throw bException("Can't write image scans");
			return;
		}
    }
    lseek(fdes, 0, SEEK_SET);
    if (write_image_header(fdes, &hdr)) {
		throw bException("Can't write output image header");
		return;
    }

    close(fdes);
    return;
}

