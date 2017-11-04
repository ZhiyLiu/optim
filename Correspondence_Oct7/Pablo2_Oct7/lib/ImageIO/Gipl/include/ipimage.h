/*
 * image.h
 * 
 * author:  msturm 
 * created: 14 Feb 1997
 * changes: 
 */

// Modified to work under Visual C++ on MS Windows systems,
// Nov 30, 2004, A. G. Gash.  Function ipOpenImage was
// changed to open files in binary mode.


#ifndef __IP_IMAGE_H__
#define __IP_IMAGE_H__

#include <stdio.h>
#include <string.h>
#include <errno.h>
#include <signal.h>
#include <stdlib.h>
#include "macros.h"
#include "ipmatrix.h"


/* error numbers */
enum _ipimage_errors { IP_IE_UNSUPDEPTH = 256,
		       IP_IE_UNKNOWNFORMAT,
		       IP_IE_NOFILENAME,
		       IP_IE_UNSUPDIM,
		       IP_IE_UNSUPTYPE };

typedef enum _ipimage_errors ipImageError;

/* commands used to handled compressed files */

#define GzipCommand  "|gzip -cf > %s"
#define GunzipCommand  "|gzip -cdfq %s"
#define CompressCommand  "|compress -c > %s"
#define UncompressCommand  "|uncompress -c %s"


/* shortened version of ImageMagicks Image type, 
   used in all image related and file io functions 
   values correspond to AVS definitions (see avs.h) if available */

typedef struct
{
  FILE
  *file;

  ipImageError
  status;       

  ipArray
  field;

  float
  *origin,        // origin in world coordinates
    *pixdims;       // mapping to real world dimensions [mm]

  char
  *filename;

  long int
  filesize;

  int
  pipe;

} ipImage;

typedef struct
{
  FILE
  *file;

  ipImageError
  status;       /* errno or internal errors, 0 = no error */

  ipArraybyte
  field;

  float         
  *origin,        // origin in world coordinates
    *pixdims;       // mapping to real world dimensions [mm]

  char
  *filename;

  long int
  filesize;

  int
  pipe;

} ipImagebyte;

typedef struct
{
  FILE
  *file;

  ipImageError
  status;       /* errno or internal errors, 0 = no error */

  ipArrayfloat
  field;

  float         
  *origin,        // origin in world coordinates
    *pixdims;       // mapping to real world dimensions [mm]

  char
  *filename;

  long int
  filesize;

  int
  pipe;

} ipImagefloat;


/* function prototypes */

int ipAllocateImageData(ipImage *image, const ipVectorint *dimsvec, 
			const ipDataType type);
void ipClearImageData(ipImage *image);

/* image structure handling */
ipImage* ipNewImage();
ipImage* ipDeepCopyImage(const ipImage *image);

void ipDeleteImage(ipImage *image);

// set / get specific members of image structure
void ipSetImageData(ipImage *image, const void *data);
void ipSetImageType(ipImage *image, const ipDataType type);
void ipSetImageOrigin(ipImage *image, const ipVectorfloat *origin);
void ipSetImagePixDims(ipImage *image, const ipVectorfloat *pixdims);
void ipSetImageDim2Dto3D(ipImage *image, int dimz);
void ipSetImageDim3Dto2D(ipImage *image);
void ipSetImageDim3Dto2D(ipImage *image,char orient);
ipVectorint*   ipGetImageDims(const ipImage *image);
void*          ipGetImageData(const ipImage *image);
ipDataType     ipGetImageType(const ipImage *image);
ipVectorfloat* ipGetImageOrigin(const ipImage *image);
ipVectorfloat* ipGetImagePixDims(const ipImage *image);

void ipSetImageFilename(ipImage *image, const char *filename);

/* file handling */
void ipOpenImage(ipImage *image, const char * mode);
void ipCloseImage(ipImage *image);

/* data conversion utilities */
int ipConvertImageData(ipImage *image, const int norm, const ipDataType type);
ipImage* ipConvertImage(const ipImage *image, const int norm, const ipDataType type);

#endif

