/*
 *
 * image.cc - image file i/o functions
 *
 * author:  msturm
 * created: 17 Feb 1997
 * changes:
 *          adapted by Martin Andreas Styner 3. Jul 97
 * 
 * based on ftools.c implementation by O. Henricsson & G. Danuser
 */

#include "ipimage.h"
#ifndef WIN32
#include <unistd.h>
#else
#include <io.h>
#endif

/*------------------------------------------------------------------------------*/
/* image structure handling routines */
/*------------------------------------------------------------------------------*/

/* ipImage constructor */
ipImage* ipNewImage()
{
  ipImage* image = NULL;

  if (!(image = (ipImage *) malloc(sizeof(ipImage)))) {
    fprintf(stderr,"Error: ipNewImage [%s, line %d]: memory allocation failed:",
	    __FILE__, __LINE__);
    perror("");
    exit(errno);
  }
  memset(image, 0, sizeof(ipImage));

  return(image);
}


/* ipImage deep copy constructor */
ipImage* ipDeepCopyImage(const ipImage *src) 
{
  ipImage* dest = ipNewImage();

  dest->file = NULL;
  dest->filename = NULL;

  if (ipCopyArrayData(&(dest->field), &(src->field))) {
    ipDeleteImage(dest);

    return(NULL);
  }

  if (!(dest->origin = (float *) malloc(src->field.ndim * sizeof(float)))) {
    fprintf(stderr,
		"Error: ipDeepCopyImage [%s, line %d]: memory allocation failed:",
	    __FILE__, __LINE__);
    perror("");
    ipDeleteImage(dest);

    exit(errno);
  }
  memcpy(dest->origin, src->origin, src->field.ndim * sizeof(float));

  if (!(dest->pixdims = (float *) malloc(src->field.ndim * sizeof(float)))) {
    fprintf(stderr,
		"Error: ipDeepCopyImage [%s, line %d]: memory allocation failed:",
	    __FILE__, __LINE__);
    perror("");
    ipDeleteImage(dest);

    exit(errno);
  }
  memcpy(dest->pixdims, src->pixdims, src->field.ndim * sizeof(float));

  ipSetImageFilename(dest, src->filename);

  return(dest);
}


/* allocates memory for the image data and initializes it to 0
 *
 * params: ipImage *image  pointer to initialised ipImage structure
 *         ipImageType type  type of data elements
 *         int        ndim   number of dimensions
 *         int         ...   respective dimensions [ndim]
 */
int ipAllocateImageData(ipImage *image, 
			const ipVectorint *dimsvec,
			const ipDataType type) 
{
  int retval = -1;

  if (!image) {
#ifdef DEBUG
    fprintf(stderr,
		"Warning: ipAllocateImageData: initialized ipImage structure required\n");
#endif
    return(-1);
  }

  if (retval = ipAllocateArrayData(&(image->field), dimsvec, type)) {
    image->status = IP_IE_UNSUPTYPE;

    return(retval);
  }

  if (!(image->origin = (float *) malloc(dimsvec->dims[0] * sizeof(float)))) {
    fprintf(stderr,
		"Error: ipAllocateImageData [%s, line %d]: memory allocation failed:",
	    __FILE__, __LINE__);
    perror("");
    ipDeleteImage(image);

    exit(errno);
  }
  // init image origin to origin of world coordinate system
  memset(image->origin, 0, dimsvec->dims[0] * sizeof(float));

  if (!(image->pixdims = (float *) malloc(dimsvec->dims[0] * sizeof(float)))) {
    fprintf(stderr,
		"Error: ipAllocateImageData [%s, line %d]: memory allocation failed:",
	    __FILE__, __LINE__);
    perror("");
    ipDeleteImage(image);

    exit(errno);
  }
  // init voxel dimensions to (hyper)cubic voxels with side length 1mm
  for (int i = 0; i < image->field.ndim; i++)
    image->pixdims[i] = 1.0;

  return(0);
}

/*
 * Utility functions to set specific members of image structure
 */

void ipSetImageData(ipImage *image, const void *data)
{
  ipSetArrayData(&(image->field), data);
}

void ipSetImageType(ipImage *image, const ipDataType type)
{
  image->field.datatype = type;
}

/* sets the data of the image to 0 */
void ipClearImageData(ipImage *image) 
{
  memset(image->field.data._void, 0, 
	  image->field.datasize * image->field.elementsize);
}

// increase dimensions from two dimension to three dimension
void ipSetImageDim2Dto3D(ipImage *image, int dimz)
{
  int *a = new int[3];
  image->field.ndim = 3;
  a[0] = image->field.dims[0];
  a[1] = image->field.dims[1];
  a[2]  = dimz;
  image->field.dims = a;
  a = NULL;
  image->field.datasize = 
	  image->field.dims[0] * image->field.dims[1] * image->field.dims[2];
}

// decrease dimensions from three dimension to two dimension
void ipSetImageDim3Dto2D(ipImage *image)
{
  image->field.dims[2] = 1;
  image->field.ndim = 2;
  image->field.datasize = image->field.dims[0] * image->field.dims[1];
}

void ipSetImageDim3Dto2D(ipImage *image,char orient)
{
  if (orient == 'Z') {
    image->field.dims[2] = 1;
  } else if (orient == 'Y') {
    image->field.dims[1] = image->field.dims[2];
  } else {
    image->field.dims[0] = image->field.dims[1];
    image->field.dims[1] = image->field.dims[2];
  }
  image->field.dims[2] = 1;
  image->field.ndim = 2;
  image->field.datasize = image->field.dims[0] * image->field.dims[1];
}

/* set origin in real world coordinates [mm] */
void ipSetImageOrigin(ipImage *image, const ipVectorfloat *origin)
{
#ifdef BOUNDS_CHECK
  if (image->field.ndim != origin->dims[0]) {
    fprintf(stderr,
		"Error: ipSetImageOrigin [%s, line %d]: dimensions mismatch: \n\
       (image->ndims = %d vs. origin->ndims = %d)\n",
	    __FILE__, __LINE__, image->field.ndim, origin->dims[0]);
    exit(-1);
  }
#endif

  memcpy(image->origin, origin->data, image->field.ndim * sizeof(float));
}


/* set pixel dimensions in real world units [mm] */
void ipSetImagePixDims(ipImage *image, const ipVectorfloat *pixdims)
{
#ifdef BOUNDS_CHECK
  if (image->field.ndim != pixdims->dims[0]) {
    fprintf(stderr,"Error: ipSetImagePixDims [%s, line %d]: dimensions mismatch: \n       (image->ndims = %d vs. pixdims->ndims = %d)\n",
	    __FILE__, __LINE__, image->field.ndim, pixdims->dims[0]);
    exit(-1);
  }
#endif

  memcpy(image->pixdims, pixdims->data, image->field.ndim * sizeof(float));
}


/* sets filename entry in ipImage structure */
void ipSetImageFilename(ipImage *image, const char *filename)
{
  ipCloseImage(image);

  if (image->filename) free(image->filename); /* ev. reuse of old struct */

  if (!(image->filename = strdup(filename))) {
    fprintf(stderr,"Error: ipSetImageFilename [%s, line %d]: strdup failed:",
	    __FILE__, __LINE__);
    perror("");
    exit(errno);
  }

}

/*
 * Utility functions to set specific members of image structure
 */

ipVectorint* ipGetImageDims(const ipImage *image)
{
  ipVectorint *tmp_veci = ipNewVectorint(image->field.ndim);
  ipSetVectorData((ipVector *) tmp_veci, image->field.dims);

  return(tmp_veci);
}

void* ipGetImageData(const ipImage *image)
{
  return(ipGetArrayData(&(image->field)));
}

ipDataType ipGetImageType(const ipImage *image)
{
  return(image->field.datatype);
}

ipVectorfloat* ipGetImageOrigin(const ipImage *image)
{
  ipVectorfloat *tmp_vecf = ipNewVectorfloat(image->field.ndim);
  ipSetVectorData((ipVector *) tmp_vecf, image->origin);

  return(tmp_vecf);
}

ipVectorfloat* ipGetImagePixDims(const ipImage *image)
{
  ipVectorfloat *tmp_vecf = ipNewVectorfloat(image->field.ndim);
  ipSetVectorData((ipVector *) tmp_vecf, image->pixdims);

  return(tmp_vecf);
}

/*
 * Destructors
 */


/* ipImage destructor */
void ipDeleteImage(ipImage *image)
{
  if (image) {
    if (image->file)
      /* close file first ! */
      ipCloseImage(image);

    ipDeleteArrayData(&(image->field));

    free(image->origin);
    free(image->pixdims);
    free(image->filename);
    free(image);
  }
}


/*------------------------------------------------------------------------------*/
/* file handling routines */
/*------------------------------------------------------------------------------*/

/* opens an image file for reading or writing (depending on type = ['r'|'w'])
 * and returns a FILE* handle
 *
 * the image structure has to be initalized by a call to ipNewImage
 *
 * handles gzipped or compressed files
 * adapted from ImageMagick's OpenImage
 */

void ipOpenImage(ipImage* image, const char * mode)
{
  char *act_name = NULL;

  if (!image->filename) {
#ifdef DEBUG
    fprintf(stderr, "Warning: ipOpenImage: empty filename\n");
#endif
    image->file = NULL;
    image->status = IP_IE_NOFILENAME;
    return;
  }

  if (((int) strlen(image->filename) > 3) &&
      (strcmp(image->filename+strlen(image->filename)-3,".gz") == 0))
    {
      /*
       * Uncompress/compress image file with GNU compress utilities.
       */
      if (mode[0] == 'r') {
	act_name = 
		new char[strlen(GunzipCommand) + strlen(image->filename) + 1];
	(void) sprintf(act_name,GunzipCommand,image->filename);
      } else {
	act_name = 
		new char[strlen(GzipCommand) + strlen(image->filename) + 1];
	(void) sprintf(act_name,GzipCommand,image->filename);
      }
    }
  else
    if (((int) strlen(image->filename) > 2) &&
	(strcmp(image->filename+strlen(image->filename)-2,".Z") == 0)) {
      /*
       * Uncompress/compress image file with UNIX compress utilities.
       */
      if (mode[0] == 'r') {
	act_name = new char[strlen(UncompressCommand) + 
			    strlen(image->filename) + 1];	  
	(void) sprintf(act_name,UncompressCommand,image->filename);
      }	else {
	act_name = new char[strlen(CompressCommand) + 
			    strlen(image->filename) + 1];	  
	(void) sprintf(act_name,CompressCommand,image->filename);
      }
    }

  image->pipe   = 0;
  image->status = (ipImageError) 0;
  image->file   = (FILE *) NULL;
  if (act_name) {
    /*
      Pipe image to or from a system command.
    */
	  //We won't allow this on windows
#ifndef WIN32
    if (mode[0] == 'w')
      (void) signal(SIGPIPE,SIG_IGN);

// the mode passed to popen on some linux systems needs to be "r" or "w"
    char singleCharMode[2];
	singleCharMode[0] = mode[0];
	singleCharMode[1] = '\0';
//

    if (!(image->file = (FILE *) popen(act_name+1, singleCharMode))) {
//#ifdef DEBUG
      fprintf(stderr,
		  "Warning: ipOpenImage [%s, line %d]: popen of \"%s\" failed:",
	      __FILE__, __LINE__, image->filename);
      perror("");
//#endif	
      image->status = (ipImageError) errno;
      return;
    }
    image->pipe = 1;
#endif

  } else {

    if (!(image->file = (FILE *) fopen(image->filename, mode))) {
#ifdef DEBUG
      fprintf(stderr,
		  "Warning: ipOpenImage [%s, line %d]: fopen of \"%s\" failed:",
	      __FILE__, __LINE__, image->filename);
      perror("");
#endif	
      image->status = (ipImageError) errno;
      return;
    } else {
      (void) fseek(image->file,0L,2);
      image->filesize=ftell(image->file);
      (void) fseek(image->file,0L,0);
    }

  }
  delete [] act_name;
}


/* Close image file previously opened by ipOpenImage */

void ipCloseImage(ipImage *image)
{
  if (image == (ipImage *) NULL)
    return;
  if (image->file == (FILE *) NULL)
    return;
  image->status = (ipImageError) ferror(image->file);
  //not allowed for windows

  if (image->pipe){
#ifndef WIN32
    (void) pclose(image->file);
#endif
  }else{
    if ((image->file != stdin) && (image->file != stdout))
      (void) fclose(image->file);
  }

  image->file=(FILE *) NULL;
}


/*------------------------------------------------------------------------------*/
/* data conversion utilities */
/*------------------------------------------------------------------------------*/

/* converts data of image to ipDataType "type"
 * (only IP_BYTE and IP_FLOAT and IP_SHORT supported)
 *
 * if norm = 1, data is scaled to the range of 0..255
 * no error is given when original-type and type are equal
 */
int ipConvertImageData(ipImage *image, const int norm, 
					   const ipDataType type)
{

  if ((image->field.datatype == type) && !norm)
    return(0);

  switch (type) {
  case IP_BYTE:
    if (ipArrayDatatoByte(&(image->field), norm)) {
      ipDeleteArrayData(&(image->field));

      return(-1);
    }
    break;
  case IP_SHORT:
    if (ipArrayDatatoShort(&(image->field), norm)) {
      ipDeleteArrayData(&(image->field));

      return(-1);
    }
    break;
  case IP_FLOAT:
    if (ipArrayDatatoFloat(&(image->field), norm)) {
      ipDeleteArrayData(&(image->field));

      return(-1);
    }
    break;
  default:
    return(-1);
  }

  return(0);
}

ipImage* ipConvertImage(const ipImage *src, const int norm, 
			const ipDataType type)
{
  ipImage* dest = ipDeepCopyImage(src);

  dest->file = NULL;

  if (ipConvertImageData(dest, norm, type)) {
    ipDeleteImage(dest);

    return(NULL);
  }

  return(dest);
}


