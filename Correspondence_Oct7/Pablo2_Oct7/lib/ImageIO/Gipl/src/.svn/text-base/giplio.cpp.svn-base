/*
 * giplio.c - gipl image file handling functions
 * (Guys Hospital London, Image Processing Library File format)
 *
 * author:  msturm
 * created: implemented as avs read module, 05 Nov 1996 
 * changes: included in ip lib, 18 Feb 1997 
 *
 */

/* Modified 20 May 2003
 * changes : handle byte order (byte swapping) using BinaryIO.h
 * author : FX Kuhn
 */


#include <iostream> 
#include "giplio.h"
#include <BinaryIO.h>

using namespace std;


GIPL_HD *ipNewGiplHeader(void)
{
  GIPL_HD* header = NULL;
  const char     *dummy_text = "Generated header";

  if (!(header = (GIPL_HD *) malloc(sizeof(GIPL_HD)))) {
    fprintf(stderr,"Error: ipNewGiplHeader [%s, line %d]: memory allocation failed:",
	    __FILE__, __LINE__);
    perror("");
    exit(errno);
  }
  memset(header, 0, sizeof(GIPL_HD));
  if (sizeof(dummy_text) < sizeof(header->line1))
    strcpy(header->line1, dummy_text);
  else
    strncpy(header->line1, dummy_text, sizeof(header->line1));
  header->magic_number = GIPL_MAGIC_NUMBER;

  return(header);
}

ipGiplImage* ipNewGipl(void)
{
  ipGiplImage* image = NULL;

  if (!(image = (ipGiplImage *) malloc(sizeof(ipGiplImage)))) {
    fprintf(stderr,"Error: ipNewGipl [%s, line %d]: memory allocation failed:",
	    __FILE__, __LINE__);
    perror("");
    exit(errno);
  }
  memset(image, 0, sizeof(ipGiplImage));

  return(image);
}

void ipDeleteGipl(ipGiplImage *in)
{
  if (in->header)
    free(in->header);
  if (in->image)
    ipDeleteImage(in->image);
  free(in);
}

ipGiplImage* ipReadGipl(const char *filename)
{
  int i, tmpi;
  unsigned char   *colormap = NULL, 
    *rgb = NULL;

  // Flag for byte swapping
  bool ByteSwapNeeded = false ;


  ipVectorint   *dimsvec = NULL;
  ipVectorfloat *floatvec;
  ipGiplImage   *gipl = ipNewGipl();
  ipDataType     imaget = IP_SHORT;

  gipl->image = ipNewImage();
  gipl->header = ipNewGiplHeader();

  ipImage       *image = gipl->image;
  GIPL_HD       *header = gipl->header;

  ipSetImageFilename(image, filename);

  ipOpenImage(image, "rb");

  if ((!image->file) || image->status) {

    ipDeleteImage(image);

    return(NULL);
  }

  /* thanks to the incredibly stupid design of the gipl header
     we have to parse the file header entry by entry.... */

  fread(header->dim,(size_t) sizeof(header->dim), (size_t) 1, image->file);

  BinaryIO binIO;

  // Here we check if the dimensions are valid, i.e, if
  // byte swapping is necessary. In that case we switch the flag
  // ByteSwapNeeded to true.

  if ( (header->dim[0] <= 0 || header->dim[0] > 1000) ||
       (header->dim[1] <= 0 || header->dim[1] > 1000) ||
       (header->dim[2] <= 0 || header->dim[2] > 1000) ){

    ByteSwapNeeded = true;
    // We "byte swap" the dimensions
    for (int i = 0 ; i < 4 ; i ++){
      header->dim[i]=binIO.swabShort(header->dim[i]);
    }  
  }
//  cout << "ByteSwapNeeded = " << ByteSwapNeeded << endl;

  if (header->dim[2] < 2) tmpi = 2;
  else if (header->dim[3] < 2) tmpi = 3;
  else tmpi = 4;

  dimsvec = ipNewVectorint(tmpi);
  floatvec = ipNewVectorfloat(tmpi);

  for (i = 0; i < dimsvec->dims[0]; i++)
    dimsvec->data[i] = (int) header->dim[i];

  fread(&header->image_type,(size_t) sizeof(header->image_type), (size_t) 1,image->file);

  // We "byte swap" the image_type
  if (ByteSwapNeeded  ){
    header->image_type = binIO.swabShort( header->image_type );
  }

  switch(header->image_type) {
  case IT_BINARY:
    break;
  case IT_CHAR:
  case IT_U_CHAR:
    imaget = IP_BYTE;
    break;
  case IT_SHORT:
  case IT_U_SHORT:
    imaget = IP_SHORT;
    break;
  case IT_INT:
  case IT_U_INT:
    imaget = IP_INT;
    break;
  case IT_FLOAT:
    imaget = IP_FLOAT;
    break;
  case IT_DOUBLE:
    imaget = IP_DOUBLE;
    break;
  default:
#ifdef DEBUG
    fprintf(stderr, "Warning: ipReadGipl: unsuported data type: %d\n",
	    header->image_type);
#endif
    cout<<"Warning: ipReadGipl: unsuported data type:"<<endl;
    ipDeleteImage(image);

    return(NULL);
  }
  ipAllocateImageData(image, dimsvec, imaget);


  // pixel dimensions
  fread(header->pixdim,(size_t) sizeof(header->pixdim), (size_t) 1,image->file);

  // We "byte swap" the pixel dimensions
  if (ByteSwapNeeded){
    for (int i = 0 ; i < 4 ; i ++ ){
      header->pixdim[i]=binIO.swabFloat(header->pixdim[i]);
    }
  }

  floatvec->dims[0] = dimsvec->dims[0];
  for (i = 0; i < floatvec->dims[0]; i++)
    floatvec->data[i] = (float) header->pixdim[i];
  image->field.ndim = dimsvec->dims[0];
  ipSetImagePixDims(image, floatvec);

  fread(header->line1,(size_t) sizeof(header->line1), (size_t) 1,image->file);
  fread(header->matrix,(size_t) sizeof(header->matrix), (size_t) 1,image->file);
  fread(&header->flag1,(size_t) sizeof(header->flag1), (size_t) 1,image->file);
  fread(&header->flag2,(size_t) sizeof(header->flag2), (size_t) 1,image->file);
  fread(&header->min,(size_t) sizeof(header->min), (size_t) 1,image->file);
  fread(&header->max,(size_t) sizeof(header->max), (size_t) 1,image->file);

  // Byte Swapping
  if (ByteSwapNeeded){
    for (int i = 0 ; i < 20 ; i ++ ){
      header->matrix[i]=binIO.swabFloat(header->matrix[i]);
    }
     header->min=binIO.swabDouble(header->min);
     header->max=binIO.swabDouble(header->max);
  }

  // image origin
  fread(header->origin,(size_t) sizeof(header->origin), (size_t) 1,image->file);
  // We "byte swap" the origin
  if (ByteSwapNeeded){
    for (int i = 0 ; i < 4 ; i ++ ){
      header->origin[i]=binIO.swabDouble(header->origin[i]);
    }
  }
  for (i = 0; i < floatvec->dims[0]; i++)
    floatvec->data[i] = (float) header->origin[i];
  ipSetImageOrigin(image, floatvec);

  fread(&header->pixval_offset,(size_t) sizeof(header->pixval_offset), (size_t) 1,
	image->file);
  fread(&header->pixval_cal,(size_t) sizeof(header->pixval_cal), (size_t) 1,image->file);
  fread(&header->user_def1,(size_t) sizeof(header->user_def1), (size_t) 1,image->file);
  fread(&header->user_def2,(size_t) sizeof(header->user_def2), (size_t) 1,image->file);
  fread(&header->magic_number,(size_t) sizeof(header->magic_number), (size_t) 1,
	image->file);

  // Byte Swapping
  if (ByteSwapNeeded){
    header->pixval_offset=binIO.swabFloat(header->pixval_offset);
    header->pixval_cal=binIO.swabFloat(header->pixval_cal);
    header->user_def1=binIO.swabFloat(header->user_def1);
    header->user_def2=binIO.swabFloat(header->user_def2);
    header->magic_number=binIO.swabInt(header->magic_number);
  }
/*
   cout<<"header->pixdim[0] "<<header->pixdim[0]<<endl;
   cout<<"header->pixdim[1] "<<header->pixdim[1]<<endl;
   cout<<"header->pixdim[2] "<<header->pixdim[2]<<endl;
   cout<<"header->pixdim[3] "<<header->pixdim[3]<<endl;
   cout<<"header->line1 "<<header->line1<<endl;
   cout<<"header->flag1 "<<header->flag1<<endl;
   cout<<"header->flag2 "<<header->flag2<<endl;
   cout<<"header->min "<<header->min<<endl;
   cout<<"header->max "<<header->max<<endl;
   cout<<"header->origin[0] "<<header->origin[0]<<endl;
   cout<<"header->origin[1] "<<header->origin[1]<<endl;
   cout<<"header->origin[2] "<<header->origin[2]<<endl;
   cout<<"header->origin[3] "<<header->origin[3]<<endl;
   cout<<"header->pixval_offset "<<header->pixval_offset<<endl;
   cout<<"header->pixval_cal "<<header->pixval_cal<<endl;
   cout<<"header->user_def1 "<<header->user_def1<<endl;
   cout<<"header->user_def2 "<<header->user_def2<<endl;
   cout<<"header->magic_number "<<header->magic_number<<endl;
   cout<<"image->field.datasize "<<image->field.datasize<<endl;
   cout<<"image->field.elementsize "<<image->field.elementsize<<endl;
*/

  // image data
  fread(image->field.data._void, image->field.elementsize, image->field.datasize, 
	image->file);

  // We "byte swap" the image data
  if (ByteSwapNeeded){
    short *newTypeS;
    int *newTypeI;
    float *newTypeF;
    double *newTypeD;

    // Depending on the image type, we cast the _void array to an array 
    // of image_type.
    // With the previous fread, all the data is stored in an array of void
    // (no type) but to do byte swapping we need this array of data to be 
    // an array of type image_type.
    // For instance, newTypeS is a pointer on short of _void. Using newTypeS,
    // we access the _void array as if it was an array of short.
    // NB : when we use swabShort on newTypeS we actually change the byte
    // order of _void because newTypeS is just a pointer to _void.

	int i =0;
    switch(header->image_type) {
    case IT_BINARY:
      break;
    case IT_CHAR:
      break;
    case IT_U_CHAR:
      break;
    case IT_SHORT:
    case IT_U_SHORT:
      newTypeS = static_cast<short *> (image->field.data._void);
      for (i = 0 ; i < image->field.datasize ; i ++ ){
	newTypeS[i]=binIO.swabShort(newTypeS[i]);
      }
      break;
    case IT_INT:
    case IT_U_INT:
    newTypeI = static_cast<int *> (image->field.data._void);
    for (i = 0 ; i < image->field.datasize ; i ++ ){
      newTypeI[i]=binIO.swabInt( newTypeI[i]);
    }
    break;
    case IT_FLOAT:
      newTypeF = static_cast<float *> (image->field.data._void);
      for (i = 0 ; i < image->field.datasize ; i ++ ){
	newTypeF[i]=binIO.swabFloat(newTypeF[i]);
      }
      break;
    case IT_DOUBLE:
      newTypeD = static_cast<double *> (image->field.data._void);
      for (i = 0 ; i < image->field.datasize ; i ++ ){
	newTypeD[i]=binIO.swabDouble(newTypeD[i]);
      }
      break;
    }
  }

  ipDeleteVectorint(dimsvec);
  ipDeleteVectorfloat(floatvec);

  ipCloseImage(image);


  return (gipl);
}

GIPL_HD * ipCreateGiplHeader(ipGiplImage *gipl)
{
  if (gipl->header)
    {
      free(gipl->header);
      gipl->header = NULL;
    }

  GIPL_HD *header = ipNewGiplHeader();

  int i = 0;
  for (i = 0; i < 4; i++) header->dim[i] = 1;
  for (i = 0; i < 4 && i < gipl->image->field.ndim; i++)
    header->dim[i] = (short) gipl->image->field.dims[i];

  switch(gipl->image->field.datatype) {
  case IP_BYTE:
    header->image_type = IT_CHAR;
    break;
  case IP_SHORT:
    header->image_type = IT_SHORT;
    break;
  case IP_INT:
    header->image_type = IT_INT;
    break;
  case IP_FLOAT:
    header->image_type = IT_FLOAT;
    break;
  case IP_DOUBLE:
    header->image_type = IT_DOUBLE;
    break;
  default:
#ifdef DEBUG
    fprintf(stderr, "Warning: ipCreateGiplHeader: unsuported data type: %d\n",
	    gipl->image->field.datatype);
#endif
    return NULL;
  }

  for (i = 0; i < gipl->image->field.ndim; i++)
    header->pixdim[i] = (double) gipl->image->pixdims[i];

  for (i = 0; i < gipl->image->field.ndim; i++)
    header->origin[i] = (double) gipl->image->origin[i];

  gipl->header = header;

  return header;
}


int ipWriteGipl(ipGiplImage *gipl, const int normalize)
{
  int            i;
  ipImage        *tmp_image = NULL;
  GIPL_HD        *header = gipl->header;
  ipImage        *image = gipl->image;

  if (!image) {
    fprintf(stderr,"Error: ipWriteGipl: no image");
    return -1;
  }

  if (!header)
    header = ipNewGiplHeader();

  ipOpenImage(image, "wb");
  if (!image->file) return(-1);

  for (i = 0; i < 4; i++) header->dim[i] = 1;
  for (i = 0; i < 4 && i < image->field.ndim; i++)
    header->dim[i] = (short) image->field.dims[i];

  fwrite(header->dim,(size_t) sizeof(header->dim), (size_t) 1, image->file);

  switch(image->field.datatype) {
  case IP_BYTE:
    header->image_type = IT_CHAR;
    break;
  case IP_SHORT:
    header->image_type = IT_SHORT;
    break;
  case IP_INT:
    header->image_type = IT_INT;
    break;
  case IP_FLOAT:
    header->image_type = IT_FLOAT;
    break;
  case IP_DOUBLE:
    header->image_type = IT_DOUBLE;
    break;
  default:
#ifdef DEBUG
    fprintf(stderr, "Warning: ipWriteGipl: unsuported data type: %d\n",
	    image->field.datatype);
#endif
    return(-1);
  }

  fwrite(&header->image_type,(size_t) sizeof(header->image_type), (size_t) 1,image->file);

  for (i = 0; i < image->field.ndim; i++)
    header->pixdim[i] = (double) image->pixdims[i];
  fwrite(header->pixdim,(size_t) sizeof(header->pixdim), (size_t) 1,image->file);
  fwrite(header->line1,(size_t) sizeof(header->line1), (size_t) 1,image->file);

  fwrite(header->matrix,(size_t) sizeof(header->matrix), (size_t) 1,image->file);
  fwrite(&header->flag1,(size_t) sizeof(header->flag1), (size_t) 1,image->file);
  fwrite(&header->flag2,(size_t) sizeof(header->flag2), (size_t) 1,image->file);
  // don't calculate explizitly maxima
  fwrite(&header->min,(size_t) sizeof(header->min), (size_t) 1,image->file);
  fwrite(&header->max,(size_t) sizeof(header->max), (size_t) 1,image->file);

  for (i = 0; i < image->field.ndim; i++)
    header->origin[i] = (double) image->origin[i];
  fwrite(header->origin,(size_t) sizeof(header->origin), (size_t) 1,image->file);

  fwrite(&header->pixval_offset,(size_t) sizeof(header->pixval_offset), (size_t) 1,
	 image->file);
  fwrite(&header->pixval_cal,(size_t) sizeof(header->pixval_cal), (size_t) 1,image->file);
  fwrite(&header->user_def1,(size_t) sizeof(header->user_def1), (size_t) 1,image->file);
  fwrite(&header->user_def2,(size_t) sizeof(header->user_def2), (size_t) 1,image->file);

  fwrite(&header->magic_number,(size_t) sizeof(header->magic_number), (size_t) 1,
	 image->file);
/*
   cout<<"header->pixdim[0] "<<header->pixdim[0]<<endl;
   cout<<"header->pixdim[1] "<<header->pixdim[1]<<endl;
   cout<<"header->pixdim[2] "<<header->pixdim[2]<<endl;
   cout<<"header->pixdim[3] "<<header->pixdim[3]<<endl;
   cout<<"header->line1 "<<header->line1<<endl;
   cout<<"header->flag1 "<<header->flag1<<endl;
   cout<<"header->flag2 "<<header->flag2<<endl;
   cout<<"header->min "<<header->min<<endl;
   cout<<"header->max "<<header->max<<endl;
   cout<<"header->origin[0] "<<header->origin[0]<<endl;
   cout<<"header->origin[1] "<<header->origin[1]<<endl;
   cout<<"header->origin[2] "<<header->origin[2]<<endl;
   cout<<"header->origin[3] "<<header->origin[3]<<endl;
   cout<<"header->pixval_offset "<<header->pixval_offset<<endl;
   cout<<"header->pixval_cal "<<header->pixval_cal<<endl;
   cout<<"header->user_def1 "<<header->user_def1<<endl;
   cout<<"header->user_def2 "<<header->user_def2<<endl;
   cout<<"header->magic_number "<<header->magic_number<<endl;
   cout<<"image->field.datasize "<<image->field.datasize<<endl;
   cout<<"image->field.elementsize "<<image->field.elementsize<<endl;
*/

  // image data
  if (normalize) {
    if (!(tmp_image = ipConvertImage(image, 1, image->field.datatype))) return(-1);
    fwrite(tmp_image->field.data._void, image->field.elementsize, 
	   image->field.datasize, image->file);
    ipDeleteImage(tmp_image);
  }
  else
    fwrite(image->field.data._void, image->field.elementsize, 
	   image->field.datasize, image->file);

  ipCloseImage(image);

  return(0);
}

int ipWriteGiplDirect(ipGiplImage *gipl)
{
  ipImage        *tmp_image = NULL;
  GIPL_HD        *header = gipl->header;
  ipImage        *image = gipl->image;

  if (!image) {
    fprintf(stderr,"Error: ipWriteGiplDirect: no image");
    return -1;
  }

  if (!header) {
    fprintf(stderr,"Error: ipWriteGiplDirect: no header");
    return -1;
  }


  ipOpenImage(image, "wb");
  if (!image->file) return(-1);

  fwrite(header->dim,(size_t) sizeof(header->dim), (size_t) 1, image->file);

  fwrite(&header->image_type,(size_t) sizeof(header->image_type), (size_t) 1,image->file);

  fwrite(header->pixdim,(size_t) sizeof(header->pixdim), (size_t) 1,image->file);
  fwrite(header->line1,(size_t) sizeof(header->line1), (size_t) 1,image->file);

  fwrite(header->matrix,(size_t) sizeof(header->matrix), (size_t) 1,image->file);
  fwrite(&header->flag1,(size_t) sizeof(header->flag1), (size_t) 1,image->file);
  fwrite(&header->flag2,(size_t) sizeof(header->flag2), (size_t) 1,image->file);
  // don't calculate explizitly maxima
  fwrite(&header->min,(size_t) sizeof(header->min), (size_t) 1,image->file);
  fwrite(&header->max,(size_t) sizeof(header->max), (size_t) 1,image->file);

  fwrite(header->origin,(size_t) sizeof(header->origin), (size_t) 1,image->file);

  fwrite(&header->pixval_offset,(size_t) sizeof(header->pixval_offset), (size_t) 1,
	 image->file);
  fwrite(&header->pixval_cal,(size_t) sizeof(header->pixval_cal), (size_t) 1,image->file);
  fwrite(&header->user_def1,(size_t) sizeof(header->user_def1), (size_t) 1,image->file);
  fwrite(&header->user_def2,(size_t) sizeof(header->user_def2), (size_t) 1,image->file);

  fwrite(&header->magic_number,(size_t) sizeof(header->magic_number), (size_t) 1,
	 image->file);

  // image data
  fwrite(image->field.data._void, image->field.elementsize, 
	 image->field.datasize, image->file);

  ipCloseImage(image);

  return(0);
}


