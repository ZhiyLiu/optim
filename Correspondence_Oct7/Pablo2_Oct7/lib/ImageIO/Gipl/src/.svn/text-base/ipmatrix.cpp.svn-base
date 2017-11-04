/*
 * matrix.cc
 *
 * author:    Manuel Sturm
 * created:   24 Jan 1997
 *
 * some matrix basics
 * adapted by Martin Andreas Styner 03 Jul 97
 */

#include "ipmatrix.h"


/*------------------------------------------------------------------------------*/
/* Internally used utility functions */
/*------------------------------------------------------------------------------*/

int ipAllocateArrayData(ipArray *array, const ipVectorint *dimsvec,
			const ipDataType type)
{
  int i;

  array->ndim = dimsvec->dims[0];
  array->datatype = type;

  if (!(array->dims = (int *) realloc(array->dims, array->ndim * sizeof(int)))) {
    fprintf(stderr,"Error: ipAllocateArrayData [%s, line %d]: memory allocation failed:",
	    __FILE__, __LINE__);
    perror("");
    exit(errno); 
  }

  memcpy((void *) array->dims, (void *) dimsvec->data, array->ndim * sizeof(int));

  if (!(array->pos = (int *) realloc(array->pos, array->ndim * sizeof(int)))) {
    fprintf(stderr,"Error: ipAllocateArrayData [%s, line %d]: memory allocation failed:",
	    __FILE__, __LINE__);
    perror("");
    exit(errno);
  }

  memset((void *) array->pos, 0, array->ndim * sizeof(int));

  for (i = 1, array->datasize = array->dims[0]; i < array->ndim; i++)
    array->datasize *= array->dims[i];

  if ((array->elementsize = ipGetDataSize(type)) < 0) {
    free(array->dims);
    free(array->pos);

    return(array->elementsize);
  }

  array->data._void = ipAllocateData(array->datasize, array->elementsize);

  return(0);
}


ipArray* ipNewArray()
{
  ipArray* array = NULL;

  if (!(array = (ipArray *) malloc(sizeof(ipArray)))) {
    fprintf(stderr,"Error: ipNewArray [%s, line %d]: memory allocation failed:",
	    __FILE__, __LINE__);
    perror("");
    exit(errno);
  }

  memset(array, 0, sizeof(ipArray));

  return(array);

}


ipVector* ipNewVector(const int cols, const ipDataType type)
{
  ipVector* vec = NULL;

  if (!(vec = (ipVector *) malloc(sizeof(ipVector)))) {
    fprintf(stderr,"Error: ipNewVector [%s, line %d]: memory allocation failed:",
	    __FILE__, __LINE__);
    perror("");
    exit(errno);
  }

  memset(vec, 0, sizeof(ipVector));

  if ((vec->elementsize = ipGetDataSize(type)) < 0) {
    ipDeleteVector(vec);

    return(NULL);
  }

  vec->data._void = ipAllocateData(cols, vec->elementsize);

  vec->ndim = 1;
  vec->dims[0] = vec->datasize = cols;
  vec->pos[0] = 0;
  vec->datatype = type;

  return(vec);
}


ipMatrix* ipNewMatrix(const int cols, const int rows,
		      const ipDataType type)
{
  ipMatrix* mat = NULL;


  if (!(mat = (ipMatrix *) malloc(sizeof(ipMatrix)))) {
    fprintf(stderr,"Error: ipNewMatrix [%s, line %d]: memory allocation failed:",
	    __FILE__, __LINE__);
    perror("");
    exit(errno);
  }

  memset(mat, 0, sizeof(ipMatrix));

  if ((mat->elementsize = ipGetDataSize(type)) < 0) {
    ipDeleteMatrix(mat);

    return(NULL);
  }

  mat->data._void = ipAllocateData(cols * rows, mat->elementsize);

  mat->ndim = 2;
  mat->dims[0] = cols;
  mat->dims[1] = rows;
  mat->datasize = cols * rows;
  mat->pos[0] = mat->pos[1] = 0;
  mat->datatype = type;

  return(mat);
}

/*------------------------------------------------------------------------------*/
/* Utility functions for direct data access                                     */
/*------------------------------------------------------------------------------*/

void ipSetVectorData(ipVector *vec, const void *data)
{
  if (vec->data._void)
    memcpy(vec->data._void, data, vec->dims[0] * vec->elementsize);
  else {
    fprintf(stderr,"Error: ipSetVectorData [%s, line %d]: write attempt to unallocated memory!\n",
	    __FILE__, __LINE__);
    exit(-1);
  }
}


void ipSetMatrixData(ipMatrix *mat, const void *data)
{
  if (mat->data._void)
    memcpy(mat->data._void, data, mat->datasize * mat->elementsize);
  else {
    fprintf(stderr,"Error: ipSetMatrixData [%s, line %d]: write attempt to unallocated memory!\n",
	    __FILE__, __LINE__);
    exit(-1);
  }
}


void ipSetArrayData(ipArray *arr, const void *data)
{
  if (arr->data._void)
    memcpy(arr->data._void, data, arr->datasize * arr->elementsize);
  else {
    fprintf(stderr,"Error: ipSetArrayData [%s, line %d]: write attempt to unallocated memory!\n",
	    __FILE__, __LINE__);
    exit(-1);
  }
}


void* ipGetVectorData(const ipVector *vec)
{
  void *tmpp = NULL;

  if (!(tmpp = malloc(vec->dims[0] * vec->elementsize))) {
    fprintf(stderr,"Error: ipGetVectorData [%s, line %d]: memory allocation failed:",
	    __FILE__, __LINE__);
    perror("");
    exit(errno);
  }

  memcpy(tmpp, vec->data._void, vec->dims[0] * vec->elementsize);

  return(tmpp);
}


void* ipGetMatrixData(const ipMatrix *mat)
{
  void *tmpp = NULL;

  if (!(tmpp = malloc(mat->datasize * mat->elementsize))) {
    fprintf(stderr,"Error: ipGetMatrixData [%s, line %d]: memory allocation failed:",
	    __FILE__, __LINE__);
    perror("");
    exit(errno);
  }

  memcpy(tmpp, mat->data._void, mat->datasize * mat->elementsize);

  return(tmpp);
}


void* ipGetArrayData(const ipArray *arr)
{
  void *tmpp = NULL;

  if (!(tmpp = malloc(arr->datasize * arr->elementsize))) {
    fprintf(stderr,"Error: ipGetArrayData [%s, line %d]: memory allocation failed:",
	    __FILE__, __LINE__);
    perror("");
    exit(errno);
  }

  memcpy(tmpp, arr->data._void, arr->datasize * arr->elementsize);

  return(tmpp);
}



/*------------------------------------------------------------------------------*/
/* Deep Copy Constructors                                                       */
/*------------------------------------------------------------------------------*/

int ipCopyArrayData(ipArray *dest, const ipArray *src)
{
  ipVectorint tmp_vec;

  tmp_vec.dims[0] = src->ndim; 
  tmp_vec.data = src->dims;

  if (ipAllocateArrayData(dest, &tmp_vec, src->datatype))
    return(-1);

  memcpy(dest->data._void, src->data._void, src->datasize * src->elementsize);

  return(0);
}

ipArray* ipDeepCopyArray(const ipArray *src)
{
  ipArray *dest = ipNewArray();

  if (ipCopyArrayData(dest, src)) {
    ipDeleteArray(dest);

    return(NULL);
  }

  return(dest);
}

ipVector* ipCopyVector(const ipVector *src)
{
  ipVector* dest = ipNewVector(src->dims[0], src->datatype);

  if (dest)
    memcpy(dest->data._void, src->data._void, src->dims[0] * src->elementsize);

  return(dest);
}

ipMatrix* ipCopyMatrix(const ipMatrix *src)
{
  ipMatrix* dest = ipNewMatrix(src->dims[0], src->dims[1], src->datatype);

  if (dest)
    memcpy(dest->data._void, src->data._void, src->dims[0] * src->dims[1] * 
	   src->elementsize);

  return(dest);
}


/*------------------------------------------------------------------------------*/
/* Converters                                                                   */
/*------------------------------------------------------------------------------*/

/* converts data of array to IP_BYTE
 * 
 * if norm = 1, data is scaled to the range of 0..255
 */

int ipArrayDatatoByte(ipArray *arr, const int norm)
{
  int            i;
  double         dmax = 1.0e-38, dmin = 1.0e38, scale = 1.0;
  unsigned char *ucptr = NULL, *tmp_data = NULL;

  /* data is scaled to 0..255 if norm=1 or minvalue < 0 or maxvalue > 255 */

  ipDataUnion  dptr;

  if ((arr->datatype ==  IP_BYTE) && !norm) {
#ifdef DEBUG
    fprintf(stderr, "Info: ipArrayDatatoByte: use ipCopyArrayData() instead\n");
#endif
    return(-1);
  }

  tmp_data = (unsigned char*) ipAllocateData(arr->datasize, ipGetDataSize(IP_BYTE));

  switch (arr->datatype) {
  case IP_BYTE:
    for (i = 0, ucptr = arr->data._byte; i < arr->datasize; 
	 i++, ucptr++){
      dmax = MAX(dmax, (double) *ucptr);
      dmin = MIN(dmin, (double) *ucptr);
    }

    if ((dmax - dmin) > 0.0) {
      scale = 255.0 / (dmax - dmin);
      for (i = 0, dptr._byte = arr->data._byte, 
	     ucptr = tmp_data; 
	   i < arr->datasize; i++, (dptr._byte)++, ucptr++)
	*ucptr =(unsigned char) (((double) *(dptr._byte) - dmin) * scale);
    }
    break;
  case IP_SHORT:
    for (i = 0, dptr._short = arr->data._short; i < arr->datasize; 
	 i++, (dptr._short)++){
      dmax = MAX(dmax, (double) *(dptr._short));
      dmin = MIN(dmin, (double) *(dptr._short));
    }

    if(!norm && (dmin >= 0.0) && (dmax <= 255.0))
      for (i = 0, dptr._short = arr->data._short,
	     ucptr = tmp_data; 
	   i < arr->datasize; i++, (dptr._short)++, ucptr++)
	*ucptr = (unsigned char) *(dptr._short);
    else 
      if ((dmax - dmin) > 0.0) {
	scale = 255.0 / (dmax - dmin);
	for (i = 0, dptr._short = arr->data._short,
	       ucptr = tmp_data;
	     i < arr->datasize; i++, (dptr._short)++, ucptr++)
	  *ucptr = (unsigned char) (((double) *(dptr._short) - dmin) * scale);
      }
    break;
  case IP_INT:
    for (i = 0, dptr._int = arr->data._int; i < arr->datasize; 
	 i++, (dptr._int)++){
      dmax = MAX(dmax, (double) *(dptr._int));
      dmin = MIN(dmin, (double) *(dptr._int));
    }

    if(!norm && (dmin >= 0.0) && (dmax <= 255.0))
      for (i = 0, dptr._int = arr->data._int,
	     ucptr = tmp_data; 
	   i < arr->datasize; i++, (dptr._int)++, ucptr++)
	*ucptr = (unsigned char) *(dptr._int);
    else {
      if ((dmax - dmin) > 0.0)
	scale = 255.0 / (dmax - dmin);
      for (i = 0, dptr._int = arr->data._int,
	     ucptr = tmp_data;
	   i < arr->datasize; i++, (dptr._int)++, ucptr++)
	*ucptr = (unsigned char) (((double) *(dptr._int) - dmin) * scale);
    }
    break;
  case IP_FLOAT:
    for (i = 0, dptr._float = arr->data._float; i < arr->datasize; 
	 i++, (dptr._float)++){
      dmax = MAX(dmax, (double) *(dptr._float));
      dmin = MIN(dmin, (double) *(dptr._float));
    }

    if(!norm && (dmin >= 0.0) && (dmax <= 255.0))
      for (i = 0, dptr._float = arr->data._float,
	     ucptr = tmp_data; 
	   i < arr->datasize; i++, (dptr._float)++, ucptr++)
	*ucptr = (unsigned char) *(dptr._float);
    else 
      if ((dmax - dmin) > 0.0) {
	scale = 255.0 / (dmax - dmin);
	for (i = 0, dptr._float = arr->data._float,
	       ucptr = tmp_data;
	     i < arr->datasize; i++, (dptr._float)++, ucptr++)
	  *ucptr = (unsigned char) (((double) *(dptr._float) - dmin) * scale);
      }
    break;
  case IP_DOUBLE:
    for (i = 0, dptr._double = arr->data._double; i < arr->datasize; 
	 i++, (dptr._double)++){
      dmax = MAX(dmax, *(dptr._double));
      dmin = MIN(dmin, *(dptr._double));
    }

    if(!norm && (dmin >= 0.0) && (dmax <= 255.0))
      for (i = 0, dptr._double = arr->data._double,
	     ucptr = tmp_data; 
	   i < arr->datasize; i++, (dptr._double)++, ucptr++)
	*ucptr = (unsigned char) *(dptr._double);
    else
      if ((dmax - dmin) > 0.0) {
	scale = 255.0 / (dmax - dmin);
	for (i = 0, dptr._double = arr->data._double,
	       ucptr = tmp_data;
	     i < arr->datasize; i++, (dptr._double)++, ucptr++)
	  *ucptr = (unsigned char) ((*(dptr._double) - dmin) * scale);
      }
    break;
  default:
#ifdef DEBUG
    fprintf(stderr, "Warning: ipArraytoByte: conversion from datatype [%d] not yet supported\n",
	    arr->datatype);
#endif

    ipDeleteArrayData(arr);
    free(tmp_data);

    return(-1);
  }

  free(arr->data._void);

  arr->elementsize = sizeof(unsigned char);
  arr->datatype = IP_BYTE;
  arr->data._byte = tmp_data;

  return(0);
}

/* converts data of array to IP_FLOAT */

int ipArrayDatatoFloat(ipArray *arr, const int norm)
{
  int            i;
  double         dmax = 1.0e-38, dmin = 1.0e38, scale = 1.0;
  float         *ucptr = NULL, *tmp_data = NULL;

  ipDataUnion  dptr;

  /* scaling to 0..255 is performed when norm = 1 */

  if ((arr->datatype ==  IP_FLOAT) && !norm) {
#ifdef DEBUG
    fprintf(stderr, "Info: ipArrayDatatoFloat: use ipCopyArrayData() instead\n");
#endif
    return(-1);
  }

  tmp_data = (float *) ipAllocateData(arr->datasize, ipGetDataSize(IP_FLOAT));

  switch (arr->datatype) {
  case IP_BYTE:
    for (i = 0, dptr._byte = arr->data._byte; i < arr->datasize; 
	 i++, (dptr._byte)++){
      dmax = MAX(dmax, (double) *(dptr._byte));
      dmin = MIN(dmin, (double) *(dptr._byte));
    }

    if(!norm)
      for (i = 0, dptr._byte = arr->data._byte,
	     ucptr = tmp_data; 
	   i < arr->datasize; i++, (dptr._byte)++, ucptr++)
	*ucptr = (float) *(dptr._byte);
    else
      if ((dmax - dmin) > 0.0) {
	scale = 255.0 / (dmax - dmin);
	for (i = 0, dptr._byte = arr->data._byte,
	       ucptr = tmp_data;
	     i < arr->datasize; i++, (dptr._byte)++, ucptr++)
	  *ucptr = (float) (((double) *(dptr._byte) - dmin) * scale);
      }
    break;
  case IP_SHORT:
    for (i = 0, dptr._short = arr->data._short; i < arr->datasize; 
	 i++, (dptr._short)++){
      dmax = MAX(dmax, (double) *(dptr._short));
      dmin = MIN(dmin, (double) *(dptr._short));
    }

    if(!norm)
      for (i = 0, dptr._short = arr->data._short,
	     ucptr = tmp_data; 
	   i < arr->datasize; i++, (dptr._short)++, ucptr++)
	*ucptr = (float) *(dptr._short);
    else
      if ((dmax - dmin) > 0.0) {
	scale = 255.0 / (dmax - dmin);
	for (i = 0, dptr._short = arr->data._short,
	       ucptr = tmp_data;
	     i < arr->datasize; i++, (dptr._short)++, ucptr++)
	  *ucptr = (float) (((double) *(dptr._short) - dmin) * scale);
      }
    break;
  case IP_INT:
    for (i = 0, dptr._int = arr->data._int; i < arr->datasize; 
	 i++, (dptr._int)++){
      dmax = MAX(dmax, (double) *(dptr._int));
      dmin = MIN(dmin, (double) *(dptr._int));
    }

    if(!norm)
      for (i = 0, dptr._int = arr->data._int,
	     ucptr = tmp_data; 
	   i < arr->datasize; i++, (dptr._int)++, ucptr++)
	*ucptr = (float) *(dptr._int);
    else
      if ((dmax - dmin) > 0.0) {
	scale = 255.0 / (dmax - dmin);
	for (i = 0, dptr._int = arr->data._int,
	       ucptr = tmp_data;
	     i < arr->datasize; i++, (dptr._int)++, ucptr++)
	  *ucptr = (float) (((double) *(dptr._int) - dmin) * scale);
      }
    break;
  case IP_FLOAT:
    for (i = 0, dptr._float = arr->data._float; i < arr->datasize; 
	 i++, (dptr._float)++){
      dmax = MAX(dmax, (double) *(dptr._float));
      dmin = MIN(dmin, (double) *(dptr._float));
    }

    if ((dmax - dmin) > 0.0) {
      scale = 255.0 / (dmax - dmin);
      for (i = 0, dptr._float = arr->data._float,
	     ucptr = tmp_data;
	   i < arr->datasize; i++, (dptr._float)++, ucptr++)
	*ucptr = (float) (((double) *(dptr._float) - dmin) * scale);
    }
    break;
  case IP_DOUBLE:
    for (i = 0, dptr._double = arr->data._double; i < arr->datasize; 
	 i++, (dptr._double)++){
      dmax = MAX(dmax, *(dptr._double));
      dmin = MIN(dmin, *(dptr._double));
    }

    if(!norm)
      for (i = 0, dptr._double = arr->data._double,
	     ucptr = tmp_data; 
	   i < arr->datasize; i++, (dptr._double)++, ucptr++)
	*ucptr = (float) *(dptr._double);
    else
      if ((dmax - dmin) > 0.0) {
	scale = 255.0 / (dmax - dmin);
	for (i = 0, dptr._double = arr->data._double,
	       ucptr = tmp_data;
	     i < arr->datasize; i++, (dptr._double)++, ucptr++)
	  *ucptr = (float) ((*(dptr._double) - dmin) * scale);
      }
    break;
  default:
#ifdef DEBUG
    fprintf(stderr, "Warning: ipArraytoFloat: conversion from datatype [%d] not yet supported\n",
	    arr->datatype);
#endif
    ipDeleteArrayData(arr);
    free(tmp_data);

    return(-1);
  }

  free(arr->data._void);

  arr->elementsize = sizeof(float);
  arr->datatype = IP_FLOAT;
  arr->data._float = tmp_data;

  return(0);
}

// Convert Data to Short

int ipArrayDatatoShort(ipArray *arr, const int norm)
{
  int            i;
  double         dmax = 0, dmin = 255, scale = 1.0;
  short         *ucptr = NULL, *tmp_data = NULL;

  ipDataUnion  dptr;

  /* scaling to 0..255 is performed when norm = 1 */

  if ((arr->datatype ==  IP_SHORT) && !norm) {
#ifdef DEBUG
    fprintf(stderr, "Info: ipArrayDatatoShort: use ipCopyArrayData() instead\n");
#endif
    return(-1);
  }

  tmp_data = (short *) ipAllocateData(arr->datasize, ipGetDataSize(IP_SHORT));

  switch (arr->datatype) {
  case IP_BYTE:
    for (i = 0, dptr._byte = arr->data._byte; i < arr->datasize; 
	 i++, (dptr._byte)++){
      dmax = MAX(dmax, (double) *(dptr._byte));
      dmin = MIN(dmin, (double) *(dptr._byte));
    }

    if(!norm)
      for (i = 0, dptr._byte = arr->data._byte,
	     ucptr = tmp_data; 
	   i < arr->datasize; i++, (dptr._byte)++, ucptr++)
	*ucptr = (short) *(dptr._byte);
    else
      if ((dmax - dmin) > 0.0) {
	scale = 255.0 / (dmax - dmin);
	for (i = 0, dptr._byte = arr->data._byte,
	       ucptr = tmp_data;
	     i < arr->datasize; i++, (dptr._byte)++, ucptr++)
	  *ucptr = (short) (((double) *(dptr._byte) - dmin) * scale);
      }
    break;
  case IP_SHORT:
    for (i = 0, dptr._short = arr->data._short; i < arr->datasize; 
	 i++, (dptr._short)++){
      dmax = MAX(dmax, (double) *(dptr._short));
      dmin = MIN(dmin, (double) *(dptr._short));
    }

    if ((dmax - dmin) > 0.0) {
      scale = 255.0 / (dmax - dmin);
      for (i = 0, dptr._short = arr->data._short,
	     ucptr = tmp_data;
	   i < arr->datasize; i++, (dptr._short)++, ucptr++)
	*ucptr = (short) (((double) *(dptr._short) - dmin) * scale);
    }
    break;
  case IP_INT:
    for (i = 0, dptr._int = arr->data._int; i < arr->datasize; 
	 i++, (dptr._int)++){
      dmax = MAX(dmax, (double) *(dptr._int));
      dmin = MIN(dmin, (double) *(dptr._int));
    }

    if(!norm)
      for (i = 0, dptr._int = arr->data._int,
	     ucptr = tmp_data; 
	   i < arr->datasize; i++, (dptr._int)++, ucptr++)
	*ucptr = (short) *(dptr._int);
    else
      if ((dmax - dmin) > 0.0) {
	scale = 255.0 / (dmax - dmin);
	for (i = 0, dptr._int = arr->data._int,
	       ucptr = tmp_data;
	     i < arr->datasize; i++, (dptr._int)++, ucptr++)
	  *ucptr = (short) (((double) *(dptr._int) - dmin) * scale);
      }
    break;
  case IP_FLOAT:
    for (i = 0, dptr._float = arr->data._float; i < arr->datasize; 
	 i++, (dptr._float)++){
      dmax = MAX(dmax, (double) *(dptr._float));
      dmin = MIN(dmin, (double) *(dptr._float));
    }

    if(!norm)
      for (i = 0, dptr._float = arr->data._float,
	     ucptr = tmp_data; 
	   i < arr->datasize; i++, (dptr._float)++, ucptr++)
	*ucptr = (short) *(dptr._float);
    else 
      if ((dmax - dmin) > 0.0) {
	scale = 255.0 / (dmax - dmin);
	for (i = 0, dptr._float = arr->data._float,
	       ucptr = tmp_data;
	     i < arr->datasize; i++, (dptr._float)++, ucptr++)
	  *ucptr = (short) (((double) *(dptr._float) - dmin) * scale);
      }
    break;
  case IP_DOUBLE:
    for (i = 0, dptr._double = arr->data._double; i < arr->datasize; 
	 i++, (dptr._double)++){
      dmax = MAX(dmax, *(dptr._double));
      dmin = MIN(dmin, *(dptr._double));
    }

    if(!norm)
      for (i = 0, dptr._double = arr->data._double,
	     ucptr = tmp_data; 
	   i < arr->datasize; i++, (dptr._double)++, ucptr++)
	*ucptr = (short) *(dptr._double);
    else
      if ((dmax - dmin) > 0.0) {
	scale = 255.0 / (dmax - dmin);
	for (i = 0, dptr._double = arr->data._double,
	       ucptr = tmp_data;
	     i < arr->datasize; i++, (dptr._double)++, ucptr++)
	  *ucptr = (short) ((*(dptr._double) - dmin) * scale);
      }
    break;
  default:
#ifdef DEBUG
    fprintf(stderr, "Warning: ipArraytoShort: conversion from datatype [%d] not yet supported\n",
	    arr->datatype);
#endif
    ipDeleteArrayData(arr);
    free(tmp_data);

    return(-1);
  }

  free(arr->data._void);

  arr->elementsize = sizeof(short);
  arr->datatype = IP_SHORT;
  arr->data._short = tmp_data;

  return(0);
}

int ipConvertArrayData(ipArray *arr, const int norm, const ipDataType type)
{

  if ((arr->datatype ==  type) && !norm)
    return(0);

  switch (type) {
  case IP_BYTE:
    if (ipArrayDatatoByte(arr, norm)) {
      ipDeleteArrayData(arr);

      return(-1);
    }
    break;
  case IP_SHORT:
    if (ipArrayDatatoShort(arr, norm)) {
      ipDeleteArrayData(arr);

      return(-1);
    }
    break;
  case IP_FLOAT:
    if (ipArrayDatatoFloat(arr, norm)) {
      ipDeleteArrayData(arr);

      return(-1);
    }
    break;
  default:
#ifdef DEBUG
    fprintf(stderr, "Warning: ipConvertArrayData: conversion to datatype [%d] not yet supported\n",
	    type);
#endif

    return(-1);
  }

  return(0);
}

/* these procedure creates a new array and then calls the conversion routine */

ipArray* ipConvertArray(const ipArray *src, const int norm, const ipDataType type)
{
  ipArray *dest = ipDeepCopyArray(src);

  if (ipConvertArrayData(dest, norm, type)) {
    ipDeleteArray(dest);

    return(NULL);
  }

  return(dest);
}


/*------------------------------------------------------------------------------*/
/* Destructors                                                                  */
/*------------------------------------------------------------------------------*/

void ipDeleteArrayData(ipArray *array)
{
  if (array) {
    free(array->dims);
    free(array->pos);
    free(array->data._void);

    array->dims = array->pos = NULL;
    array->data._void = NULL;
  }
}

void ipDeleteArray(ipArray *array)
{
  if (array) {
    ipDeleteArrayData(array);

    free(array);
  }
}

void ipDeleteVector(ipVector *vec)
{
  if (vec) {
    free(vec->data._void);
    free(vec);
  }
}

void ipDeleteMatrix(ipMatrix *mat)
{
  if (mat) {
    free(mat->data._void);
    free(mat);
  }
}


