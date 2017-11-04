/******************************************************
 general declarations of routines in matrix.c
 ******************************************************/

#ifndef __IP_MATRIX_H__
#define __IP_MATRIX_H__

#include <stdlib.h>
#include <stdio.h>
#include <errno.h>
#include <string.h>
#include "macros.h"
#include "misc.h"

typedef struct {
  int           ndim,
                dims[1],
                pos[1],
                datasize;
  ipDataType  datatype;
  size_t        elementsize;
  ipDataUnion data;
} ipVector;

typedef struct {
  int           ndim,
                dims[1],
                pos[1],
                datasize;
  ipDataType  datatype;
  size_t        elementsize;
  int*          data;
} ipVectorint;

typedef struct {
  int           ndim,
                dims[1],
                pos[1],
                datasize;
  ipDataType  datatype;
  size_t        elementsize;
  float*        data;
} ipVectorfloat;

typedef struct {
  int           ndim,
                dims[2],
                pos[2],
                datasize;
  ipDataType  datatype;
  size_t        elementsize;
  ipDataUnion data;
} ipMatrix;

typedef struct {
  int           ndim,
                dims[2],
                pos[2],
                datasize;
  ipDataType  datatype;
  size_t        elementsize;
  float*        data;
} ipMatrixfloat;

typedef struct {
  int           ndim,           /* nr of dimensions */                   
                *dims,		/* dims[ndim], size in each dimension */ 
                *pos,		/* act position [ndim] in data for iterator */
                datasize;	/* product of all dims */
  ipDataType  datatype;	/* image data type */                        
  size_t        elementsize;	/* nr of bytes per element */ 
  ipDataUnion data;		/* pointer to image data */                  
} ipArray;

typedef struct {
  int           ndim,           /* nr of dimensions */                   
                *dims,		/* dims[ndim], size in each dimension */ 
                *pos,		/* act position [ndim] in data for iterator */
                datasize;	/* product of all dims */
  ipDataType  datatype;	/* image data type */                        
  size_t        elementsize;	/* nr of bytes per element */ 
  unsigned char *data;
} ipArraybyte;

typedef struct {
  int           ndim,           /* nr of dimensions */                   
                *dims,		/* dims[ndim], size in each dimension */ 
                *pos,		/* act position [ndim] in data for iterator */
                datasize;	/* product of all dims */
  ipDataType  datatype;	/* image data type */                        
  size_t        elementsize;	/* nr of bytes per element */ 
  float         *data;
} ipArrayfloat;


/* function prototypes */
int ipAllocateArrayData(ipArray *array, const ipVectorint *dimsvec,
			  const ipDataType type);

ipArray* ipNewArray();
ipVector* ipNewVector(const int cols, const ipDataType type);
ipMatrix* ipNewMatrix(const int cols, const int rows, const ipDataType type);

// initializes data; sizes must match! (no bounds check)
void ipSetVectorData(ipVector *vec, const void *data);
void ipSetMatrixData(ipMatrix *mat, const void *data);
void ipSetArrayData(ipArray *arr, const void *data);

// direct access to data members
void* ipGetVectorData(const ipVector *vec);
void* ipGetMatrixData(const ipMatrix *mat);
void* ipGetArrayData(const ipArray *arr);

int ipCopyArrayData(ipArray *dest, const ipArray *src);
ipArray* ipDeepCopyArray(const ipArray *arr);

ipVector* ipCopyVector(const ipVector *vec);
ipMatrix* ipCopyMatrix(const ipMatrix *mat);

int ipArrayDatatoByte(ipArray *arr, const int norm);
int ipArrayDatatoShort(ipArray *arr, const int norm);
int ipArrayDatatoFloat(ipArray *arr, const int norm);
int ipConvertArrayData(ipArray *arr, const int norm, const ipDataType type);
ipArray* ipConvertArray(const ipArray *src, const int norm, const ipDataType type);

void ipDeleteArrayData(ipArray *arr);
void ipDeleteArray(ipArray *arr);
void ipDeleteVector(ipVector *vec);
void ipDeleteMatrix(ipMatrix *mat);

#define ipNewVectorfloat(cols) ((ipVectorfloat *) ipNewVector(cols, IP_FLOAT))
#define ipNewVectorint(cols) ((ipVectorint *) ipNewVector(cols, IP_INT))
#define ipNewMatrixfloat(cols, rows) ((ipMatrixfloat *) ipNewMatrix(cols, rows, IP_FLOAT))
#define ipCopyVectorfloat(vec) ((ipVectorfloat *) ipCopyVector((ipVector *) vec))
#define ipCopyMatrixfloat(mat) ((ipMatrixfloat *) ipCopyMatrix((ipMatrix *) mat))

#define ipDeleteVectorfloat(vec) (ipDeleteVector((ipVector *) vec))
#define ipDeleteVectorint(vec) (ipDeleteVector((ipVector *) vec))
#define ipDeleteMatrixfloat(mat) (ipDeleteMatrix((ipMatrix *) mat))


#endif

