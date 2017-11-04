/*
 * giplio.h
 * (Guys Hospital London, Image Processing Library File format)
 *
 * author:  msturm
 * created: implemented as avs read module, 05 Nov 1996 
 * changes: included in ip lib, 18 Feb 1997 
 *
 */

/*
 * giplio.h
 * 
 * author:  msturm 
 * created: 14 Feb 1997
 * changes: 
 */

#ifndef __IP_GIPLIO_H__
#define __IP_GIPLIO_H__

#include <stdio.h>
#include <string.h>
#include "macros.h"
#include "ipimage.h"
#include "gipl_header.h"

typedef struct
{ 
   ipImage *image;
   GIPL_HD *header;

} ipGiplImage;

// Constructor, members are initialized to zero, call their constructors
// if necessary
ipGiplImage* ipNewGipl();
// Destructor 
void ipDeleteGipl(ipGiplImage *in);

GIPL_HD * ipCreateGiplHeader(ipGiplImage *gipl);
// Read gipl file into a newly created image-struct
ipGiplImage* ipReadGipl(const char *filename);

// if struct-member 'header' of 'gipl' is == NULL then
// a default header will be written out
int ipWriteGipl(ipGiplImage *gipl, const int normalize);
int ipWriteGiplDirect(ipGiplImage *gipl);

#endif

