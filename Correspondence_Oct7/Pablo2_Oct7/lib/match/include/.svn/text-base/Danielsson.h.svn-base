#ifndef DANIELSSON_H
#define DANIELSSON_H

/*  Danielsson Distance Map Algorithm
*/

#ifdef BINARY

int edt3ddan(char * bdryIm, int xdim, int ydim, int zdim,
			 double spacingx, double spacingy, double spacingz,
			 unsigned char metric,	// 6 or 26: # of neighbor distances to check
			 short ** odx,	// output distance maps, per direction
			 short ** ody,
			 short ** odz);
#else
int edt3ddan(char * bdryIm, int xdim, int ydim, int zdim,
			 unsigned char metric,
			 short ** odx, short ** ody, short ** odz);
#endif


#endif

