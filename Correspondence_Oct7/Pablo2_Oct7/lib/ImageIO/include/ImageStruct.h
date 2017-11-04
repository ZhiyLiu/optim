#ifndef IMAGE_STRUCT_H
#define IMAGE_STRUCT_H


// This structure is used for conversion between image files
// other than .raw3 and pablo's internal variables.
struct ImageStruct
{
	int dims[3];
	unsigned long len;	// dims[0]*dims[1]*dims[2]; not always set
	float spacing[3];                                       // AGG: use double?
	float origin[3];

	// The intensity range should never exceed the range of
	// either a short or unsigned short.  It is stored here
	// as an int, only because no 17-bit quantity is available.
	int min;		// Not always set
	int max;		// Not always set

	// This pointer may be to data stored as either short or
	// unsigned short.  If headerOnly is true when loading, then
	// this will be NULL.
	unsigned short * voxels;
	bool dataIsShort;	// True if voxels contains shorts.

	// File Image3D.h defines enum modality_t.  This enumeration may
	// be used for all Image I/O, including .raw3, so it is important
	// that the following integer have a value that is identical to
	// the enumeration.
	int modality;	// Must be set to correspond to enum modality_t
};



#endif

