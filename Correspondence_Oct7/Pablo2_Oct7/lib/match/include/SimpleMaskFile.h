#ifndef SIMPLE_MASK_FILE_H
#define SIMPLE_MASK_FILE_H

/*
	Joshua Stough, 3/04
	This is a very simplified version of MaskFile, which I didn't 
	write and don't completely understand.  This is merely to use
	the target image of a training case for the testing on a new 
	image.  So a user first runs a getprofiles program with their
	model and image it is fit to.  The result is read by this 
	and the profiles replace the mask profiles in mask.cpp.
*/

class SimpleMaskFile
{
public:
    SimpleMaskFile() {}

    void readBinary(const char * filename);
	void read(const char * filename);

	bool getValue(int figureId, int point_ind, int sample, double & value);
	double getTau(int figureId, int which);
	int getNumSamples(int figureId);
	int getNumPoints(int figureId) { return numPoints[figureId]; }

private:
    int numFigures;
	double ** values;

	int * numPoints;
	int * numSamples;
	char ** figNames;
	int * figIds;
	double *taus;

};

#endif

