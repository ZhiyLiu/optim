#ifndef _PCA_H_
#define _PCA_H_

#include <vector>

#ifdef PCA_VISUALIZER

class BYU;

class PCA
{

public:

	PCA();
	~PCA();
	void clearModes();
	bool status() const { return modes != NULL; }

	int getNumModes() { return nModes; }
	void setMeanPCA(const BYU * byu);

	bool readPCAfile(const char * filename);
	bool pcaDeform(BYU * byu, std::vector<float> & vals);

private:

	PCA(const PCA&);				// Not implemented
	PCA & operator=(const PCA&);	// Not implemented

	int nModes;
	int nVerts;
	int numMeanVertices;
	float * meanVertices;
	float * modes;
};
#endif


#endif 

