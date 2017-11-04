#include "globalBuildParms.h"
#ifdef PCA_VISUALIZER


#ifdef WIN32
#include <windows.h>
#endif
#include <stdio.h>
#include <iostream>
#include <sys/stat.h>
#include "BYU.h"
#include "PCA.h"


#ifndef S_ISREG
#define S_ISREG(m)        (((S_IFREG & m) != 0) ? true : false)
#endif

using namespace std;


PCA::PCA()
{
	nModes = 0;
	nVerts = 0;
	meanVertices = NULL;
	modes = NULL;
}

PCA::~PCA()
{
	if (meanVertices != NULL)
		delete meanVertices;
	clearModes();
}

void PCA::clearModes()
{
	if (modes != NULL) 
		delete [] modes;
	modes = NULL;
}

void PCA::setMeanPCA(const BYU * byu)
{
	numMeanVertices = byu->vertexCount();
	meanVertices = new float[3*numMeanVertices];
	for (int i = 0; i < 3*numMeanVertices; i++)
			meanVertices[i] = byu->vertexList()[i];
}

bool PCA::readPCAfile(const char *filename)
{
	struct stat buf;
	int i, j;
	FILE *fp;

	if (filename == NULL)
		return false;
	else {
		if (stat(filename, &buf) != 0 || S_ISREG(buf.st_mode) == false) {
			cout << "Not a valid file.\n";
			return false;
		}
	}

	fp = fopen(filename, "r");
    if (fp == NULL)
        return false;

	clearModes();

	cout << "Loading PCA file: " << filename <<'\n';

    // Read the PCA residue file
    fscanf(fp, "%d %d\n", &nVerts, &nModes);

	if (nVerts == 0 || nModes == 0) {
		cout << "File contained either no vertices or no number of modes. \n";
		fclose(fp);
		return false;
	}

	modes = new float[(3*nVerts)*nModes];

	for (i = 0; i < (3*nVerts)*nModes ; i++)
		modes[i] = 0.0;

	for (i = 0; i < nModes ; i++) {
		for (j = i*nVerts; j < (i+1)*nVerts ; j++) {
			fscanf(fp, "%f %f %f\n", &modes[3*j], &modes[3*j+1], &modes[3*j+2]);
		}
	}

	if (fp != NULL) {
		fclose(fp);
		return true;
	}
	else {
		cout << "Not a valid file.\n";
		return false;
	}
}

bool PCA::pcaDeform(BYU * byu, std::vector<float> & vals)
{
	int i,j;

	if (byu == NULL) {
		cout << "NULL byu object is passed.\n";
		return false;
	}

	if (meanVertices == NULL) {
		cout << "PCA file is not loaded.\n";
		return false;
	}

	if (nModes > vals.size()) {
		for (i = vals.size() ; i < nModes ; i++) 
			vals.push_back(0.0);
	}

	float *vec = new float[3*nVerts];
	for (i = 0 ; i < 3*nVerts ; i++) {
		vec[i] = 0.0; 
	}

	for (i = 0 ; i < nModes ; i++) {
		for (j =0 ; j < 3*nVerts ; j++) 
			vec[j] += vals[i]*modes[j+i*3*nVerts];
	}

	if (nVerts != byu->vertexCount() || nVerts != numMeanVertices) {
		cout << "BYU object and meanBYU do not have the same number of vertices. \n";
		if (vec)
			delete [] vec;
		return false; 
	}

	float * verts = (float *) byu->vertexList();
	for (i = 0; i <3*nVerts ; i++) {
		verts[i] = meanVertices[i] + vec[i];
	}

	//cout<< vec[11] <<'\n';
	//cout<< mean[11] <<'\n';

	if (vec)
		delete [] vec;
	return true;
}


#endif	/* PCA_VISUALIZER */

