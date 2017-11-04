#ifndef DQFMATCH_H
#define DQFMATCH_H

#include <math.h>
#include <time.h>
#include <vector>
#include <iostream>
#include "M3DObject.h"
#include "M3DObjectFile.h"
#include "M3DQuadFigure.h"
#include "Image3D.h"
#include "DQFImage.h"
#include "DistanceMap3D.h"
#include "RAWImageFile.h"
#include "Registry.h"
#include "Shapedepend.h"

/*
Joshua Stough
MIDAG, Pizer/Chaney
April 2007

  This class is used in training and segmentation using distance quantile functions
  (DQF).  If Eli's SurfacePatchEnsemble was not littered with problem specific 
  improvements to histogram matching, one could possibly use that code for this
  (that is, until this code gets similarly littered).

  I use a DQuantileFunction implementation that does some of the generic things for 
  a quantile function, like store data and compute that data's qf.

  DQFTraining and DQFMatch use DQuantileFunction to perform their tasks.


  */



typedef struct
{
	double element;
	double weight;
}DDataElement;

//The class DQuantileFunction can store data and return the quantile function of it.
class DQuantileFunction
{
public:
	DQuantileFunction();
	DQuantileFunction(int size);
	~DQuantileFunction();

	void addData(DDataElement element);
	//add a data element to the list to consider.

	void addData(double element);
	//if the user doesn't wants DQuantileFunction to specify the element weight.

	void makeQF(int size = 0); 
	//Compute the quantile function for the current rawData.

	void setQFSize(int size);
	//sets the qf size and declares it needs to be remade;

	void clearData(){rawData.clear();}

	void clearQF(){delete []qf;}

	void useDefaultWeighting(double umu = 2.0, double uvar = 1.0);

	void useCustomWeighting(){useMyWeighting = false;}

	bool usingDefaultWeighting(){return useMyWeighting;}

	double *getQF();

	double getMean(); //common functionality;


private:
	double *qf;
	//the quantile function.

	std::vector<DDataElement> rawData;
	//The data from which the qf is computed.

	bool dirty;
	//Basically, is the qf older than rawData (as in, should the qf be recomputed).

	int qfSize;
	//the size of the expected quantile function.

	bool useMyWeighting;
	//In that I'm coding this for the purposes of distance quantile functions, 
	//the default weighting will be my 1 til mu then falloff wrt var.

	double mu, var;

};

//This class allows for the distance quantile training, whether that's
//writing out the dqf's for fit models in greyscales, or writing the 
//dqf image of a greyscale.
class DQFTraining
{
public:
	DQFTraining();
	DQFTraining(Registry & reg);

	~DQFTraining();

	void setDQFImageBoundingBox(Vector3D uXmin, Vector3D uXmax);
	void setDQFImageBoundingBox(Image3D *image, const M3DObject * object, int figureID);
	bool createDQFImage(Image3D *image, const char * outputFilename, const M3DObject * object = NULL, int figureID = 0);
	void setBoundingBoxBuffer(double buffer);

	bool doModelDQFTraining(Image3D * image, const M3DObject * object, int figureID, const char *saveFilename, const char *dqfFilename = NULL);


private:

	void floodFill(Image3D *image, int x, int y, int z, GreyValue &targetGrey, GreyValue &replacementGrey);
	int conCompSize(Image3D *target, Image3D *seen, int x, int y, int z, int & maxComponentSize);
	void deleteBoneBits(Image3D *target, Image3D *seen, int slice, int maxComponentSize);

	void loadConfig(Registry & reg);

	Vector3D Xmin, Xmax;
	double boundingBoxBuffer;

	DQFImage * dqfImage;

	double weightingMu, weightingVar;
	int sWidth;
	int modelSubDivLevel;

	double imageDelta, imageMultiplier;

};

//Finally, the match class itself, the part that Match.cpp will know. 
//Somewhere in here is the computeDQFPenalty(image, figure) that 
//other modular matches include.

//Just thought of something.  There are weighting params that affect the 
//quantile function.  Those parameters are not stored in the dqf Image.  
//Uh-oh...Well, we'll have to keep that in mind, maybe ensure the parameters
//are the same in training the dqf images and later (here) where 
//we compute a match.

//The parameters will be apparent in the dqf model training.  When in the match
//you sample from the dqf image, the dimension needs to match the dimension for 
//the stats you have.  As for the weighting parameters and so on, one needs to
//ensure that the params in the ImageMaker and ModelMaker config files are the 
//same.

//First I also want a DQFpca class to read and store the stats for a dqf model,
//as well as compute probabilities...just for a single point.  
class DQFpca
{
public:
	DQFpca();
	~DQFpca();

	double computeDQFMahal(double * qf);
	void readDQFpca(FILE * fin, int uwidth);
	int getNEigs();

private:
	double *mu;
	double *var;
	double **eigs;

	int neigs;
	int dim;

	//mu is d x 1, var is d x 1 with the elements being the variances of the 
	//eigenmodes, eigs is d x (<d), with a column vector for every non-zero
	//element of var minus one.  That last non-zero element of var is the 
	//residual variance.

};

class DQFMatch
{
public:
	DQFMatch();
	DQFMatch(Registry &reg);
	~DQFMatch();

	void initializeDQFMatch(Image3D *uimage, const M3DObject * object, const char ** dqfImageFilenames, const char ** dqfStatFilenames);
	void loadDQFpcasForFigure(const char *dqfStatFilename, int figureID);
	void loadConfig(Registry & reg);
	
	double computeDQFMatch(Image3D * image, const M3DObject * targetObject, int figureID);

	void setWeightingParams(double umu = 2.0, double uvar = 1.0);


private:

	double mu, var; //weighting the distances in the computed match.
	//Should be the same as those used to compute the dqfImage;

	Image3D * image;

	int surfaceLevel;

	int numFigs;

	DQFpca ** dqfPca;
	//each element of dqfPca is an array of DQFpca's, for the points on the surface of the figure.
	int * numPoints;
	//to be an array of ints, describing the number of points for each figure (the length of each
	//of the elements of dqfPca.

	int * width; 
	//the number of bins for the qf's for each figure.
	
	DQFImage ** dqfImage;
	//somewhere I'll say dqfImage = new DQFImage*[numFigs] or something, so one image per figure.
};


#endif

