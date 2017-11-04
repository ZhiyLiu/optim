#ifndef SURFACE_PATCH_ENSEMBLE_H
#define SURFACE_PATCH_ENSEMBLE_H

#include <math.h>
#include <time.h>
#include <vector>
#include <iostream>
#include "M3DObject.h"
#include "M3DObjectFile.h"
#include "M3DQuadFigure.h"
#include "Image3D.h"
#include "DistanceMap3D.h"
#include "RAWImageFile.h"
#include "Registry.h"
#include "Shapedepend.h"


enum HistogramPatchType
{
	GLOBAL_HPT,
	LOCAL_HPT
};

enum HistogramMatchType
{
	GAUSSIAN_HMT,
	WEIGHTED_GAUSSIAN_HMT
};

typedef struct QuantileComponentStruct
{
	double value;
	double weight;
} QuantileComponent;

typedef struct ProfileComponentStruct
{
	int index;
	double weight;
} ProfileComponent;

bool lessProfileComponent(ProfileComponent c1, ProfileComponent c2);
bool lessQuantileComponent(QuantileComponent c1, QuantileComponent c2);


typedef struct QuantileMixtureStruct
{
	double lowFrequency, midFrequency, highFrequency;
	std::vector<QuantileComponent>  lowComponents;
	std::vector<QuantileComponent>  midComponents;
	std::vector<QuantileComponent>  highComponents;
} QuantileMixture;

class SurfacePatch;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
// SurfacePatchEnsemble
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
class SurfacePatchEnsemble
{
public:

	SurfacePatchEnsemble(HistogramMatchType hmtIn, const char * filename, const char * QFDataFilename);
	void initialize(const M3DObject * object, Image3D * refImage);

	void setMean(double *newMean);  // replace mean with same-size mean
	double *getMean() {return mean;}

	void print();
	bool save(const char * filename);
	~SurfacePatchEnsemble();
	void initializeLocalRegions(const M3DObject *object, Image3D *refImage);

	bool update(Image3D * refImage, const M3DObject * oldObject, int figureID,
		int atomID = -1, const char * binaryImageFilename = NULL, const char * patchImageFilename = NULL);
	double getMatchValue(int figureID, int atomID);

	bool getBPoints(const M3DObject * object, int figureID, int inSurfaceLevel,
									int & numPoints, Bpoint *& Bpoints, bool copy = false);

	M3DObject * getSingleFigureObject(const M3DObject * oldObject, int figureID);
	static int histogramMatchTypeToInt(HistogramMatchType hmt);
	static HistogramMatchType intToHistogramMatchType(int hmtInt);

  // For recording local image match in the credibility application
	void setLocalLogFilename(const char* c) { localLogFilename = c; }

	HistogramMatchType hmt;

private:
	int histogramPatchTypeToInt(HistogramPatchType hpt);
	HistogramPatchType intToHistogramPatchType(int hptInt); 

	bool checkQFData(int * slices);
	double * getVoxelFeatures(Image3D *refImage, int index);
	void addElement(int patchID, double *features, double profileWeight, double pixelWeight, bool isInside);

	///////////////////////////////////////////////////////////////////////////
	// CONFIGURATION VARIABLES
	///////////////////////////////////////////////////////////////////////////

	Registry reg;

	// Variables that determine the regions
	HistogramPatchType histogramPatchType;
	double inDistance, outDistance,	inGaussianSTD, outGaussianSTD;

	// Sampling Variables
	int numInProfileSteps, numOutProfileSteps, voxelProfileSurfaceLevel;
	int usePixelNeighbors;

	// Local region definition
	int localPatchSurfaceLevel;
	double localPatchRadius, localPatchSTD;

	// Gas and Bone defaults
	int useLowHigh, lowThreshold, highThreshold;

	// Feature collection
	int numLowBins, numMidBins, numHighBins;
	int numFeatures; //neighborhoodSize
	//double **featureVectors;
	//int featureVectorLength;

	double *QFData;
	int QFDataDim[3];
	int QFDataMargin[3];

	///////////////////////////////////////////////////////////////////////////
	// DATA VARIABLES
	///////////////////////////////////////////////////////////////////////////
	int numFigures, numPatches;
	std::vector<SurfacePatch *> patches;
	ThallCode::Pointlist_server2 **pList; //One for each figure

	int numValues;
	double *values; //format: [lowW midW highW quantiles] conversion to/from W's is in Matlab code
	double *mean;
	double *scaleVector;

	int numVectors;
	double **vectors;
	int *vectorLengths;
	int *vectorOffsets;
	double *vectorVariances;
	int *vectorResidualID;
	int *vectorAtomID;
	int *vectorFigureID;

	int numResiduals;
	double *residualVariances;
	bool ** residualMasks; // computed not stored
	int * residualMaskLengths; // computed not stored

	bool returnWorstMatch;
	double * steps;
	Image3D * checkImage;

	bool initialized;

  // For recording local image match in the credibility application
	const char* localLogFilename;
};

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
// SurfacePatch
//////////////////////////////////////////////////////////////////////////////////////////////////////////////

class SurfacePatch
{
public:
	SurfacePatch(int inFigureID, int inAtomID, int inOffset, int numFeatures);
	~SurfacePatch();
	void clear();
	bool addElement(double *features, double weight, bool isInside, bool isLow, bool isHigh);
	bool simplify(double *values, double *mean, int numLowBins, int numMidBins, int numHighBins);

	//set once
	int figureID, atomID, offset, numFeatures;
	std::vector<ProfileComponent> profiles;

	//cleared and set with each call to update()
	QuantileMixture *mixtures;
private:
	bool simplifyHelper(double *values, double *mean, int numBins, int offset, double count, 
											double total, std::vector<QuantileComponent> &components, bool useSort);
};

#endif

