#ifndef OBJECTRELATIVESAMPLING_H
#define OBJECTRELATIVESAMPLING_H

// MSVC6: disable compiler warning about debugging symbol names too long
#pragma warning (disable:4786)

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
#include "SurfacePatchEnsemble.h"

#define MIN(A,B) (((A) < (B)) ? (A) : (B))
#define MAX(A,B) (((A) > (B)) ? (A) : (B))

/*
Joshua Stough
MIDAG, Pizer/Chaney
October 2007

A base class for flexible model-relative image sampling. 

Using this, you create a patch sampler, where parameters include
surface level of the samples, in/out distance in cm or radius-relative,
gaussian falloff in the weighting, basically very much like Eli's
stuff, but more flexible in that inside and outside patches are independent,
patches can be combined, etc.

  I'm going to recreate a lot of the code for the sampling, this is like thrid
  or fourth generation code for sampling the collar region using profiles.

  I'm going to remake the IQuantileFunction class from the DQuantileFunction
  in DQFMatch that I made, as the means of returning the signature of a patch.
  This is going to be huge though.  I want to be able to take qf's of patches, 
  then take qf's of combined patches and all that.  So I need to be store the
  voxel index that any given sample came from... whew...


  




  */



typedef struct
{
	double element;
	double weight;
	Vector3D coord;
}ImageElement;


struct compareImageElements
{
	bool operator() (const ImageElement &a, const ImageElement &b) const
	{
		if (a.coord.getZ() < b.coord.getZ())
			return true;
		else
		{
			if (a.coord.getZ() == b.coord.getZ())
			{
				if (a.coord.getY() < b.coord.getY())
					return true;
				else
				{
					if (a.coord.getY() == b.coord.getY())
					{
						if (a.coord.getX() < b.coord.getX())
							return true;
					}
				}
			}
			return false;
		}
	}
};


//The class QuantileFunction can store data and return the quantile function of it.
class QuantileFunction
{
public:
	QuantileFunction();
	QuantileFunction(int size);
	QuantileFunction(QuantileFunction q, int size);
	~QuantileFunction();

	void addData(ImageElement element);
	//add a data element to the list to consider.

	void addData(double element, Vector3D coord);
	//if the user doesn't wants QuantileFunction to specify the element weight.

	void addQF(const QuantileFunction &q);
	//Add the data from another quantile function to this one.

	void makeQF(int size = 0); 
	//Compute the quantile function for the current rawData.

	void setCutOffs(double inLowCutOff = -1.0e20, double inHighCutOff = 1.0e20); 
	//Set low and high cutoffs.

	void setQFSize(int size);
	//sets the qf size and declares it needs to be remade;

	int getQFSize() const {return qfSize;}

	void clearData(){rawData.clear();}

	//void clearDataStrong();
	//was to delete the elements too, but I'm not sure, with the structs, if I need it.

	void clearQF(){delete []qf;}

	void useDefaultWeighting(double umu = 2.0, double uvar = 1.0);

	void useCustomWeighting(){useMyWeighting = false;}

	bool usingDefaultWeighting(){return useMyWeighting;}

	double *getQF(int size = 0);
	const double *getQF(int size = 0) const;

	double getMean(); //common functionality;

	int getSize() 
	{
		return rawData.size();
	}

	//high and low cutoffs on the data, this is the qf's responsibility, but it doesn't need to know why.
	//so when you make the qf, you can send high and low cutoff's and it will set the proportion of the weight
	//in those extremes.

	//void cutOffExtremes(double lowCutoff, double highCutOff);
	double * getWeights();
	double getMiddleWeight() const 
	{
		return (totWeight - highWeight - lowWeight);
	}


private:
	double *qf;
	//the quantile function.

	//std::set<ImageElement, compareImageElements> rawData;
	std::vector<ImageElement> rawData;
	//The data from which the qf is computed.

	bool dirty;
	//Basically, is the qf older than rawData (as in, should the qf be recomputed).

	int qfSize;
	//the size of the expected quantile function.

	bool useMyWeighting;
	//In that I'm coding this for the purposes of distance quantile functions, 
	//the default weighting will be my 1 til mu then falloff wrt var.

	double lowWeight, highWeight, totWeight;

	double lowCutOff, highCutOff;

	double mu, var;

};




/************************************************************************************
 
   Now is the PatchPCA class. This class stores the appearance model for a patch. 
   It, as the QuanTileFunction, should be a data member of a Patch object (below).

  That way, higher up (lower down in the the file), we can ask a figure to compute a 
  match, and the figure in turn asks each of its Patches for its match, and the patches
  in turn ask their patchPCA to compute a match given the qf and gas/bone/total weights.

  NO WAIT!

  It seems that only after I get into the matlab training can I see which local 
  patches I'm going to want to put together, and during optimization adding patches
  and comparing to the stats.

  There is a patch pca, which is simply pca on the qf's of a patch defined by the anchor
  points and inside/outside (but then I can also imagine grouping an inside and outside 
  region).  A cluster type appearance from training can be represented as such a 
  PatchPCA.  Then there is a separate structure that will store the variability in the 
  cluster delineation on the surface, a regional flexibility object.  

  I'm trying to write this such that local pca is a limiting case, where a patchPCA has 
  one anchor point and the stats for it are for that local region, and there is no 
  regional flexibility object...


************************************************************************************/


class PatchPCA
{
public:
	PatchPCA();
	~PatchPCA();

	double computeMahal(double * qf, bool normalize = false);
	//void readPatchPCA(FILE * fin, int uwidth);
	double computeMatch(double * qf, double *w = NULL, bool normalize = false);
	double computeGasMatch(double *w = NULL, bool normalize = false);
	double computeBoneMatch(double *w = NULL, bool normalize = false);

	void setPca(double *inPca, int inDim, int length);
	void setWeights(double *inWeights, int length);

	void setGasBoneMatters(double inGM, double inBM)
	{
		gasMatters = inGM;
		boneMatters = inBM;
	}

	int getGasBoneCount()
	{
		int count = 0;
		if (gasMatters) count ++;
		if (boneMatters) count ++;
		return count;
	}

	double getGasBoneMatters()
	{
		return gasMatters + boneMatters;
	}

	int getNEigs()
	{
		return neigs;
	}

	int getDim()
	{
		return dim;
	}

	double ** projectData(double ** qfs, int num);


private:
	double *mu;
	double *var;
	double **eigs;

	double *weights;
	
	//double *priors;
	//int type;
	//These two are for MATCH_SHIFTING.

	int neigs;
	int dim;

	//mu is d x 1, var is d x 1 with the elements being the variances of the 
	//eigenmodes, eigs is d x (<d), with a column vector for every non-zero
	//element of var minus one.  That last non-zero element of var is the 
	//residual variance.

	double gasMatters;
	double boneMatters;
	//Cause I want some flexibility here.

};




/************************************************************************************
 

   Now is the Patch class.  A patch object should know how to sample itself.  It assumes
   you give it the appropriate figure for it to ask profiles of. 

  It's strange, I'm doing this because I cannot get Eli's surfacepatchensemble stuff to 
  do what I want. So this is simplified over the many optimizations he had.  

  I have a vision of a piece of code where it says for every patch for this object, ask
  the patch to sample itself. then move on having the patches compute their quantile 
  functions (which is having the patch ask its data member qf to compute itself)...you 
  know, just as many black boxes as possible, or as I think of at least.


************************************************************************************/

enum PatchType
{
	LOCAL_PATCHES,
	DEFINED_PATCHES,
	FLEXIBLE_PATCHES
};
//These are sampling types.
//Local is simple enough, a pca for every patch at some surface level.  Classic 
//clustered is implemented with this, but in training each patch is assigned a pca
//from a set of cluster types, or possible pcas.
//Defined is a pca per patch, but the patches are specially defined by the user.
//Flexible is a placeholder.

enum StatType
{
	MATCH_DEFAULT,
	MATCH_MLE,
	MATCH_SHIFTING,
	MATCH_UNKNOWN
};
//There are also statTypes, like
//MATCH_DEFAULT, which assumes a pca at each patch. Classic clustered is this case
//	where, as described above, there is a pca per patch given.
//MATCH_MLE takes the pca models provided and determines during the getYourMatch
//which one gives the best answer. There is no prior on the model used at a patch.
	//at least for now, since I can add more blocks to the pca's later.




typedef struct
{
	double samplingDistance, distanceSTD;
	int numSteps, samplingProfileSurfaceLevel, localPatchSurfaceLevel;
	int emSize;

	double localPatchRadius, localPatchSTD;
	double lowCutOff, highCutOff; 
	bool isInside, isRRelative;//see below:, gasMatters, boneMatters;
	PatchType patchType;
	StatType statType;
	
	int numPatches;
	double **delins; //for the DEFINED_PATCHES option.
	double *gasMatters;
	double *boneMatters;
	//gas and bone-Matters are now arrays, either a single element long, for all
	//the patches, or numPatches long, for each of the DEFINED_PATCHES.

	std::vector<double> taus; //the steps to sample a profile.
	std::vector<double> tauWeights; //the weights of the samples of a profile.

}SamplingConfiguration; //for a figure and side.


struct compareProfileComponents
{
	bool operator() (const ProfileComponent &a, const ProfileComponent &b) const
	{
		return (a.index < b.index);
	}
};

class Patch
{
public:
	Patch(int inFigureID, int inSpokeID, Bpoint bp, SamplingConfiguration *inSc);
	Patch();
	~Patch();

	void clear();

	void sampleYourself(Image3D *image, Bpoint *Bpoints, Image3D *patchImage, GreyValue g, Image3D *binImage = NULL);
	
	double * getQF() {return data.getQF();}

	void addPatch(const Patch &p);

	//Simple pass throughs.  The outside doesn't need to interact with a qf, it 
	//can just assume the Patch does all this stuff.

	void setCutOffs(double inLowCutOff = -1.0e20, double inHighCutOff = 1.0e20)
	{ data.setCutOffs(inLowCutOff, inHighCutOff); }

	double * getWeights()
	{ return data.getWeights(); }

	double getMiddleWeight() const 
	{
		return data.getMiddleWeight();
	}

	void addProfile(const ProfileComponent pc)
	{
		profiles.insert(profiles.begin(), pc);
	}

	double getProfileWeight(int index);

	void addAnchorAndSpoke(const Bpoint *bp, int inSpokeID)
	{
		anchorPoints.push_back(*bp);
		spokeIDs.push_back(inSpokeID);
	}

	void setSC(SamplingConfiguration *inSc)
	{
		sc = inSc;
		data.setQFSize(sc->emSize);
		data.setCutOffs(sc->lowCutOff, sc->highCutOff);
	}

	SamplingConfiguration * getSC()
	{return sc;}

	void setPCA(PatchPCA *inPca)
	{
		pca = inPca;
	}

	void setFigureID(int inFigureID)
	{ 
		figureID = inFigureID;
	}

	int numSpokeIDs()
	{
		return spokeIDs.size();
	}

	int * getSpokeIDs()
	{
		int * arr = new int[spokeIDs.size()];
		for (int i = 0; i < spokeIDs.size(); i ++)
			arr[i] = spokeIDs[i];

		return arr;
	}

	int getFirstSpokeID()
	{
		return spokeIDs[0];
	}

	Bpoint getFirstAnchorPoint()
	{
		return *anchorPoints.begin();
	}

	int getEMSize()
	{
		return sc->emSize;
	}

	int getDataSize()
	{
		return data.getSize();
	}

	double getYourMahal(bool normalize = false)
	{
		return pca->computeMahal(getQF(), normalize);
	}

	double getYourGas(bool normalize = false)
	{
		return pca->computeGasMatch(getWeights(), normalize);
	}

	double getYourBone(bool normalize = false)
	{
		return pca->computeBoneMatch(getWeights(), normalize);
	}
	
	double getYourMatch(bool normalize = false)
	{
		return pca->computeMatch(getQF(), getWeights(), normalize);
	}

	int getGasBoneCount()
	{
		return pca->getGasBoneCount();
	}

	double getGasBoneMatters()
	{
		return pca->getGasBoneMatters();
	}

	int getNEigs()
	{
		return pca->getNEigs();
	}

private:

	QuantileFunction data;

	SamplingConfiguration *sc;

	PatchPCA * pca;

	std::vector<Bpoint> anchorPoints; //the anchor point for this patch.
	std::vector<int> spokeIDs; //The ID's of the anchor points for this patch.
	int figureID;
	std::set< ProfileComponent, compareProfileComponents > profiles; //profiles to sample the patch.

};


//This object responsible for shifting on one side of one figure, a bit easier to think about.
class RegionShift
{
public:
	RegionShift(PatchPCA *inProjPCA, std::vector< PatchPCA * > inTypePCAs, double *inP, 
		int priorSize, SamplingConfiguration *inSc);
	~RegionShift();

	void buildAdjacency(const M3DObject * object, int f, int inSurfL);

	double getMatch(double **qfs, double &totModes);

	double getCombinedMatch(double **qfs, double &totModes);
	double getCombinedGasBoneMatch(double &totModes);

	void setPatches ( std::vector< Patch *> inPatches)
	{
		//patches.insert(patches.begin(), inPatches);
		patches = inPatches;

	}

	int getNumTypes() {return numTypes;}

private:

	void doShifting(double ** cX);
	void doShiftingCombined(double ** cX);
	double computeShiftedMatch(const int *typePerPoint, double ** cX, bool normalize = false);
	double computeShiftedMatchCombined(const int *typePerPoint, double ** qfs, bool normalize = false);
	void computeWeights();
	void Permute(int n = 200);

	PatchPCA *projPCA; //the projection pca.
	std::vector< PatchPCA * > typePCAs; //to store the pcas of the clusters.
	int *initModel;
	int *finalModel;
	double *priors;
	int numTypes;
	int **adj;
	int surfL;
	int numPts;
	int *permutation;
	int numModes;
	SamplingConfiguration *sc;

	double *weights;

	std::vector< Patch *> patches;
};


class MultiPatch
{
public:
	MultiPatch() : pList(0), patchImage(0), loggingFileName(0) {} ;
	~MultiPatch();
	MultiPatch(const char * filename);
	MultiPatch(const M3DObject * object, Image3D * refImage, const char * filename);
	//Will construct the Patch structures needed for sampling each figure of the object.


	void initializePatches(const M3DObject * object, Image3D * refImage);

	bool saveMP(const char * saveFilename);

	void sampleYourself(const M3DObject * object, Image3D * image, int f, Image3D * binImage = NULL, bool debug = false);

	double getYourMatch(const M3DObject * object, Image3D * image, int f, Image3D * binImage = NULL);

	double getYourMatch(const M3DObject * object, Image3D * image, int f, 
		std::vector<std::vector< int> > insides, std::vector<std::vector< int> > outsides, Image3D * binImage = NULL);

	void flexibleMatch(const M3DObject * object, Image3D * image, int f, int whichSide, 
		double &theMatch, double &total, Image3D * binImage = NULL);

	void MLEMatch(const M3DObject * object, Image3D * image, int f, int whichSide, 
		double &theMatch, double &total, Image3D * binImage = NULL);

    // This is so that we can instrument the match, and record the value
	// for each patch
	char *loggingFileName;


private: 

	bool getBPoints(const M3DObject *object, int figureID, int whichSide,
							int inSurfaceLevel, int &numPoints, Bpoint * &Bpoints, bool copy);


	PatchType intToPatchType(int n);
	int patchTypeToInt(PatchType p);

	StatType intToStatType(int n);
	int statTypeToInt(StatType s);

	void updateTauAndWeights(int figureID, int whichSide, Image3D *refImage);

	int whichDelin(Bpoint bp, SamplingConfiguration *sc, int *matching, bool isTube = false);

	void readSamplingConfiguration(const char * filename);

	void initializeFigurePatches(const M3DObject * object, Image3D * refImage, 
								int figureID, int whichSide);





	Registry reg;

	std::vector<std::vector< Patch *> > patches;  //one inside and one outside vector of patches for each figure.
	//std::vector<std::vector< Patch *> > outsidePatches;  //one outside vector of patches for each figure.

	std::vector<std::vector< PatchPCA *> > patchPCAs; 
	//I'm thinking of this as a parallel vector of pca data. These aren't used through MultiPatch but rather 
	//through the patches.  I want the patch to know it's match. If during the readSamplingConfiguration there
	//are stats, I will pick them up and put them in patchPCAs.  Then, during initializePatches, place their 
	//references in the patches.  Then, finally, this MultiPatch can ask patch->getYourMatch() at runtime.

	//configuration variables, which describe how the patches are sampled.  Can be different 
	//for each figure inside and outside (not currently for each patch), it's all picked up 
	//from the reg file.
	//SamplingConfiguration ** samplingConfigs; 
	//std::vector< SamplingConfiguration > samplingConfigs;
	SamplingConfiguration * samplingConfigs;
	//so, an inside and outside sampling config for each figure.
	
	char scFilename[2048];

	//more variables
	int numFigures; 

	ThallCode::Pointlist_server2 **pList;  //one for each figure.

	Image3D *patchImage;		//so that once it's first used, no need to re alloc it.
	GreyValue ** randColors;     //same as patchImage.

	//all the below used in the MATCH_SHIFTING.
	std::vector< RegionShift * > regionShifts;


};















#endif
