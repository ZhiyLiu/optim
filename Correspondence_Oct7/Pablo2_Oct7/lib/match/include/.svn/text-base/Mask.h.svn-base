#ifndef MASK_H
#define MASK_H

#include "M3DObject.h"
#include "M3DFigure.h"
#include "Image3D.h"
#include "TemplateProfiles.h"
#include "SimpleMaskFile.h"
#include "renderDefinitions.h"

#define D_XFERLIST
#define D_POINTLIST_SERVERB
#define D_POINTLIST_SERVER2
#define D_SUBDIVSURF
#define D_DIATOM
#define NUM_MATCH_TYPES 7
#include "Shapedepend.h"
#include <vector>

extern const int MASK_MEDIAL_SUBDIVISIONS;

#define MASK_NUM_STEPS	5


typedef struct
{
    double x,y,z;
    double weight;
    float dist;
    int figureId;
    double u,v,t;
    double volelem;
    double tempVal;
    bool inBlendRegion;
} MaskElement;

// The commented-out terms were removed June 16, 2006
enum MatchType
{
    GAUSSIAN_DERIVATIVE_MATCH,
//	ABS_GAUSSIAN_DERIVATIVE_MATCH,
	NEG_GAUSSIAN_DERIVATIVE_MATCH,
//	NOTCH_MATCH,
//	PROFILE_MATCH,
//    TRAINING_IMAGE_WITH_INTERIOR_MATCH,
//    TRAINING_IMAGE_WITHOUT_INTERIOR_MATCH,
    SIMPLE_MASK_MATCH,
	BINARY_IMAGE_MATCH,
	HISTO_IMAGE_MATCH,
	MULTIPATCH_MATCH,
	NO_MATCH_TYPE
};

class Mask
{
public:
    Mask() {}
    Mask(M3DObject * referenceObject, int _figureId, float _cutoff,
         MatchType _type = GAUSSIAN_DERIVATIVE_MATCH,
         Image3D * traningImage = NULL, bool useWindowing = true,
		 int surfaceLevel = 2, const char * templateFilename = NULL,
		 const char * profileFilename = NULL,
		 SimpleMaskFile * simpleMask = NULL);
	bool status() { return ok; }	// Check to see if the constructor failed
    ~Mask();

    int getSize() { return maskElement.size(); }
    double getCutoff() { return cutoff; }
    double getSigma() { return sigma; }
    double getTotalVolume() { return totalVolume; }
    int getFigureId() { return figureId; }
	int getNumberOfPoints() { return NumberOfPoints; }
	int getSamplesPerPoint() { return SamplesPerPoint; }

    M3DFigure * getFigurePtr() const {
        return referenceObjectPtr->getFigurePtr(figureId);
    }
	M3DObject * getReferenceObjectPtr() {return referenceObjectPtr;}

	unsigned short getStackedImageMask() { return stackedImageMask; }

    bool interiorIsIncluded() { return includeInterior; }

    MaskElement * getElement(int i) { return maskElement[i]; }

    // Sets all the members of the mask, this is used for reading a mask
    // from a file. CAUTION: values are just set blindly.
    void setVals(int _figureId, double _cutoff, double _sigma,
        bool _includeInterior, const std::vector<MaskElement *> & elements);

    void calcPositions(M3DObject * object);
	double getRms(int i);
	double getMeanOffset(int i);
	double getMeanStd(int i);
	double getPMatchVar();
	double getMMatchVar();
#ifndef BINARY
	double getCenterIntensity();
#endif
	int getTemplateType(int i) { return mt.getProfileType(i); }

private:
	bool ok;
	TemplateProfiles mt;
	M3DFigureStats * figureStats;

	M3DObject * referenceObjectPtr;
	unsigned short stackedImageMask;

    std::vector<MaskElement *> maskElement;
	MatchType * profile;

	int numstep;
    double cutoff;
    double sigma;
	double notch_sigma;
    double totalVolume;

	int NumberOfPoints;
	int SamplesPerPoint;

    int figureId;

    bool includeInterior;
};

typedef struct 
{
	double area;
} SurfaceElement;

typedef struct
{
	double dist;
	double tempVal;
	double dTempVal;
} DerivMaskElement;



/* Derivative mask class used in the boundary optimization stage. 
*/
class DMask
{
public:
	double cutoffValues[2];


    DMask() {};
    DMask(M3DObject * object, int figureId, double _cutoff,
         MatchType _type = GAUSSIAN_DERIVATIVE_MATCH,
		 int subdivLevel = MATCH_POINT_CLOUD_SUBDIVISIONS,
		 const char * templateFilename = NULL, const char * profileFilename = NULL);
    ~DMask();

	MatchType getDmaskType() { return type;}
	void setDmaskType(MatchType _type) {type = _type;}
    int getMaskSize() { return dmaskElements.size(); }
	int getSurfaceSize(int level) { return surfaceElements[level].size(); }
    double getCutoff() { return cutoff; }
	double getCutoffValue(int j, int i);
    double getSigma() { return sigma; }
	double getTotalVolume() { return totalVolume; }

	void initializeSurfaceEls(M3DFigure * figure, int level);
	SurfaceElement * getSurfaceElement(int level, int i) 
		{ return surfaceElements[level][i]; }
	DerivMaskElement * getMaskElement(int i, int j);
	int getTemplateType(int i);

private:
	int figureId;
	double cutoff;
	double sigma;
	double notch_sigma;
	double totalVolume;

	int numstep;

	MatchType type;

	TemplateProfiles mt;

	std::vector<SurfaceElement *> surfaceElements[MAX_SUBDIV_LEVEL];

	std::vector<DerivMaskElement *> dmaskElements;

	// Another data structure must be included for profile matching,
	// where every type of filter must be available for the 
	// getMaskElement() function.

	DerivMaskElement ** dmaskElementsAll;
	double *cutoffValuesAll[2];
	// Profile will be a list (length equals the number of points on the
	// surface) of values 0, 1, 2, 3, ..., corresponding to the rows of 
	// dmaskElementsAll (containing all possible derivative filters).
};

const int MAX_NEIGHBOR_COUNT = 8;


#endif

