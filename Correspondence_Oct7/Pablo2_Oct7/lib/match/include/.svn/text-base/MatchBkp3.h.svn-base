#ifndef MATCH_H
#define MATCH_H

#include "SimilarityTransform3D.h"
#include "M3DObject.h"
#include "M3DQuadFigure.h"
#include "Image3D.h"
#include "Mask.h"
#include "SimpleMaskFile.h"
#include "renderDefinitions.h"
#include "DistanceMap3D.h"
#include "ImageDistanceMap.h"
#ifdef BINARY
#include "MatchUtility.h"
#endif
#include "matrix.h"
#include "optima.h"

#define D_XFERLIST
#define D_POINTLIST_SERVER2
#define D_SUBDIVSURF
#define D_DIATOM
#include "Shapedepend.h"

struct SurfacePointList
{
    Bpoint * points;
    int numPoints;
};

struct referencePoint
{
	double u, v, t;
    double pu, pv, pt;
	double distance;
};


#ifndef BINARY
extern const double DEFAULT_NORMAL_PENALTY_WEIGHT;
extern const double DEFAULT_SLIDING_PENALTY_WEIGHT;
extern const double DEFAULT_PENETRATION_PENALTY_WEIGHT;
extern const double DEFAULT_CREASING_PENALTY_WEIGHT;

#define MAX_NUM_FIG_MATCH_RESULTS	8
#define MAX_NUM_ATOM_MATCH_RESULTS	7
#define MAX_NUM_SPOKE_MATCH_RESULTS	6

#else	/* BINARY */

#define MAX_NUM_FIG_MATCH_RESULTS	11
#define MAX_NUM_ATOM_MATCH_RESULTS	12
#define MAX_NUM_SPOKE_MATCH_RESULTS	7

class DistanceVectorList;
class M3DAtomPredictor;

#endif

static const double ERR = 0.0000000001;


class SurfacePatchEnsemble;
class DQFMatch;
class MultiPatch;
class DistanceToPointSetFunction;

#ifdef SWIG
//
// Workaround for the lack of swig nested class support.
//

// Suppress SWIG warning about nested structures/classes
#pragma SWIG nowarn=SWIGWARN_PARSE_NESTED_CLASS

// Redefine nested class in global scope in order for SWIG to generate
// a proxy class. Only SWIG parses this definition.
struct matchResult {
	int tuningParm;		// Either bpTune_t or gpTune_t
	double value;		// Computed penalty (unweighted)
	double weight;		// Current penalty weight
	bool useIt;			// Active flag
	Components::FunctionClass_t fnClass;	// rough or smooth function.
	const char * name;	// Name to print for super-verbose output
};

%{
// SWIG thinks that Inner is a global class, so we need to trick the C++
// compiler into understanding this so called global type.
typedef Match::matchResult matchResult;
%}
#endif

#if defined(BINARY) && defined(SWIG)
// Someone forgot to define this function for binary.
%ignore Match::computeModelPenalty(const SimilarityTransform3D&);
#endif

class Match
{
public:
	struct matchResult {
		int tuningParm;		// Either bpTune_t or gpTune_t
		double value;		// Computed penalty (unweighted)
		double weight;		// Current penalty weight
		bool useIt;			// Active flag
		Components::FunctionClass_t fnClass;	// rough or smooth function.
		const char * name;	// Name to print for super-verbose output
	};

#ifndef BINARY
    Match(bool adaptive = false);
#else
    Match(M3DObject * landmarkObj);
#endif
    ~Match();

#ifndef BINARY
	void setHistogram(SurfacePatchEnsemble * h) { spe = h; }
	void setDQFMatch(DQFMatch * d) { dqf = d; }
	void setMultiPatch(MultiPatch * m) { mp = m; }
#endif
	void definePLists(M3DObject * object, int subdivLevel = -1);

#ifdef BINARY
	DistanceVectorList * getDistanceVectorList() {
		return distanceVectorList;
	}
#else
	bool initialize(M3DObject * _referenceObject, Image3D * _targetImage,
        bool buildMasks = true, Image3D * _trainingImage = NULL,
		bool useWindowing = true, MatchType _type = GAUSSIAN_DERIVATIVE_MATCH,
		int surfaceLevel = 2, const char * templateFilename = NULL,
		const char * profileFilename = NULL);

	void setWeight(int optimizerType, int tuningParm, double weight);
#endif
//    bool initialize(M3DObject * _referenceObject, Image3D * _targetImage,
//        const char * maskFilename, bool useWindowing = true);

#ifdef BINARY
	// Changed to Image3D to ImageDistanceMap in initialize
	bool initialize(M3DObject * _referenceObject, Image3D * _targetImage,
		ImageDistanceMap * _binaryDistanceMap, bool buildMasks = true,
		Image3D * _trainingImage = NULL, bool useWindowing = true,
        MatchType _type = GAUSSIAN_DERIVATIVE_MATCH, int surfaceLevel = 3,
		const char * templateFilename = NULL, const char * profileFilename = NULL);

	ImageDistanceMap * getImageDistMapPointer() { return binaryDistanceMap; }
    void setDistMap(ImageDistanceMap * _distMap) { binaryDistanceMap = _distMap; }

	// Dibyendu
	void setRefObjectAfterScale( M3DObject* _refObjectAfterScale ) { refObjectAfterScale = _refObjectAfterScale; }

#endif

	bool initContourMatch(const char* filename);

    bool buildMask(int figureId, MatchType _type, Image3D * trainingImage = NULL,
		int surfaceLevel = 2, const char * templateFilename = NULL, 
		const char * profileFilename = NULL);

#ifdef BINARY
	void setUpAtomPredictor(int figureId);
#endif

    // Access Functions
    MatchType getType() { return type; }
    M3DObject * getReferenceObject() { return referenceObject; }

	// Dibyendu
	M3DObject * getRefObjectAfterScale() { return refObjectAfterScale; }

    Image3D * getTargetImage() { return targetImage; }
    SimilarityTransform3D & getInitialTransform() { return initialTransform; }
    Mask * getMask(int figureId) { return masks[figureId]; }
	DMask * getDMask(int figureId) { return dmasks[figureId]; }
    std::vector<referencePoint> ** getConstrainedPointLists() { return constrainedPointLists; }
#ifndef BINARY
    double getNormalPenaltyWeight() { return normalPenaltyWeight; }
    double getSlidingPenaltyWeight() { return slidingPenaltyWeight; }
    double getPenetrationPenaltyWeight() { return penetrationPenaltyWeight; }
    double getCreasingPenaltyWeight() { return creasingPenaltyWeight; }
#endif

    // Set Functions
    void setNumMedialSubdivisions(int n) { numMedialSubdivisions = n; }

    // For the recorded transformation done before similarity registration
    void setInitialTransform(const SimilarityTransform3D * initXform);
    void setTargetImage(Image3D * _targetImage) { targetImage = _targetImage; }
    void setConstrainedPointLists(std::vector<referencePoint> ** _constrainedPointLists);

    // Resets the reference object; penalties are computed relative to it
    // (this should be called at the start of a new optimizer stage)
    void resetReferenceObject(M3DObject * _referenceObject);

#ifndef BINARY

    // Penalty weights for multi-object constraints
    void setNormalPenaltyWeight(double w) { normalPenaltyWeight = w; }
    void setSlidingPenaltyWeight(double w) { slidingPenaltyWeight = w; }
    void setPenetrationPenaltyWeight(double w) { penetrationPenaltyWeight = w; }

    // Boundary creasing penalty weight
    void setCreasingPenaltyWeight(double w) { creasingPenaltyWeight = w; }

    // Match Functions
	double computeModelMatch(M3DObject * targetObject);
    double computeModelMatch(const SimilarityTransform3D & transformation);
    double computeMainFigureMatch(M3DObject * targetObject, int figureId);
    double computeSubfigureMatch(M3DObject * object, int figureId);

    // Penalty Functions
    double computeModelPenalty(const SimilarityTransform3D & transformation);
    double computeMainFigurePenalty(M3DObject * targetObject, int figureId);
    double computeSubfigurePenalty(M3DObject * object, int figureId);
    double computeAtomPenalty(const M3DObject & targetObject, int figureId,
		int atomId, double neighborPenaltyWeight);

    double computeConstraintsPenalty(const M3DObject & targetObject, int figureId,
                                     bool useInterpenetration = true);

	double getWeightedPenalty(int optimizerType, int tuningParm);

#else	/* BINARY */

	double * getMaxDistance() { return maxDist; }

	double computeBinaryImageMatch(M3DObject * targetObject, int figureId,
		int enableDMApproximation, double angleTolerance, int subdivLevel);
    
	//Xiaojie
	double computeBinaryImageNormalMatch(M3DObject * targetObject, int figureId,
									  int subdivLevel);

	double computeFigureBinaryImagePlaneOrienMatchAtEnd(M3DObject * targetObject, int figureId);

	double computeAtomBinaryImagePlaneOrienMatchAtEnd(M3DObject * targetObject, int figureId, int atomId);

	double computeFigureBinaryImageVertexMatch(M3DObject * targetObject, int figureId);

	double computeAtomBinaryImageVertexMatch(M3DObject * targetObject, int figureId, int atomId);

    // Penalty Functions
    double computeModelPenalty(const SimilarityTransform3D & transformation);
	// Compares how this figure penetrates image regions it should
	// avoid (see "inter-penetration mask")
	double computeInterpenetrationPenalty(M3DObject * targetObject, int figureId);
	// These are the only functions useful in Binary Pablo.  They fold
	// together the match and penalty calculations.  Both return
	// 0 as the best match and higher values for worse matches.

#endif	/* BINARY */
    double computeSpokeMahamPenalties(const Vector & x, int figureId, int atomId);
    double computeAtomMahamPenalties(const Vector & x, int figureId = -1,
		int atomId = -1);
    Components computeFigureMatchAndPenaltiesAsComponents(M3DObject * targetObject,
		int figureId, M3DObject * similarityObject, const std::vector<double> & x,
		double & geomPenalty);

	Components computeAtomMatchAndPenaltiesAsComponents(M3DObject * targetObject,
		int figureId, int atomId, const Vector & x, double & geomPenalty);
	Components computeSpokeMatchAndPenaltiesAsComponents(M3DObject * targetObject,
		int figureId, int atomId, const Vector & x, double & geomPenalty);
    double computeFigureMatchAndPenalties(M3DObject * targetObject, int figureId,
		M3DObject * similarityObject, const std::vector<double> & x,
		double & geomPenalty)
	{
		const Components penalty =
			computeFigureMatchAndPenaltiesAsComponents(targetObject, figureId,
				similarityObject, x, geomPenalty);
		return penalty.rough + penalty.smooth;
	}

	double computeAtomMatchAndPenalties(M3DObject * targetObject, int figureId,
		int atomId, const Vector & x, double & geomPenalty)
	{
		const Components penalty =
			computeAtomMatchAndPenaltiesAsComponents(targetObject, figureId,
				atomId, x, geomPenalty);
		return penalty.rough + penalty.smooth;
	}
	double computeSpokeMatchAndPenalties(M3DObject * targetObject, int figureId,
		int atomId, const Vector & x, double & geomPenalty)
	{
		const Components penalty =
			computeSpokeMatchAndPenaltiesAsComponents(targetObject, figureId,
				atomId, x, geomPenalty);
		return penalty.rough + penalty.smooth;
	}

#ifdef BINARY
	// this function is a hack that can have serious side effects
	// we are assuming that the targetSurfacePoints computed by computeBinaryImageMatch
	// are at a high enough subdivLevel that they make sense here.
	// we make no guarantee that this targetSurfacePoints computed here (if needed)
	// are suitable for any other purpose -- levy & his attorney
    double computeLandmarkPenalty(M3DObject * targetObject, int figureId,
		M3DObject * liObject, double dilationFactorInModelUnits, int subdivLevel);
#endif

    // Initialize before deformation using inter-figural constraint
    bool initializeDistanceLists(M3DObject * object, int level);
	bool updateDistanceList(M3DObject * object, int figureId);
    bool initializeConstraintPointsLists(M3DObject * object, bool calcFigCoord = false);

    void testSeurat(M3DFigure * figure, int surf_level) const;

#ifndef BINARY
	bool isAdaptive() { return adapt; }
#endif

	const matchResult * getFigureStageResults() const { return figureResults; }
	const matchResult * getAtomStageResults() const { return atomResults; }
	const matchResult * getSpokeStageResults() const { return spokeResults; }

private:

	bool windowIntensities;
#ifndef BINARY
	bool adapt;
#endif

	ThallCode::Pointlist_server2 ** pList;
	int pListCount;

    M3DObject * referenceObject;

	// Dibyendu
	M3DObject * refObjectAfterScale ;	// saving a copy of the object right after scaling
										// because of the parameter FigureDoScale
										// To be used for extracting radius

    Image3D * targetImage;

	// Training mask, maybe
	SimpleMaskFile * simpleMask;

    // Mask for each figure
    std::vector<Mask  *> masks;

	std::vector<DMask *> dmasks;

    MatchType type;
	matchResult figureResults[MAX_NUM_FIG_MATCH_RESULTS];
	matchResult atomResults[MAX_NUM_ATOM_MATCH_RESULTS];
	matchResult spokeResults[MAX_NUM_SPOKE_MATCH_RESULTS];

	// Distance map for fitting to a binary image
	ImageDistanceMap * binaryDistanceMap;

    // Defines the fineness of the interpolation
    int numMedialSubdivisions;

    std::vector<SurfacePointList> surfacePoints; //this is for reference model in gray pablo
    SimilarityTransform3D initialTransform;
	SurfacePointList targetSurfacePoints; // this is used to lookup landmarks (all pablii) and binary image match

    // Inter-figural constraints
    std::vector<DistanceMap3D *> distanceMaps;
    int constraintsMapSize[3];
    std::vector<referencePoint> ** constrainedPointLists;
    int numPointLists;
	int surface_level_for_match;
    int constraints_surface_level;

#ifndef BINARY

    // Penalty weights for multi-object constraints
    double normalPenaltyWeight,
           slidingPenaltyWeight,
           penetrationPenaltyWeight;

    // Boundary creasing penalty weight
    double creasingPenaltyWeight;

	SurfacePatchEnsemble * spe;
	MultiPatch * mp;

	DQFMatch * dqf;

#else
	M3DAtomPredictor * atomPredictor;

	M3DObject * landmarkObject; // This contains landmark figures
	// LATER: allow match to compute LM_METHOD_OBJ penalty

	double maxDist[3];
	DistanceVectorList * distanceVectorList;	// Diagnostic distance vector list

#endif

	DistanceToPointSetFunction * contourMatch;

	void initializeAtomStage();
	double * atomMahamPenalty;
	int nAtomsMax;


	// this is to allow landmarks to point at any bpoint vertex
	// this depends on the figure being optimized
	int lmBpointIndices[100];
};



#endif

