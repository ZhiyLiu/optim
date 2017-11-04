#include <math.h>
#include <assert.h>


// dibyendu
#include <vector>

#include <typeinfo>
#include <vtksrepinterpolatecrestspokesquartic.h>
#include <vtksrepinterpolatemedialspokeshermite.h>
#include <vtkQuad.h>
#include <vtkCellArray.h>
#include <vtkLine.h>

#include "globalBuildParms.h"
#include "ImageDistanceMap.h"
//#include "snapshot.h"
#include "M3DPrimitiveCorrector.h"
#include "ccl.h"
#include "Mask.h"
//#include "MaskFile.h"
#include "SurfacePatchEnsemble.h"
#include "DQFMatch.h"
#include "ObjectRelativeSampling.h"

#define D_PLIST_SUBDIVCOMP
#include "Shapedepend.h"

#define _BP_MATCH

#include "ImageResample3D.h"
#include "DistanceVectorList.h"
#include "M3DMainFigureProblem.h"
#include "M3DAtomPredictor.h"
#include "Tuning.h"
#include "OptimizerBase.h"
#include "Match.h"
#include "P3DControl.h"
#include "DistanceToPointSetFunction.h"
#ifdef SURFACE_VISUALIZER
#include "M3DObjectSurfaceVisualizer.h"
#endif

#define WEIGHTED

//#define DEBUG_GREY_ATOM_TUBE
//#define DEBUG

// Uncomment to eliminate use of the match penalty in regular Pablo
//#define NO_MATCH_PENALTY


using namespace std;

const double FAILED = 1.0e100;
const double WORST_MATCH = -1.0e100;

const double MASK_CUTOFF = 0.3;

#ifndef BINARY
const double DEFAULT_NORMAL_PENALTY_WEIGHT = 0.5;
const double DEFAULT_SLIDING_PENALTY_WEIGHT = 0.5;
const double DEFAULT_PENETRATION_PENALTY_WEIGHT = 1.0;
const double DEFAULT_CREASING_PENALTY_WEIGHT = 0.0005;
const double TUBE_GREY_ATOM_SRAD_PENALTY_WEIGHT = 2;
const double TUBE_GREY_ATOM_CURVINESS_PENALTY_WEIGHT = 10;
const double TUBE_GREY_ATOM_MODEL_PENALTY_WEIGHT = 1000;
#endif

extern int globalVerbosity;


#ifdef BINARY
Match::Match(M3DObject * landmarkObj) : landmarkObject(landmarkObj)
#else
Match::Match(bool adaptive)
#endif
{
	int i;

#ifndef BINARY
	adapt = adaptive;
#endif
	referenceObject = NULL;

	// Dibyendu
	refObjectAfterScale = NULL ;

	constrainedPointLists = NULL;
	numPointLists = 0;

	targetImage = NULL;
	contourMatch = NULL;

#ifdef BINARY
	atomPredictor = NULL;
	surface_level_for_match = 1;

	numMedialSubdivisions = MASK_MEDIAL_SUBDIVISIONS;
	maxDist[0] = 0;
	maxDist[1] = 0;
	maxDist[2] = 0;

	distanceVectorList = new DistanceVectorList;

#else
	constraints_surface_level = 0;
	surface_level_for_match = MATCH_POINT_CLOUD_SUBDIVISIONS;

	numMedialSubdivisions = MASK_MEDIAL_SUBDIVISIONS;

	normalPenaltyWeight = DEFAULT_NORMAL_PENALTY_WEIGHT;
	slidingPenaltyWeight = DEFAULT_SLIDING_PENALTY_WEIGHT;
	penetrationPenaltyWeight = DEFAULT_PENETRATION_PENALTY_WEIGHT;
	creasingPenaltyWeight = DEFAULT_CREASING_PENALTY_WEIGHT;

	spe = NULL;
#endif
	atomMahamPenalty = NULL;
	nAtomsMax = 0;
	pList = NULL;
	pListCount = 0;

	binaryDistanceMap = NULL;

#ifdef BINARY

	assert( 11 == MAX_NUM_FIG_MATCH_RESULTS );
	figureResults[0].tuningParm = BpFigureImageMatch;
	figureResults[0].fnClass	= Components::ROUGH;
	figureResults[0].name = "Image";
	figureResults[1].tuningParm = BpFigureLandmarkMatch;
	figureResults[1].fnClass	= Components::SMOOTH;
	figureResults[1].name = "Landmark";
	figureResults[2].tuningParm = BpFigureModelMatch;
	figureResults[2].fnClass	= Components::SMOOTH;
	figureResults[2].name = "Model";	
	figureResults[3].tuningParm = BpFigureInterpenetrationMatch;
	figureResults[3].fnClass	= Components::ROUGH;
	figureResults[3].name = "Penetration";
	figureResults[4].tuningParm = FigureMahalanobisMatch;
	figureResults[4].fnClass	= Components::SMOOTH;
	figureResults[4].name = "Mahalanobis";
	figureResults[5].tuningParm = BpFigureRSradPenalty;
	figureResults[5].fnClass	= Components::SMOOTH;
	figureResults[5].name = "R*Srad";
	figureResults[6].tuningParm	= BpFigureCurvinessPenalty;
	figureResults[6].fnClass	= Components::SMOOTH;
	figureResults[6].name = "Curviness";
	figureResults[7].tuningParm = FigureContourMatch;
	figureResults[7].fnClass	= Components::SMOOTH;
	figureResults[7].name = "Contours";
	//Xiaojie
	//fnClass not sure
	figureResults[8].tuningParm = BpFigureImageNormalMatch;
	figureResults[8].fnClass	= Components::ROUGH;
	figureResults[8].name = "Normal";
	figureResults[9].tuningParm = BpFigureImagePlaneOrienAtEnd;
	figureResults[9].fnClass	= Components::ROUGH;
	figureResults[9].name = "Plane";
	figureResults[10].tuningParm = BpFigureImageVertexMatch;
	figureResults[10].fnClass	= Components::ROUGH;
	figureResults[10].name = "Vertex";


	assert( 12 == MAX_NUM_ATOM_MATCH_RESULTS );
	atomResults[0].tuningParm	= BpAtomAverageNeighborMatch;
	atomResults[0].fnClass		= Components::SMOOTH;
	atomResults[1].tuningParm	= BpAtomImageMatch;
	atomResults[1].fnClass		= Components::ROUGH;
	atomResults[2].tuningParm	= AtomMahalanobisMatch;
	atomResults[2].fnClass		= Components::SMOOTH;
	atomResults[3].tuningParm	= BpAtomLandmarkMatch;
	atomResults[3].fnClass		= Components::SMOOTH;
	atomResults[4].tuningParm	= BpAtomModelMatch;
	atomResults[4].fnClass		= Components::SMOOTH;
	atomResults[5].tuningParm	= BpAtomInterpenetrationMatch;
	atomResults[5].fnClass		= Components::ROUGH;
	atomResults[6].tuningParm	= BpAtomRSradPenalty;
	atomResults[6].fnClass		= Components::SMOOTH;
	atomResults[7].tuningParm	= BpAtomDeltaNeighborMatch;
	atomResults[7].fnClass		= Components::SMOOTH;
	atomResults[8].tuningParm	= BpAtomCurvinessPenalty;
	atomResults[8].fnClass		= Components::SMOOTH;
	//Xiaojie
	//fnClass not sure
	atomResults[9].tuningParm = BpAtomImageNormalMatch;
	atomResults[9].fnClass	= Components::ROUGH;
	atomResults[10].tuningParm = BpAtomImagePlaneOrienAtEnd;
	atomResults[10].fnClass	= Components::ROUGH;
	atomResults[11].tuningParm = BpAtomImageVertexMatch;
	atomResults[11].fnClass	= Components::ROUGH;

	// dibyendu 
	//srep stage results

	assert( 12 == MAX_NUM_SREP_MATCH_RESULTS );
	sRepResults[0].tuningParm	= BpSRepAverageNeighborMatch;
	sRepResults[0].fnClass		= Components::SMOOTH;
	sRepResults[1].tuningParm	= BpSRepImageMatch;
	sRepResults[1].fnClass		= Components::ROUGH;
	sRepResults[2].tuningParm	= SRepRadiusMatch;
	sRepResults[2].fnClass		= Components::ROUGH;
	sRepResults[3].tuningParm	= BpSRepLandmarkMatch;
	sRepResults[3].fnClass		= Components::SMOOTH;
	sRepResults[4].tuningParm	= BpSRepModelMatch;
	sRepResults[4].fnClass		= Components::SMOOTH;
	sRepResults[5].tuningParm	= BpSRepInterpenetrationMatch;
	sRepResults[5].fnClass		= Components::ROUGH;
	sRepResults[6].tuningParm	= BpSRepRSradPenalty;
	sRepResults[6].fnClass		= Components::SMOOTH;
	sRepResults[7].tuningParm	= BpSRepDeltaNeighborMatch;
	sRepResults[7].fnClass		= Components::SMOOTH;
	sRepResults[8].tuningParm	= BpSRepCurvinessPenalty;
	sRepResults[8].fnClass		= Components::SMOOTH;
	//Xiaojie
	//fnClass not sure
	sRepResults[9].tuningParm = BpSRepImageNormalMatch;
	sRepResults[9].fnClass	= Components::ROUGH;
	sRepResults[10].tuningParm = BpSRepImagePlaneOrienAtEnd;
	sRepResults[10].fnClass	= Components::ROUGH;
	sRepResults[11].tuningParm = BpSRepImageVertexMatch;
	sRepResults[11].fnClass	= Components::ROUGH;

	assert( 7 == MAX_NUM_SPOKE_MATCH_RESULTS );
	spokeResults[0].tuningParm	= BpSpokeAverageNeighborMatch;
	spokeResults[0].fnClass		= Components::SMOOTH;
	spokeResults[1].tuningParm	= BpSpokeImageMatch;
	spokeResults[1].fnClass		= Components::ROUGH;
	spokeResults[2].tuningParm	= SpokeMahalanobisMatch;
	spokeResults[2].fnClass		= Components::SMOOTH;
	spokeResults[3].tuningParm	= BpSpokeLandmarkMatch;
	spokeResults[3].fnClass		= Components::SMOOTH;
	spokeResults[4].tuningParm	= BpSpokeModelMatch;
	spokeResults[4].fnClass		= Components::SMOOTH;
	spokeResults[5].tuningParm	= BpSpokeInterpenetrationMatch;
	spokeResults[5].fnClass		= Components::ROUGH;
	spokeResults[6].tuningParm	= BpSpokeRSradPenalty;
	spokeResults[6].fnClass		= Components::SMOOTH;

#else

	// The order or parameters should be (ideally) from smallest to largest,
	// after the weights are applied, so that arithmetic errors are minimized.

	assert( 8 == MAX_NUM_FIG_MATCH_RESULTS );
	// Set the default values of penalties for figure and mainFigure stages
	figureResults[0].tuningParm = GpFigureInterObjectWeight;
	figureResults[0].fnClass	= Components::SMOOTH;
	figureResults[0].weight = tuningWt(GpFigureInterObjectWeight);
	figureResults[0].name = "cnstr";

	figureResults[1].tuningParm = GpFigureGeometryWeight;
	figureResults[1].fnClass	= Components::SMOOTH;
	figureResults[1].weight = tuningWt(GpFigureGeometryWeight);
	figureResults[1].name = "dsq";

	figureResults[2].tuningParm = FigureMahalanobisMatch;
	figureResults[2].fnClass	= Components::SMOOTH;
	figureResults[2].weight = tuningWt(FigureMahalanobisMatch);
	figureResults[2].name = "mahal";

	figureResults[3].tuningParm = FigureContourMatch;
	figureResults[3].fnClass = Components::SMOOTH;
	figureResults[3].weight = tuningWt(FigureContourMatch);
	figureResults[3].name = "contours";

	// The image match penalties go below here
#define FIRST_IMAGE_MATCH_FIGURE_TERM	4

	figureResults[4].tuningParm = GpFigureHistogramMatch;
	figureResults[4].fnClass	= Components::ROUGH;
	figureResults[4].weight = tuningWt(GpFigureHistogramMatch);
	figureResults[4].name = "hist";

	// Note: figureResults[4].tuningParm is never used.
	// GpFigureProfileMatch is only defined for calls to getWeightedPenalty().
	figureResults[5].tuningParm = GpFigureProfileMatch;	// Mask-based image match
	figureResults[5].fnClass	= Components::ROUGH;
	figureResults[5].weight = 1.0;	// See comment above
	figureResults[5].name = "prof";

	figureResults[6].tuningParm = FigureDQFMatch;
	figureResults[6].fnClass	= Components::ROUGH;
	figureResults[6].weight = tuningWt(FigureDQFMatch);
	figureResults[6].name = "dqf";

	figureResults[7].tuningParm = GpFigureMultiPatchMatch;
	figureResults[7].fnClass	= Components::ROUGH;
	figureResults[7].weight = tuningWt(GpFigureMultiPatchMatch);
	figureResults[7].name = "mp";



	// --------------------------------------
	// Set the default values of penalties for atom (deformation) stage

	assert( 7 == MAX_NUM_ATOM_MATCH_RESULTS );
	// Note: atomResults[4].tuningParm is never used.
	// GpFigureProfileMatch is only defined for calls to getWeightedPenalty().
	atomResults[0].tuningParm = GpAtomProfileMatch;	// Mask-based image match
	atomResults[0].fnClass	= Components::ROUGH;
	atomResults[0].weight = 1.0;	// See comment above
	atomResults[0].name = "prof";

	atomResults[1].tuningParm = GpAtomHistogramMatch;
	atomResults[1].fnClass	= Components::ROUGH;
	atomResults[1].weight = tuningWt(GpAtomHistogramMatch);
	atomResults[1].name = "hist";

	atomResults[2].tuningParm = GpAtomMultiPatchMatch;
	atomResults[2].fnClass	= Components::ROUGH;
	atomResults[2].weight = tuningWt(GpAtomMultiPatchMatch);
	atomResults[2].name = "mp";

	// The image match penalties go above here
#define LAST_IMAGE_MATCH_ATOM_TERM	2

	atomResults[3].tuningParm = AtomMahalanobisMatch;
	atomResults[3].fnClass	= Components::SMOOTH;
	atomResults[3].weight = tuningWt(AtomMahalanobisMatch);
	atomResults[3].name = "mahal";

	atomResults[4].tuningParm = GpAtomGeometryWeight;
	atomResults[4].fnClass	= Components::SMOOTH;
	atomResults[4].weight = tuningWt(GpAtomGeometryWeight);
	atomResults[4].name = "dsq";

	// Note: atomResults[4] is not an ordinary matchResult entry.
	// Its weight is passed to computeAtomPenalty(), which is used
	// compute the GpAtomGeometryWeight, and atomResults[4].useIt
	// is never set to true.
	atomResults[5].tuningParm = GpAtomNeighborWeight;
	atomResults[5].fnClass	= Components::SMOOTH;
	atomResults[5].weight = tuningWt(GpAtomNeighborWeight);
	atomResults[5].name = "neigh";

	atomResults[6].tuningParm = GpAtomInterWeight;
	atomResults[6].fnClass	= Components::SMOOTH;
	atomResults[6].weight = tuningWt(GpAtomInterWeight);
	atomResults[6].name = "cnstr";

	// --------------------------------------
	// Set the default values of penalties for spoke (deformation) stage

	assert( 6 == MAX_NUM_SPOKE_MATCH_RESULTS );
	// Note: spokeResults[4].tuningParm is never used.
	// GpFigureProfileMatch is only defined for calls to getWeightedPenalty().
	spokeResults[0].tuningParm = GpSpokeProfileMatch;	// Mask-based image match
	spokeResults[0].fnClass	= Components::ROUGH;
	spokeResults[0].weight = 1.0;	// See comment above
	spokeResults[0].name = "prof";

	spokeResults[1].tuningParm = GpSpokeHistogramMatch;
	spokeResults[1].fnClass	= Components::ROUGH;
	spokeResults[1].weight = tuningWt(GpSpokeHistogramMatch);
	spokeResults[1].name = "hist";

	// The image match penalties go above here
#define LAST_IMAGE_MATCH_SPOKE_TERM	1

	spokeResults[2].tuningParm = SpokeMahalanobisMatch;
	spokeResults[2].fnClass	= Components::SMOOTH;
	spokeResults[2].weight = tuningWt(SpokeMahalanobisMatch);
	spokeResults[2].name = "mahal";

	spokeResults[3].tuningParm = GpSpokeGeometryWeight;
	spokeResults[3].fnClass	= Components::SMOOTH;
	spokeResults[3].weight = tuningWt(GpSpokeGeometryWeight);
	spokeResults[3].name = "dsq";

	// Note: spokeResults[4] is not an ordinary matchResult entry.
	// Its weight is passed to computeAtomPenalty(), which is used
	// compute the GpSpokeGeometryWeight, and spokeResults[4].useIt
	// is never set to true.
	spokeResults[4].tuningParm = GpSpokeNeighborWeight;
	spokeResults[4].fnClass	= Components::SMOOTH;
	spokeResults[4].weight = tuningWt(GpSpokeNeighborWeight);
	spokeResults[4].name = "neigh";

	spokeResults[5].tuningParm = GpSpokeInterWeight;
	spokeResults[5].fnClass	= Components::SMOOTH;
	spokeResults[5].weight = tuningWt(GpSpokeInterWeight);
	spokeResults[5].name = "cnstr";

#endif

	for (i = 0; i < MAX_NUM_FIG_MATCH_RESULTS; i++) {
		figureResults[i].value = 0.0;
		figureResults[i].useIt = false;
	}
	for (i = 0; i < MAX_NUM_ATOM_MATCH_RESULTS; i++) {
		atomResults[i].value = 0.0;
		atomResults[i].useIt = false;
	}

	// dibyendu

	for (i = 0; i < MAX_NUM_SREP_MATCH_RESULTS; i++) {
		sRepResults[i].value = 0.0;
		sRepResults[i].useIt = false;
	}

	for (i = 0; i < MAX_NUM_SPOKE_MATCH_RESULTS; i++) {
		spokeResults[i].value = 0.0;
		spokeResults[i].useIt = false;
	}

	targetSurfacePoints.numPoints = 0;
	targetSurfacePoints.points = 0;
	lmBpointIndices[0] = -1;
}

Match::~Match()
{
	int i;

	if (constrainedPointLists != NULL) {
		for (i = 0; i < numPointLists; i++) {
			if (constrainedPointLists[i] != NULL) {
				for (int j = 0; j < numPointLists; j++)
					constrainedPointLists[i][j].clear();        // Not strictly needed
			}
			delete [] constrainedPointLists[i];
		}
	}
	delete [] constrainedPointLists;

	for (i = 0; i < distanceMaps.size(); i++)
		delete distanceMaps[i];

	for (i = 0; i < masks.size(); i++) {
		if (masks[i] != NULL)
			delete masks[i];
	}
	masks.clear();			// Not strictly needed

	for (i = 0; i < dmasks.size(); i++) {
		if (dmasks[i] != NULL)
			delete dmasks[i];
	}
	dmasks.clear();			// Not strictly needed

	if (referenceObject != NULL)
		delete referenceObject;

	// Dibyendu
	if ( refObjectAfterScale != NULL )
		delete refObjectAfterScale ;

	for (i = 0; i < surfacePoints.size(); i++) {
		if (surfacePoints[i].points != NULL)
			delete[] surfacePoints[i].points;
	}

	for (i = 0; i < pListCount; i++) {
		if (pList[i] != NULL)
			delete pList[i];
	}
	delete [] pList;

	if (binaryDistanceMap != NULL)
		delete binaryDistanceMap;

#ifdef BINARY
	if (distanceVectorList != NULL)
		delete distanceVectorList;

	if (atomPredictor != NULL)
		delete atomPredictor;
#endif

	if (atomMahamPenalty != NULL)
		delete [] atomMahamPenalty;

	// spe doesn't belong to Match

	// cleaning up list of surface bpoints
	if (targetSurfacePoints.points) {
		delete [] targetSurfacePoints.points;
		targetSurfacePoints.points = 0;
		targetSurfacePoints.numPoints = 0;
	}

	if (contourMatch)
		delete contourMatch;
}

#ifndef BINARY
void Match::setWeight(int optimizerType, int tuningParm, double weight)
{
	OptimizerBase::optimizer_t type = (OptimizerBase::optimizer_t) optimizerType;
	if (type == OptimizerBase::MainFigure) {
		for (int i = 0; i < MAX_NUM_FIG_MATCH_RESULTS; i++)
			if (tuningParm == figureResults[i].tuningParm)
				figureResults[i].weight = weight;
	}
	else if (type == OptimizerBase::Deformation) {
		for (int i = 0; i < MAX_NUM_ATOM_MATCH_RESULTS; i++)
			if (tuningParm == atomResults[i].tuningParm)
				atomResults[i].weight = weight;
	}
	else if (type == OptimizerBase::Spoke) {
		for (int i = 0; i < MAX_NUM_SPOKE_MATCH_RESULTS; i++)
			if (tuningParm == spokeResults[i].tuningParm)
				spokeResults[i].weight = weight;
	}
	else
		cout << "Warning: ignoring attempt to set weights for optimizer type"
		<< optimizerType << endl;
}

double Match::getWeightedPenalty(int optimizerType, int tuningParm)
{
	OptimizerBase::optimizer_t type = (OptimizerBase::optimizer_t) optimizerType;
	if (type == OptimizerBase::MainFigure) {
		for (int i = 0; i < MAX_NUM_FIG_MATCH_RESULTS; i++)
			if (tuningParm == figureResults[i].tuningParm)
				return figureResults[i].value*figureResults[i].weight;
	}
	else if (type == OptimizerBase::Deformation) {
		for (int i = 0; i < MAX_NUM_ATOM_MATCH_RESULTS; i++)
			if (tuningParm == figureResults[i].tuningParm)
				return figureResults[i].value*figureResults[i].weight;
	}
	return -1.0;
}
#endif

void Match::setInitialTransform(const SimilarityTransform3D * initXform)
{
	if (initXform)
		initialTransform = *initXform;
	else	// P3DControl::objectSimTransform can be NULL
		initialTransform.setToIdentity();
}

void Match::definePLists(M3DObject * object, int subdivLevel)
{
	int i;
	Xferlist * xferList;
	M3DFigure *figure;

#ifdef DEBUG
	cout << "Match::definePLists()" << endl;
#endif

	if (pList != NULL) {
		for (i = 0; i < pListCount; i++) {
			if (pList[i] != NULL)
				delete pList[i];
		}
		delete [] pList;
	}

	int figureCount = object->getFigureCount();
	pListCount = figureCount + figureCount;

	pList = new ThallCode::Pointlist_server2*[pListCount];

	if (subdivLevel < 0)
		subdivLevel = surface_level_for_match;

	for (i = 0; i < figureCount; i++) {
		figure = object->getFigurePtr(i);
		// Rohit: Why should figure ever be NULL?
		assert( figure != NULL );
		//if (figure == NULL)
		//	continue;
		xferList = convertM3DtoXfer(figure);

		pList[i] = new ThallCode::Pointlist_server2;
		pList[i]->init(xferList);
		pList[i]->ComputeSubdivPointCloud(subdivLevel);

		pList[figureCount + i] = new ThallCode::Pointlist_server2;
		pList[figureCount + i]->init(xferList);
		pList[figureCount + i]->ComputeSubdivPointCloud(subdivLevel);

		delete [] (xferList->atomlist);
		delete xferList;
	}
}

#ifdef BINARY

bool Match::initialize(M3DObject * _referenceObject, Image3D * _targetImage,
					   ImageDistanceMap * _binaryDistanceMap, bool buildMasks, Image3D * trainingImage,
					   bool useWindowing, MatchType _type, int surfaceLevel, const char * templateFilename, 
					   const char * profileFilename)
{
	int figureCount, i;
	DMask * newDMask;

#ifdef DEBUG
	cout << "Match::initialize() called" << endl;
#endif

	windowIntensities = useWindowing;
	initialTransform = *(_referenceObject->getTransformation());
	resetReferenceObject(_referenceObject);
	targetImage = _targetImage;
	type = _type;
	//surface_level_for_match = surfaceLevel;

	for (i = 0; i < masks.size(); i++) {
		if (masks[i] != NULL)
			delete masks[i];
		if (dmasks[i] != NULL)
			delete dmasks[i];
	}

	masks.erase(masks.begin(), masks.end());
	dmasks.clear();

	if (type == BINARY_IMAGE_MATCH) {
		if (referenceObject->getFigureCount() <= tuningWt(BpFigureId)) {
			cout << "Binary image match was requested for figure " << tuningWt(BpFigureId) 
				<< " but the model has only " << referenceObject->getFigureCount() << " figures" << endl;
			return false;
		}
		if (!_binaryDistanceMap) {
			cout << "Warning: Binary image matching requires you to load or calc a distance map; remember to disable image match!" << endl;
			//return false;	// NO! image match might be disabled
		}

		binaryDistanceMap = _binaryDistanceMap;    // Take ownership	
		// Figure stage is responsible for method of moments initialization
	}
	else if (type == SIMPLE_MASK_MATCH) {
		simpleMask = new SimpleMaskFile;
		simpleMask->readBinary(profileFilename);
	}

	figureCount = referenceObject->getFigureCount();
	for (i = 0; i < figureCount; i++) {
		masks.push_back(NULL);
		if (buildMasks)
			buildMask(i, type, trainingImage, surfaceLevel, templateFilename,
			profileFilename);

		newDMask = new DMask(referenceObject, i, MASK_CUTOFF, type, surfaceLevel, 
			templateFilename, profileFilename);
		dmasks.push_back(newDMask);
	}

	return true;
}

bool Match::initContourMatch(const char * filename)
{
	// Initialize contour match
	if (contourMatch) {
		delete contourMatch;
		contourMatch = NULL;
	}

	if (tuningWt(FigureContourMatch) == 0.0)
		return false;

	if (filename) {
		if (globalVerbosity >= 1) {
			if (targetImage)
				cout << "Setting image for contourMatch." << endl;
			else
				cout << "No image for contourMatch." << endl;
		}
		contourMatch = new DistanceToPointSetFunction(filename,
			DistanceToPointSetFunction::Full,
			(int) tuningWt(FigureContourLevel),
			(int) tuningWt(FigureContourFigureId), targetImage,
			(int) tuningWt(BpFigureContourPadding)); 
		contourMatch->setClosenessLimit(tuningWt(FigureContourCloseness));
		contourMatch->setSizeWeight(tuningWt(BpFigureContourSizeWeight));
	}
	return true;
}

#else	/* BINARY */

bool Match::initContourMatch(const char * filename)
{
	Image3D * image;

	// Initialize contour match
	if (contourMatch) {
		delete contourMatch;
		contourMatch = NULL;
	}
	if (! filename)
		return false;

	if (tuningWt(FigureContourMatch) == 0.0)
		return false;

	if (0 == (int) tuningWt(GpFigureContourCoords))
		image = NULL;
	else
		image = targetImage;

	if (globalVerbosity >= 1) {
		cout << "Using contour penalty for figure "
			<< (int) tuningWt(FigureContourFigureId) << '\n';
		cout << "Contour points are in ";
		if (image)
			cout << "world";
		else
			cout << "model";
		cout << " coordinates" << endl;
	}

	contourMatch = new DistanceToPointSetFunction(filename,
		DistanceToPointSetFunction::Full,
		(int) tuningWt(FigureContourLevel),
		(int) tuningWt(FigureContourFigureId), image);
	contourMatch->setClosenessLimit(tuningWt(FigureContourCloseness));
	return true;
}

bool Match::initialize(M3DObject * _referenceObject, Image3D * _targetImage,
					   bool buildMasks, Image3D * trainingImage, bool useWindowing,
					   MatchType _type, int surfaceLevel, const char * templateFilename, 
					   const char * profileFilename)
{
	int figureCount,
		i;
	DMask * newDMask;

#ifdef DEBUG
	cout << "Match::initialize() called" << endl;
#endif

	windowIntensities = useWindowing;
	resetReferenceObject(_referenceObject);

	targetImage = _targetImage;
	type = _type;
	surface_level_for_match = surfaceLevel;

	if (_type == HISTO_IMAGE_MATCH)
		return true;

	if (_type == MULTIPATCH_MATCH)
	{
		mp->initializePatches(_referenceObject, _targetImage);
		return true;
	}

	for (i = 0; i < masks.size(); i++) {
		if(masks[i] != NULL)
			delete masks[i];
	}
	masks.clear();
	for (i = 0; i < dmasks.size(); i++) {
		if (dmasks[i] != NULL)
			delete dmasks[i];
	}
	dmasks.clear();

	if (type == BINARY_IMAGE_MATCH) {

		if (_referenceObject->getFigureCount() != 1) {
			cout << "Binary image matching only works on a single figure" << endl;
			return false;
		}

		// Create a distance map from the image
		binaryDistanceMap = new ImageDistanceMap(targetImage);
		(void) binaryDistanceMap->createMap();
		const int * map_size = binaryDistanceMap->getMapSize();
		if (map_size[0] <= 0) {
			cout << "Image distance map_size of 0!" << endl;
			return false;
		}
		if (binaryDistanceMap->status() == false) {
			cout << "Warning: Binary image matching requires you to load or calc a distance map; remember to disable image match!" << endl;
			return false;
		}

		/*		{	// Create snapshot image of the distance map extrapolated over the full input image
		GreyValue * img = new GreyValue[map_size[0]*map_size[1]*map_size[2]];
		for (int k = 0; k < map_size[2]; k++)
		for (int j = 0; j < map_size[1]; j++)
		for (i = 0; i < map_size[0]; i++) {
		float d = binaryDistanceMap->getWorldDistance(i, j, k);
		img[i + map_size[0]*(j + map_size[1]*k)] = (GreyValue) (0.5 + d);
		}
		Image3D im;
		im.setVoxels(img, map_size[0], map_size[1], map_size[2]);
		im.setSpacing(1.0, 1.0, 1.0);
		snapshot("../../../data/images/extrap_map.raw3", im);
		}
		*/
	}
	else if (type == SIMPLE_MASK_MATCH) {
		simpleMask = new SimpleMaskFile;
		if (profileFilename == NULL)
			return false;
		simpleMask->readBinary(profileFilename);
	}

	figureCount = referenceObject->getFigureCount();
	for (i = 0; i < figureCount; i++) {
		masks.push_back(NULL);
		if (buildMasks)
			if (! buildMask(i, type, trainingImage, surfaceLevel, templateFilename,
				profileFilename))
			{
				cout << "Mask constructor failed" << endl;
				return false;
			}

			newDMask = new DMask(referenceObject, i, MASK_CUTOFF, type, surfaceLevel, 
				templateFilename, profileFilename);
			dmasks.push_back(newDMask);
	}
	initializeAtomStage();
	return true;
}

#endif	/* BINARY */

// This must be used after initialize() is called and before each atom stage
void Match::initializeAtomStage()
{
	int figureCount, i, nAtoms;
	M3DQuadFigure * figure;

	figureCount = referenceObject->getFigureCount();
	nAtoms = 0;
	for (i = 0; i < figureCount; i++) {
		figure = (M3DQuadFigure *) referenceObject->getFigurePtr(i);
		int n = figure->getPrimitiveCount();
		if (n > nAtoms)
			nAtoms = n;
	}

	if (atomMahamPenalty != NULL && nAtoms > nAtomsMax) {
		delete [] atomMahamPenalty;
		atomMahamPenalty = NULL;
	}
	if (atomMahamPenalty == NULL) {
		// Allocate one time for all atom stages
		nAtomsMax = nAtoms;
		atomMahamPenalty = new double[nAtomsMax];
	}

	// Initialize for later recording the maham distance for each atom
	memset(atomMahamPenalty, 0, nAtomsMax*sizeof(double));

	lmBpointIndices[0] = -1;
}

/*  Resets the reference object; penalties are computed relative to it
(this should be called at the start of a new optimizer stage).
Atoms to be used in the optimization stages must be selected in
the reference object before this function is called.  This is also
true if this function is called from one of the initialize() functions
above.

This function builds the surfacePoints object, which contains one
SurfacePointList entry per figure.  The surface points are computed
using a transformed copy of the mrep model.
*/
void Match::resetReferenceObject(M3DObject * _referenceObject)
{
	M3DFigure * figure;
	Xferlist * xferList;
	SurfacePointList surfacePointList;    
	int i, figureCount;
	M3DObject * tempObject;


#ifdef DEBUG
	cout << "Match resetReferenceObject called" << endl;
#endif
	if (_referenceObject == NULL)
		return;

	if (referenceObject != NULL)
		delete referenceObject;

	for (i = 0; i < surfacePoints.size(); i++) {
		if(surfacePoints[i].points != NULL)
			delete[] surfacePoints[i].points;
	}
	surfacePoints.clear();

	referenceObject = _referenceObject->assign();
	tempObject = referenceObject->assign();
	tempObject->select();	// For applySimilarityAboutCOG()

	// Transform object by initial transform
#ifdef UNFLIPPED
	tempObject->applySimilarity(initialTransform, Vector3D(0.0, 1.0, 0.0));
#else
	tempObject->applySimilarity(initialTransform, Vector3D(0.0, 0.0, 0.0));
#endif

	figureCount = referenceObject->getFigureCount();
	for (i = 0; i < figureCount; i++) {
		figure = tempObject->getFigurePtr(i);

		if (figure == NULL)
			continue;

		xferList = convertM3DtoXfer(figure);

#ifdef BINARY
		int surface_level_for_match = (int) tuningWt(BpAtomSurfaceLevel);					// AGG: This probably should be put in the constructor
		// AGG: There is an error here: the above definition hides Match::surface_level_for_match
#endif
		pList[i]->UpdateSubdivPointCloud(surface_level_for_match, xferList);

		pList[i]->copySubdivBoundaryInfo(&(surfacePointList.numPoints),
			&(surfacePointList.points));		
		surfacePoints.insert(surfacePoints.end(), surfacePointList);

		delete [] (xferList->atomlist);
		delete xferList;
	}
	delete tempObject;

	initializeAtomStage();
}

bool Match::buildMask(int figureId, MatchType _type, Image3D * trainingImage, 
					  int surfaceLevel, const char * templateFilename, const char * profileFilename)
{
	if (figureId < 0 || figureId >= masks.size())
		return false;

	if (masks[figureId] != NULL)
		delete masks[figureId];

	// Types NOTCH_MATCH, PROFILE_MATCH, TRAINING_IMAGE_WITH_INTERIOR_MATCH.
	// and TRAINING_IMAGE_WITHOUT_INTERIOR_MATCH were deleted June 16, 2006.
	switch(_type)
	{
	case GAUSSIAN_DERIVATIVE_MATCH:
		masks[figureId] = new Mask(referenceObject, figureId, MASK_CUTOFF,
			GAUSSIAN_DERIVATIVE_MATCH, NULL, windowIntensities,
			surfaceLevel, templateFilename);
		break;

	case NEG_GAUSSIAN_DERIVATIVE_MATCH:
		masks[figureId] = new Mask(referenceObject, figureId, MASK_CUTOFF,
			NEG_GAUSSIAN_DERIVATIVE_MATCH, NULL, windowIntensities,
			surfaceLevel, templateFilename);
		break;

	case SIMPLE_MASK_MATCH:
		masks[figureId] = new Mask(referenceObject, figureId, MASK_CUTOFF,
			SIMPLE_MASK_MATCH, NULL, windowIntensities, surfaceLevel,
			templateFilename, profileFilename, simpleMask);
		break;

	default:
		masks[figureId] = NULL;
		return false;
	};
	return masks[figureId]->status();
}

void Match::setConstrainedPointLists(std::vector<referencePoint> **
									 _constrainedPointLists)
{
	int i;

	if (referenceObject == NULL)
		return;

	if (constrainedPointLists != NULL) {
		for (i = 0; i < numPointLists; i++) {
			if (constrainedPointLists[i] != NULL)
				delete [] constrainedPointLists[i];
		}
	}
	delete [] constrainedPointLists;

	constrainedPointLists = _constrainedPointLists;
	numPointLists = referenceObject->getFigureCount();
}

#ifndef BINARY

// A function to calculate the new image match given a new target object
double Match::computeModelMatch(M3DObject * targetObject)
{
	int i, figureCount;
	Vector3D coord, newCoord;
	double figuralMatchValue;
	double matchValue;

	if (targetImage == NULL) {
		cout << "Error: No target image\n" << endl;
		return FAILED;
	}

	// The match of each object(figure) will contribute to the overall image
	// match by the same ratio they take in the sense of volume, but
	// still, this is to be proved to be better than the simple average
	figureCount = masks.size();
	double * maskVolumes = new double[figureCount];
	double * maskVolumeRatios = new double [figureCount];
	double maskVolumeTotal = 0.0;
	Mask * matchMask;

	for (i = 0; i < figureCount; i++) {
		matchMask = masks[i];
		if (matchMask == NULL) {
			maskVolumes[i] = 0;
			continue;
		}
		maskVolumes[i] = matchMask->getTotalVolume();
		maskVolumeTotal += maskVolumes[i];
	}

	for (i = 0; i < figureCount; i++)
		maskVolumeRatios[i] = maskVolumes[i]/maskVolumeTotal;

	matchValue = 0.0;
	for (i = 0; i < figureCount; i++) {
		// If the target image is stacked, we want to set the mask here,
		// according to which image we want.  Since no other information
		// seems to be available here right now, we make the simplifying
		// assumptions that all objects are single figures, and thus the
		// the figure index (or id I guess) tells which image to use...

		matchMask = masks[i];

		if (targetImage->getIsImageStacked())
			targetImage->setStackedMask((GreyValue) matchMask->getStackedImageMask());

		// CL: assume each object is single figure
		// TO DO: needs to loop through figure tree to deal with
		// multi-object, multi-figure models
		if (i == 0 || isAdaptive())
			figuralMatchValue = computeMainFigureMatch(targetObject, i);
		else
			figuralMatchValue = computeSubfigureMatch(targetObject, i);

		matchValue += maskVolumeRatios[i] * figuralMatchValue;
	}

#ifdef DEBUG
	printf("Match Value: %f\n", matchValue);
#endif

	delete [] maskVolumes;
	delete [] maskVolumeRatios;

	return -matchValue;
}

double Match::computeModelMatch(const SimilarityTransform3D & newTransformation)
{
	double x, y, z;
	int i, j, figureCount;
	SimilarityTransform3D transformation;
	MaskElement * maskElement;
	Vector3D coord, newCoord;
	double val1, val2;
	double figuralMatchValue;
	double matchValue;
	double mu, totval2Sqr, rms;
	double maskVolume;

	if(targetImage == NULL)
		return WORST_MATCH;

	transformation = initialTransform;
	transformation.multiplyBy(newTransformation);

	matchValue = 0.0;

	figureCount = masks.size();
	for (i = 0; i < figureCount; i++) {
		Mask * matchMask = masks[i];
		if (matchMask == NULL)
			continue;

		// If the target image is stacked, we want to set the mask here,
		// according to which image we want.  Since no other information
		// seems to be available here right now, we make the simplifying
		// assumptions that all objects are single figures, and thus the
		// the figure index (or id I guess) tells which image to use...

		if (targetImage->getIsImageStacked())
			targetImage->setStackedMask((GreyValue) matchMask->getStackedImageMask());

		maskVolume = matchMask->getTotalVolume();

		figuralMatchValue = 0.0;

		mu = 0.0;
		totval2Sqr = 0.0;

		for (j = 0; j < matchMask->getSize(); j++) {
			maskElement = matchMask->getElement(j);

			if (maskElement == NULL) {
				printf("NULL mask element.\n");
				continue;
			}
			if (maskElement->inBlendRegion)
				continue;

			x = maskElement->x;
			y = maskElement->y;
			z = maskElement->z;

			newCoord.set(x, y, z);

			transformation.transformVector(newCoord);

			targetImage->modelToImageCoordinates(newCoord);

			if (windowIntensities)
				val2 = targetImage->mapDisplayToRelative(
				targetImage->getWindowedInterpolatedVoxelValue(
				newCoord.getX(), newCoord.getY(), newCoord.getZ()));
			else
				val2 = targetImage->mapDisplayToRelative(
				targetImage->getInterpolatedVoxelValue(
				newCoord.getX(), newCoord.getY(), newCoord.getZ()));

			mu += val2 * maskElement->volelem;
			totval2Sqr += val2 * val2 * maskElement->volelem;
		}

		mu /= maskVolume;
		rms = (totval2Sqr / maskVolume) - (mu * mu);
		if (rms <= R_SMALL_TOLERANCE)
			rms = 0.0;
		else
			rms = sqrt(rms);

		if (rms <= R_SMALL_TOLERANCE)
			continue;

		for (j = 0; j < matchMask->getSize(); j++) {
			maskElement = matchMask->getElement(j);
			if (maskElement == NULL) {
				printf("NULL mask element\n");
				continue;
			}
			if (maskElement->inBlendRegion)
				continue;

			x = maskElement->x;
			y = maskElement->y;
			z = maskElement->z;

			newCoord.set(x, y, z);

			transformation.transformVector(newCoord);

			targetImage->modelToImageCoordinates(newCoord);

			val1 = maskElement->tempVal;
			if (windowIntensities)
				val2 = targetImage->mapDisplayToRelative(
				targetImage->getWindowedInterpolatedVoxelValue(
				newCoord.getX(), newCoord.getY(), newCoord.getZ()));
			else
				val2 = targetImage->mapDisplayToRelative(
				targetImage->getInterpolatedVoxelValue(
				newCoord.getX(), newCoord.getY(), newCoord.getZ()));

			figuralMatchValue += val1*((val2 - mu) / rms)*maskElement->volelem;
		}

		matchValue += (figuralMatchValue / maskVolume);
	}

	if (figureCount > 0)
		matchValue /= ((double) figureCount);

#ifdef DEBUG
	printf("Match Value mod: %f\n", matchValue);
#endif
	return matchValue;
}

#else	/* BINARY */


// RET: unit cube distance^2 summed over all landmarks.
//   Each Euclidian distance is from (first atom of
//   figure named "LM-*") to (tip of b-vector of
//   same-named atom in targetmodeObject).
// Larger values mean a worse fit.
// Values [0..1] mean all landmarks are reasonably well fit.
// [1..] means worse fits.
double Match::computeLandmarkPenalty(M3DObject * targetObject, int figureId,
									 M3DObject * liObject, const double dilationFactorInModelUnits, 
									 int subdivLevel)
{
	double penalty = 0;
	if (!targetObject || ! liObject)
		return 0;

	static int complainedAlready = 0;	// don't be too chatty next call
	int complaining = 0;				// was chatty this call

	// foreach landmark figure ilm,
	//  - mlm = find model-space position of corresponding landmark in targetObject
	//  - pen += ((ilm - mlm) / mlm.tolerance) ^ 2;

	M3DFigure * targetFigure;

	if (figureId < 0 || targetObject == NULL || liObject == NULL)
		return 0;		// nothing to match

	targetFigure = targetObject->getFigurePtr(figureId);

	if (!targetSurfacePoints.points) {

		Xferlist * xferList = convertM3DtoXfer(targetFigure);

		pList[figureId]->UpdateSubdivPointCloud(subdivLevel, xferList);
		pList[figureId]->copySubdivBoundaryInfo(&targetSurfacePoints.numPoints, &targetSurfacePoints.points);

		delete [] (xferList->atomlist);
		delete xferList;

	}




	// "li"-style values are the position of the landmark in the image, which is fixed.
	// "lm" are the position of the model landmark, which is moved during optimization.

	// Loop over model lm
	// Each model lm is associated with a Bpoint that comes from seurat.
	// The indices into that list of Bpoints are stored in lmBpointIndices
	// so that the linear search (to match a bpoint to a landmark)
	// only gets performed once, rather than for every model
	int nlm = targetFigure->getLandmarkCount();
	int nilm = liObject->getFigureCount();

	int foundNumLM = 0;

	bool needToComputeIdx = (lmBpointIndices[0] == -1);
	for (int lm = 0; lm < nlm; lm++) {
		char * lmName = targetFigure->getLandmarkName(lm);

		// HACK!! negating the LM T value to match the T's in Bpoint
		// which seem to be inverted from the ones in M3DPrimitiveRenderer
		double lmT = -targetFigure->getLandmarkAtomT(lm);

		double lmU = targetFigure->getLandmarkAtomU(lm);
		double lmV = targetFigure->getLandmarkAtomV(lm);
		if (needToComputeIdx) {
			//need to compute u,v from lmIndex
			// to be in sync with what happens in quad mesh,
			// u column, v is row
			lmBpointIndices[lm] = -1;

			/*  BUG: This is too slow... a faster method to get the bpoint is 
			Pointlist_server2::subdivBPosition */
			Bpoint bpnt; 
			lmBpointIndices[lm] = pList[figureId]->subdivNearestBPosition(&bpnt, lmU, lmV, lmT );
			if (lmBpointIndices[lm] != -1)
				cout << "Info: bpoint for landmark:  " << lmName << " is: " << lmBpointIndices[lm] << endl;
			else
				cout<<"Info: bpoint for landmark:  " << lmName << " not found!" << endl;
			/*for (int idx = 0; idx < targetSurfacePoints.numPoints; idx++) {
			if ( 
			(targetSurfacePoints.points[idx].u == lmU) && 
			(targetSurfacePoints.points[idx].v == lmV) &&
			(targetSurfacePoints.points[idx].t == lmT) 
			) {

			lmBpointIndices[lm] = idx;
			cout << "Info: bpoint for landmark:  " << lmName << " is: " << lmBpointIndices[lm] << endl;
			break;
			}
			}*/
		}

		// li has 1 landmark per figure (use only the 1st atom of each figure);
		// find that figure by matching it's name against model lm's name
		int gotOne = 0;
		M3DFigure * liFigure;
		for (int li = 0; li < nilm && !gotOne; li++) {
			liFigure = liObject->getFigurePtr(li);
			// HACK: skip optional "LM-" part of name
			int offset = 0;
			if (! strncmp(liFigure->getName(), "LM-", 3))
				offset = 3;
			if (! strcmp(liFigure->getName() + offset, lmName))
				gotOne = 1;
		}

		// Warn user that the corresponding image lm was not found
		if (! gotOne) {
			complaining = 1;
			if (! complainedAlready) {
				cout << "Warning: cannot find image landmark for model landmark "
					<< lmName << " -- ignoring landmark" << endl;
			}
			continue;		// Skip to next model lm
		}

		M3DPrimitive * liAtom = liFigure->getPrimitivePtr(0);
		Vector3D liPosition = liAtom->getX();
		double tolerance = liAtom->getR();


		if (lmBpointIndices[lm] == -1) {
			cout << "Error: could not find bpoint for landmark:  " << lmName << endl;
			continue;
		}
		else {
			foundNumLM++;
		}

		/*
		b allows us to put the landmark anywhere on the subdivision surface.
		the correction by dilationFactor * the unit surface normal
		is to support thin objects (where the distance map is actually measuring
		distance to the boundary of a dilated version of the model)
		*/
		Bpoint b = targetSurfacePoints.points[lmBpointIndices[lm]];
		Vector3D lmPosition(b.pnt);
		if (dilationFactorInModelUnits != 0.0) {
			lmPosition -= (dilationFactorInModelUnits * Vector3D(b.norm));
		}   

		// Update penalty
		double weightedDist = sqrt(liPosition.distSquare(lmPosition)) / tolerance;
		double weightedDistSquared = weightedDist * weightedDist;
		if (globalVerbosity >= 1)
			cout << " atom [col: " << lmU << ", row: " << lmV << "] (" << lmName << ")^2 = " << weightedDistSquared;
		penalty += weightedDistSquared;
	}

	if (complaining) complainedAlready = 1;	// silent on next call

	if (foundNumLM) penalty /= foundNumLM;		// Average over all landmarks

	if (globalVerbosity >= 1)
		cout << "\n     lmAvg = " << penalty << endl;
	return penalty;

}

void Match::setUpAtomPredictor(int figureId)
{
	if (referenceObject == NULL)
		return;

	M3DFigure * figure = referenceObject->getFigurePtr(figureId);
	// FIXME: tubes: need setInitialTubeFigure and supporting code -GST/RSS
	delete atomPredictor;
	if (typeid(*figure) == typeid(M3DQuadFigure)) {
		M3DAtomPredictorQuad* atomPredictor = new M3DAtomPredictorQuad;
		atomPredictor->setInitialFigures(referenceObject->getFigurePtr(figureId), 
			referenceObject->getFigurePtr(figureId));
		this->atomPredictor	= atomPredictor;
	}
	else if(typeid(*figure) == typeid(M3DTubeFigure)) {
		M3DAtomPredictorTube* atomPredictor	= new M3DAtomPredictorTube;
		this->atomPredictor	= atomPredictor;
	}
}

#endif	/* BINARY */

Components Match::computeAtomMatchAndPenaltiesAsComponents(M3DObject * targetObject, int figureId,
														   int atomId, const Vector& x, double & geomPenalty)
{
#ifdef BINARY
	DistanceType distType;
#endif
    double penalty = 0;
	Components totalPenalty;
	totalPenalty.smooth		= 0.0;
	totalPenalty.rough		= 0.0;
	totalPenalty.smoothSet	= false;
	totalPenalty.roughSet	= false;
	int nTerms, i;

	nTerms = 0;

	// cleaning up list of surface bpoints
	if (targetSurfacePoints.points) {
		delete [] targetSurfacePoints.points;
		targetSurfacePoints.points = 0;
		targetSurfacePoints.numPoints = 0;
	}

	if (type == BINARY_IMAGE_MATCH) {
#ifdef BINARY
		geomPenalty	= 0.0;
		if (globalVerbosity >= 1)
			cout << "   ";

		distType = (enum DistanceType) (int) tuningWt(BpAtomDistanceType);
		// Convert units for dilationFactor.
		// (Distance map needs to have uniform imageToModelScale)
		// FIXME: Check for this condition.
		double dilationFactorInModelUnits = 0;
		if(tuningWt(BpDilationFactor) > 0) {
			const Vector3D imageToModelScale(
				1.0/binaryDistanceMap->getModelToImageScale()[0],
				1.0/binaryDistanceMap->getModelToImageScale()[1],
				1.0/binaryDistanceMap->getModelToImageScale()[2] );
			dilationFactorInModelUnits = tuningWt(BpDilationFactor) * imageToModelScale.getX();
		}
		if (tuningWt(BpAtomAverageNeighborMatch) != 0.0) {
			PrimNeighborhoodDefn neiType = (enum PrimNeighborhoodDefn)
				(int) tuningWt(BpAtomNeighborhood);
			penalty = targetObject->dist2FromAveOfNeighbors(neiType,
				(int) tuningWt(BpFigureId), distType);	//MultiObject
			atomResults[0].value = penalty;
			atomResults[0].useIt = true;
			geomPenalty	+= tuningWt(atomResults[0].tuningParm) * penalty;
		}
		if (tuningWt(BpAtomImageMatch) != 0.0) {
			penalty = computeBinaryImageMatch(targetObject, figureId,
				int(tuningWt(BpDDMEnable)), tuningWt(BpDDMAngleTolerance),
				int(tuningWt(BpAtomSurfaceLevel)));
			atomResults[1].value = penalty;
			atomResults[1].useIt = true;
		}
		if (tuningWt(AtomMahalanobisMatch) != 0.0) {
			penalty = computeAtomMahamPenalties(x, figureId, atomId);
			atomResults[2].value = penalty;
			atomResults[2].useIt = true;
			geomPenalty	+= tuningWt(atomResults[2].tuningParm) * penalty;
		}
		if (tuningWt(BpAtomLandmarkMatch) != 0.0) {
			penalty = computeLandmarkPenalty(targetObject, figureId, landmarkObject, 
				dilationFactorInModelUnits, int(tuningWt(BpAtomSurfaceLevel)));
			atomResults[3].value = penalty;
			atomResults[3].useIt = true;
			geomPenalty	+= tuningWt(atomResults[3].tuningParm) * penalty;
		}
		if (tuningWt(BpAtomModelMatch) != 0.0) {
			penalty = targetObject->dist2FromObject(targetObject->loadedObject(),
				(int) tuningWt(BpFigureId), distType);
			atomResults[4].value = penalty;
			atomResults[4].useIt = true;
			geomPenalty	+= tuningWt(atomResults[4].tuningParm) * penalty;
		}
		if (tuningWt(BpAtomInterpenetrationMatch) != 0.0) {
			penalty = computeInterpenetrationPenalty(targetObject, figureId);
			atomResults[5].value = penalty;
			atomResults[5].useIt = true;
			geomPenalty	+= tuningWt(atomResults[5].tuningParm) * penalty;
		}
		if (tuningWt(BpAtomRSradPenalty) != 0.0) {
			M3DFigure * figure = targetObject->getFigurePtr(figureId);
            if (tuningWt(BpAtomSumMatch) != 0.0){
				// Sum over all atoms.
                if(atomPredictor){
                    penalty = atomPredictor->getFigureRSradPenalty(figure, -1, PEL_BOTH, tuningWt(RSradPenaltyNorm), tuningWt(RSradPenaltyThreshold) );
                }else{
                    cout<<"atomPredictor is Null: Penalty value set to 0;"<<endl;
                }
            }else{
                if(atomPredictor){
                    penalty = atomPredictor->getFigureRSradPenalty(figure, atomId, PEL_BOTH, tuningWt(RSradPenaltyNorm), tuningWt(RSradPenaltyThreshold) );
                }else{
                    cout<<"atomPredictor is Null: Penalty value set to 0;"<<endl;
                }
            }
			atomResults[6].value = penalty;
			atomResults[6].useIt = true;
			geomPenalty	+= tuningWt(atomResults[6].tuningParm) * penalty;
		} 
		if (tuningWt(BpAtomDeltaNeighborMatch) != 0.0) {
			M3DQuadFigure * figure = dynamic_cast<M3DQuadFigure*>( targetObject->getFigurePtr(figureId));
			if (tuningWt(BpAtomSumMatch) != 0.0) {
				// Sum over all atoms
				for (int i = 0; i < figure->getPrimitiveCount(); i++)
					penalty = dynamic_cast<M3DAtomPredictorQuad*>(atomPredictor)->getAtomPenaltyWithoutOpt(figure, i);
				penalty /= figure->getPrimitiveCount();
			}
			else 
				penalty = dynamic_cast<M3DAtomPredictorQuad*>(atomPredictor)->getAtomPenaltyWithoutOpt(figure, atomId);
			atomResults[7].value = penalty;
			atomResults[7].useIt = true;
			geomPenalty	+= tuningWt(atomResults[7].tuningParm) * penalty;
		}
		if ( tuningWt(BpAtomCurvinessPenalty) != 0.0) {
			penalty	= targetObject->getFigurePtr(figureId)->curviness();
			atomResults[8].value = penalty;
			atomResults[8].useIt = true;
			geomPenalty	+= tuningWt(atomResults[8].tuningParm) * penalty;
		}

		//Xiaojie
		if (tuningWt(BpAtomImageNormalMatch) != 0.0) {
            penalty = computeAtomBinaryImageNormalMatch(targetObject, figureId, atomId, int(tuningWt(BpAtomSurfaceLevel)));
			atomResults[9].value = penalty;
			atomResults[9].useIt = true;
		}
		if (tuningWt(BpAtomImagePlaneOrienAtEnd) != 0.0) {
			penalty = computeAtomBinaryImagePlaneOrienMatchAtEnd(targetObject, figureId, atomId);
			atomResults[10].value = penalty;
			atomResults[10].useIt = true;
		}
		if (tuningWt(BpAtomImageVertexMatch) != 0.0) {
			penalty = computeAtomBinaryImageVertexMatch(targetObject, figureId, atomId);
			atomResults[11].value = penalty;
			atomResults[11].useIt = true;
		}

		for (i = 0; i < MAX_NUM_ATOM_MATCH_RESULTS; i++) {
			double wtMatch;

			if (atomResults[i].useIt) {
				atomResults[i].useIt = false;
				wtMatch = tuningWt(atomResults[i].tuningParm)*atomResults[i].value;
				switch(atomResults[i].fnClass) {
					case Components::ROUGH:
						totalPenalty.rough		+= wtMatch;
						totalPenalty.roughSet	= true;
						break;
					case Components::SMOOTH:
						totalPenalty.smooth		+= wtMatch;
						totalPenalty.smoothSet	= true;
						break;
					default:
						assert(false);
				}
				if (globalVerbosity >= 1) {
					cout << tuningShortName(atomResults[i].tuningParm) << ':'
						<< tuningWt(atomResults[i].tuningParm) << '*'
						<< atomResults[i].value << ' ';
				}
			}
		}

		if (globalVerbosity >= 1)
			cout << " = " << totalPenalty << endl;
#else
		cout << "Binary image match is not implemented\n";
		return Components(FAILED);
#endif
	}
#ifdef BINARY
	else {
		cout << "Error: computeAtomMatchAndPenalties failed\n";
		return Components(FAILED);
	}
#else	/* BINARY */
	else {	// Greyscale Pablo
		if (type != HISTO_IMAGE_MATCH && type != MULTIPATCH_MATCH) {
			// Mask-based image match
#ifdef NO_MATCH_PENALTY
			penalty = 0.0;
#else
			penalty = computeMainFigureMatch(targetObject, figureId);
#endif
			atomResults[0].value = penalty;
			atomResults[0].useIt = true;
		}
		else if (type == HISTO_IMAGE_MATCH) {
			// Histogram image match
			if (atomResults[1].weight != 0.0) {
				if (figureId < 0 || targetImage == NULL || targetObject == NULL) {
					cout << "Warning: worst histogram match being returned\n";
					return Components(FAILED);
				}
				if (spe == NULL) {
					cout << "Error: no spe\n";
					return Components(FAILED);
				}
				if (! spe->update(targetImage, targetObject, figureId, -1))
				{
					//	cout << "Error: surface patch (histogram) image match failed\n";
					return Components(FAILED);
				}				
				penalty = spe->getMatchValue(figureId, -1);


				//printf("Histogram match value %.4lf\n", penalty);
				atomResults[1].value = penalty;
				atomResults[1].useIt = true;

			}
		}
		else if (type == MULTIPATCH_MATCH)
		{
			if (atomResults[2].weight != 0.0) {
				if (figureId < 0 || targetImage == NULL || targetObject == NULL) {
					cout << "Warning: worst MutliPatch match being returned\n";
					return Components(FAILED);
				}
				if (mp == NULL) {
					cout << "Error: no MultiPatch mp \n";
					return Components(FAILED);
				}

				penalty = mp->getYourMatch(targetObject, targetImage, figureId, NULL);

				//printf("MutliPatch match value %.4lf\n", penalty);
				//I'm piggybacking on Eli's histo match, cause it looks like a hassle to put in my own.
				atomResults[2].value = penalty;
				atomResults[2].useIt = true;

			}
		}
		else
		{
			cout << "ERROR: not sure what match you want...";
			return Components(FAILED);
		}

		if (atomResults[3].weight != 0.0) {
			// Mahalanobis penalty
			penalty = computeAtomMahamPenalties(x, figureId, atomId);
			atomResults[3].value = penalty;
			atomResults[3].useIt = true;
		}
		else {
			// Geometric penalty
			if (atomResults[4].weight != 0.0) {
				penalty = computeAtomPenalty(*targetObject, figureId,
					atomId, atomResults[4].weight);
				atomResults[4].value = penalty;
				atomResults[4].useIt = true;
			}
		}

		// Constraints penalty
		if (atomResults[6].weight > 0.0) {
			penalty = computeConstraintsPenalty(*targetObject, figureId);
			atomResults[6].value = penalty;
			atomResults[6].useIt = true;
		}

		// Now collect the terms
		if (globalVerbosity >= 2)
			cout << "  ";
		geomPenalty = 0.0;
		for (i = 0; i < MAX_NUM_ATOM_MATCH_RESULTS; i++) {
			double wtMatch;

			if (atomResults[i].useIt) {
#ifndef OPTIMIZATION_VISUALIZER
				atomResults[i].useIt = false;
#endif
				wtMatch = atomResults[i].value*atomResults[i].weight;
				switch(atomResults[i].fnClass) {
					case Components::ROUGH:
						totalPenalty.rough		+= wtMatch;
						totalPenalty.roughSet	= true;
						break;
					case Components::SMOOTH:
						totalPenalty.smooth		+= wtMatch;
						totalPenalty.smoothSet	= true;
						break;
					default:
						assert(false);
				}
				// The image match entries do not contribute to geomPenalty
				if (i > LAST_IMAGE_MATCH_ATOM_TERM)
					geomPenalty += wtMatch;

				if (globalVerbosity >= 1) {
					cout << atomResults[i].name << ':' << atomResults[i].weight
						<< '*' << atomResults[i].value << ' ';
				}
			}
		}

		if (globalVerbosity >= 1)
			cout << "= " << totalPenalty << endl;
	}
#endif	/* BINARY */

	return totalPenalty;
}


Components Match::computeSRepMatchAndPenaltiesAsComponents(M3DObject * targetObject, int figureId,
														   int atomId, int spokeId, const Vector& x, double & geomPenalty)
{
#ifdef BINARY
	DistanceType distType;
#endif
	double penalty;
	Components totalPenalty;
	totalPenalty.smooth		= 0.0;
	totalPenalty.rough		= 0.0;
	totalPenalty.smoothSet	= false;
	totalPenalty.roughSet	= false;
	int nTerms, i;

	nTerms = 0;

	// cleaning up list of surface bpoints
	if (targetSurfacePoints.points) {
		delete [] targetSurfacePoints.points;
		targetSurfacePoints.points = 0;
		targetSurfacePoints.numPoints = 0;
	}

	if (type == BINARY_IMAGE_MATCH) {
#ifdef BINARY
		geomPenalty	= 0.0;
		if (globalVerbosity >= 1)
			cout << "   ";

		distType = (enum DistanceType) (int) tuningWt(BpSRepDistanceType);
		// Convert units for dilationFactor.
		// (Distance map needs to have uniform imageToModelScale)
		// FIXME: Check for this condition.
		double dilationFactorInModelUnits = 0;
		if(tuningWt(BpDilationFactor) > 0) {
			const Vector3D imageToModelScale(
				1.0/binaryDistanceMap->getModelToImageScale()[0],
				1.0/binaryDistanceMap->getModelToImageScale()[1],
				1.0/binaryDistanceMap->getModelToImageScale()[2] );
			dilationFactorInModelUnits = tuningWt(BpDilationFactor) * imageToModelScale.getX();
		}


		// No need of neighborhood penalty if the atom does not move in the SRep stage
		if (tuningWt(SRepAtomRestrictMovement) == 1) {
			sRepResults[0].value = 0;
			sRepResults[0].useIt = true;			
		}
		else {
			if (tuningWt(BpSRepAverageNeighborMatch) != 0.0) {
				PrimNeighborhoodDefn neiType = (enum PrimNeighborhoodDefn)
					(int) tuningWt(BpSRepNeighborhood);
				penalty = targetObject->dist2FromAveOfNeighbors(neiType,
					(int) tuningWt(BpFigureId), distType);	//MultiObject
				sRepResults[0].value = penalty;
				sRepResults[0].useIt = true;
				geomPenalty	+= tuningWt(sRepResults[0].tuningParm) * penalty;
			}
		}

		if (tuningWt(BpSRepImageMatch) != 0.0) {
			M3DQuadFigure* figure = dynamic_cast<M3DQuadFigure*>( targetObject->getFigurePtr( figureId ) );
			//interpolateSpoke(figure, 2, 4.5, 0);
			penalty = computeSpokeBinaryImageMatch(targetObject, figureId, atomId, spokeId, int(tuningWt(BpSRepSurfaceLevel))) ;
			sRepResults[1].value = penalty;
			sRepResults[1].useIt = true;
		}
		if (tuningWt(SRepRadiusMatch) != 0.0) {
			penalty = computeAtomRadiusMatch(targetObject, figureId, atomId);
			sRepResults[2].value = penalty;
			sRepResults[2].useIt = true;
			geomPenalty	+= tuningWt(sRepResults[2].tuningParm) * penalty;
		}
		if (tuningWt(BpSRepLandmarkMatch) != 0.0) {
			penalty = computeLandmarkPenalty(targetObject, figureId, landmarkObject, 
				dilationFactorInModelUnits, int(tuningWt(BpSRepSurfaceLevel)));
			sRepResults[3].value = penalty;
			sRepResults[3].useIt = true;
			geomPenalty	+= tuningWt(sRepResults[3].tuningParm) * penalty;
		}
		if (tuningWt(BpSRepModelMatch) != 0.0) {
			penalty = computeSpokeDistFromRef( targetObject, figureId, atomId, spokeId, distType ) ;			
			sRepResults[4].value = penalty;
			sRepResults[4].useIt = true;
			geomPenalty	+= tuningWt(sRepResults[4].tuningParm) * penalty;
		}
		if (tuningWt(BpSRepInterpenetrationMatch) != 0.0) {
			penalty = computeInterpenetrationPenalty(targetObject, figureId);
			sRepResults[5].value = penalty;
			sRepResults[5].useIt = true;
			geomPenalty	+= tuningWt(sRepResults[5].tuningParm) * penalty;
		}
		if (tuningWt(BpSRepRSradPenalty) != 0.0) {
			M3DFigure * figure = targetObject->getFigurePtr(figureId);
			if (tuningWt(BpSRepSumMatch) != 0.0)
				// Sum over all atoms.
				penalty = atomPredictor->getFigureRSradPenalty(
				figure, -1, PEL_BOTH, tuningWt(RSradPenaltyNorm), tuningWt(RSradPenaltyThreshold) );
			else
				penalty = atomPredictor->getFigureRSradPenalty(
				figure, atomId, PEL_BOTH, tuningWt(RSradPenaltyNorm), tuningWt(RSradPenaltyThreshold) );
			sRepResults[6].value = penalty;
			sRepResults[6].useIt = true;
			geomPenalty	+= tuningWt(sRepResults[6].tuningParm) * penalty;
		} 
		if (tuningWt(BpSRepDeltaNeighborMatch) != 0.0) {
			M3DQuadFigure * figure = dynamic_cast<M3DQuadFigure*>( targetObject->getFigurePtr(figureId));
			if (tuningWt(BpSRepSumMatch) != 0.0) {
				// Sum over all atoms
				for (int i = 0; i < figure->getPrimitiveCount(); i++)
					penalty = dynamic_cast<M3DAtomPredictorQuad*>(atomPredictor)->getAtomPenaltyWithoutOpt(figure, i);
				penalty /= figure->getPrimitiveCount();
			}
			else 
				penalty = dynamic_cast<M3DAtomPredictorQuad*>(atomPredictor)->getAtomPenaltyWithoutOpt(figure, atomId);
			sRepResults[7].value = penalty;
			sRepResults[7].useIt = true;
			geomPenalty	+= tuningWt(sRepResults[7].tuningParm) * penalty;
		}
		if ( tuningWt(BpSRepCurvinessPenalty) != 0.0) {
			penalty	= targetObject->getFigurePtr(figureId)->curviness();
			sRepResults[8].value = penalty;
			sRepResults[8].useIt = true;
			geomPenalty	+= tuningWt(sRepResults[8].tuningParm) * penalty;
		}

		//Xiaojie
		if (tuningWt(BpSRepImageNormalMatch) != 0.0) {
            penalty = computeSpokeBinaryImageNormalMatch(targetObject, figureId, atomId, spokeId, int(tuningWt(BpSRepSurfaceLevel)));
			sRepResults[9].value = penalty;
			sRepResults[9].useIt = true;
		}
		if (tuningWt(BpSRepImagePlaneOrienAtEnd) != 0.0) {
			penalty = computeSpokeBinaryImagePlaneOrienMatchAtEnd(targetObject, figureId, spokeId, atomId);
			sRepResults[10].value = penalty;
			sRepResults[10].useIt = true;
		}
		if (tuningWt(BpSRepImageVertexMatch) != 0.0) {
			penalty = computeSpokeBinaryImageVertexMatch(targetObject, figureId, atomId, spokeId);
			sRepResults[11].value = penalty;
			sRepResults[11].useIt = true;
		}

		for (i = 0; i < MAX_NUM_SREP_MATCH_RESULTS; i++) {
			double wtMatch;

			if (sRepResults[i].useIt) {
				sRepResults[i].useIt = false;
				wtMatch = tuningWt(sRepResults[i].tuningParm)*sRepResults[i].value;
				switch(sRepResults[i].fnClass) {
					case Components::ROUGH:
						totalPenalty.rough		+= wtMatch;
						totalPenalty.roughSet	= true;
						break;
					case Components::SMOOTH:
						totalPenalty.smooth		+= wtMatch;
						totalPenalty.smoothSet	= true;
						break;
					default:
						assert(false);
				}
				if (globalVerbosity >= 1) {
					cout << tuningShortName(sRepResults[i].tuningParm) << ':'
						<< tuningWt(sRepResults[i].tuningParm) << '*'
						<< sRepResults[i].value << ' ';
				}
			}
		}

		

		if (globalVerbosity >= 1)
			cout << " = " << totalPenalty << endl;
#else
		cout << "Binary image match is not implemented\n";
		return Components(FAILED);
#endif
	}
#ifdef BINARY
	else {
		cout << "Error: computeSRepMatchAndPenalties failed\n";
		return Components(FAILED);
	}
#else	/* BINARY */
	else {	// Greyscale Pablo
		if (type != HISTO_IMAGE_MATCH && type != MULTIPATCH_MATCH) {
			// Mask-based image match
#ifdef NO_MATCH_PENALTY
			penalty = 0.0;
#else
			penalty = computeMainFigureMatch(targetObject, figureId);
#endif
			sRepResults[0].value = penalty;
			sRepResults[0].useIt = true;
		}
		else if (type == HISTO_IMAGE_MATCH) {
			// Histogram image match
			if (sRepResults[1].weight != 0.0) {
				if (figureId < 0 || targetImage == NULL || targetObject == NULL) {
					cout << "Warning: worst histogram match being returned\n";
					return Components(FAILED);
				}
				if (spe == NULL) {
					cout << "Error: no spe\n";
					return Components(FAILED);
				}
				if (! spe->update(targetImage, targetObject, figureId, -1))
				{
					//	cout << "Error: surface patch (histogram) image match failed\n";
					return Components(FAILED);
				}				
				penalty = spe->getMatchValue(figureId, -1);


				//printf("Histogram match value %.4lf\n", penalty);
				sRepResults[1].value = penalty;
				sRepResults[1].useIt = true;

			}
		}
		else if (type == MULTIPATCH_MATCH)
		{
			if (sRepResults[2].weight != 0.0) {
				if (figureId < 0 || targetImage == NULL || targetObject == NULL) {
					cout << "Warning: worst MutliPatch match being returned\n";
					return Components(FAILED);
				}
				if (mp == NULL) {
					cout << "Error: no MultiPatch mp \n";
					return Components(FAILED);
				}

				penalty = mp->getYourMatch(targetObject, targetImage, figureId, NULL);

				//printf("MutliPatch match value %.4lf\n", penalty);
				//I'm piggybacking on Eli's histo match, cause it looks like a hassle to put in my own.
				sRepResults[2].value = penalty;
				sRepResults[2].useIt = true;

			}
		}
		else
		{
			cout << "ERROR: not sure what match you want...";
			return Components(FAILED);
		}

		if (sRepResults[3].weight != 0.0) {
			// Mahalanobis penalty
			penalty = computeAtomMahamPenalties(x, figureId, atomId);
			sRepResults[3].value = penalty;
			sRepResults[3].useIt = true;
		}
		else {
			// Geometric penalty
			if (sRepResults[4].weight != 0.0) {
				penalty = computeAtomPenalty(*targetObject, figureId,
					atomId, sRepResults[4].weight);
				sRepResults[4].value = penalty;
				sRepResults[4].useIt = true;
			}
		}

		// Constraints penalty
		if (sRepResults[6].weight > 0.0) {
			penalty = computeConstraintsPenalty(*targetObject, figureId);
			sRepResults[6].value = penalty;
			sRepResults[6].useIt = true;
		}

		// Now collect the terms
		if (globalVerbosity >= 2)
			cout << "  ";
		geomPenalty = 0.0;
		for (i = 0; i < MAX_NUM_ATOM_MATCH_RESULTS; i++) {
			double wtMatch;

			if (sRepResults[i].useIt) {
#ifndef OPTIMIZATION_VISUALIZER
				sRepResults[i].useIt = false;
#endif
				wtMatch = sRepResults[i].value*sRepResults[i].weight;
				switch(sRepResults[i].fnClass) {
					case Components::ROUGH:
						totalPenalty.rough		+= wtMatch;
						totalPenalty.roughSet	= true;
						break;
					case Components::SMOOTH:
						totalPenalty.smooth		+= wtMatch;
						totalPenalty.smoothSet	= true;
						break;
					default:
						assert(false);
				}
				// The image match entries do not contribute to geomPenalty
				if (i > LAST_IMAGE_MATCH_ATOM_TERM)
					geomPenalty += wtMatch;

				if (globalVerbosity >= 1) {
					cout << sRepResults[i].name << ':' << sRepResults[i].weight
						<< '*' << sRepResults[i].value << ' ';
				}
			}
		}

		if (globalVerbosity >= 1)
			cout << "= " << totalPenalty << endl;
	}
#endif	/* BINARY */

	return totalPenalty;
}

double Match::computeAtomMahamPenalties(const Vector & x, int atomId, int figureId)
{
	double mahamPenaltyVal = 0.0;	
	int i;

	for (i = 0; i < x.size(); i++)
		mahamPenaltyVal += x(i) * x(i);	

#ifndef BINARY
#ifdef WEIGHTED
	// IF NORMALIZE TO ZERO MEAN UNIT STD/VAR 
	mahamPenaltyVal -= x.size() ;
	mahamPenaltyVal /= sqrt((double)(2*x.size()) );			
#endif 
#endif			
	if (atomId >= 0 )
		atomMahamPenalty[atomId] =  mahamPenaltyVal;

	return mahamPenaltyVal;
}


/**
* Spoke stage penalty computation
*/
Components Match::computeSpokeMatchAndPenaltiesAsComponents(M3DObject * targetObject, int figureId,
															int atomId, const Vector& x, double & geomPenalty)
{
#ifdef BINARY
	DistanceType distType;
#endif
	double penalty;
	Components totalPenalty;
	totalPenalty.smooth		= 0.0;
	totalPenalty.rough		= 0.0;
	totalPenalty.smoothSet	= false;
	totalPenalty.roughSet	= false;
	int nTerms, i;

	nTerms = 0;

	// cleaning up list of surface bpoints
	if (targetSurfacePoints.points) {
		delete [] targetSurfacePoints.points;
		targetSurfacePoints.points = 0;
		targetSurfacePoints.numPoints = 0;
	}

	if (type == BINARY_IMAGE_MATCH) {
#ifdef BINARY
		if (globalVerbosity >= 1)
			cout << "   ";

		distType = (enum DistanceType) (int) tuningWt(BpSpokeDistanceType);
		// Convert units for dilationFactor.
		// (Distance map needs to have uniform imageToModelScale)
		// FIXME: Check for this condition.
		double dilationFactorInModelUnits = 0;
		if(tuningWt(BpDilationFactor) > 0) {
			const Vector3D imageToModelScale(
				1.0/binaryDistanceMap->getModelToImageScale()[0],
				1.0/binaryDistanceMap->getModelToImageScale()[1],
				1.0/binaryDistanceMap->getModelToImageScale()[2] );
			dilationFactorInModelUnits = tuningWt(BpDilationFactor) * imageToModelScale.getX();
		}
		if (tuningWt(BpSpokeAverageNeighborMatch) != 0.0) {
			PrimNeighborhoodDefn neiType = (enum PrimNeighborhoodDefn)
				(int) tuningWt(BpSpokeNeighborhood);
			penalty = targetObject->dist2FromAveOfNeighbors(neiType,
				(int) tuningWt(BpFigureId), distType);	//MultiObject
			spokeResults[0].value = penalty;
			spokeResults[0].useIt = true;
		}
		if (tuningWt(BpSpokeImageMatch) != 0.0) {
			penalty = computeBinaryImageMatch(targetObject, figureId,
				int(tuningWt(BpDDMEnable)), tuningWt(BpDDMAngleTolerance),
				int(tuningWt(BpSpokeSurfaceLevel)));
			spokeResults[1].value = penalty;
			spokeResults[1].useIt = true;
		}
		if (tuningWt(SpokeMahalanobisMatch) != 0.0) {
			penalty = computeSpokeMahamPenalties(x, figureId, atomId);
			spokeResults[2].value = penalty;
			spokeResults[2].useIt = true;
		}
		if (tuningWt(BpSpokeLandmarkMatch) != 0.0) {
			penalty = computeLandmarkPenalty(targetObject, figureId, landmarkObject, 
				dilationFactorInModelUnits, int(tuningWt(BpSpokeSurfaceLevel)));
			spokeResults[3].value = penalty;
			spokeResults[3].useIt = true;
		}
		if (tuningWt(BpSpokeModelMatch) != 0.0) {
			penalty = targetObject->dist2FromObject(targetObject->loadedObject(),
				(int) tuningWt(BpFigureId), distType);
			spokeResults[4].value = penalty;
			spokeResults[4].useIt = true;
		}
		if (tuningWt(BpSpokeInterpenetrationMatch) != 0.0) {
			penalty = computeInterpenetrationPenalty(targetObject, figureId);
			spokeResults[5].value = penalty;
			spokeResults[5].useIt = true;
		}
		if (tuningWt(BpSpokeRSradPenalty) != 0.0) {
			M3DFigure * figure = targetObject->getFigurePtr(figureId);
			if (tuningWt(BpSpokeSumMatch) != 0.0)
				// Sum over all atoms.
				penalty = atomPredictor->getFigureRSradPenalty(
				figure, -1, PEL_BOTH, tuningWt(RSradPenaltyNorm), tuningWt(RSradPenaltyThreshold) );
			else
				penalty = atomPredictor->getFigureRSradPenalty(
				figure, atomId, PEL_BOTH, tuningWt(RSradPenaltyNorm), tuningWt(RSradPenaltyThreshold) );
			spokeResults[6].value = penalty;
			spokeResults[6].useIt = true;
		} 

		for (i = 0; i < MAX_NUM_SPOKE_MATCH_RESULTS; i++) {
			double wtMatch;

			if (spokeResults[i].useIt) {
				spokeResults[i].useIt = false;
				wtMatch = tuningWt(spokeResults[i].tuningParm)*spokeResults[i].value;
				switch(spokeResults[i].fnClass) {
					case Components::ROUGH:
						totalPenalty.rough		+= wtMatch;
						totalPenalty.roughSet	= true;
						break;
					case Components::SMOOTH:
						totalPenalty.smooth		+= wtMatch;
						totalPenalty.smoothSet	= true;
						break;
					default:
						assert(false);
				}
				if (globalVerbosity >= 1) {
					cout << tuningShortName(spokeResults[i].tuningParm) << ':'
						<< tuningWt(spokeResults[i].tuningParm) << '*'
						<< spokeResults[i].value << ' ';
				}
			}
		}

		if (globalVerbosity >= 1)
			cout << " = " << totalPenalty << endl;
#else
		cout << "Binary image match is not implemented\n";
		return Components(FAILED);
#endif
	}
#ifdef BINARY
	else {
		cout << "Error: computeSpokeMatchAndPenalties failed\n";
		return Components(FAILED);
	}
#else	/* BINARY */
	else {	// Greyscale Pablo
		if (type != HISTO_IMAGE_MATCH && type != MULTIPATCH_MATCH) {
			// Mask-based image match
#ifdef NO_MATCH_PENALTY
			penalty = 0.0;
#else
			penalty = computeMainFigureMatch(targetObject, figureId);
#endif
			spokeResults[0].value = penalty;
			spokeResults[0].useIt = true;
		}
		else if (type == HISTO_IMAGE_MATCH) {
			// Histogram image match
			if (spokeResults[1].weight != 0.0) {
				if (figureId < 0 || targetImage == NULL || targetObject == NULL) {
					cout << "Warning: worst histogram match being returned\n";
					return Components(FAILED);
				}
				if (spe == NULL) {
					cout << "Error: no spe\n";
					return Components(FAILED);
				}
				if (! spe->update(targetImage, targetObject, figureId, -1))
				{
					//	cout << "Error: surface patch (histogram) image match failed\n";
					return Components(FAILED);
				}
				penalty = spe->getMatchValue(figureId, -1);
				//printf("Histogram match value %.4lf\n", penalty);
				spokeResults[1].value = penalty;
				spokeResults[1].useIt = true;
			}
		}
		else if (type == MULTIPATCH_MATCH)
		{
			// MultiPatch image match.  Again, I ride Eli's histo match so as til someone helps me 
			// put in a new spokeResult.
			if (spokeResults[1].weight != 0.0) {
				if (figureId < 0 || targetImage == NULL || targetObject == NULL) {
					cout << "Warning: worst multipatch match being returned\n";
					return Components(FAILED);
				}
				if (mp == NULL) {
					cout << "Error: no MultiPatch mp\n";
					return Components(FAILED);
				}

				penalty = mp->getYourMatch(targetObject, targetImage, figureId, NULL);

				//printf("MultiPatch match value %.4lf\n", penalty);
				spokeResults[1].value = penalty;
				spokeResults[1].useIt = true;
			}
		}
		else
		{
			cout << "ERROR: Hmmm, no spoke match defined.\n";
		}

		if (spokeResults[2].weight != 0.0) {
			// Mahalanobis penalty
			penalty = computeSpokeMahamPenalties(x, figureId, atomId);
			spokeResults[2].value = penalty;
			spokeResults[2].useIt = true;
		}
		else {
			// Geometric penalty
			if (spokeResults[3].weight != 0.0) {
				penalty = computeAtomPenalty(*targetObject, figureId,
					atomId, spokeResults[4].weight);
				spokeResults[3].value = penalty;
				spokeResults[3].useIt = true;
			}
		}

		// Constraints penalty
		if (spokeResults[5].weight > 0.0) {
			penalty = computeConstraintsPenalty(*targetObject, figureId);
			spokeResults[5].value = penalty;
			spokeResults[5].useIt = true;
		}

		// Now collect the terms
		if (globalVerbosity >= 2)
			cout << "  ";
		geomPenalty = 0.0;
		for (i = 0; i < MAX_NUM_SPOKE_MATCH_RESULTS; i++) {
			double wtMatch;

			if (spokeResults[i].useIt) {
#ifndef OPTIMIZATION_VISUALIZER
				spokeResults[i].useIt = false;
#endif
				wtMatch = spokeResults[i].value*spokeResults[i].weight;
				switch(spokeResults[i].fnClass) {
					case Components::ROUGH:
						totalPenalty.rough		+= wtMatch;
						totalPenalty.roughSet	= true;
						break;
					case Components::SMOOTH:
						totalPenalty.smooth		+= wtMatch;
						totalPenalty.smoothSet	= true;
						break;
					default:
						assert(false);
				}
				// The image match entries do not contribute to geomPenalty
				if (i > LAST_IMAGE_MATCH_SPOKE_TERM)
					geomPenalty += wtMatch;
				if (globalVerbosity >= 1) {
					cout << spokeResults[i].name << ':' << spokeResults[i].weight
						<< '*' << spokeResults[i].value << ' ';
				}
			}
		}

		//		if (globalVerbosity >= 2)
		cout << "= " << totalPenalty << endl;
	}
#endif	/* BINARY */

	return totalPenalty;
}

double Match::computeSpokeMahamPenalties(const Vector & x, int figureId, int atomId)
{
	double mahamPenaltyVal;
	double penaltyVal =0.0;
	M3DFigure * figure;
	int i, nAtoms;

	if (atomId < 0)
		return FAILED;

	figure = referenceObject->getFigurePtr(figureId);
	nAtoms = figure->getPrimitiveCount();

	mahamPenaltyVal = 0.0;
	for (i = 0; i < x.size(); i++)
		mahamPenaltyVal += x(i) * x(i);	

	return mahamPenaltyVal;
}



Components Match::computeFigureMatchAndPenaltiesAsComponents(M3DObject * targetObject, int figureId,
							     M3DObject * similarityObject, const std::vector<double>& x, double & geomPenalty )
{
#ifdef BINARY
	DistanceType distType;
#endif
	double penalty;
	Components totalPenalty;
	totalPenalty.smooth		= 0.0;
	totalPenalty.rough		= 0.0;
	totalPenalty.smoothSet	= false;
	totalPenalty.roughSet	= false;
	int nTerms, i;

	nTerms = 0;

	// cleaning up list of surface bpoints
	if (targetSurfacePoints.points) {
		delete [] targetSurfacePoints.points;
		targetSurfacePoints.points = 0;
		targetSurfacePoints.numPoints = 0;
	}

#ifdef OPTIMIZATION_VISUALIZER
	for (i = 0; i < MAX_NUM_FIG_MATCH_RESULTS; i++)
		figureResults[i].useIt = false;
#endif

	if (type == BINARY_IMAGE_MATCH) {
#ifdef BINARY
		if (globalVerbosity >= 1)
			cout << "  ";
		// Convert units for dilationFactor.
		// (Distance map needs to have uniform imageToModelScale)
		// FIXME: Check for this condition.
		double dilationFactorInModelUnits = 0;
		if(tuningWt(BpDilationFactor) > 0){
			const Vector3D imageToModelScale(
				1.0/binaryDistanceMap->getModelToImageScale()[0],
				1.0/binaryDistanceMap->getModelToImageScale()[1],
				1.0/binaryDistanceMap->getModelToImageScale()[2] );

			dilationFactorInModelUnits	= tuningWt(BpDilationFactor) * imageToModelScale.getX();
		}

		geomPenalty	= 0.0;
		if (tuningWt(BpFigureImageMatch) != 0.0) {
			penalty = computeBinaryImageMatch(targetObject, figureId,
				int(tuningWt(BpDDMEnable)), tuningWt(BpDDMAngleTolerance),
				int(tuningWt(BpFigureSurfaceLevel)));
			figureResults[0].value = penalty;
			figureResults[0].useIt = true;
		}

		if (tuningWt(BpFigureLandmarkMatch) != 0.0) {
			penalty = computeLandmarkPenalty(targetObject, figureId, landmarkObject,
				dilationFactorInModelUnits, int(tuningWt(BpFigureSurfaceLevel)));
			figureResults[1].value = penalty;
			figureResults[1].useIt = true;
			geomPenalty	+= tuningWt(figureResults[1].tuningParm) * penalty;
		}
		if (tuningWt(BpFigureModelMatch) != 0.0) {
			distType = (enum DistanceType) (int) tuningWt(BpFigureDistanceType);
			penalty = targetObject->dist2FromObject(similarityObject, 
				(int) tuningWt(BpFigureId), distType);
			figureResults[2].value = penalty;
			figureResults[2].useIt = true;
			geomPenalty	+= tuningWt(figureResults[2].tuningParm) * penalty;
		}
		if (tuningWt(BpFigureInterpenetrationMatch) != 0.0) {
			penalty = computeInterpenetrationPenalty(targetObject, figureId);
			figureResults[3].value = penalty;
			figureResults[3].useIt = true;
			geomPenalty	+= tuningWt(figureResults[3].tuningParm) * penalty;
		}
		if (tuningWt(FigureMahalanobisMatch) != 0.0) {
			penalty = 0.0;
			for (int i = 0; i < x.size(); i++)						
				penalty += x[i] * x[i];
			figureResults[4].value = penalty;
			figureResults[4].useIt = true;
			geomPenalty	+= tuningWt(figureResults[4].tuningParm) * penalty;
		}

		if (tuningWt(BpFigureRSradPenalty) != 0.0) {
			M3DFigure * figure = targetObject->getFigurePtr(figureId);
			penalty = atomPredictor->getFigureRSradPenalty(
				figure, -1, PEL_BOTH, tuningWt(RSradPenaltyNorm), tuningWt(RSradPenaltyThreshold) );
			figureResults[5].value = penalty;
			figureResults[5].useIt = true;
			geomPenalty	+= tuningWt(figureResults[5].tuningParm) * penalty;
		}
		if (tuningWt(BpFigureCurvinessPenalty) != 0.0) {
			penalty	= targetObject->getFigurePtr(figureId)->curviness();
			figureResults[6].value = penalty;
			figureResults[6].useIt = true;
			geomPenalty	+= tuningWt(figureResults[6].tuningParm) * penalty;
		}
		if (tuningWt(FigureContourMatch) != 0.0) {
			penalty = (contourMatch ? contourMatch->computeDistanceSqr(targetObject) : 0.0);
			figureResults[7].value = penalty;
			figureResults[7].useIt = true;
			geomPenalty	+= tuningWt(figureResults[7].tuningParm) * penalty;
		}

		//Xiaojie
		if (tuningWt(BpFigureImageNormalMatch) != 0.0) {
			penalty = computeFigureBinaryImageNormalMatch(targetObject, figureId);
			figureResults[8].value = penalty;
			figureResults[8].useIt = true;
		}
		if (tuningWt(BpFigureImagePlaneOrienAtEnd) != 0.0) {
			penalty = computeFigureBinaryImagePlaneOrienMatchAtEnd(targetObject, figureId);
			figureResults[9].value = penalty;
			figureResults[9].useIt = true;
		}
		if (tuningWt(BpFigureImageVertexMatch) != 0.0) {
			penalty = computeFigureBinaryImageVertexMatch(targetObject, figureId);
			figureResults[10].value = penalty;
			figureResults[10].useIt = true;
		}

		// Now collect the terms
		if (globalVerbosity >= 2)
			cout << "   ";
		for (i = 0; i < MAX_NUM_FIG_MATCH_RESULTS; i++) {
			double wtMatch;

			if (figureResults[i].useIt) {
#ifndef OPTIMIZATION_VISUALIZER
				figureResults[i].useIt = false;
#endif
				wtMatch = tuningWt(figureResults[i].tuningParm)*figureResults[i].value;
				switch(figureResults[i].fnClass) {
					case Components::ROUGH:
						totalPenalty.rough		+= wtMatch;
						totalPenalty.roughSet	= true;
						break;
					case Components::SMOOTH:
						totalPenalty.smooth		+= wtMatch;
						totalPenalty.smoothSet	= true;
						break;
					default:
						assert(false);
				}
				if (globalVerbosity >= 1) {
					cout << tuningShortName(figureResults[i].tuningParm) << ':'
						<< tuningWt(figureResults[i].tuningParm) << '*'
						<< figureResults[i].value << ' ';
				}
			}
		}
		if (globalVerbosity >= 1)
			cout << "= " << totalPenalty << endl;
#else
		cout << "Binary image match is not implemented\n";
		return Components(FAILED);
#endif
	}
#ifdef BINARY
	else {
		cout << "Error: computeFigureMatchAndPenalties failed\n";
		return Components(WORST_MATCH);
	}
#else	/* BINARY */
	else {	// Greyscale Pablo
		bool useInterpenetration;

		if (type != HISTO_IMAGE_MATCH && type != MULTIPATCH_MATCH) {
			// TOM: turn off match if constraint penalty is all the way up
			useInterpenetration = true;
			if (figureResults[0].weight == M3DMainFigureProblem::CONSTRAINT_PENALTY_SCALE_FACTOR)
				useInterpenetration = false;
			else {
				// Mask-based image match
#ifdef NO_MATCH_PENALTY
				penalty = 0.0;
#else
				penalty = computeMainFigureMatch(targetObject, figureId);
#endif
				figureResults[5].value = penalty;
				figureResults[5].useIt = true;
			}
		}
		else if (type == HISTO_IMAGE_MATCH) {
			// Histogram image match
			if (figureResults[4].weight != 0.0) {
				if (figureId < 0 || targetImage == NULL || targetObject == NULL) {
					cout << "Warning: worst histogram match being returned\n";
					return Components(FAILED);
				}
				if (spe == NULL) {
					cout << "Error: no spe\n";
					return Components(FAILED);
				}
				if (! spe->update(targetImage, targetObject, figureId, -1))
				{
					cout << "Error: surface patch (histogram) image match failed\n";
					return Components(FAILED);
				}
				penalty = spe->getMatchValue(figureId, -1);
				figureResults[4].value = penalty;
				figureResults[4].useIt = true;
			}
		}
		else if (type == MULTIPATCH_MATCH)
		{
			// Histogram image match
			if (figureResults[7].weight != 0.0) {
				if (figureId < 0 || targetImage == NULL || targetObject == NULL) {
					cout << "Warning: worst MultiPatch match being returned\n";
					return Components(FAILED);
				}
				if (mp == NULL) {
					cout << "Error: no MultiPatch mp\n";
					return Components(FAILED);
				}

				penalty = mp->getYourMatch(targetObject, targetImage, figureId, NULL);
				figureResults[7].value = penalty;
				figureResults[7].useIt = true;
			}
		}
		else
		{
			cout << "ERROR: Hmmm...., no figure match defined?\n";
		}

		//This is the DQF Match term.
		if (figureResults[6].weight != 0.0)
		{
			if (figureId < 0 || targetImage == NULL || targetObject == NULL) {
				cout << "Warning: worst dqf match being returned\n";
				return Components(FAILED);
			}
			if (dqf == NULL)
			{
				cout << "Error: no dqf\n";
				return Components(FAILED);
			}
			penalty = dqf->computeDQFMatch(targetImage, targetObject, figureId);
			figureResults[6].value = penalty;
			figureResults[6].useIt = true;
		}

		if (figureResults[2].weight != 0.0) {
			// Mahalanobis penalty
			penalty = 0.0;
			for (int i = 0; i < x.size(); i++)
				penalty += x[i] * x[i];

			if(spe && (spe->hmt == WEIGHTED_GAUSSIAN_HMT)) { // IF NORMALIZE TO ZERO MEAN UNIT STD/VAR 
				penalty -= x.size();
				penalty /= sqrt(2.0*x.size());
				//printf("Geom penalty: %lf, after weighting of %d\n", penalty, sqrt(2*(x.size()  - start)));
			}
			figureResults[2].value = penalty;
			figureResults[2].useIt = true;
		}
		else {
			// Geometric penalty
			if (figureResults[1].weight != 0.0) {
				penalty = computeMainFigurePenalty(targetObject, figureId);
				figureResults[1].value = penalty;
				figureResults[1].useIt = true;
			}
		}

		// Constraints penalty
		if (figureResults[0].weight > 0.0) {
			penalty = computeConstraintsPenalty(*targetObject, figureId,
				useInterpenetration);
			figureResults[0].value = penalty;
			figureResults[0].useIt = true;
		}

		// Contour penalty
		if (figureResults[3].weight >  0.0) {
			penalty = contourMatch->computeDistanceSqr(targetObject);
			// To make sense of what is happening, get the avg dist between the contour points
			// and boundary points of the current model in units of voxels.  computeDistanceSqr
			// returns the squared avg distance, so first get the total squared dist, take the
			// root, and then divide by number of points.
			//cout << "Avg distance to contour points in voxels: ";
			//cout << targetImage->getXDim()*sqrt(penalty*contourMatch->getNumPoints())/contourMatch->getNumPoints() << endl;
			penalty = pow(targetImage->getXDim()*sqrt(penalty*contourMatch->getNumPoints())
				/ contourMatch->getNumPoints(),4);

			figureResults[3].value = penalty;
			figureResults[3].useIt = true;
		}

		// Now collect the terms
		if (globalVerbosity >= 2)
			cout << "  ";
		geomPenalty = 0.0;
		for (i = 0; i < MAX_NUM_FIG_MATCH_RESULTS; i++) {
			double wtMatch;

			if (figureResults[i].useIt) {
#ifndef OPTIMIZATION_VISUALIZER
				figureResults[i].useIt = false;
#endif
				wtMatch = figureResults[i].value*figureResults[i].weight;
				switch(figureResults[i].fnClass) {
					case Components::ROUGH:
						totalPenalty.rough += wtMatch;
						totalPenalty.roughSet = true;
						break;
					case Components::SMOOTH:
						totalPenalty.smooth += wtMatch;
						totalPenalty.smoothSet = true;
						break;
					default:
						assert(false);
				}
				// The image match entries do not contribute to geomPenalty
				if (i < FIRST_IMAGE_MATCH_FIGURE_TERM)
					geomPenalty += wtMatch;
				if (globalVerbosity >= 2) {
					cout << figureResults[i].name << ':' << figureResults[i].weight
						<< '*' << figureResults[i].value << ' ';
				}
			}
		}

		if (globalVerbosity >= 2)
			cout << "= " << totalPenalty << endl;
	}
#endif	/* BINARY */

	return totalPenalty;
}

#ifndef BINARY

double Match::computeMainFigureMatch(M3DObject * targetObject, int figureId)
{
	int i;
	MaskElement * maskElement;
	M3DFigure * figure;
	Vector3D coord;

	double val1, val2;
	double matchValue;
	double mu, totval2Sqr, rms;
	double maskVolume;

	Mask * matchMask;

	if (figureId < 0 || figureId >= masks.size())
		return FAILED;

	matchMask = masks[figureId];
	if (targetImage == NULL || matchMask == NULL)
		return FAILED;

	if (matchMask->getSize() == 0)
		return FAILED;

	maskVolume = matchMask->getTotalVolume();

	matchValue = 0.0;

	figure = targetObject->getFigurePtr(figureId);

	// If the target image is stacked, we want to set the mask here,
	// according to which image we want.  Since no other information
	// seems to be available here right now, we make the simplifying
	// assumptions that all objects are single figures, and thus the
	// the figure index (or id I guess) tells which image to use...

	if (targetImage->getIsImageStacked())
		targetImage->setStackedMask((GreyValue) matchMask->getStackedImageMask());

	Xferlist * xferList = convertM3DtoXfer(figure);

	pList[figureId]->UpdateSubdivPointCloud(surface_level_for_match, xferList);

	int j;
	int numPoints;
	Bpoint * surfacePointList2;
	pList[figureId]->subdivboundaryinfo(&numPoints, &surfacePointList2);

	delete [] (xferList->atomlist);
	delete xferList;

	//DEBUGGING, I'm going to output the current target image and weights.
	//	FILE *f1, *f2;

	//f1 = fopen("current_target.txt", "wb");
	//f2 = fopen("current_volumes.txt", "wb");

#ifdef SURFACE_VISUALIZER
	objectSurfaceVisualizer->initializeSurfaceMap(0, "image intensities",
		targetObject, figureId, surface_level_for_match, numPoints, pList, false);
	objectSurfaceVisualizer->initializeSurfaceMap(1, "objective function",
		targetObject, figureId, surface_level_for_match, numPoints, pList, true);
	objectSurfaceVisualizer->initializeSurfaceMap(2, "vertex volume",
		targetObject, figureId, surface_level_for_match, numPoints, pList, true);
#endif

	mu = 0.0;
	totval2Sqr = 0.0;
	for (i = 0; i < matchMask->getSize(); i++) {
		maskElement = matchMask->getElement(i);
		if (maskElement == NULL) {
			printf("NULL mask element\n");
			continue;
		}
		if (maskElement->inBlendRegion)
			continue;

		// Add the distance along the normal.
		// Note that the distance is r Normalized.
		j = i/(matchMask->getSamplesPerPoint());
		coord.setX(surfacePointList2[j].pnt[0] + surfacePointList2[j].norm[0] *
			maskElement->dist * surfacePointList2[j].rad);
		coord.setY(surfacePointList2[j].pnt[1] + surfacePointList2[j].norm[1] *
			maskElement->dist * surfacePointList2[j].rad);
		coord.setZ(surfacePointList2[j].pnt[2] + surfacePointList2[j].norm[2] *
			maskElement->dist * surfacePointList2[j].rad);

		targetImage->modelToImageCoordinates(coord);

		if (windowIntensities)
			val2 = targetImage->mapDisplayToRelative(
			targetImage->getWindowedInterpolatedVoxelValue(
			coord.getX(), coord.getY(), coord.getZ()));
		else
			val2 = targetImage->mapDisplayToRelative(
			targetImage->getInterpolatedVoxelValue(
			coord.getX(), coord.getY(), coord.getZ()));

		//fprintf(f1,"%lf ", val2);
		//fprintf(f2,"%e ", maskElement->volelem);

		mu += val2 * maskElement->volelem;
		totval2Sqr += val2 * val2 * maskElement->volelem;
	}

	//fclose(f1);
	//fclose(f2);

	mu /= maskVolume;
	rms = (totval2Sqr / maskVolume) - (mu * mu);
	if(rms <= R_SMALL_TOLERANCE)
		rms = 0.0;
	else
		rms = sqrt(rms);

	if(rms <= R_SMALL_TOLERANCE)
		return 0.0;

	for (i = 0; i < matchMask->getSize(); i++) {
		maskElement = matchMask->getElement(i);
		if (maskElement == NULL) {
			printf("NULL mask element\n");
			continue;
		}
		if (maskElement->inBlendRegion)
			continue;

		// Add the distance along the normal.
		// Note that the distance is r Normalized.
		j=i/(matchMask->getSamplesPerPoint());
		coord.setX(surfacePointList2[j].pnt[0] + surfacePointList2[j].norm[0] *
			maskElement->dist * surfacePointList2[j].rad);
		coord.setY(surfacePointList2[j].pnt[1] + surfacePointList2[j].norm[1] *
			maskElement->dist * surfacePointList2[j].rad);
		coord.setZ(surfacePointList2[j].pnt[2] + surfacePointList2[j].norm[2] *
			maskElement->dist * surfacePointList2[j].rad);

		targetImage->modelToImageCoordinates(coord);

		val1 = maskElement->tempVal;
		if (windowIntensities)
			val2 = targetImage->mapDisplayToRelative(
			targetImage->getWindowedInterpolatedVoxelValue(
			coord.getX(), coord.getY(), coord.getZ()));
		else
			val2 = targetImage->mapDisplayToRelative(
			targetImage->getInterpolatedVoxelValue(
			coord.getX(), coord.getY(), coord.getZ()));

#ifdef SURFACE_VISUALIZER
		if ((i + 6) % 11 == 0)	// True for the 6th (middle) element of each profile
		{
			double intensity = targetImage->mapDisplayToRelative(
				targetImage->getVoxelValue(
				coord.getX(), coord.getY(), coord.getZ()));
			objectSurfaceVisualizer->setSurfaceColor(figureId, j, intensity);

			objectSurfaceVisualizer->setSurfaceColor(figureId, j,
				val1 * ((val2 - mu) / rms) * maskElement->volelem, 1);

			objectSurfaceVisualizer->setSurfaceColor(figureId, j,
				maskElement->volelem, 2);
		}
#endif
		matchValue += val1 * ((val2 - mu) / rms) * maskElement->volelem;
	}
	matchValue /= maskVolume;

	//#ifdef DEBUG
	printf("Match Value mfig: %lf\n", matchValue);
	//#endif

	// This negation is so that the objective function's terms may be added in
	// computeFigureMatchAndPenalties() and computeAtomMatchAndPenalties().
	// The image match is maximized, while other terms added there are
	// minimized.
	return -matchValue;		// This returns values in [-1.0, 1.0]
}

double Match::computeSubfigureMatch(M3DObject * targetObject, int figureId)
{
	int testFigureId;
	int i;

	MaskElement * maskElement;
	M3DFigure * figure;
	Vector3D coord;
	double val1, val2;
	double totval2Sqr, rms, mu;
	double maskVolume;
	double matchValue; 

	Mask * matchMask;

	if (figureId < 0 || figureId >= masks.size())
		return FAILED;

	matchMask = masks[figureId];

	if (targetImage == NULL || matchMask == NULL)
		return FAILED;

	if (matchMask->getSize() == 0)
		return FAILED;

	matchValue = 0.0;

	maskVolume = matchMask->getTotalVolume();

	figure = targetObject->getFigurePtr(figureId);

	Xferlist * xferList = convertM3DtoXfer(figure);

	pList[figureId]->UpdateSubdivPointCloud(surface_level_for_match, xferList);

	int j;
	int numPoints;
	Bpoint * surfacePointList2;
	pList[figureId]->subdivboundaryinfo(&numPoints, &surfacePointList2);

	mu = 0.0;
	totval2Sqr = 0.0;
	for (i = 0; i < matchMask->getSize(); i++) {
		maskElement = matchMask->getElement(i);
		if (maskElement == NULL) {
			printf("NULL mask element\n");
			continue;
		}
		if (maskElement->inBlendRegion)
			continue;

		// Add the distance along the normal.
		// Note that the distance is r Normalized.
		j=i/(matchMask->getSamplesPerPoint());
		coord.setX(surfacePointList2[j].pnt[0] + surfacePointList2[j].norm[0] *
			maskElement->dist * surfacePointList2[j].rad);
		coord.setY(surfacePointList2[j].pnt[1] + surfacePointList2[j].norm[1] *
			maskElement->dist * surfacePointList2[j].rad);
		coord.setZ(surfacePointList2[j].pnt[2] + surfacePointList2[j].norm[2] *
			maskElement->dist * surfacePointList2[j].rad);

		targetImage->modelToImageCoordinates(coord);

		if (windowIntensities)
			val2 = targetImage->mapDisplayToRelative(
			targetImage->getWindowedInterpolatedVoxelValue(
			coord.getX(), coord.getY(), coord.getZ()));
		else
			val2 = targetImage->mapDisplayToRelative(
			targetImage->getInterpolatedVoxelValue(
			coord.getX(), coord.getY(), coord.getZ()));

		mu += val2 * maskElement->volelem;
		totval2Sqr += val2 * val2 * maskElement->volelem;
	}

	mu /= maskVolume;
	rms = (totval2Sqr / maskVolume) - (mu * mu);
	if(rms <= R_SMALL_TOLERANCE)
		rms = 0.0;
	else
		rms = sqrt(rms);

	if(rms <= R_SMALL_TOLERANCE)
		return 0.0;

	int count = 0;
	for (i = 0; i < matchMask->getSize(); i++) {
		maskElement = matchMask->getElement(i);
		if (maskElement == NULL) {
			printf("NULL mask element\n");
			continue;
		}
		if (maskElement->inBlendRegion)
			continue;

		testFigureId = maskElement->figureId;
		if (testFigureId != figureId)
			continue;

		// Add the distance along the normal.
		// Note that the distance is r Normalized.
		count++;
		j=i/(matchMask->getSamplesPerPoint());
		coord.setX(surfacePointList2[j].pnt[0] + surfacePointList2[j].norm[0] *
			maskElement->dist * surfacePointList2[j].rad);
		coord.setY(surfacePointList2[j].pnt[1] + surfacePointList2[j].norm[1] *
			maskElement->dist * surfacePointList2[j].rad);
		coord.setZ(surfacePointList2[j].pnt[2] + surfacePointList2[j].norm[2] *
			maskElement->dist * surfacePointList2[j].rad);

		targetImage->modelToImageCoordinates(coord);

		val1 = maskElement->tempVal;
		if (windowIntensities)
			val2 = targetImage->mapDisplayToRelative(
			targetImage->getWindowedInterpolatedVoxelValue(
			coord.getX(), coord.getY(), coord.getZ()));
		else
			val2 = targetImage->mapDisplayToRelative(
			targetImage->getInterpolatedVoxelValue(
			coord.getX(), coord.getY(), coord.getZ()));

		matchValue += val1 * ((val2 - mu) / rms) * maskElement->volelem;
	}
	matchValue /= maskVolume;

	delete [] (xferList->atomlist);
	delete xferList;

#ifdef DEBUG
	printf("Match Value sub: %f\n", matchValue);
#endif
	return -matchValue;
}

#endif	/* ! BINARY */


#ifdef BINARY

// Dibyendu
// function to compute distance between a spoke and corresponding reference spoke in the s-rep stage

double Match::computeSpokeDistFromRef(M3DObject *targetObject, int figureId, int atomId, int spokeId, DistanceType DT) {

	if( figureId < 0 || atomId < 0 || spokeId < 0 || spokeId > 2 ) {
		cout << "Wrong input parameters to Match::computeSpokeDistFromRef() !" << endl ;
		return WORST_MATCH;
	}	

    if(dynamic_cast<M3DQuadFigure*>(this->getReferenceObject()->getFigurePtr(figureId))){
        M3DPrimitive * refAtom = this->getReferenceObject()->getFigurePtr(figureId)->getPrimitivePtr(atomId) ;
        M3DPrimitive * thisAtom = targetObject->getFigurePtr(figureId)->getPrimitivePtr(atomId) ;

        M3DQuadPrimitive * refQuadAtom = dynamic_cast<M3DQuadPrimitive*>( refAtom ) ;
        M3DQuadPrimitive * thisQuadAtom = dynamic_cast<M3DQuadPrimitive*>( thisAtom ) ;

        Vector3D refU ;
        Vector3D thisU ;

        double refR, thisR, effectiveR ;
        double penalty = 0.0 ;

        if( this->refObjectAfterScale != NULL )
            effectiveR = this->refObjectAfterScale->getFigurePtr(figureId)->getPrimitivePtr(atomId)->getR() ;
        else
            effectiveR = refAtom->getR() ;

        switch( spokeId ) {

            case 0 :	refR = refQuadAtom->getR0() ;
                        thisR = thisQuadAtom->getR0() ;

                        refU = refQuadAtom->getU0() ;
                        thisU = thisQuadAtom->getU0() ;
                        break ;

            case 1 :	refR = refQuadAtom->getR1() ;
                        thisR = thisQuadAtom->getR1() ;

                        refU = refQuadAtom->getU1() ;
                        thisU = thisQuadAtom->getU1() ;
                        break ;

            case 2 :    if( thisQuadAtom->type() == M3D_END_PRIMITIVE ) {

                        refR = (dynamic_cast<M3DEndPrimitive*>(refQuadAtom))->getREnd() ;
                        thisR = (dynamic_cast<M3DEndPrimitive*>(thisQuadAtom))->getREnd() ;

                        refU = refQuadAtom->getUEnd() ;
                        thisU = thisQuadAtom->getUEnd() ;
                        }
                        else {
                            cout << "Wrong spokeId ! This is not an end atom." << endl ;
                            return WORST_MATCH ;
                        }
                        break ;

            default :	cout << "Wrong spokeId ! Spoke Id has to be between 0 and 2." << endl ;
                        return WORST_MATCH ;
        }

        double delR = effectiveR * (log(thisR) - log(refR)) ;

        Quat q = ShapeSpace::S2::rotationToOrigin( refU ) ;
        q.normalize();
        q.rotateVector( thisU );

        Vector2D logMapDelU = ShapeSpace::S2::Log( thisU ) * effectiveR ;

        double delU = logMapDelU.norm() ;

        penalty = delR * delR + delU * delU  ;

        return( penalty ) ;
    }else{


        M3DPrimitive * refAtom = this->getReferenceObject()->getFigurePtr(figureId)->getPrimitivePtr(atomId) ;
        M3DPrimitive * thisAtom = targetObject->getFigurePtr(figureId)->getPrimitivePtr(atomId) ;

        M3DTubePrimitive * reftubeAtom = dynamic_cast<M3DTubePrimitive*>( refAtom ) ;
        M3DTubePrimitive * thistubeAtom = dynamic_cast<M3DTubePrimitive*>( thisAtom ) ;

        Vector3D refU ;
        Vector3D thisU ;

        double refR, thisR, effectiveR ;
        double penalty = 0.0 ;

        if( this->refObjectAfterScale != NULL )
            effectiveR = this->refObjectAfterScale->getFigurePtr(figureId)->getPrimitivePtr(atomId)->getR() ;
        else
            effectiveR = refAtom->getR() ;

        refR = reftubeAtom->getRN(spokeId) ;
        thisR = thistubeAtom->getRN(spokeId) ;

        refU = reftubeAtom->getNormalizedYN(spokeId) ;
        thisU = thistubeAtom->getNormalizedYN(spokeId) ;

        double delR = effectiveR * (log(thisR) - log(refR)) ;

        Quat q = ShapeSpace::S2::rotationToOrigin( refU ) ;
        q.normalize();
        q.rotateVector( thisU );

        Vector2D logMapDelU = ShapeSpace::S2::Log( thisU ) * effectiveR ;

        double delU = logMapDelU.norm() ;

        penalty = delR * delR + delU * delU  ;

        return( penalty ) ;
    }



}


// FEATURES:
// - use the distance map, the "search along normal" distance, or a combination of both. The search
//   is about 4x slower than the distance map, per surface evaluation. Of course the distance map
//   takes some time to compute, initially.
// - save distance vectors: a visualization of the distance map vectors and the "search along normal"
//   methods which can show, for each surface point, the direction and distance and other parameters
//   of the last method used. The data is saved in this routine, then displayed in the model display
//   routines. The data can be saved to a file for latter display. Rename the file to disable it's
//   display.

double Match::computeBinaryImageMatch(M3DObject * targetObject, int figureId,
									  int enableDDM, double angleTolerance, int subdivLevel)
{
	int i;
	M3DFigure * figure;

	if (figureId < 0 || targetImage == NULL || binaryDistanceMap == NULL)
		return WORST_MATCH;

	figure = targetObject->getFigurePtr(figureId);

    /*if (!targetSurfacePoints.points) {

		Xferlist * xferList = convertM3DtoXfer(figure);

		pList[figureId]->UpdateSubdivPointCloud(subdivLevel, xferList);
		pList[figureId]->copySubdivBoundaryInfo(&targetSurfacePoints.numPoints, &targetSurfacePoints.points);

		delete [] (xferList->atomlist);
		delete xferList;

    }*/

	std::vector <M3DSpoke*> relevantSpokes = getRelevantSpokes(targetObject, subdivLevel, figureId);
	int numSpokes = relevantSpokes.size();

	Bpoint * newPoints = new Bpoint[numSpokes];

	for (i = 0; i < numSpokes; i++)
	{
		newPoints[i].setvals(relevantSpokes[i]->getB().getX(), relevantSpokes[i]->getB().getY(), 
			relevantSpokes[i]->getB().getZ(), relevantSpokes[i]->getU().getX(), 
			relevantSpokes[i]->getU().getY(), relevantSpokes[i]->getU().getZ(), 
			relevantSpokes[i]->getR(), i, i, 0);
	}

	SurfacePointList targetSurfacePoints2;
	targetSurfacePoints2.points = newPoints;
	targetSurfacePoints2.numPoints = numSpokes;

	double match	= 0.0;
	double avg_dist	= 0.0;

    bool isOutside = false;

	double max_dist	= 0.0, min_dist	= +1e20, dist;
	// LATER: change this from voxels to cm

	const bool save_dvec = (tuningWt(SaveDistanceVectors) != 0);
	// Tells if "distance vectors" info should be saved

	static int infoMsg = 0;	// print info message(s) once only, then set this flag
	static int thinInfoMsg = 0;	// print info message(s) once only, then set this flag

    const double dilationFactor = tuningWt(BpDilationFactor);
	if (!thinInfoMsg && globalVerbosity > -1) {
		if( dilationFactor != 0.0 )
			cout << "INFO: Special handling for thin objects enabled."
			<< "Dilation Factor = " << dilationFactor << " voxels"
			<< endl;
		thinInfoMsg = 1;
	}

	const double endWeight	= tuningWt(WeightMatchAtEnds);

	double numPoints	= 0;

    if (enableDDM == 1) {

		// Sum distance map values (only)
		if (! infoMsg && globalVerbosity > -1) {
			cout << "Image-to-object distances enabled" << endl;
			infoMsg = 1;
		}

		for (i = 0; i < targetSurfacePoints2.numPoints; i++) {
			const double weight	= (pList[figureId]->thisMesh->isEndMedcoord(targetSurfacePoints2.points[i]))
				? endWeight : 1.0;
			numPoints	+= weight;
			dist = fabs(binaryDistanceMap->getDistance(targetSurfacePoints2.points[i].pnt)
				- dilationFactor );

			if (dist > max_dist)
				max_dist	= dist;
			if (dist < min_dist)
				min_dist = dist;
			match += weight * dist*dist;
			avg_dist += weight * dist;	// average unsigned distance
		}
	}
	else {
		// Sum either "search normal" distance or distance map; decision is per surface-point
		if (! infoMsg && globalVerbosity > -1) {
			if (enableDDM == 0)
				cout << "Object-to-image distances enabled"
				<< ", distance threshold = " << tuningWt(BpDDMSearchThreshold) << " voxels"
				<< endl;
			else
				cout << "Conditional object-to-image distances enabled"
				<< ", angle tolerance = " << angleTolerance << " degrees"
				<< ", distance threshold = " << tuningWt(BpDDMSearchThreshold) << " voxels"
				<< endl;
			infoMsg = 1;
		}
		int pointsHit = 0;	// Count # of surface points for either method

		DistanceVectorList::DistVectorPList & pointlist = distanceVectorList->plist;
		if (save_dvec)
            pointlist.resize(targetSurfacePoints2.numPoints);

		const Vector3D imageToModelScale(
			1.0/binaryDistanceMap->getModelToImageScale()[0],
			1.0/binaryDistanceMap->getModelToImageScale()[1],
			1.0/binaryDistanceMap->getModelToImageScale()[2] );

		// dp is delta voxel in unit co-ordinates.
		Vector3D dp = imageToModelScale;
		dp = dp / 2.;	// halve the step size as an initial guess
		double dpx = dp.getX();		// faster forms of dp
		double dpy = dp.getY();
		double dpz = dp.getZ();

		// The length of the projection of a vector at angleTolerance onto the unit vector
		const double toleranceLen = cos(angleTolerance*M_PI/180.0);

		for (i = 0; i < targetSurfacePoints2.numPoints; i++) {
			const double weight	= (pList[figureId]->thisMesh->isEndMedcoord(targetSurfacePoints2.points[i]))
				? endWeight : 1.0;
			numPoints	+= weight;
			const double x	= targetSurfacePoints2.points[i].pnt[0];
			const double y	= targetSurfacePoints2.points[i].pnt[1];
			const double z	= targetSurfacePoints2.points[i].pnt[2];
			Vector3D norm = targetSurfacePoints2.points[i].norm;
			norm.normalize();
			// Read distance map and adjust for dilation
            double dist = binaryDistanceMap->getDistance(x, y, z)
                - dilationFactor;	// start search here.

			Vector3D grad;

			// Save distance vector info to display later
			DistanceVector & p = (save_dvec) ? pointlist[i] : *((DistanceVector *) NULL);

			bool doNormalSearch = false;	// Should we search the normals?
			if (enableDDM == 0)
				doNormalSearch = true;
			else if (enableDDM == 2) {
				// The gradient is derived from the distance map, a piecewise-linear map which is tri-lerp'd by
				// getDistance() at positions slightly around (x,y,z). The second term (weighted by .0532) helps
				// decide cases in which getDistance(+dp) == getDistance(-dp).
				// FIXME: Check for bounds
				// (No need to adjust for dilation here, as we are only interested in differences)
                grad = 0.9679 * Vector3D(
					binaryDistanceMap->getDistance(x+dpx,y,z) - binaryDistanceMap->getDistance(x-dpx,y,z),
					binaryDistanceMap->getDistance(x,y+dpy,z) - binaryDistanceMap->getDistance(x,y-dpy,z),
					binaryDistanceMap->getDistance(x,y,z+dpz) - binaryDistanceMap->getDistance(x,y,z-dpz) )
					+ 0.0532 * Vector3D(
					binaryDistanceMap->getDistance(x+2*dpx,y,z) - binaryDistanceMap->getDistance(x-2*dpx,y,z),
					binaryDistanceMap->getDistance(x,y+2*dpy,z) - binaryDistanceMap->getDistance(x,y-2*dpy,z),
                    binaryDistanceMap->getDistance(x,y,z+2*dpz) - binaryDistanceMap->getDistance(x,y,z-2*dpz) );
                grad.normalize();

				doNormalSearch = (grad*norm <= toleranceLen);
			}

			bool isOutside;	// Is surface point outside binary image?
			if (dist > 0.0) {
				// If (x,y,z) is outside the binary image boundary, then invert the normal & the gradient
				norm = -norm;
				grad = -grad;
				isOutside = true;
			}
			else {
				dist = -dist;	// Make distance positive
				isOutside = false;
			}

			if (save_dvec) {
				p.n	= Vector3D(norm).vprod(imageToModelScale);
				p.p.set(x,y,z);
				p.grad	= Vector3D(grad).vprod(imageToModelScale * dist);
				p.dist	= dist;
				p.method= 1 + doNormalSearch;
			}

			if (doNormalSearch) {
				// Use "search normal" method for this surface point
				pointsHit++;

				const double ddm_dist = dist;		// ddm_dist should be positive

				Vector3D from = Vector3D(x,y,z);
				norm = norm.vprod(imageToModelScale); // a unit vector in image coords.
				// dist's direction is now *along grad direction*
				double step_size;
#define MAX_STEPS 10	// number of times to sample along the normal;

				Vector3D check	= from + dist*norm;

				// Start looking for the boundary from the known distance for 10 voxels
				int steps;
				bool done	= false;
				for (steps = 0; steps != MAX_STEPS; steps++ ) {
					// read distance map and adjust for dilation :
					const double intensity	= binaryDistanceMap->getDistance(check.getX(),check.getY(),check.getZ())
						- dilationFactor;

					// Ran outside the binary image -> stop early (pretend we ran out of steps) & use distance map
					if (intensity > 1000) {
						steps = MAX_STEPS;
						break;
					}
					// If we're "close enough" to the boundary, trust the dist map &
					// jump the rest of the way, then stop searching
					if (fabs(intensity) < tuningWt(BpDDMSearchThreshold)) {
						dist += (isOutside) ? intensity : -intensity;
						done	= true;
						break;
					}

					// If we passed the boundary, back up to the last place before the
					// boundary and repeat the search with a finer sampling
					if ((isOutside && intensity <= 0.0) || (! isOutside && intensity > 0.0)) {
						dist += (isOutside) ? intensity : -intensity;	// adjust for overstep
						break;
					}

					check = check + norm;
					dist = dist + 1.0;	// step size is 1 voxel
				}

				if (steps == MAX_STEPS)
					dist = ddm_dist;	// if we ran out of steps, just accept the dist map
				else if (dist <= 3.0 && ! done) {
					// earch small distances again, with a finer increment this time
					step_size	= dist/(MAX_STEPS-2);
					// (MAX_STEPS-2) has to be smaller than the number of iterations to
					// guarantee that we search past the initial stopping point
					check = from;
					dist = 0.0;
					for (steps = 0; steps != MAX_STEPS; steps++) {
						// read distance map and adjust for dilation :
						const double intensity	= binaryDistanceMap->getDistance(check.getX(),check.getY(),check.getZ())
							- dilationFactor;
						if(( isOutside && intensity <= 0.0)
							|| (! isOutside && intensity > 0.0)
							|| fabs(intensity) < tuningWt(BpDDMSearchThreshold)) {
								dist += (isOutside) ? intensity : -intensity;	// trust the dist map & jump the rest of the way
								break;
						}
						check = check + step_size*norm;
						dist = dist + step_size;
					}
				}

				// Now 'dist' is the furthest distance we searched, which might be short
				if (save_dvec) {
					p.n	= p.n*dist;
					p.dist	= dist;	// Overwrite ddm distance
				}
			}

			dist = fabs(dist);		// Might be negative if it was thresholded
			if (dist > max_dist)
				max_dist	= dist;
			if (dist < min_dist)
				min_dist = dist;
			match += weight * dist*dist;
			avg_dist += weight * dist;
		}
	}

	if( numPoints > 0 ) {
		match /= numPoints;
		avg_dist /= numPoints;
	}
	else {
		match = 0.0 ;
		avg_dist = 0.0 ;
	}

	MatchUtility::max_dist = max_dist;
	MatchUtility::min_dist = min_dist;
	MatchUtility::nPoints = numPoints;
	MatchUtility::avg_dist = avg_dist;
	MatchUtility::avg_dist_squared = match;

	return match;
}

// Dibyendu 
// this function computes the binary image match at the spoke ends and its interpolated neighbors only

double Match::computeSpokeBinaryImageMatch(M3DObject * targetObject, int figureId, int atomId, int spokeId, int subdivLevel)
{

	if( figureId < 0 || atomId < 0 || spokeId < 0 || spokeId > 2 || binaryDistanceMap == NULL ) {
		cout << "Wrong input parameters to Match::computeSpokeBinaryImageMatch() !" << endl ;
		return WORST_MATCH;
	}

	M3DFigure * figure = targetObject->getFigurePtr( figureId ) ;	
	M3DPrimitive * atom = figure->getPrimitivePtr( atomId ) ; 
	M3DQuadPrimitive* quadAtom = dynamic_cast<M3DQuadPrimitive*>( atom ) ;

	if( quadAtom == NULL ) {
        cout << "Match::computeSpokeBinaryImageMatch() is only for quad figures and quad primitives" << endl ;
        return WORST_MATCH ;
	}

	// dibyendu
	// There are several vector classes in use in Pablo. I shall use the std::vector wherever possible 
	// and use Vector3D to represent points

	std::vector <Vector3D> relevantSurfacePoints ;

	//vector<M3DSpoke*> relevantSpokes = getRelevantSpokes(targetObject, 1, figureId, 5, 0);

	//cout << relevantSpokes.size() << endl;
	
	/*for (int blah = 0; blah < relevantSpokes.size(); blah++)
	{
		cout << "{" << relevantSpokes[blah]->getB().getX() << "," << relevantSpokes[blah]->getB().getY() <<
			"," << relevantSpokes[blah]->getB().getZ() << "},";
	}

    cout << endl << endl;*/
	
	// Load the surface points - create a separate function when Jared's interpolation is in place. 
	// Till then, just take the tip of the spoke which is in question	

    Vector3D tempPoint = atom->getX() ;
	vector<M3DSpoke*> relevantSpokes;
	int numSpokes, i;

	switch( spokeId ) {

		case 0 :	//tempPoint = tempPoint + quadAtom->getY0() ;
					relevantSpokes = getRelevantSpokes(targetObject, subdivLevel, figureId, atomId, 0);
					numSpokes = relevantSpokes.size();

					//int i;
					for (i = 0; i < numSpokes; i++)
					{
						relevantSurfacePoints.push_back( relevantSpokes[i]->getB() );
					}


					break ;

		case 1 :	//tempPoint = tempPoint + quadAtom->getY1() ;
					relevantSpokes = getRelevantSpokes(targetObject, subdivLevel, figureId, atomId, 1);
					numSpokes = relevantSpokes.size();

					//int i;
					for (i = 0; i < numSpokes; i++)
					{
						relevantSurfacePoints.push_back( relevantSpokes[i]->getB() );
					}

					break ;

        case 2 :    if( atom->type() == M3D_END_PRIMITIVE )
                        //tempPoint = tempPoint + (dynamic_cast <M3DQuadEndPrimitive *> (quadAtom))->getYEnd() ;
                        //relevantSurfacePoints.push_back( tempPoint );

                    relevantSpokes = getRelevantSpokes(targetObject, subdivLevel, figureId, atomId, 2);
                    numSpokes = relevantSpokes.size();

                    for (i = 0; i < numSpokes; i++)
                    {
                        relevantSurfacePoints.push_back( relevantSpokes[i]->getB() );
                    }

					break ;

		default :	cout << "Wrong spokeId ! Spoke Id has to be between 0 and 2." << endl ;
					return WORST_MATCH ;					
	}

	//relevantSurfacePoints.push_back( tempPoint ) ;

	// Calculate the match by taking values from the distance map

	double match	= 0.0;
	double avg_dist	= 0.0;	

	double max_dist	= 0.0, min_dist	= +1e20, dist;
	

	static int infoMsg = 0;		// print info message(s) once only, then set this flag
	static int thinInfoMsg = 0;	// print info message(s) once only, then set this flag

	const double dilationFactor = tuningWt(BpDilationFactor);

	if (!thinInfoMsg && globalVerbosity > -1) {
		if( dilationFactor != 0.0 )
			cout << "INFO: Special handling for thin objects enabled."
			<< "Dilation Factor = " << dilationFactor << " voxels"
			<< endl;
		thinInfoMsg = 1;
	}

	double numPoints = relevantSurfacePoints.size() ;

	for( int i = 0; i < relevantSurfacePoints.size(); i++ ) {
		
		dist = fabs( binaryDistanceMap->getDistance(relevantSurfacePoints[i]) - dilationFactor );

		if( dist > max_dist )
			max_dist = dist ;
		if( dist < min_dist )
			min_dist = dist ;

		match += dist * dist ;

		avg_dist += dist ;	// average unsigned distance
	}	

	if( numPoints ) {
		match /= numPoints;
		avg_dist /= numPoints;
	}

	MatchUtility::max_dist = max_dist;
	MatchUtility::min_dist = min_dist;
	MatchUtility::nPoints = numPoints;
	MatchUtility::avg_dist = avg_dist;
	MatchUtility::avg_dist_squared = match;

	relevantSurfacePoints.clear() ;

	return match;
}
									  
// Xiaojie

double Match::computeFigureBinaryImageNormalMatch( M3DObject * targetObject, int figureId )
{

	M3DFigure * figure;

	if (figureId < 0 || targetImage == NULL || binaryDistanceMap == NULL)
		return WORST_MATCH;

	figure = targetObject->getFigurePtr(figureId);

	double match	= 0.0;

	double numPoints	= 0;

	//TOFIX
	//should take care link atom case
	for (int atomId = 0; atomId != figure->getPrimitiveCount(); atomId++){

		//TOCHECK
		//use the BpAtomImageVertexMatch indirectly
        match += computeAtomBinaryImageNormalMatch(targetObject, figureId,atomId, 0 );

	}

	match /= figure->getPrimitiveCount();

	return match;

}


// Xiaojie

double Match::computeAtomBinaryImageNormalMatch( M3DObject * targetObject, int figureId,int atomId, int subdivLevel  )
{


	if (figureId < 0 || targetImage == NULL || binaryDistanceMap == NULL)
		return WORST_MATCH;



    M3DFigure * figure = targetObject->getFigurePtr( figureId ) ;
    M3DQuadFigure* quadfig = dynamic_cast<M3DQuadFigure*>( figure );

    if(quadfig){

        //atom
        M3DPrimitive * atom;
        M3DPrimitive * atom_afterScale;

        atom = figure->getPrimitivePtr(atomId);

        // dibyendu - check for NULL condition for reference object after scale

        if( refObjectAfterScale != NULL )
            atom_afterScale = refObjectAfterScale->getPrimitivePtr(atomId);
        else
            atom_afterScale = atom ;


        double match	 = 0.0;

        double numPoints = 0;

        const double endMidWeight	= tuningWt(WeightNormalMatchAtEndMidSpoke);//ADD this weight

        //TOFIX:should take care of the hinge atom case
        //if not end atom
        //TOCHECK
        if(atom->type()==M3D_STANDARD_PRIMITIVE){
            //do not consider tube model now
            M3DQuadPrimitive* quadAtom = dynamic_cast<M3DQuadPrimitive*>(atom);
            M3DQuadPrimitive* quadAtom_afterScale = dynamic_cast<M3DQuadPrimitive*>(atom_afterScale);

            vector<M3DSpoke*> relevantSpokes = getRelevantSpokes(targetObject, subdivLevel, figureId, atomId);

            int i;

            vector <int>::size_type numSpokes = relevantSpokes.size();

            for (i = 0; i < numSpokes; i++)
            {

                //two spoke ends
                //Vector3D SpokeEnd0 = quadAtom->getX()+quadAtom->getY0();
                //Vector3D SpokeEnd1 = quadAtom->getX()+quadAtom->getY1();

                Vector3D SpokeEnd0 = relevantSpokes[i]->getB();
                i++;
                Vector3D SpokeEnd1 = relevantSpokes[i]->getB();

                //get_graddistance???sign???
                //gradient of distancemap at two spoke ends
                Vector3D grad0 = binaryDistanceMap->getGradDistance(SpokeEnd0);
                grad0.normalize();
                Vector3D grad1 = binaryDistanceMap->getGradDistance(SpokeEnd1);
                grad1.normalize();

                //TOCHECK

                double r_afterScale = quadAtom_afterScale->getR();

                double dot_product ;

                dot_product = min( 1.0, quadAtom->getU0()*grad0 ) ;
                dot_product = max(-1.0, dot_product ) ;

                double dist0 = acos( dot_product )* r_afterScale;

                dot_product = min( 1.0, quadAtom->getU1()*grad1 ) ;
                dot_product = max(-1.0, dot_product ) ;

                double dist1 = acos( dot_product )* r_afterScale;//  ;

                match += (dist0 + dist1)/2;
            }

            match /= (numSpokes/2);
        }
        //if end atom
        else
        {

            M3DQuadEndPrimitive* endAtom = dynamic_cast<M3DQuadEndPrimitive*>(atom);
            M3DQuadEndPrimitive* endAtom_afterScale = dynamic_cast<M3DQuadEndPrimitive*>(atom_afterScale);

            //middle spoke end position
            Vector3D midSpokeEnd = endAtom->getX()+endAtom->getYEnd();
            //other two spoke ends
            Vector3D SpokeEnd0 = endAtom->getX()+endAtom->getY0();
            Vector3D SpokeEnd1 = endAtom->getX()+endAtom->getY1();
            //get_graddistance???sign???
            //gradient of distancemap at three spoke ends
            Vector3D grad0 = binaryDistanceMap->getGradDistance(SpokeEnd0);grad0.normalize();
            Vector3D grad1 = binaryDistanceMap->getGradDistance(SpokeEnd1);grad1.normalize();
            Vector3D gradMid = binaryDistanceMap->getGradDistance(midSpokeEnd);gradMid.normalize();

            //TOCHECK
            //time r_init

            double r_afterScale = endAtom_afterScale->getR();
            double rEnd_afterScale = endAtom_afterScale->getREnd();

            double dot_product ;

            dot_product = min( 1.0, endAtom->getU0()*grad0 ) ;
            dot_product = max( -1.0, dot_product ) ;

            double dist0 = acos( dot_product ) * r_afterScale;//;

            dot_product = min( 1.0, endAtom->getU1()*grad1 ) ;
            dot_product = max( -1.0, dot_product ) ;

            double dist1 = acos( dot_product )* r_afterScale;// ;

            dot_product = min( 1.0, endAtom->getB()*gradMid ) ;
            dot_product = max( -1.0, dot_product ) ;

            double distMid = acos( dot_product )* rEnd_afterScale;// ;

            match = (dist0 + dist1 + distMid * endMidWeight)/( 2 + endMidWeight );
        }

        return match ;
    }else{

        M3DTubeFigure* tube = dynamic_cast<M3DTubeFigure*>( figure );


        //atom
        M3DPrimitive * atom;
        M3DPrimitive * atom_afterScale;

        atom = figure->getPrimitivePtr(atomId);

        // dibyendu - check for NULL condition for reference object after scale

        if( refObjectAfterScale != NULL )
            atom_afterScale = refObjectAfterScale->getPrimitivePtr(atomId);
        else
            atom_afterScale = atom ;


        double match	 = 0.0;

        double numPoints = 0;

        const double endMidWeight	= tuningWt(WeightNormalMatchAtEndMidSpoke);//ADD this weight


        M3DTubePrimitive* tubeAtom = dynamic_cast<M3DTubePrimitive*>(atom);
        M3DTubePrimitive* tubeAtom_afterScale = dynamic_cast<M3DTubePrimitive*>(atom_afterScale);

        vector<M3DSpoke*> relevantSpokes = getRelevantSpokes(targetObject, subdivLevel, figureId, atomId);

        int i;

        vector <int>::size_type numSpokes = relevantSpokes.size();

        for (i = 0; i < numSpokes; i++)
        {

            Vector3D SpokeEnd0 = relevantSpokes[i]->getB();
            Vector3D grad0 = binaryDistanceMap->getGradDistance(SpokeEnd0);
            grad0.normalize();

            double r_afterScale = tubeAtom_afterScale->getR();

            double dot_product ;

            dot_product = min( 1.0, relevantSpokes[i]->getU()*grad0 ) ;
            dot_product = max( -1.0, dot_product ) ;

            double dist0 = acos( dot_product ) * r_afterScale;

            match += dist0;
        }

        match /= numSpokes;

        return match ;

    }

}

// dibyendu -  the binaryImageNormalMatch penalty, only for the spoke and its interpolated neighbors 

double Match::computeSpokeBinaryImageNormalMatch( M3DObject * targetObject, int figureId, int atomId, int spokeId, int subdivLevel  )
{
	int i;

	if( figureId < 0 || atomId < 0 || spokeId < 0 || spokeId > 2 || binaryDistanceMap == NULL ) {
		cout << "Wrong input parameters to Match::computeSpokeBinaryImageNormalMatch() !" << endl ;
		return WORST_MATCH;
	}

	M3DFigure * figure = targetObject->getFigurePtr(figureId);

	//atom

	M3DPrimitive * atom = figure->getPrimitivePtr(atomId);

	M3DPrimitive * atom_afterScale;	

	// dibyendu - check for NULL condition for reference object after scale

	if( refObjectAfterScale != NULL )
		atom_afterScale = refObjectAfterScale->getPrimitivePtr(atomId);
	else
		atom_afterScale = atom ;

	M3DQuadPrimitive * quadAtom = dynamic_cast <M3DQuadPrimitive *> (atom) ;
	M3DQuadPrimitive * quadAtom_afterScale = dynamic_cast <M3DQuadPrimitive *> (atom_afterScale) ;

	if( quadAtom == NULL || quadAtom_afterScale == NULL ) {
		cout << "Match::computeSpokeBinaryImageNormalMatch is only for quad figure" << endl ;
		return WORST_MATCH ;
	}

	double match	 = 0.0 ;

	std::vector <Vector3D> spokeEnd ;		// end of the spoke
	std::vector <Vector3D> spokeDir ;		// direction of the spoke
	Vector3D gradDir ;		// direction of distance map gradient at the spoke end

	//std::vector <Vector3D> relevantSurfacePoints;
	std::vector <M3DSpoke*> relevantSpokes;
	
	double rMultiplier;
	int numSpokes;

	
	switch( spokeId ) {

		case 0 :	
			//spokeEnd = quadAtom->getX() + quadAtom->getY0() ;					
			//spokeDir = quadAtom->getU0() ;

            relevantSpokes = getRelevantSpokes( targetObject, subdivLevel, figureId, atomId, 0);
			numSpokes = relevantSpokes.size();

			for (i = 0; i < numSpokes; i++)
			{
				spokeEnd.push_back( relevantSpokes[i]->getB() );
				spokeDir.push_back( relevantSpokes[i]->getU() );
			}
			
			
			rMultiplier = quadAtom_afterScale->getR0() ;
			break ;

		case 1 :

			//spokeEnd = quadAtom->getX() + quadAtom->getY1() ;
			//spokeDir = quadAtom->getU1() ;

            relevantSpokes = getRelevantSpokes( targetObject, subdivLevel, figureId, atomId, 1);
			numSpokes = relevantSpokes.size();

			for (i = 0; i < numSpokes; i++)
			{
				spokeEnd.push_back( relevantSpokes[i]->getB() );
				spokeDir.push_back( relevantSpokes[i]->getU() );
			}


			rMultiplier = quadAtom_afterScale->getR1() ;
			break ;

            // This could be cleaned up now - Jared
		case 2 :

            if( ( atom->type() != M3D_END_PRIMITIVE ) || ( atom_afterScale->type() != M3D_END_PRIMITIVE ) ) {
				cout << "Wrong spokeId in computeSpokeBinaryImageNormalMatch. This is not an end atom" << endl ;
				return WORST_MATCH ;
			}

            //numSpokes = 1;

            relevantSpokes = getRelevantSpokes (targetObject, subdivLevel, figureId, atomId, 2);
            numSpokes = relevantSpokes.size();

            for (i = 0; i < numSpokes; i++)
            {
                spokeEnd.push_back( relevantSpokes[i]->getB() );
                spokeDir.push_back( relevantSpokes[i]->getU() );
            }

            //spokeEnd.push_back( quadAtom->getX() + (dynamic_cast <M3DQuadEndPrimitive *> (quadAtom))->getYEnd() );
            //spokeDir.push_back( (dynamic_cast <M3DQuadEndPrimitive *> (quadAtom))->getUEnd() );

			rMultiplier = (dynamic_cast <M3DQuadEndPrimitive *> (quadAtom_afterScale))->getREnd() ;
	}
	for (i = 0; i < numSpokes; i++)
	{
		gradDir = binaryDistanceMap->getGradDistance( spokeEnd[i] ) ;
		
		spokeDir[i].normalize() ;
		gradDir.normalize() ;

		double dotProduct = spokeDir[i] * gradDir ;

		dotProduct = (dotProduct < -1) ? -1 : dotProduct ;
		dotProduct = (dotProduct > 1 ) ? 1  : dotProduct ;

		match += acos( dotProduct ) * rMultiplier ;
	}

	match /= numSpokes;

	return match ;

}
//Xiaojie:
//in Figure Stage
//computeFigureBinaryImagePlaneOrienMatchAtEnd:
//To minimize the difference between the end atom spoke plane orientation and the plane orientations implied by the principal direction with larger principal 
//(using the the level surface of ddm the spoke end is on as the binary image) 

double Match::computeFigureBinaryImagePlaneOrienMatchAtEnd(M3DObject * targetObject, int figureId)
{
	M3DFigure * figure;

	if (figureId < 0 || targetImage == NULL || binaryDistanceMap == NULL)
		return WORST_MATCH;

	figure = targetObject->getFigurePtr(figureId);

	double numPoints = 0;

	double match	= 0.0;

	//TOFIX
	//should take care link atom case
	for (int atomId = 0; atomId != figure->getPrimitiveCount(); atomId++){

		M3DPrimitive * atom;
		atom = figure->getPrimitivePtr(atomId);
		if(atom->type() == M3D_END_PRIMITIVE)
			numPoints++;

		//TOCHECK
		//use the BpAtomImagePlaneOrienAtEnd indirectly
		match += computeAtomBinaryImagePlaneOrienMatchAtEnd(targetObject, figureId, atomId );

	}

	match /= numPoints;

	//cout<<"match"<<match<<endl;

	return match;

}
//Xiaojie:
//in Atom Stage
//computeAtomBinaryImagePlaneOrienMatchAtEnd:
//To minimize the difference between the end atom spoke plane orientation and the plane orientations implied by the principal direction with larger principal 
//(using the the level surface of ddm the spoke end is on as the binary image) 

double Match::computeAtomBinaryImagePlaneOrienMatchAtEnd(M3DObject * targetObject, int figureId, int atomId)
{

	M3DFigure * figure;

	if (figureId < 0 || targetImage == NULL || binaryDistanceMap == NULL)
		return WORST_MATCH;

	figure = targetObject->getFigurePtr(figureId);

	//atom
	M3DPrimitive * atom;
	M3DPrimitive * atom_afterScale;

	atom = figure->getPrimitivePtr(atomId);

	// dibyendu - ensure that refObjectAfterScale is not NULL

	if( refObjectAfterScale != NULL )
		atom_afterScale = refObjectAfterScale->getPrimitivePtr(atomId);
	else
		atom_afterScale = atom ;

	double match = 0;
	//TOFIX:should take care of the hinge atom case
	//if not end atom
	//TOCHECK
	if(atom->type()==M3D_STANDARD_PRIMITIVE)
		return 0;
	//if end atom
	else 
	{

        //cout<<atom->type()<<endl;


        if(dynamic_cast<M3DQuadEndPrimitive*>(atom)){
            M3DQuadEndPrimitive* endAtom = dynamic_cast<M3DQuadEndPrimitive*>(atom);
            M3DQuadEndPrimitive* endAtom_afterScale = dynamic_cast<M3DQuadEndPrimitive*>(atom_afterScale);

            //p1_spoke, p2_spoke
            Vector3D p1_spoke = endAtom->getN();
            Vector3D p2_spoke = endAtom->getBPerp();

            //TOFIX M3DEndPrimitive abstract??
            //middle spoke end position
            Vector3D midSpokeEnd = endAtom->getX()+endAtom->getYEnd();

            //end atom middle spoke length
            double rEnd_afterScale = endAtom_afterScale->getREnd();

            //p1_dist
            Vector3D p1_dist = binaryDistanceMap->getCurvature(midSpokeEnd).p1;
            Vector3D p2_dist = binaryDistanceMap->getCurvature(midSpokeEnd).p2;


            p1_spoke.normalize() ;
            p1_dist.normalize() ;
            p2_dist.normalize() ;

            //0417
            //double dot_product = fabs(p1_spoke * p1_dist) ;
            double dot_product = fabs(p2_spoke * p2_dist);

            //TOCHECK: umbilic
            if(dot_product == 0)
                dot_product = 1;

            if(dot_product > 1.0)
                dot_product = 1.0 ;
            if(dot_product < -1.0)
                dot_product = -1.0 ;

            match = acos(dot_product) * rEnd_afterScale;
        }else{

            M3DTubeEndPrimitive* endAtom = dynamic_cast<M3DTubeEndPrimitive*>(atom);
            M3DTubeEndPrimitive* endAtom_afterScale = dynamic_cast<M3DTubeEndPrimitive*>(atom_afterScale);

            //p1_spoke, p2_spoke
            Vector3D p1_spoke = endAtom->getN();
            Vector3D p2_spoke = endAtom->getBPerp();

            //TOFIX M3DEndPrimitive abstract??
            //middle spoke end position
            Vector3D midSpokeEnd = endAtom->getX()+endAtom->getU0();

            //end atom middle spoke length
            double rEnd_afterScale = endAtom_afterScale->getREnd();

            //p1_dist
            Vector3D p1_dist = binaryDistanceMap->getCurvature(midSpokeEnd).p1;
            Vector3D p2_dist = binaryDistanceMap->getCurvature(midSpokeEnd).p2;


            p1_spoke.normalize() ;
            p1_dist.normalize() ;
            p2_dist.normalize() ;

            //0417
            //double dot_product = fabs(p1_spoke * p1_dist) ;
            double dot_product = fabs(p2_spoke * p2_dist);

            //TOCHECK: umbilic
            if(dot_product == 0)
                dot_product = 1;

            if(dot_product > 1.0)
                dot_product = 1.0 ;
            if(dot_product < -1.0)
                dot_product = -1.0 ;

            match = acos(dot_product) * rEnd_afterScale;


        }

	}

	return match;

}


// dibyendu
// the binary image plane orientation penalty function only for the spoke and its interpolated neighbors

double Match::computeSpokeBinaryImagePlaneOrienMatchAtEnd(M3DObject *targetObject, int figureId, int atomId, int spokeId) {

	M3DFigure * figure;

	if (figureId < 0 || targetImage == NULL || binaryDistanceMap == NULL)
		return WORST_MATCH;

	figure = targetObject->getFigurePtr(figureId);

	//atom
	M3DPrimitive * atom;
	M3DPrimitive * atom_afterScale;

	atom = figure->getPrimitivePtr(atomId);

	// dibyendu - ensure that refObjectAfterScale is not NULL

	if( refObjectAfterScale != NULL )
		atom_afterScale = refObjectAfterScale->getPrimitivePtr(atomId);
	else
		atom_afterScale = atom ;

	double match = 0;
	//TOFIX:should take care of the hinge atom case
	//if not end atom
	//TOCHECK
	if(atom->type()==M3D_STANDARD_PRIMITIVE)
		return 0;
	//if end atom
	else 
	{

		M3DQuadEndPrimitive* endAtom = dynamic_cast<M3DQuadEndPrimitive*>(atom);
		M3DQuadEndPrimitive* endAtom_afterScale = dynamic_cast<M3DQuadEndPrimitive*>(atom_afterScale);

		//p1_spoke, p2_spoke
		Vector3D p1_spoke = endAtom->getN();
		Vector3D p2_spoke = endAtom->getBPerp();

		//TOFIX M3DEndPrimitive abstract??
		//middle spoke end position
		Vector3D midSpokeEnd = endAtom->getX()+endAtom->getYEnd();

		//end atom middle spoke length
		double rEnd_afterScale = endAtom_afterScale->getREnd();

		//p1_dist
		Vector3D p1_dist = binaryDistanceMap->getCurvature(midSpokeEnd).p1;
		Vector3D p2_dist = binaryDistanceMap->getCurvature(midSpokeEnd).p2;


		p1_spoke.normalize() ;
		p1_dist.normalize() ;		
		p2_dist.normalize() ;

		//0417
	    double dot_product = fabs(p1_spoke * p1_dist) ;	

		//TOCHECK: umbilic
		if(dot_product == 0)
			dot_product = 1;

		if(dot_product > 1.0)
			dot_product = 1.0 ;
		if(dot_product < -1.0)
			dot_product = -1.0 ;

		match = acos(dot_product) * rEnd_afterScale;

	}

	return match;

}

//Xiaojie:
//Figure stage
//computeBinaryVertexMatch:
//To minimize the difference between the middle spoke end of the end atom and crest point of the image 
//(using the the level surface of ddm the spoke end is on as the binary image) 
double Match::computeFigureBinaryImageVertexMatch(M3DObject * targetObject, int figureId)
{
	M3DFigure * figure;

	if (figureId < 0 || targetImage == NULL || binaryDistanceMap == NULL)
		return WORST_MATCH;

	figure = targetObject->getFigurePtr(figureId);

	double match	= 0.0;

	double numPoints = 0;

	//TOFIX
	//should take care link atom case
	for (int atomId = 0; atomId != figure->getPrimitiveCount(); atomId++){

		M3DPrimitive * atom;
		atom = figure->getPrimitivePtr(atomId);
		if(atom->type() == M3D_END_PRIMITIVE)
			numPoints++;

		//TOCHECK
		//use the BpAtomImageVertexMatch indirectly
		match += computeAtomBinaryImageVertexMatch(targetObject, figureId,atomId );

	}

	match /= numPoints;

	return match;

}
//Xiaojie:
//Atom stage
//computeBinaryVertexMatch:
//To minimize the difference between the middle spoke end of the end atom and crest point of the image 
//(using the the level surface of ddm the spoke end is on as the binary image) 

double Match::computeAtomBinaryImageVertexMatch(M3DObject * targetObject, int figureId, int atomId)
{	
	
	M3DFigure * figure;

	if (figureId < 0 || targetImage == NULL || binaryDistanceMap == NULL)
		return WORST_MATCH;

	figure = targetObject->getFigurePtr(figureId);

	//atom
	M3DPrimitive * atom;
	M3DPrimitive * atom_afterScale;

	atom = figure->getPrimitivePtr(atomId);

	// dibyendu - ensure that refObjectAfterScale is not NULL

	if( refObjectAfterScale != NULL )
		atom_afterScale = refObjectAfterScale->getPrimitivePtr(atomId);
	else
		atom_afterScale = atom ;

	double match = 0;

	//TOFIX:should take care of the hinge atom case
	//if not end atom
	//TOCHECK
	if(atom->type()==M3D_STANDARD_PRIMITIVE)
		return 0;
	//if end atom
	else 
	{

		M3DQuadEndPrimitive* endAtom = dynamic_cast<M3DQuadEndPrimitive*>(atom);
		M3DQuadEndPrimitive* endAtom_afterScale = dynamic_cast<M3DQuadEndPrimitive*>(atom_afterScale);

        if(endAtom && endAtom_afterScale){

            //p1_spoke, p2_spoke
            Vector3D p1_spoke = endAtom->getN();
            Vector3D p2_spoke = endAtom->getBPerp();

            //TOFIX M3DEndPrimitive abstract??
            //middle spoke end position
            Vector3D midSpokeEnd = endAtom->getX()+endAtom->getYEnd();

            //end atom middle spoke length
            double rEnd_afterScale = endAtom_afterScale->getREnd();

            //p1_dist
            Vector3D p1_dist = binaryDistanceMap->getCurvature(midSpokeEnd).p1;
            Vector3D p2_dist = binaryDistanceMap->getCurvature(midSpokeEnd).p2;

            //gradient of kappa1 at middle spoke end
            Vector3D grad_kappa1, grad_kappa2;
            binaryDistanceMap->getGradCurvature(midSpokeEnd, grad_kappa1, grad_kappa2);

            grad_kappa1.normalize() ;
            p1_spoke.normalize() ;

            match = fabs( p1_spoke * grad_kappa1) * rEnd_afterScale * rEnd_afterScale * rEnd_afterScale;
        }

	}

	return match;

}


// dibyendu
// the previous penalty function only for the spoke and its interpolated neighbors

double Match::computeSpokeBinaryImageVertexMatch(M3DObject * targetObject, int figureId, int atomId, int spokeId ) {
	
	M3DFigure * figure;

	if( figureId < 0 || atomId < 0 || spokeId < 0 || spokeId > 2 || binaryDistanceMap == NULL ) {
		cout << "Wrong input parameters to Match::computeSpokeBinaryImageVertexMatch() !" << endl ;
		return WORST_MATCH;
	}

	figure = targetObject->getFigurePtr(figureId);

	//atom
	M3DPrimitive * atom;
	M3DPrimitive * atom_afterScale;

	atom = figure->getPrimitivePtr(atomId);

	// dibyendu - ensure that refObjectAfterScale is not NULL

	if( refObjectAfterScale != NULL )
		atom_afterScale = refObjectAfterScale->getPrimitivePtr(atomId);
	else
		atom_afterScale = atom ;

	double match = 0;

	//TOFIX:should take care of the hinge atom case
	//if not end atom
	//TOCHECK
	if( ( atom->type() == M3D_STANDARD_PRIMITIVE ) || ( spokeId < 2 ) )
		return 0;
	//if end atom
	else 
	{

		M3DQuadEndPrimitive* endAtom = dynamic_cast<M3DQuadEndPrimitive*>(atom);
		M3DQuadEndPrimitive* endAtom_afterScale = dynamic_cast<M3DQuadEndPrimitive*>(atom_afterScale);
		//p1_spoke, p2_spoke
		Vector3D p1_spoke = endAtom->getN();
		Vector3D p2_spoke = endAtom->getBPerp();

		//TOFIX M3DEndPrimitive abstract??
		//middle spoke end position
		Vector3D midSpokeEnd = endAtom->getX()+endAtom->getYEnd();

		//end atom middle spoke length
		double rEnd_afterScale = endAtom_afterScale->getREnd();

		//p1_dist
		Vector3D p1_dist = binaryDistanceMap->getCurvature(midSpokeEnd).p1;
		Vector3D p2_dist = binaryDistanceMap->getCurvature(midSpokeEnd).p2;

		//gradient of kappa1 at middle spoke end
		Vector3D grad_kappa1, grad_kappa2;
		binaryDistanceMap->getGradCurvature(midSpokeEnd, grad_kappa1, grad_kappa2);

		grad_kappa1.normalize() ;
		p1_spoke.normalize() ;

		match = fabs( p1_spoke * grad_kappa1) * rEnd_afterScale * rEnd_afterScale * rEnd_afterScale;

	}

	return match;

}


// dibyendu
// penalize mismatch between the lengths of the two spokes (non-end) of an atom

double Match::computeFigureRadiusMatch( M3DObject* targetObject, int figureId ) 
{
	if (figureId < 0 || targetImage == NULL || binaryDistanceMap == NULL)
		return WORST_MATCH;

	double match = 0.0 ;

	M3DFigure * figure = targetObject->getFigurePtr(figureId);
	
	for (int atomId = 0; atomId != figure->getPrimitiveCount(); atomId++){

		match += computeAtomRadiusMatch(targetObject, figureId, atomId );

	}

	match /= figure->getPrimitiveCount();

	return match ;
}

// dibyendu
// penalize mismatch between the lengths of the two spokes (non-end) of an atom

double Match::computeAtomRadiusMatch( M3DObject* targetObject, int figureId, int atomId ) 
{
	if ( figureId < 0 || atomId <0 || targetImage == NULL || binaryDistanceMap == NULL )
		return WORST_MATCH;

	double match = 0.0 ;

	double r0 = 0.0 , r1 = 0.0 ;

	// get figure

	M3DFigure * figure = targetObject->getFigurePtr( figureId ) ;

	// get atom

	M3DPrimitive * atom = figure->getPrimitivePtr( atomId ) ;

	// pointers to store the actual derived class objects

	M3DQuadPrimitive * quadAtom = NULL ;

	M3DQuadEndPrimitive * quadEndAtom = NULL ;

	// dynamic cast depending on atom type

	switch( atom->type() ) {

		case M3D_STANDARD_PRIMITIVE:	

			quadAtom = dynamic_cast <M3DQuadPrimitive*>( atom );

            if(quadAtom){

                r0 = quadAtom->getR0() ;
                r1 = quadAtom->getR1() ;
            }else{

                M3DTubePrimitive* tubeAtom = dynamic_cast<M3DTubePrimitive*>(atom);

                if(tubeAtom){

                    cout<<"TODO double Match::computeAtomRadiusMatch( M3DObject* targetObject, int figureId, int atomId ) "<<endl;
                    return 0;

                }

            }

			break ;

		case M3D_END_PRIMITIVE:

			quadEndAtom = dynamic_cast <M3DQuadEndPrimitive*>( atom );

            if(quadEndAtom){

                r0 = quadEndAtom->getR0() ;
                r1 = quadEndAtom->getR1() ;


            }else{

                M3DTubeEndPrimitive* tubeAtom = dynamic_cast<M3DTubeEndPrimitive*>(atom);

                if(tubeAtom){

                    cout<<"TODO double Match::computeAtomRadiusMatch( M3DObject* targetObject, int figureId, int atomId ) "<<endl;
                    return 0;

                }

            }



			break ;
	}

	match = fabs( ( r0 - r1 ) / sqrt( r0 * r1 ) ) ;

	return match ;
}

#endif	/* BINARY */

#ifndef BINARY

// Penalty computation of similarity transformation of full model
double Match::computeModelPenalty(const SimilarityTransform3D & newTransformation)
{
	int figureCount;
	int numSurfacePoints;
	int i, j;

	double x, y, z;

	Bpoint * surfacePointPtr;
	double penalty;

	Vector3D coord,
		newCoord;

	SimilarityTransform3D transformation;


	transformation = initialTransform;
	transformation.multiplyBy(newTransformation);

	figureCount = surfacePoints.size();
	penalty = 0.0;
	numSurfacePoints = 0;
	for(i = 0; i < figureCount; i++)
	{
		surfacePointPtr = surfacePoints[i].points;
		for(j = 0; j < surfacePoints[i].numPoints; j++)
		{
			x = (surfacePointPtr->pnt[0]);
			y = (surfacePointPtr->pnt[1]);
			z = (surfacePointPtr->pnt[2]);
			surfacePointPtr++;

			coord.set((double) x, (double) y, (double) z);
			newCoord = coord;

			transformation.transformVector(newCoord);

			penalty += (newCoord - coord) * (newCoord - coord);
		}

		numSurfacePoints += surfacePoints[i].numPoints;
	}

	penalty = penalty / (double) numSurfacePoints;

	return penalty;
}
/*
// Penalty computation of similarity transform of object in multi-object case
double Match::computeObjectPenalty(const SimilarityTransform3D & newTransformation,
int treeIndex)
{
M3DFigureTreeNode * treeRoot;
int figureId;
int figureCount;
int numSurfacePoints;
int i;

double x, y, z;

Bpoint * surfacePointPtr;
double penalty;

Vector3D coord,
newCoord;

SimilarityTransform3D transformation;


if(referenceObject == NULL)
return FAILED;

treeRoot = referenceObject->getFigureTreeRoot(treeIndex);
if(treeRoot == NULL)
return FAILED;

figureId = treeRoot->getFigureId();
figureCount = surfacePoints.size();

if(figureId >= figureCount || figureId < 0)
return FAILED;

penalty = 0.0;
numSurfacePoints = surfacePoints[figureId].numPoints;

if(numSurfacePoints == 0)
return 0.0;

transformation = initialTransform;
transformation.multiplyBy(newTransformation);

surfacePointPtr = surfacePoints[figureId].points;
for(i = 0; i < numSurfacePoints; i++)
{
x = (surfacePointPtr->pnt[0]);
y = (surfacePointPtr->pnt[1]);
z = (surfacePointPtr->pnt[2]);
surfacePointPtr++;

coord.set((double) x, (double) y, (double) z);
newCoord = coord;

transformation.transformVector(newCoord);

penalty += (newCoord - coord) * (newCoord - coord);
}

penalty = penalty / (double) numSurfacePoints;

return penalty;
}
*/
// Penalty computation for main figure of object
double Match::computeMainFigurePenalty(M3DObject * targetObject, int figureId)
{
	M3DFigure * figure;
	int numSurfacePoints;
	int j;

	double x, y, z;

	Xferlist * xferList;
	Bpoint * surfacePointPtr;

	double penalty;

	Vector3D coord,
		newCoord;


	if(figureId >= surfacePoints.size())
		return 0.0;

	if(figureId < 0 || figureId >= surfacePoints.size())
	{
		printf("computePenalty - figure not saved in pointlists\n");
		return FAILED;
	}

	figure = dynamic_cast<M3DQuadFigure*>( targetObject->getFigurePtr(figureId));
	if(figure == NULL)
	{
		printf("computePenalty - figure == NULL.\n");
		return FAILED;
	}

	xferList = convertM3DtoXfer(figure);
	pList[figureId]->UpdateSubdivPointCloud(surface_level_for_match, xferList);

	int numPoints;
	Bpoint *surfacePointList2;
	pList[figureId]->subdivboundaryinfo(&numPoints, &surfacePointList2);

	penalty = 0.0;
	numSurfacePoints = 0;
	surfacePointPtr = surfacePoints[figureId].points;
	for(j = 0; j < surfacePoints[figureId].numPoints; j++)
	{
		x = surfacePointPtr->pnt[0];
		y = surfacePointPtr->pnt[1];
		z = surfacePointPtr->pnt[2];

		coord.set((double) x, (double) y, (double) z);

		newCoord.setX(surfacePointList2[j].pnt[0]);
		newCoord.setY(surfacePointList2[j].pnt[1]);
		newCoord.setZ(surfacePointList2[j].pnt[2]);

		penalty += (newCoord - coord) * (newCoord - coord);

		surfacePointPtr++;
	}

	penalty = penalty / (double) surfacePoints[figureId].numPoints;

	delete [] (xferList->atomlist);
	delete xferList;

	return penalty;
}

// Penalty computation for subfigure of object
double Match::computeSubfigurePenalty(M3DObject * targetObject, int figureId)
{
	M3DFigure * figure;
	Xferlist * xferList;    
	SurfacePointList surfaceList;
	int numPoints,
		i;

	Bpoint * savedPointsPtr;
	Bpoint * targetPointsPtr;

	double x1, y1, z1;
	double x2, y2, z2;
	double penalty;


	if(figureId < 0 || figureId >= surfacePoints.size())
	{
		printf("computePenalty - figure not saved in pointlists\n");
		return FAILED;
	}

	figure = targetObject->getFigurePtr(figureId);
	if(figure == NULL)
	{
		printf("computePenalty - figure == NULL.\n");
		return FAILED;
	}

	xferList = convertM3DtoXfer(figure);
	pList[figureId]->UpdateSubdivPointCloud(surface_level_for_match, xferList);
	pList[figureId]->subdivboundaryinfo(&numPoints, &targetPointsPtr);
	surfaceList = surfacePoints[figureId];
	if(surfaceList.numPoints != numPoints)
	{
		printf("computePenalty - number of surface points don't match.\n");
		return FAILED;
	}

	savedPointsPtr = surfaceList.points;
	penalty = 0.0;
	for(i = 0; i < numPoints; i++)
	{
		x1 = savedPointsPtr->pnt[0];
		y1 = savedPointsPtr->pnt[1];
		z1 = savedPointsPtr->pnt[2];
		savedPointsPtr++;

		x2 = (targetPointsPtr->pnt[0]);
		y2 = (targetPointsPtr->pnt[1]);
		z2 = (targetPointsPtr->pnt[2]);
		targetPointsPtr++;

		penalty += (x2 - x1) * (x2 - x1) +
			(y2 - y1) * (y2 - y1) +
			(z2 - z1) * (z2 - z1);
	}

	if(numPoints != 0)
		penalty /= (double) numPoints;

	delete [] (xferList->atomlist);
	delete xferList;

	return penalty;
}


// Classic atom penalty computation

// For an atom, the penalty is computed as the sum of the squares of the distances
// from the position predicted by each of the atom's neighbors to the proposed
// position of the atom.  Minimum penalty is attained if the atom moves to the
// average of the predicted positions of the atom's neighbors.
double Match::computeAtomPenalty(const M3DObject & targetObject,int figureId,
								 int atomId, double neighborPenaltyWeight)
{
	M3DQuadFigure * targetFigure;
	M3DQuadFigure * referenceFigure;
	M3DQuadPrimitive referencePredictor;
	M3DQuadPrimitive newPredictor;
	M3DQuadPrimitive newPrimitive;
	M3DQuadPrimitive referencePrimitive;
	M3DQuadPrimitive PredictedPrimitive;
	M3DQuadPrimitive averagePrimitive;
	int i,j;
	double penalty;
	double penalty1, penalty2;
	double numOfPred;
	Xferlist * xferList;
	Xferlist * xferList2;
	Xferlist * xferList3;
	int row, column;

	// Interpolation doesn't work below level 2
	int level = 3;

	targetFigure = dynamic_cast<M3DQuadFigure*>( targetObject.getFigurePtr(figureId));
	referenceFigure = dynamic_cast<M3DQuadFigure*>( referenceObject->getFigurePtr(figureId));
	if( typeid(*targetObject.getFigurePtr(figureId)) == typeid(M3DTubeFigure) ) {
		// FIXME:
		//	cout << "computeAtomPenalty() only works for quads.\n";
		//	assert(false);
		//sRad penalty
		M3DAtomPredictorTube atomPredictor;
		penalty = TUBE_GREY_ATOM_SRAD_PENALTY_WEIGHT * atomPredictor.getFigureRSradPenalty(
			targetObject.getFigurePtr(figureId), -1, PEL_BOTH, tuningWt(RSradPenaltyNorm), tuningWt(RSradPenaltyThreshold) );

		penalty1 = TUBE_GREY_ATOM_MODEL_PENALTY_WEIGHT * targetObject.getFigurePtr(figureId)->dist2FromFigure(
			referenceObject->getFigurePtr(figureId));

		penalty2 = TUBE_GREY_ATOM_CURVINESS_PENALTY_WEIGHT*  targetObject.getFigurePtr(figureId)->curviness();
#ifdef DEBUG_GREY_ATOM_TUBE
		cout <<"srad:"<<penalty<<"   model:"<<penalty1<<"   curviness:"<<penalty2<<"\n";
#endif

		return (penalty+ penalty1 +penalty2);
	}

	penalty = 0;
	numOfPred = 0;
	xferList = convertM3DtoXfer(referenceFigure);
	if (xferList == NULL) {
		cout << "Error: NULL xferList in computePenalty\n";
		return FAILED;
	}
	xferList2 = convertM3DtoXfer(targetFigure);
	if (xferList2 == NULL) {
		cout << "Error: NULL xferList(2) in computePenalty\n";
		return FAILED;
	}
	xferList3 = NULL;

	int figureCount = targetObject.getFigureCount();
	int numRows = targetFigure->getRowCount();
	int numColumns = targetFigure->getColumnCount();
	targetFigure->figuralCoordinates(atomId, row, column);

	pList[figureId]->UpdateSubdivBoundaryTiles(level, xferList);
	pList[figureCount + figureId]->UpdateSubdivBoundaryTiles(level,
		xferList2);

	ThallCode::Plist_subdivcomp Pcomp;
	Pcomp.init(pList[figureId]);

	// penalty1 is the d/r distance between target and reference objects
	// Stough: looking in comparemesh, it looks more like d^2/r.   Also, it seems to have 
	// virtually no effect on the geometric deviation, because the creasing penalty
	// dominates.  I'm changing the comparemesh to be (d/r)^2.
	penalty1 = Pcomp.comparemesh(pList[figureCount + figureId], column, row);

	// Now the neighbor penalty (penalty2)
	M3DQuadFigure neighborFigure = (*targetFigure);
	ThallCode::Plist_subdivcomp Pcomp2;
	M3DPrimitive * ap;

	penalty2 = 0.0;
	int numNeighbors = 0;

	referencePrimitive = *dynamic_cast<M3DQuadPrimitive*>(referenceFigure->getPrimitivePtr(atomId));
	newPrimitive = *dynamic_cast<M3DQuadPrimitive*>(targetFigure->getPrimitivePtr(atomId));
	for (i = -1; i <= 1; i = i + 2) {
		if (((row + i) < 0) || ((row + i) >= numRows))
			continue;
		referencePredictor = *dynamic_cast<M3DQuadPrimitive*>(referenceFigure->getPrimitivePtr(row + i, column));
		newPredictor = *dynamic_cast<M3DQuadPrimitive*>(targetFigure->getPrimitivePtr(row + i, column));
		M3DPrimitiveCorrector::predictPrimitive(PredictedPrimitive,referencePredictor,
			referencePrimitive, newPredictor);

		ap = new M3DQuadPrimitive(PredictedPrimitive);
		neighborFigure.setPrimitivePtr(row, column, ap);
		if (xferList3) {
			delete [] xferList3->atomlist;
			delete xferList3;
		}
		xferList3 = convertM3DtoXfer(&neighborFigure);
		if(xferList3 == NULL)
		{
			cout << "Error: xferList(3) in computePenalty\n";
			return FAILED;
		}

		pList[figureId]->UpdateSubdivBoundaryTiles(level, xferList3);
		Pcomp2.init(pList[figureId]);
		penalty2 += Pcomp2.comparemesh(pList[figureId], column, row);

		numNeighbors++;
	}

	for (j = -1; j <= 1; j = j + 2) {
		if (((column + j) < 0) || ((column + j) >= numColumns))
			continue;
		referencePredictor = *dynamic_cast<M3DQuadPrimitive*>(referenceFigure->getPrimitivePtr(row, column + j));
		newPredictor = *dynamic_cast<M3DQuadPrimitive*>(targetFigure->getPrimitivePtr(row, column + j));
		M3DPrimitiveCorrector::predictPrimitive(PredictedPrimitive,referencePredictor,
			referencePrimitive, newPredictor);

		ap = new M3DQuadPrimitive(PredictedPrimitive);
		neighborFigure.setPrimitivePtr(row, column, ap);
		if (xferList3) {
			delete [] xferList3->atomlist;
			delete xferList3;
		}
		xferList3 = convertM3DtoXfer(&neighborFigure);
		if(xferList3 == NULL)
		{
			cout << "Error: NULL xferList(3) in computePenalty\n";
			return FAILED;
		}

		pList[figureId]->UpdateSubdivBoundaryTiles(level, xferList3);
		Pcomp2.init(pList[figureCount + figureId]);
		penalty2 += Pcomp2.comparemesh(pList[figureId], column, row);

		numNeighbors++;
	}

	if(numNeighbors != 0)
		penalty2 /= numNeighbors;

	// Finally, penalty3 is the surface creasing (i.e., high curvature) penalty
	double penalty3 = 0.0;
	int count = 0;

	int numPoints;
	Bpoint * bPoints;
	ThallCode::Tileneighbors * neighbors;
	int index1, index2;
	Vector3D p1, p2, v1, v2, norm1, norm2;

	const static int MAX_NEIGHBOR_COUNT = 8;

	pList[figureId]->UpdateSubdivPointCloud(surface_level_for_match, xferList2);

	pList[figureId]->subdivboundaryinfo(&numPoints, &bPoints);
	pList[figureId]->subdivvertexneighbors(&numNeighbors, &neighbors);
	// qiong working
	for(i = 0; i < numPoints; i++)
	{
		p1.set(bPoints[i].pnt);

		index1 = neighbors[i][0];
		index2 = neighbors[i][1];
		if(index1 < 0 || index1>=numPoints || index2 < 0 || index2>=numPoints)//if(index1 < 0 || index2 < 0)
		{
			if(index1>=numPoints || index2>=numPoints)
				cout << "Invalid vertex neighbor index " << index1 << "," << index2 \
				<< " in Match::computeAtomPenalty(...)" << endl;
			continue;
		}   // 2003/07/07

		p2.set(bPoints[index1].pnt);
		v1 = p2 - p1;

		p2.set(bPoints[index2].pnt);
		v2 = p2 - p1;

		norm1 = v1.cross(v2);
		norm1.normalize();

		v1 = v2;

		for(j = 2; j < MAX_NEIGHBOR_COUNT; j++)
		{
			index2 = neighbors[i][j];

			if(index2 < 0 || index2>=numPoints)
			{
				if(index2>=numPoints)
					cout << "Invalid vertex neighbor index " << index2
					<< " in Match::computeAtomPenalty(...)" << endl;
				continue;
			}   // 2003/07/07

			p2.set(bPoints[index2].pnt);
			v2 = p2 - p1;

			norm2 = v1.cross(v2);
			norm2.normalize();
			penalty3 -= norm1 * norm2;


			norm1 = norm2;
			v1 = v2;

			count++;
		}
	}

	penalty = penalty1 * (1 - neighborPenaltyWeight) + 
		penalty2 * neighborPenaltyWeight;

	penalty += penalty3 * creasingPenaltyWeight;

#ifdef DEBUG
	cout << "Penalty1 = " << penalty1 << endl <<
		"Penalty2 = " << penalty2 << endl <<
		"Penalty3 = " << penalty3 << endl << 
		"TOTAL = " << penalty << endl;
#endif

	delete [] (xferList->atomlist);
	delete xferList;
	delete [] (xferList2->atomlist);
	delete xferList2;
	delete [] (xferList3->atomlist);
	delete xferList3; 

	return penalty;
} 

// Compute the geometric constraint penalty for figure FigureID in targetObject.
double Match::computeConstraintsPenalty(const M3DObject & targetObject,
										int figureId, bool useInterpenetration)
{
	double penalty;
	M3DFigure * targetFigure;
	M3DFigure * referenceFigure;
	bool validPoint;

	int numInvalidPoints = 0;

	referenceFigure = referenceObject->getFigurePtr(figureId);
	targetFigure = targetObject.getFigurePtr(figureId);

	using namespace ThallCode;

	Xferlist * xferList = convertM3DtoXfer(targetFigure);
	pList[figureId]->UpdateSubdivBoundaryTiles(CONSTRAINTS_RESOLUTION_BIAS +
		constraints_surface_level, xferList);

	delete [] (xferList->atomlist);
	delete xferList;

	penalty = 0.0;
	if (constrainedPointLists[figureId] == NULL)
		return penalty;		// This figure has no governors

	double normalPenalty = 0.0,
		penetrationPenalty = 0.0,
		slidingPenalty = 0.0;
	InterfiguralConstraints & governors = referenceFigure->inverseConstraints();
	for (int j = 0; j < governors.size(); j++) {

		int numConstrainedPoints;
		double refDist, trgDist;
		Vector3D targetPos, governorPos;
		referencePoint * point;
		double figuralPenalty;

		int govId = governors.figure(j);
		M3DFigure * govFigure = referenceObject->getFigurePtr(govId);

		// Compute boundary points for the governor figure
		xferList = convertM3DtoXfer(govFigure);

		pList[govId]->UpdateSubdivBoundaryTiles(CONSTRAINTS_RESOLUTION_BIAS +
			constraints_surface_level, xferList);

		delete [] (xferList->atomlist);
		delete xferList;

		figuralPenalty = 0.0;
		numConstrainedPoints = constrainedPointLists[figureId][govId].size();

		for(int n = 0; n < numConstrainedPoints; n++)
		{
			// Distance of a point on the reference figure from the governor
			point = &constrainedPointLists[figureId][govId][n];
			refDist = point->distance;

			Bpoint targetPoint, governorPoint;
			validPoint = pList[figureId]->subdivBposition(&targetPoint,
				point->u, point->v, point->t);
			if(!validPoint)
			{
				numInvalidPoints++;
				continue;
			}

			validPoint = pList[govId]->subdivBposition(&governorPoint,
				point->pu, point->pv, point->pt);
			if(!validPoint)
			{
				numInvalidPoints++;
				continue;
			}

			targetPos.set(targetPoint.pnt);
			Vector3D predPos;
			Vector3D actualPos;
			Vector3D normal;

			predPos.set(governorPoint.pnt);
			normal.set(governorPoint.norm);

			predPos += (normal * governorPoint.rad * point->distance);

			actualPos.set(targetPoint.pnt);

			// Normal component of difference
			Vector3D diff = predPos - actualPos;
			double normalComponent = diff * normal;
			double rSqr = governorPoint.rad * governorPoint.rad;
			normalPenalty += (normalComponent * normalComponent) / rSqr;

			// Tangential component of difference
			diff = (diff - normalComponent * normal);
			slidingPenalty += (diff * diff) / rSqr;

			//bool inBounds;
			trgDist = distanceMaps[govId]->getDistance(targetPos.getX(),
				targetPos.getY(), targetPos.getZ()/*, inBounds*/);
			double x = trgDist;
			// Penetration penalty is exponential for negative distance
			if(useInterpenetration && /*inBounds &&*/ x < 0.0)
				penetrationPenalty += exp(x * x);
		}
		if (numConstrainedPoints - numInvalidPoints > 0)
		{
			normalPenalty /= (numConstrainedPoints - numInvalidPoints);
			slidingPenalty /= (numConstrainedPoints - numInvalidPoints);
		}

		penalty += normalPenalty * normalPenaltyWeight +
			penetrationPenalty * penetrationPenaltyWeight +
			slidingPenalty * slidingPenaltyWeight;

		if (numInvalidPoints)
			cout << "Match::computeConstraintsPenalty() detected " << numInvalidPoints 
			<< " out of " << numConstrainedPoints << " invalid points" << endl;
	}

#ifdef DEBUG
	cout << "Geometric constraints penalty = " << penalty << endl;
	cout << "Normal penalty = " << normalPenalty << endl;
	cout << "Sliding penalty = " << slidingPenalty<< endl;
	cout << "Penetration penalty = " << penetrationPenalty << endl << endl;
#endif

	return penalty;
}

#endif	/* BINARY */

// This function computes a distance map for each governor used in computing
//   the geometric constraints penalty.
bool Match::initializeDistanceLists(M3DObject * object, int level)
{
	int figureId;
	//	bool firstTime;
	M3DFigure * figure;

#ifdef DEBUG
	cout << "Match::initializeDistanceLists() called" << endl;
#endif
	if (object == NULL)
		return false;

	if (level < 0) constraints_surface_level = 0;
	else constraints_surface_level = level;

	//	firstTime = true;
	for (int i = 0; i < distanceMaps.size(); i++) {
		if(distanceMaps[i] != NULL)
			delete distanceMaps[i];
		distanceMaps[i] = NULL;
	}
	distanceMaps.clear();

	int num_figures = object->getFigureCount();
	for (figureId = 0; figureId < num_figures; figureId++) {

		figure = object->getFigurePtr(figureId);

		InterfiguralConstraints & ifc = figure->constraints();
		if (ifc.size() == 0) {
			// Not a governor
			distanceMaps.push_back(NULL);
			continue;
		}

		// Create an object containing a copy of the governor
		M3DObject * tmp_object = new M3DObject;
		M3DFigure * tmp_figure = figure->assign();
		tmp_object->addFigure(tmp_figure);

		// Create a distance map from that object's boundary surface
		DistanceMap3D * map;
		BoundingBox bbox;
		int mapSize[3];

		// TOM: SET THE CUTOFF FOR THE DISTANCE MAP		AGG: Why set it to 2.0?
		map = new DistanceMap3D;
		map->defineBBox(tmp_object, bbox, mapSize,
			CONSTRAINTS_RESOLUTION_BIAS + constraints_surface_level, 2.0); // HARDCODED
		// AGG: if the cutoff could be given as a Euclidean distance, I could use
		// the max distance stored in ifc, instead of 2.0*r.
		map->initializeWithBBox(tmp_object, &bbox, mapSize,
			CONSTRAINTS_RESOLUTION_BIAS + constraints_surface_level);
		map->createMap(true);	// Produce an r-normalized map
		const int * map_size_ptr = map->getMapSize();
		if (map_size_ptr[0] <= 0) {
			delete tmp_object;
			cout << "map_size of 0!" << endl;
			return false;
		}
		if (map->status() == false) {
			delete tmp_object;
			cout << "No distance map!" << endl;
			return false;
		}
		else {
			/*
			if (firstTime) {
			constraintsMapSize[0] = map_size_ptr[0];
			constraintsMapSize[1] = map_size_ptr[1];
			constraintsMapSize[2] = map_size_ptr[2];
			firstTime = false;
			}
			else {
			if (constraintsMapSize[0] != map_size_ptr[0] || constraintsMapSize[1] != map_size_ptr[1]
			|| constraintsMapSize[2] != map_size_ptr[2])
			cout << "Warning: distance maps are not commensurate!" << endl;
			}
			*/
			distanceMaps.push_back(map);
		}
		delete tmp_object;
	}
	return true;
}

bool Match::updateDistanceList(M3DObject * object, int figureId) {
	M3DFigure * figure;

#ifdef DEBUG
	cout << "Match::updateDistanceList() called" << endl;
#endif
	figure = object->getFigurePtr(figureId);

	InterfiguralConstraints & ifc = figure->constraints();
	if (ifc.size() == 0) {
		// Not a governor
		distanceMaps.push_back(NULL);
		return true;
	}

	// Make sure we have enough distance maps
	if(distanceMaps.size() <= figureId)
		distanceMaps.resize(figureId + 1, NULL);

	// Discard the old map
	if(distanceMaps[figureId] != NULL)
	{
		delete distanceMaps[figureId];
		distanceMaps[figureId] = NULL;
	}

	// Create an object containing a copy of the governor
	M3DObject * tmp_object = new M3DObject;
	M3DFigure * tmp_figure = figure->assign();
	tmp_object->addFigure(tmp_figure);

	// Create a distance map from that object's boundary surface
	DistanceMap3D * map;
	BoundingBox bbox;
	int mapSize[3];

	// TOM: SET THE CUTOFF FOR THE DISTANCE MAP		AGG: Why set it to 2.0?
	map = new DistanceMap3D;
	map->defineBBox(tmp_object, bbox, mapSize,
		CONSTRAINTS_RESOLUTION_BIAS + constraints_surface_level, 2.0); // HARDCODED
	// AGG: if the cutoff could be given as a Euclidean distance, I could use
	// the max distance stored in ifc, instead of 2.0*r.
	map->initializeWithBBox(tmp_object, &bbox, mapSize,
		CONSTRAINTS_RESOLUTION_BIAS + constraints_surface_level);
	map->createMap(true);	// Produce an r-normalized map
	const int * map_size_ptr = map->getMapSize();
	if (map_size_ptr[0] <= 0) {
		delete tmp_object;
		cout << "map_size of 0!" << endl;
		return false;
	}
	if (map->status() == false) {
		delete tmp_object;
		cout << "No distance map!" << endl;
		return false;
	}

	// Store the new map
	distanceMaps[figureId] = map;

	delete tmp_object;
	return true;
}

// This function computes a list of points on the figural surface used in
//   computing the geometric constraints penalty.  Member function
//   initializeDistanceLists() must be called before using this function.
bool Match::initializeConstraintPointsLists(M3DObject * object, bool calcFiguralCoord)
{
	int figureId;
	M3DFigure * figure;
	int i, j;

#ifdef DEBUG
	cout << "Match::initializeConstraintPointsLists() called" << endl;
#endif
	if (object == NULL)
		return false;

	int num_figures = object->getFigureCount();

	//    if (constrainedPointLists == NULL) {
	//    if (num_figures != numPointLists) {
	if (constrainedPointLists != NULL) {
		for (i = 0; i < numPointLists; i++) {
			if (constrainedPointLists[i] != NULL) {
				for (int j = 0; j < numPointLists; j++)
					constrainedPointLists[i][j].clear();        // AGG: needed?
			}
			delete [] constrainedPointLists[i];
		}
		delete [] constrainedPointLists;
	}

	constrainedPointLists = new std::vector<referencePoint> *[num_figures];
	numPointLists = num_figures;
	for (figureId = 0; figureId < num_figures; figureId++)
		constrainedPointLists[figureId] =
		new std::vector<referencePoint> [num_figures];
	//	}
	//	else {
	//		for (figureId = 0; figureId < numPointLists; figureId++) {
	//			if (constrainedPointLists[figureId] != NULL) {
	//				for (i = 0; i < numPointLists; i++)
	//					constrainedPointLists[figureId][i].clear();
	//			}
	//		}
	//	}

	for (figureId = 0; figureId < num_figures; figureId++) {

		figure = object->getFigurePtr(figureId);

		InterfiguralConstraints & governors = figure->inverseConstraints();
#ifdef DEBUG
		cout << "Figure " << figureId << " has " << governors.size() << " governors\n";
#endif
		if (governors.size() == 0)
			continue;	// This figure is not governed

		using namespace ThallCode;

		int numBPoints, govId;
		Bpoint * bpoints;

		// Compute boundary points for the governed figure.
		// A dense point cloud is not used - just verticies of tiles.
		Xferlist * xferList = convertM3DtoXfer(figure);

		pList[figureId]->UpdateSubdivPointCloud(CONSTRAINTS_RESOLUTION_BIAS +
			constraints_surface_level, xferList);
		pList[figureId]->subdivboundaryinfo(&numBPoints, &bpoints);

#ifdef DEBUG
		cout << "Figure " << figureId << " has " << numBPoints
			<< " boundary points." << endl;
#endif

		// Create a list per governor of distances from boundary points on the
		// governed figure to the governor, that are inside the cut-off distance.
		// These distances are the reference distances in computeConstraintsPenalty().
		if (constrainedPointLists[figureId] == NULL)
			constrainedPointLists[figureId] =
			new std::vector<referencePoint>[numPointLists];
		for (j = 0; j < governors.size(); j++) {
			govId = governors.figure(j);

			M3DFigure * govFigure = object->getFigurePtr(govId);

			// Compute boundary points for the governor figure
			Xferlist * xferList2 = convertM3DtoXfer(govFigure);

			pList[govId]->UpdateSubdivBoundaryTiles(CONSTRAINTS_RESOLUTION_BIAS
				+ constraints_surface_level, xferList2);

			delete [] (xferList2->atomlist);
			delete xferList2;

			for (i = 0; i < numBPoints; i++) {
				double x = bpoints[i].pnt[0];
				double y = bpoints[i].pnt[1];
				double z = bpoints[i].pnt[2];

				Vector3D pos(x, y, z);

				// For each boundary point that is close enough, store its distance to
				// the governor and the u,v,t coordinates of the boundary point.
				referencePoint point;
				FiguralCoordinate figuralCoord;
				bool inBounds;
				double distance =
					distanceMaps[govId]->getDistance(x, y, z, inBounds);
				if (inBounds && distance < governors.distance(j)) {
					point.u = bpoints[i].u;
					point.v = bpoints[i].v;
					point.t = bpoints[i].t;
					if(calcFiguralCoord)
					{
						if (! distanceMaps[govId]->figuralCoordinates(figuralCoord, pos))
							cout << "Bad figural coordinate\n";
						point.pu = figuralCoord.u;
						point.pv = figuralCoord.v;
						point.pt = figuralCoord.t;

						Bpoint governorPoint;

						pList[govId]->subdivBposition(&governorPoint, point.pu,
							point.pv, point.pt);

						Vector3D govPos(governorPoint.pnt);

						point.distance = (pos - govPos).norm()/governorPoint.rad;
					}
					constrainedPointLists[figureId][govId].push_back(point);
				}
			}
		}
		delete [] (xferList->atomlist);
		delete xferList;
#ifdef DEBUG
		for (j = 0; j < governors.size(); j++) {
			govId = governors.figure(j);
			if (constrainedPointLists[figureId][govId].size() > 0)
				cout << "Governor " << govId << " affects " <<
				constrainedPointLists[figureId][govId].size()
				<< " points of figure " << figureId << endl;
		}
#endif
	}
#ifdef DEBUG
	for (i = 0; i < num_figures; i++) {
		if (constrainedPointLists[i] != NULL)
			for (int j = 0; j < num_figures; j++) {
				if (constrainedPointLists[i][j].size() > 0)
					cout << "constrainedPointLists(" << i << ", " << j << ") has " <<
					constrainedPointLists[i][j].size() << " bpoints" << endl;
			}
	}
#endif
	return true;
}

void Match::testSeurat(M3DFigure * figure, int surf_level) const {

	using namespace ThallCode;

	Pointlist_server2 pL;
	int numBPoints;
	Bpoint * bpoints;
	Bpoint bp;

	Xferlist * xferList = convertM3DtoXfer(figure);
	pL.init(xferList);
	pL.ComputeSubdivBoundaryTiles(surf_level);
	pL.subdivtileinfo(&numBPoints, &bpoints);

	int count = 0;
	cout << "numBPoints = " << numBPoints << endl;
	for (int i = 0; i < numBPoints; i++) {
		if (! pL.subdivBposition(&bp, bpoints[i].u, bpoints[i].v, bpoints[i].t)) {
			bpoints[i].print();
			cout << "^^^^^ Invalid ^^^^^" << endl;
			count++;
		}
	}

	Pointlist_server2 pList2;
	int numBPoints2;
	Bpoint * bpoints2;

	pList2.init(xferList);
	pList2.ComputeSubdivBoundaryTiles(surf_level);
	pList2.subdivtileinfo(&numBPoints2, &bpoints2);

	if (numBPoints != numBPoints2) {
		if (count)
			cout << "Pointlist contained " << count << " invalid results." << endl;
		cout << "Pointlist computation is indeterminate: bpoints received = "
			<< numBPoints << " and " << numBPoints2 << endl;
	}
	else {
		int match = 0;
		for (int i = 0; i < numBPoints; i++) {
			if (bpoints[i].u != bpoints2[i].u && bpoints[i].v != bpoints2[i].v
				&& bpoints[i].t != bpoints2[i].t && bpoints[i].pnt[0] != bpoints2[i].pnt[0]
			&& bpoints[i].pnt[1] != bpoints2[i].pnt[1] && bpoints[i].pnt[2] != bpoints2[i].pnt[2]
			&& bpoints[i].rad != bpoints2[i].rad && bpoints[i].norm[0] != bpoints2[i].norm[0]
			&& bpoints[i].norm[1] != bpoints2[i].norm[1] && bpoints[i].norm[2] != bpoints2[i].norm[2])
				match++;
		}
		if (match)
			cout << "Identically computed Pointlists yielded " << match
			<< " different Bpoints." << endl;
		else {
			if (count)
				cout << "Pointlist contained " << count << " invalid results." << endl;
			else
				cout << "Passed." << endl;
		}
	}

	delete [] (xferList->atomlist);
	delete xferList;
}

#ifdef BINARY

// Return sum(squared thickness) that penetrates inside image
double Match::computeInterpenetrationPenalty(M3DObject * targetObject, int figureId) {
	GreyValue mask = (GreyValue) tuningWt(RepulsionMask);
	if (mask == 0 || targetImage == NULL)
		return 0.0; 
	targetImage->pushImageIsStacked(true, mask);

	double penalty = 0;
	int numPoints;
	Bpoint * surfacePointList;
	pList[figureId]->subdivboundaryinfo(&numPoints, &surfacePointList);

	for (int i = 0; i < numPoints; i++) {
		Vector3D v(surfacePointList[i].norm);
		v.normalize();

		// Find thickness of object inside model: walk back
		// towards medial sheet integrating image intensities
		// until outside of object is found
		Vector3D boundary(surfacePointList[i].pnt);
		targetImage->modelToImageCoordinates(boundary);
		double x = boundary.getX();
		double y = boundary.getY();
		double z = boundary.getZ();
		double localPenalty = 0;
		double interpolatedValue = 0;
		while ((interpolatedValue = targetImage->getInterpolatedVoxelValue(x, y, z)) != 0) {
			// TODO: design penalty function
			localPenalty += interpolatedValue;

			// Probably need other bounds checks here
			// Are we out of the image?
			// Have we gotten back to the medial sheet of the object?
			x -= v.getX();
			y -= v.getY();
			z -= v.getZ();
		}
		penalty += (localPenalty * localPenalty);
	}
	//use an average value so that the same weights can be used
	//on different shapes and at different subdiv levels

	if( numPoints )
		penalty /= numPoints; 		

	targetImage->popImageIsStacked();
	return penalty;	
}


// Jared's Spoke Interpolation Code
M3DSpoke* Match::interpolateSpoke(M3DFigure *figure, double u, double v, int side)
{
	int ubase = (int)floor(u);
	int vbase = (int)floor(v);

	u = u - ubase;
	v = v - vbase;

	M3DQuadFigure *tempFigure = dynamic_cast<M3DQuadFigure*>(figure);

	// Get four corner atoms

	M3DPrimitive *atom11 = tempFigure->getPrimitivePtr(ubase,vbase);
	M3DQuadPrimitive* quadAtom11 = dynamic_cast<M3DQuadPrimitive*>( atom11 );

	M3DPrimitive *atom21 = tempFigure->getPrimitivePtr(ubase+1,vbase);
	M3DQuadPrimitive* quadAtom21 = dynamic_cast<M3DQuadPrimitive*>( atom21 );

	M3DPrimitive *atom12 = tempFigure->getPrimitivePtr(ubase,vbase+1);
	M3DQuadPrimitive* quadAtom12 = dynamic_cast<M3DQuadPrimitive*>( atom12 );

	M3DPrimitive *atom22 = tempFigure->getPrimitivePtr(ubase+1,vbase+1);
	M3DQuadPrimitive* quadAtom22 = dynamic_cast<M3DQuadPrimitive*>( atom22 );

//	cout << "Positions: " << endl;
//	cout << quadAtom11->getX().getX() << ", " << quadAtom11->getX().getY() << ", " << quadAtom11->getX().getZ() << endl;
//	cout << quadAtom21->getX().getX() << ", " << quadAtom21->getX().getY() << ", " << quadAtom21->getX().getZ() << endl;
//	cout << quadAtom12->getX().getX() << ", " << quadAtom12->getX().getY() << ", " << quadAtom12->getX().getZ() << endl;
//	cout << quadAtom22->getX().getX() << ", " << quadAtom22->getX().getY() << ", " << quadAtom22->getX().getZ() << endl;

//	cout << "Spokes: " << endl;
//	cout << quadAtom11->getU0().getX() << ", " << quadAtom11->getU0().getY() << ", " << quadAtom11->getU0().getZ() << endl;
//	cout << quadAtom21->getU0().getX() << ", " << quadAtom21->getU0().getY() << ", " << quadAtom21->getU0().getZ() << endl;
//	cout << quadAtom12->getU0().getX() << ", " << quadAtom12->getU0().getY() << ", " << quadAtom12->getU0().getZ() << endl;
//	cout << quadAtom22->getU0().getX() << ", " << quadAtom22->getU0().getY() << ", " << quadAtom22->getU0().getZ() << endl;

	// Get neighboring atoms for use in derivatives
	// If ubase/vbase is 0 or ubase+1/vbase+1 is max, cant use neighboring, use current instead
	// Currently using non-central differences for ease, need to fix later
	/*M3DPrimitive *atom10 = tempFigure->getPrimitivePtr(1,0);
	M3DQuadPrimitive* quadAtom10 = dynamic_cast<M3DQuadPrimitive*>( atom10 );

	M3DPrimitive *atom20 = tempFigure->getPrimitivePtr(2,0);
	M3DQuadPrimitive* quadAtom20 = dynamic_cast<M3DQuadPrimitive*>( atom20 );

	M3DPrimitive *atom31 = tempFigure->getPrimitivePtr(3,1);
	M3DQuadPrimitive* quadAtom31 = dynamic_cast<M3DQuadPrimitive*>( atom31 );

	M3DPrimitive *atom32 = tempFigure->getPrimitivePtr(3,2);
	M3DQuadPrimitive* quadAtom32 = dynamic_cast<M3DQuadPrimitive*>( atom32 );

	M3DPrimitive *atom23 = tempFigure->getPrimitivePtr(2,3);
	M3DQuadPrimitive* quadAtom23 = dynamic_cast<M3DQuadPrimitive*>( atom23 );

	M3DPrimitive *atom13 = tempFigure->getPrimitivePtr(1,3);
	M3DQuadPrimitive* quadAtom13 = dynamic_cast<M3DQuadPrimitive*>( atom13 );

	M3DPrimitive *atom02 = tempFigure->getPrimitivePtr(0,2);
	M3DQuadPrimitive* quadAtom02 = dynamic_cast<M3DQuadPrimitive*>( atom02 );

	M3DPrimitive *atom01 = tempFigure->getPrimitivePtr(0,1);
	M3DQuadPrimitive* quadAtom01 = dynamic_cast<M3DQuadPrimitive*>( atom01 );*/

	// Take the positions of the atoms and do simple finite difference derivative calculations
	/*Vector3D x10 = quadAtom10->getX();
	Vector3D x20 = quadAtom20->getX();
	Vector3D x31 = quadAtom31->getX();
	Vector3D x32 = quadAtom32->getX();
	Vector3D x23 = quadAtom23->getX();
	Vector3D x13 = quadAtom13->getX();
	Vector3D x02 = quadAtom02->getX();
	Vector3D x01 = quadAtom01->getX();*/

	Vector3D x11 = quadAtom11->getX();
	Vector3D x21 = quadAtom21->getX();
	Vector3D x12 = quadAtom12->getX();
	Vector3D x22 = quadAtom22->getX();

	Vector3D U0_11, U0_12, U0_21, U0_22;

	if (side == 0)
	{
		U0_11 = quadAtom11->getU0();
		U0_21 = quadAtom21->getU0();
		U0_12 = quadAtom12->getU0();
		U0_22 = quadAtom22->getU0();
	}
	else
	{
		U0_11 = quadAtom11->getU1();
		U0_21 = quadAtom21->getU1();
		U0_12 = quadAtom12->getU1();
		U0_22 = quadAtom22->getU1();
	}

	Vector3D du11, du21, du12, du22, dv11, dv21, dv12, dv22;

	Vector3D dU0du11, dU0dv11, dU0du21, dU0dv21, dU0du12, dU0dv12, dU0du22, dU0dv22;

	// Get derivatives for each corner, starting at 11, then 21, 12, and 22

	// Corner 11
	// du

	if (ubase == 0)
	{
		du11 = x21 - x11;
		dU0du11 = U0_21 - U0_11;
	}
	else
	{
		M3DQuadPrimitive* tempAtom = dynamic_cast<M3DQuadPrimitive*>(tempFigure->getPrimitivePtr(ubase-1, vbase));
		du11 = (x21 - tempAtom->getX()) / 2;

		if (side == 0)
			dU0du11 = (U0_21 - tempAtom->getU0()) / 2;
		else
			dU0du11 = (U0_21 - tempAtom->getU1()) / 2;
	}

	//dv

	if (vbase == 0)
	{
		dv11 = x12 - x11;
		dU0dv11 = U0_12 - U0_11;
	}
	else
	{
		M3DQuadPrimitive* tempAtom = dynamic_cast<M3DQuadPrimitive*>(tempFigure->getPrimitivePtr(ubase, vbase-1));
		dv11 = (x12 - tempAtom->getX()) / 2;

		if (side == 0)
			dU0dv11 = (U0_12 - tempAtom->getU0()) / 2;
		else
			dU0dv11 = (U0_12 - tempAtom->getU1()) / 2;
	}

	// Corner 21

	// du

	if (ubase + 1 == tempFigure->getRowCount() - 1)
	{
		du21 = x21 - x11;
		dU0du21 = U0_21 - U0_11;
	}
	else
	{
		M3DQuadPrimitive* tempAtom = dynamic_cast<M3DQuadPrimitive*>(tempFigure->getPrimitivePtr(ubase+2, vbase));
		du21 = (tempAtom->getX() - x11) / 2;

		if (side == 0)
			dU0du21 = (tempAtom->getU0() - U0_11) / 2;
		else
			dU0du21 = (tempAtom->getU1() - U0_11) / 2;
	}

	//dv

	if (vbase == 0)
	{
		dv21 = x22 - x21;
		dU0dv21 = U0_22 - U0_21;
	}
	else
	{
		M3DQuadPrimitive* tempAtom = dynamic_cast<M3DQuadPrimitive*>(tempFigure->getPrimitivePtr(ubase+1, vbase-1));
		dv21 = (x22 - tempAtom->getX()) / 2;

		if (side == 0)
			dU0dv21 = (U0_22 - tempAtom->getU0()) / 2;
		else
			dU0dv21 = (U0_22 - tempAtom->getU1()) / 2;
	}

	// Corner 12

	// du

	if (ubase == 0)
	{
		du12 = x22 - x12;
		dU0du12 = U0_22 - U0_12;
	}
	else
	{
		M3DQuadPrimitive* tempAtom = dynamic_cast<M3DQuadPrimitive*>(tempFigure->getPrimitivePtr(ubase-1,vbase+1));
		du12 = (x22 - tempAtom->getX()) / 2;

		if (side == 0)
			dU0du12 = (U0_22 - tempAtom->getU0()) / 2;
		else
			dU0du12 = (U0_22 - tempAtom->getU1()) / 2;
	}

	if (vbase + 1 == tempFigure->getColumnCount() - 1)
	{
		dv12 = x12 - x11;
		dU0dv12 = U0_12 - U0_11;
	}
	else
	{
		M3DQuadPrimitive* tempAtom = dynamic_cast<M3DQuadPrimitive*>(tempFigure->getPrimitivePtr(ubase, vbase+2));
		dv12 = (tempAtom->getX() - x11) / 2;

		if (side == 0)
			dU0dv12 = (tempAtom->getU0() - U0_11) / 2;
		else
			dU0dv12 = (tempAtom->getU1() - U0_11) / 2;
	}

	// Corner 22
	
	if (ubase + 1 == tempFigure->getRowCount() - 1)
	{
		du22 = x22 - x12;
		dU0du22 = U0_22 - U0_12;
	}
	else
	{
		M3DQuadPrimitive* tempAtom = dynamic_cast<M3DQuadPrimitive*>(tempFigure->getPrimitivePtr(ubase+2,vbase+1));
		du22 = (tempAtom->getX() - x12) / 2;

		if (side == 0)
			dU0du22 = (tempAtom->getU0() - U0_12) / 2;
		else
			dU0du22 = (tempAtom->getU1() - U0_12) / 2;
	}

	if (vbase + 1 == tempFigure->getColumnCount() - 1)
	{
		dv22 = x22 - x21;
		dU0dv22 = U0_22 - U0_21;
	}
	else
	{
		M3DQuadPrimitive* tempAtom = dynamic_cast<M3DQuadPrimitive*>(tempFigure->getPrimitivePtr(ubase+1,vbase+2));
		dv22 = (tempAtom->getX() - x21) / 2;

		if (side == 0)
			dU0dv22 = (tempAtom->getU0() - U0_21) / 2;
		else
			dU0dv22 = (tempAtom->getU1() - U0_21) / 2;
	}

//	cout << endl << "U derivatives: " << endl;
//	cout << du11.getX() << ", " << du11.getY() << ", " << du11.getZ() << endl;
//	cout << du21.getX() << ", " << du21.getY() << ", " << du21.getZ() << endl;
//	cout << du12.getX() << ", " << du12.getY() << ", " << du12.getZ() << endl;
//	cout << du22.getX() << ", " << du22.getY() << ", " << du22.getZ() << endl;

//	cout << endl << "V derivatives: " << endl;
//	cout << dv11.getX() << ", " << dv11.getY() << ", " << dv11.getZ() << endl;
//	cout << dv21.getX() << ", " << dv21.getY() << ", " << dv21.getZ() << endl;
//	cout << dv12.getX() << ", " << dv12.getY() << ", " << dv12.getZ() << endl;
//	cout << dv22.getX() << ", " << dv22.getY() << ", " << dv22.getZ() << endl;


	/*Vector3D du11 = x21 - x11;
	Vector3D du21 = x21 - x11;
	Vector3D du12 = x22 - x12;
	Vector3D du22 = x22 - x12;

	Vector3D dv11 = x12 - x11;
	Vector3D dv21 = x22 - x21;
	Vector3D dv12 = x12 - x11;
	Vector3D dv22 = x22 - x21;*/

	

	/*Vector3D du11 = (x21 - x01) / 2;
	Vector3D du21 = (x31 - x11) / 2;
	Vector3D du12 = (x22 - x02) / 2;
	Vector3D du22 = (x32 - x12) / 2;

	Vector3D dv11 = (x12 - x10) / 2;
	Vector3D dv21 = (x22 - x20) / 2;
	Vector3D dv12 = (x13 - x11) / 2;
	Vector3D dv22 = (x23 - x21) / 2;*/

	// Get unit normals for the four corners and project derivatives on to the tangent plane
	Vector3D n11 = quadAtom11->getU0() - quadAtom11->getU1();
	n11.normalize();
	Vector3D n21 = quadAtom21->getU0() - quadAtom21->getU1();
	n21.normalize();
	Vector3D n12 = quadAtom12->getU0() - quadAtom22->getU1();
	n12.normalize();
	Vector3D n22 = quadAtom22->getU0() - quadAtom22->getU1();
	n22.normalize();

	Vector3D du11t = du11 - (du11 * n11) * n11;
	Vector3D du21t = du21 - (du21 * n21) * n21;
	Vector3D du12t = du12 - (du12 * n12) * n12;
	Vector3D du22t = du22 - (du22 * n22) * n22;

	Vector3D dv11t = dv11 - (dv11 * n11) * n11;
	Vector3D dv21t = dv21 - (dv21 * n21) * n21;
	Vector3D dv12t = dv12 - (dv12 * n12) * n12;
	Vector3D dv22t = dv22 - (dv22 * n22) * n22;

	// Build matrices for hermite interpolation of medial sheet
	double hx[16] = { x11.getX(), x21.getX(), du11t.getX(), du21t.getX(), 
		x12.getX(), x22.getX(), du12t.getX(), du22t.getX(), dv11t.getX(), 
		dv21t.getX(), 0, 0, dv21t.getX(), dv22t.getX(), 0, 0};

	double hy[16] = { x11.getY(), x21.getY(), du11t.getY(), du21t.getY(), 
		x12.getY(), x22.getY(), du12t.getY(), du22t.getY(), dv11t.getY(), 
		dv21t.getY(), 0, 0, dv21t.getY(), dv22t.getY(), 0, 0};

	double hz[16] = { x11.getZ(), x21.getZ(), du11t.getZ(), du21t.getZ(), 
		x12.getZ(), x22.getZ(),du12t.getZ(), du22t.getZ(), dv11t.getZ(), 
		dv21t.getZ(), 0, 0, dv21t.getZ(), dv22t.getZ(), 0, 0};

	Matrix hxmat = Matrix(4, 4, hx, true);
	Matrix hymat = Matrix(4, 4, hy, true);
	Matrix hzmat = Matrix(4, 4, hz, true);

	double hu[4] = { h1(u), h2(u), h3(u), h4(u) };
	double hv[4] = { h1(v), h2(v), h3(v), h4(v) };
	Matrix humat = Matrix(1, 4, hu, true);
	Matrix hvmat = Matrix(4, 1, hv, true);

	Matrix xn = humat * hxmat * hvmat;
	Matrix yn = humat * hymat * hvmat;
	Matrix zn = humat * hzmat * hvmat;

	// Calculate sRad matrices using finite differences
	/*Vector3D dU0du11 = (quadAtom21->getU0() - quadAtom01->getU0()) / 2;
	Vector3D dU0dv11 = (quadAtom12->getU0() - quadAtom10->getU0()) / 2;

	Vector3D dU0du21 = (quadAtom31->getU0() - quadAtom11->getU0()) / 2;
	Vector3D dU0dv21 = (quadAtom22->getU0() - quadAtom20->getU0()) / 2;

	Vector3D dU0du12 = (quadAtom22->getU0() - quadAtom02->getU0()) / 2;
	Vector3D dU0dv12 = (quadAtom12->getU0() - quadAtom10->getU0()) / 2;

	Vector3D dU0du22 = (quadAtom32->getU0() - quadAtom12->getU0()) / 2;
	Vector3D dU0dv22 = (quadAtom23->getU0() - quadAtom21->getU0()) / 2;*/

	/*Vector3D dU0du11, dU0dv11, dU0du21, dU0dv21, dU0du12, dU0dv12, dU0du22, dU0dv22;

	if (side == 0)
	{
		dU0du11 = quadAtom21->getU0() - quadAtom11->getU0();
		dU0dv11 = quadAtom12->getU0() - quadAtom11->getU0();

		dU0du21 = quadAtom21->getU0() - quadAtom11->getU0();
		dU0dv21 = quadAtom22->getU0() - quadAtom21->getU0();

		dU0du12 = quadAtom22->getU0() - quadAtom12->getU0();
		dU0dv12 = quadAtom12->getU0() - quadAtom11->getU0();

		dU0du22 = quadAtom22->getU0() - quadAtom12->getU0();
		dU0dv22 = quadAtom22->getU0() - quadAtom21->getU0();
	}
	else
	{
		dU0du11 = quadAtom21->getU1() - quadAtom11->getU1();
		dU0dv11 = quadAtom12->getU1() - quadAtom11->getU1();

		dU0du21 = quadAtom21->getU1() - quadAtom11->getU1();
		dU0dv21 = quadAtom22->getU1() - quadAtom21->getU1();

		dU0du12 = quadAtom22->getU1() - quadAtom12->getU1();
		dU0dv12 = quadAtom12->getU1() - quadAtom11->getU1();

		dU0du22 = quadAtom22->getU1() - quadAtom12->getU1();
		dU0dv22 = quadAtom22->getU1() - quadAtom21->getU1();
	}*/

	// These form the non-orthogonal medial coordinate system for the spoke derivatives

	//Vector3D U0_11, U0_12, U0_21, U0_22;

	//if (side == 0)
	//{
	//	U0_11 = quadAtom11->getU0();
	//	U0_21 = quadAtom21->getU0();
	//	U0_12 = quadAtom12->getU0();
	//	U0_22 = quadAtom22->getU0();
	//}
	//else
	//{
	//	U0_11 = quadAtom11->getU1();
	//	U0_21 = quadAtom21->getU1();
	//	U0_12 = quadAtom12->getU1();
	//	U0_22 = quadAtom22->getU1();
	//}

	double basis11[9] = { du11.getX(), du11.getY(), du11.getZ(), dv11.getX(), dv11.getY(), dv11.getZ(), 
		U0_11.getX(), U0_11.getY(), U0_11.getZ() };

	double basis21[9] = { du21.getX(), du21.getY(), du21.getZ(), dv21.getX(), dv21.getY(), dv21.getZ(), 
		U0_21.getX(), U0_21.getY(), U0_21.getZ() };

	double basis12[9] = { du12.getX(), du12.getY(), du12.getZ(), dv12.getX(), dv12.getY(), dv12.getZ(), 
		U0_12.getX(), U0_12.getY(), U0_12.getZ() };

	double basis22[9] = { du22.getX(), du22.getY(), du22.getZ(), dv22.getX(), dv22.getY(), dv22.getZ(), 
		U0_22.getX(), U0_22.getY(), U0_22.getZ() };
	


	Matrix C11 = Matrix(3, 3, basis11, true);
	Matrix C21 = Matrix(3, 3, basis21, true);
	Matrix C12 = Matrix(3, 3, basis12, true);
	Matrix C22 = Matrix(3, 3, basis22, true);

	Matrix Ci11, Ci21, Ci12, Ci22;
	
	C11.inverse(Ci11);
	C21.inverse(Ci21);
	C12.inverse(Ci12);
	C22.inverse(Ci22);

	double dU0du11_coeffs[3] = { dU0du11.getX(), dU0du11.getY(), dU0du11.getZ() };
	double dU0dv11_coeffs[3] = { dU0dv11.getX(), dU0dv11.getY(), dU0dv11.getZ() };
	double dU0du21_coeffs[3] = { dU0du21.getX(), dU0du21.getY(), dU0du21.getZ() };
	double dU0dv21_coeffs[3] = { dU0dv21.getX(), dU0dv21.getY(), dU0dv21.getZ() };
	double dU0du12_coeffs[3] = { dU0du12.getX(), dU0du12.getY(), dU0du12.getZ() };
	double dU0dv12_coeffs[3] = { dU0dv12.getX(), dU0dv12.getY(), dU0dv12.getZ() };
	double dU0du22_coeffs[3] = { dU0du22.getX(), dU0du22.getY(), dU0du22.getZ() };
	double dU0dv22_coeffs[3] = { dU0dv22.getX(), dU0dv22.getY(), dU0dv22.getZ() };

	Matrix Au11 = Matrix(3,1,dU0du11_coeffs,true);
	Matrix Av11 = Matrix(3,1,dU0dv11_coeffs,true);
	Matrix Au21 = Matrix(3,1,dU0dv21_coeffs,true);
	Matrix Av21 = Matrix(3,1,dU0dv21_coeffs,true);
	Matrix Au12 = Matrix(3,1,dU0du12_coeffs,true);
	Matrix Av12 = Matrix(3,1,dU0dv12_coeffs,true);
	Matrix Au22 = Matrix(3,1,dU0du22_coeffs,true);
	Matrix Av22 = Matrix(3,1,dU0dv22_coeffs,true);

	Matrix Bu11 = Ci11 * Au11;
	Matrix Bv11 = Ci11 * Av11;
	Matrix Bu21 = Ci21 * Au21;
	Matrix Bv21 = Ci21 * Av21;
	Matrix Bu12 = Ci12 * Au12;
	Matrix Bv12 = Ci12 * Av12;
	Matrix Bu22 = Ci22 * Au22;
	Matrix Bv22 = Ci22 * Av22;

	double b11[4] = { -1*Bu11(0,0), -1*Bu11(1,0), -1*Bv11(0,0), -1*Bv11(1,0) };
	double b21[4] = { -1*Bu21(0,0), -1*Bu21(1,0), -1*Bv21(0,0), -1*Bv21(1,0) };
	double b12[4] = { -1*Bu12(0,0), -1*Bu12(1,0), -1*Bv12(0,0), -1*Bv12(1,0) };
	double b22[4] = { -1*Bu22(0,0), -1*Bu22(1,0), -1*Bv22(0,0), -1*Bv22(1,0) };

	Matrix Srad11 = Matrix(2,2,b11,true);
	Matrix Srad21 = Matrix(2,2,b21,true);
	Matrix Srad12 = Matrix(2,2,b12,true);
	Matrix Srad22 = Matrix(2,2,b22,true);

	Matrix rSrad11 = Srad11 * quadAtom11->getR();
	Matrix rSrad21 = Srad21 * quadAtom21->getR();
	Matrix rSrad12 = Srad12 * quadAtom12->getR();
	Matrix rSrad22 = Srad22 * quadAtom22->getR();

//	cout << endl << "rSrad: " << endl;
//	cout << rSrad11(0,0) << ", " << rSrad11(0,1) << ", " << rSrad11(1,0) << ", " << rSrad11(1,1) << endl;
//	cout << rSrad21(0,0) << ", " << rSrad21(0,1) << ", " << rSrad21(1,0) << ", " << rSrad21(1,1) << endl;
//	cout << rSrad12(0,0) << ", " << rSrad12(0,1) << ", " << rSrad12(1,0) << ", " << rSrad12(1,1) << endl;
//	cout << rSrad22(0,0) << ", " << rSrad22(0,1) << ", " << rSrad22(1,0) << ", " << rSrad22(1,1) << endl;

	Vector L11, L21, L12, L22;
	Matrix V11, V21, V12, V22;

	rSrad11.factorEV(L11, V11, NON_SYM);
	rSrad21.factorEV(L21, V21, NON_SYM);
	rSrad12.factorEV(L12, V12, NON_SYM);
	rSrad22.factorEV(L22, V22, NON_SYM);

	//cout << "HERE" << endl;

	double Lam1_11, Lam2_11, Lam1_21, Lam2_21, Lam1_12, Lam2_12, Lam1_22, Lam2_22;
	Vector e1_11, e2_11, e1_21, e2_21, e1_12, e2_12, e1_22, e2_22;

	// Corner 1
	if (L11(0) > L11(1))
	{
		Lam1_11 = L11(0);
		Lam2_11 = L11(1);

		e1_11 = V11.getColumn(0);
		e2_11 = V11.getColumn(1);
	}
	else
	{
		Lam1_11 = L11(1);
		Lam2_11 = L11(0);

		e1_11 = V11.getColumn(1);
		e2_11 = V11.getColumn(0);
	}

	// Corner 2
	if (L21(0) > L21(1))
	{
		Lam1_21 = L21(0);
		Lam2_21 = L21(1);

		e1_21 = V21.getColumn(0);
		e2_21 = V21.getColumn(1);
	}
	else
	{
		Lam1_21 = L21(1);
		Lam2_21 = L21(0);

		e1_21 = V21.getColumn(1);
		e2_21 = V21.getColumn(0);
	}

	//Corner 3
	if (L12(0) > L12(1))
	{
		Lam1_12 = L12(0);
		Lam2_12 = L12(1);

		e1_12 = V12.getColumn(0);
		e2_12 = V12.getColumn(1);
	}
	else
	{
		Lam1_12 = L12(1);
		Lam2_12 = L12(0);

		e1_12 = V12.getColumn(1);
		e2_12 = V12.getColumn(0);
	}

	//Corner 4
	if (L22(0) > L22(1))
	{
		Lam1_22 = L22(0);
		Lam2_22 = L22(1);

		e1_22 = V22.getColumn(0);
		e2_22 = V22.getColumn(1);
	}
	else
	{
		Lam1_22 = L22(1);
		Lam2_22 = L22(0);

		e1_22 = V22.getColumn(1);
		e2_22 = V22.getColumn(0);
	}

//	cout << endl << "Eigenvectors: " << endl;
//	cout << "(" << e1_11(0) << ", " << e1_11(1) << "), (" << e2_11(0) << ", " << e2_11(1) << ")" << endl;
//	cout << "(" << e1_21(0) << ", " << e1_21(1) << "), (" << e2_21(0) << ", " << e2_21(1) << ")" << endl;
//	cout << "(" << e1_12(0) << ", " << e1_12(1) << "), (" << e2_12(0) << ", " << e2_12(1) << ")" << endl;
//	cout << "(" << e1_22(0) << ", " << e1_22(1) << "), (" << e2_22(0) << ", " << e2_22(1) << ")" << endl;

	double testlambda[4] = { Lam1_11, 0, 0, Lam2_11 };

	
	
	Matrix testL = Matrix(2,2,testlambda,true);
	Matrix testinv;

	Matrix testsrad = V11 * testL * V11.inverse(testinv);

//	cout << endl << endl;
//	cout << testsrad(0,0) << ", " << testsrad(0,1) << ", " << testsrad(1,0) << ", " << testsrad(1,1) << endl;
//	cout << endl << endl;

//	cout << endl << "Eigenvalues: " << endl;
//	cout << Lam1_11 << ", " << Lam2_11 << endl;
//	cout << Lam1_21 << ", " << Lam2_21 << endl;
//	cout << Lam1_12 << ", " << Lam2_12 << endl;
//	cout << Lam1_22 << ", " << Lam2_22 << endl;


	double logLam1_11 = log(1 - Lam1_11);
	double logLam2_11 = log(1 - Lam2_11);
	double logLam1_21 = log(1 - Lam1_21);
	double logLam2_21 = log(1 - Lam2_21);
	double logLam1_12 = log(1 - Lam1_12);
	double logLam2_12 = log(1 - Lam2_12);
	double logLam1_22 = log(1 - Lam1_22);
	double logLam2_22 = log(1 - Lam2_22);

	double logAvg1 = (1-u)*(1-v)*logLam1_11 + (u)*(1-v)*logLam1_21 + (1-u)*(v)*logLam1_12 + (u)*(v)*logLam1_22;
	double logAvg2 = (1-u)*(1-v)*logLam2_11 + (u)*(1-v)*logLam2_21 + (1-u)*(v)*logLam2_12 + (u)*(v)*logLam2_22;

	double Lam1 = 1 - exp(logAvg1);
	double Lam2 = 1 - exp(logAvg2);

	// Uses a form of bilinear interpolation to get the eigenvectors: This should probably be a Frechet mean of the thetas 
	// on the unit circle in the future

	/*double theta1_11 = atan2(e1_11(1), e1_11(0));
	double theta2_11 = atan2(e2_11(1), e2_11(0));
	double theta1_21 = atan2(e1_21(1), e1_21(0));
	double theta2_21 = atan2(e2_21(1), e2_21(0));
	double theta1_12 = atan2(e1_12(1), e1_12(0));
	double theta2_12 = atan2(e2_12(1), e2_12(0));
	double theta1_22 = atan2(e1_22(1), e1_22(0));
	double theta2_22 = atan2(e2_22(1), e2_22(0));*/

	double PI = 3.14159265;

	double avgx1, avgy1, avgx2, avgy2, avgtheta1, avgtheta2;//, at1_11_21, at2_11_21, at1_12_22, at2_12_22;

	avgx1 = (1-u)*(1-v)*e1_11(0) + (u)*(1-v)*e1_21(0) + (1-u)*(v)*e1_12(0) + (u)*(v)*e1_22(0);
	avgy1 = (1-u)*(1-v)*e1_11(1) + (u)*(1-v)*e1_21(1) + (1-u)*(v)*e1_12(1) + (u)*(v)*e1_22(1);
	avgtheta1 = atan2(avgy1, avgx1);

	avgx2 = (1-u)*(1-v)*e2_11(0) + (u)*(1-v)*e2_21(0) + (1-u)*(v)*e2_12(0) + (u)*(v)*e2_22(0);
	avgy2 = (1-u)*(1-v)*e2_11(1) + (u)*(1-v)*e2_21(1) + (1-u)*(v)*e2_12(1) + (u)*(v)*e2_22(1);
	avgtheta2 = atan2(avgy2, avgx2);

	// Go from 11 to 21 (u direction)
	//double thetadist1_11_21 = theta1_21 - theta1_11;
	//if (thetadist1_11_21 > PI)
	//{
	//	theta1_21 -= 2*PI;
	//	thetadist1_11_21 = theta1_21 - theta1_11;
	//}
	//else if (thetadist1_11_21 < -PI)
	//{
	//	theta1_21 += 2*PI;
	//	thetadist1_11_21 = theta1_21 - theta1_11;
	//}
	//at1_11_21 = theta1_11 + u*thetadist1_11_21;

	//double thetadist2_11_21 = theta2_21 - theta2_11;
	//if (thetadist2_11_21 > PI)
	//{
	//	theta2_21 -= 2*PI;
	//	thetadist2_11_21 = theta2_21 - theta2_11;
	//}
	//else if (thetadist2_11_21 < -PI)
	//{
	//	theta2_21 += 2*PI;
	//	thetadist2_11_21 = theta2_21 - theta2_11;
	//}
	//at2_11_21 = theta2_11 + u*thetadist2_11_21;

	//// Go from 12 to 22 (u direction)
	//double thetadist1_12_22 = theta1_22 - theta1_12;
	//if (thetadist1_12_22 > PI)
	//{
	//	theta1_22 -= 2*PI;
	//	thetadist1_12_22 = theta1_22 - theta1_12;
	//}
	//else if (thetadist1_12_22 < -PI)
	//{
	//	theta1_22 += 2*PI;
	//	thetadist1_12_22 = theta1_22 - theta1_12;
	//}
	//at1_12_22 = theta1_12 + u*thetadist1_12_22;

	//double thetadist2_12_22 = theta2_22 - theta2_12;
	//if (thetadist2_12_22 > PI)
	//{
	//	theta2_22 -= 2*PI;
	//	thetadist2_12_22 = theta2_22 - theta2_12;
	//}
	//else if (thetadist2_12_22 < -PI)
	//{
	//	theta2_22 += 2*PI;
	//	thetadist2_12_22 = theta2_22 - theta2_12;
	//}
	//at2_12_22 = theta2_12 + u*thetadist2_12_22;

	//// Now go v direction
	//double thetadist1_v = at1_12_22 - at1_11_21;
	//if (thetadist1_v > PI)
	//{
	//	at1_12_22 -= 2*PI;
	//	thetadist1_v = at1_12_22 - at1_11_21;
	//}
	//else if (thetadist1_v < -PI)
	//{
	//	at1_12_22 += 2*PI;
	//	thetadist1_v = at1_12_22 - at1_11_21;
	//}
	//avgtheta1 = at1_11_21 + v*thetadist1_v;

	//double thetadist2_v = at2_12_22 - at2_11_21;
	//if (thetadist2_v > PI)
	//{
	//	at2_12_22 -= 2*PI;
	//	thetadist2_v = at2_12_22 - at2_11_21;
	//}
	//else if (thetadist2_v < -PI)
	//{
	//	at2_12_22 += 2*PI;
	//	thetadist2_v = at2_12_22 - at2_11_21;
	//}
	//avgtheta2 = at2_11_21 + v*thetadist2_v;

	double neweigenv[4] = { cos(avgtheta1), sin(avgtheta1), cos(avgtheta2), sin(avgtheta2) };

	Matrix NewV = Matrix(2,2,neweigenv,true);

	double newlambda[4] = { Lam1, 0, 0, Lam2 };

	Matrix NewL = Matrix(2,2,newlambda,true);
	
	/*cout << V11(0,0) << ", " << V11(1,0)  << ", " << V11(0,1) << ", " << V11(1,1) << endl;
	cout << V21(0,0) << ", " << V21(1,0)  << ", " << V21(0,1) << ", " << V21(1,1) << endl;
	cout << V12(0,0) << ", " << V12(1,0)  << ", " << V12(0,1) << ", " << V12(1,1) << endl;
	cout << V22(0,0) << ", " << V22(1,0)  << ", " << V22(0,1) << ", " << V22(1,1) << endl;
	cout << NewV(0,0) << ", " << NewV(1,0) << ", " << NewV(0,1) << ", " << NewV(1,1) << endl;*/
	
	Matrix NewVi;
	NewV.inverse(NewVi);

	Matrix NewrSrad = NewV * NewL * NewVi;

	double curru = 0;
	double currv = 0;
	double du = u / 10;
	double dv = v / 10;

	if (side == 0)
		U0_11 = quadAtom11->getU0();
	else
		U0_11 = quadAtom11->getU1();
	
	double Pdata[6] = { du11.getX(), dv11.getX(), du11.getY(), dv11.getY(), du11.getZ(), dv11.getZ() };
	Matrix P = Matrix(2,3,Pdata,true);

	double Udata[3] = { U0_11.getX(), U0_11.getY(), U0_11.getZ() };
	Matrix U = Matrix(1,3,Udata,true);

	double ru11 = du11 * U0_11;
	double rv11 = dv11 * U0_11;

	double ru21 = du21 * U0_21;
	double rv21 = dv21 * U0_21;

	double ru12 = du12 * U0_12;
	double rv12 = dv12 * U0_12;

	double ru22 = du22 * U0_22;
	double rv22 = dv22 * U0_22;

//	cout << endl << "Ru & Rv: " << endl;
//	cout << ru11 << ", " << rv11 << endl;
//	cout << ru21 << ", " << rv21 << endl;
//	cout << ru12 << ", " << rv12 << endl;
//	cout << ru22 << ", " << rv22 << endl;

	double Idata[9] = {1,0,0,0,1,0,0,0,1};
	Matrix I = Matrix(3,3,Idata,true);

	Matrix Q = P * ( (U.t() * U) - I );

	Matrix R = -1 * P * U.t();

	Matrix dS = ( rSrad11.t() * Q ) + ( R * U );

	Vector dSdu = dS.getRow(0);
	Vector dSdv = dS.getRow(1);

	Vector3D dSu = Vector3D(dSdu(0), dSdu(1), dSdu(2));
	Vector3D dSv = Vector3D(dSdv(0), dSdv(1), dSdv(2));

	Vector3D S0;

	if (side == 0)
		S0 = U0_11 * quadAtom11->getR0();
	else 
		S0 = U0_11 * quadAtom11->getR1();

	Vector3D NewSpoke = S0 + (du * dSu) + (dv * dSv);

	curru = curru + du;
	currv = currv + dv;

	hu[0] = h1(curru); hu[1] = h2(curru); hu[2] = h3(curru); hu[3] = h4(curru); 
	hv[0] = h1(currv); hv[1] = h2(currv); hv[2] = h3(currv); hv[3] = h4(currv);	

	humat = Matrix(1,4,hu,true);
	hvmat = Matrix(4,1,hv,true);

	xn = humat * hxmat * hvmat;
	yn = humat * hymat * hvmat;
	zn = humat * hzmat * hvmat;

	Vector3D newpos = Vector3D(xn(0,0), yn(0,0), zn(0,0));

	for (int i = 0; i < 10; i++)
	{
		// To find rSrad, need to calculate lambdas and eigenvectors
		logAvg1 = (1-curru)*(1-currv)*logLam1_11 + (curru)*(1-currv)*logLam1_21 + (1-curru)*(currv)*logLam1_12 + (curru)*(currv)*logLam1_22;
		logAvg2 = (1-curru)*(1-currv)*logLam2_11 + (curru)*(1-currv)*logLam2_21 + (1-curru)*(currv)*logLam2_12 + (curru)*(currv)*logLam2_22;
	
		Lam1 = 1 - exp(logAvg1);
		Lam2 = 1 - exp(logAvg2);


		avgx1 = (1-curru)*(1-currv)*e1_11(0) + (curru)*(1-currv)*e1_21(0) + (1-curru)*(currv)*e1_12(0) + (curru)*(currv)*e1_22(0);
		avgy1 = (1-curru)*(1-currv)*e1_11(1) + (curru)*(1-currv)*e1_21(1) + (1-curru)*(currv)*e1_12(1) + (curru)*(currv)*e1_22(1);
		avgtheta1 = atan2(avgy1, avgx1);

		avgx2 = (1-u)*(1-v)*e2_11(0) + (u)*(1-v)*e2_21(0) + (1-u)*(v)*e2_12(0) + (u)*(v)*e2_22(0);
		avgy2 = (1-u)*(1-v)*e2_11(1) + (u)*(1-v)*e2_21(1) + (1-u)*(v)*e2_12(1) + (u)*(v)*e2_22(1);
		avgtheta2 = atan2(avgy2, avgx2);

		neweigenv[0] = cos(avgtheta1); neweigenv[1] = sin(avgtheta1); neweigenv[2] = cos(avgtheta2); neweigenv[3] = sin(avgtheta2);
		NewV = Matrix(2,2,neweigenv,true);
		NewV.inverse(NewVi);

		newlambda[0] = Lam1; newlambda[1] = 0; newlambda[2] = 0; newlambda[3] = Lam2;
		NewL = Matrix(2,2,newlambda,true);

		NewrSrad = NewV * NewL * NewVi;

		double hu[4] = { h1(curru), h2(curru), h3(curru), h4(curru) };
		double hv[4] = { h1(currv), h2(currv), h3(currv), h4(currv) };
		Matrix humat = Matrix(1,4,hu,true);
		Matrix hvmat = Matrix(4,1,hv,true);

		double hup[4] = { h1p(curru), h2p(curru), h3p(curru), h4p(curru) };
		double hvp[4] = { h1p(currv), h2p(currv), h3p(currv), h4p(currv) };
		Matrix hupmat = Matrix(1,4,hup,true);
		Matrix hvpmat = Matrix(4,1,hvp,true);

		Matrix pux = hupmat * hxmat * hvmat;
		Matrix puy = hupmat * hymat * hvmat;
		Matrix puz = hupmat * hzmat * hvmat;

		Matrix pvx = humat * hxmat * hvpmat;
		Matrix pvy = humat * hymat * hvpmat;
		Matrix pvz = humat * hzmat * hvpmat;
	
		double pdata[6] = { pux(0,0), pvx(0,0), puy(0,0), pvy(0,0), puz(0,0), pvz(0,0) };
		P = Matrix(2,3,pdata,true);

		double r = NewSpoke.normalize();

		double udata[3] = { NewSpoke.getX(), NewSpoke.getY(), NewSpoke.getZ() };
		U = Matrix(1,3,udata,true);

		Q = P * ( (U.t() * U) - I );
		
		R = -1 * P * U.t();
		
		dS = ( rSrad11.t() * Q ) + ( R * U );

		dSdu = dS.getRow(0);
		dSdv = dS.getRow(1);

		dSu = Vector3D(dSdu(0), dSdu(1), dSdu(2));
		dSv = Vector3D(dSdv(0), dSdv(1), dSdv(2));

		Vector3D TempSpoke = r * NewSpoke;

		NewSpoke = TempSpoke + (du * dSu) + (dv * dSv);

		curru = curru + du;
		currv = currv + dv;

	}

	hu[0] = h1(u); hu[1] = h2(u); hu[2] = h3(u); hu[3] = h4(u); 
	hv[0] = h1(v); hv[1] = h2(v); hv[2] = h3(v); hv[3] = h4(v);	

	humat = Matrix(1,4,hu,true);
	hvmat = Matrix(4,1,hv,true);

	xn = humat * hxmat * hvmat;
	yn = humat * hymat * hvmat;
	zn = humat * hzmat * hvmat;

	newpos = Vector3D(xn(0,0), yn(0,0), zn(0,0));

	//cout << "{" << quadAtom11->getX().getX() + (quadAtom11->getR0() * quadAtom11->getU0().getX()) << ", " << quadAtom11->getX().getY() + (quadAtom11->getR0() * quadAtom11->getU0().getY()) << ", " << quadAtom11->getX().getZ() + (quadAtom11->getR0() * quadAtom11->getU0().getZ()) << "}" << endl;
	//cout << "{" << quadAtom21->getX().getX() + (quadAtom21->getR0() * quadAtom21->getU0().getX()) << ", " << quadAtom21->getX().getY() + (quadAtom21->getR0() * quadAtom21->getU0().getY()) << ", " << quadAtom21->getX().getZ() + (quadAtom21->getR0() * quadAtom21->getU0().getZ()) << "}" << endl;
	//cout << "{" << quadAtom12->getX().getX() + (quadAtom12->getR0() * quadAtom12->getU0().getX()) << ", " << quadAtom12->getX().getY() + (quadAtom12->getR0() * quadAtom12->getU0().getY()) << ", " << quadAtom12->getX().getZ() + (quadAtom12->getR0() * quadAtom12->getU0().getZ()) << "}" << endl;
	//cout << "{" << quadAtom22->getX().getX() + (quadAtom22->getR0() * quadAtom22->getU0().getX()) << ", " << quadAtom22->getX().getY() + (quadAtom22->getR0() * quadAtom22->getU0().getY()) << ", " << quadAtom22->getX().getZ() + (quadAtom22->getR0() * quadAtom22->getU0().getZ()) << "}" << endl;
	//cout << "{" << newpos.getX() + NewSpoke.getX() << ", " << newpos.getY() + NewSpoke.getY() << ", " << newpos.getZ() + NewSpoke.getZ() << endl;

	//cout << quadAtom11->getR0() << endl;
	//cout << quadAtom21->getR0() << endl;
	//cout << quadAtom12->getR0() << endl;
	//cout << quadAtom22->getR0() << endl;
	//cout << NewSpoke.norm() << endl;

	//cout << quadAtom11->getX().getX() << ", " << quadAtom11->getX().getY() << ", " << quadAtom11->getX().getZ() << endl;
	//cout << quadAtom21->getX().getX() << ", " << quadAtom21->getX().getY() << ", " << quadAtom21->getX().getZ() << endl;
	//cout << quadAtom12->getX().getX() << ", " << quadAtom12->getX().getY() << ", " << quadAtom12->getX().getZ() << endl;
	//cout << quadAtom22->getX().getX() << ", " << quadAtom22->getX().getY() << ", " << quadAtom22->getX().getZ() << endl;

	double NewR = NewSpoke.normalize();

	//cout << newpos.getX() << ", " << newpos.getY() << ", " << newpos.getZ() << endl;

	//cout << newpos.getX() << ", " << NewR << ", " << NewSpoke.getX() << endl;

	M3DSpoke* spoke = new M3DSpoke(newpos, NewSpoke, NewR);
	

	return spoke;
}

vector<M3DSpoke*> Match::getRelevantSpokes(M3DObject* targetObject, int level, int figureId)
{

    /*vector<M3DSpoke*> spokes;
        vtkSmartPointer<vtkSRep> srepfig = vtkSmartPointer<vtkSRep>::New();
        GetVtkSrepFig(targetObject, figureId, srepfig);
        GetInternalSpokes(srepfig, level, spokes);
        GetCrestSpokes(srepfig, level, spokes);

        return spokes;*/
    M3DFigure * figure = targetObject->getFigurePtr( figureId ) ;
    M3DQuadFigure* quadFig = dynamic_cast<M3DQuadFigure*>( figure );

    if(quadFig){

        int umax = quadFig->getRowCount() - 1;
        int vmax = quadFig->getColumnCount() - 1;

        vector<double> steps;

        switch (level)
        {
        case 3:
            steps.push_back(.125);
            steps.push_back(.375);
            steps.push_back(.625);
            steps.push_back(.875);
        case 2:
            steps.push_back(.25);
            steps.push_back(.75);
        case 1:
            steps.push_back(.5);
            break;
        }

        vector<M3DSpoke*> spokes;


        int u,v,i, j;
        for (u = 0; u < umax; u++)
        {

            for (v = 0; v < vmax; v++)
            {
                M3DQuadPrimitive* quadAtom = dynamic_cast<M3DQuadPrimitive*>( quadFig->getPrimitivePtr(u, v) );
                Vector3D X = quadAtom->getX();
                Vector3D U0 = quadAtom->getU0();
                Vector3D U1 = quadAtom->getU1();
                double R0 = quadAtom->getR0();
                double R1 = quadAtom->getR1();

                M3DSpoke* top = new M3DSpoke(X, U0, R0);
                M3DSpoke* bottom = new M3DSpoke(X, U1, R1);

                spokes.push_back(top);
                spokes.push_back(bottom);

                for (i = 0; i < steps.size(); i++)
                {
                    for(j = 0; j < steps.size(); j++)
                    {
                        spokes.push_back(interpolateSpoke(figure, u + steps[i], v + steps[j], 0));
                        spokes.push_back(interpolateSpoke(figure, u + steps[i], v + steps[j], 1));

                        spokes.push_back(interpolateSpoke(figure, u, v + steps[j], 0));
                        spokes.push_back(interpolateSpoke(figure, u, v + steps[j], 0));

                        spokes.push_back(interpolateSpoke(figure, u + steps[i], v, 1));
                    }
                }
            }
        }



        vtkSmartPointer<vtkSRep> srepfig = vtkSmartPointer<vtkSRep>::New();
        GetVtkSrepFig(targetObject, figureId, srepfig);
        GetCrestSpokes(srepfig, level, spokes);

        return spokes;
    }else{

        vector<M3DSpoke*> spokes;
        vtkSmartPointer<vtkSRep> srepfig = vtkSmartPointer<vtkSRep>::New();
        GetVtkSrepFig(targetObject, figureId, srepfig);
        GetCrestSpokes(srepfig, level, spokes, true);

        return spokes;

    }
}

vector<M3DSpoke*> Match::getRelevantSpokes(M3DObject* targetObject, int level, int figureId, int atomId)
{
    /*vector<M3DSpoke*> spokes;
        vtkSmartPointer<vtkSRep> srepfig = vtkSmartPointer<vtkSRep>::New();
        vtkIdType vtkAtomId = GetVtkSrepFig(targetObject, figureId, srepfig, atomId);
        GetInternalSpokes(srepfig, level, spokes, vtkAtomId);
        GetCrestSpokes(srepfig, level, spokes, vtkAtomId);

        return spokes;*/
    M3DQuadFigure* figure = dynamic_cast<M3DQuadFigure*>( targetObject->getFigurePtr( figureId ) );
    if(figure){
        M3DQuadPrimitive* quadAtom = dynamic_cast<M3DQuadPrimitive*>( figure->getPrimitivePtr( atomId ) );

        int umax = figure->getRowCount() - 1;
        int vmax = figure->getColumnCount() - 1;


        int numCols = figure->getColumnCount();

        int u = atomId / numCols;
        int v = atomId % numCols;

        vector<double> steps;

        switch (level)
        {
        case 3:
            steps.push_back(.125);
            steps.push_back(.375);
            steps.push_back(.625);
            steps.push_back(.875);
        case 2:
            steps.push_back(.25);
            steps.push_back(.75);
        case 1:
            steps.push_back(.5);
            break;
        }

        vector<M3DSpoke*> spokes;

        Vector3D X = quadAtom->getX();
        Vector3D U0 = quadAtom->getU0();
        Vector3D U1 = quadAtom->getU1();
        double R0 = quadAtom->getR0();
        double R1 = quadAtom->getR1();

        M3DSpoke* top = new M3DSpoke(X, U0, R0);
        M3DSpoke* bottom = new M3DSpoke(X, U1, R1);

        spokes.push_back(top);
        spokes.push_back(bottom);

        int i, j;
        for (i = 0; i < steps.size(); i++)
        {

            for(j = 0; j < steps.size(); j++)
            {
                if (u == 0)
                {
                    if (v == 0)
                    {
                        spokes.push_back(interpolateSpoke(figure, u + steps[i], v + steps[j], 0));
                        spokes.push_back(interpolateSpoke(figure, u + steps[i], v + steps[j], 1));

                        spokes.push_back(interpolateSpoke(figure, u + steps[i], v, 0));
                        spokes.push_back(interpolateSpoke(figure, u + steps[i], v, 1));

                        spokes.push_back(interpolateSpoke(figure, u, v + steps[i], 0));
                        spokes.push_back(interpolateSpoke(figure, u, v + steps[i], 1));
                    }
                    else if (v == vmax)
                    {
                        spokes.push_back(interpolateSpoke(figure, u + steps[i], v - steps[j], 0));
                        spokes.push_back(interpolateSpoke(figure, u + steps[i], v - steps[j], 1));


                        //spokes.push_back(interpolateSpoke(figure, u + steps[i], v, 0));
                        //spokes.push_back(interpolateSpoke(figure, u + steps[i], v, 1));


                        spokes.push_back(interpolateSpoke(figure, u, v - steps[i], 0));
                        spokes.push_back(interpolateSpoke(figure, u, v - steps[i], 1));
                    }
                    else
                    {
                        spokes.push_back(interpolateSpoke(figure, u + steps[i], v + steps[j], 0));
                        spokes.push_back(interpolateSpoke(figure, u + steps[i], v + steps[j], 1));

                        spokes.push_back(interpolateSpoke(figure, u + steps[i], v - steps[j], 0));
                        spokes.push_back(interpolateSpoke(figure, u + steps[i], v - steps[j], 1));

                        spokes.push_back(interpolateSpoke(figure, u + steps[i], v, 0));
                        spokes.push_back(interpolateSpoke(figure, u + steps[i], v, 1));

                        spokes.push_back(interpolateSpoke(figure, u, v + steps[j], 0));
                        spokes.push_back(interpolateSpoke(figure, u, v + steps[j], 1));

                        spokes.push_back(interpolateSpoke(figure, u, v - steps[j], 0));
                        spokes.push_back(interpolateSpoke(figure, u, v - steps[j], 1));

                    }
                }
                else if (u == umax)
                {
                    if (v == 0)
                    {
                        spokes.push_back(interpolateSpoke(figure, u - steps[i], v + steps[j], 0));
                        spokes.push_back(interpolateSpoke(figure, u - steps[i], v + steps[j], 1));

                        spokes.push_back(interpolateSpoke(figure, u - steps[i], v, 0));
                        spokes.push_back(interpolateSpoke(figure, u - steps[i], v, 1));

                        //spokes.push_back(interpolateSpoke(figure, u, v + steps[i], 0));
                        //spokes.push_back(interpolateSpoke(figure, u, v + steps[i], 1));

                    }
                    else if (v == vmax)
                    {
                        spokes.push_back(interpolateSpoke(figure, u - steps[i], v - steps[j], 0));
                        spokes.push_back(interpolateSpoke(figure, u - steps[i], v - steps[j], 1));

                        //spokes.push_back(interpolateSpoke(figure, u - steps[i], v, 0));
                        //spokes.push_back(interpolateSpoke(figure, u - steps[i], v, 1));

                        //spokes.push_back(interpolateSpoke(figure, u, v - steps[i], 0));
                        //spokes.push_back(interpolateSpoke(figure, u, v - steps[i], 1));
                    }
                    else
                    {
                        spokes.push_back(interpolateSpoke(figure, u - steps[i], v + steps[j], 1));
                        spokes.push_back(interpolateSpoke(figure, u - steps[i], v + steps[j], 0));

                        spokes.push_back(interpolateSpoke(figure, u - steps[i], v - steps[j], 0));
                        spokes.push_back(interpolateSpoke(figure, u - steps[i], v - steps[j], 1));

                        spokes.push_back(interpolateSpoke(figure, u - steps[i], v, 0));
                        spokes.push_back(interpolateSpoke(figure, u - steps[i], v, 1));

                        //spokes.push_back(interpolateSpoke(figure, u, v + steps[j], 0));
                        //spokes.push_back(interpolateSpoke(figure, u, v + steps[j], 1));

                        //spokes.push_back(interpolateSpoke(figure, u, v - steps[j], 0));
                        //spokes.push_back(interpolateSpoke(figure, u, v - steps[j], 1));
                    }
                }
                else
                {
                    if (v==0)
                    {
                        spokes.push_back(interpolateSpoke(figure, u + steps[i], v + steps[j], 0));
                        spokes.push_back(interpolateSpoke(figure, u + steps[i], v + steps[j], 1));

                        spokes.push_back(interpolateSpoke(figure, u - steps[i], v + steps[j], 0));
                        spokes.push_back(interpolateSpoke(figure, u - steps[i], v + steps[j], 1));

                        spokes.push_back(interpolateSpoke(figure, u + steps[i], v, 0));
                        spokes.push_back(interpolateSpoke(figure, u + steps[i], v, 1));

                        spokes.push_back(interpolateSpoke(figure, u - steps[i], v, 0));
                        spokes.push_back(interpolateSpoke(figure, u - steps[i], v, 1));

                        spokes.push_back(interpolateSpoke(figure, u, v + steps[i], 0));
                        spokes.push_back(interpolateSpoke(figure, u, v + steps[i], 1));

                    }
                    else if (v == vmax)
                    {
                        spokes.push_back(interpolateSpoke(figure, u - steps[i], v - steps[j], 0));
                        spokes.push_back(interpolateSpoke(figure, u - steps[i], v - steps[j], 1));

                        spokes.push_back(interpolateSpoke(figure, u + steps[i], v - steps[j], 0));
                        spokes.push_back(interpolateSpoke(figure, u + steps[i], v - steps[j], 1));

                        //spokes.push_back(interpolateSpoke(figure, u + steps[i], v, 0));
                        //spokes.push_back(interpolateSpoke(figure, u + steps[i], v, 1));

                        //spokes.push_back(interpolateSpoke(figure, u - steps[i], v, 0));
                        //spokes.push_back(interpolateSpoke(figure, u - steps[i], v, 1));

                        spokes.push_back(interpolateSpoke(figure, u, v - steps[i], 0));
                        spokes.push_back(interpolateSpoke(figure, u, v - steps[i], 1));

                    }
                    else
                    {
                        spokes.push_back(interpolateSpoke(figure, u + steps[i], v + steps[j], 0));
                        spokes.push_back(interpolateSpoke(figure, u + steps[i], v + steps[j], 1));

                        spokes.push_back(interpolateSpoke(figure, u - steps[i], v - steps[j], 0));
                        spokes.push_back(interpolateSpoke(figure, u - steps[i], v - steps[j], 1));

                        spokes.push_back(interpolateSpoke(figure, u + steps[i], v - steps[j], 0));
                        spokes.push_back(interpolateSpoke(figure, u + steps[i], v - steps[j], 1));

                        spokes.push_back(interpolateSpoke(figure, u - steps[i], v + steps[j], 0));
                        spokes.push_back(interpolateSpoke(figure, u - steps[i], v + steps[j], 1));

                        spokes.push_back(interpolateSpoke(figure, u + steps[i], v, 0));
                        spokes.push_back(interpolateSpoke(figure, u + steps[i], v, 1));

                        spokes.push_back(interpolateSpoke(figure, u - steps[i], v, 0));
                        spokes.push_back(interpolateSpoke(figure, u - steps[i], v, 1));

                        spokes.push_back(interpolateSpoke(figure, u, v + steps[i], 0));
                        spokes.push_back(interpolateSpoke(figure, u, v + steps[j], 1));

                        spokes.push_back(interpolateSpoke(figure, u, v - steps[i], 0));
                        spokes.push_back(interpolateSpoke(figure, u, v - steps[j], 1));

                    }
                }
            }
        }

        if(dynamic_cast<M3DQuadEndPrimitive*>(quadAtom) != 0){
            vtkSmartPointer<vtkSRep> srepfig = vtkSmartPointer<vtkSRep>::New();
            vtkIdType vtkAtomId = GetVtkSrepFig(targetObject, figureId, srepfig);
            GetCrestSpokes(srepfig, level, spokes,false, vtkAtomId);
        }

        return spokes;
    }else{

        vector<M3DSpoke*> spokes;
        vtkSmartPointer<vtkSRep> srepfig = vtkSmartPointer<vtkSRep>::New();
        vtkIdType vtkAtomId = GetVtkSrepFig(targetObject, figureId, srepfig);
        GetCrestSpokes(srepfig, level, spokes,true, vtkAtomId);

        return spokes;

    }
}

vector<M3DSpoke*> Match::getRelevantSpokes(M3DObject* targetObject, int level, int figureId, int atomId, int spokeId)//spoke id is top or bottom
{

    /*vector<M3DSpoke*> spokes;
        vtkSmartPointer<vtkSRep> srepfig = vtkSmartPointer<vtkSRep>::New();
        vtkIdType vtkAtomId = GetVtkSrepFig(targetObject, figureId, srepfig, atomId);
        GetInternalSpokes(srepfig, level, spokes, vtkAtomId, spokeId);
        GetCrestSpokes(srepfig, level, spokes, vtkAtomId, spokeId);

        return spokes;*/

    M3DQuadFigure* figure = dynamic_cast<M3DQuadFigure*>( targetObject->getFigurePtr( figureId ) );
    if(figure){
        M3DQuadPrimitive* quadAtom = dynamic_cast<M3DQuadPrimitive*>( figure->getPrimitivePtr( atomId ) );

        int umax = figure->getRowCount() - 1;
        int vmax = figure->getColumnCount() - 1;
        int numCols = figure->getColumnCount();

        int u = atomId / numCols;
        int v = atomId % numCols;

        vector<double> steps;

        switch (level)
        {
        case 3:
            steps.push_back(.125);
            steps.push_back(.375);
            steps.push_back(.625);
            steps.push_back(.875);
        case 2:
            steps.push_back(.25);
            steps.push_back(.75);
        case 1:
            steps.push_back(.5);
            break;
        }

        vector<M3DSpoke*> spokes;

        Vector3D X = quadAtom->getX();
        Vector3D U;
        double R;
        if (spokeId == 0)
        {
            U = quadAtom->getU0();
            R = quadAtom->getR0();
        }
        else if (spokeId == 1)
        {
            U = quadAtom->getU1();
            R = quadAtom->getR1();
        }
        else
        {
            U = dynamic_cast<M3DQuadEndPrimitive*>(quadAtom)->getUEnd();
            R = dynamic_cast<M3DQuadEndPrimitive*>(quadAtom)->getREnd();
        }

        M3DSpoke* thisSpoke = new M3DSpoke(X, U, R);
        spokes.push_back(thisSpoke);


        int i, j;

        if (spokeId == 0 || spokeId == 1)
        {

            for (i = 0; i < steps.size(); i++)
            {
                for(j = 0; j < steps.size(); j++)
                {
                    if (u == 0)
                    {
                        if (v == 0)
                        {
                            spokes.push_back(interpolateSpoke(figure, u + steps[i], v + steps[j], spokeId));
                            spokes.push_back(interpolateSpoke(figure, u + steps[i], v, spokeId));
                            spokes.push_back(interpolateSpoke(figure, u, v + steps[i], spokeId));
                        }
                        else if (v == vmax)
                        {
                            spokes.push_back(interpolateSpoke(figure, u + steps[i], v - steps[j], spokeId));
                            //spokes.push_back(interpolateSpoke(figure, u + steps[i], v, spokeId));
                            spokes.push_back(interpolateSpoke(figure, u, v - steps[i], spokeId));
                        }
                        else
                        {
                            spokes.push_back(interpolateSpoke(figure, u + steps[i], v + steps[j], spokeId));
                            spokes.push_back(interpolateSpoke(figure, u + steps[i], v - steps[j], spokeId));
                            spokes.push_back(interpolateSpoke(figure, u + steps[i], v, spokeId));
                            spokes.push_back(interpolateSpoke(figure, u, v + steps[j], spokeId));
                            spokes.push_back(interpolateSpoke(figure, u, v - steps[j], spokeId));
                        }
                    }
                    else if (u == umax)
                    {
                        if (v == 0)
                        {
                            spokes.push_back(interpolateSpoke(figure, u - steps[i], v + steps[j], spokeId));
                            spokes.push_back(interpolateSpoke(figure, u - steps[i], v, spokeId));
                            //spokes.push_back(interpolateSpoke(figure, u, v + steps[i], spokeId));

                        }
                        else if (v == vmax)
                        {
                            spokes.push_back(interpolateSpoke(figure, u - steps[i], v - steps[j], spokeId));
                            //spokes.push_back(interpolateSpoke(figure, u - steps[i], v, spokeId));
                            //spokes.push_back(interpolateSpoke(figure, u, v - steps[i], spokeId));
                        }
                        else
                        {
                            spokes.push_back(interpolateSpoke(figure, u - steps[i], v + steps[j], spokeId));
                            spokes.push_back(interpolateSpoke(figure, u - steps[i], v - steps[j], spokeId));
                            spokes.push_back(interpolateSpoke(figure, u - steps[i], v, spokeId));
                            //spokes.push_back(interpolateSpoke(figure, u, v + steps[j], spokeId));
                            //spokes.push_back(interpolateSpoke(figure, u, v - steps[j], spokeId));
                        }
                    }
                    else
                    {
                        if (v==0)
                        {
                            spokes.push_back(interpolateSpoke(figure, u + steps[i], v + steps[j], spokeId));
                            spokes.push_back(interpolateSpoke(figure, u - steps[i], v + steps[j], spokeId));
                            spokes.push_back(interpolateSpoke(figure, u + steps[i], v, spokeId));
                            spokes.push_back(interpolateSpoke(figure, u - steps[i], v, spokeId));
                            spokes.push_back(interpolateSpoke(figure, u, v + steps[i], spokeId));

                        }
                        else if (v == vmax)
                        {
                            spokes.push_back(interpolateSpoke(figure, u - steps[i], v - steps[j], spokeId));
                            spokes.push_back(interpolateSpoke(figure, u + steps[i], v - steps[j], spokeId));
                            //spokes.push_back(interpolateSpoke(figure, u + steps[i], v, spokeId));
                            //spokes.push_back(interpolateSpoke(figure, u - steps[i], v, spokeId));
                            spokes.push_back(interpolateSpoke(figure, u, v - steps[i], spokeId));

                        }
                        else
                        {
                            spokes.push_back(interpolateSpoke(figure, u + steps[i], v + steps[j], spokeId));
                            spokes.push_back(interpolateSpoke(figure, u - steps[i], v - steps[j], spokeId));
                            spokes.push_back(interpolateSpoke(figure, u + steps[i], v - steps[j], spokeId));
                            spokes.push_back(interpolateSpoke(figure, u - steps[i], v + steps[j], spokeId));
                            spokes.push_back(interpolateSpoke(figure, u + steps[i], v, spokeId));
                            spokes.push_back(interpolateSpoke(figure, u - steps[i], v, spokeId));
                            spokes.push_back(interpolateSpoke(figure, u, v + steps[i], spokeId));
                            spokes.push_back(interpolateSpoke(figure, u, v - steps[i], spokeId));
                        }
                    }
                }
            }
        }

        if(dynamic_cast<M3DQuadEndPrimitive*>(quadAtom) != 0){
            vtkSmartPointer<vtkSRep> srepfig = vtkSmartPointer<vtkSRep>::New();
            vtkIdType vtkAtomId = GetVtkSrepFig(targetObject, figureId, srepfig);
            GetCrestSpokes(srepfig, level, spokes,false, vtkAtomId, spokeId);
        }

        return spokes;

    }else{

        vector<M3DSpoke*> spokes;
        vtkSmartPointer<vtkSRep> srepfig = vtkSmartPointer<vtkSRep>::New();
        vtkIdType vtkAtomId = GetVtkSrepFig(targetObject, figureId, srepfig);
        GetCrestSpokes(srepfig, level, spokes,true, vtkAtomId, spokeId);

        return spokes;

    }
}


void Match::GetCrestSpokes(vtkSmartPointer<vtkSRep> srepfig, int level, vector<M3DSpoke*>& spokes,bool istube, vtkIdType vtkAtomId, int spokeId){

    if(false){//Set this to true to disable the crest interpolation, without crest interpolation, it should return only the base crest spokes
        level = 0;
    }

    vtkSmartPointer<vtkSRepInterpolateCrestSpokesQuartic> interpolatecrestspokes = vtkSmartPointer<vtkSRepInterpolateCrestSpokesQuartic>::New();
    interpolatecrestspokes->SetInput(srepfig);
    interpolatecrestspokes->SetInterpolationLevel(level);
    interpolatecrestspokes->SetAtomId(vtkAtomId);
    if(vtkAtomId != -1){
        interpolatecrestspokes->SetGamma_t(0.5);
        interpolatecrestspokes->SetGamma_theta(0.5);
    }else{
        interpolatecrestspokes->SetGamma_t(1);
        interpolatecrestspokes->SetGamma_theta(1);
    }
    if(istube){        
        interpolatecrestspokes->SetCyclicCurve(false);
        interpolatecrestspokes->SetCyclicSpokes(true);
        interpolatecrestspokes->SetSpokeType(spokeId);
    }else{
        interpolatecrestspokes->SetCyclicCurve(true);
        interpolatecrestspokes->SetCyclicSpokes(false);
        if(spokeId == 0){
            interpolatecrestspokes->SetSpokeType(vtkSRep::TOP_SPOKE);
        }else if(spokeId == 1){
            interpolatecrestspokes->SetSpokeType(vtkSRep::BOTTOM_SPOKE);
        }else if(spokeId == 2){
            interpolatecrestspokes->SetSpokeType(vtkSRep::CREST_SPOKE);
        }
    }

    interpolatecrestspokes->Update();


    vtkSRep* srepcrest = interpolatecrestspokes->GetSRepOutput();


    for(unsigned i = 0; i < srepcrest->GetNumberOfPoints(); i++){

        double point[3];

        srepcrest->GetPoint(i, point);

        vtkSRep::VectorVNLType currentspokes = srepcrest->GetSpokes(i);        

        for(unsigned j = 0; j < currentspokes.size(); j++){
            vtkSRep::VNLType s = currentspokes[j];
            double radius = s.magnitude();
            s.normalize();
            M3DSpoke* spoke = new M3DSpoke(point[0], point[1], point[2], s[0], s[1], s[2], radius);

            spokes.push_back(spoke);
        }
    }

}
void Match::GetInternalSpokes(vtkSmartPointer<vtkSRep> srepfig, int level, vector<M3DSpoke*>& spokes, vtkIdType atomId, int spokeId){
    vtkSmartPointer< vtkSRepInterpolateMedialSpokesHermite > medialspokesinterpolator = vtkSmartPointer< vtkSRepInterpolateMedialSpokesHermite >::New();

    medialspokesinterpolator->SetInput(srepfig);
    medialspokesinterpolator->SetInterpolationLevel(level);
    if(atomId != -1){
        medialspokesinterpolator->SetAtomId(atomId);
        medialspokesinterpolator->SetGamma_u(0.5);
        medialspokesinterpolator->SetGamma_v(0.5);
    }
    if(spokeId == 0){
        medialspokesinterpolator->SetSpokeType(vtkSRep::TOP_SPOKE);
    }else if(spokeId == 1){
        medialspokesinterpolator->SetSpokeType(vtkSRep::BOTTOM_SPOKE);
    }

    medialspokesinterpolator->Update();

    vtkSRep* srepcrest = medialspokesinterpolator->GetSRepOutput();


    for(unsigned i = 0; i < srepcrest->GetNumberOfPoints(); i++){

        double point[3];

        srepcrest->GetPoint(i, point);

        vtkSRep::VectorVNLType currentspokes = srepcrest->GetSpokes(i);
        vtkSRep::VectorDoubleType radius = srepcrest->GetSpokesRadius(i);

        for(unsigned j = 0; j < currentspokes.size(); j++){
            M3DSpoke* spoke = new M3DSpoke(point[0], point[1], point[2], currentspokes[j][0], currentspokes[j][1], currentspokes[j][2], radius[j]);

            spokes.push_back(spoke);
        }
    }
}

vtkIdType Match::GetVtkSrepFig(M3DObject* targetObject, int figureId, vtkSmartPointer<vtkSRep>& srepfig, int atomId){

    M3DFigure * figure = targetObject->getFigurePtr( figureId ) ;
    M3DQuadFigure* quadfig = dynamic_cast<M3DQuadFigure*>( figure );

    vtkIdType vtkAtomId = -1;

    if(quadfig){

        vtkSmartPointer< vtkPoints > hubpos = vtkSmartPointer< vtkPoints >::New();
        vtkSmartPointer<vtkCellArray> cellarray = vtkSmartPointer<vtkCellArray>::New();

        vtkSRep::VectorSRepIdsType pointsIds;

        vtkSRep::RadiusVectorType allradius;
        vtkSRep::SpokesVectorType allspokes;

        for(int u = 0; u < quadfig->getRowCount(); u++){
            pointsIds.push_back(vtkSRep::VectorIdsType());
            for(int v = 0; v < quadfig->getColumnCount(); v++){

                M3DQuadPrimitive* prim0 = dynamic_cast<M3DQuadPrimitive*>(quadfig->getPrimitivePtr(u, v));

                Vector3D x = prim0->getX();
                Vector3D u0 = prim0->getU0();
                Vector3D u1 = prim0->getU1();

                vtkSRep::VectorVNLType vnlspokes;
                vtkSRep::VNLType s(3);
                s[0] = u0.getX();
                s[1] = u0.getY();
                s[2] = u0.getZ();
                vnlspokes.push_back(s);

                s[0] = u1.getX();
                s[1] = u1.getY();
                s[2] = u1.getZ();
                vnlspokes.push_back(s);

                vtkSRep::VectorDoubleType radius;
                radius.push_back(prim0->getR0());
                radius.push_back(prim0->getR1());

                if(u == 0 || u == quadfig->getRowCount() - 1 || v == 0 || v == quadfig->getColumnCount() - 1){

                    M3DQuadEndPrimitive* prim0 = dynamic_cast<M3DQuadEndPrimitive*>(quadfig->getPrimitivePtr(u, v));
                    Vector3D uend = prim0->getUEnd();

                    s[0] = uend.getX();
                    s[1] = uend.getY();
                    s[2] = uend.getZ();

                    vnlspokes.push_back(s);

                    radius.push_back(prim0->getREnd());

                }

                int numCols = quadfig->getColumnCount();
                if(atomId != -1 && atomId / numCols == u && atomId % numCols == v){
                    vtkAtomId = hubpos->InsertNextPoint(x.getX(), x.getY(), x.getZ());
                    pointsIds[u].push_back(vtkAtomId);
                }else{
                    pointsIds[u].push_back(hubpos->InsertNextPoint(x.getX(), x.getY(), x.getZ()));
                }


                allspokes.push_back(vnlspokes);
                allradius.push_back(radius);
            }
        }


        for(unsigned i = 0; i < pointsIds.size() - 1; i++){
             for(unsigned j = 0; j < pointsIds[i].size() - 1; j++){

                 vtkSmartPointer<vtkQuad> quad = vtkSmartPointer<vtkQuad>::New();
                 quad->GetPointIds()->SetId(0, pointsIds[i][j]);
                 quad->GetPointIds()->SetId(1, pointsIds[i+1][j]);
                 quad->GetPointIds()->SetId(2, pointsIds[i+1][j+1]);
                 quad->GetPointIds()->SetId(3, pointsIds[i][j+1]);

                 //quad->Print(cout);

                 cellarray->InsertNextCell(quad);

             }
         }

        srepfig->SetPoints(hubpos);
        srepfig->SetPolys(cellarray);
        srepfig->SetAllSpokes(allspokes);
        srepfig->SetAllRadius(allradius);
        //srepfig->SetGridTopolgyIds(pointsIds);

        for(unsigned i = 0; i < pointsIds.size(); i++){
            pointsIds[i].clear();
        }
        pointsIds.clear();
        for(unsigned i = 0; i < allradius.size(); i++){
            allradius[i].clear();
            allspokes[i].clear();
        }
        allradius.clear();
        allspokes.clear();

    }else{

        M3DTubeFigure* tubefig = dynamic_cast<M3DTubeFigure*>(figure);

            vtkSmartPointer< vtkPoints > hubpos = vtkSmartPointer< vtkPoints >::New();
            vtkSmartPointer<vtkCellArray> cellarray = vtkSmartPointer<vtkCellArray>::New();

            vtkSRep::VectorSRepIdsType pointsIds;

            vtkSRep::RadiusVectorType allradius;
            vtkSRep::SpokesVectorType allspokes;



            pointsIds.push_back(vtkSRep::VectorIdsType());

            for(int u = 0; u < (int)tubefig->getPrimitiveCount(); u++){

                M3DTubePrimitive* prim0 = dynamic_cast<M3DTubePrimitive*>(tubefig->getPrimitivePtr(u));
                Vector3D x = prim0->getX();

                vtkIdType id = hubpos->InsertNextPoint(x.getX(), x.getY(), x.getZ());

                if(atomId == u){
                    vtkAtomId = id;
                }

                pointsIds[0].push_back(id);

                vtkSRep::VectorVNLType vnlspokes;
                vtkSRep::VectorDoubleType radius;

                for(int v = 0; v < (int)tubefig->getNumberOfSpokes(); v++){


                    Vector3D u0 = prim0->getYN(v);


                    vtkSRep::VNLType s(3);
                    s[0] = u0.getX();
                    s[1] = u0.getY();
                    s[2] = u0.getZ();

                    radius.push_back(s.magnitude());
                    vnlspokes.push_back(s.normalize());


                }

                allspokes.push_back(vnlspokes);
                allradius.push_back(radius);
            }


            for(unsigned i = 0; i < pointsIds.size(); i++){
                 for(unsigned j = 0; j < pointsIds[i].size() - 1; j++){

                     vtkSmartPointer<vtkLine> line = vtkSmartPointer<vtkLine>::New();
                     line->GetPointIds()->SetId(0, pointsIds[i][j]);
                     line->GetPointIds()->SetId(1, pointsIds[i][j+1]);

                     //quad->Print(cout);

                     cellarray->InsertNextCell(line);

                 }
             }

            srepfig->SetPoints(hubpos);
            srepfig->SetLines(cellarray);
            srepfig->SetAllSpokes(allspokes);
            srepfig->SetAllRadius(allradius);

    }

    return vtkAtomId;
}


#endif

