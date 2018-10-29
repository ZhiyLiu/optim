// dibyendu - for memory leak detection

//#define _CRTDBG_MAP_ALLOC
//#include <stdlib.h>
//#include <crtdbg.h>
// -----------------------------------
// dibyendu
#include <iomanip>

#include <stdio.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#ifdef _WIN32
#include <direct.h>r
#else
#include <limits.h>
#include <unistd.h>
#include <errno.h>
#endif
#include <iostream>
#include <iterator>
#include <string.h>

#include <sstream>

#include <typeinfo>
#include <memory>

#include "globalBuildParms.h"
#ifdef BINARY
#include <fstream>
#endif
#include "ControlParms.h"
#include "SimilarityTransform3D.h"
#include "M3DObject.h"
#include "M3DObjectFile.h"
#include "WorldSystem.h"
#include "ImageResample3D.h"
#include "GaussianBlur3D.h"
#include "DQFImage.h"
#include "RAWImageFile.h"
#include "AllImageIO.h"
#include "LandmarkDeformation.h"
#include "SurfacePatchEnsemble.h"
#include "DQFMatch.h"
#include "ObjectRelativeSampling.h"
#include "M3DQuadInterpolator.h"
#include "M3DSubfigureTransformation.h"
#include "M3DFigureTreeNode.h"

#include "M3DSpokeLengthOptimizer.h"
#include "M3DSpokeAngleOptimizer.h"
#include "M3DNEWUOAOptimizer.h"

#include "DistanceMap3D.h"
#include "utility.h"
#include "BYU.h"
#ifdef BINARY
#include "DistanceVectorList.h"
#endif
#include "M3DAtomPredictor.h"
#include "Geodesic.h"
#include "OptimizerBase.h"
#include "M3DAdaptiveRegistrationPGAOptimizer.h"
#ifndef BINARY
#include "M3DMainFigResiduePGAOptimizer.h"
#include "M3DPGAStats.h"
#else
#include "M3DPGAPrimitiveStats.h"
#include "utility.h"
#endif
#ifdef PCA_VISUALIZER
#include "PCA.h"
#endif
#include "P3DControl.h"
#ifdef OPTIMIZATION_VISUALIZER
#include "LogManager.h"
#endif
#include "FitUnlabeledPoints.h"
#include "Tuning.h"
//#include "snapshot.h"

// dibyendu - to make optimization video
#include "MakeOptimizationVideo.h"

//#define DEBUG

//#define FAST_ALIGNMENT  /* Strictly for debugging */


#ifndef S_ISREG
#define S_ISREG(m)        (((S_IFREG & m) != 0) ? true : false)
#endif

using namespace std;


M3DObject * readOldObjectFile(const char * filename);

const imageLandmark_t emptyImageLandmark = { 0.0, 0.0, 0.0, -1, -1, NULL };

const char mappingMsg[] =
"Transforming the model to the world coordinate system of the loaded image\n";

// The following global variable is declared in Tuning.h
#ifdef BINARY
class bpTuning tuneVals;		// Instantiation of tuning parameters for Binary Pablo
#else
class gpTuning tuneVals;		// Instantiation of tuning parameters for regular Pablo
#endif

bool P3DControl::verifyApplication(application_t id) const
{
	bool appDefined = false;

#ifndef BINARY
#ifdef AE2_BUILD
	appDefined = true;
	if (id == AE2)
		return true;
#endif

	/*  Add other non-BINARY applications here, in this manner:
	#ifdef XYZ_BUILD
	appDefined = true;
	if (id == xyz)
	return true;
	#endif
	*/
#endif	/* BINARY */

#ifdef BINARY
#ifdef AE2_BUILD
	appDefined = true;
	if (id == AE2Binary)
		return true;
#endif
	//  Add other BINARY applications here, as above
#endif	/* BINARY */

	if (! appDefined) {
		if (id == Pablo)
			return true;
	}
	return false;
}

bool P3DControl::verifyOrientation() const
{
#ifdef UNFLIPPED
	return true;
#else
	return false;
#endif
}

P3DControl::P3DControl(int undoLength) : undoList(undoLength)
{
#ifdef DEBUG
	cout << "P3DControl::P3DControl()  this = 0x" << hex << this << dec << '\n';
#endif
	object = NULL;
#ifdef AE2_BUILD
	anaList = NULL;
#endif
	objectSimTransform = NULL;
	recordedTransformation = NULL;
#ifndef BINARY
	handPlacementTransform = new SimilarityTransform3D;
#endif

	imagePlanes.setSmoothing(globalControl->readBool(SmoothImages));
	AllImageIO::setImageScaling(globalControl->readBool(ScaleImages));
	imagePlanesChanged = true;
	// Indicates user changed atom/figure selection in an optimization stage
	selectionChanged = false;
	worldWasApplied = false;

	tileSet = NULL;
	tileSetMobile = false;
	for (int f=0; f < 6; f++)
		tileSetMobileFigure[f] = true;
	auxTileSet = NULL;

	objectCopyBuffer = NULL;

	setEditingMode(ANY_EDITING_MODE);
	allowSelection(true);

	inOperation = false;
	modelIsSaved();

	modelFileName = NULL;
	pgaModelFileName = NULL;

#ifdef BINARY
	//	lmrkFileName = NULL;
	dist_map = NULL;	// See setDistMapPtr()
#endif
	imageFileName = NULL;

	markedPrimitiveId = INVALID_PRIMITIVE_ID;

	interpolatedPrimitive = NULL;

	match = NULL;

	drawing_partial_surface = false;
	drawing_full_surface = false;

#ifndef BINARY
	optimizeWithWindowing = true;	// Same as in P3DUserInterface.fl
#else
	optimizeWithWindowing = false;
#endif
	registrationOptimizer	= NULL;
	mainFigureOptimizer		= NULL;
	subfigureOptimizer		= NULL;
	deformationOptimizer	= NULL;
	//dibyendu
	sRepOptimizer			= NULL;
	spokeOptimizer			= NULL;
	voxelStageOptimizer		= NULL;

#ifndef BINARY
	penaltyWeights = NULL;
	penaltyNames = NULL;
	numPenaltyWeights = 0;
#endif

	referenceObject = NULL;
	recordingTransformation = false;

	elongater = NULL;

	currentOptimizerStage = NO_OPTIMIZER_STAGE;

	numImageLandmarks = 0;

	treeFigures = NULL;
	slideShowSavedObject = NULL;
	slideShowSavedEditingMode = ANY_EDITING_MODE;
	inSlideShow = false;

#ifdef LM_METHOD_OBJ
	landmarkObject = NULL;
#endif

#ifndef BINARY
	spe = NULL;
	dqf = NULL;
	mp  = NULL;
#endif

	matchReferenceObject = NULL;

	registrationPenaltyWeight = 0.5;
	pgaStats = NULL;	// Owned by *object
	atomPgaStats = NULL;	// Owned by *object
	attachedStats = false;
	residuePGAStage = 0;

	// dibyendu - cpns
	cpnsStats = NULL ;	// Owned by *object

#ifdef PCA_VISUALIZER
	pca = NULL;
#endif

	updateGui = NULL;
	updateGuiStage = NULL;

	// default to some useful clip planes
	clips[0].mask = 0;
	clips[0].equation[0] =  1.0;
	clips[0].equation[1] =  0.0;
	clips[0].equation[2] =  0.0;
	clips[0].equation[3] =  0.5;

	clips[1].mask = 0;
	clips[1].equation[0] = -1.0;
	clips[1].equation[1] =  0.0;
	clips[1].equation[2] =  0.0;
	clips[1].equation[3] =  0.5;

	clips[2].mask = 0;
	clips[2].equation[0] =  0.0;
	clips[2].equation[1] =  1.0;
	clips[2].equation[2] =  0.0;
	clips[2].equation[3] =  0.5;

	clips[3].mask = 0;
	clips[3].equation[0] =  0.0;
	clips[3].equation[1] = -1.0;
	clips[3].equation[2] =  0.0;
	clips[3].equation[3] =  0.5;

	clips[4].mask = 0;
	clips[4].equation[0] =  0.0;
	clips[4].equation[1] =  0.0;
	clips[4].equation[2] =  1.0;
	clips[4].equation[3] =  0.5;

	clips[5].mask = 0;
	clips[5].equation[0] =  0.0;
	clips[5].equation[1] =  0.0;
	clips[5].equation[2] = -1.0;
	clips[5].equation[3] =  0.5;
}

P3DControl::~P3DControl()
{	

#ifdef DEBUG
	cout << "P3DControl::~P3DControl()  this = 0x" << hex << this << dec << '\n';
#endif
	if (inSlideShow)
		quitModelSlideShow();	// side-effect: quits tile slide show, too

	delete object;
	delete objectSimTransform;
	delete recordedTransformation;
#ifdef AE2_BUILD
	if (anaList != NULL)
		delete anaList;
#endif

#ifndef BINARY
	delete handPlacementTransform;
#endif

	delete tileSet;

#ifndef BINARY
	if (auxTileSet != NULL) {
		BYU * auxBYU;	// This mirrors auxTileSet
		auxBYU = (BYU *) auxTileSet;
		delete auxBYU;
	}
#else
	if (dist_map)
	{
		delete dist_map;
		dist_map = NULL;
	}
#endif

	if (objectCopyBuffer != NULL)
		delete objectCopyBuffer;

	if (modelFileName != NULL)
		delete [] modelFileName;

	if (pgaModelFileName != NULL)
		delete [] pgaModelFileName;

	//if (imageFileName != NULL)
	//	delete [] imageFileName;		

	if (interpolatedPrimitive != NULL)
		delete interpolatedPrimitive;

	if (match != NULL)
		delete match;

	if (registrationOptimizer != NULL)
		delete registrationOptimizer;

	if (mainFigureOptimizer != NULL)
		delete mainFigureOptimizer;

	if (subfigureOptimizer != NULL)
		delete subfigureOptimizer;

	if (deformationOptimizer != NULL)
		delete deformationOptimizer;

	if (sRepOptimizer != NULL)
		delete sRepOptimizer;

	if (spokeOptimizer != NULL)
		delete spokeOptimizer;

	if (voxelStageOptimizer != NULL)
		delete voxelStageOptimizer;

#ifndef BINARY
	if (penaltyWeights != NULL) {
		delete [] penaltyWeights;
		delete [] penaltyNames;
	}

	if (spe != NULL)
		delete spe;
#else
	if (dist_map)
	{
		delete dist_map;
		dist_map = NULL;
	}
#endif

	if (referenceObject != NULL)
		delete referenceObject;

	delete elongater;

	if (treeFigures != NULL)
		delete [] treeFigures;

#ifdef LM_METHOD_OBJ
	delete landmarkObject;
#endif

	if (slideShowSavedObject != NULL)
		delete slideShowSavedObject;

	for (int i = 0; i < modelSlideShow.size(); i++)
		delete modelSlideShow[i];
	modelSlideShow.clear();

	for (int i = 0; i < tileSlideShow.size(); i++)
		delete tileSlideShow[i];
	tileSlideShow.clear();

#ifdef PCA_VISUALIZER
	if (pca != NULL)
		delete pca;
#endif

	// dibyendu - to detect memory leaks

	//_CrtDumpMemoryLeaks();

	// -----------------------------------
}

void P3DControl::setGuiCallbacks(void (* guiCB)(int), void (* guiStageCB)(OptimizerStage))
{
	updateGui = guiCB;
	updateGuiStage = guiStageCB;
}

#ifdef LM_METHOD_OBJ

/*	Returns a position for the named landmark, either from a normal model
(in which case the landmark names point to atoms) or in an image landmark
model (in which case the figure's name is the landmark's name and the
position of it's first atom is the landmark's position).
The objType is either 1 for a normal model or 2 for a landmark model.
FigureId is only used for normal models.
*/
int P3DControl::findLandmarkPositByName(Vector3D &lmPosit, M3DObject * obj, int figureId,
										const char * landmarkName, int objType)
{
	if (obj == NULL) {
		cout << "WARNING: internal error in findLandmarkPositByName: invalid obj" << endl;
		return 0;		// nothing to match
	}

	M3DPrimitive *lmAtom;
	M3DFigure *fig;
	int lmAtomIndex, li, nli, lmAtomT, lmIndex;

	switch (objType)
	{
	case 1:
		// "lm"-style values are the position of the model landmark,
		//	which is moved during optimization.
		if (figureId < 0) {
			cout << "WARNING: internal error in findLandmarkPositByName: figureId < 0"
				<< endl;
			return 0;		// nothing to match
		}

		// find the atom whose spoke is the landmark
		fig = obj->getFigurePtr(figureId);
		lmAtomIndex = fig->findLandmarkAtomIndexByName(landmarkName);
		if (lmAtomIndex < 0)
		{
			cout << "WARNING: cannot find model landmark \""
				<< landmarkName << "\"" << endl;
			return 0;
		}
		lmAtom = fig->getPrimitivePtr(lmAtomIndex);

		// choose which spoke to find tip of
		lmIndex = fig->findLandmarkIndexByName(landmarkName);
		lmAtomT = (int) fig->landmarkAtomTs[lmIndex];
		switch (lmAtomT)
		{
		case -1:
			lmPosit = lmAtom->getX() + lmAtom->getY0();
			break;
		case 0:
			lmPosit = lmAtom->getX() + lmAtom->getExtendedB();
			break;
		case 1:
			lmPosit = lmAtom->getX() + lmAtom->getY1();
			break;
		}
		break;
	case 2:			
		// "li"-style values are the position of the landmark in the image,
		//	never moved after loading
		nli = obj->getFigureCount();
		for (li = 0; li < nli; li++)
		{
			fig = obj->getFigurePtr(li);
			if (! strcmp(landmarkName, fig->getName()+3))
			{
				lmPosit = fig->getPrimitivePtr(0)->getX();
				return 1;
			}
		}
		if (globalVerbosity >= 0)
			cout << "WARNING: cannot find image landmark \""
			<< landmarkName << "\"" << endl;
		return 0;
		break;
	default:
		cout << "WARNING: internal error in findLandmarkPositByName: "
			<< "invalid objType (" << objType << ")" << endl;
		return 0;
		break;
	}
	return 1;
}

#endif

#ifdef BINARY

// Resample the binary image if it's not available as a file;
// write the resampled file for next time
Image3D * P3DControl::getResampedBinaryScript(Registry & script,
											  char * binImageFileKey, char * resBinImageFileKey)
{
	const char * filename;
	Image3D * image;
	RAWImageFile imageFile;

	filename = script.getStringValue(resBinImageFileKey, NULL);
	image = imageFile.read(filename);

	if (! image) {
		filename = script.getStringValue(binImageFileKey, NULL);
		image = imageFile.read(filename);

		cout << "Loaded Binary Image: " << filename << endl;

		//resample the first binary
		ImageResample3D resampler;
		resampler.isotropicSample(*image);

		cout << "Resampled Binary Image.\n";

		filename = script.getStringValue(resBinImageFileKey, NULL);
		if(imageFile.write(filename, *image))
			cout << "Saved resampled Image: " << filename << endl;
		else
			cout << "Resampled Image failed to save: " << filename << endl;
	}
	else
		cout << "Loaded Resampled Binary Image: " << filename << endl;
	return image;
}

ImageDistanceMap * P3DControl::getDDMScript(Registry & script,
											Image3D * image, char * filenameIn)
{
	const char * filename;
	RAWImageFile imageFile;
	ImageDistanceMap * binaryDistanceMap;
	Image3D * distMapImage;
	Image3D * imageOut;

	filename = script.getStringValue(filenameIn, NULL);
	distMapImage = imageFile.read(filename);

	if (distMapImage != NULL) {
		cout << "Loading DistanceMap from file.\n";
		binaryDistanceMap = new ImageDistanceMap;
		binaryDistanceMap->fromImage3D(distMapImage);
	}
	else {        
		cout << "Computing DistanceMap.\n";
		binaryDistanceMap = new ImageDistanceMap(image);
		(void) binaryDistanceMap->createMap();
		filename = script.getStringValue(filenameIn, NULL);
		imageOut = binaryDistanceMap->operator Image3D *();
		saveImage(filename, imageOut);
		delete imageOut;
	}

	return binaryDistanceMap;
}

#endif

#ifdef LM_METHOD_OBJ

void P3DControl::binaryPabloListLandmarks()
{
	if (object == NULL)
		return;

	if (globalVerbosity > 0) {
		cout << " Found " << object->getLandmarkCount() << " landmarks" << endl;

		// Print landmarks in each figure
		int nf = object->getFigureCount();
		for (int f = 0; f < nf; f++) {
			M3DFigure *fig = object->getFigurePtr(f);

			// Print the position of each landmark in figure
			int nlm = fig->getLandmarkCount();
			for (int lm = 0; lm < nlm; lm++) {					
				Vector3D lmPosit;
				char * lmName = (char*) fig->getLandmarkName(lm);
				findLandmarkPositByName(lmPosit, object, f, lmName, 1);

				cout << "  lm #" << lm << ": name=\"" << lmName << "\"" << endl;
				cout << "        position=<ucs>";
				lmPosit.print();
				Image3D * image = getImagePtr();
				if (image) {
					image->modelToImageCoordinates(lmPosit);
					cout << "        position=<ics>";
					lmPosit.print();
				}
			}
		}
	}
}

#endif

// ----------------------------- Start of Scripting Code --------------------------------

#ifdef LM_METHOD_OBJ

bool P3DControl::binaryPabloLoadLandmarkModel(const char * filename)
{
	if (filename && *filename) {
		// Load the first model
		if (! readLandmarkModel(filename)) {
			if (filename && *filename)
				cout << "Error: unable to load landmark model: " << filename << endl;
			else
				cout << "Error: no landmark model filename supplied" << endl;
			return false;
		}

		if (globalVerbosity > -1)
			cout << "Loaded landmark model " << filename << endl;
		if (globalVerbosity > 0) {
			cout << " Found " << landmarkObject->getFigureCount() << " landmarks" << endl;

			// print the position of each landmark
			int nli = landmarkObject->getFigureCount();
			for (int li = 0; li < nli; li++) {
				const char * lmName = landmarkObject->getFigurePtr(li)->getName();

				Vector3D lmPosit;
				findLandmarkPositByName(lmPosit, landmarkObject, 0, (char*) lmName+3, 2);

				cout << "  lm #" << li << ": name=\"" << lmName << "\" position=<ucs>";
				lmPosit.print();
			}
		}

		return true;
	}
	return false;	// null filename
}

#endif

#ifdef BINARY

// This only gets called by Eli's scripts - this isn't the generic GUI-based figure stage
bool P3DControl::binaryPabloFigureStage(int maxNumIterations,
										int doScale, // bitFlags 1=scale
										const char * pgaFilename)
{
	cout << "binaryPabloFigureStage" << endl;
	if (! initializeMainFigureOptimizer(doScale, pgaFilename))
		return false;
	cout << "Optimizer Initialized" << endl;

	// Dibyendu

	// set the reference object after scale M3DObject inside the Match class

	//if ( match != NULL ) {
	//	M3DObject* scaledObject = object->assign() ;
	//	match->setRefObjectAfterScale( scaledObject ) ;
	//}

	if (! maxNumIterations)		// nothing to do: used to do 65535 iterations
		return true;

	if (updateGuiStage != NULL)
		updateGuiStage(currentOptimizerStage);

	for (int i = 0; i < maxNumIterations; i++) {
		doMainFigureIterations(1);
		if (mainFigureStageSkipped())
			return true;

		if (updateGui != NULL)
			updateGui(i);

		// BINARY PABLO
		// Set output values, added appropriate accessors to control
		if (globalVerbosity > -1) {
			const Match::matchResult * mr = match->getFigureStageResults();
			cout << " " << i << ": ";
			cout << "ObjFunc:" << getLastMainFigureObjectiveFunctionValue();
			for (int j = 0; j < MAX_NUM_FIG_MATCH_RESULTS; j++) {
				bpTune_t index = (bpTune_t) mr[j].tuningParm;
				switch (index) {	// Why not just print the "short name"?
					case BpFigureImageMatch:			cout << " Match";
						break;
					case BpFigureLandmarkMatch:			cout << " Landmark";
						break;
					case BpFigureModelMatch:			cout << " Ref";
						break;
					case FigureMahalanobisMatch:		cout << " Maha";
						break;
					case BpFigureInterpenetrationMatch:	cout << " Pen";
						break;
					case BpFigureRSradPenalty:			cout << " RSrad";
						break;
					case BpFigureCurvinessPenalty:		cout << " Curviness";
						break;
					case FigureContourMatch:			cout << " Contour";
						break;
					case BpFigureImageNormalMatch:		cout << " Boundary Normal";
						break ;
					case BpFigureImagePlaneOrienAtEnd:	cout << " Crest principal direction";
						break ;
					case BpFigureImageVertexMatch:		cout << " Crest vertex" ;
						break ;

					default:		cout << " Unknown match type";
						break;
				}
				cout << ':' << mr[j].value;
			}
			cout << '\n';
			cout.flush();
		}		

		if( MakeOptimizationVideo::makeVideoOn )			
			MakeOptimizationVideo::writeCurrentModel( object ) ;	// dibyendu - save current model		


	} // end of figure stage iterations

	if (globalVerbosity > -1) {
		cout << "Match after Figure Stage:" << endl;
		int figureId = (int) tuningWt(BpFigureId);
		MatchUtility::calculateAverageDistance(object, object->loadedObject(),
			figureId, getImageDistMapPointer(), (int) tuningWt(BpFigureSurfaceLevel),
			(enum PrimNeighborhoodDefn) (int) tuningWt(BpFigureNeighborhood),
			(enum DistanceType) (int) tuningWt(BpFigureDistanceType),
			globalVerbosity);
	}

	return true;
}

bool P3DControl::binaryPabloAtomStage(int maxNumIterations)
{
	int i;

	bool optimizerStopRequested = false;
	testAndResetModelChanged();
	initializeDeformationOptimizer();

	if (updateGuiStage != NULL)
		updateGuiStage(currentOptimizerStage);

	// Optimize the figure
	for (i = 0; i < maxNumIterations; i++) {
		// Essentially a while(1) loop
		int figureId = (int) tuningWt(BpFigureId); //MultiObject

		// Optimize selected atoms:
		// this next call IS the atom loop - it returns true when enough atoms
		// are at their best, which includes skipping atoms which
		// seem to stop improving.
		optimizerStopRequested = doDeformationIterations(figureId);
		if (updateGui != NULL)
			updateGui(i);

		// Set output values
		if (globalVerbosity > -1) {
			const Match::matchResult * mr = match->getAtomStageResults();
			cout << " " << i << ": ";
			cout << "ObjFunc:" << getLastDeformationObjectiveFunctionValue(false);
			for (int j = 0; j < MAX_NUM_ATOM_MATCH_RESULTS; j++) {
				bpTune_t index = (bpTune_t) mr[j].tuningParm;
				switch (index) {	// Why not just print the "short name"?
					case BpAtomImageMatch:				cout << " Match";
						break;
					case BpAtomLandmarkMatch:			cout << " Landmark";
						break;
					case BpAtomAverageNeighborMatch:	cout << " Reg";
						break;
					case BpAtomModelMatch:				cout << " Ref";
						break;
					case AtomMahalanobisMatch:			cout << " Maha";
						break;
					case BpAtomInterpenetrationMatch:	cout << " Pen";
						break;
					case BpAtomRSradPenalty:			cout << " Srad";
						break;
					case BpAtomDeltaNeighborMatch:		cout << " DelN";
						break;
					case BpAtomCurvinessPenalty:		cout << " Curviness";
						break;
					case BpAtomImageNormalMatch:		cout << " Boundary Normal";
						break ;
					case BpAtomImagePlaneOrienAtEnd:	cout << " Crest principal direction";
						break ;
					case BpAtomImageVertexMatch:		cout << " Crest vertex" ;
						break ;
				}
				cout << ':' << mr[j].value;
			}
			cout << '\n';
			cout.flush();
		}

		if( MakeOptimizationVideo::makeVideoOn )			
			MakeOptimizationVideo::writeCurrentModel( object ) ;	// dibyendu - save current model		

		if (optimizerStopRequested == true)
			break;

	} // end of iterations

	// IS THIS DOING ANYTHING? THE PURPOSE WAS TO TELL THE USER THAT MORE
	// OPTIMIZATION IS REQUIRED THAN WAS ALLOWED, EG, THE OBJFUNC DID NOT
	// LEVEL OUT YET
	if (i == (maxNumIterations - 1))
		if (globalVerbosity > -1)
			cout << "\nOptimizer could not reach an optimum in "
			<< maxNumIterations << " iterations" << endl;

	if (globalVerbosity > 0) {
		const Match::matchResult * mr = match->getAtomStageResults();
		cout << "\nMatch after Atom Stage:\n";
		// Note: BP omits printing of the penetration penalty
		for (int j = 0; j < MAX_NUM_ATOM_MATCH_RESULTS - 1; j++) {
			bpTune_t index = (bpTune_t) mr[j].tuningParm;
			switch (index) {	// Why not just print the "short name"?
				case BpAtomImageMatch:				cout << " Match Distance";
					break;
				case BpAtomLandmarkMatch:			cout << " Landmark Distance";
					break;
				case BpAtomAverageNeighborMatch:	cout << " Match Regularity";
					break;
				case BpAtomModelMatch:				cout << " Match Reference";
					break;
				case AtomMahalanobisMatch:			cout << " Mahalanobis Distance";
					break;					
				case BpAtomInterpenetrationMatch:	cout << " Pen";
					break;
				case BpAtomRSradPenalty:			cout << " Srad";
					break;
				case BpAtomDeltaNeighborMatch:		cout << " Delta Neighbor";
					break;
			}
			cout << ": " << mr[j].value << '\n';
		}
		cout << flush;
	}

	int figureId = (int) tuningWt(BpFigureId);
	(void) MatchUtility::calculateAverageDistance(object, object->loadedObject(),
		figureId, getImageDistMapPointer(), (int) tuningWt(BpAtomSurfaceLevel),
		(enum PrimNeighborhoodDefn) (int) tuningWt(BpAtomNeighborhood),
		(enum DistanceType) (int) tuningWt(BpAtomDistanceType),
		globalVerbosity);

	return true;
}



bool P3DControl::binaryPabloSRepStage(int maxNumIterations)
{
	int i;

	bool optimizerStopRequested = false;
	testAndResetModelChanged();
	initializeSRepOptimizer();

	if (updateGuiStage != NULL)
		updateGuiStage(currentOptimizerStage);

	// Optimize the figure
	for (i = 0; i < maxNumIterations; i++) {
		// Essentially a while(1) loop
		int figureId = (int) tuningWt(BpFigureId); //MultiObject

		// Optimize selected atoms:
		// this next call IS the atom loop - it returns true when enough atoms
		// are at their best, which includes skipping atoms which
		// seem to stop improving.
		optimizerStopRequested = doSRepIterations(figureId);
		if (updateGui != NULL)
			updateGui(i);

		// Set output values
		if (globalVerbosity > -1) {
			const Match::matchResult * mr = match->getSRepStageResults();
			cout << " " << i << ": ";
			cout << "ObjFunc:" << getLastSRepObjectiveFunctionValue(false);
			for (int j = 0; j < MAX_NUM_SREP_MATCH_RESULTS; j++) {
				bpTune_t index = (bpTune_t) mr[j].tuningParm;
				switch (index) {	// Why not just print the "short name"?
					case BpSRepAverageNeighborMatch:	cout << " Reg";
						break;
					case BpSRepImageMatch:				cout << " Match";
						break;
					case SRepRadiusMatch:			cout << " Atom Radius";
						break;
					case BpSRepLandmarkMatch:			cout << " Landmark";
						break;					
					case BpSRepModelMatch:				cout << " Ref";
						break;
					case BpSRepInterpenetrationMatch:	cout << " Pen";
						break;
					case BpSRepRSradPenalty:			cout << " Srad";
						break;
					case BpSRepDeltaNeighborMatch:		cout << " DelN";
						break;
					case BpSRepCurvinessPenalty:		cout << " Curviness";
						break;
					case BpSRepImageNormalMatch:		cout << " Boundary Normal";
						break ;
					case BpSRepImagePlaneOrienAtEnd:	cout << " Crest principal direction";
						break ;
					case BpSRepImageVertexMatch:		cout << " Crest vertex" ;
						break ;
				}
				cout << ':' << mr[j].value;
			}
			cout << '\n';
			cout.flush();
		}

		if( MakeOptimizationVideo::makeVideoOn )			
			MakeOptimizationVideo::writeCurrentModel( object ) ;	// dibyendu - save current model		

		if (optimizerStopRequested == true)
			break;
	} // end of iterations

	// IS THIS DOING ANYTHING? THE PURPOSE WAS TO TELL THE USER THAT MORE
	// OPTIMIZATION IS REQUIRED THAN WAS ALLOWED, EG, THE OBJFUNC DID NOT
	// LEVEL OUT YET
	if (i == (maxNumIterations - 1))
		if (globalVerbosity > -1)
			cout << "\nOptimizer could not reach an optimum in "
			<< maxNumIterations << " iterations" << endl;

	if (globalVerbosity > 0) {
		const Match::matchResult * mr = match->getSRepStageResults();
		cout << "\nMatch after Atom Stage:\n";
		// Note: BP omits printing of the penetration penalty
		for (int j = 0; j < MAX_NUM_ATOM_MATCH_RESULTS - 1; j++) {
			bpTune_t index = (bpTune_t) mr[j].tuningParm;
			switch (index) {	// Why not just print the "short name"?
				case BpSRepImageMatch:				cout << " Match Distance";
					break;
				case BpSRepLandmarkMatch:			cout << " Landmark Distance";
					break;
				case BpSRepAverageNeighborMatch:	cout << " Match Regularity";
					break;
				case BpSRepModelMatch:				cout << " Match Reference";
					break;
				case SRepRadiusMatch:				cout << " Radius mismatch";
					break;					
				case BpSRepInterpenetrationMatch:	cout << " Pen";
					break;
				case BpSRepRSradPenalty:			cout << " Srad";
					break;
				case BpSRepDeltaNeighborMatch:		cout << " Delta Neighbor";
					break;
			}
			cout << ": " << mr[j].value << '\n';
		}
		cout << flush;
	}

	int figureId = (int) tuningWt(BpFigureId);
	(void) MatchUtility::calculateAverageDistance(object, object->loadedObject(),
		figureId, getImageDistMapPointer(), (int) tuningWt(BpSRepSurfaceLevel),
		(enum PrimNeighborhoodDefn) (int) tuningWt(BpSRepNeighborhood),
		(enum DistanceType) (int) tuningWt(BpSRepDistanceType),
		globalVerbosity);

	return true;
}


bool P3DControl::binaryPabloSpokeStage(int maxNumIterations)
{
	int i;

	bool optimizerStopRequested = false;
	testAndResetModelChanged();
	initializeSpokeOptimizer();

	if (updateGuiStage != NULL)
		updateGuiStage(currentOptimizerStage);

	// Optimize the figure
	for (i = 0; i < maxNumIterations; i++) {
		// Essentially a while(1) loop
		int figureId = (int) tuningWt(BpFigureId); //MultiObject

		// Optimize spokes of selected atoms:
		// this next call IS the atom loop - it returns true when enough atoms
		// are at their best, which includes skipping atoms which
		// seem to stop improving.
		optimizerStopRequested = doSpokeIterations(figureId);
		if (updateGui != NULL)
			updateGui(i);

		// Set output values
		if (globalVerbosity > -1) {
			const Match::matchResult * mr = match->getSpokeStageResults();
			cout << " " << i << ": ";
			cout << "ObjFunc:" << getLastSpokeObjectiveFunctionValue(false);
			for (int j = 0; j < MAX_NUM_SPOKE_MATCH_RESULTS; j++) {
				bpTune_t index = (bpTune_t) mr[j].tuningParm;
				switch (index) {	// Why not just print the "short name"?
case BpSpokeImageMatch:				cout << " Match";
	break;
case BpSpokeLandmarkMatch:			cout << " Landmark";
	break;
case BpSpokeAverageNeighborMatch:	cout << " Reg";
	break;
case BpSpokeModelMatch:				cout << " Ref";
	break;
case SpokeMahalanobisMatch:			cout << " Maha";
	break;
case BpSpokeInterpenetrationMatch:	cout << " Pen";
	break;
case BpSpokeRSradPenalty:			cout << " Srad";
	break;
				}
				cout << ':' << mr[j].value;
			}
			cout << '\n';
		}
		if (optimizerStopRequested == true)
			break;
	}

	// IS THIS DOING ANYTHING? THE PURPOSE WAS TO TELL THE USER THAT MORE
	// OPTIMIZATION IS REQUIRED THAN WAS ALLOWED, EG, THE OBJFUNC DID NOT
	// LEVEL OUT YET
	if (i == (maxNumIterations - 1))
		if (globalVerbosity > -1)
			cout << "\nOptimizer could not reach an optimum in "
			<< maxNumIterations << " iterations" << endl;

	if (globalVerbosity > 0) {
		const Match::matchResult * mr = match->getSpokeStageResults();
		cout << "\nMatch after Spoke Stage:\n";
		// Note: BP omits printing of the penetration penalty
		for (int j = 0; j < MAX_NUM_SPOKE_MATCH_RESULTS - 1; j++) {
			bpTune_t index = (bpTune_t) mr[j].tuningParm;
			switch (index) {	// Why not just print the "short name"?
case BpSpokeImageMatch:				cout << " Match Distance";
	break;
case BpSpokeLandmarkMatch:			cout << " Landmark Distance";
	break;
case BpSpokeAverageNeighborMatch:	cout << " Match Regularity";
	break;
case BpSpokeModelMatch:				cout << " Match Reference";
	break;
case SpokeMahalanobisMatch:			cout << " Mahalanobis Distance";
	break;					
case BpSpokeInterpenetrationMatch:	cout << " Pen";
	break;
case BpSpokeRSradPenalty:			cout << " Srad";
	break;
			}
			cout << ": " << mr[j].value << '\n';
		}
		cout << flush;
	}

	int figureId = (int) tuningWt(BpFigureId);
	(void) MatchUtility::calculateAverageDistance(object, object->loadedObject(),
		figureId, getImageDistMapPointer(), (int) tuningWt(BpSpokeSurfaceLevel),
		(enum PrimNeighborhoodDefn) (int) tuningWt(BpSpokeNeighborhood),
		(enum DistanceType) (int) tuningWt(BpSpokeDistanceType),
		globalVerbosity);

	return true;
}


void P3DControl::binaryPabloNextOptimizationStage(int surfaceLevel)
{
	// Make certain we have an object
	if (! object) {
		cerr << "A model must first be loaded" << endl;
		return;
	}

	// Make certain we have an image
	if (! haveImage()) {
		cerr << "An image must first be loaded" << endl;
		return;
	}

	// First handle the stage we are leaving
	switch (currentOptimizerStage)
	{
	case NO_OPTIMIZER_STAGE:
	case OPTIMIZER_SETUP_STAGE:
		if (! initializeMatch(BINARY_IMAGE_MATCH, surfaceLevel))
			return;
		break;

	case OPTIMIZER_OBJECT_STAGE:
		break;

	case OPTIMIZER_ENSEMBLE_STAGE:
		break;

	case OPTIMIZER_FIGURE_STAGE:
		break;

	case OPTIMIZER_MAIN_FIGURE_STAGE:
		break;

	case OPTIMIZER_SUBFIGURE_STAGE:
		break;

	case OPTIMIZER_PRIMITIVE_STAGE:
		break;

	case OPTIMIZER_SREP_STAGE:
		break;


	case OPTIMIZER_SPOKE_STAGE:
		break;

	case OPTIMIZER_BOUNDARY_DISPLACEMENT_STAGE:
		break;

	case OPTIMIZER_END_STAGE:
		break;

	default:
		cout << "WARNING: Current optimizer stage (" << currentOptimizerStage << ") is not defined:" << endl;
		break;
	}

	// Get the next stage we are entering (depends on model)
	currentOptimizerStage = nextOptimizerStage();
	//cout << "Entering stage " << currentOptimizerStage << endl;
}



// xiaojie dibyendu

char* calcAntialiasedDMFilename(char* _imageFilename) {

	// if image file name is image.raw3
	// DM file name is image-Antialiased.raw3

	int len = strlen(_imageFilename);

	// input file must end in ".raw3", ".raw3z", or one of the following formats
    //  derived file will have ".raw3"

    char * suffix = NULL;

    if ((strstr(_imageFilename, ".raw3" ) == (_imageFilename + len - strlen(".raw3" ))))
        suffix = ".raw3";

    if ((strstr(_imageFilename, ".raw3z" ) == (_imageFilename + len - strlen(".raw3z" ))))
        suffix = ".raw3z";

    if ((strstr(_imageFilename, ".gipl" ) == (_imageFilename + len - strlen(".gipl" )))) 
        suffix = ".gipl";

    if ((strstr(_imageFilename, ".mhd" ) == (_imageFilename + len - strlen(".mhd" ))))
        suffix = ".mhd";

    if ((strstr(_imageFilename, ".mha" ) == (_imageFilename + len - strlen(".mha" ))))
        suffix = ".mha";

    if (!suffix)
        return NULL;

    int suffixLen = strlen(suffix);    

	const char* ddmFileExtension = ".raw3";	

	char* ddmFilename = new char[ len + 20 + 1 ];

	strcpy(ddmFilename, _imageFilename);

	strcpy(ddmFilename + len - suffixLen, "-antiAliased");

	strcat(ddmFilename, ddmFileExtension);

	return ddmFilename;
}

// dibyendu

char* calcOutputDMFilename(char* _imageFilename) {

	// if image file name is image.raw3
	// DM file name is image-Antialiased.raw3

	int len = strlen(_imageFilename);

	// input file must end in ".raw3", ".raw3z", or one of the following formats
    //  derived file will have ".raw3"

    char * suffix = NULL;

    if ((strstr(_imageFilename, ".raw3" ) == (_imageFilename + len - strlen(".raw3" ))))
        suffix = ".raw3";

    if ((strstr(_imageFilename, ".raw3z" ) == (_imageFilename + len - strlen(".raw3z" ))))
        suffix = ".raw3z";

    if ((strstr(_imageFilename, ".gipl" ) == (_imageFilename + len - strlen(".gipl" )))) 
        suffix = ".gipl";

    if ((strstr(_imageFilename, ".mhd" ) == (_imageFilename + len - strlen(".mhd" ))))
        suffix = ".mhd";

    if ((strstr(_imageFilename, ".mha" ) == (_imageFilename + len - strlen(".mha" ))))
        suffix = ".mha";

    if (!suffix)
        return NULL;

    int suffixLen = strlen(suffix);    

	const char* ddmFileExtension = ".raw3";	

	char* ddmFilename = new char[ len + 20 + 1 ];

	strcpy(ddmFilename, _imageFilename);

	strcpy(ddmFilename + len - suffixLen, "-outputDM");

	strcat(ddmFilename, ddmFileExtension);

	return ddmFilename;
}


// derive the ddm image filename from the image filename, eg,
// return "I-ddm.raw3" if imageFilename is "I.raw3", eg, in
// cshell, $a:r-ddm.$a:e


// DEBUG THIS LATER, IF TIME

char * calcDDMFilename(char * _imageFilename, unsigned short mask = 0xFFFF)
{
	char * ddmFilename;		// the derived filename
	int len = strlen(_imageFilename);
	const char* ddmFileExtension = ".raw3";

	// input file must end in ".raw3" or ".raw3z"; derived file will have raw3
	char * suffix = NULL;
	if ((strstr(_imageFilename, ".raw3" ) == (_imageFilename + len - strlen(".raw3" ))))
		suffix = ".raw3";
	if ((strstr(_imageFilename, ".raw3z" ) == (_imageFilename + len - strlen(".raw3z" ))))
		suffix = ".raw3z";
	if ((strstr(_imageFilename, ".gipl" ) == (_imageFilename + len - strlen(".gipl" ))))
		suffix = ".gipl";
	if ((strstr(_imageFilename, ".mhd" ) == (_imageFilename + len - strlen(".mhd" ))))
		suffix = ".mhd";
	if ((strstr(_imageFilename, ".mha" ) == (_imageFilename + len - strlen(".mha" ))))
		suffix = ".mha";
	if (!suffix)
		return NULL;

	int suffixLen = strlen(suffix);
	ddmFilename = new char[len+suffixLen+1+7];		// insert '-ddm' -0xABCD
	strcpy(ddmFilename, _imageFilename);
	strcpy(ddmFilename + len - suffixLen, "-ddm");
	if (mask != 0xFFFF) {
		sprintf(ddmFilename + len - suffixLen + 4, "-0x%0.4x", mask);
	}
	strcat(ddmFilename, ddmFileExtension);
	return ddmFilename;
}


#ifdef INTEGRATE_WHEN_OBJ_HISTORY_STRING_DONE
// Append nameUpdate to object's history,
// if nameUpdate is not already in the parent
void appendHistoy(M3DObject * obj, const char * nameUpdate)
{
	char newName[1000];
	char * token;		// is newName already in name?

	if (obj && nameUpdate && *nameUpdate) {
		const char *objName = obj->getName();
		token = strstr(objName, nameUpdate);
		if (! token) {
			if (objName)
				strcpy(newName, objName);
			else
				newName[0] = '\0';
			strcat(newName, ":");
			strcat(newName, nameUpdate);
			obj->setName(newName);
			if (globalVerbosity > 0)
				cout << " renaming model object to " << newName << endl;
		}

		M3DFigure * fig = obj->getFigurePtr(tuningWt(BpFigureId)); //Multiobject
		if (fig) {
			const char * figName = fig->getName();
			token = strstr(figName, nameUpdate);
			if (! token) {
				if (figName)
					strcpy(newName, figName);
				else
					newName[0] = '\0';

				strcat(newName, ":");
				strcat(newName, nameUpdate);
				fig->setName(newName);
				if (globalVerbosity > 0)
					cout << " renaming model's first figure to " << newName << endl;
			}
		}
	}
}
#endif



#include <algorithm>

struct TPSPoint {
	double x,y,z;	// The point in xyz coords
	int m;			// The closest model
	double dist;	// The distance to the closest model
	double u,v,t;	// The uvt coords of the closest point on the surface.
	TPSPoint() : dist(HUGE_VAL), m(-1) {}

	TPSPoint& operator = (const TPSPoint& a) {
		x	= a.x;
		y	= a.y;
		z	= a.z;
		m	= a.m;
		dist	= a.dist;
		u	= a.u;
		v	= a.v;
		t	= a.t;
		return (*this);
	}
};


bool operator < (const TPSPoint& a, const TPSPoint& b) {
	return (a.m < b.m) || ( (a.m == b.m) &&
		((a.u < b.u) || ( (a.u == b.u)  &&
		((a.v < b.v) || ( (a.v == b.v) && (a.t < b.t))))));
}

bool operator == (const TPSPoint& a, const TPSPoint& b) {
	return (a.m == b.m) && (a.u == b.u) && (a.v == b.v) && (a.t == b.t);
}

int P3DControl::runThinPlateSpline(Registry &registry)
{
	const int subdivLevel		= registry.getIntValue("Subdivisionlevel",3);
	const char* inPointsfName	= registry.getStringValue("Inpoints","");
	const char* outPointsfName	= registry.getStringValue("Outpoints","");
	// these are fixed to data<stepnum>.csv
	const int numModels			= registry.getIntValue("NumModels",0);
	const int numSteps			= registry.getIntValue("NumSteps", 2);
	// Image is needed to convert from model to image coords and vice versa.
	const char* imgName			= registry.getStringValue("Image",NULL);
	const bool saveInWorldCoordinates
		= registry.getBooleanValue("SaveInWorldCoordinates",false);
	const char* targetImageName	= registry.getStringValue("TargetImage",NULL);

	std::vector<const char*> srcModelNames(numModels), dstModelNames(numModels);
	std::vector<M3DObject*> srcObjects(numModels), dstObjects(numModels), interpObjects(numModels);
	std::vector<ThallCode::Pointlist_server2*> pList(numModels);

#if 1
	// Keeps a count of total number of projections on each model
	std::vector<int> projectionCount(numModels,0);
#endif

	M3DObjectFile m3dfile;
	Xferlist* xfer;
	char buffer[1000];
	int p,m;

	int nSrc, nDst;

	if(!imgName || !*imgName) {
		cerr << "ERROR: Loading image " << endl;
		return 7;
	}

	if( globalVerbosity > -1 ) {
		cout << "Loading image " << imgName << " ... " << endl;
	}
	Image3D* image	= RAWImageFile().read(imgName);
	if( image == NULL ) {
		cerr << "ERROR: Loading image " << imgName << ", reason: " << strerror(errno) << endl;
		return 8;
	}

	Image3D* targetImage	= NULL;
	if(saveInWorldCoordinates) {
		if(!targetImageName || !*targetImageName) {
			cerr << "ERROR: Loading target image " << endl;
			return 7;
		}

		if( globalVerbosity > -1 ) {
			cout << "Loading target image " << targetImageName << " ... " << endl;
		}
		targetImage	= RAWImageFile().read(targetImageName);
		if( targetImage == NULL ) {
			cerr << "ERROR: Loading target image " << targetImageName << ", reason: " << strerror(errno) << endl;
			return 8;
		}
	}

	for(m = 0; m != numModels; ++m ) {
		// Get Model names
		srcModelNames[m]		= registry.getStringValue("Srcmodel[%d]", NULL,m);
		// Load the models
		if( globalVerbosity > -1 ) {
			cout << "Loading src model " << srcModelNames[m] << " ... ";
		}
		srcObjects[m]			= m3dfile.read(srcModelNames[m]);
		if(srcObjects[m] == NULL) {
			cerr << "ERROR: Could not load src model " << srcModelNames[m] << endl;
			return 1;
		}
		// Compute the surface point cloud.
		if( globalVerbosity > -1 ) {
			cout << "computing point cloud ... ";
		}
		xfer					= convertM3DtoXfer(srcObjects[m]->getFigurePtr(0));
		pList[m]				= new ThallCode::Pointlist_server2();
		pList[m]->init(xfer);
		pList[m]->ComputeSubdivPointCloud(subdivLevel);
		delete[] xfer->atomlist;
		delete xfer;
		if( globalVerbosity > -1 ) {
			cout << "done" << endl;
		}

		interpObjects[m]		= srcObjects[m]->clone();

		// Get Model names
		dstModelNames[m]		= registry.getStringValue("Dstmodel[%d]", NULL,m);
		// Load the models
		if( globalVerbosity > -1 ) {
			cout << "Loading dst model " << dstModelNames[m] << " ... ";
		}
		dstObjects[m]			= m3dfile.read(dstModelNames[m]);
		if(dstObjects[m] == NULL) {
			cerr << "ERROR: Could not load dst model " << dstModelNames[m] << endl;
			return 2;
		}
		if( globalVerbosity > -1 ) {
			cout << "done" << endl;
		}

		// Check if primitive counts match ...
		nSrc	= srcObjects[m]->getPrimitiveCount();
		nDst	= dstObjects[m]->getPrimitiveCount();
		if( nSrc != nDst ) {
			cerr << "ERROR: Primitive counts do not match: " << srcModelNames[m] << "(" << nSrc << ") != "
				<< dstModelNames[m] << "(" << nDst << ")" << endl;
			return 3;
		}
	}

	// Now load the points which have to be registered - (x,y,z) triplets
	vector<TPSPoint> points;
	std::ifstream inPoints;
	std::ofstream outPoints;
	if (inPointsfName[0] != '\0') {
		inPoints.open(inPointsfName);
		if( inPoints.bad() ) {
			cerr << "ERROR: Could not load inpoints file " << inPointsfName << ", reason: " << strerror(errno) << endl;
			// Should cleanup here ... oh well ...
			return 4;
		}
		outPoints.open(outPointsfName);
		if( outPoints.bad() ) {
			cerr << "ERROR: Could not load outpoints file " << outPointsfName << ", reason: " << strerror(errno) << endl;
			return 5;
		}
	}
	double dist;
	int multiplicity;
	ThallCode::DbVector3 nearpointXYZ;
	ThallCode::DbVector3 nearpointUVT;
	Vector3D coord;

	Bpoint bp;
	std::vector< std::vector<DistanceVector> > pointlist(numModels);

	if( registry.getIntValue("BackProject", 0) == 0 ) {
		if( inPointsfName[0] != '\0') {
			for( p = 0; !inPoints.eof(); ++p ) {
				if( p % 1000 == 0) {
					cout << "\rTransforming ... " << p;
					cout.flush();
				}
				double x,y,z;
				inPoints >> x >> y >> z;
				if( inPoints.eof() ) {
					break;
				}
				// Pablo internally uses a left-handed coordinate system.
				// These byu and image files are all in a right-handed system. When such a system
				// is given to pablo (the width in some direction being -ve), it moves the origin
				// to the other side of the axis.
				// need to do this world to model coords transform differently.
				// FIXME: Never works across two images :'(
				coord.set(x,y,z);
				image->imageToModelCoordinates(coord);

				// Converting back to world scale for target image:
				if(saveInWorldCoordinates) 
					targetImage->modelToWorldCoordinates(coord);
				outPoints << coord.getX() << '\t' << coord.getY() << '\t' << coord.getZ() << endl;
			}
			cout << endl;
		}
		for( int m = 0; m != numModels; ++m) {
			int numPoints;
			Bpoint *surfacePointList;
			pList[m]->subdivboundaryinfo(&numPoints, &surfacePointList);
			projectionCount[m]	= numPoints;
			pointlist[m].resize(projectionCount[m]);
			projectionCount[m]	= 0;
			for( p = 0; p != numPoints; ++p ) {
				if( !pList[m]->thisMesh->isEndMedcoord(ThallCode::Medcoord(surfacePointList[p])) ) {
					DistanceVector& dv	= pointlist[m][p];
					dv.p.set(surfacePointList[p].pnt[0], surfacePointList[p].pnt[1], surfacePointList[p].pnt[2]);
					dv.grad.set(0.001,0.000,0.000);
					dv.n.set(0.0,0.001,0.0);
					dv.method	= 1;
					dv.dist		= 1.0;
					TPSPoint pnt;
					pnt.dist	= 1.0;
					pnt.m		= m;
					pnt.u		= surfacePointList[p].u;
					pnt.v		= surfacePointList[p].v;
					pnt.t		= surfacePointList[p].t;
					pnt.x		= surfacePointList[p].pnt[0];
					pnt.y		= surfacePointList[p].pnt[1];
					pnt.z		= surfacePointList[p].pnt[2];
					points.push_back( pnt );
					projectionCount[m]++;
				}
			}
			pointlist[m].resize(projectionCount[m]);
		}
	}
	else {
		for( p = 0; !inPoints.eof(); ++p ) {
			if( p % 100 == 0) {
				cout << "\rBack projecting ... " << p;
				cout.flush();
			}
			double x,y,z;
			inPoints >> x >> y >> z;
			if( inPoints.eof() ) {
				break;
			}
			// Pablo internally uses a left-handed coordinate system.
			// These byu and image files are all in a right-handed system. When such a system
			// is given to pablo (the width in some direction being -ve), it moves the origin
			// to the other side of the axis.
			// need to do this world to model coords transform differently.
			// FIXME: Never works across two images :'(
			coord.set(x,y,z);
			image->imageToModelCoordinates(coord);
			// Back-project the point onto the models (find the closest model in the process too)
			points.push_back( TPSPoint() );
			for( m = 0; m != numModels; ++m ) {
				dist = pList[m]->subdivNearpoint(coord.getX(), coord.getY(), coord.getZ(), &nearpointXYZ, &nearpointUVT, &multiplicity);
				if( dist < points[p].dist || points[p].m == -1 ) {
					points[p].dist	= dist;
					points[p].m		= m;
					points[p].u		= nearpointUVT.x();
					points[p].v		= nearpointUVT.y();
					points[p].t		= nearpointUVT.z();
					points[p].x		= nearpointXYZ.x();
					points[p].y		= nearpointXYZ.y();
					points[p].z		= nearpointXYZ.z();
				}
			}
#if 1
			projectionCount[points[p].m]++;
#endif
			// Converting back to world scale for target image:
			if(saveInWorldCoordinates)
				targetImage->modelToWorldCoordinates(coord);
			outPoints << coord.getX() << '\t' << coord.getY() << '\t' << coord.getZ() << endl;
		}
		std::sort( points.begin(), points.end() );
		points.erase( std::unique( points.begin(), points.end() ), points.end() );

#if 1
		for( m = 0; m != numModels; ++m ) {
			pointlist[m].resize(projectionCount[m]);
			projectionCount[m]	= 0;
		}
		for( p = 0; p != points.size(); ++p ) {
			DistanceVector& dv	= pointlist[points[p].m][projectionCount[points[p].m]++];
			/*	if(!pList[points[p].m]->subdivBposition( &bp, points[p].u, points[p].v, points[p].t)) {
			cerr << "ERROR: " << __FILE__ << ':' << __LINE__ << " p=" << p << " m=" << m << endl
			<< "uvt=" << points[p].u << ',' << points[p].v << ',' << points[p].t << endl;
			return 6;
			}
			outPoints << bp.pnt[0] << ", " << bp.pnt[1] << ", " << bp.pnt[2] << ", " << points[p].dist
			<< ", " << points[p].m << endl;
			dv.p.set(bp.pnt[0], bp.pnt[1],	bp.pnt[2]);
			*/	
			dv.p.set(points[p].x,points[p].y,points[p].z);
			dv.grad.set(0.001,0.000,0.000);
			dv.n.set(0.0,0.001,0.0);
			dv.method	= 1;
			dv.dist		= 1.0;
		}
		for( m = 0; m != numModels; ++m ) {
			pointlist[m].resize(projectionCount[m]);
		}
		cout << endl;
	}
	for( m = 0; m != numModels; ++m ) {
		sprintf(buffer,"%s.list",srcModelNames[m]);
		std::ofstream plistout(buffer, ios::out|ios::trunc);
		std::copy(pointlist[m].begin(), pointlist[m].begin() + projectionCount[m],
			std::ostream_iterator<DistanceVector>(plistout));
		plistout.close();
	}
#endif
	// Now do the geodesic stepping for generating all the landmarks ...
	M3DPrimitive *prim, *srcP, *dstP;
	for( int step_num = 1; step_num <= numSteps; ++step_num ) {
		sprintf(buffer, "csv/Data%d.csv", step_num );
		std::ofstream outPoints(buffer);
		if( outPoints.bad() ) {
			cerr << "ERROR: Could not load outpoints file " << buffer << ", reason: " << strerror(errno) << endl;
			// Should cleanup here ... oh well ...
			return 5;
		}
		// step over models ...
		cout << "\r" << "Step " << step_num << " of " << numSteps;
		cout.flush();
		for( m = 0; m != numModels; ++m ) {
			nSrc	= srcObjects[m]->getPrimitiveCount();
			// step over primitives ...
			for( p = 0; p != nSrc; ++p ) {
				srcP	= srcObjects[m]->getPrimitivePtr(p);
				dstP	= dstObjects[m]->getPrimitivePtr(p);
				prim	= srcP->copyPtr();
				prim->atomInterp((step_num - 1.0)/(numSteps - 1.0), srcP, dstP);
				*(interpObjects[m]->getPrimitivePtr(p))	= *prim;
			}
#if 1
			// Save the intermediate models
			const char* filename	= strrchr(srcModelNames[m], '/');
			if( filename == NULL ) {
				filename	= strrchr(srcModelNames[m], '\\');
			}
			if( filename == NULL ) {
				filename	= srcModelNames[m];
			}
			sprintf(buffer,"debug/%02d-%s",step_num, filename);
			//			m3dfile.write(buffer, *interpObjects[m] );
			projectionCount[m]	= 0;
#endif
			xfer		= convertM3DtoXfer(interpObjects[m]->getFigurePtr(0));
			delete pList[m];
			pList[m]	= new ThallCode::Pointlist_server2();
			pList[m]->init(xfer);
			pList[m]->ComputeSubdivPointCloud(subdivLevel);
			delete[] xfer->atomlist;
			delete xfer;
		}
		// Now step over all points and find the new xyz for the same uvt
		for( p = 0; p != points.size(); ++p ) {
			if(!pList[points[p].m]->subdivBposition( &bp, points[p].u, points[p].v, points[p].t)) {
				cerr << "ERROR: " << __FILE__ << ':' << __LINE__ << " p=" << p << " m=" << m << endl
					<< "uvt=" << points[p].u << ',' << points[p].v << ',' << points[p].t << endl;
				return 6;
			}
			//			outPoints << bp.pnt[0] << ", " << bp.pnt[1] << ", " << bp.pnt[2] << ", " << points[p].dist
			//				<< ", " << points[p].m << endl;
#if 1
			DistanceVector& dv	= pointlist[points[p].m][projectionCount[points[p].m]++];
			dv.p.set(bp.pnt[0],bp.pnt[1],bp.pnt[2]);
			dv.grad.set(0.001,0.000,0.000);
			dv.n.set(0.0,0.001,0.0);
			dv.method	= 1;
			dv.dist		= 1.0;
#endif
			if(saveInWorldCoordinates) {
#if 1
				//bp.pnt[1]	= -bp.pnt[1];
#endif
				bp.pnt[0] = targetImage->modelXToWorld(bp.pnt[0]);
				bp.pnt[1] = targetImage->modelYToWorld(bp.pnt[1]);
				bp.pnt[2] = targetImage->modelZToWorld(bp.pnt[2]);
			}
			outPoints << bp.pnt[0] << ", " << bp.pnt[1] << ", " << bp.pnt[2] << endl;
		}
#if 1
		for( m = 0; m != numModels; ++m ) {
			const char* filename	= strrchr(srcModelNames[m], '/');
			if( filename == NULL ) {
				filename	= strrchr(srcModelNames[m], '\\');
			}
			if( filename == NULL ) {
				filename	= srcModelNames[m];
			}
			sprintf(buffer,"debug/%02d-%s.list",step_num, filename);
			std::ofstream plistout(buffer, ios::out|ios::trunc);
			std::copy(pointlist[m].begin(), pointlist[m].begin() + projectionCount[m],
				std::ostream_iterator<DistanceVector>(plistout));
			plistout.close();
		}
#endif
	}

	for( m = 0; m != numModels; ++m ) {
		delete srcObjects[m];
		delete dstObjects[m];
		delete interpObjects[m];
		delete pList[m];
	}
	return 0;
}


bool P3DControl::runBinaryPablo(Registry & script)
{
	bool status;
	const char * filename;
	const char * modelFilename;
	const char * xformFilename;
	const char * pgaFilename;

	// dibyendu
	string spokeDistFileName ;

	// dibyendu - making a video of the optimization process

	const char * videoFolderName = script.getStringValue("VideoFolder", NULL) ;

	if( videoFolderName ) {
		// MakeOptimizationVideo::turnVideoOn( (char *) videoFolderName ) ;
		MakeOptimizationVideo::makeVideoOn = 1 ;
		MakeOptimizationVideo::videoFolder = (char *) videoFolderName  ;
		// MakeOptimizationVideo::frameCount = 0 ; 
	} 
	else {
		//MakeOptimizationVideo::turnVideoOff() ;
		MakeOptimizationVideo::makeVideoOn = 0 ;
		MakeOptimizationVideo::videoFolder = NULL ;
	}

	videoFolderName = NULL ;
		


	if (globalVerbosity > -1)
		cout << "---- LOADING FILES" << endl;

	// Optionally, load a PGA mean model File (command-line option -ipr)
	filename = script.getStringValue("PgaResidueModel", NULL);
	if (filename != NULL) {
		pgaModelFileName = new char[strlen(filename) + 1];
		strcpy(pgaModelFileName, filename);
		if (globalVerbosity > -1)
			cout << "PGA mean model: " << pgaModelFileName << '\n';
	}

	// Get the replacement PGA model file (command-line option -ip)
	pgaFilename = script.getStringValue("PgaModel", NULL);

	// Load the model.  If pgaModelFileName contains a path, the
	// PgaResidueModel is also loaded.  If the PgaModel is specified,
	// its statistics will replace those in the model.
	modelFilename = script.getStringValue("Model", NULL);
	xformFilename = script.getStringValue("Transform", NULL);
	// Scripts cannot map to world
	if (! read(modelFilename, false, xformFilename, pgaFilename)) {
		if (globalVerbosity > -1) {
			if (modelFilename == NULL) 
				cout << "Error: No model file specified\n";
			else
				cout << "Error: Unable to load model: " << modelFilename << '\n';
		}
		return false;
	}
	else
		if (globalVerbosity > -1)
			cout << "Loaded model " << modelFilename << '\n';
	binaryPabloListLandmarks();
	selectAll();

	// Load the Landmark Model (optional)
	filename = script.getStringValue("LandmarkModel", NULL);
	if (filename && *filename) {
		// Landmark model is optional
		status = binaryPabloLoadLandmarkModel(filename);
		if (! status)
			return false;
	}

	// Load the image
	filename = script.getStringValue("Image", NULL);

	char * _imageFilename;	// Save to calculated distance map's filename

	if (filename && *filename) {
		_imageFilename = strdup(filename);
		if (globalVerbosity > -1)
			cout << "Loading image " << filename << endl;

		status = loadImage(filename, false);	// Scripts cannot map to world
	}
	else {
		if (globalVerbosity > -1)
			cout << "Error: An image file is required for optimization" << endl;
		status = false;
	}

	cout << "Image loaded" << endl;

	// If outModelImageFile is defined, export to image
	filename = script.getStringValue("OutModelImage", NULL);
	if (filename && *filename) {
		exportToImage(filename);
		return true;    // Successful
	}

	if (! status)
		return false;

	if( filename != NULL )
		delete [] filename ;
	filename  = NULL ;

	Image3D * image = getImagePtr();
	unsigned short attMask = (unsigned short) tuningWt(AttractionMask);
	image->pushImageIsStacked(attMask != 0xFFFF, attMask);


	// LOAD OR COMPUTE THE DISTANCE MAP
	// - look for the file specified in the registry
	// - if absent, look for I-ddm.raw3 (if the image was I.raw3)
	// - if the image is stacked, the filename to look for is
	//   I-ddm-0xABCD.raw3 where ABCD is hex code for the attraction mask
	// - DDMEnable==3 disables reading or calc'ing the distance map, on the
	//   assumption that it's not needed.	

	if (((int) tuningWt(BpDDMEnable) != 3))
	{
		char * deallocate_filename	= NULL;

		// Load the filename from the command line input or config script
		filename = script.getStringValue("DistanceMap", NULL);	

		if (! filename) {
			if( (bool) tuningWt(BpDDMAntialiased) ) {	
				// distance map filename = "<filename>-antiAliased.raw3"
				deallocate_filename = calcAntialiasedDMFilename( _imageFilename ) ;
			}
			else {
				// distance map filename = "<filename>-ddm.raw3"
				deallocate_filename = calcDDMFilename(_imageFilename, attMask);
			}
			
			filename = deallocate_filename;

			if (globalVerbosity > -1)
				cout << "Distance map image: " << filename << endl;
		}
		
		


		// deepak: 
		if( _imageFilename )
		{
			free(_imageFilename);		
			_imageFilename = NULL;
		}

		// We can read and write this signed raw3 image, but users should
		// NOT try to display them...
#define SIGNED_RAW3
#ifdef SIGNED_RAW3
		//Xiaojie		

		if (((int) tuningWt(BpIgnoreDistFile) == 1) || ( ! loadDistMap(filename) ) ) {
#else
			cout << "WARNING: file ignored; signed distance map files are not supported in this version" << endl;
#endif
			distMap(filename);	//load OR compute & save the distance map

			cout << "Loaded distance map " << filename << endl ;
		}

		if( deallocate_filename )
			delete[] deallocate_filename;
		deallocate_filename = NULL ;
		
		filename = NULL ;

	}
	else
		cout << "Distance Map disabled by DDMEnable" << endl;

#define DEBUG_DIBYENDU_DIST_MAP
#ifdef DEBUG_DIBYENDU_DIST_MAP

	ImageDistanceMap * donotDeleteDMap = NULL ;

	const int * mapSize = NULL ;

	if( match && match->getImageDistMapPointer() ) {
		mapSize = match->getImageDistMapPointer()->getMapSize() ;
		donotDeleteDMap = match->getImageDistMapPointer() ;
	}
	else if( dist_map ) {
		mapSize = dist_map->getMapSize() ;
		donotDeleteDMap = dist_map ;
	}
	else
		cout << "NULL Image distance map pointer !" << endl ;

	if( mapSize && donotDeleteDMap ) {

		const char* debugFilename = "../../bin/test_dist_map.txt" ;

		std::ofstream f ;
		f.open( (const char*) debugFilename,  ios::trunc);	

		if( f ) {
			f << setiosflags(ios::fixed) ;
			cout << debugFilename << " has been opened for writing distance map" << endl ;

			f << setw(12) << setprecision(8) << mapSize[0] << " " << mapSize[1] << " " << mapSize[2] << endl ;

			for( int k = 0 ; k < mapSize[2] ; k++ ) {
				for( int i = 0 ; i < mapSize[0] ; i++ ) {
					for( int j = 0 ; j < mapSize[1] ; j++ ) 
							f << setw(12) << setprecision(8) << donotDeleteDMap->getDistance(i, j, k) ;		

					f << endl ;					
				}
				f << endl ;
			}
		}
		f.close() ;
	}

#endif 

	// Seed the random number generator before optimizing
	srand((int) tuningWt(OptimizerRandomSeed));

	// INITIALIZATION STAGE: use method of moments to CALC A TRANSFORM
	// THAT WOULD place the figure in the image (if it had been applied),
	// but do NOT apply the transform here...
	

	if (globalVerbosity > -1)
		cout << "---- INITIALIZATION STAGE" << endl;

	GLOBALMOMTransform = MatchUtility::initializeObjectToImage(object, getImagePtr(),
		(int) tuningWt(BpDoMethodOfMoments), globalVerbosity);
		

	const char* outSimTransModelFile = script.getStringValue("OutSimTransModel", NULL);

	if (outSimTransModelFile && *outSimTransModelFile) {
		cout << "Writing similarity transform to " << outSimTransModelFile << endl;
		GLOBALMOMTransform.writeSimilarity(outSimTransModelFile, false);
	}

	if( MakeOptimizationVideo::makeVideoOn )			
		MakeOptimizationVideo::writeCurrentModel( object ) ;	// dibyendu - save current model		

	
	// Now dilate the object, if requested
	
	// Convert units for dilationFactor.
	// (Distance map needs to have uniform imageToModelScale)
	// FIXME: Check for this condition.

	// it sure would be nice if *distance* image access methods were available instead of just the *coordinate* methods
	double dilationFactorInModelUnits	= 0;
	if (tuningWt(BpDilationFactor) != 0) {
		double imageToModelScale = 1.0 / getDistMapPtr()->getModelToImageScale()[0];
		dilationFactorInModelUnits = tuningWt(BpDilationFactor) * imageToModelScale;
		getObjectPtr()->dilate(dilationFactorInModelUnits);
		
	}

	if (int(tuningWt(BpSubdivide)) != 0) {
		cout << "Subdividing model " << int(tuningWt(BpSubdivide)) << " times" << endl;
		for( int i = 0; i < int(tuningWt(BpSubdivide)); ++i) {
			getObjectPtr()->subdivide();
			getObjectPtr()->loadedObject()->subdivide();
		}
	}

	// FIGURE STAGE: reset loaded model and optimize
	if (globalVerbosity > -1)
		cout << "---- FIGURE STAGE" << endl;

	M3DQuadFigure * qFigThis = dynamic_cast <M3DQuadFigure *>(object->getFigurePtr(0)) ;
	
	double * sRepVector = qFigThis->vectorize() ;	

	int numAtoms = qFigThis->getPrimitiveCount() ;

	// separate containers for positions, spoke directions and spoke radii

	double * positions = new double[ 3*numAtoms ] ;
	for( int i = 0 ; i < 3*numAtoms ; i++ )
		positions[i] = 0.0 ;

	for( int i = 0 ; i < numAtoms ; i++ ) {
		positions[ 3*i + 0 ] = sRepVector[ 3*i + 0 ] ;
		positions[ 3*i + 1 ] = sRepVector[ 3*i + 1 ] ;
		positions[ 3*i + 2 ] = sRepVector[ 3*i + 2 ] ;
	}

	if( sRepVector != NULL ) 
		delete [] sRepVector ;
	sRepVector = NULL ;

	double meanX = 0.0, meanY =0.0 , meanZ =0.0 ;

	for( int i = 0 ; i < numAtoms ; i++ ) {
		meanX += positions[ 3*i + 0 ] ;
		meanY += positions[ 3*i + 1 ] ;
		meanZ += positions[ 3*i + 2 ] ;
	}
	
	meanX /= numAtoms ;
	meanY /= numAtoms ;
	meanZ /= numAtoms ;

	double *newMean = new double[3];
	newMean[0] = meanX;
	newMean[1] = meanY;
	newMean[2] = meanZ;

	if (object->getCPNSStats() != NULL)
		object->getCPNSStats()->setMeanShape(newMean);

	bool windowIntens = script.getBooleanValue("WindowIntensities", false);
	if (globalVerbosity > 0 && ! windowIntens)
		cout << "Intensities will be windowed in optimization\n";

	setOptimizerIntensityWindowing(windowIntens);	// Must precede call of initializeMatch()	

	binaryPabloNextOptimizationStage((int) tuningWt(BpFigureSurfaceLevel));	
	
		
	

	if (tuningWt(BpFigureResetModelBefore) == 1.0) {
		M3DFigure *fig = object->getFigurePtr((int) tuningWt(BpFigureId));

#ifdef UNFLIPPED
		Vector3D cog = fig->getCOG();
		cog.setY(1.0 - cog.getY());
		fig->applySimilarity(GLOBALMOMTransform, cog);
#else
		fig->applySimilarityAboutCOG(&GLOBALMOMTransform);
#endif

		// Also apply the transform to the original object
		M3DObject * lobj = object->loadedObject();
		fig = lobj->getFigurePtr((int) tuningWt(BpFigureId));
#ifdef UNFLIPPED
		cog = fig->getCOG();
		cog.setY(1.0 - cog.getY());
		fig->applySimilarity(GLOBALMOMTransform, cog);
#else
		fig->applySimilarityAboutCOG(&GLOBALMOMTransform);
#endif

		GLOBALMOMTransform.setToIdentity();
	}
	if (globalVerbosity >= 0) {
		if (tuningWt(BpFigureResetModelBefore) == 1.0)
			cout << "Reset reference model via MOM transform" << endl;
		else
			cout << "Did not reset reference model via MOM transform" << endl;
	}

	bool have_PGA = havePGA();
	cpnsStats = object->getCPNSStats();
	bool have_CPNS = (cpnsStats != NULL);
	bool wantPGA = (tuningWt(FigureMahalanobisMatch) != 0);

	if (wantPGA && !have_PGA && !have_CPNS) {
		cout << "Error: PGA/CPNS match weight set (" << tuning(FigureMahalanobisMatch).longName
			<< "=" << tuningWt(FigureMahalanobisMatch)
			<< ") but PGA/CPNS missing from model " << modelFilename << endl;
		return false;	// Wrong optimizer would be invoked, so fail
	}
	else if (! wantPGA && have_PGA) {
		cout << "Warning: model contains PGA stat's but not enabled with " <<
			tuning(FigureMahalanobisMatch).longName << endl;
	}

#define DEBUG_DIBYENDU_WRITE_DIST

#ifdef DEBUG_DIBYENDU_WRITE_DIST

	// dibyendu
	// writing out the distances at the spoke ends of all atoms

	if( modelFilename != NULL ) {

		spokeDistFileName = string( modelFilename ) + ".dist.txt" ;

		this->writeDistancesAtSpokeEnds( (const char *) spokeDistFileName.c_str() ) ;
	}

#endif // DEBUG_DIBYENDU_WRITE_DIST


	if (tuningWt(FigureIterations) > 0)
	{
		// Match has already been initialized, so it is safe to pass it the 
		// contour filename
		const char * filename = script.getStringValue("Contour", NULL);
		if (filename) {
			if (! match->initContourMatch(filename))
				cout <<
				"Warning: Contour penalty not specified; ignoring contour points file"
				<< endl;
			else
				cout << "Using contour file " << filename << endl;
		}
		else {
			if (tuningWt(FigureContourMatch) != 0.0) {
				cout << "Error: A contour file must be specified to use contour matching"
					<< endl;
				return false;
			}
		}

		binaryPabloFigureStage((int) tuningWt(FigureIterations),
			(int) tuningWt(BpFigureDoScale), modelFilename);


		if (outSimTransModelFile && *outSimTransModelFile) {
			cout << "Writing similarity transform to " << outSimTransModelFile << endl;
			GLOBALMOMTransform.writeSimilarity(outSimTransModelFile, false);
		}
	}

	// ATOM STAGE: reset reference model and optimize
	if (globalVerbosity > -1)
		cout << "---- ATOM STAGE" << endl;

	binaryPabloNextOptimizationStage((int) tuningWt(BpAtomSurfaceLevel));
	if (tuningWt(BpAtomResetModelBefore) == 1.0) {
		// Use fig stage's model for atomNeighbor and modelMatch distances

		// Replace the reference object with the transformed object.
		// Because the "loaded object" is never supposed to differ from
		// the model in the input file, this operation is bad practice.
		// However, it will suffice until a better solution is determined.
		M3DObject * lobj = object->loadedObject();
		*lobj = *object;
	}
	if (globalVerbosity >= 0) {
		if (tuningWt(BpAtomResetModelBefore) == 1.0)
			cout << "Reset reference model" << endl;
		else
			cout << "Did NOT reset reference model" << endl;
	}

	int n = (int) tuningWt(AtomIterations);
	if (n < 0)
		n = MAX_OPTIMIZER_ITERATIONS;
	if (n != 0)
		binaryPabloAtomStage((int) tuningWt(AtomIterations));

	// dibyendu

	// ----------------------------------------------------------------

	// SREP STAGE: reset reference model and optimize
	if (globalVerbosity > -1)
		cout << "---- SREP STAGE" << endl;

	binaryPabloNextOptimizationStage((int) tuningWt(BpSRepSurfaceLevel));
	if (tuningWt(BpSRepResetModelBefore) == 1.0) {
		// Use atom stage's model for atomNeighbor and modelMatch distances

		// Replace the reference object with the transformed object.
		// Because the "loaded object" is never supposed to differ from
		// the model in the input file, this operation is bad practice.
		// However, it will suffice until a better solution is determined.
		M3DObject * lobj = object->loadedObject();
		*lobj = *object;
	}
	if (globalVerbosity >= 0) {
		if (tuningWt(BpSRepResetModelBefore) == 1.0)
			cout << "Reset reference model" << endl;
		else
			cout << "Did NOT reset reference model" << endl;
	}

	int nS = (int) tuningWt(SRepIterations);
	if (nS < 0)
		n = MAX_OPTIMIZER_ITERATIONS;
	if (nS != 0)
		binaryPabloSRepStage((int) tuningWt(SRepIterations));

	// ----------------------------------------------------------------


	//
	// SPOKE STAGE: reset reference model and optimize
	if (globalVerbosity > -1)
		cout << "---- SPOKE STAGE" << endl;

	binaryPabloNextOptimizationStage((int) tuningWt(BpSpokeSurfaceLevel));
	if (tuningWt(BpSpokeResetModelBefore) == 1.0) {
		// Use fig/atom stage's model for atomNeighbor and modelMatch distances

		// Replace the reference object with the transformed object.
		// Because the "loaded object" is never supposed to differ from
		// the model in the input file, this operation is bad practice.
		// However, it will suffice until a better solution is determined.
		M3DObject * lobj = object->loadedObject();
		*lobj = *object;
	}
	if (globalVerbosity >= 0) {
		if (tuningWt(BpSpokeResetModelBefore) == 1.0)
			cout << "Reset reference model" << endl;
		else
			cout << "Did NOT reset reference model" << endl;
	}

	n = (int) tuningWt(SpokeIterations);
	if (n < 0)
		n = MAX_OPTIMIZER_ITERATIONS;
	if (n != 0)
		binaryPabloSpokeStage((int) tuningWt(SpokeIterations));

#ifndef OPTIMIZATION_VISUALIZER
	cancelOptimization();	// Delete optimizers
#endif	/* OPTIMIZATION_VISUALIZER */

	//
	// Finally shrink(erode) the fitted model if the image was dilated:
	//
	if (tuningWt(BpDilationFactor) != 0)
		getObjectPtr()->erode(dilationFactorInModelUnits);


	// SAVE THE MODEL
	if (globalVerbosity > -1)
		cout << "---- SAVE" << endl;

	filename = script.getStringValue("OutModel", NULL);
	if (filename && *filename) {
#ifdef INTEGRATE_WHEN_OBJ_HISTORY_STRING_DONE
		// Append filename to figure/model name, so
		// users can distinguish models in Pablo's vis window
		if ((int) tuningWt(BP_TUNE_MRENAME) == 1)
			renameObjectAndFigure(object, filename);
#endif

		if (globalVerbosity > -1)
			cout << "Writing model to " << filename << endl;
		write(filename);	// No status returned, so assume it worked

#ifdef DEBUG_DIBYENDU_WRITE_DIST

		// dibyendu
		// writing out the distances at the spoke ends of all atoms

		if( filename != NULL ) {

			spokeDistFileName = string( filename ) + ".dist.txt" ;		

			this->writeDistancesAtSpokeEnds( (const char *) spokeDistFileName.c_str() ) ;
		}

#endif // DEBUG_DIBYENDU_WRITE_DIST


		if (tuningWt(SaveDistanceVectors) != 0) {
			// Save distance vector list for diagnostic display
			char * dvec_fname = new char[strlen(filename) + strlen(".list") + 1];
			sprintf(dvec_fname,"%s.list",filename);

			if (match != NULL) {
				DistanceVectorList * dvlist = match->getDistanceVectorList();
				if (globalVerbosity > -1)
					cout << "Writing distance vectors to " << dvec_fname << endl;
				std::ofstream outf(dvec_fname, ios::out|ios::trunc);
				std::copy(dvlist->plist.begin(), dvlist->plist.end(),
					ostream_iterator<DistanceVector>(outf));
				outf.close();
			}
			delete [] dvec_fname;
		}
	}
	else
		if (globalVerbosity > -1)
			cout << "OutModel filename undefined.  Model not saved." << endl;

	// Compute and save the tiles of the final model
	filename = script.getStringValue("OutTile", NULL);
	if (filename && *filename)
		exportToBYUTiles(filename, (tuningWt(TileQuads) == 1.0),
		(tuningWt(TileUnit) == 0.0),
		(int) tuningWt(TileSurfaceLevel));
	else
		if (globalVerbosity > -1)
			cout << "OutTile filename undefined.  Tiles not saved." << endl;

		// Chong's modification:
	//
	// The job of the following piece of code is to print out the sum of the objective function values of all the 
	//  atoms on the figure, in order to decide the overall performance of the fitting result.
	cout << "Printing the overall fitting information:" << endl;

	// assume that there is only one figure, get the reference of that figure
	M3DObject * lobj = object->loadedObject();
	*lobj = *object;
	M3DFigure * fig = lobj->getFigurePtr(0);

	// get the number of atoms on that figure
	int numOfAtoms = fig->getPrimitiveCount();
	double penalty = 0.0;

	// for every atom, get the value of the objective function, print it out....
	double sumPenalty = 0.0; 
	for (int i = 0; i < numOfAtoms; i++) {
		Components p = match->computeAtomMatchAndPenaltiesAsComponents(lobj , 0, i, NULL, penalty);
		cout << p.rough << endl;
		sumPenalty += p.rough; 
	//	sumPenalty += p.smooth; <---- do not take smooth value into account
	}
	cout << "Sum of the objective function value: " << sumPenalty << endl;

	// if it runs a figure stage: <---- should output the result for the last iteration of the figure stage, 
	//  since the eigenmodes values are lost in the s-rep stage
	if (tuningWt(FigureIterations) > 0) {
		cout << "Figure stage optimization was performed." << endl;	
	}
	else {
		cout << "No figure stage." << endl;	
	}


	return true;
}

#else	/* BINARY */

/*	Process a script in greyscale Pablo.

The options for loading a model and an image and to optimize must have been
specified.  The registry that is passed should contain valid model and image
pathnames.  At present only a single model and a single image can be loaded.
Also, the image cannot be a stacked image.

The text output of the various operations here is designed to match what
is seen when running interactively.
*/
bool P3DControl::runPablo(Registry & script)
{
	const char * filename;
	const char * outFilename;
	const char * xformFilename;
	const char * pgaFilename;

	// Check for conflicts in the script before running any stages
	if (script.getDoubleValue("FigureGeometryWeight", -1.0) > 0.0 &&
		script.getDoubleValue("FigureMahalanobisMatch", -1.0) > 0.0)
	{
		cout << "Mahalanobis and D/R penalties cannot both be used" << endl;
		return false;
	}

	// Optionally, load a PGA mean model File (command-line option -ipr)
	filename = script.getStringValue("PgaResidueModel", NULL);
	if (filename != NULL) {
		pgaModelFileName = new char[strlen(filename) + 1];
		strcpy(pgaModelFileName, filename);
		if (globalVerbosity > -1)
			cout << "Using PGA mean model " << pgaModelFileName << '\n';
	}

	// Get the replacement PGA model file (command-line option -ip)
	pgaFilename = script.getStringValue("PgaModel", NULL);
	if (pgaFilename != NULL) {
		if (globalVerbosity > -1)
			cout << "PGA statistics will be taken from model "
			<< pgaFilename << "\n";
	}

	// Load the model.  Optionally, in ConStruct, the model may
	// already be loaded.  If pgaModelFileName contains a path, the
	// PgaResidueModel is also loaded.  If the PgaModel is specified,
	// its statistics will replace those in the model.
	filename = script.getStringValue("Model", NULL);
	xformFilename = script.getStringValue("Transform", NULL);
	// Scripts cannot map to world
	if (! read(filename, false, xformFilename, pgaFilename)) {
		if (filename == NULL) {
#ifdef AE2_BUILD
			if (0 != (int) script.getDoubleValue("LoadedModel", 0.0))
				cout << "Using loaded model" << endl;
			else
#endif
			{
				if (globalVerbosity > -1)
					cout << "Error: No model file specified\n";
				return false;
			}
		}
		else {
			if (globalVerbosity > -1) {
				cout << "Error: Unable to load " << filename;
				if (pgaFilename != NULL)
					cout << " or " << pgaFilename;
				cout << '\n';
			}
			return false;
		}
	}
	else
		if (globalVerbosity > -1)
			cout << "Loaded model " << filename << '\n';

	// Load the image
	filename = script.getStringValue("Image", NULL);
	if (filename) {
		if (! loadImage(filename, false)) {		// Scripts cannot map to world
			if (globalVerbosity > -1) {
				if (filename)
					cout << "Error: Unable to load image " << filename << '\n';
				else
					cout << "Error: An image file is required for optimization\n";
			}
			return false;
		}
		else {
			if (globalVerbosity > -1) {
				cout << "Loaded";
				//		if (stacked)
				//			cout << " stacked";
				cout << " image " << filename << '\n';
			}
		}
	}
	else {
#ifdef AE2_BUILD
		filename = script.getStringValue("LoadedImage", NULL);
		if (filename)
			cout << "Using image " << filename << '\n';
		else
#endif
			cout << "Cannot optimize without an input image\n";
	}

	if (script.getBooleanValue("CmdHistTraining", false)) {
		doCommandHistogramTraining(script);
		return true;
	}

	if (script.getBooleanValue("CmdMultiPatchTraining", false)) {
		doCommandMultiPatchTraining(script);
		return true;
	}

	if (script.getBooleanValue("CmdDQFTraining", false)) {
		doCommandDistanceQuantileFunctions(script);
		return true;
	}

	if (script.getBooleanValue("CmdDQFImageMaker", false)) {
		doCommandDistanceQuantileFunctionImageMaker(script);
		return true;
	}

	if (script.getBooleanValue("CmdDQFModelMaker", false)) {
		doCommandDistanceQuantileFunctionModelMaker(script);
		return true;
	}


	// Make sure an output model file name is available
	// NOT EVERYTHING requires an output model, like histogram training
	outFilename = script.getStringValue("OutModel", NULL);
	if (outFilename == NULL) {
		cout << "Optimizing via script requires an output model filename\n";
		return false;
	}

	long rns = script.getIntValue("OptimizerRandomSeed", -1);
	if (rns >= 0) {
		// Seed the random number generator before beginning work
		srand(rns);
		if (globalVerbosity > -1)
			cout << "Set random number seed to " << rns << '\n';
	}

	bool windowIntens = script.getBooleanValue("WindowIntensities", true);
	if (globalVerbosity > 0 && ! windowIntens)
		cout << "Intensities will not be windowed in optimization\n";
	setOptimizerIntensityWindowing(windowIntens);	// Must precede call of initializeMatch()

	if (! initializeOptimizer(script))		// Calls initializeMatch()
		return false;

	// Match has already been initialized, so it is safe to pass it the
	// contour filename
	filename = script.getStringValue("Contour", NULL);
	if (filename) {
		if (! match->initContourMatch(filename))
			cout <<
			"Warning: Contour penalty not specified; ignoring contour points file"
			<< endl;
		else
			cout << "Using contour file " << filename << endl;
	}
	else {
		if (tuningWt(FigureContourMatch) != 0.0) {
			cout << "Error: A contour file must be specified to use contour matching"
				<< endl;
			return false;
		}
	}

	currentOptimizerStage = OPTIMIZER_SETUP_STAGE;
	while (currentOptimizerStage != OPTIMIZER_END_STAGE) {
		OptimizerStage stage = nextOptimizerStage();
		switch (stage) {
case OPTIMIZER_ENSEMBLE_STAGE:
	if (globalVerbosity > -1)
		cout << "---- Ensemble Stage ----" << endl;
	modelStage(true, script);
	break;

case OPTIMIZER_OBJECT_STAGE:
	if (globalVerbosity > -1)
		cout << "---- Object Stage ----" << endl;
	modelStage(false, script);
	break;

case OPTIMIZER_FIGURE_STAGE:
	if (globalVerbosity > -1)
		cout << "---- Figure Stage ----" << endl;
	figureStage(false, script);
	break;

case OPTIMIZER_MAIN_FIGURE_STAGE:
	if (globalVerbosity > -1)
		cout << "---- Main Figure Stage ----" << endl;
	figureStage(true, script);
	break;

case OPTIMIZER_SUBFIGURE_STAGE:
	if (globalVerbosity > -1)
		cout << "---- Subfigure Stage ----" << endl;
	subfigureStage(script);
	break;

case OPTIMIZER_PRIMITIVE_STAGE:
	if (globalVerbosity > -1)
		cout << "---- Atom Stage ----" << endl;
	atomStage(script);
	break;

case OPTIMIZER_SPOKE_STAGE:
	if (globalVerbosity > -1)
		cout << "---- Spoke Stage ----" << endl;
	spokeStage(script);
	break;

case OPTIMIZER_BOUNDARY_DISPLACEMENT_STAGE:
	if (globalVerbosity > -1)
		cout << "---- Boundary Displacement Stage ----" << endl;
	voxelStage(script);
	break;

case OPTIMIZER_END_STAGE:
	if (updateGuiStage != NULL)
		updateGuiStage(OPTIMIZER_END_STAGE);
	break;
		}
		if (globalVerbosity > -1) {
			if (stage == OPTIMIZER_END_STAGE)
				cout << "Optimization completed\n";
			else
				cout << "Stage completed\n";
		}
	}

	bool ret = true;
	if (! write(outFilename)) {
		cout << "Error: no output model was written\n";
		ret = false;
	}
	else
		cout << "Saved model " << outFilename << endl;

	// Compute and save the tiles of the final model
	filename = script.getStringValue("OutTile", NULL);
	if (filename && *filename)
		exportToBYUTiles(filename, (tuningWt(TileQuads) == 1.0),
		(tuningWt(TileUnit) == 0.0),
		(int) tuningWt(TileSurfaceLevel),
		(unsigned long) tuningWt(TileFiguresMask),
		(bool) (tuningWt(TileSeparateFiles) != 0.0));

#ifdef CLUELESS
	// save the final similarity transform
	filename = script.getStringValue("OutSimTransModel", NULL);
	if (filename && *filename)
	{
		// I have no clue as to what to write here! -GT
		//  one cannot write a simTrans from the registration optimizer 
		//  that does not export the simTrans to P3DControl...
		// Use the simpler binary pablo version if you need access to
		//  the simTrans    
	}
#endif

	return ret;
}

bool P3DControl::initializeOptimizer(Registry & script)
{
	const char * filename;
	bool ret;

	int surfaceLevel = script.getIntValue("SurfaceLevel", 2);	// AGG: This should be done per stage

	if (globalVerbosity > 0)
		cout << "Surface level = " << surfaceLevel << endl;
	if (object == NULL) {
		// This should not happen
		cout << "No object!\n";
		return false;
	}
	if (! object->verifyInBounds())
		cout << "Warning: the model is larger than or outside the unit cube" << endl;

	int method = script.getIntValue("OptimizerMatchMethod", -1);

	// Stough's DQF stuff
	const char * dqfConfigFilename = script.getStringValue("DqfConfig", NULL);

	Registry dqfReg;
	bool didDQF = false;

	if(dqfConfigFilename != NULL) {
		bool regError = false;		// Error reading registry; set in 'catch' clause

		// Load the script (a registry)
		try {
			dqfReg.readFromFile(dqfConfigFilename);
		}
		catch (RException excp) {
			excp.print(std::cout);
			cout << "ERR: script file " << dqfConfigFilename << " not found (failed to load).\n";
			regError = true;	// CAN I JUST return HERE?
		}
		if (regError)
			return false;

		didDQF = true;
	}

	if (didDQF)
	{
		dqf = new DQFMatch(dqfReg);

		int numFigs = object->getFigureCount();
		char key[256];
		const char ** dqfStatFilenames = new const char * [numFigs];
		const char ** dqfImageFilenames = new const char * [numFigs];
		for (int f = 0; f < numFigs; f ++)
		{
			sprintf(key, "Dqf%d", f);
			dqfStatFilenames[f]      = script.getStringValue(key, NULL);

			//I also need the images for the dqf.
			sprintf(key, "DqfImage%d", f);
			dqfImageFilenames[f]       = script.getStringValue(key, NULL);
		}

		dqf->initializeDQFMatch(getImagePtr(), object, dqfImageFilenames, dqfStatFilenames);

	}

	ret = false;
	if (globalVerbosity > 0)
		cout << "Match method is ";
	switch (method) {
case -1:	// No match method; fall through for default method
case 0:		// Gaussian Derivative Match
	if (globalVerbosity > 0)
		cout << "light-to-dark Gaussian derivative\n";
	ret = initializeMatch(GAUSSIAN_DERIVATIVE_MATCH, surfaceLevel);
	break;
case 1:		// Negative Gaussian Derivative Match
	if (globalVerbosity > 0)
		cout << "dark-to-light Gaussian derivative\n";
	ret = initializeMatch(NEG_GAUSSIAN_DERIVATIVE_MATCH, surfaceLevel);
	break;
case 6:		// Simple Mask Match
	if (globalVerbosity > 0)
		cout << "simple mask\n";
	filename = script.getStringValue("SimpleMaskFile", NULL);
	if (filename == NULL) {
		cout << "Error: No simple mask file was specified\n";
		return false;
	}
	if (globalVerbosity > 0)
		cout << "Simple mask file = " << filename << endl;
	ret = initializeMatch(SIMPLE_MASK_MATCH, surfaceLevel, NULL, filename);
	break;
case 7:		// Binary Match
	cout << "Binary match cannot be used in a script: run Binary Pablo instead\n";
	return false;
case 8:		// Histogram Match
	if (globalVerbosity > 0)
		cout << "histogram\n";
	// Load the SPE file
	filename = script.getStringValue("Hist", NULL);
	if (! filename)	{	
		if (globalVerbosity > -1) {
			cout << "Error: No surface patch ensemble histogram supplied\n";
			return false;
		}
	}
	else {
		loadSPEObject(filename, script);	// Deletes any old spe object

		if (globalVerbosity > -1)
			cout << "Surface patch ensemble file = " << filename << ".\n";

		// override spe's mean vector with a vector from another SPE file
		const char *speMeanFilename = script.getStringValue("HistMean", NULL);
		if (speMeanFilename && *speMeanFilename) {
			HistogramMatchType hmt = SurfacePatchEnsemble::intToHistogramMatchType(
				script.getIntValue("HistogramMatchType", 0));
			SurfacePatchEnsemble *speMean = new SurfacePatchEnsemble(hmt, speMeanFilename, NULL);
			spe->setMean(speMean->getMean());
			delete speMean;
		}
	}
	ret = initializeMatch(HISTO_IMAGE_MATCH, surfaceLevel);
	break;
case 9:
	//Do the MultiPatch match.
	if (globalVerbosity > 0)
		cout << "MultiPatch match...\n";

	filename = script.getStringValue("MpConfig", NULL);
	if (! filename)	{	
		if (globalVerbosity > -1) {
			cout << "Error: No MultiPatch file supplied\n";
			return false;
		}
	}
	else
	{
		if (! loadMPObject(filename, script))
		{
			cout << "Error: MultiPatch Training failed, MP Object failed to load\n";
			return false;
		}

		if (globalVerbosity > -1)
			cout << "MultiPatch file = " << filename << ".\n";
	}

	ret = initializeMatch(MULTIPATCH_MATCH, surfaceLevel, (char *) script.getStringValue("MPFDebug", NULL));
	break;
default:
	cout << "Unknown match method specified\n";
	return false;
	}

	if (didDQF)
	{
		ret = initializeMatch(HISTO_IMAGE_MATCH, surfaceLevel);
		return true;
	}


	return ret;
}

const char sp4[] = "    ";
const char sp5[] = "     ";
const char sp6[] = "      ";

void P3DControl::modelStage(bool modelStage, Registry & script)
{
	if (0 != script.getIntValue("SkipModelStage", 0)) {
		if (modelStage)
			cout << "Ensemble";
		else
			cout << "Model";
		cout << " stage skipped" << endl;
		return;
	}

	double wt = script.getDoubleValue("ModelPenaltyWeight", 0.5);
	if (globalVerbosity > 0)
		cout << "Penalty weight = " << wt << '\n';
	setRegistrationPenaltyWeight(wt);

	int surfaceLevel = script.getIntValue("SurfaceLevel", 2);

	bool usePGA = script.getBooleanValue("UsePGA", true);			// AGG: The usePGA variable is being ignored below: probably bug
	if (globalVerbosity > 0 && usePGA)								//		Function initializeRegistrationOptimizer() does use PGA,
		cout << "Using PGA data\n";									//		but it apparently cannot be turned off.

	initializeRegistrationOptimizer();
	if (updateGuiStage != NULL)
		updateGuiStage(currentOptimizerStage);

	const char * sp = sp4;
	for (int i = 1; i < MAX_OPTIMIZER_ITERATIONS; i++) {
		doRegistrationIterations(1);
		if (updateGui != NULL)
			updateGui(i);
		if (globalVerbosity > 0) {
			cout << i << ".  Current penalty: "
				<< getLastRegistrationPenaltyValue() << '\n';;
			if ((i  + 1)% 10 == 0) {
				if ((i + 1) % 100 == 0)
					sp = sp6;
				else
					sp = sp5;
			}
			cout << sp << "Current match: " << -getLastRegistrationMatchValue()
				<< '\n';
		}

		if (registrationOptimizerIsFinished())
			break;
	}
}

bool writeModelSimply(const char * filename, M3DObject * object)
{
	if (!filename || !*filename || !object)
		return false;

	M3DObjectFile objectFile;
	return objectFile.write(filename, *object);
}

void P3DControl::figureStage(bool mainFigureStage, Registry & script)
{
	double wt;
	int n;

	bool usePGA = script.getBooleanValue("UsePGA", true)  && (script.getDoubleValue("FigureMahalanobisMatch",-1.0) >0.0);

	if (usePGA) {
		if (! haveFigurePGA()) {
			cout << "No PGA statistics available; ignoring request to use PGA data" << endl;
			usePGA = false;
		}
		else if (globalVerbosity > 0)
			cout << "Using PGA data\n";
	}

	clearPenaltyWeights();

	// Set the penalty weights
	wt = script.getDoubleValue("FigureMahalanobisMatch", -1.0);
	if (wt <= 0.0) {
		wt = script.getDoubleValue("FigureGeometryWeight", 0.5);
		if (globalVerbosity > 0)
			cout << "Geometry weight = " << wt << '\n';
		setMainFigurePenaltyWeight(GpFigureGeometryWeight, wt);
	}
	else {
		if (globalVerbosity > 0)
			cout << "Mahalanobis weight = " << wt << '\n';
		setMainFigurePenaltyWeight(FigureMahalanobisMatch, wt);
	}

	if (object->getFigureCount() > 1) {
		wt = script.getDoubleValue("FigureInterObjectWeight", 0.5);
		if (globalVerbosity > 0)
			cout << "Inter-object weight = " << wt << '\n';
	}
	else
		wt = 0.0;	// Single figure models don't have constraint penalties
	setMainFigurePenaltyWeight(GpFigureInterObjectWeight, wt);


	if (! initializeMainFigureOptimizer(false, usePGA)) {
		cout << "Figure stage could not be run" << endl;
		return;
	}

	// why is this not a tuning value? --GST
	n = script.getIntValue("FigureIterations", MAX_OPTIMIZER_ITERATIONS - 1);
	if (globalVerbosity > 0) {
		cout << "Number of iterations ";
		if (n == MAX_OPTIMIZER_ITERATIONS - 1)
			cout << "is unlimited\n";
		else
			cout << "= " << n << '\n';
	}
	n++;

	if (updateGuiStage != NULL)
		updateGuiStage(currentOptimizerStage);


	const char * sp = sp4;
	if( n == 0 ) {
		// Even if iterations are 0, run the optimizer just to set the starting position if given.
		doMainFigureIterations(0, true);

		if (globalVerbosity > -1) {
			int i = 0;
			cout << i << ".  Geometric deviation: "
				<< getLastMainFigureGeometricMatchValue() << '\n';
			if (i % 10 == 0) {
				if (i % 100 == 0)
					sp = sp6;
				else
					sp = sp5;
			}
			cout << sp << "Image match: " << (getLastMainFigureObjectiveFunctionValue()
				- getLastMainFigureGeometricMatchValue()) << '\n';
			cout << sp << "Objective function: "
				<< getLastMainFigureObjectiveFunctionValue() << '\n';
		}
	}

	// If this string is defined, it's used both to trigger writing a model
	// after each iteration and to specify the printf-style format string.
	const char *writeModel = tuningStr(FigureModelPerIteration);

	for (int i = 1; i < n; i++) {
		doMainFigureIterations(1, i == 1);
		if (mainFigureStageSkipped())	// Was this figure skipped?
			break;

		if (updateGui != NULL)
			updateGui(i);

		if (globalVerbosity > 0) {
			cout << i << ".  Geometric deviation: "
				<< getLastMainFigureGeometricMatchValue() << '\n';
			if (i % 10 == 0) {
				if (i % 100 == 0)
					sp = sp6;
				else
					sp = sp5;
			}
			cout << sp << "Image match: " << (getLastMainFigureObjectiveFunctionValue()
				- getLastMainFigureGeometricMatchValue()) << '\n';
			cout << sp << "Objective function: "
				<< getLastMainFigureObjectiveFunctionValue() << '\n';
			cout.flush();
		}

		if (writeModel && *writeModel) {
			char modelFileName[1000];
			sprintf(modelFileName, writeModel, i);
			cout << "  Writing model figPerIter to " << modelFileName << endl;
			writeModelSimply(modelFileName, mainFigureOptimizer->getLastBestObject());
		}

		if (mainFigureOptimizerIsFinished())
			break;
	}

	writeModel = tuningStr(FigureModelResult);
	if (writeModel && *writeModel) {
		cout << "  Writing model result to " << writeModel << endl;
		writeModelSimply(writeModel, object);
	}	
}

void P3DControl::atomStage(Registry & script)
{
	M3DFigure * currFigure;
	int i, n, m, figureId, numFigures, len;
	double wt;

	clearPenaltyWeights();

	// Set the penalty weights
	wt = script.getDoubleValue("AtomGeometryWeight", 0.5);
	if (globalVerbosity > 0)
		cout << "Geometry weight = " << wt << '\n';
	setDeformationPenaltyWeight(GpAtomGeometryWeight, wt);

	wt = script.getDoubleValue("AtomInterWeight", 0.5);
	if (globalVerbosity > 0)
		cout << "Inter-atom weight = " << wt << '\n';
	setDeformationPenaltyWeight(GpAtomInterWeight, wt);

	wt = script.getDoubleValue("AtomNumberWeight", 0.5);
	if (globalVerbosity > 0)
		cout << "Neighbor penalty weight = " << wt << '\n';
	setDeformationPenaltyWeight(GpAtomNeighborWeight, wt);

	n = script.getIntValue("AtomIterations", -1);
	if (n < 0)
		n = MAX_OPTIMIZER_ITERATIONS;
	m = script.getIntValue("AtomSubIterations", -1);	// Useful for debugging
	if (m < 0)
		m = MAX_OPTIMIZER_ITERATIONS;
	if (globalVerbosity > 0) {
		cout << "Number of iterations ";
		if (n >= MAX_OPTIMIZER_ITERATIONS)
			cout << "is unlimited\n";
		else
			cout << "= " << n << '\n';
		if (m < MAX_OPTIMIZER_ITERATIONS)
			cout << "Number of atom fits = " << m << '\n';
	}

	initializeDeformationOptimizer();
	if (updateGuiStage != NULL)
		updateGuiStage(currentOptimizerStage);

	// Find the (selected) figure to be optimized
	M3DDeformationOptimizer * deformationOptimizer = getDeformationOptimizer();
	M3DObject * candidateObject = deformationOptimizer->getCandidateObject();
	numFigures = candidateObject->getFigureCount();
	for (figureId = 0; figureId < numFigures; figureId++) {
		currFigure = candidateObject->getFigurePtr(figureId);
		if (currFigure != NULL && currFigure->isAnySelected())
			break;
	}

	const int * skipFigs = tuningList(AtomSkipFigures, len);

	// Optimize the figure
	const char * sp = sp4;
	if (figureId < numFigures) {
		// Check to see if this figure is to be skipped
		for (i = 0; i < len; i++)
			if (skipFigs[i] == figureId) {
				cout << "Figure number " << figureId;
				M3DFigure * figure = object->getFigurePtr(figureId);
				const char * name = figure->getName();
				if (name != NULL)
					cout << " (" << name << ')';
				cout << " skipped" << endl;
				return;
			}

			if (globalVerbosity >= 0) {
				cout << "Figure number " << figureId;
				M3DFigure * figure = object->getFigurePtr(figureId);
				const char * name = figure->getName();
				if (name != NULL)
					cout << " (" << name << ")\n";
				else
					cout << '\n';
			}

			// Atom PGA 
			bool wantPGA = (tuningWt(AtomMahalanobisMatch) != 0.0);
			if (wantPGA) {
				if (atomPgaStats == NULL) {
					cout << "Error: no atom PGA data available; ignoring AtomMahalanobisMatch"
						<< endl;
					wantPGA = false;
				}
			}

			if (wantPGA) {	
				deformationOptimizer->initializeAtomPGAStage(atomPgaStats,figureId);
			}

			// An iteration in the atom stage is a cycle through all figural atoms
			n *= (object->getFigurePtr(figureId))->getPrimitiveCount();
			if (m < n) {	// Use the lesser of the AtomIterations and AtomSubIterations settings
				n = m;
				if (globalVerbosity >= 1)
					cout << "AtomSubIterations supersedes atomIterations\n";
			}

			for (i = 0; i < n; i++) {
				bool optimizerStopRequested = doDeformationIterations(figureId, false, (globalVerbosity <= 0));
				if (updateGui != NULL)
					updateGui(i);		

				cout << i + 1 << ".  Geometric deviation: ";
				if (tuningWt(AtomMahalanobisMatch)>0)
					cout << getLastDeformationGeometricMatchValue()/tuningWt(AtomMahalanobisMatch) << '\n';
				else 
					cout <<  getLastDeformationGeometricMatchValue()<< '\n';

				if ((i  + 1)% 10 == 0) {
					if ((i + 1) % 100 == 0)
						sp = sp6;
					else
						sp = sp5;
				}
				cout << sp << "Image match: " << (getLastDeformationObjectiveFunctionValue(false)
					-getLastDeformationGeometricMatchValue()) << '\n';
				cout << sp << "Objective function: " <<
					getLastDeformationObjectiveFunctionValue(false) << endl;
				//	}

				if (optimizerStopRequested)
					break;

				if (optimizerStopRequested)
					break;
			}

			// Apply the initial transformation to the final atom stage result, which is deformed
			// in the feature space but "evaluat"ed in image space each time in order to use the
			// statistics
			if (atomPgaStats != NULL){
				SimilarityTransform3D * xform = candidateObject->getTransformation();
				if (xform!= NULL){
					Vector3D zero(0.0, 0.0, 0.0);
					object->applySimilarity( *xform, zero );                			
				}
			}
	}

	if (globalVerbosity <= 0)
		cout << endl;	// Needed because the deformation optimizer prints dots as it runs.
}

void P3DControl::spokeStage(Registry & script)
{
	M3DFigure * currFigure;
	int i, n, figureId, numFigures, len;
	double wt;

	clearPenaltyWeights();

	// Set the penalty weights
	wt = script.getDoubleValue("SpokeGeometryWeight", 0.5);
	if (globalVerbosity > 0)
		cout << "Geometry weight = " << wt << '\n';
	setSpokePenaltyWeight(GpSpokeGeometryWeight, wt);

	wt = script.getDoubleValue("SpokeInterWeight", 0.5);
	if (globalVerbosity > 0)
		cout << "Inter-atom weight = " << wt << '\n';
	setSpokePenaltyWeight(GpSpokeInterWeight, wt);

	wt = script.getDoubleValue("SpokeNumberWeight", 0.5);
	if (globalVerbosity > 0)
		cout << "Neighbor penalty weight = " << wt << '\n';
	setSpokePenaltyWeight(GpSpokeNeighborWeight, wt);

	n = int(tuningWt(SpokeIterations));
	if (n < 0)
		n = MAX_OPTIMIZER_ITERATIONS;
	if (globalVerbosity > 0) {
		cout << "Number of iterations ";
		if (n >= MAX_OPTIMIZER_ITERATIONS)
			cout << "is unlimited\n";
		else
			cout << "= " << n << '\n';
	}

	initializeSpokeOptimizer();
	if (updateGuiStage != NULL)
		updateGuiStage(currentOptimizerStage);

	// Find the (selected) figure to be optimized
	M3DSpokeOptimizer * spokeOptimizer = getSpokeOptimizer();
	M3DObject * candidateObject = spokeOptimizer->getCandidateObject();
	numFigures = candidateObject->getFigureCount();
	for (figureId = 0; figureId < numFigures; figureId++) {
		currFigure = candidateObject->getFigurePtr(figureId);
		if (currFigure != NULL && currFigure->isAnySelected())
			break;
	}

	const int * skipFigs = tuningList(SpokeSkipFigures, len);

	// Optimize the figure
	const char * sp = sp4;
	if (figureId < numFigures) {
		// Check to see if this figure is to be skipped
		for (i = 0; i < len; i++)
			if (skipFigs[i] == figureId) {
				cout << "Figure number " << figureId;
				M3DFigure * figure = object->getFigurePtr(figureId);
				const char * name = figure->getName();
				if (name != NULL)
					cout << " (" << name << ')';
				cout << " skipped" << endl;
				return;
			}

			if (globalVerbosity >= 0) {
				cout << "Figure number " << figureId;
				M3DFigure * figure = object->getFigurePtr(figureId);
				const char * name = figure->getName();
				if (name != NULL)
					cout << " (" << name << ")\n";
				else
					cout << '\n';
			}

			// Spoke PGA 
			bool wantPGA = (tuningWt(SpokeMahalanobisMatch) != 0.0);
			if (wantPGA) {
				if (spokePgaStats == NULL) {
					cout << "Error: no spoke PGA data available; ignoring SpokeMahalanobisMatch"
						<< endl;
					wantPGA = false;
				}
			}

			if (wantPGA) {
				M3DObject * meanObj;

				// The meanObj is used for alignment before applying the spoke difference PGA
				meanObj = NULL;
				if (pgaStats != NULL)
					meanObj = pgaStats->getMeanObjPtr();

				if (meanObj != NULL)
					spokeOptimizer->initializeSpokePGAStage(spokePgaStats,meanObj,figureId);
				else {
					cout << "Error: meanObj is not available!\n" <<endl;
					return;
				}
			}

			// An iteration in the spoke stage is a cycle through all figural atoms
			n *= (object->getFigurePtr(figureId))->getPrimitiveCount();

			for (i = 0; i < n; i++) {
				bool optimizerStopRequested = doSpokeIterations(figureId, false, (globalVerbosity <= 0));
				if (updateGui != NULL)
					updateGui(i);
				if (globalVerbosity > 0) {

					cout << i + 1 << ".  Geometric deviation: ";
					if (tuningWt(SpokeMahalanobisMatch)>0)
						cout << -getLastSpokeGeometricMatchValue()/tuningWt(SpokeMahalanobisMatch) << '\n';
					else 
						cout << -getLastSpokeGeometricMatchValue()<< '\n';

					if ((i  + 1)% 10 == 0) {
						if ((i + 1) % 100 == 0)
							sp = sp6;
						else
							sp = sp5;
					}
					cout << sp << "Image match: " << -(getLastSpokeObjectiveFunctionValue(false)
						-getLastSpokeGeometricMatchValue()) << '\n';
					cout << sp << "Objective function: " <<
						-getLastSpokeObjectiveFunctionValue(false) << endl;
				}
				if (optimizerStopRequested)
					break;
			}
	}

	if (globalVerbosity <= 0)
		cout << endl;	// Needed because the spoke optimizer prints dots as it runs.
}


void P3DControl::subfigureStage(Registry & script)
{
	cout << "Not implemented" << endl;

	if (updateGuiStage != NULL)
		updateGuiStage(currentOptimizerStage);
}

void P3DControl::voxelStage(Registry & script)
{
	cout << "Not implemented" << endl;

	if (updateGuiStage != NULL)
		updateGuiStage(currentOptimizerStage);
}

//As below from doCommandHistogramTraining
bool P3DControl::doCommandMultiPatchTraining(Registry & script)
{
	Image3D * image = getImagePtr();

	const char * mpFilename       = script.getStringValue("MpConfig", NULL);

	if (! loadMPObject(mpFilename, script)){
		cout << "MultiPatch Training failed, MP Object failed to load\n";
		return false;
	}


	const char * saveFilename = script.getStringValue("MpOut", NULL);

	if (saveFilename == NULL){
		cout << "MultiPatch Training failed, no save filename found.\n";
		return false;
	}

	mp->initializePatches(object, image);


	int f, numFigures = object->getFigureCount();
	char * key = new char[1024];
	AllImageIO imageFile;
	Image3D *binaryImage = NULL;

	for (f = 0; f < numFigures; f++)
	{
		sprintf(key, "MPBinaryImage%d", f);
		const char * binaryImageFilename = script.getStringValue(key, NULL);

		if (binaryImageFilename != NULL)
		{
			binaryImage = imageFile.read(binaryImageFilename);

			if (binaryImage == NULL) {
				cout << "Binary image " << binaryImageFilename << " FAILED to load.\n";
				return false;
			}
		}

		mp->sampleYourself(object, image, f, binaryImage, true);

		if (binaryImage != NULL)
			delete binaryImage;
		binaryImage = NULL;

	}

	mp->saveMP(saveFilename);

	cout << "MutliPatch file saved " << saveFilename << ".\n";


	return true;
}

bool P3DControl::doCommandHistogramTraining(Registry & script)
{
	// Model and image have been loaded
	Image3D * image = getImagePtr();

	const char * speFilename       = script.getStringValue("Hist", NULL);

	if (! loadSPEObject(speFilename, script)) {	// Deletes any old spe object
		cout << "Histogram Training failed, SPE Object failed to load\n";
		return false;
	}

	const char * saveFilename = script.getStringValue("HistOut", NULL);

	if (saveFilename == NULL){
		cout << "Histogram Training failed, no save filename found.\n";
		return false;
	}

	int f, numFigures = object->getFigureCount();
	char key[256];

	for (f = 0; f < numFigures; f++){
		sprintf(key, "HistPatchImage%d", f);
		const char * patchImageFilename = script.getStringValue(key, NULL);

		sprintf(key, "HistBinaryImage%d", f);
		const char * binaryImageFilename = script.getStringValue(key, NULL);

		spe->update(image, object, f, -1, binaryImageFilename, patchImageFilename);
	}

	spe->save(saveFilename);
	cout << "Surface Patch Ensemble File Saved: " << saveFilename << ".\n";

	return true;
}

const char * P3DControl::calcDQFFilename(char * imageFilename, int figID)
{
	char * dqfFilename;		// the derived filename
	int len = strlen(imageFilename);

	// input file must end in ".raw3" or ".raw3z"; derived file will have raw3
	char * suffix = NULL;
	if ((strstr(imageFilename, ".raw3" ) == (imageFilename + len - strlen(".raw3" ))))
		suffix = ".raw3";
	if ((strstr(imageFilename, ".raw3z" ) == (imageFilename + len - strlen(".raw3z" ))))
		suffix = ".raw3z";
	if ((strstr(imageFilename, ".gipl" ) == (imageFilename + len - strlen(".gipl" ))))
		suffix = ".gipl";
	if ((strstr(imageFilename, ".mhd" ) == (imageFilename + len - strlen(".mhd" ))))
		suffix = ".mhd";
	if ((strstr(imageFilename, ".mha" ) == (imageFilename + len - strlen(".mha" ))))
		suffix = ".mha";
	if (!suffix)
		return NULL;

	int suffixLen = strlen(suffix);
	dqfFilename = new char[len+suffixLen+1+7];		// insert '-ddm' -0xABCD
	strcpy(dqfFilename, imageFilename);
	strcpy(dqfFilename + len - suffixLen, "-dqf");

	sprintf(dqfFilename + len - suffixLen + 4, "-%i", figID);

	strcat(dqfFilename, suffix);

	//so asdf.raw3 should be asdf-dqf-1.raw3 or whatever figID is.
	return dqfFilename;
}

bool P3DControl::doCommandDistanceQuantileFunctions(Registry & script)
{

	Image3D * image = getImagePtr();

	const char * filename = script.getStringValue("Image", NULL);
	char * imageFilename = strdup(filename);

	const char * dqfFilename = calcDQFFilename(imageFilename,0);

	const char * dqfConfigFilename = script.getStringValue("DqfConfig", NULL);

	Registry dqfReg;

	if(dqfConfigFilename != NULL) {
		bool regError = false;		// Error reading registry; set in 'catch' clause

		// Load the script (a registry)
		try {
			dqfReg.readFromFile(dqfConfigFilename);
		}
		catch (RException excp) {
			excp.print(std::cout);
			cout << "ERR: script file " << dqfConfigFilename << " not found (failed to load).\n";
			regError = true;	// CAN I JUST return HERE?
		}
		if (regError)
			return false;
	}

	DQFTraining *myTrainer = new DQFTraining(dqfReg);

	//if (myTrainer->createDQFImage(image, dqfFilename, object, 0))
	//{
	//	cout << "DQF Image maker succeeded." << endl;
	//}


	if (myTrainer->doModelDQFTraining(image, object, 0, NULL, dqfFilename))
	{ 
		cout << "made dqf model file.\n";
	}

	return true;
}

bool P3DControl::doCommandDistanceQuantileFunctionImageMaker(Registry & script)
{

	Image3D * image = getImagePtr();

	const char * filename = script.getStringValue("Image", NULL);
	char * imageFilename = strdup(filename);

	int figID = (int) tuningWt(GP_TUNE_FIGURE_ID);

	const char * dqfFilename = calcDQFFilename(imageFilename, figID);

	const char * dqfConfigFilename = script.getStringValue("DqfConfig", NULL);

	Registry dqfReg;

	if(dqfConfigFilename != NULL) {
		bool regError = false;		// Error reading registry; set in 'catch' clause

		// Load the script (a registry)
		try {
			dqfReg.readFromFile(dqfConfigFilename);
		}
		catch (RException excp) {
			excp.print(std::cout);
			cout << "ERR: script file " << dqfConfigFilename << " not found (failed to load).\n";
			regError = true;	// CAN I JUST return HERE?
		}
		if (regError)
			return false;
	}

	DQFTraining *myTrainer = new DQFTraining(dqfReg);

	if (myTrainer->createDQFImage(image, dqfFilename, object, figID))
	{
		cout << "DQF Image maker succeeded." << endl;
	}

	return true;
}

bool P3DControl::doCommandDistanceQuantileFunctionModelMaker(Registry & script)
{

	Image3D * image = getImagePtr();

	const char * filename = script.getStringValue("Image", NULL);
	char * imageFilename = strdup(filename);

	int figID = (int) tuningWt(GP_TUNE_FIGURE_ID);

	//Rather than compute the dqf filename from the image filename, since it may have
	//been moved, I'm going to let the -dqfImage1, -dqfImage2 etc refer to the dqf images in the 
	//modelMaker instance.  In the match instance, dqf's refer to dqf statistics files.  
	//I may also make it that in the imageMaker instance, they may refer to the output dqf
	//image names...  
	//Basically, I will send DQFTraining a const char ** with the filenames I pick up.

	const char ** dqfImageFilenames = new const char * [object->getFigureCount()];
	char key[256];
	//const char * dqfImageFilename = calcDQFFilename(imageFilename, figID);
	for (int i = 0; i < object->getFigureCount(); i ++)
	{
		sprintf(key, "DqfImage%d", i);
		dqfImageFilenames[i] = script.getStringValue(key, NULL);
	}
	//Actually, I forgot, I only need the one dqf image file now, since training will 
	//happen per figure. It's in the match where I need to fix stuff...  Crap, in dqfmatch
	//I'm calling the calcDQF stuff, that's stupid...

	const char * dqfConfigFilename = script.getStringValue("DqfConfig", NULL);

	Registry dqfReg;

	if(dqfConfigFilename != NULL) {
		bool regError = false;		// Error reading registry; set in 'catch' clause

		// Load the script (a registry)
		try {
			dqfReg.readFromFile(dqfConfigFilename);
		}
		catch (RException excp) {
			excp.print(std::cout);
			cout << "ERR: script file " << dqfConfigFilename << " not found (failed to load).\n";
			regError = true;	// CAN I JUST return HERE?
		}
		if (regError)
			return false;
	}

	const char * dqfOutFilename = script.getStringValue("DqfOut", NULL);

	DQFTraining *myTrainer = new DQFTraining(dqfReg);

	//Here's where we could do multiple figure training, with a for loop.  But I'll keep
	//it single figure for now. (one could conceiveably want a different config file for each
	//figure...
	if (myTrainer->doModelDQFTraining(image, object, figID, dqfOutFilename, dqfImageFilenames[figID]))
	{ 
		cout << "made dqf model file.\n";
	}

	return true;
}

bool P3DControl::loadMPObject(const char *mpFilename, Registry &script)
{
	if (mp != NULL){
		delete mp;
		mp = NULL;
	}

	Registry mpReg;
	//At this point there is only one match type (as in mahal, euclidean, etc.) only gaussian.
	//and as for the patch type, that is is the mpconfig file.

	if(mpFilename == NULL) {
		cout << "ERR: MultiPatch (MP) script file not supplied - required for MP training.\n";
		return false;
	}

	mp = new MultiPatch(mpFilename);


	return true;
}

bool P3DControl::loadSPEObject(const char *speFilename, Registry &script)
{
	if (spe != NULL){
		delete spe;
		spe = NULL;
	}

	Registry speReg;
	HistogramMatchType hmt = SurfacePatchEnsemble::intToHistogramMatchType(script.getIntValue("HistogramMatchType", 0));

	if(speFilename == NULL) {
		cout << "ERR: SPE script file non supplied - required for SPE training.\n";
		return false;
	}

	const char * QFDataFilename = script.getStringValue("QFData", NULL);

	spe = new SurfacePatchEnsemble(hmt, speFilename, QFDataFilename);

	return true;
}

#endif	/* BINARY */

// ----------------------------- End of Scripting Code --------------------------------

// Function called as the final step when loading an image for Pablo's use.
// The main thing does is install the image into the ImagePlanes object.  If
// external is true, the image will not become owned by ImagePlanes.  This
// capability is provided for programs other than Pablo.
bool P3DControl::finishImageLoading(Image3D * image, bool stacked, bool mapToWorld,
									bool external)
{
	if (image == NULL)
		return false;

	if (mapToWorld) {
		if (object != NULL && worldWasApplied) {
			Image3D * im = imagePlanes.getImagePtr();
			object->unApplyWorld(im);
		}
	}

	image->setIsImageStacked(stacked);
	imagePlanes.setImagePtr(image, 512, external);
	imagePlanesChanged = true;

	if (match != NULL)
		match->setTargetImage(image);

	if (mapToWorld) {
		applyWorld(image);

		if (referenceObject != NULL) {
			delete referenceObject;
			referenceObject = NULL;
		}
	}

	blurred = false;
	return true;
}

void P3DControl::discardImage(bool external)
{
	imagePlanes.setImagePtr(NULL, 0, external);
	imagePlanesChanged = true;

	if (match != NULL)
		match->setTargetImage(NULL);

	if (imageFileName != NULL) {
		delete [] imageFileName;
		imageFileName = NULL;
	}
}

// With only the filename, i.e. with the second argument NULL, the main
// image (stored in ImagePlanes) is saved.  For special cases, a different
// image may be provided by the second argument.
bool P3DControl::saveImage(const char * filename, Image3D * image,
						   GreyValue lowIntens, GreyValue highIntens)      // AGG: These args should be ints
{
	AllImageIO file;

	if (image == NULL)
		image = imagePlanes.getImagePtr();

	if (image == NULL)
		return false;

	if (filename == NULL || filename[0] == '\0')
		return false;

	if (lowIntens > highIntens)
		image->range(lowIntens, highIntens);
	return file.write(filename, *image, lowIntens, highIntens);
}

bool P3DControl::saveDQF(const char * filename, DQFImage * dqf) const
{
	AllImageIO file;

	if (dqf == NULL)
		return false;

	if (filename == NULL || filename[0] == '\0')
		return false;

	// Range is meaningless for DQF images, and scaling and shifting
	// when reading must be avoided.
	return file.write(filename, *dqf, MIN_GREY_VALUE, MAX_GREY_VALUE,
		AllImageIO::Raw3);
}

// Lowest-level function to read and return an image.
DQFImage * P3DControl::readDQF(const char * filename) const
{
	AllImageIO file;

	if (filename == NULL)
		return NULL;

	return (DQFImage *) file.read(filename, false, false);
}

// Lowest-level function to read and return an image.
Image3D * P3DControl::readImage(const char * filename, bool stacked,
								bool headerOnly) const
{
	AllImageIO file;

	if (filename == NULL)
		return NULL;

	return file.read(filename, stacked, headerOnly);
}

bool P3DControl::loadImage(const char * filename, bool stacked,
						   bool mapToWorld, bool headerOnly, bool external)
{
	int len;
	Image3D * image;

	if (filename == NULL)
		return false;

	len = strlen(filename);
	if (0 == len)
		return false;

	image = readImage(filename, stacked, headerOnly);
	if (image == NULL)
		return false;

	if (imageFileName != NULL) 
	{
		delete [] imageFileName;
		imageFileName = NULL;
	}
#ifdef _WIN32
	if (filename[1] != ':' || (filename[2] != '\\' && filename[2] != '/'))
#else
	if (filename[0] != '/')
#endif
	{
		// Filename is a relative path
#ifdef _WIN32
#define MAX_PATH_LEN 512
#else
#define MAX_PATH_LEN PATH_MAX
#endif
		char cwd[MAX_PATH_LEN + 1];
#ifdef _WIN32
		if (_getcwd(cwd, MAX_PATH_LEN) == NULL)
#else
		if (getcwd(cwd, MAX_PATH_LEN) == NULL)
#endif
			strcpy(cwd, "...");
		imageFileName = new char[strlen(cwd) + len + 200];
		//imageFileName = (char*)malloc(sizeof(char) * strlen(cwd) + len + 200);
		strcpy(imageFileName, cwd);
		strcat(imageFileName, "/");
		strcat(imageFileName, filename);
	}
	else {
		// Filename is an absolute path
		imageFileName = new char[len + 1];
		//imageFileName = (char*)malloc(sizeof(char) * len + 1);
		strcpy(imageFileName, filename);
	}
#ifdef _WIN32
	reverseBackSlashes(imageFileName);
#endif

	return finishImageLoading(image, stacked, mapToWorld, external);
}

bool P3DControl::getIsImageStacked()
{
	Image3D * image;

	image = imagePlanes.getImagePtr();
	if (image == NULL)
		return false;
	return image->getIsImageStacked();
}

void P3DControl::setInputImageScaling(bool scaleInput)
{
	AllImageIO::setImageScaling(scaleInput);
}

#ifdef BINARY

void P3DControl::setDistMapPtr(ImageDistanceMap * im)
{
	dist_map = im ;

	if (match == NULL)
		return;

	// Distance maps can be created either in P3DControl or Match.
	// Ownership of the distance map is passed over to Match,
	// because Match deletes the distance map; so, this class
	// can only delete it if match never exists.

	match->setDistMap(dist_map);

	dist_map = NULL;

};

// Access to the distance map must be through this function, because
// is variously owned by classes P3DControl and Match.
ImageDistanceMap * P3DControl::getDistMapPtr()
{
	if (dist_map)
		return dist_map;
	if (match)
		return match->getImageDistMapPointer();
	return NULL;
}

bool P3DControl::loadDistMap(const char * filename)
{
	AllImageIO imageFile;
	Image3D * distMapImage;	

	if (! filename || ! *filename)
		return false;

	cout << filename << endl;
	distMapImage = imageFile.read(filename);

	if (distMapImage == NULL)
		return false;

	ImageDistanceMap * dm = new ImageDistanceMap;
	dm->fromImage3D(distMapImage);
	setDistMapPtr(dm);
	return true;

# define WRITE_DISTANCE_MAPS

#ifdef WRITE_DISTANCE_MAPS

	char* outputDistanceMapDataFileName = calcOutputDMFilename( (char *) filename ) ;	

	GreyValue minI, maxI ;

	distMapImage->range( minI, maxI ) ;

	imageFile.write( outputDistanceMapDataFileName, *distMapImage, minI, maxI ) ; 

#endif

}

bool P3DControl::loadAntialiasedDistMap(const char * filename )
{
	int len = strlen(filename);
	const char* ddmFileExtension = ".mhd";
	int suffixLen = strlen(".mhd");

	char* filename_gradX = NULL ;
	char* filename_gradY = NULL ;
	char* filename_gradZ = NULL ;

	// dibyendu

	// test to see if the gradient map loading works for the pablo-generated distance map file

	bool USE_PABLO_DMAP = 1 ;

	if( USE_PABLO_DMAP ) {

		filename_gradX = new char[len] ;
		strcpy( filename_gradX, filename ) ;

		filename_gradY = new char[len] ;
		strcpy( filename_gradY, filename ) ;

		filename_gradZ = new char[len] ;
		strcpy( filename_gradZ, filename ) ;

	}
	else {

		filename_gradX = new char[len+suffixLen+6];
		strcpy(filename_gradX, filename);
		strcpy(filename_gradX + len - suffixLen, "_gradX");
		filename_gradX = strcat(filename_gradX, ddmFileExtension);

		filename_gradY = new char[len+suffixLen+6];
		strcpy(filename_gradY, filename);
		strcpy(filename_gradY + len - suffixLen, "_gradY");
		filename_gradY = strcat(filename_gradY, ddmFileExtension);

		filename_gradZ = new char[len+suffixLen+6];
		strcpy(filename_gradZ, filename);
		strcpy(filename_gradZ + len - suffixLen, "_gradZ");
		filename_gradZ = strcat(filename_gradZ, ddmFileExtension);

	}

//char * filename_gradX = calcAntialisedGradFilename(filename,'X');
//char * filename_gradY = calcAntialisedGradFilename(filename,'Y');
//char * filename_gradZ = calcAntialisedGradFilename(filename,'Z');

	AllImageIO imageIOHandler ;			//xiaojie ??

	Image3D * antialiasedDistMapImage ;

	Image3D * antialiasedDistGradXImage ;
	Image3D * antialiasedDistGradYImage ;
	Image3D * antialiasedDistGradZImage ;

	if (! filename || ! *filename||! filename_gradX || ! *filename_gradX||! filename_gradY 
		|| ! *filename_gradY||! *filename_gradZ||! filename_gradZ)
		return false;


	antialiasedDistMapImage = imageIOHandler.read(filename);

	if( antialiasedDistMapImage == NULL )
		return( false ) ;

	// -------------------------------------------------------------------

    antialiasedDistGradXImage = imageIOHandler.read(filename_gradX);

    antialiasedDistGradYImage = imageIOHandler.read(filename_gradY);

    antialiasedDistGradZImage = imageIOHandler.read(filename_gradZ);

	if ( (antialiasedDistMapImage == NULL) || (antialiasedDistGradXImage == NULL) ||
		 (antialiasedDistGradYImage == NULL) || (antialiasedDistGradZImage == NULL) ) {
		return false;
	}

	// -------------------------------------------------------------------

	ImageDistanceMap * dm = new ImageDistanceMap ;

	dm->fromImage3D(antialiasedDistMapImage);

// -------------------------------------------------------------------

	dm->createGradDistMapFromGradXYZ(antialiasedDistGradXImage,antialiasedDistGradYImage,antialiasedDistGradZImage);

	// dm->createGradDistMapFromGradXYZ(antialiasedDistMapImage,antialiasedDistMapImage,antialiasedDistMapImage);//for test
// -------------------------------------------------------------------

	setDistMapPtr(dm);
		
	//if( filename_gradX ) 
	//		delete[] filename_gradX ;
	//if( filename_gradY ) 
	//		delete[] filename_gradY ;
	//if( filename_gradZ ) 
	//		delete[] filename_gradZ ;

	return true;
}
#endif

bool P3DControl::loadTileSet(const char * filename, bool split, bool useWorldCoords)
{
	if (filename == NULL || 0 == strlen(filename))
		return false;

	Image3D * image;
	if (useWorldCoords) {
		image = imagePlanes.getImagePtr();
		if (image == NULL) {
			cout << "An image must exist to load BYU tiles in world coordinates"
				<< endl;
			beep();
			image = NULL;
		}
	}
	else {
		image = NULL;
	}

	if (tileSet != NULL) {
		if (tileSet->loadFile(filename, split, image) == false)
			return false;
	}
	else {
		tileSet = new TileSet(filename, split, image);
		if (tileSet->tileCount() == 0)
			return false;
	}
	tileSet->calcCOG(tileSetMobileFigure);
	return true;
}

bool P3DControl::appendTileSet(const char * filename, bool split, bool useWorldCoords)
{
	if (tileSet == NULL)
		return loadTileSet(filename, split, useWorldCoords);

	if (filename == NULL || 0 == strlen(filename))
		return false;

	if (! tileSet->byuFormat()) {
		cout << "First tile set loaded must have BYU format to merge tile sets"
			<< endl;
		return false;
	}

	Image3D * image;
	if (useWorldCoords) {
		image = imagePlanes.getImagePtr();
		if (image == NULL) {
			cout << "An image must exist to load BYU tiles in world coordinates"
				<< endl;
			beep();
			image = NULL;
		}
	}
	else {
		image = NULL;
	}

	if (tileSet->loadAppendFile(filename, split, image) == false)
		return false;

	return true;
}

//	Verify the specified files are loadable tilesets and determine
//	their range of types (min and max are set to 3 or 4).
/*bool P3DControl::evaluateTileSetTypes(int numFiles, const char ** filenames,
int & min, int & max)
{
BYU byu;
BYU::BYUPoly_t type1;
BYU::BYUPoly_t type2;
bool done;

type1 = byu.type(filenames[0]);
if (type1 == BYU::NoTiles)
return false;
cout << "File " << filenames[0] << " is type ";
if (type1 == BYU::Tris) {
min = 3;
max = 3;
cout << "3\n";
}
else {
min = 4;
max = 4;
cout << "4\n";
}

done = false;
for (int i = 1; i < numFiles; i++) {
type2 = byu.type(filenames[i]);
if (type2 == BYU::NoTiles)
return false;

cout << "File " << filenames[i] << " is type ";
if (type2 == BYU::Tris)
cout << "3\n";
else
cout << "4\n";

if (! done && type1 != type2) {
min = 3;
max = 4;
done = true;
// Continue looping to make sure all file names are valid paths
}
}
return true;
}*/
bool P3DControl::evaluateTileSetTypes(int numFiles, const char ** filenames,
									  int & min, int & max)
{
	TileSet ts;
	int type1, type2;
	bool isBYU, done;

	type1 = ts.type(filenames[0], isBYU);
	// Presently, .til files cannot be merged, so fail here
	if (! isBYU) {
		cout << "Old-style .til files cannot be loaded with other tile sets" << endl;
		return false;
	}
	min = type1;
	max = type1;
	if (type1 == 0)
		return false;

	done = false;
	for (int i = 1; i < numFiles; i++) {
		type2 = ts.type(filenames[i], isBYU);
		// Presently, .til files cannot be merged, so fail here
		if (! isBYU) {
			cout << "Old-style .til files cannot be loaded with other tile sets" << endl;
			return false;
		}
		if (type2 == 0)
			return false;

		if (! done && type1 != type2) {
			min = 3;
			max = 4;
			done = true;
			// Continue looping to make sure all file names are valid paths
		}
	}
	return true;
}

// This is presently used for loading PCA data
bool P3DControl::loadAuxTileSet(const char * filename)
{
	BYU * auxBYU;	// This mirrors auxTileSet

	if(filename == NULL || 0 == strlen(filename))
		return false;

	auxBYU = (BYU *) auxTileSet;
	if(auxBYU != NULL) {
		auxBYU->clear();
		if (auxBYU->read(filename) == false)
			return false;
	}
	else {
		auxBYU = new BYU(filename);
		auxTileSet = auxBYU;
		if (auxBYU->tileListLength() == 0)
			return false;
	}

	return true;
}

void P3DControl::clearTileSet()
{
	if (tileSet)
	{
		delete tileSet;
		tileSet = NULL;
	}
	if (!tileSet)
		tileSet = new TileSet;
}

// Tile the image, discarding any existing TileSet
bool P3DControl::generateTileSet()
{
	Image3D * image;

	image = imagePlanes.getImagePtr();
	if (image == NULL)
		return false;

#ifdef BINARY
	if (globalVerbosity > -1)
		cout << "Replacing tiles by generated tiles..." << flush;
#else
	cout << "Generating tileset " << flush;
#endif
	if (! tileSet)
		tileSet = new TileSet;

	if (! tileSet->simpleMarchingCubes(image)) {	// Should not fail
		beep();
		return false;
	}
	//  tileSet->print(4);
#ifdef BINARY
	if (globalVerbosity > -1)
#endif
		cout << " done" << endl;

	return true;
}

// When useWorldCoords is true, an image is expected to be present
bool P3DControl::saveTileSet(const char * filename, bool useWorldCoords)
{
	Image3D * image;

	if (filename == NULL || 0 == strlen(filename))
		return false;

	if (tileSet == NULL)
		return false;
	if (tileSet->tileCount() == 0)
		return false;

	if (useWorldCoords) {
		image = imagePlanes.getImagePtr();
		if (image == NULL) {
			cout << "An image must exist to save BYU tiles in world coordinates"
				<< endl;
			beep();
			return false;
		}
	}
	else
		image = NULL;

	bool quads = globalControl->readBool(BYUOutputType);
	return tileSet->saveFile(filename, quads, image);
}

bool P3DControl::haveTileSet() const
{
	if (tileSet == NULL)
		return false;
	if (tileSet->tileCount() == 0)
		return false;
	return true;
}

int P3DControl::addQuadFigure(int numRows, int numColumns, char * figureName,
							  float * color)
{
	M3DFigure * figure;

	int makeLandmark = 0;
	// GSTHACK: if name contains "LM-", make a special landmark figure
	int offset = 0;				// char offset from "LM-foo" to "foo"
	if (strstr(figureName, "LM-"))  // == figureName

	{

		offset = 3;
		makeLandmark = 1;
	}
	// end of GSTHACK


	const float white[] = { 1.0, 1.0, 1.0 };

	if(editingMode != ANY_EDITING_MODE)
		return -1;

	if(object == NULL)
		object = new M3DObject;

	startOperation();

	if( numRows == 1 ) {
		M3DTubeFigure* tube = new M3DTubeFigure;
		if (makeLandmark)
			tube->makeLandmarkStockFigure(2);
		else 
			tube->makeEllipticalStockFigure(numColumns);
		figure	= tube;
	}
	else {
		M3DQuadFigure* quad = new M3DQuadFigure;
		if (makeLandmark)
			quad->makeLandmarkStockFigure(2,2);
		else
			quad->makeEllipticalStockFigure(numRows, numColumns);
		figure	= quad;
	}

	if (figureName && *figureName)
		figure->setName(figureName + offset);	// GSTHACK: offset
	if (color[0] != 0.0 || color[1] != 0.0 || color[2] != 0.0)
		figure->setColor(color);
	else
		figure->setColor(white);
	object->deselect();
	figure->select();

	object->addFigure((M3DFigure*) figure);

	endOperation();

	modelNeedsSaving();
	selectionChanged = true;

	return object->getFigureCount() - 1;    // Number of the added figure   
}

// Removes all selected figures
void P3DControl::removeSelectedFigures()
{
	M3DFigure * figurePtr;
	int numFigures;
	int i;
	int markedFigure, primIndex;

	if(editingMode != ANY_EDITING_MODE)
		return;

	if(object == NULL)
		return;

	startOperation();
	getRelativeMarkedPrimitiveId(markedFigure, primIndex);

	numFigures = object->getFigureCount();
	for(i = numFigures - 1; i >= 0; i--)
	{
		figurePtr = object->getFigurePtr(i);
		if(figurePtr != NULL && figurePtr->isSelected()) {
			if (i == markedFigure)
				markedPrimitiveId = INVALID_PRIMITIVE_ID;
			else if (i < markedFigure)
				markedPrimitiveId -= figurePtr->getPrimitiveCount();
			object->deleteFigure(i);
		}
	}
	if (markedPrimitiveId < 0)
		markedPrimitiveId = INVALID_PRIMITIVE_ID;

	numFigures = object->getFigureCount();
	if (numFigures == 0)
	{
		// Empty objects are not allowed
		delete object;
		object = NULL;
		if (! attachedStats) {
			delete pgaStats;
			delete atomPgaStats;
		}
		pgaStats = NULL;
		atomPgaStats = NULL;			
		cpnsStats = NULL ;		// cpns stats attached to the object by default

		if (globalVerbosity >= 0)
			cout << "All figures removed: model discarded" << endl;

		if (referenceObject != NULL) {
			delete referenceObject;
			referenceObject = NULL;
		}
		markedPrimitiveId = INVALID_PRIMITIVE_ID;

		modelIsSaved();
		shiftOnPaste = false;
	}
	else {
		modelNeedsSaving(true);
		selectionChanged = true;
	}

	endOperation();
}


// Check whether or not the entire sub-figure tree has been selected
bool P3DControl::isFigureSubTreeSelected(M3DFigureTreeNode * node)
{
	if(!object->getFigurePtr(node->getFigureId())->isSelected())
		return false;

	M3DFigureTreeNode * child;
	int childCount, i;

	childCount = node->getChildCount();
	for(i = 0; i < childCount; i++)
	{
		child = node->getChild(i);
		if (! isFigureSubTreeSelected(child))
			return false;
	}

	return true;
}

// Absolute translate, rotate & scale
void P3DControl::translate(const Vector3D &trans)
{
	if(editingMode != ANY_EDITING_MODE &&
		editingMode != SIMILARITY_ONLY_EDITING_MODE)
		return;

	if(object)
	{

#ifdef BINARY
		// we are NOT clear as to what is going on here. why is figure #0 the only
		// figure selected, instead of any figure in the tree? -gst/jl

		int figureId = (int) tuningWt(BpFigureId); // Multiobject: figureId comes from tuning vals
#endif

		M3DFigure *figure;
		M3DFigure *oldFigure;
#ifdef BINARY
		figure = object->getFigurePtr(figureId); // WHY? who calls this routine, anyway?
#else
		figure = object->getFigurePtr(0);	// figureId=0 is hard-coded right now
#endif
		oldFigure = figure->assign();

		startOperation();

		object->translateBy(trans);
#ifndef BINARY
		handPlacementTransform->translateBy(trans);
#endif

		if(recordingTransformation)
			recordedTransformation->translateBy(trans);

#ifndef BINARY
		// Multiobject / Multifigure question
		// This is correct for the single object-single figure case
		// figureID will be 0 and figure tree will have no children
		// 
		// This is correct for the multiple-single object case
		// each figureID is also the id for a tree with no children
		// 
		// This is correct for the single-object, multiple figure case
		// figureId will be 0 and it will have children
		//
		// In the future, if we have multiple objects with multiple figures
		// I'm not sure if this will be right or not.  I don't know how the figures 
		// will be ordered.
#endif
		M3DFigureTreeNode *figureTreeNode, *subFigNode;
		int i;
#ifndef BINARY
		int figureId = 0;				// hard coded
#endif

		// Han: this operation is not necessary the entire object is translated or rotated at the same time.
		M3DSubfigureTransformation subfigure_xform;
		// Han: work to do: multi-level compatibility
		figureTreeNode = object->getFigureTreeNode(figureId);
		if (figure->isAnySelected() && ! isFigureSubTreeSelected(figureTreeNode))
		{
			figureTreeNode = object->getFigureTreeNode(figureId);

			i = 0;
			subFigNode = figureTreeNode->getChild(i);
			while (subFigNode != NULL)
			{
				subfigure_xform.init(object, subFigNode);
				subfigure_xform.updateSubfigure(dynamic_cast<M3DQuadFigure*>(oldFigure));
				i++;
				subFigNode = figureTreeNode->getChild(i);
			}
		}

		modelNeedsSaving(true);
		selectionChanged = true;
	}

	if (tileSet && tileSetMobile) {
		tileSet->translateBy(trans, tileSetMobileFigure);
		tileSet->setGenerated();
	}
}

void P3DControl::rotate(Quat &q)
{
	if(editingMode != ANY_EDITING_MODE &&
		editingMode != SIMILARITY_ONLY_EDITING_MODE)
		return;

	if(object)
	{
		M3DFigure *figure;
		M3DFigure *oldFigure;
#ifdef BINARY
		int figureId = (int) tuningWt(BpFigureId); // Multiobject: figureId comes from tuning vals
		figure = object->getFigurePtr(figureId);
#else
		figure = object->getFigurePtr(0);	// figureId=0 is hard-coded right now
#endif
		oldFigure = figure->assign();

		startOperation();

		object->rotateBy(q);

#ifndef BINARY
		// Set handplacementTransform
		Vector3D gravcenter = object->getCOG();
		handPlacementTransform->translateBy(-gravcenter);
		handPlacementTransform->rotateBy(q);
		handPlacementTransform->translateBy(gravcenter);
#endif

		if(recordingTransformation)
		{
			Vector3D cog = object->getCOG();
			recordedTransformation->translateBy(-cog);
			recordedTransformation->rotateBy(q);
			recordedTransformation->translateBy(cog);
		}

		M3DFigureTreeNode *figureTreeNode, *subFigNode;
		int i;
#ifndef BINARY
		int figureId = 0; // hard coded - multiple figures are not yet being handled
#else
		// MultiObject: See the extended comment in translate
#endif

		// This operation is not necessary; the entire object is translated or rotated at the same time.
		M3DSubfigureTransformation subfigure_xform;
		// work to do: multi-level compatibility
		figureTreeNode = object->getFigureTreeNode(figureId);
		if (figure->isAnySelected() && ! isFigureSubTreeSelected(figureTreeNode))
		{
			figureTreeNode = object->getFigureTreeNode(figureId);

			i = 0;
			subFigNode = figureTreeNode->getChild(i);
			while (subFigNode != NULL)
			{
				subfigure_xform.init(object, subFigNode);
				subfigure_xform.updateSubfigure(dynamic_cast<M3DQuadFigure*>(oldFigure));
				i++;
				subFigNode = figureTreeNode->getChild(i);
			}
		}

		modelNeedsSaving(true);
		selectionChanged = true;
	}
}

void P3DControl::rotateTileSet(Quat &q)
{
	if (tileSet && tileSetMobile) {
		tileSet->rotateBy(q, tileSetMobileFigure);
		tileSet->setGenerated();
	}
}

void P3DControl::pirouetteTileSet(Quat &q)
{
	if(editingMode != ANY_EDITING_MODE &&
		editingMode != SIMILARITY_ONLY_EDITING_MODE)
		return;

	if (tileSet && tileSetMobile) {
		tileSet->pirouetteBy(q, tileSetMobileFigure);
		tileSet->setGenerated();
	}
}


// Scale about center of gravity
void P3DControl::scale(double val)
{
	if(editingMode != ANY_EDITING_MODE &&
		editingMode != SIMILARITY_ONLY_EDITING_MODE)
		return;

	if(object)
	{
		M3DFigure *figure;
		M3DFigure *oldFigure;
#ifdef BINARY
		int figureId = (int) tuningWt(BpFigureId); // Multiobject: figureId comes from tuning vals
		figure = object->getFigurePtr(figureId);
#else
		figure = object->getFigurePtr(0);	// figureId=0 is hard-coded right now
#endif
		oldFigure = figure->assign();

		startOperation();

		object->scaleBy(val);

#ifndef BINARY
		Vector3D gravcenter = object->getCOG();
		handPlacementTransform->translateBy(-gravcenter);
		handPlacementTransform->scaleBy(val);
		handPlacementTransform->translateBy(gravcenter);
#endif

		if(recordingTransformation)
		{
			Vector3D cog = object->getCOG();
			recordedTransformation->translateBy(-cog);
			recordedTransformation->scaleBy(val);
			recordedTransformation->translateBy(cog);
		}

		M3DFigureTreeNode *figureTreeNode, *subFigNode;
		int i;
#ifdef BINARY
		// MultiObject: See the extended comment in translate
#else
		int figureId = 0; // hard coded
#endif

		// this operation is not necessary the entire object is translated or rotated at the same time.
		M3DSubfigureTransformation subfigure_xform;
		// work to do: multi-level compatibility
		figureTreeNode = object->getFigureTreeNode(figureId);
		if (figure->isAnySelected() && ! isFigureSubTreeSelected(figureTreeNode))
		{
			figureTreeNode = object->getFigureTreeNode(figureId);

			i = 0;
			subFigNode = figureTreeNode->getChild(i);
			while (subFigNode != NULL)
			{
				subfigure_xform.init(object, subFigNode);
				subfigure_xform.updateSubfigure(dynamic_cast<M3DQuadFigure*>(oldFigure));
				i++;
				subFigNode = figureTreeNode->getChild(i);
			}
		}

		modelNeedsSaving(true);
		selectionChanged = true;
	}

	if (tileSet && tileSetMobile) {
		tileSet->scaleBy(val, tileSetMobileFigure);
		tileSet->setGenerated();
	}

}


void P3DControl::scaleWidth(double val)
{
	if(editingMode != ANY_EDITING_MODE)
		return;

	if(object == NULL)
		return;

	startOperation();

	object->scaleWidth(val);

	modelNeedsSaving(true);
	selectionChanged = true;
}

void P3DControl::elongateRow(double val)
{
	if(editingMode != ANY_EDITING_MODE)
		return;

	if (! initializeElongater())
		return;

	if(elongater == NULL)
		return;

	startOperation();

	//	cout << "P3DControl::elongateRow(" << val << ")\n";
	elongater->setElongationMethod(M3D_ELONGATE_ROW);
	elongater->elongate(val);
	endOperation();
	modelNeedsSaving(true);
	selectionChanged = true;
	elongater->update();
}

void P3DControl::elongateColumn(double val)
{
	if(editingMode != ANY_EDITING_MODE)
		return;

	if (! initializeElongater())
		return;

	if (elongater == NULL)
		return;

	startOperation();

	//	cout << "P3DControl::elongateColumn(" << val << ")\n";
	elongater->setElongationMethod(M3D_ELONGATE_COLUMN);
	elongater->elongate(val);
	endOperation();
	modelNeedsSaving(true);
	selectionChanged = true;
	elongater->update();
}

bool P3DControl::initializeElongater()
{
	M3DFigure * figure;
	int figureId, atomId;

	if(elongater != NULL) {
		delete elongater;
		elongater = NULL;
	}

	if (object == NULL)
		return false;

	getRelativeMarkedPrimitiveId(figureId, atomId);
	if (figureId < 0)
		return false;

	// FIXME: What should we do in case of tubes? - gst/rrs
	figure = object->getFigurePtr(figureId);
	elongater = new M3DFigureElongater(dynamic_cast<M3DQuadFigure*>(figure));
	return true;
}

void P3DControl::setTheta(double val)
{
	M3DPrimitive * markedPrimitive;

	if(editingMode != ANY_EDITING_MODE)
		return;

	if(object == NULL)
		return;

	markedPrimitive = getMarkedPrimitive();
	if(markedPrimitive == NULL)
		return;

	// Commented out by AGG.  Note that P3DUserInterfaceCallback::updateMarkedPrimitive()
	// only requires there to be a marked primitive for the Atom Editor to be active.
	// If this is restored, that function should be made consistent with this one.
	//    if(markedPrimitive->isSelected() == false)
	//        return;

	startOperation();

	markedPrimitive->setTheta(val);

	endOperation();
	modelNeedsSaving(true);
	selectionChanged = true;
}

void P3DControl::setElongation(double val)
{
	M3DPrimitive * markedPrimitive;

	if(editingMode != ANY_EDITING_MODE)
		return;

	if(object == NULL)
		return;

	markedPrimitive = getMarkedPrimitive();
	if(markedPrimitive == NULL || markedPrimitive->type() != M3D_END_PRIMITIVE)
		return;

	// Commented out by AGG.  Note that P3DUserInterfaceCallback::updateMarkedPrimitive()
	// only requires there to be a marked primitive for the Atom Editor to be active.
	// If this is restored, that function should be made consistent with this one.
	//    if(markedPrimitive->isSelected() == false)
	//        return;

	startOperation();

	dynamic_cast<M3DEndPrimitive*>(markedPrimitive)->setElongation(val);

	endOperation();
	modelNeedsSaving(true);
	selectionChanged = true;
}

void P3DControl::toggleFigure(int figId)
{
	M3DFigure * figurePtr;

	if(object == NULL || ! selectionIsAllowed())
		return;

	figurePtr = object->getFigurePtr(figId);
	if (figurePtr != NULL) {
		if (figurePtr->isSelected())
			figurePtr->deselect();
		else
			figurePtr->select();
	}
	selectionChanged = true;
}

void P3DControl::togglePrimitive(int primId)
{
	M3DPrimitive * primitivePtr;

	if(object == NULL || ! selectionIsAllowed())
		return;

	primitivePtr = object->getPrimitivePtr(primId);
	if(primitivePtr != NULL)
		primitivePtr->toggleSelected();
	selectionChanged = true;
}

void P3DControl::selectAll()
{
	if(object == NULL || ! selectionIsAllowed())
		return;

	object->select();
	selectionChanged = true;
}

void P3DControl::deselectAll()
{
	if(object == NULL || ! selectionIsAllowed())
		return;

	object->deselect();
	selectionChanged = true;
}

void P3DControl::toggleAll()
{
	if(object == NULL || ! selectionIsAllowed())
		return;

	object->toggleSelected();
	selectionChanged = true;
}

// Copies the selected figures into the copy buffer
void P3DControl::copy()
{
	M3DFigure * figurePtr;
	M3DFigure * newFigurePtr;
	int numFigures, numTrees;
	M3DFigureTreeNode * node;
	M3DFigureTreeNode * parent;
	M3DFigureTreeNode * copiedNode;
	M3DPrimitiveLinkInfo * newLink;
	M3DPrimitiveLinkInfo * oldLink;
	int * copiedFigures;	// Source & target figure Id's of copied figures
	int numCopiedFigures;	// Twice the number of copied figures
	int childId;
	int newChildId;


	if(object == NULL)
		return;

	if(objectCopyBuffer != NULL)
	{
		delete objectCopyBuffer;
		objectCopyBuffer = NULL;
	}

	objectCopyBuffer = new M3DObject;
	shiftOnPaste = true;

	numFigures = object->getFigureCount();
	copiedFigures = new int[numFigures + numFigures];	// Source and target figure Id's
	numCopiedFigures = 0;

	// Copy the selected figures to the buffer
	newChildId = 0;
	for (childId = 0; childId < numFigures; childId++)
	{
		figurePtr = object->getFigurePtr(childId);

		if (figurePtr == NULL || !(figurePtr->isSelected()))
			continue;

		newFigurePtr = figurePtr->assign();  // Doesn't duplicate the figure name (see paste())
		objectCopyBuffer->addFigure(newFigurePtr);	// Creates one tree node per figure
		copiedFigures[numCopiedFigures++] = childId;
		copiedFigures[numCopiedFigures++] = newChildId;
		newChildId++;
	}

	// Copy any existing connectivity between the selected figures
	numTrees = object->getFigureTreeCount();
	if (numTrees == numFigures)
		return;		// All figures are singletons

	for (childId = 0; childId < numFigures; childId++)
	{
		int parentId;
		int newParentId;
		int i;
		bool copy_attachment;

		figurePtr = object->getFigurePtr(childId);

		if (figurePtr == NULL || !(figurePtr->isSelected()))
			continue;

		// Figure childId was copied in the loop above
		node = object->getFigureTreeNode(childId);
		parent = node->getParent();
		if (parent == NULL)
			continue;

		// If the copied figure is a child of a copied parent, copy the attachment
		parentId = parent->getFigureId();
		copy_attachment = false;
		for (i = 0; i < numCopiedFigures; i += 2)
			if (copiedFigures[i] == parentId) {
				// The parent was copied, so must copy the attachment
				newParentId = copiedFigures[i + 1];
				copy_attachment = true;
				break;
			}
			if (copy_attachment) {
				// Look up the new Id of the copied child
				for (i = 0; i < numCopiedFigures; i += 2)
					if (copiedFigures[i] == childId) {
						newChildId = copiedFigures[i + 1];
						break;
					}
					// Fill in and attach the new child's node
					copiedNode = objectCopyBuffer->getFigureTreeNode(newChildId);
					copiedNode->setBlendAmount(node->getBlendAmount());
					copiedNode->setBlendExtent(node->getBlendExtent());
					for (i = 0; i < node->getLinkCount(); i++) {
						newLink = new M3DPrimitiveLinkInfo;
						oldLink = node->getLink(i);
						*newLink = *oldLink;
						copiedNode->addLink(newLink);
					}
					objectCopyBuffer->attachFigureTreeNode(copiedNode, newParentId);
			}
	}
}

// Pastes the contents of the copy buffer into the model
void P3DControl::paste()
{
	M3DFigure * figurePtr;
	M3DFigure * newFigurePtr;
	int numFigures, numTrees;
	M3DFigureTreeNode * node;
	M3DFigureTreeNode * parent;
	M3DFigureTreeNode * copiedNode;
	M3DPrimitiveLinkInfo * newLink;
	M3DPrimitiveLinkInfo * oldLink;
	int * copiedFigures;	// Source & target figure Id's of copied figures
	int numCopiedFigures;	// Twice the number of copied figures
	int childId;
	int newChildId;


	if(editingMode != ANY_EDITING_MODE)
		return;

	if(object == NULL || objectCopyBuffer == NULL)
		return;

	startOperation();

	object->deselect();
	numFigures = objectCopyBuffer->getFigureCount();

	copiedFigures = new int[numFigures + numFigures];	// Source and target figure Id's
	numCopiedFigures = 0;

	// Use R of the first primitive of the object for an arbitrary shift amount
	Vector3D trans;
	if (shiftOnPaste && numFigures > 0) {
		figurePtr = objectCopyBuffer->getFigurePtr(0);
		M3DPrimitive * prim = figurePtr->getPrimitivePtr(0);
		double t = prim->getR()/3.0;
		trans = Vector3D(t, t, t);
	}

	// Paste the selected figures from the buffer
	newChildId = object->getFigureCount();
	for (childId = 0; childId < numFigures; childId++)
	{
		figurePtr = objectCopyBuffer->getFigurePtr(childId);

		if (figurePtr == NULL)
			continue;

		newFigurePtr = figurePtr->clone();  // Duplicates the figure's name (see copy())
		newFigurePtr->setModified(true);
		if (shiftOnPaste)
			newFigurePtr->translateBy(trans);
		object->addFigure(newFigurePtr);
		copiedFigures[numCopiedFigures++] = childId;
		copiedFigures[numCopiedFigures++] = newChildId;
		newChildId++;
	}

	// Copy any connectivity from the buffer to the copied figures
	numTrees = objectCopyBuffer->getFigureTreeCount();
	if (numTrees != numFigures) {
		for (childId = 0; childId < numFigures; childId++)
		{
			int parentId;
			int newParentId;
			int i;
			bool copy_attachment;

			figurePtr = objectCopyBuffer->getFigurePtr(childId);

			if (figurePtr == NULL)
				continue;

			// Figure childId was copied in the loop above
			node = objectCopyBuffer->getFigureTreeNode(childId);
			parent = node->getParent();
			if (parent == NULL)
				continue;

			// If the copied figure is a child of a copied parent, copy the attachment
			parentId = parent->getFigureId();
			copy_attachment = false;
			for (i = 0; i < numCopiedFigures; i += 2)
				if (copiedFigures[i] == parentId) {
					// The parent was copied, so must copy the attachment
					newParentId = copiedFigures[i + 1];
					copy_attachment = true;
					break;
				}
				if (copy_attachment) {
					// Look up the new Id of the copied child
					for (i = 0; i < numCopiedFigures; i += 2)
						if (copiedFigures[i] == childId) {
							newChildId = copiedFigures[i + 1];
							break;
						}
						// Fill in and attach the new child's node
						copiedNode = object->getFigureTreeNode(newChildId);
						copiedNode->setBlendAmount(node->getBlendAmount());
						copiedNode->setBlendExtent(node->getBlendExtent());
						for (i = 0; i < node->getLinkCount(); i++) {
							newLink = new M3DPrimitiveLinkInfo;
							oldLink = node->getLink(i);
							*newLink = *oldLink;
							copiedNode->addLink(newLink);
						}
						object->attachFigureTreeNode(copiedNode, newParentId);
				}
		}
	}

	endOperation();

	modelNeedsSaving(true);
	shiftOnPaste = true;
}

// Combines mergeObject into the existing object
void P3DControl::merge(M3DObject * mergeObject)
{
	int numAddedTrees;
	M3DFigureTreeNode * node;
	int * correspondence;
	int clen;
	M3DFigure * newFigure;
	int numfigs, orig_numfigs;


	if (mergeObject == NULL)
		return;

	if (object == NULL) {
		*object = *mergeObject;
		return;
	}

	// To adjust constraints, we must have a map of the figure ID changes
	numfigs = mergeObject->getFigureCount();
	correspondence = new int[numfigs];
	clen = numfigs;
	orig_numfigs = object->getFigureCount();

	numAddedTrees = mergeObject->getFigureTreeCount();
	for (int treeId = 0; treeId < numAddedTrees; treeId++)
	{
		node = mergeObject->getFigureTreeRoot(treeId);
		object->addTree(node, mergeObject, correspondence);
	}

	// Correct the figureID's of the constraints of added figures.  This could
	// not have been done earlier, because the order of figures may differ in
	// the *mergeObject from what it now is in *object.
	numfigs = object->getFigureCount();
	for (int newFigureId = orig_numfigs; newFigureId < numfigs; newFigureId++) {
		newFigure = object->getFigurePtr(newFigureId);
		InterfiguralConstraints & constraints = newFigure->constraints();
		constraints.remapFigureIds(correspondence, clen);
	}
	delete [] correspondence;

	// Make object->that the same as object since restore() is later called
	M3DObject * backup = object->loadedObject();
	*backup = *object;

	modelNeedsSaving(true);
}

#ifdef LM_METHOD_OBJ

// Combines mergeObject into the existing object
void P3DControl::mergeLandmarkModel(M3DObject * mergeObject)
{
	int numAddedTrees;
	M3DFigureTreeNode * node;
	int * correspondence;
	int clen;
	M3DFigure * newFigure;
	int numfigs, orig_numfigs;


	if (mergeObject == NULL)
		return;

	if (landmarkObject == NULL) {
		*landmarkObject = *mergeObject;
		return;
	}

	// To adjust constraints, we must have a map of the figure ID changes
	numfigs = mergeObject->getFigureCount();
	correspondence = new int[numfigs];
	clen = numfigs;
	orig_numfigs = object->getFigureCount();

	numAddedTrees = mergeObject->getFigureTreeCount();
	for (int treeId = 0; treeId < numAddedTrees; treeId++)
	{
		node = mergeObject->getFigureTreeRoot(treeId);
		landmarkObject->addTree(node, mergeObject, correspondence);
	}

	// Correct the figureID's of the constraints of added figures.  This could
	// not have been done earlier, because the order of figures may differ in
	// the *mergeObject from what it now is in *object.
	numfigs = landmarkObject->getFigureCount();
	for (int newFigureId = orig_numfigs; newFigureId < numfigs; newFigureId++) {
		newFigure = landmarkObject->getFigurePtr(newFigureId);
	}
	delete [] correspondence;
}

#endif

// Mirrors the selected figures about the y-z plane
void P3DControl::mirror()
{
	M3DFigure    * figurePtr;
	M3DPrimitive * primitivePtr;

	int numFigures;
	int numPrimitives;

	Vector3D translateVector(0.0, 0.0, 0.0),
		posVector,
		newBVector;

	Quat flipAboutY;
	Quat q;

	double xCOG;
	int count,
		i,
		j;


	if(editingMode != ANY_EDITING_MODE)
		return;

	if(object == NULL)
		return;

	startOperation();

	flipAboutY.setAxisAngle(Vector3D(0.0, 1.0, 0.0), R_PI);

	xCOG = 0.0;
	count = 0;

	numFigures = object->getFigureCount();
	for(i = 0; i < numFigures; i++)
	{
		figurePtr = object->getFigurePtr(i);
		if(figurePtr == NULL || !(figurePtr->isSelected()))
			continue;

		xCOG += figurePtr->getCOG().getX();
		count++;
	}

	if(count != 0)
		xCOG = xCOG / (double) count;
	else
		return;

	for(i = 0; i < numFigures; i++)
	{
		figurePtr = object->getFigurePtr(i);
		if(figurePtr == NULL || !(figurePtr->isSelected()))
			continue;

		numPrimitives = figurePtr->getPrimitiveCount();
		for(j = 0; j < numPrimitives; j++)
		{
			primitivePtr = figurePtr->getPrimitivePtr(j);
			if(primitivePtr == NULL)
				continue;

			posVector = primitivePtr->getX();

			q = primitivePtr->getQ();
			q.setX(-q.getX());
			q.setW(-q.getW());

			primitivePtr->setQ(q * flipAboutY);

			newBVector = primitivePtr->getB();
			q.setAxisAngle(newBVector, R_PI);
			primitivePtr->rotateBy(q);

			// Translate primitive to opposite side
			translateVector.setX(-2.0 * posVector.getX());
			primitivePtr->translateBy(translateVector);
		}

		translateVector.setX(2.0 * xCOG);
		figurePtr->translateBy(translateVector);
	}

	endOperation();

	modelNeedsSaving(true);
	selectionChanged = true;
}

// Called when starting a transformation operation (adds to undo list)
int P3DControl::startOperation()
{
	int n = undoList.size();
	if (! inOperation) {
		undoList.addUndo(object, markedPrimitiveId);
		inOperation = true;
	}
	return n;
}

// Called when a transformation operation is completed
int P3DControl::endOperation()
{
	inOperation = false;
	return undoList.size();
}

// Undoes previous transformation operation
bool P3DControl::undo()
{
	if(editingMode != ANY_EDITING_MODE)
		return false;

	if(object == NULL)
		return false;

	M3DObject * obj = undoList.undoObject(object, markedPrimitiveId);
	if (obj == object)
		return false;	// There was nothing to undo

	object = obj;
	object->setModified(true);
	selectionChanged = true;

	//updateArbPlane(false);
	return true;
}

// Undoes previous operations without possibility of redo
void P3DControl::resetOperations(int npushes)
{
	int l = undoList.size();
	if (npushes > l)
		npushes = l;
	for (int i = 0; i < npushes; i++)
		object = undoList.deleteUndo(markedPrimitiveId);

	if(object != NULL)
		object->setModified(true);
}

// Redoes the previous undone transformation operation
bool P3DControl::redo()
{
	if(editingMode != ANY_EDITING_MODE)
		return false;

	if(object == NULL)
		return false;

	M3DObject * obj = undoList.redoObject(object, markedPrimitiveId);
	if (obj == object)
		return false;	// There was nothing to redo

	object = obj;
	object->setModified(true);
	selectionChanged = true;

	//updateArbPlane(false);
	return true;
}

void P3DControl::reorderSubModel(M3DFigureTreeNode * node, int * newFigureNums, int & n)
{
	newFigureNums[n++] = node->getFigureId();
	int nChildren = node->getChildCount();
	for (int c = 0; c < nChildren; c++) {
		M3DFigureTreeNode * child = node->getChild(c);
		if (child == NULL)
			continue;
		reorderSubModel(child, newFigureNums, n);
	}
}

// Generate new figure numbers in depth-first order
bool P3DControl::reorderModel()
{
	int * newFigureNums;
	int numFigures;
	int numTrees, t, n;
	M3DFigureTreeNode * node;

	if (object == NULL)
		return false;

	numFigures = object->getFigureCount();
	if (numFigures == 0)
		return false;

	newFigureNums = new int[numFigures];
	numTrees = object->getFigureTreeCount();
	n = 0;
	// Get figure numbers of the tree's nodes in depth-first order
	for (t = 0; t < numTrees; t++) {
		node = object->getFigureTreeRoot(t);
		newFigureNums[n++] = node->getFigureId();
		int nChildren = node->getChildCount();
		for (int c = 0; c < nChildren; c++) {
			M3DFigureTreeNode * child = node->getChild(c);
			if (child == NULL)
				continue;
			reorderSubModel(child, newFigureNums, n);
		}
	}
	//	for (int i = 0; i < numFigures; i++)
	//		cout << newFigureNums[i] << ' ';
	//	cout << endl;

	undoList.addUndo(object, markedPrimitiveId);	// Save original
	object->renumber(newFigureNums);

	return true;
}

bool P3DControl::testReorderSubModel(M3DFigureTreeNode * node, int & n)
{
	if (n++ != node->getFigureId())
		return true;
	int nChildren = node->getChildCount();
	for (int c = 0; c < nChildren; c++) {
		M3DFigureTreeNode * child = node->getChild(c);
		if (child == NULL)
			continue;
		if (testReorderSubModel(child, n))
			return true;
	}
	return false;
}

// Determine if figure numbers need to be put into depth-first order
bool P3DControl::testReorderModel()
{
	int numFigures;
	int numTrees, t, n;
	M3DFigureTreeNode * node;

	if (object == NULL)
		return false;

	numFigures = object->getFigureCount();
	if (numFigures == 0)
		return false;

	numTrees = object->getFigureTreeCount();
	n = 0;
	// Get figure numbers of the tree's nodes in depth-first order
	for (t = 0; t < numTrees; t++) {
		node = object->getFigureTreeRoot(t);
		if (n++ != node->getFigureId())
			return true;
		int nChildren = node->getChildCount();
		for (int c = 0; c < nChildren; c++) {
			M3DFigureTreeNode * child = node->getChild(c);
			if (child == NULL)
				continue;
			if (testReorderSubModel(child, n))
				return true;
		}
	}
	return false;
}

char * P3DControl::getModelFileName()
{
	char * name;

	if(modelFileName == NULL)
		return NULL;

	name = new char[strlen(modelFileName) + 1];
	strcpy(name, modelFileName);

	return name;
}

// Clear current model and start with a blank object
void P3DControl::newModel()
{
	if(editingMode != ANY_EDITING_MODE)
		return;

	if(modelFileName != NULL)
	{
		delete [] modelFileName;
		modelFileName = NULL;
	}

	if(object != NULL)
	{
		delete object;
		object = NULL;

		if (globalVerbosity >= 0) {
			cout << "Started a new model";
			if (pgaStats != NULL || atomPgaStats != NULL)
				cout << "; discarded geometric statistics";
			cout << endl;
		}

		if (! attachedStats) {
			delete pgaStats;
			delete atomPgaStats;
		}
		pgaStats = NULL;
		atomPgaStats = NULL;
		cpnsStats = NULL ;		// cpns stats attached to the object by default
	}
	// Note: figureless objects are not allowed
	if(referenceObject != NULL) {
		delete referenceObject;
		referenceObject = NULL;
	}
	markedPrimitiveId = INVALID_PRIMITIVE_ID;

	modelIsSaved();
	shiftOnPaste = false;
}

#ifdef LM_METHOD_OBJ

// Clear current model and start with a blank object
void P3DControl::newLandmarkModel()
{
	if(landmarkObject != NULL)
	{
		delete landmarkObject;
		landmarkObject = NULL;

		if (globalVerbosity >= 0)
			cout << "Started a new LANDMARK model" << endl;
	}

	modelIsSaved();
	shiftOnPaste = false;
	worldWasApplied = false;						// AGG: this is a bug: need a separate variable for landmark models
}

#endif

void P3DControl::newModel(M3DObject *newModel)
{
	if (newModel == NULL)
		return;

	if(object != NULL)
	{
		delete object;
		object = NULL;
		if (! attachedStats) {
			delete pgaStats;
			delete atomPgaStats;
		}
		pgaStats = NULL;
		atomPgaStats = NULL;
		cpnsStats = NULL ;
	}

	object = newModel;
	if(referenceObject != NULL) {
		delete referenceObject;
		referenceObject = NULL;
	}
	markedPrimitiveId = INVALID_PRIMITIVE_ID;

	modelIsSaved();
	shiftOnPaste = false;
	worldWasApplied = false;

	if (globalVerbosity >= 0)
		cout << "Replaced current model with another" << endl;
}

void P3DControl::applyWorld(Image3D * image)
{
	if (image == NULL)
		return;
	if (object == NULL)
		return;

	startOperation();
	if (object->applyWorld(image))
		cout << mappingMsg;
	else {
		endOperation();
		resetOperations(1);
		return;
	}
	endOperation();
	object->setModified(true);
	modelNeedsSaving(true);
	worldWasApplied = true;
}

#ifdef OPTIMIZATION_VISUALIZER

/*	Note: this function is intended only for temporary display
of models.  It does not delete objectSimTransform,
handPlacementTransform, or referenceObject.  The user
must save the object pointer and restore it after display
of the new model is finished.
*/
M3DObject * P3DControl::replaceObjectPtr(M3DObject * obj)
{
	M3DObject * o;

	o = object;
	object = obj; 
	if (object == NULL)
		return o;
	object->setModified(true);
	modelNeedsSaving(true); 
	return o;
}

#endif

// I/O routines to read/write a model from/to a .m3d file

void P3DControl::verifyOrientationMatches(M3DObject * model)
{
	if (verifyOrientation()) {
		if (model->orientation() == false) {
			cout << "Warning: this model has the Y-axis oriented positively;\n";
			cout << "    this model should not be fit with this version of ";
#ifdef AE2_BUILD
			cout << "ConStruct." << endl;
#else
			cout << "Pablo." << endl;
#endif
		}
	}
	else {
		if (model->orientation() == true) {
			cout << "Warning: this model has the Y-axis oriented negatively;\n";
			cout << "    this model should not be fit with this version of ";
#ifdef AE2_BUILD
			cout << "ConStruct." << endl;
#else
			cout << "Pablo." << endl;
#endif
		}
	}
}

/*  Read a model file.

If the xformFilename is provided, the transform will be read from
that file.  If the xformFilename is NULL, then the similarity
transform will be looked for in the model file.  If it is not
found there, it will be looked for in a separate .sim file with
the same path and name as the model file.

If a pgaFilename is provided, the PGA data will be read from that
file, replacing the PGA data that would normally be read from 
filename.  Then the model given by filename is read.  After it has
been read, if pgaModelFileName is set, PGA data will be read from
that file and the pgaStats and atomPgaStats pointers set to it.
Note that these two ways of getting auxiliary PGA data are not
the same, because the PGA operations in the M3DObjectFile::read()
function are only done for the first of them.
*/
bool P3DControl::read(const char * filename, bool mapToWorld,
					  const char * xformFilename, const char * pgaFilename)
{
    //cout << "Reading..." << endl;
	M3DObjectFile objectFile;
	SimilarityTransform3D * xform;
	M3DObject * obj;
	WorldSystem * world;
	Image3D * image;
	char * str;
	M3DObjectFile pgaObjectFile;
	M3DObject * pgaObj;

	if (filename == NULL)
		return false;

	if (editingMode != ANY_EDITING_MODE)
		return false;

	if (modelFileName != NULL)
		delete [] modelFileName;
	modelFileName = new char[strlen(filename) + 1];
	if (modelFileName == NULL)
		return false;
	strcpy(modelFileName, filename);

	xform = new SimilarityTransform3D;
	world = new WorldSystem;
	obj = objectFile.read(modelFileName, markedPrimitiveId,
		(xformFilename != NULL ? NULL : xform), world, pgaFilename);

    // cout << "Model file read successfully" << endl ;

	if (obj != NULL) {
		if (object != NULL) {
			delete object;
			if (! attachedStats) {
				delete pgaStats;
				delete atomPgaStats;
			}
			pgaStats = NULL;
			atomPgaStats = NULL;
			cpnsStats = NULL ;		// cpnsStats are tied to the object by default
			// A model must exist for there to be a similarity transform
			if (objectSimTransform != NULL)
				delete objectSimTransform;
			objectSimTransform = NULL;

			
			
		}
		object = obj;
		object->setWorld(world, true);	// Sets world in object and object->that
	}
	else {
		delete world;
		return false;
	}

	verifyOrientationMatches(object);

	if (pgaModelFileName != NULL) {   // pgaModelFileName is set by -ipr in pablo.cpp
		pgaObj = pgaObjectFile.read(pgaModelFileName);

		// Notice that in this case, the statistics in object are not replaced
		pgaStats = pgaObj->getPGAStats();
		atomPgaStats = pgaObj->getAtomPGAStats();
		// The pgaStats and atomPgaStats objects must be separately deleted
		attachedStats = false;
		pgaObj->setPGAStats(NULL, true);
		pgaObj->setAtomPGAStats(NULL, true);
		// Note: at this point, pgaStats and atomPgaStats, are not attached to
		// any instance of M3DObject
		delete pgaObj;
	}
	else {
		pgaStats = object->getPGAStats();
		atomPgaStats = object->getAtomPGAStats();
		// The pgaStats and atomPgaStats objects will be deleted with object
		attachedStats = true;
	}

	// dibyendu - cpns
	cpnsStats = object->getCPNSStats() ;
	if( cpnsStats )
		cout << endl << "CPNS Statistics loaded into the object" << endl ;

	str = NULL;
	if (xformFilename != NULL)
		xform->readSimilarity(xformFilename);
	else if (! xform->wasRead()) {
		int len = strlen(modelFileName);
		str = new char[len + 1];
		strncpy(str, modelFileName, len - 3);
		str[len - 3] = '\0';
		strcat(str, "sim");
		xform->readSimilarity(str);
		xformFilename = str;
	}

	if (xform->wasRead()) {
#ifdef DEBUG
		cout << "Similarity transformation loaded" << endl;
#endif
		if (xformFilename != NULL)
			cout << "Loaded similarity transformation " << xformFilename << endl;
	}
	else {			// The transformation in object is set to identity
		cout << "Similarity transformation not loaded" << endl;
	}
	objectSimTransform = xform;
	object->setTransformation(objectSimTransform);
	object->loadedObject()->setTransformation(objectSimTransform);

	if (referenceObject != NULL) {
		delete referenceObject;
		referenceObject = NULL;
	}

	modelIsSaved();
	shiftOnPaste = false;
	delete [] str;

	image = imagePlanes.getImagePtr();
	if (image == NULL)
		return true;

	if (mapToWorld) {
		if (! world->equals(image->getModelOrigin(), image->getModelBound(), image->getSpacing()))
			applyWorld(image);
	}
	else
		worldWasApplied = false;

	return true;
}

bool P3DControl::readAppend(const char * filename, bool mapToWorld, bool & ignoredPGA)
{
	M3DObjectFile objectFile;
	SimilarityTransform3D * xform;
	WorldSystem * world;
	Image3D * image;
	int m;

	if(object == NULL)
		return read(filename);

	if(filename == NULL)
		return false;

	if(editingMode != ANY_EDITING_MODE)
		return false;

	xform = NULL;
	world = new WorldSystem;
	M3DObject * newObject = objectFile.read(filename, m, xform, world);
	if (newObject == NULL)
		return false;
	if (markedPrimitiveId == INVALID_PRIMITIVE_ID)
		markedPrimitiveId = m;

	verifyOrientationMatches(newObject);

	// AGG: This function needs to handle similarity transforms; temporarily it merely ignores them

	// AGG: I've disabled this because it is no longer correct, in that it
	// doesn't work properly with the attachedStats flag, and it is also
	// probably not needed.
#if 0
	ignoredPGA = false;
	if (pgaStats == NULL)
		pgaStats = newObject->getPGAStats();
	else {
		if (newObject->getPGAStats() != NULL)
			ignoredPGA = true;
	}

	if (atomPgaStats == NULL)
		atomPgaStats = newObject->getAtomPGAStats();
	else {
		if (newObject->getAtomPGAStats() != NULL)
			ignoredPGA = true;
	}
#else
	ignoredPGA = true;
#endif

	merge(newObject);
	delete newObject;

	if (referenceObject != NULL) {
		delete referenceObject;
		referenceObject = NULL;
	}

	modelIsSaved();

	image = imagePlanes.getImagePtr();
	if (image == NULL)
		return true;

	bool worldWasAppliedTo1stModel = worldWasApplied;
	if (mapToWorld) {
		if (! world->equals(image->getModelOrigin(), image->getModelBound(), image->getSpacing()))
			applyWorld(image);
	}
	else
		worldWasApplied = false;
	worldWasApplied |= worldWasAppliedTo1stModel;
	/* Here, worldWasApplied is true, if either model was mapped.  This is
	equivalent to both models being mapped, since a model will not be
	mapped, if the worlds are identical: mapping from a world to the
	same world results in no change.  This is all needed so that
	unApplyWorld() can be applied if a second image is loaded after
	the model(s) are loaded.
	*/

	return true;
}

bool P3DControl::loadSimTransform(const char * xformFilename)
{
	SimilarityTransform3D * xform;
	bool ret;

	// The model must be loaded before any similarity transform
	if (object == NULL)
		return false;

	if (xformFilename != NULL) {
		xform = new SimilarityTransform3D;
		ret = xform->readSimilarity(xformFilename);
	}
	else
		return false;

	if (ret) {
		cout << "Loaded similarity transformation " << xformFilename << endl;
		if (objectSimTransform != NULL)
			delete objectSimTransform;
		objectSimTransform = xform;
		object->setTransformation(objectSimTransform);
		object->loadedObject()->setTransformation(objectSimTransform);
	}
	else
		cout << "Similarity transformation " << xformFilename << " not loaded" << endl;

	return ret;
}

#ifdef AE2_BUILD
void P3DControl::applySimTransform(SimilarityTransform3D * transform)
{
	if (transform != NULL)
		*objectSimTransform = *transform;
#else
void P3DControl::applySimTransform()
{
#endif
	if (objectSimTransform == NULL) {
		// Should never occur
		cout << "Error: applySimTransform: no transformation present\n";
		return;
	}

	if (object == NULL) {
		// Should never occur
		cout << "Error: applySimTransform: no model present\n";
		return;
	}

#ifdef UNFLIPPED
	object->applySimilarity(*objectSimTransform, Vector3D(0.0, 1.0, 0.0));
#else
	object->applySimilarity(*objectSimTransform, Vector3D(0.0, 0.0, 0.0));
#endif
	object->setModified(true);
	modelNeedsSaving(true);
	objectSimTransform->setToIdentity();
}

#ifdef LM_METHOD_OBJ

bool P3DControl::readLandmarkModel(const char* filename)
{
	M3DObjectFile objectFile;
	int lmarkedPrimitiveId;

	if (filename == NULL)
		return false;

	if (editingMode != ANY_EDITING_MODE)
		return false;

	/*    if (lmrkFileName != NULL)
	delete [] lmrkFileName;
	lmrkFileName = new char[strlen(filename) + 1];
	if (lmrkFileName == NULL)
	return false;
	strcpy(lmrkFileName, filename);*/

	if (landmarkObject != NULL)
		delete landmarkObject;

	//    landmarkObject = objectFile.read(lmrkFileName, markedPrimitiveId);
	landmarkObject = objectFile.read(filename, lmarkedPrimitiveId);
	if (landmarkObject == NULL)
		return false;

	shiftOnPaste = false;
	return true;
}

bool P3DControl::readAppendLandmarkModel(const char* filename)
{
	M3DObjectFile objectFile;

	if(landmarkObject == NULL)
		return read(filename);

	if(filename == NULL)
		return false;

	if(editingMode != ANY_EDITING_MODE)
		return false;

	M3DObject * newObject = objectFile.read(filename);
	if (newObject == NULL)
		return false;

	mergeLandmarkModel(newObject);

	return true;
}

#endif

// dibyendu
// Function to write out the distances at the spoke ends for all spokes in the current model to a text file
// This gives us a metric to compare model fittings quantitatively

bool P3DControl::writeDistancesAtSpokeEnds(const char *filename) {

	if( ( filename == NULL ) || ( object == NULL ) || ( match == NULL ) || ( match->getImageDistMapPointer() == NULL ) ) {
		cout << endl << "NULL arguments passed to P3DControl::writeDistancesAtSpokeEnds() !" << endl ;
		return( 0 ) ;
	}

	// Collect the distance data into a vector first	

	vector <double> distVector ;

	vector <double> rowNumVector ;
	vector <double> colNumVector ;
	vector <double> spokeNumVector ;

	// take the first figure, get the distances and put them in the vector

	M3DFigure * figThis = object->getFigurePtr(0) ;

	M3DQuadFigure * quadFigThis = dynamic_cast <M3DQuadFigure *> (figThis) ;

	if( quadFigThis == NULL ) {
		cout << endl << "P3DControl::writeDistancesAtSpokeEnds() does not work with Tube figures" << endl ;
		return (0) ;
	}

	double dilationFactor = object->getDilationFactorInModelUnits() ;

	int nPrims = quadFigThis->getPrimitiveCount() ;	

	Vector3D thisSpokeEnd ;

	int nRows = quadFigThis->getRowCount() ;
	int nCols = quadFigThis->getColumnCount() ;

	for( int ii = 0 ; ii < quadFigThis->getRowCount() ; ii ++ ) {

		for( int jj = 0 ; jj < quadFigThis->getColumnCount() ; jj ++ ) {			

			M3DPrimitive * primThis = quadFigThis->getPrimitivePtr( ii, jj ) ;			

			// spoke 0

			thisSpokeEnd = primThis->getX() + primThis->getY0() ;

			distVector.push_back( fabs( match->getImageDistMapPointer()->getDistance(thisSpokeEnd) - dilationFactor ) ) ;

			rowNumVector.push_back( ii ) ;	colNumVector.push_back( jj ) ; spokeNumVector.push_back( 0 ) ;

			// debug
			double d0 = fabs( match->getImageDistMapPointer()->getDistance(thisSpokeEnd) - dilationFactor ) ;

			// spoke 1

			thisSpokeEnd = primThis->getX() + primThis->getY1() ;

			distVector.push_back( fabs( match->getImageDistMapPointer()->getDistance(thisSpokeEnd) - dilationFactor ) ) ;

			rowNumVector.push_back( ii ) ;	colNumVector.push_back( jj ) ; spokeNumVector.push_back( 1 ) ;

			// debug
			double d1 = fabs( match->getImageDistMapPointer()->getDistance(thisSpokeEnd) - dilationFactor ) ;

			// spoke End

			if( primThis->type() == M3D_END_PRIMITIVE ) {		

				thisSpokeEnd = primThis->getX() + (dynamic_cast <M3DQuadEndPrimitive*> (primThis))->getYEnd() ;

				distVector.push_back( fabs( match->getImageDistMapPointer()->getDistance(thisSpokeEnd) - dilationFactor ) ) ;

				// debug
				double d2 = fabs( match->getImageDistMapPointer()->getDistance(thisSpokeEnd) - dilationFactor ) ;

				rowNumVector.push_back( ii ) ;	colNumVector.push_back( jj ) ; spokeNumVector.push_back( 2 ) ;				

			}

		}
	}

	// open the text file for writing 

	std::ofstream f;

	f.open( filename,  ios::trunc);

	if (!f) {
		string s("Unable to open " + string(filename) + " for writing spoke end distances!" ) ;
		cout << s << endl ;
		return(0) ;
	}	

	f << setw(10) << "Row No." << setw(10) << "Col. No." << setw(10) << "Spoke No." << setw(20) << "Distance" << endl ;

	for( int ii = 0 ; ii < distVector.size() ; ii ++ ) {

		f << setw(10) << rowNumVector[ii] ;
		f << setw(10) << colNumVector[ii] ;
		f << setw(10) << spokeNumVector[ii] ;
		f << setw(20) << distVector[ii] ;
		f << endl ;
	}

	f.close() ;

	distVector.clear() ;

	cout << endl << "Writing distances at spoke ends to " << filename << endl ;

	return( 1 ) ;	

}

/*  Write the model to a file.

If a filename is provided, it will become the new modelFileName.  If the
filename is NULL, then the model will be written to modelFileName.  If
modelFileName is NULL and no filename is provided, nothing will happen.

If there is no similarity transform, then no similarity transform output
will be produced.  Otherwise, if the xformFilename is NULL, then the
similarity transform will either be written into the model file or a
separate file, depending on the SimTransformSeparate preference.  When
writing to a separate file, if an xformFilename is provided, it will be
the name of the output file.  If not, then the modelFileName will be used,
but the extension will be ".sim".  The format of the written transform
is determined by the SimTransformMatrix preference.

This function returns false if no file is written.
*/

bool P3DControl::write(const char * filename, const char * xformFilename)
{
	M3DObjectFile objectFile;
	Image3D * image;
	char * str;
	char * xfptr;


	if (filename != NULL)
	{
		if(modelFileName != NULL)
			delete [] modelFileName;

		modelFileName = new char[strlen(filename) + 1];
		strcpy(modelFileName, filename);
	}
	else if (modelFileName == NULL)
		return false;

	bool writeMatrix = globalControl->readBool(SimTransformMatrix);
	bool separate = globalControl->readBool(SimTransformSeparate);

	str = NULL;
	if (objectSimTransform) {
		const char dot = '.';
		int len;

		if (xformFilename == NULL) {
			if (separate) {
				char * ext = strrchr(modelFileName, dot);
				if (ext == NULL)
					len = strlen(modelFileName);
				else
					len = ext - modelFileName;
				str = new char[len + 5];
				strncpy(str, modelFileName, len);
				str[len] = '\0';
				strcat(str, ".sim");
			}
			xfptr = str;
		}
		else
			xfptr = (char *) xformFilename;
	}
	else {
		xfptr = NULL;
		separate = true;
	}

	if (object != NULL) {
		//
		// Some info display for tube figures.
		//
		for (int i = 0 ; i != object->getFigureCount(); ++i) {
			M3DFigure* figure			= object->getFigurePtr(i);
			figure->fixGlobalConsistency();
			if (globalVerbosity >= 1) {
				cout << "Curviness     : " << figure->curviness() << endl;
				M3DTubeFigure* tube_fig	= dynamic_cast<M3DTubeFigure*>(figure);
				if (tube_fig != NULL) {
					cout << "rSrad penalty : " <<
						M3DAtomPredictorTube().getFigureRSradPenalty(tube_fig, -1, PEL_BOTH, tuningWt(RSradPenaltyNorm), tuningWt(RSradPenaltyThreshold) )
						<< endl;
				}
				M3DQuadFigure* quad_fig	= dynamic_cast<M3DQuadFigure*>(figure);
				if (quad_fig != NULL) {
					cout << "rSrad penalty : " <<
						M3DAtomPredictorQuad().getFigureRSradPenalty(quad_fig, -1, PEL_BOTH, tuningWt(RSradPenaltyNorm), tuningWt(RSradPenaltyThreshold) )
						<< endl;
				}
			}
		}
		image = imagePlanes.getImagePtr();
		//
		// FIXME: write PGA should be a switch?
		//
		//		if (objectFile.write(modelFileName, *object,
		//			(separate ? NULL : objectSimTransform), writeMatrix,
		//			markedPrimitiveId, image, imageFileName, ! modelIsModified()))
		if (objectFile.write(modelFileName, *object,
			(separate ? NULL : objectSimTransform), writeMatrix,
			markedPrimitiveId, image, imageFileName, true))
			modelIsSaved();
	}

	if (separate && xfptr != NULL) {
		if (objectSimTransform->writeSimilarity(xfptr, writeMatrix))
			cout << "Saved similarity transformation " << xfptr << endl;
		else
			cout << "Unable to save similarity transformation" << endl;
	}

	if (str != NULL)
		delete [] str;
	return true;
}

bool P3DControl::readOldFile(const char * filename)
{
	if(editingMode != ANY_EDITING_MODE)
		return false;

	if(filename == NULL || 0 == strlen(filename))
		return false;

	if(modelFileName != NULL)
		delete [] modelFileName;

	modelFileName = new char[strlen(filename) + 5];
	if(modelFileName == NULL)
		return false;

	strcpy(modelFileName, filename);
	strcat(modelFileName, ".m3d");

	if(object != NULL) {
		delete object;
		if (! attachedStats) {
			delete pgaStats;
			delete atomPgaStats;
		}
		pgaStats = NULL;
		atomPgaStats = NULL;
		cpnsStats = NULL ;
	}

	object = readOldObjectFile(filename);
	if (object == NULL)
		return false;

	modelNeedsSaving();
	shiftOnPaste = false;
	worldWasApplied = false;
	return true;
}

void P3DControl::exportToImage(const char * filename)
{

	if(object == NULL)
		return;
	if (! object->verifyInBounds()) {
		cout << "This model is larger than or outside the unit cube --\n"
			<< "    cannot export a correct image" << endl;
		beep();
		return;
	}

	DistanceMap3D distanceMap;
	BoundingBox bbox;
	bbox.origin = Vector3D(0.0, 0.0, 0.0);
	bbox.size = Vector3D(1.0, 1.0, 1.0);
	int map_size[3];
	map_size[0] = (int) tuningWt(OutModelImageX);
	map_size[1] = (int) tuningWt(OutModelImageY);
	map_size[2] = (int) tuningWt(OutModelImageZ);

	Image3D* sourceImage = this->getImagePtr();
	if (sourceImage != NULL) {
		tuneVals.setWeight(OutModelImageFixedSize, 1);
		bbox.size.setX(sourceImage->getXModelExtent());
		bbox.size.setY(sourceImage->getYModelExtent());
		bbox.size.setZ(sourceImage->getZModelExtent());
		map_size[0] = sourceImage->getXDim();
		map_size[1] = sourceImage->getYDim();
		map_size[2] = sourceImage->getZDim();
	}

	if (tuningWt(OutModelImageFixedSize)) {
		cout << "map size starts as "
			<< map_size[0] << " " 
			<< map_size[1] << " " 
			<< map_size[2] << " " 
			<< endl;
	}

	int sizeChanged = distanceMap.initializeWithBBox(
		object, &bbox, map_size,
		(int) tuningWt(OutModelImageSurfaceLevel),
		0.6,
		tuningWt(OutModelImageFixedSize) > 0 ? true : false);

	if (!tuningWt(OutModelImageFixedSize) && sizeChanged > 0) {
		cout << "map size changed from "
			<< map_size[0] << " " 
			<< map_size[1] << " " 
			<< map_size[2];
		map_size[0] = distanceMap.getMapSize()[0];
		map_size[1] = distanceMap.getMapSize()[1];
		map_size[2] = distanceMap.getMapSize()[2];
		cout << " to "
			<< map_size[0] << " " 
			<< map_size[1] << " " 
			<< map_size[2] << " " 
			<< endl;
	}

	distanceMap.createMap(true);
	if(distanceMap.status() == false) {
		cout << "Error: could not export an image" << endl;
		beep();
		return;
	}

	Image3D* image = new Image3D(map_size[0], map_size[1], map_size[2]);
	// FIXME: Make spacing and origin also an option
	//image->setSpacingAndOrigin(0.05, 0.05, 0.05 );
	if (sourceImage == NULL) {
		image->setSpacingAndOrigin(
			1.0 / (double)map_size[0],
			1.0 / (double)map_size[1],
			1.0 / (double)map_size[2]);
	} else {
		Vector3D origin = sourceImage->getWorldOrigin();
		image->setSpacingAndOrigin(
			sourceImage->getXSpacing(),
			sourceImage->getYSpacing(),
			sourceImage->getZSpacing(),
			&origin);
	}


	GreyValue* voxels = image->getVoxels();
	int i;
	int size = map_size[0] * map_size[1] * map_size[2];

	for(i = 0; i < size; i++)
	{
		if(distanceMap.getDistance(i) <= 0.0f)
#ifdef BINARY
			voxels[i] = 1;
#else
			voxels[i] = MAX_GREY_VALUE_WRITTEN;
#endif
		else
			voxels[i] = MIN_GREY_VALUE;
	}


	AllImageIO file;
	if (file.write(filename, *image, MIN_GREY_VALUE, MAX_GREY_VALUE_WRITTEN)) {
		if (globalVerbosity >= 0)
			cout << "Exported model to image " << filename << endl;
	}
	else {
		cout << "Unable to export model to image: write failed" << endl;
		beep();
	}

	delete image;
}

#ifndef PRODUCTION_VERSION

void P3DControl::exportDistanceMap(const char * filename)
{
	GreyValue * voxels;
	float min, max;
	BoundingBox bbox;

	if(object == NULL)
		return;
	if (! object->verifyInBounds()) {
		cout << "This model is larger than or outside the unit cube --\n"
			<< "    cannot export a correct distance map" << endl;
		beep();
		return;
	}

	static const int MAP_SIZE[3] = {128, 128, 128};

	DistanceMap3D distanceMap;
	bbox.origin = Vector3D(0.0, 0.0, 0.0);
	bbox.size = Vector3D(1.0, 1.0, 1.0);

	(void) distanceMap.initializeWithBBox(object, &bbox, MAP_SIZE, 6, 0.6, true);

	distanceMap.createMap(true);
	if (distanceMap.status() == false)
		return;

	distanceMap.range(min, max);
	if (globalVerbosity >= 1)
		cout << "Distance map range = (" << min << ", " << max << ")\n";
	const int * mapSize = distanceMap.getMapSize();
	if (globalVerbosity >= 1)
		cout << "Image size: " << mapSize[0] << ", " << mapSize[1] << ", "
		<< mapSize[2] << endl;
	Image3D image(mapSize[0], mapSize[1], mapSize[2]);
	image.setSpacingAndOrigin(1.0 / (double)mapSize[0], 1.0 / (double)mapSize[1],
		1.0 / (double)mapSize[2]);

	voxels = image.getVoxels();
	int i;
	int size = mapSize[0] * mapSize[1] * mapSize[2];

	for(i = 0; i < size; i++)
	{
		float d = distanceMap.getDistance(i);
		if (d < 0)
			voxels[i] = MIN_GREY_VALUE;
		else
			voxels[i] = (GreyValue) (MAX_GREY_VALUE_WRITTEN*d/max);
	}

	AllImageIO file;
	if (! file.write(filename, image, MIN_GREY_VALUE, MAX_GREY_VALUE_WRITTEN)) {
		cout << "Unable to write distance map image" << endl;
		beep();
	}
	if (globalVerbosity >= 0)
		cout << "Exported distance map to " << filename << endl;
}

#endif	/* PRODUCTION_VERSION */

/*
// This version of exportDistanceMap() extends significantly beyond the bounding box,
// to show that the "city block" interpolation code is working.
void P3DControl::exportDistanceMap(const char * filename)
{
GreyValue * voxels;
float min, max;
BoundingBox bbox;

if(object == NULL)
return;
if (! object->verifyInBounds()) {
cout << "This model is larger than or outside the unit cube --\n"
<< "    cannot export a correct distance map" << endl;
beep();
return;
}

static const int MAP_SIZE[3] = {128, 128, 128};

DistanceMap3D distanceMap(object, 6, 0.6);
distanceMap.initialize();

distanceMap.createMap(true);
if (distanceMap.status() == false)
return;

distanceMap.range(min, max);
cout << "Distance map range = (" << min << ", " << max << ")\n";
int * mapSize = distanceMap.getMapSize();
cout << "Distance map size: " << mapSize[0] << ", " << mapSize[1] << ", "
<< mapSize[2] << endl;
mapSize = (int *) &MAP_SIZE;
cout << "Image size: " << mapSize[0] << ", " << mapSize[1] << ", "
<< mapSize[2] << endl;
Image3D image(mapSize[0], mapSize[1], mapSize[2]);
image.setSpacing(1.0 / (double)mapSize[0], 1.0 / (double)mapSize[1],
1.0 / (double)mapSize[2]);

voxels = image.getVoxels();
int i;
int size = MAP_SIZE[0] * MAP_SIZE[1] * MAP_SIZE[2];

double x = 0.0;
double y = 0.0;
double z = 0.0;
i = 0;
for(z = 0.0; z < 1.0; z += 1.0/MAP_SIZE[2])
for(y = 0.0; y < 1.0; y += 1.0/MAP_SIZE[1])
for(x = 0.0; x < 1.0; x += 1.0/MAP_SIZE[0])
{
float d = distanceMap.getDistance(x, y, z);
voxels[i++] = 1000*(GreyValue) (0.5 + d);	// Arbitrary factor
}

AllImageIO file;
if (! file.write(filename, image)) {
cout << "Unable to write distance map image" << endl;
beep();
}
if (globalVerbosity >= 0)
cout << "Exported distance map to " << filename << endl;
}
*/

/*void P3DControl::exportToBYUTiles(const char * filename,
Image3D * imageForWorldCoords,		// may be NULL
bool quads, bool useWorldCoords,
int level)
{
ThallCode::Pointlist_server2 * pList;
M3DFigure * figure;
FILE * fp;
Bpoint * tileset;
int * tileCount;
int numtiles, numFigures;
int i, figureNo;

if(object == NULL)
return;

numFigures = object->getFigureCount();
tileCount = new int[numFigures];
pList = new ThallCode::Pointlist_server2[numFigures];

// Generate the surfaces
numtiles = 0;
for (figureNo = 0; figureNo < numFigures; figureNo++)
{
figure = object->getFigurePtr(figureNo);
if(figure == NULL)
continue;

Xferlist * xferList = convertM3DtoXfer(figure);
pList[figureNo].init(xferList);
pList[figureNo].ComputeSubdivBoundaryTiles(level);
delete [] (xferList->atomlist);
delete xferList;

// Get boundary displacement information
SubdivBoundary * boundary = figure->getBoundaryPtr();
if (boundary != NULL)
{
Displacements * disp = boundary->getDisplacements(level);
if (disp != NULL)
{
pList[figureNo].ComputeSubdivPointCloud(level);
pList[figureNo].applyDisplacements(disp->getNumPts(), disp->getVals());
}
}

pList[figureNo].subdivtileinfo(tileCount + figureNo, &tileset);
numtiles += tileCount[figureNo];
}

// Output the header lines
fp = fopen(filename, "wb");
// Numbers of: figures, vertices, quads/triangles, and edges.
if (quads)
fprintf(fp, "%7d %7d %7d %7d\n", numFigures, 4*numtiles, numtiles, 4*numtiles);
else
fprintf(fp, "%7d %7d %7d %7d\n", numFigures, 4*numtiles, 2*numtiles, 6*numtiles);
i = 1;
for (figureNo = 0; figureNo < numFigures; figureNo++) {
if (quads) {
fprintf(fp, "%7d %7d ", i, tileCount[figureNo] + i - 1);
i += tileCount[figureNo];
}
else {
fprintf(fp, "%7d %7d ", i, 2*tileCount[figureNo] + i - 1);
i += 2*tileCount[figureNo];
}
}
fprintf(fp, "\n");

// Output the vertex lines - scale to world coord's of image if present
if (imageForWorldCoords)
cout << "Scaling tiles to world coordinates of image" << endl;

for (figureNo = 0; figureNo < numFigures; figureNo++)
{
figure = object->getFigurePtr(figureNo);
if(figure == NULL)
continue;

pList[figureNo].subdivtileinfo(&numtiles, &tileset);

for(i = 0; i < 4*numtiles; i++)
{
Vector3D p(tileset[i].pnt[0], tileset[i].pnt[1], tileset[i].pnt[2]);

if (imageForWorldCoords)
imageForWorldCoords->modelToWorldCoordinates(p);

fprintf(fp, "%f %f %f\n", p.getX(), p.getY(), p.getZ());
}
}

// Output the polygon lines
int offset = 1;
if (quads) {	// quads
for (figureNo = 0; figureNo < numFigures; figureNo++)
{
figure = object->getFigurePtr(figureNo);
if(figure == NULL)
continue;

for(i = 0; i < 4*tileCount[figureNo]; i += 4)
fprintf(fp, "%d %d %d -%d\n", offset + i + 3, offset + i,
offset + i + 1, offset + i + 2);
offset += 4*tileCount[figureNo];
}
}
else {	// triangles
for (figureNo = 0; figureNo < numFigures; figureNo++)
{
figure = object->getFigurePtr(figureNo);
if(figure == NULL)
continue;

for(i = 0; i < 4*tileCount[figureNo]; i += 4) {
fprintf(fp, "%d %d -%d\n", offset + i + 3, offset + i,
offset + i + 1);
fprintf(fp, "%d %d -%d\n", offset + i + 1,
offset + i + 2, offset + i + 3);
}
offset += 4*tileCount[figureNo];
}
}

delete [] pList;
delete [] tileCount;
fclose(fp);

if (globalVerbosity >= 0)
cout << "Writing BYU tiles to " << filename << endl;
}*/

/*  Output part or all of the model as BYU tiles.  The default behavior is
to output all figures of a multi-figural model into a single BYU file.
If separate is true, the figures specified by the mask (all if it is
zero), will be written to files with names consisting of the supplied
filename followed by the figure number.
*/
void P3DControl::exportToBYUTiles(const char * filename, bool quads,
								  bool useWorldCoords, int level, unsigned long figureMask, bool separate)
{
	Image3D * image;
	BYU::BYUPoly_t type;

	if(object == NULL)
		return;

	if (useWorldCoords) {
		image = imagePlanes.getImagePtr();
		if (image == NULL) {
			cout << "Error: An image must exist to export BYU tiles in world coordinates"
				<< endl;
			beep();
			return;
		}
#ifdef BINARY
		cout << "Scaling tiles to world coordinates of image" << endl;
#endif
	}
	else
		image = NULL;

	BYU byu;
	if (quads)
		type = BYU::Quads;
	else
		type = BYU::Tris;

#ifdef BINARY
	const bool discardEnds = (tuningWt(WeightMatchAtEnds) < 1e-5);
#else
	const bool discardEnds = (tuningWt(IgnoreMatchAtEnds) != 0.0);
#endif
	if ( !byu.exportModelToBYU(filename, type, level, object,
		discardEnds, figureMask, separate, image, BYU::XYZ) )
	{
		cout << "Error: BYU file not saved" << endl;
		beep();
		return;
	}

	int fname_len	= strlen(filename);
	auto_ptr<char> fname( new char[fname_len+13] );
	if(strstr(filename, ".byu") + 4 - fname_len == filename) {
		fname_len	-= 4;
	}
	memcpy(fname.get(), filename, sizeof(char)*(fname_len+1));
	memcpy(fname.get() + fname_len, "-uvt.byu", 9*sizeof(char));
	if ( !byu.exportModelToBYU(fname.get(), type, level, object,
		discardEnds, figureMask, separate, image, BYU::UVT) )
	{
		cout << "Error: BYU file (uvt data) not saved" << endl;
		beep();
		return;
	}

	memcpy(fname.get() + fname_len, "-normals.byu", 13*sizeof(char));
	if ( !byu.exportModelToBYU(fname.get(), type, level, object,
		discardEnds, figureMask, separate, image, BYU::NORMALS) )
	{
		cout << "Error: BYU file (normal data) not saved" << endl;
		beep();
		return;
	}


	if (globalVerbosity >= 0) {
		cout << "Exported BYU tiles as " << (quads ? "quads" : "triangles")
			<< " at level " << level << " to ";
		if (separate)
			cout << "files named \n  " << byu.baseName() << "NN" << byu.extension()
			<< ", where 'NN' is the figure number" << endl;
		else
			cout << filename << endl;
	}
}

void P3DControl::getRelativeMarkedPrimitiveId(int & figIndex, int & primIndex)
{
	M3DFigure * figurePtr;
	int numFigures;
	int numPrimitives;
	int i;

	primIndex = markedPrimitiveId;
	figIndex = -1;
	if(object == NULL || markedPrimitiveId == INVALID_PRIMITIVE_ID)
		return;

	numFigures = object->getFigureCount();
	for(i = 0; i < numFigures; i++)
	{
		figurePtr = object->getFigurePtr(i);
		if(figurePtr == NULL)
			continue;

		numPrimitives = figurePtr->getPrimitiveCount();
		if(primIndex >= numPrimitives)
			primIndex -= numPrimitives;

		else
		{
			figIndex = i;
			return;
		}
	}
}

void P3DControl::primitiveFiguralCoordinates(int figIndex, int primIndex, int & u, int & v)
{
	M3DFigure * figurePtr;
	int numFigures;
	int numPrimitives;

	u = -1;
	v = -1;
	if(object == NULL || figIndex < 0 || primIndex == INVALID_PRIMITIVE_ID)
		return;

	numFigures = object->getFigureCount();
	if (figIndex >= numFigures)
		return;

	figurePtr = object->getFigurePtr(figIndex);
	if (figurePtr == NULL)
		return;

	numPrimitives = figurePtr->getPrimitiveCount();
	if (primIndex >= numPrimitives)
		return;

	if(typeid(*figurePtr) == typeid(M3DQuadFigure)) {
		static_cast<M3DQuadFigure*>(figurePtr)->figuralCoordinates(primIndex, u, v);
	}
	else {
		u	= primIndex;
		v	= 0;
	}
	return;
}

void P3DControl::toggleMarkedPrimitive(int id)
{
	if(object == NULL)
		return;

	if(id == markedPrimitiveId)
		markedPrimitiveId = INVALID_PRIMITIVE_ID;
	else
		markedPrimitiveId = id;
}

M3DPrimitive * P3DControl::getMarkedPrimitive()
{
	if(object == NULL || markedPrimitiveId == INVALID_PRIMITIVE_ID)
		return NULL;

	return object->getPrimitivePtr(markedPrimitiveId);
}

void P3DControl::clearInterpolatedPrimitive() {
	if(interpolatedPrimitive != NULL)
	{
		delete interpolatedPrimitive;
		interpolatedPrimitive = NULL;
	}
}

void P3DControl::setInterpolatedPrimitive(double u, double v)
{
	int figureId, primitiveId;
	M3DFigure * figurePtr;

	clearInterpolatedPrimitive();

	if(object == NULL)
		return;

	getRelativeMarkedPrimitiveId(figureId, primitiveId);
	if (figureId < 0) {
		cout << "You must first mark an atom of the figure you want to interpolate"
			<< endl;
		beep();
		return;
	}

	figurePtr = object->getFigurePtr(figureId);
	if(figurePtr == NULL)
		return;

	M3DSubfigureTransformation transform;

	interpolatedPrimitive = new M3DQuadPrimitive;
	cout << "Interpolated primitive at u = " << u << ", v = " << v;
	// FIXME: What should we do in case of tubes? - gst/rrs
	if (transform.getInterpolatedPrimitive(*interpolatedPrimitive, dynamic_cast<M3DQuadFigure*>(figurePtr), u, v, 0)) {
		cout << ":\n";
		interpolatedPrimitive->print();
	}
	else {
		cout << " unknown" << endl;
		clearInterpolatedPrimitive();
	}
}

// Regularize the model by adjusting the atom positions
void P3DControl::regularize(double stepsize, int iterations)
{
	if (stepsize <= 0.0 || iterations < 0)
		return;

	if (object == NULL)
		return;

	startOperation();

	for( int i = 0; i != object->getFigureCount(); ++i ) {
		M3DFigure* fig	= object->getFigurePtr(i);
		if( typeid(*fig) == typeid(M3DTubeFigure)) {
			fig->resampleForRegularSpacing();
		}
		else {
			//fig->regularize(stepsize, iterations, (globalVerbosity > 0));
		}
	}

	object->setModified(true);
	modelNeedsSaving(true);

	endOperation();
	selectionChanged = true;

	if (globalVerbosity >= 0)
		cout << "Regularized model" << endl;
}

// ----------------------------- Start of Landmarks Code --------------------------------

void P3DControl::allocateImageLandmarks(int count)
{
	imageLandmarks.clear();
	for (int i = 0; i < count; i++)
		imageLandmarks.push_back(emptyImageLandmark);
	numImageLandmarks = 0;
}

void P3DControl::insertImageLandmark(int at)
{
	imageLandmarks.push_back(emptyImageLandmark);	// Add a slot at the end
	for (int i = imageLandmarks.size() - 1; i > at; i--)	// Make room in the array
		imageLandmarks[i] = imageLandmarks[i - 1];
	imageLandmarks[at] = emptyImageLandmark;	// Clear the new slot
}

void P3DControl::appendImageLandmark()
{
	imageLandmarks.push_back(emptyImageLandmark);	// Add a slot at the end
}

bool P3DControl::setImageLandmark(int landmarkId, double x, double y, double z,
								  int figureId, int number, const char * name)
{
	if (landmarkId < 0 || landmarkId >= imageLandmarks.size())
		return false;

	if (imageLandmarks[landmarkId].number < 0)	// Allows replacement or filling an empty slot
		numImageLandmarks++;
	imageLandmarks[landmarkId].coords[0] = x;
	imageLandmarks[landmarkId].coords[1] = y;
	imageLandmarks[landmarkId].coords[2] = z;
	imageLandmarks[landmarkId].figureId = figureId;
	imageLandmarks[landmarkId].number = number;
	imageLandmarks[landmarkId].name = name;
	return true;
}

imageLandmark_t P3DControl::getImageLandmark(int landmarkId) const
{
	if (landmarkId < 0 || landmarkId >= imageLandmarks.size())
		return emptyImageLandmark;
	else
		return imageLandmarks[landmarkId];
}

/*
bool P3DControl::anyLandmarkIsPaired(int figureId)
{
if (numImageLandmarks == 0)
return false;

for (int l = 0; l < imageLandmarks.size(); l++) {
if (imageLandmarks[l].figureId == figureId && imageLandmarks[l].number >= 0)
return true;
}
return false;
}
*/

bool P3DControl::landmarkIsPaired(int figureId, int landmark)
{
	if (numImageLandmarks == 0)
		return false;

	for (int l = 0; l < imageLandmarks.size(); l++) {
		if (imageLandmarks[l].figureId == figureId && imageLandmarks[l].number == landmark)
			return true;
	}

	return false;
}

/*
Vector3D P3DControl::landmarkCoordinates(int figureId, int landmark)
{
Vector3D v(-1.0, -1.0, -1.0);

if (numImageLandmarks == 0)
return v;

for (int l = 0; l < imageLandmarks.size(); l++) {
if (imageLandmarks[l].figureId == figureId && imageLandmarks[l].number == landmark) {
v.setX(imageLandmarks[l].coords[0]);
v.setY(imageLandmarks[l].coords[1]);
v.setZ(imageLandmarks[l].coords[2]);
return v;
}
}
return v;
}
*/

int P3DControl::figureOfLandmark(int index, int & landmark)
{
	if (index < imageLandmarks.size()) {
		landmark = imageLandmarks[index].number;
		return imageLandmarks[index].figureId;
	}
	else {
		landmark = -1;
		return -1;
	}
}

// Discards all image landmarks and the space to store them.
// Called when all surface landmarks are also being deleted
// or replaced (e.g. when a new model file is loaded).
void P3DControl::deleteImageLandmarks() {
#ifdef DEBUG
	cout << "P3DControl::deleteImageLandmarks()" << endl;
#endif
	imageLandmarks.clear();
	numImageLandmarks = 0;
}

// Discards all image landmarks, but retains the space for storing
// them.
void P3DControl::clearImageLandmarks() {
	for (int i = 0; i < imageLandmarks.size(); i++)
		if (imageLandmarks[i].figureId >= 0)
			imageLandmarks[i] = emptyImageLandmark;

	numImageLandmarks = 0;
}

/*
// Clears a single image landmark by setting its relative landmark
// number and coordinates to -1.  The figureId must not be corrupted
// because it will be needed if the surface landmark is then deleted.
void P3DControl::clearImageLandmark(int figureId, int landmark)
{
int i;

#ifdef DEBUG
cout << "P3DControl::clearImageLandmark(" << figureId << ", " << landmark << ')' << endl;
#endif
if (figureId < 0 || landmark < 0)
return;

for (i = 0; i < imageLandmarks.size(); i++)
if (imageLandmarks[i].figureId == figureId && imageLandmarks[i].number == landmark)
{
imageLandmarks[i] = emptyImageLandmark;
imageLandmarks[i].figureId = figureId;		// restoration

numImageLandmarks--;
callback->setLandmarksImageNumber(numImageLandmarks);

return;
}
cout << "There is no landmark to delete\n";
}
*/

// Clears a single image landmark by setting its relative landmark
// number and coordinates to -1.  The figureId must not be corrupted
// because it will be needed if the surface landmark is then deleted.
bool P3DControl::clearImageLandmark(int index)
{
	int figureId;

#ifdef DEBUG
	cout << "P3DControl::clearImageLandmark(" << index << ')' << endl;
#endif
	if (index < 0 || index >= imageLandmarks.size())
		return false;

	/*	cout << "Landmarks (before):\n";
	for (int ii = 0; ii < imageLandmarks.size(); ii++)
	cout << '(' << imageLandmarks[ii].figureId << ", " << imageLandmarks[ii].number << ")\n";
	*/
	if (imageLandmarks[index].number >= 0)
	{
		figureId = imageLandmarks[index].figureId;
		imageLandmarks[index] = emptyImageLandmark;
		imageLandmarks[index].figureId = figureId;		// restoration

		numImageLandmarks--;

		/*		cout << "NumImageLandmarks = " << numImageLandmarks << endl;
		cout << "Landmarks (after):\n";
		for (int ii = 0; ii < imageLandmarks.size(); ii++)
		cout << '(' << imageLandmarks[ii].figureId << ", " << imageLandmarks[ii].number << ")\n";
		*/
		return true;
	}

	cout << "There is no landmark to delete\n";
	return false;
}

// Discards a single image landmark and deletes the space used.
// This assumes that the surface landmark has or will be deleted.
// This assumes the number of figures in the model does not change.
void P3DControl::deleteImageLandmark(int index)
{
	int i;
	int figureId;
	bool paired;

#ifdef DEBUG
	cout << "P3DControl::deleteImageLandmark(" << index << ')' << endl;
#endif
	if (index < 0 || index >= imageLandmarks.size())
		return;

	/*	cout << "Landmarks (before):\n";
	for (i = 0; i < imageLandmarks.size(); i++)
	cout << '(' << imageLandmarks[i].figureId << ", " << imageLandmarks[i].number << ")\n";
	*/
	if (imageLandmarks[index].number >= 0)
		paired = true;
	else
		paired = false;
	figureId = imageLandmarks[index].figureId;
	// Remove the old landmark from the list
	for (i = index + 1; i < imageLandmarks.size(); i++) {
		imageLandmarks[i - 1] = imageLandmarks[i];
		// Must adjust the relative index of the landmark on the figure of deletion
		if (figureId >= 0 && imageLandmarks[i - 1].figureId == figureId)
			imageLandmarks[i - 1].number--;
	}
	imageLandmarks.pop_back();

	if (paired)
		numImageLandmarks--;

	/*	cout << "NumImageLandmarks = " << numImageLandmarks << endl;
	cout << "Landmarks (after):\n";
	for (i = 0; i < imageLandmarks.size(); i++)
	cout << '(' << imageLandmarks[i].figureId << ", " << imageLandmarks[i].number << ")\n";
	*/
}

#ifndef BINARY

void P3DControl::applyLandmarkTransformation()
{
	M3DFigure * figure;
	int numFigures;
	int numLandmarks;

	if (object == NULL)
		return;

	numLandmarks = object->getLandmarkCount();
	if (numLandmarks < 2) {
		cout << "An insufficient number of landmarks has been set" << endl;
		beep();
		return;
	}
	if (numLandmarks != numImageLandmarks) {
		cout << "All landmarks have not been paired" << endl;
		beep();
		return;
	}

	using namespace pauly;

	Matrix surfaceLandmarks(numLandmarks, 3, NULL, false);
	Matrix targetLandmarks(numLandmarks, 3, NULL, false);

	int col = 0;
	numFigures = object->getFigureCount();
	for (int figureId = 0; figureId < numFigures; figureId++) {
		figure = object->getFigurePtr(figureId);
		int n = figure->getLandmarkCount();
		for (int i = 0; i < n; i++) {
			Vector3D l = figure->getLandmark(i);
			surfaceLandmarks(col, 0) = l.getX();
			surfaceLandmarks(col, 1) = l.getY();
			surfaceLandmarks(col, 2) = l.getZ();
			if (imageLandmarks[col].coords[0] < 0)	// Landmark is not paired; this should never happen
				continue;
			targetLandmarks(col, 0) = imageLandmarks[col].coords[0];
			targetLandmarks(col, 1) = imageLandmarks[col].coords[1];
			targetLandmarks(col, 2) = imageLandmarks[col].coords[2];
			col++;
		}
	}

	LandmarkDeformation procrustes(surfaceLandmarks, targetLandmarks, col);
	//	procrustes.print_results();

	startOperation();

	selectAll();

	// Scale object
	object->scaleBy(procrustes.getBeta());

	// get rotation quat
	Matrix rotMat(procrustes.getGamma());

	double mat[4][4];

	for(int i = 0; i < 3; i++)
	{
		for(int j = 0; j < 3; j++)
		{
			mat[i][j] = rotMat(i, j);
		}
	}

	Quat qRot;
	qRot.matrixToQuat(mat);

	// Rotate object
	object->rotateBy(qRot);

	// Translate object
	Vector x = procrustes.getOffset();

	Vector3D trans(x(0), x(1), x(2));
	object->translateBy(trans);

	deselectAll();

	endOperation();

	modelNeedsSaving(true);
	selectionChanged = true;
}

#endif

// ----------------------------- End of Landmarks Code --------------------------------


// --------------------------- Start of Optimization Code -----------------------------

// Produces a vector of the IDs of figures in a single tree, a level at a time.
// The output, array treeFigures, contains the refinement order for the tree.
// This is a breadth-first version of the function below.
/*void P3DControl::listTreeFigures(M3DFigureTreeNode * node)
{
M3DFigureTreeNode * child;
M3DFigureTreeNode ** nodes;
int numFigures, childCount, i, lo, hi;

if (treeFigures != NULL)
delete [] treeFigures;
numFigures = node->countFigures();
treeFigures = new int[numFigures];

treeFiguresCount = 0;
treeFigures[0] = node->getFigureId();
if (treeFigures[0] < 0)
return;

treeFiguresCount++;
nodes = new M3DFigureTreeNode *[numFigures];	// queue
nodes[0] = node;
lo = 0;
hi = 0;
numFigures--;

// Iterative breadth-first search of the tree
while (lo <= hi) {
childCount = nodes[lo]->getChildCount();
for (i = 0; i < childCount; i++)
{
child = nodes[lo]->getChild(i);
if (child != NULL) {
treeFigures[treeFiguresCount++] = child->getFigureId();
nodes[++hi] = child;
}
}
lo++;
}
delete [] nodes;
}*/

// Produces a vector of the IDs of figures in a single tree rooted at node,
// which cannot be NULL.  The figures are determined by preorder depth-first search.
// The output, array treeFigures, contains the refinement order for the tree.
void P3DControl::listTreeFigures(M3DFigureTreeNode * root)
{
	M3DFigureTreeNode * child;
	M3DFigureTreeNode * node;

	if (treeFigures != NULL)
		delete [] treeFigures;
	treeFigures = new int[root->countFigures()];
	treeFiguresCount = 0;

	if (root->getFigureId() < 0)
		return;

	if (root->getChildCount() == 0) {
		// The tree is a single node
		treeFigures[treeFiguresCount++] = root->getFigureId();
		return;
	}

	// Iterative depth-first search of the tree
	node = root;
	while (node != NULL)
	{
		// Preorder recording of figure IDs
		treeFigures[treeFiguresCount++] = node->getFigureId();

		// Try going deeper in the tree first
		child = node->getChild(0);
		if (child == NULL) {
			// Now try going to the next sibling
			child = node->getSibling();
			while (child == NULL) {
				node = node->getParent();	// Nothing worked - go up the tree
				if (node == root) {
					child = NULL;	// Parent of root
					break;
				}
				child = node->getSibling();	// Parent already has been recorded
			}
		}
		node = child;
	}
}

// Called to start over with the first figure of the first tree
M3DFigureTreeNode * P3DControl::initializeOptimizerTree()
{
	M3DFigureTreeNode * root;

	optimizerTreeIndex = 0;	// Start with first tree of forest
	root = object->getFigureTreeRoot(optimizerTreeIndex);
	if (root == NULL) {
		cout << "Figure tree " << optimizerTreeIndex << " contains no figures\n";
		currentOptimizerStage = OPTIMIZER_END_STAGE;
		return NULL;
	}

	// Do a breadth-first search of this tree to deterimine refinement order
	listTreeFigures(root);
	if (treeFiguresCount <= 0) {
		// This should never happen
		cout << "Figure tree " << optimizerTreeIndex << " contains no figures\n";
		currentOptimizerStage = OPTIMIZER_END_STAGE;
		return NULL;
	}
	optimizerFigureId = treeFigures[0];	// Figure ID of root of tree

	// CL: reset the first PGA optimization step
	currentPGAOptimizationStep = 0;

	return root;
}

void P3DControl::selectTree(M3DFigureTreeNode * node)
{
	M3DFigureTreeNode * child;
	M3DFigure * figure;
	int childCount, i;

	i = node->getFigureId();
	figure = object->getFigurePtr(i);
	figure->select();

	childCount = node->getChildCount();
	for(i = 0; i < childCount; i++)
	{
		child = node->getChild(i);
		if (child != NULL)
			selectTree(child);
	}
}

void P3DControl::deselectTree(M3DFigureTreeNode * node)
{
	M3DFigureTreeNode * child;
	M3DFigure * figure;
	int childCount, i;

	i = node->getFigureId();
	figure = object->getFigurePtr(i);
	figure->deselect();

	childCount = node->getChildCount();
	for(i = 0; i < childCount; i++)
	{
		child = node->getChild(i);
		if (child != NULL)
			deselectTree(child);
	}
}

void P3DControl::toggleTree(M3DFigureTreeNode * node)
{
	M3DFigureTreeNode * child;
	M3DFigure * figure;
	int childCount, i;

	i = node->getFigureId();
	figure = object->getFigurePtr(i);
	figure->toggleSelected();

	childCount = node->getChildCount();
	for(i = 0; i < childCount; i++)
	{
		child = node->getChild(i);
		if (child != NULL)
			toggleTree(child);
	}
}

void P3DControl::setOptimizerIntensityWindowing(bool toggle)
{
	optimizeWithWindowing = toggle;
}

OptimizerStage P3DControl::nextOptimizerStage()
{
	int numTrees;
	M3DFigureTreeNode * figureTreeRoot;


	if (object == NULL)
		return currentOptimizerStage;

	numTrees = object->getFigureTreeCount();

	switch(currentOptimizerStage)
	{
	case NO_OPTIMIZER_STAGE:
	case OPTIMIZER_SETUP_STAGE:
		// Next object (tree) to be optimized
		figureTreeRoot = initializeOptimizerTree();
		if (figureTreeRoot == NULL)
			break;
		currentFigureIndex = 0;

		if (numTrees == 1) {
			deselectAll();
			if (treeFiguresCount > 1) {
				currentOptimizerStage = OPTIMIZER_OBJECT_STAGE;
				toggleTree(figureTreeRoot);
			}
			else {
				currentOptimizerStage = OPTIMIZER_FIGURE_STAGE;
				toggleFigure(optimizerFigureId);
			}
		}
		else {
			selectAll();
			currentOptimizerStage = OPTIMIZER_ENSEMBLE_STAGE;
		}
		currentPGAOptimizationStep = 0;  // CL: set to the first step
		break;

	case OPTIMIZER_ENSEMBLE_STAGE:
		deselectAll();
		if (treeFiguresCount > 1) {
			// Next object (tree) to be optimized
			figureTreeRoot = initializeOptimizerTree();
			if (figureTreeRoot == NULL)
				break;
			toggleTree(figureTreeRoot);
			currentOptimizerStage = OPTIMIZER_OBJECT_STAGE;
		}
		else {
			// Next figure to be optimized
#ifdef BINARY
			// This doesn't work with multiple single
			// figure because we need to select the
			// figure specified in the tuning param's
			optimizerFigureId = (int) tuningWt(BpFigureId);
#else
			optimizerFigureId = treeFigures[currentFigureIndex];
#endif
			toggleFigure(optimizerFigureId);
			currentOptimizerStage = OPTIMIZER_FIGURE_STAGE;
		}
		currentPGAOptimizationStep++; // CL: move to the next step
		break;

	case OPTIMIZER_OBJECT_STAGE:
		// Next figure to be optimized
		optimizerFigureId = treeFigures[currentFigureIndex];
		deselectAll();
		toggleFigure(optimizerFigureId);
		if (treeFiguresCount > 1)
			currentOptimizerStage = OPTIMIZER_MAIN_FIGURE_STAGE;
		else
			currentOptimizerStage = OPTIMIZER_FIGURE_STAGE;
		currentPGAOptimizationStep++; // CL: move to the next step
		break;

	case OPTIMIZER_FIGURE_STAGE:
	case OPTIMIZER_MAIN_FIGURE_STAGE:
	case OPTIMIZER_SUBFIGURE_STAGE:
		deselectAll();
		toggleFigure(optimizerFigureId);
		currentOptimizerStage = OPTIMIZER_PRIMITIVE_STAGE;
		currentPGAOptimizationStep++;
		break;

	case OPTIMIZER_PRIMITIVE_STAGE:

		deselectAll();

		// Recompute distance map of moved governor for contraints penalty
		// FIXME: I've no idea what this does - rrs
		(void) match->updateDistanceList(object, optimizerFigureId);

		toggleFigure(optimizerFigureId);
		currentOptimizerStage = OPTIMIZER_SPOKE_STAGE;
		break;

	case OPTIMIZER_SPOKE_STAGE:

		deselectAll();

		// Recompute distance map of moved governor for contraints penalty
		// FIXME: I've no idea what this does - rrs
		(void) match->updateDistanceList(object, optimizerFigureId);

		// Take next figure going depth-first through the figure tree
		currentFigureIndex++;
		if (currentFigureIndex < treeFiguresCount) {

			// Next figure (has be a subfigure) to be optimized
			optimizerFigureId = treeFigures[currentFigureIndex];
			toggleFigure(optimizerFigureId);
			currentOptimizerStage = OPTIMIZER_SUBFIGURE_STAGE;
		}
		else {
			// Finished tree; check the next tree in the forest
			optimizerTreeIndex++;
			if (optimizerTreeIndex < object->getFigureTreeCount())
			{
				// Process the next tree
				figureTreeRoot = object->getFigureTreeRoot(optimizerTreeIndex);
				if (figureTreeRoot == NULL) {
					// This should never happen
					cout << "Figure tree " << optimizerTreeIndex << " contains no figures\n";
					currentOptimizerStage = OPTIMIZER_END_STAGE;
					break;
				}
				else {
					listTreeFigures(figureTreeRoot);	// Compute new figure IDs
					if (treeFiguresCount <= 0) {
						// This should never happen
						cout << "Figure tree " << optimizerTreeIndex << " contains no figures\n";
						currentOptimizerStage = OPTIMIZER_END_STAGE;
						break;
					}
					optimizerFigureId = treeFigures[0];	// Figure ID of root of tree
					currentFigureIndex = 0;
				}

				if (treeFiguresCount > 1) {
					toggleTree(figureTreeRoot);
					currentOptimizerStage = OPTIMIZER_OBJECT_STAGE;
				}
				else {
					toggleFigure(optimizerFigureId);
					currentOptimizerStage = OPTIMIZER_FIGURE_STAGE;
				}
			}
			else {
				// Finished with all trees; reset to 1st figure of 1st tree
				figureTreeRoot = initializeOptimizerTree();
				if (figureTreeRoot == NULL)
					break;
				currentFigureIndex = 0;
				toggleFigure(optimizerFigureId);
				currentOptimizerStage = OPTIMIZER_BOUNDARY_DISPLACEMENT_STAGE;                
			}
		}
		break;


	case OPTIMIZER_BOUNDARY_DISPLACEMENT_STAGE:

		deselectAll();

		// Recompute distance map of moved governor for contraints penalty
		(void) match->updateDistanceList(object, optimizerFigureId);

		// Take next figure going depth-first through the figure tree
		currentFigureIndex++;
		if (currentFigureIndex < treeFiguresCount) {

			// Next figure (has be a subfigure) to be optimized
			optimizerFigureId = treeFigures[currentFigureIndex];
			toggleFigure(optimizerFigureId);
		}
		else {
			// Finished tree; check the next tree in the forest
			optimizerTreeIndex++;
			if (optimizerTreeIndex < object->getFigureTreeCount())
			{
				// Process the next tree
				figureTreeRoot = object->getFigureTreeRoot(optimizerTreeIndex);
				if (figureTreeRoot == NULL) {
					// This should never happen
					cout << "Figure tree " << optimizerTreeIndex << " contains no figures\n";
					currentOptimizerStage = OPTIMIZER_END_STAGE;
					break;
				}
				else {
					listTreeFigures(figureTreeRoot);	// Compute new figure IDs
					if (treeFiguresCount <= 0) {
						// This should never happen
						cout << "Figure tree " << optimizerTreeIndex << " contains no figures\n";
						currentOptimizerStage = OPTIMIZER_END_STAGE;
						break;
					}
					optimizerFigureId = treeFigures[0];	// Figure ID of root of tree
					currentFigureIndex = 0;
					toggleFigure(optimizerFigureId);
				}
			}
			else // Finished with all trees
				currentOptimizerStage = OPTIMIZER_END_STAGE;                
		}
		break;

	case OPTIMIZER_END_STAGE:
		cancelOptimization();
		break;

	default:
		currentOptimizerStage = NO_OPTIMIZER_STAGE;
		break;
	};

	return currentOptimizerStage;
}

OptimizerBase * P3DControl::getOptimzer(OptimizerBase::optimizer_t type)
{
	OptimizerBase * o;
	switch (type) {
case OptimizerBase::BoundaryDispl:
	o = voxelStageOptimizer;
	break;
case OptimizerBase::Deformation:
	o = deformationOptimizer;
	break;
case OptimizerBase::Spoke:
	o = spokeOptimizer;
	break;
case OptimizerBase::MainFigure:
	o = mainFigureOptimizer;
	break;
case OptimizerBase::Registration:
	o = registrationOptimizer;
	break;
case OptimizerBase::Subfigure:
	o = subfigureOptimizer;
	break;
	}
	return o;
}

void P3DControl::cancelOptimization()
{
#ifdef OPTIMIZATION_VISUALIZER
	globalLogManager.setOptimizer(NULL);
#endif

	if (registrationOptimizer != NULL) {
		delete registrationOptimizer;
		registrationOptimizer = NULL;
	}
	if (mainFigureOptimizer != NULL) {
		delete mainFigureOptimizer;
		mainFigureOptimizer = NULL;
	}

	if (subfigureOptimizer != NULL) {
		delete subfigureOptimizer;
		subfigureOptimizer = NULL;
	}
	if (deformationOptimizer != NULL) {
		delete deformationOptimizer;
		deformationOptimizer = NULL;
	}
	if (spokeOptimizer != NULL) {
		delete spokeOptimizer;
		spokeOptimizer = NULL;
	}
	if (voxelStageOptimizer != NULL) {
		delete voxelStageOptimizer;
		voxelStageOptimizer = NULL;
	}

	if (recordedTransformation)
		recordedTransformation->setToIdentity();

#ifndef BINARY
	handPlacementTransform->setToIdentity();
#endif
	currentOptimizerStage = NO_OPTIMIZER_STAGE;
}

bool P3DControl::initializeMatch(MatchType type, int surfaceLevel,
								 char * templateFilename, const char * profileFilename, Image3D * trainingImage)
{
	Image3D * image;

#ifdef DEBUG
	cout << "P3DControl::initializeMatch()\n";
#endif

	if (object == NULL)
		return false;

	image = imagePlanes.getImagePtr();
	if (image == NULL)
		return false;

#ifdef BINARY
	if (match == NULL) {
		match = new Match(landmarkObject);

		// Dibyendu
		//xiaojie 0428

		// set the reference object after scalind the M3DObject inside the Match class

		if ( match != NULL ) {
			M3DObject* scaledObject = object->assign() ;
			match->setRefObjectAfterScale( scaledObject ) ;
		}
	}

	match->definePLists(object, surfaceLevel);
#else
	if (match == NULL)
		match = new Match(isAdaptive());

	match->definePLists(object);

	if (spe != NULL)
		match->setHistogram(spe);
	if (dqf != NULL)
		match->setDQFMatch(dqf);
	if (mp != NULL) {
		mp->loggingFileName = templateFilename;
		match->setMultiPatch(mp);
	}
#endif

	if (trainingImage != NULL)
		trainingImage->intensityWindow(image->getMinIntensityWindow(),
		image->getMaxIntensityWindow());

	// This is to compute the inter-figural distance penalty
	object->invertConstraints();

#ifdef BINARY
	if (! match->initialize(object, image, dist_map, true, trainingImage,			// AGG: This can pass a NULL mistakenly!
		optimizeWithWindowing, type, surfaceLevel, templateFilename, profileFilename))
		return false;
	dist_map = NULL;	// See setDistMapPtr()
#else
	if (! match->initialize(object, image, true, trainingImage, optimizeWithWindowing,
		type, surfaceLevel, templateFilename, profileFilename))
		return false;
#endif

	setInitialTransformInMatch(objectSimTransform);

	if (matchReferenceObject)
		match->resetReferenceObject(matchReferenceObject);

	if (! match->initializeDistanceLists(object, surfaceLevel))
		return false;
	if (! match->initializeConstraintPointsLists(object, true))
		return false;

	return true;
}

#ifndef BINARY

void P3DControl::initializeRegistrationOptimizer()
{
#ifdef DEBUG
	cout << "P3DControl::initializeRegistrationOptimizer()\n";
#endif
	if (match == NULL || object == NULL)
		return;

	if (registrationOptimizer != NULL)
		delete registrationOptimizer;

	if (pgaStats && isAdaptive()) {
		/* CL: 08/2004

		For adaptive pablo.  The differences between this version and the
		"regular" version are:			
		- there is no similarity component in the optimizer

		JJ: 05/2005
		- currentPGAOptimizationStep == 0 at this stage
		*/
		registrationOptimizer = new M3DAdaptiveRegistrationPGAOptimizer(pgaStats,
			currentPGAOptimizationStep);
#ifdef DEBUG
		cout << "\tM3DAdaptiveRegistrationPGAOptimizer() called" << endl;
		cout << "\tcurrentPGAOptimizationStep: " << currentPGAOptimizationStep << endl;
#endif
	}
	else
		registrationOptimizer = new M3DRegistrationOptimizer;

	selectionChanged = false;

	setInitialTransformInMatch(objectSimTransform);
	if (isAdaptive())
		registrationOptimizer->setPenaltyWeight(registrationPenaltyWeight);
	registrationOptimizer->initialize(match);

	// Make sure the atoms selected by the user are the ones used for optimization
	if (matchReferenceObject != NULL)	// Used in the bootstrap
		match->resetReferenceObject(matchReferenceObject);
	else
		match->resetReferenceObject(object);

	/*JVS:	In the adaptive version, there has been no handplacement, but there may have been 
	a bone transform applied to object.  In this case, what object do we want to use in
	createTargetObject?  Problem is, if I use the PGA original Mean (in m3dAdapRegPGAProb::
	createTargetObject), then it doesn't have the bone transform applied yet.  The issue is
	when is the bone transform applied. 
	JJ:	05/2005 Applied to image right now.
	*/

	if (! isAdaptive())
		if (objectSimTransform)
			objectSimTransform->setToIdentity();

#ifdef OPTIMIZATION_VISUALIZER
	globalLogManager.setOptimizer(NULL);
#endif
}

// Object ensemble (model) stage of the optimizer
void P3DControl::doRegistrationIterations(int numIterations)
{
	double initVal;

	if (registrationOptimizer == NULL || object == NULL)
		return;

	initVal = registrationOptimizer->getLastBestVal();

	registrationOptimizer->performIterations(numIterations);

	if (registrationOptimizer->getLastBestVal() < initVal) {
		if (isAdaptive()) {
			startOperation();

			if (object != NULL)
				delete object;

			object = registrationOptimizer->getLastBestObject()->assign();

#ifdef DEBUG
			cout << "\t(Adaptive) P3DControl::doRegistrationIterations()" << endl;
			Vector x = registrationOptimizer->getLastBestX();
			int vecSize = x.size();
			for (int i = 0; i < vecSize; i++)
				cout << x(i)*M3DRegistrationProblem::ENSEMBLE_PGA_SCALE_FACTOR << "  ";
			cout<<endl;
#endif

		}
		else {
			M3DObject *savedObject;
			SimilarityTransform3D transformation;

			savedObject = registrationOptimizer->getReferenceObject();
			if (object == NULL || savedObject == NULL)
				return;

			startOperation();

			(*object) = (*savedObject);

			if (objectSimTransform)
				transformation = *objectSimTransform;
			transformation.multiplyBy(registrationOptimizer->getLastBestTransformation());
#ifdef UNFLIPPED
			object->applySimilarity(transformation, Vector3D(0.0, 1.0, 0.0));
#else
			object->applySimilarity(transformation, Vector3D(0.0, 0.0, 0.0));
#endif
		}

		object->setModified(true);
		modelNeedsSaving(true);

		endOperation();
	}
}

#endif	/* BINARY */

double P3DControl::getLastRegistrationMatchValue()
{
	if (registrationOptimizer == NULL)
		return 0.0;

	return registrationOptimizer->getLastBestVal();
}

double P3DControl::getLastRegistrationPenaltyValue()
{
	if (registrationOptimizer == NULL)
		return 0.0;

	return registrationOptimizer->getLastPenalty();
}

#ifndef BINARY

void P3DControl::setRegistrationPenaltyWeight(double w)
{
	if (registrationOptimizer != NULL)
		registrationOptimizer->setPenaltyWeight(w);
	if (isAdaptive())
		registrationPenaltyWeight = w;
}

#endif	/* BINARY */

bool P3DControl::registrationOptimizerIsFinished()
{
	if (registrationOptimizer == NULL)
		return false;

	return registrationOptimizer->isFinished();
}

#ifdef BINARY

bool P3DControl::initializeMainFigureOptimizer(int doScale, const char *pgaFilename)
{
#ifdef DEBUG
	cout << "P3DControl::initializeMainFigureOptimizer()\n";
#endif
	if (! match || ! object)
		return false;

	if (mainFigureOptimizer != NULL)
		delete mainFigureOptimizer;
	
	//bool havePGA = pgaFilename && *pgaFilename;  // if we get here, the model must have PGA stat's
	bool wantPGA = (tuningWt(FigureMahalanobisMatch) != 0.0);

	if (wantPGA /*&& havePGA*/) 
	{	// havePGA was tested by runBinaryPablo()
	
		cpnsStats = object->getCPNSStats();
		if (pgaStats == NULL && cpnsStats == NULL) {
			// This should never happen
			cout << "Internal Error: No PGA/CPNS statistics" << endl;
			return false;
		}
		if (cpnsStats != NULL)
		{
			cout << "Using CPNS Mahalanobis Stage" << endl;
			mainFigureOptimizer = new M3DSimilarityCPNSOptimizer(match, cpnsStats);
		}
		else
		{
			cout << "Using PGA Mahalanobis Stage" << endl;
			mainFigureOptimizer = new M3DSimilarityPGAOptimizer(match, pgaStats);
		}
	}
	else {
		mainFigureOptimizer =
			new M3DSimilarityElongationOptimizer(match, (doScale ? 7 : 6));
		match->setInitialTransform(object->getTransformation());	// probably unity
	}
	selectionChanged = false;

	setInitialTransformInMatch(objectSimTransform);

	// MultiObject -- the optimizerTreeIndex isn't used
	if (! mainFigureOptimizer->initialize(object, optimizerTreeIndex))
	{
		return false;
	}
	/*
	THIS PROBABLY IS NOT RIGHT FOR SIMILARITY - USE LOADED OBJECT INSTEAD?
	if(matchReferenceObject)
	match->resetReferenceObject(matchReferenceObject);
	*/
	return true;
}

#else	/* BINARY */

bool P3DControl::initializeMainFigureOptimizer(bool rigid, bool usePGA)
{
#ifdef DEBUG
	cout << "P3DControl::initializeMainFigureOptimizer()\n";
#endif
	if (match == NULL || object == NULL)
		return false;

	if (mainFigureOptimizer != NULL)
		delete mainFigureOptimizer;

	if (usePGA) {
		/* CL: HACK: if currentPGAOptimizationStep == 0, use similarity + PGA (Single figure object)
		if currentPGAOptimizationStep >= 1, use residue PGA  (Multi-object of single figure) 
		*/
		if (currentPGAOptimizationStep == 0) {	//Single figure, PGA case, regular Pablo
			if (pgaStats == NULL) {
				// This should never happen
				cout << "Internal Error: No PGA statistics" << endl;
				return false;
			}

			mainFigureOptimizer = new M3DSimilarityPGAOptimizer(match, pgaStats);

			//JJ : Where to reset HandPlacement?
			pgaStats->setHandPlacementTransform(*handPlacementTransform);
			handPlacementTransform->setToIdentity();

#ifdef DEBUG
			cout << "\tM3DSimilarityPGAOptimizer() called" << endl;
#endif
		}
		else if (currentPGAOptimizationStep >= 1) {
			if (pgaStats == NULL) {
				// This should never happen
				cout << "Internal Error: No PGA statistics" << endl;
				return false;
			}
			// Multifigure PGA case, adaptive Pablo
			mainFigureOptimizer = new M3DMainFigResiduePGAOptimizer(match, pgaStats, 
				currentPGAOptimizationStep);

#ifdef DEBUG
			cout << "\tM3DMainFigResiduePGAOptimizer() called" << endl;
#endif
		}
		else
			cerr << "	Error: unknown optimization step." << endl;
#ifdef DEBUG
		cout << "\tcurrentPGAOptimizationStep: " << currentPGAOptimizationStep << endl;
#endif
	}
	else {
		// When rigid, don't optimize scale or elongation
		mainFigureOptimizer = new M3DSimilarityElongationOptimizer(match, (rigid ? 6 : 8));

#ifdef DEBUG
		cout << "\tM3DSimilarityElongationOptimizer() called" << endl;
#endif
	}
	selectionChanged = false;

	setInitialTransformInMatch(objectSimTransform);

	mainFigureOptimizer->setPenaltyWeightsInMatch((int *) penaltyNames, penaltyWeights,
		numPenaltyWeights);

	//XX :BUG FIX?  can not do multi-figure serially
	if (usePGA) {
		// Write the object w/o similarity applied later on.  In other words, record the 
		// optimized similarity in the transform of the m3d.
		if (currentPGAOptimizationStep >= 1) {	// "M3DMainFigResiduePGAOptimizer" is called
			M3DObject * referenceObj = pgaStats->getLatestObjPtr()->assign();
			// JVS: In this case the objectSimTransform, as currently constituted, was applied
			// during the object ensemble stage
			if (! mainFigureOptimizer->initialize(referenceObj, optimizerTreeIndex))
				return false;
		}
		else {	// "M3DSimilarityPGAOptimizer" is called
			if (pgaModelFileName != NULL) {
				cout << "Warning: Considering model " << pgaModelFileName
					<< " to be a residue model" << endl;
				M3DObject * meanResidue	= pgaStats->getLatestObjPtr()->assign();
				if (! mainFigureOptimizer->setAddToResidueObject(object->assign()) )
					return false;
				if (! mainFigureOptimizer->initialize(meanResidue, optimizerTreeIndex))
					return false;
			}
			else {
				if (! mainFigureOptimizer->initialize(object, optimizerTreeIndex))
					return false;
			}
			//
			// Rohit - get a copy of the object back from the optimizer because starting
			// coefficients might have caused the object to change from what we have
			// *ONLY* if M3DSimilarityPGAOptimizer was called.
			//
			if (object != NULL)
				delete object;

			object = mainFigureOptimizer->getLastBestObject()->assign();
		}
	}
	else	// "M3DSimilarityElongationOptimizer" is called
		if (! mainFigureOptimizer->initialize(object, optimizerTreeIndex))
			return false;
	// JVS: In the case that mainFigResPGAOpt is the mainfigureoptimizer, object may be 
	// the wrong argument, since it has been modified by similarities.  Need a model in 
	// the PGA coordinates.  Basically, if usePGA is true, then the reference object should
	// be without similarities applied to it.

	if (matchReferenceObject) 
		match->resetReferenceObject(matchReferenceObject);

#ifdef OPTIMIZATION_VISUALIZER
	globalLogManager.setOptimizer(mainFigureOptimizer);
#endif
	return true;
}

#endif	/* BINARY */

// Main figure stage of the optimizer
void P3DControl::doMainFigureIterations(int numIterations, bool verbose)
{
	int treeIndex;
	int figureId, len;
	M3DFigureTreeNode * treeRoot;
	double initVal;

	if (object == NULL)
		return;

	if (mainFigureOptimizer == NULL)
		return;

	// See if the current figure is to be skipped
	const int * skipFigs = tuningList(FigureSkipFigures, len);

	for (int i = 0; i < len; i++)
		if (skipFigs[i] == optimizerFigureId) {
			cout << "Figure number " << optimizerFigureId;
			M3DFigure * figure = object->getFigurePtr(optimizerFigureId);
			const char * name = figure->getName();
			if (name != NULL)
				cout << " (" << name << ')';
			cout << " skipped" << endl;
			mainFigureOptimizer->skipStage();	// Indicate skipping figure
			return;
		}

		initVal = mainFigureOptimizer->getLastBestVal();

		mainFigureOptimizer->performIterations(numIterations);

		treeIndex = mainFigureOptimizer->getTreeIndex();
		treeRoot = object->getFigureTreeRoot(treeIndex);
		if (treeRoot == NULL)
			return;

		figureId = treeRoot->getFigureId();

		if (verbose) {
			cout << "Figure number " << figureId;
			M3DFigure * figure = object->getFigurePtr(figureId);
			const char * name = figure->getName();
			if (name != NULL)
				cout << " (" << name << ")\n";
			else
				cout << '\n';
		}
#ifdef DEBUG
		cout << "figureId in P3DControl::doMainFigureIterations() = " << figureId << endl;
#endif

		//if (mainFigureOptimizer->getLastBestVal() < initVal) {
		// LEVY, SABOO 9-20-2007 Always do this, because PGA optimization
		// with a fixed simtrans, needs to apply the simtrans even if 
		// deformation doesn't improve the objective function
		// (or if you pass in a starting PGA vector, that is optimal,
		// it should still be applied)
		{
			startOperation();

			if (object != NULL)
				delete object;

			object = mainFigureOptimizer->getLastBestObject()->assign();

#ifdef DEBUG
			cout << "\tP3DControl::doMainFigureIterations()" << endl;
			Vector x = mainFigureOptimizer->getLastBestX();
			int vecSize = x.size();
			for (int j = 0; j < vecSize; j++)
				cout << x(j)*tuningWt(FigurePGAFactor) << '\t';
			cout << endl;
#endif

			M3DFigure * figure = 
				mainFigureOptimizer->getReferenceObject()->getFigurePtr(figureId);

			M3DFigureTreeNode * figureTreeNode = object->getFigureTreeNode(figureId);

			int i = 0;
			M3DFigureTreeNode * subFigNode = figureTreeNode->getChild(i);
			while (subFigNode != NULL) {
				subfigureTransform.init(object, subFigNode);
				// FIXME: What should we do for tubes? - gst/rrs
				subfigureTransform.updateSubfigure(dynamic_cast<M3DQuadFigure*>(figure));
				i++;
				subFigNode = figureTreeNode->getChild(i);
			}

			object->setModified(true);
			modelNeedsSaving(true);

			endOperation();
		}
}

bool P3DControl::mainFigureStageSkipped()
{
	if (mainFigureOptimizer == NULL)
		return true;

	return mainFigureOptimizer->stageSkipped();
}

double P3DControl::getLastMainFigureObjectiveFunctionValue()
{
	if (mainFigureOptimizer == NULL)
		return 0.0;

	return mainFigureOptimizer->getLastBestVal();
}

double P3DControl::getLastMainFigureGeometricMatchValue()
{
	if (mainFigureOptimizer == NULL)
		return 0.0;

	return mainFigureOptimizer->getLastPenalty();
}

#ifndef BINARY

// Discard penalty weights (used at the start of optimizer stages)
void P3DControl::clearPenaltyWeights() {
	if (penaltyWeights != NULL) {
		delete [] penaltyWeights;
		penaltyWeights = NULL;
		delete [] penaltyNames;
		penaltyNames = NULL;
	}
	numPenaltyWeights = 0;
}

// This function cannot be used to change a previously set weight
void P3DControl::setMainFigurePenaltyWeight(gpTune_t tuningParam, double w)
{
	if (penaltyWeights == NULL) {
		penaltyWeights = new double[MAX_NUM_FIG_MATCH_RESULTS];
		penaltyNames = new gpTune_t[MAX_NUM_FIG_MATCH_RESULTS];
	}
	if (numPenaltyWeights >= MAX_NUM_FIG_MATCH_RESULTS) {
		cout << "Error: too many penalty weights\n";
		return;
	}
	penaltyWeights[numPenaltyWeights] = w;
	penaltyNames[numPenaltyWeights++] = tuningParam;
}

#endif	/* BINARY */

bool P3DControl::mainFigureOptimizerIsFinished()
{
	if (mainFigureOptimizer == NULL)
		return false;

	return mainFigureOptimizer->isFinished();
}

#ifndef BINARY

void P3DControl::initializeSubfigureOptimizer()
{
	if (match == NULL || object == NULL)
		return;

	if (subfigureOptimizer != NULL)
		delete subfigureOptimizer;

	subfigureOptimizer = new M3DSubfigureOptimizer;
	selectionChanged = false;

	setInitialTransformInMatch(objectSimTransform);

	subfigureOptimizer->initialize(match, object, optimizerFigureId);

	if(matchReferenceObject)
		match->resetReferenceObject(matchReferenceObject);

#ifdef OPTIMIZATION_VISUALIZER
	globalLogManager.setOptimizer(NULL);
#endif
}

// Subfigure stage of the optimizer
void P3DControl::doSubfigureIterations(int numIterations)
{
	M3DFigureTreeNode * treeNode;
	double initVal;
	Vector x;
	M3DObject * savedObject;
	M3DSubfigureTransformation transform;


	if (object == NULL)
		return;

	if (subfigureOptimizer == NULL)
		return;

	initVal = subfigureOptimizer->getLastBestVal();

	subfigureOptimizer->performIterations(numIterations);

	if (subfigureOptimizer->getLastBestVal() < initVal) {
		x = subfigureOptimizer->getLastBestX();

#ifdef DEBUG
		// Debugging information: to print the optimized parameters at this stage
		cout << "------------------------------------" << endl;
		for(int dex=0; dex<x.size(); dex++)
			if(dex>=3 && dex<=4)
				cout << pow(LOG_BASE, x(dex)) << ", ";
			else
				cout << x(dex) << ", ";
		cout << "\n------------------------------------" << endl;
#endif

		savedObject = subfigureOptimizer->getReferenceObject();
		if (savedObject == NULL)
			return;

		startOperation();

		if (object != NULL)
			delete object;

		object = savedObject->assign();

		treeNode = object->getFigureTreeNode(optimizerFigureId);
		if (treeNode == NULL)
			return;
		transform.init(object, treeNode);

		transform.translate(M3DSubfigureProblem::SUB_TRANSLATION_FACTOR*x(0), \
			M3DSubfigureProblem::SUB_TRANSLATION_FACTOR*x(1));
		transform.rotate(M3DSubfigureProblem::SUB_ROTATION_FACTOR*x(2));
		transform.widen(pow(LOG_BASE, M3DSubfigureProblem::SUB_WIDEN_FACTOR*x(3)));			//(exp(x(3)));
		transform.elongate(pow(LOG_BASE, M3DSubfigureProblem::SUB_ELONGATE_FACTOR*x(4)));	//(exp(x(4)));
		// no hinge rotation at this stage, yet

		object->setModified(true);
		modelNeedsSaving(true);

		endOperation();
	}
}

double P3DControl::getLastSubfigureMatchValue()
{
	if (subfigureOptimizer == NULL)
		return 0.0;

	return subfigureOptimizer->getLastBestVal();
}

bool P3DControl::subfigureOptimizerIsFinished()
{
	if (subfigureOptimizer == NULL)
		return false;

	return subfigureOptimizer->isFinished();
}

#endif	/* BINARY */

/**
* Atom stage deformation optimizer stuff.
*/

bool P3DControl::initializeDeformationOptimizer(int appID)
{
	if (match == NULL) {
		cout << "Error: No Match class!\n";
		return false;
	}

	deformationObjFunc = 0.0;
	deformationGeomPenalty = 0.0;
	if (deformationOptimizer == NULL)
		deformationOptimizer = new M3DDeformationOptimizer(appID);
	selectionChanged = false;

#ifndef BINARY
	if (! deformationOptimizer->initialize(match, object)) {
		cout << "Error: Unable to initialize deformation optimizer\n";
		return false;
	}

	deformationOptimizer->setPenaltyWeightsInMatch((int *) penaltyNames, penaltyWeights,
		numPenaltyWeights);
#else

	// Atom PGA 
	bool wantPGA = (tuningWt(AtomMahalanobisMatch) != 0.0);
	if (wantPGA) {
		if (atomPgaStats == NULL) {
			cout << "Error: no atom PGA data available; ignoring AtomMahalanobisMatch"
				<< endl;
			wantPGA = false;
		}
	}

	if (wantPGA) {			
		deformationOptimizer->initialize(match, object);
		int figureId = (int) tuningWt(BpFigureId);//multiobj??xiaoxiao
		deformationOptimizer->initializeAtomPGAStage(atomPgaStats,figureId);

	}
	else {
		deformationOptimizer->initialize(match, object);
	}
#endif

	if (matchReferenceObject)
		match->resetReferenceObject(matchReferenceObject);

#ifdef OPTIMIZATION_VISUALIZER
	globalLogManager.setOptimizer(NULL);
#endif

	return true;
}

void P3DControl::resetDeformationOptimizer()
{
	if (match == NULL)
		return;

	deformationObjFunc = 0.0;
	deformationGeomPenalty = 0.0;
	selectionChanged = false;

	if (deformationOptimizer == NULL)
		deformationOptimizer = new M3DDeformationOptimizer;
	deformationOptimizer->reset(match, object);

	if (matchReferenceObject)
		match->resetReferenceObject(matchReferenceObject);
}

// Atom stage of the optimizer
bool P3DControl::doDeformationIterations(int figureId, bool preview, bool verbosity)
{
	Image3D * image;
	bool ret;

	if (object == NULL)
		return false;

	image = imagePlanes.getImagePtr();
	if (image == NULL)
		return false;

	if(deformationOptimizer == NULL) {
		cout << "NULL optimizer\n";
		return false;
	}
#ifndef BINARY
	if (selectionChanged) {
		/* Note: it is not good enough just to call resetDeformationOptimizer()
		here, since the candidate and reference objects of the deformation
		optimizer are available to users of this class.  See for example
		P3DUserInterfaceCallback::startPrimitiveStage().
		*/
		cout << "The optimizer must be reset" << endl;
		return false;
	}
#endif

	if (preview) {
		return deformationOptimizer->performIterations(figureId, true, false);
	}

	startOperation();

	ret = deformationOptimizer->performIterations(figureId, false, verbosity);
	if (ret) {
		/* This code allows the deformation to be repeated by resetting the variables
		associated with the skipping of atoms in the deformation.  However, it
		also keeps track of the results of each fit and when they repeat, it stops
		resetting the variables, resulting in the next call to performIterations()
		simply returning true, indicating that refinement is complete.
		*/
		if (deformationObjFunc != deformationOptimizer->getLastBestVal()
			|| deformationGeomPenalty != deformationOptimizer->getLastPenalty())
			deformationOptimizer->readyNextRun();
		deformationObjFunc = deformationOptimizer->getLastBestVal();
		deformationGeomPenalty = deformationOptimizer->getLastPenalty();
	}

	if (deformationOptimizer->getCandidateObject() != NULL)
	{
		(*object) = (*(deformationOptimizer->getCandidateObject()));

		M3DFigure * figure = deformationOptimizer->getReferenceObject()->getFigurePtr(figureId);

		M3DFigureTreeNode * figureTreeNode = object->getFigureTreeNode(figureId);

		int i = 0;
		M3DFigureTreeNode * subFigNode = figureTreeNode->getChild(i);
		while (subFigNode != NULL)
		{
			subfigureTransform.init(object, subFigNode);
			// FIXME: What should we do in case of tubes? - gst/rrs
			subfigureTransform.updateSubfigure(dynamic_cast<M3DQuadFigure*>(figure));
			i++;
			subFigNode = figureTreeNode->getChild(i);
		}
	}

	object->setModified(true);
	modelNeedsSaving(true);
	endOperation();
	return ret;
}

double P3DControl::getLastDeformationObjectiveFunctionValue(bool preview)
{
	if (deformationOptimizer == NULL)
		return 0.0;

	return deformationOptimizer->getLastBestVal(preview);
}

double P3DControl::getLastDeformationGeometricMatchValue()
{
	if (deformationOptimizer == NULL)
		return 0.0;

	return deformationOptimizer->getLastPenalty();
}


#ifndef BINARY


// This function cannot be used to change a previously set weight
void P3DControl::setDeformationPenaltyWeight(gpTune_t tuningParam, double w)
{
	if (penaltyWeights == NULL) {
		penaltyWeights = new double[MAX_NUM_ATOM_MATCH_RESULTS];
		penaltyNames = new gpTune_t[MAX_NUM_ATOM_MATCH_RESULTS];
	}
	if (numPenaltyWeights >= MAX_NUM_ATOM_MATCH_RESULTS) {
		cout << "Error: too many penalty weights\n";
		return;
	}
	penaltyWeights[numPenaltyWeights] = w;
	penaltyNames[numPenaltyWeights++] = tuningParam;
}

#endif	/* BINARY */

// dibyendu

/**
* SRep stage deformation optimizer stuff.
*/

bool P3DControl::initializeSRepOptimizer(int appID)
{
	if (match == NULL) {
		cout << "Error: No Match class!\n";
		return false;
	}

	sRepObjFunc = 0.0;
	sRepGeomPenalty = 0.0;
	if (sRepOptimizer == NULL)
		sRepOptimizer = new M3DSRepOptimizer(appID);
	selectionChanged = false;

#ifndef BINARY
	if (! sRepOptimizer->initialize(match, object)) {
		cout << "Error: Unable to initialize SRep optimizer\n";
		return false;
	}

	sRepOptimizer->setPenaltyWeightsInMatch((int *) penaltyNames, penaltyWeights,
		numPenaltyWeights);
#else

	// SRep PGA 
	bool wantPGA = false ;

	sRepOptimizer->initialize(match, object);

#endif

	if (matchReferenceObject)
		match->resetReferenceObject(matchReferenceObject);

#ifdef OPTIMIZATION_VISUALIZER
	globalLogManager.setOptimizer(NULL);
#endif

	return true;
}

void P3DControl::resetSRepOptimizer()
{
	if (match == NULL)
		return;

	sRepObjFunc = 0.0;
	sRepGeomPenalty = 0.0;
	selectionChanged = false;

	if (sRepOptimizer == NULL)
		sRepOptimizer = new M3DSRepOptimizer;
	sRepOptimizer->reset(match, object);

	if (matchReferenceObject)
		match->resetReferenceObject(matchReferenceObject);
}

// SRep stage of the optimizer
bool P3DControl::doSRepIterations(int figureId, bool preview, bool verbosity)
{
	Image3D * image;
	bool ret;

	if (object == NULL)
		return false;

	image = imagePlanes.getImagePtr();
	if (image == NULL)
		return false;

	if(sRepOptimizer == NULL) {
		cout << "NULL optimizer\n";
		return false;
	}
#ifndef BINARY
	if (selectionChanged) {
		/* Note: it is not good enough just to call resetsRepOptimizer()
		here, since the candidate and reference objects of the sRep
		optimizer are available to users of this class.  See for example
		P3DUserInterfaceCallback::startPrimitiveStage().
		*/
		cout << "The optimizer must be reset" << endl;
		return false;
	}
#endif

	if (preview) {
		return sRepOptimizer->performIterations(figureId, true, false);
	}

	startOperation();

	ret = sRepOptimizer->performIterations(figureId, false, verbosity);
	if (ret) {
		/* This code allows the deformation to be repeated by resetting the variables
		associated with the skipping of atoms in the deformation.  However, it
		also keeps track of the results of each fit and when they repeat, it stops
		resetting the variables, resulting in the next call to performIterations()
		simply returning true, indicating that refinement is complete.
		*/
		if (sRepObjFunc != sRepOptimizer->getLastBestVal()
			|| sRepGeomPenalty != sRepOptimizer->getLastPenalty())
			sRepOptimizer->readyNextRun();
		sRepObjFunc = sRepOptimizer->getLastBestVal();
		sRepGeomPenalty = sRepOptimizer->getLastPenalty();
	}

	if (sRepOptimizer->getCandidateObject() != NULL)
	{
		(*object) = (*(sRepOptimizer->getCandidateObject()));

		M3DFigure * figure = sRepOptimizer->getReferenceObject()->getFigurePtr(figureId);

		M3DFigureTreeNode * figureTreeNode = object->getFigureTreeNode(figureId);

		int i = 0;
		M3DFigureTreeNode * subFigNode = figureTreeNode->getChild(i);
		while (subFigNode != NULL)
		{
			subfigureTransform.init(object, subFigNode);
			// FIXME: What should we do in case of tubes? - gst/rrs
			subfigureTransform.updateSubfigure(dynamic_cast<M3DQuadFigure*>(figure));
			i++;
			subFigNode = figureTreeNode->getChild(i);
		}
	}

	object->setModified(true);
	modelNeedsSaving(true);
	endOperation();
	return ret;
}

double P3DControl::getLastSRepObjectiveFunctionValue(bool preview)
{
	if (sRepOptimizer == NULL)
		return 0.0;

	return sRepOptimizer->getLastBestVal(preview);
}

double P3DControl::getLastSRepGeometricMatchValue()
{
	if (sRepOptimizer == NULL)
		return 0.0;

	return sRepOptimizer->getLastPenalty();
}





/**
* Spoke stage deformation optimizer stuff.
*/

bool P3DControl::initializeSpokeOptimizer(int appID)
{
	if (match == NULL) {
		cout << "Error: No Match class!\n";
		return false;
	}

	spokeObjFunc = 0.0;
	spokeGeomPenalty = 0.0;
	if( spokeOptimizer == NULL ) {
		spokeOptimizer	= new M3DSpokeOptimizer(appID);
	}
	selectionChanged = false;

#ifndef BINARY
	if (! spokeOptimizer->initialize(match, object)) {
		cout << "Error: Unable to initialize spoke optimizer\n";
		return false;
	}

	spokeOptimizer->setPenaltyWeightsInMatch((int *) penaltyNames, penaltyWeights,
		numPenaltyWeights);
#else

	// Spoke PGA 
	bool wantPGA = (tuningWt(SpokeMahalanobisMatch) != 0.0);
	if (wantPGA) {
		if (spokePgaStats == NULL) {
			cout << "Error: no spoke PGA data available; ignoring SpokeMahalanobisMatch"
				<< endl;
			wantPGA = false;
		}
	}

	/*
	if (wantPGA) {
	M3DObject * meanObj;

	// The meanObj is used for alignment before applying the atom difference PGA
	meanObj = NULL;
	if (pgaStats != NULL)
	meanObj = pgaStats->getMeanObjPtr();

	spokeOptimizer->initialize(match, object);
	int figureId = (int) tuningWt(BpFigureId);//multiobj??xiaoxiao
	if (meanObj != NULL)
	spokeOptimizer->initializeSRepPGAStage(SRepPgaStats, meanObj,figureId);

	}
	else
	*/
	{
		spokeOptimizer->initialize(match, object);
	}
#endif

	if (matchReferenceObject)
		match->resetReferenceObject(matchReferenceObject);

#ifdef OPTIMIZATION_VISUALIZER
	globalLogManager.setOptimizer(NULL);
#endif

	return true;
}

void P3DControl::resetSpokeOptimizer()
{
	if (match == NULL)
		return;

	spokeObjFunc = 0.0;
	spokeGeomPenalty = 0.0;
	selectionChanged = false;

	if (spokeOptimizer == NULL)
		spokeOptimizer = new M3DSpokeOptimizer;
	spokeOptimizer->reset(match, object);

	if (matchReferenceObject)
		match->resetReferenceObject(matchReferenceObject);
}

// Spoke stage of the optimizer
bool P3DControl::doSpokeIterations(int figureId, bool preview, bool verbosity)
{
	Image3D * image;
	bool ret;

	if (object == NULL)
		return false;

	image = imagePlanes.getImagePtr();
	if (image == NULL)
		return false;

	if( spokeOptimizer == NULL ) {
		cout << "NULL optimizer\n";
		return false;
	}
#ifndef BINARY
	if (selectionChanged) {
		/* Note: it is not good enough just to call resetSpokeOptimizer()
		here, since the candidate and reference objects of the spoke
		optimizer are available to users of this class.  See for example
		P3DUserInterfaceCallback::startPrimitiveStage().
		*/
		cout << "The optimizer must be reset" << endl;
		return false;
	}
#endif

	if (preview) {
		return spokeOptimizer->performIterations(figureId, true, false);
	}

	startOperation();

	ret = spokeOptimizer->performIterations(figureId, false, verbosity);
	if (ret) {
		/* This code allows the deformation to be repeated by resetting the variables
		associated with the skipping of atoms in the deformation.  However, it
		also keeps track of the results of each fit and when they repeat, it stops
		resetting the variables, resulting in the next call to performIterations()
		simply returning true, indicating that refinement is complete.
		*/
		if (spokeObjFunc != spokeOptimizer->getLastBestVal()
			|| spokeGeomPenalty != spokeOptimizer->getLastPenalty())
			spokeOptimizer->readyNextRun();
		spokeObjFunc = spokeOptimizer->getLastBestVal();
		spokeGeomPenalty = spokeOptimizer->getLastPenalty();
	}

	if (spokeOptimizer->getCandidateObject() != NULL)
	{
		(*object) = (*(spokeOptimizer->getCandidateObject()));

		M3DFigure * figure = spokeOptimizer->getReferenceObject()->getFigurePtr(figureId);

		M3DFigureTreeNode * figureTreeNode = object->getFigureTreeNode(figureId);

		int i = 0;
		M3DFigureTreeNode * subFigNode = figureTreeNode->getChild(i);
		while (subFigNode != NULL)
		{
			subfigureTransform.init(object, subFigNode);
			// FIXME: What should we do in case of tubes? - gst/rrs
			subfigureTransform.updateSubfigure(dynamic_cast<M3DQuadFigure*>(figure));
			i++;
			subFigNode = figureTreeNode->getChild(i);
		}
	}
	object->setModified(true);
	modelNeedsSaving(true);
	endOperation();
	return ret;
}

double P3DControl::getLastSpokeObjectiveFunctionValue(bool preview)
{
	if (spokeOptimizer == NULL)
		return 0.0;

	return spokeOptimizer->getLastBestVal(preview);
}

double P3DControl::getLastSpokeGeometricMatchValue()
{
	if (spokeOptimizer == NULL)
		return 0.0;

	return spokeOptimizer->getLastPenalty();
}

#ifndef BINARY


// This function cannot be used to change a previously set weight
void P3DControl::setSpokePenaltyWeight(gpTune_t tuningParam, double w)
{
	if (penaltyWeights == NULL) {
		penaltyWeights = new double[MAX_NUM_ATOM_MATCH_RESULTS];
		penaltyNames = new gpTune_t[MAX_NUM_ATOM_MATCH_RESULTS];
	}
	if (numPenaltyWeights >= MAX_NUM_ATOM_MATCH_RESULTS) {
		cout << "Error: too many penalty weights\n";
		return;
	}
	penaltyWeights[numPenaltyWeights] = w;
	penaltyNames[numPenaltyWeights++] = tuningParam;
}

#endif	/* BINARY */



/**
* Voxel optimizer.
*/
bool P3DControl::initializeVoxelOptimizer(int surfaceLevel)
{
#ifdef DEBUG
	cout << "P3DControl::initializeVoxelOptimizer()\n";
#endif
	if (match == NULL || object == NULL)
		return false;

	if (voxelStageOptimizer == NULL)
		voxelStageOptimizer = new M3DVoxelOptimizer;
	selectionChanged = false;

	return voxelStageOptimizer->
#ifdef BINARY
		initialize(match, object, optimizerFigureId, surfaceLevel);
#else
		initialize(match, object, optimizerFigureId, surfaceLevel, optimizeWithWindowing);
#endif

#ifdef OPTIMIZATION_VISUALIZER
	globalLogManager.setOptimizer(NULL);
#endif
}

// Boundary deformation stage of the optimizer
void P3DControl::doVoxelIterations(int numIterations, bool preview)
{
	Image3D * image;

	if (object == NULL) {
		cout << "NULL object." << endl;
		return;
	}

	image = imagePlanes.getImagePtr();
	if (image == NULL) {
		cout << "NULL image." << endl;
		return;
	}

	if (voxelStageOptimizer == NULL) {
		cout << "NULL optimizer." << endl;
		return;
	}

	if (! preview)
		startOperation();

	voxelStageOptimizer->performIterations(numIterations);
	if (preview)
		return;

	if (voxelStageOptimizer->getTargetObject() != NULL) {
		(*object) = (*(voxelStageOptimizer->getTargetObject()));

		object->setModified(true);
		modelNeedsSaving(true);
	}

	endOperation();
}

void P3DControl::setVoxelPenaltyWeight(double w)
{
	if (voxelStageOptimizer == NULL)
		voxelStageOptimizer = new M3DVoxelOptimizer;

	voxelStageOptimizer->setPenaltyWeight(w);
}

void P3DControl::setVoxelConstraintsPenaltyWeight(double w)
{
	if (voxelStageOptimizer == NULL)
		voxelStageOptimizer = new M3DVoxelOptimizer;

	voxelStageOptimizer->setConstraintsPenaltyWeight(w);
}

double P3DControl::getLastVoxelMatchValue()
{
	if (voxelStageOptimizer == NULL)
		return 0.0;

	return voxelStageOptimizer->getLastMatchValue();
}

double P3DControl::getLastVoxelPenaltyValue()
{
	if (voxelStageOptimizer == NULL)
		return 0.0;

	return voxelStageOptimizer->getLastPenalty();
}

double P3DControl::getLastVoxelObjectiveFunctionValue()
{
	if (voxelStageOptimizer == NULL)
		return 0.0;

	return voxelStageOptimizer->getLastObjectiveFunctionValue();
}

bool P3DControl::voxelStageOptimizerIsFinished()
{
	if (voxelStageOptimizer == NULL)
		return false;

	return voxelStageOptimizer->isFinished();
}

// ---------------------------- End of Optimization Code -------------------------------

void P3DControl::setIntensityWindow(double min, double max)
{
	Image3D * image;

	image = imagePlanes.getImagePtr();
	if(image == NULL)
		return;

	imagePlanes.setIntensityWindow(min, max);
	image->intensityWindow(min, max);
}

void P3DControl::getIntensityWindow(double & min, double & max)
{
	Image3D * image;

	image = imagePlanes.getImagePtr();
	if(image == NULL) {
		min = 1.0;
		max = -1.0;
		return;
	}

	min = image->getMinIntensityWindow();
	max = image->getMaxIntensityWindow();
}

/*  Attach the tree (or partial tree) starting at childId to figure parentId.
Note that the tree containing the subtree to be attached may be the same
as the tree to which it will be attached, resulting in moving of a subtree.
The attachment mode should be either M3DFigureTreeNode::PROTRUDE or
M3DFigureTreeNode::INDENT.
*/
bool P3DControl::attachSubfigure(int parentId, int childId,
								 M3DFigureTreeNode::SubfigureAttachment_t mode)
{
	M3DFigureTreeNode * treeNode;
	M3DFigure * parentFigure;
	M3DQuadFigure * childFigure;
	M3DPrimitive * primitive;
	M3DPrimitiveLinkInfo * link;
	Vector3D figCoord;
	double u, v, t;
	int i, row, col;

	static const int ATTACH_MEDIAL_SUBDIVISIONS = 6;


	if(editingMode != ANY_EDITING_MODE)
		return false;

	if(object == NULL)
		return false;

	// FIXME: What should we do in case of tubes? - gst/rrs
	parentFigure = object->getFigurePtr(parentId);
	childFigure = dynamic_cast<M3DQuadFigure*>(object->getFigurePtr(childId));
	if(parentFigure == NULL || childFigure == NULL)
		return false;

	int op_count = startOperation();

	treeNode = object->getFigureTreeNode(childId);
	treeNode->setAttachmentMode(mode);

	for (i = treeNode->getLinkCount() - 1; i >= 0; i--)
		treeNode->removeLink(i);
	// In M3DObject::detachFigureTreeNode() the blend settings are retained
	if (treeNode->getBlendAmount() == 0.0)
		treeNode->setBlendAmount(0.5);
	if (treeNode->getBlendExtent() == 0.0)
		treeNode->setBlendExtent(0.5);

	// Require that the attached subfigure be oriented so that U is
	// pointed away from the main figure.
	// This code assumes either a full row or a full column is selected.
	row = childFigure->rowSelected();
	if (row < 0) {
		// A column is selected
		col = childFigure->columnSelected();
		if (col < 0) {
			// wrong selection of hinge (link) atoms!
			return false;
		}

		if (col != 0) {
			// rotate the labeling for 180 degree
			childFigure->reverseRows();
			childFigure->reverseColumns();
		}
	}
	else {	// A row is selected
		// need to re-label the atoms
		int j;
		M3DQuadFigure * newChildFigure =
			new M3DQuadFigure(childFigure->getColumnCount(), childFigure->getRowCount());
		M3DQuadFigure * oldChildFigure = (M3DQuadFigure *)(childFigure->clone());

		if (row == 0) {
			//	(1 0) (1 1) (1 2) -> (2 1) (1 1) (0 1)
			//	(0 0) (0 1) (0 2)    (2 0) (1 0) (0 0)
			// this wasn't working because an m-rep has an internal "orientation"
			//childFigure->reverseColumns();
			//childFigure->reverseRows();
			for (i=0; i<newChildFigure->getRowCount(); i++)
				for (j=0; j<newChildFigure->getColumnCount(); j++)
					newChildFigure->setPrimitivePtr(i, j, 
					oldChildFigure->getPrimitivePtr(j, childFigure->getColumnCount()-1-i));
		}
		else {
			//	(0 2) (0 1) (0 0)    (2 1) (1 1) (0 1)
			//	(1 2) (1 1) (1 0) -> (2 0) (1 0) (0 0)
			for (i=0; i<newChildFigure->getRowCount(); i++)
				for (j=0; j<newChildFigure->getColumnCount(); j++)
					newChildFigure->setPrimitivePtr(i, j, 
					oldChildFigure->getPrimitivePtr(childFigure->getRowCount()-1-j, i));
		}

		*childFigure = *newChildFigure;

		delete oldChildFigure;
	}

	// Because the child has only one parent, and the child's links were discarded
	// earlier, the tree does not have to be updated to reflect the above transformation.

	Xferlist * xferList = convertM3DtoXfer(parentFigure);

	ThallCode::Pointlist_server2 pList;
	pList.init(xferList);
	pList.ComputeSubdivBoundaryTiles(ATTACH_MEDIAL_SUBDIVISIONS);

	// Create the attachment, moving the link primitives to the surface of the parent
	// and all other subfigure atoms accordingly
	int	linkId = 0, 
		linkNum = 0,
		rowNum = childFigure->getRowCount(),		// v = row
		colNum = childFigure->getColumnCount();		// u = col
	Vector3D * shifts;
	int * linkRows,
		* linkCols;
	linkNum = childFigure->numberSelected();
	shifts = new Vector3D[linkNum];
	linkCols = new int[linkNum];
	linkRows = new int[linkNum];

	// to translate all the link atoms and record the translation
	for (i = 0; i < childFigure->getPrimitiveCount(); i++)
	{
		Bpoint Bpnt;
		double dist;
		bool outside;
		Vector3D pos;

		primitive = childFigure->getPrimitivePtr(i);

		if (primitive == NULL)
			continue;

		if (primitive->isSelected()) {
			pos = primitive->getX();
			if (! pList.subdivNearSurface(pos.getX(), pos.getY(), pos.getZ(), &Bpnt, &dist, &outside)) {
				cout << "Error in computing figural coordinates: cannot perform attachment" << endl;
				// Undo the primitive shifts and model reindexing
				op_count = endOperation() - op_count;
				resetOperations(op_count);
				return false;
			}
			u = Bpnt.u;
			v = Bpnt.v;
			t = Bpnt.t;
			pos = Vector3D(Bpnt.pnt[0], Bpnt.pnt[1], Bpnt.pnt[2]);

			Vector3D shift = pos - primitive->getX();
			if (globalVerbosity >= 1) {
				cout << "Translating Primitive " << i << " at ";
				cout << "u: " << u << ", v: " << v << ", t: " << t << " by \n\t";
				shift.print();
			}
			primitive->translateBy(shift);
			primitive->toggleHinge(true);
			link = new M3DPrimitiveLinkInfo;
			link->primitiveId = i;
			link->u = u;
			link->v = v;
			link->t = t;
			treeNode->addLink(link);

			shifts[linkId] = shift;
			childFigure->figuralCoordinates(i, linkRows[linkId], linkCols[linkId]);
			linkId++;
		}
	}

	if (linkId != linkNum) {
		// should not happened unless there is mem corruption
		op_count = endOperation() - op_count;
		resetOperations(op_count);
		cout << "Translations discarded; figure not attached" << endl;
		beep();
		return false;
	}

	// to translate all the other subfigure atoms using the recorded translation
	for (i = 0; i < childFigure->getPrimitiveCount(); i++) {
		primitive = childFigure->getPrimitivePtr(i);

		if (primitive == NULL)
			continue;

		if (primitive->isSelected()) {
		}
		else {
			int rowDex, colDex;
			childFigure->figuralCoordinates(i, rowDex, colDex);

			// find which column this atom belongs to
			int j;
			for (j=0; j<linkNum; j++) {
				if (rowDex == linkRows[j])
					break;
			}

			if (j >= linkNum) {
				// should not happened unless there is mem corruption
				op_count = endOperation() - op_count;
				resetOperations(op_count);
				cout << "Translations discarded; figure not attached" << endl;
				beep();
				return false;
			}
			else {
				double ratio;
				if (linkCols[j] == 0) {
					if (globalVerbosity >= 1) {
						cout << "Translating Primitive " << i << " by \n\t";
						shifts[j].print();
					}

					ratio = 1.0 - 1.0*colDex/(colNum-1);
					primitive->translateBy(shifts[j] * ratio);
				}
				else {
					if (globalVerbosity >= 1) {
						cout << "Translating Primitive " << i << " by \n\t";
						shifts[j].print();
					}

					ratio = 1.0*colDex/(colNum-1);
					primitive->translateBy(shifts[j] * ratio);
				}
			}
		}
	}
	delete []shifts;
	delete []linkCols;
	delete []linkRows;

	if (object->attachFigureTreeNode(treeNode, parentId)) {
		endOperation();
		if (globalVerbosity >= 0)
			cout << "Attached figure " << object->getFigureIndex(childFigure) << " to figure "
			<< object->getFigureIndex(parentFigure) << endl;

		markedPrimitiveId = INVALID_PRIMITIVE_ID;
		modelNeedsSaving(true);
		selectionChanged = true;
	}
	else {
		// Undo the primitive shifts and model reindexing
		op_count = endOperation() - op_count;
		resetOperations(op_count);
		cout << "Translations discarded; figure not attached" << endl;
		beep();
		return false;
	}
	return true;
}

bool P3DControl::detachSubfigure(int figureId)
{
	M3DFigureTreeNode * treeNode;
	M3DFigure * figure;
	bool ret;

	if(editingMode != ANY_EDITING_MODE)
		return false;

	if(object == NULL)
		return false;

	figure = object->getFigurePtr(figureId);
	if(figure == NULL)
		return false;

	setEditingMode(NO_EDITING_MODE);

	(void) startOperation();
	treeNode = object->detachFigureTreeNode(figureId);
	endOperation();

	if (treeNode != NULL) {
		if (globalVerbosity >= 0)
			cout << "Detached subfigure " << figureId << endl;
		ret = true;
	}
	else {	// This should never happen
		cout << "Subfigure could not be detached" << endl;
		ret = false;
	}

	modelNeedsSaving(true);
	selectionChanged = true;
	setEditingMode(ANY_EDITING_MODE);
	return ret;
}

void P3DControl::initializeSubfigureTransformation()
{
	int figureId,
		primId;

	M3DFigureTreeNode * figureTreeNode;

	if(object == NULL)
		return;

	getRelativeMarkedPrimitiveId(figureId, primId);
	if (figureId < 0)
		return;
	figureTreeNode = object->getFigureTreeNode(figureId);
	if(figureTreeNode == NULL)
		return;

	subfigureTransform.init(object, figureTreeNode);
}

void P3DControl::translateSubfigure(double u, double v)
{
	if(editingMode != ANY_EDITING_MODE)
		return;

	subfigureTransform.translate(u, v);

	modelNeedsSaving(true);
	selectionChanged = true;
}

void P3DControl::rotateSubfigure(double angle)
{
	if(editingMode != ANY_EDITING_MODE)
		return;

	subfigureTransform.rotate(angle);

	modelNeedsSaving(true);
	selectionChanged = true;
}

void P3DControl::scaleSubfigure(double scale)
{
	if(editingMode != ANY_EDITING_MODE)
		return;

	subfigureTransform.widen(scale);

	modelNeedsSaving(true);
	selectionChanged = true;
}

void P3DControl::elongateSubfigure(double elongation)
{
	if(editingMode != ANY_EDITING_MODE)
		return;

	subfigureTransform.elongate(elongation);

	modelNeedsSaving(true);
	selectionChanged = true;
}

void P3DControl::hingeSubfigure(double angle)
{
	if(editingMode != ANY_EDITING_MODE)
		return;

	subfigureTransform.hinge(angle);

	modelNeedsSaving(true);
	selectionChanged = true;
}

bool P3DControl::haveModel()
{
	if (object == NULL)
		return false;
	else
		return true;
}

bool P3DControl::haveImage()
{
	Image3D * image;

	image = imagePlanes.getImagePtr();
	if (image == NULL)
		return false;
	else
		return true;
}

void P3DControl::gaussBlur()
{
	Image3D * image;
	ImageResample3D resampler;
	GaussianBlur3D gBlur;

	image = imagePlanes.getImagePtr();

	if (image == NULL)
		return;

	gBlur.blur(*image, 1.0, 7);

	blurred = true;
	imagePlanesChanged = true;
}

// Load image file's header and return max extent of all dimensions.
double P3DControl::calcMaxExtent(const char * imageFile)
{
	if (! loadImage(imageFile, false, false, true))
		return 0.0;

	Image3D * image = imagePlanes.getImagePtr();
	if (! image)
		return 0.0;

	return image->maxExtent();
}

#ifdef BINARY

// All dimensions are up-sampled to maxExtent cm
void P3DControl::resampleImageUnitCubic(float maxExtent)
{
	Image3D *image = imagePlanes.getImagePtr();
	if (! image)
		return;

	// Make xSpacing = ySpacing = zSpacing
	ImageResample3D resampler;
	resampler.isotropicSample(*image, tuningWt(BpVoxelSpacing));

	// Threshold? nope, not if input is binary

	// Pad to fit in a cube in WORLD coords
	resampler.cubeSample(*image, maxExtent);
	imagePlanesChanged = true;
}

// All dimensions are up-sampled to maxExtent cm
void P3DControl::resampleImage()
{
	Image3D *image = imagePlanes.getImagePtr();
	if (! image)
		return;

	// Make xSpacing = ySpacing = zSpacing
	ImageResample3D resampler;
	resampler.isotropicSample(*image, tuningWt(BpVoxelSpacing));

	imagePlanesChanged = true;
}

// Computes a distance map for later use
bool P3DControl::distMap(const char *filename)
{
	Image3D * image = imagePlanes.getImagePtr();
	if(image == NULL)
		return false;

	ImageDistanceMap * binaryDistanceMap;

#if NEW_DM_CODE		/* AGG: This code is under development */
	// Create a temporary image to avoid destroying the original image
	Image3D * tempImage = new Image3D;
	Vector3D origin = image->getWorldOrigin();
	GreyValue * voxels = new GreyValue[image->getVoxelCount()];
	memcpy(voxels, image->getVoxels(), image->getVoxelCount()*sizeof(GreyValue));
	tempImage->setVoxels(voxels, image->getXDim(), image->getYDim(), image->getZDim());
	tempImage->setSpacingAndOrigin(image->getXSpacing(), image->getYSpacing(),
		image->getZSpacing(), &origin);
	tempImage->pushImageIsStacked(image->getIsImageStacked(), image->getStackedMask());

	// Create a distance map from the image
	binaryDistanceMap = new ImageDistanceMap(tempImage);
#else

	// Create a distance map from the image
	binaryDistanceMap = new ImageDistanceMap(image);
#endif
	(void) binaryDistanceMap->createMap();
	const int * map_size = binaryDistanceMap->getMapSize();
	const double * map_spacing = binaryDistanceMap->getMapSpacing();
	if (map_size[0] <= 0) {
		cout << "Error: Image distance map_size of 0!" << endl;
		return false;
	}
	if (binaryDistanceMap->status() == false) {
		cout << "Error: No image distance map!" << endl;
		return false;
	}
#if NEW_DM_CODE
	delete tempImage;
#endif

	//binaryDistanceMap->setSpacingAndOrigin(image);	// copy flipped coords
	setDistMapPtr(binaryDistanceMap);

	Image3D * im = binaryDistanceMap->operator Image3D *();

	// Create image file of the distance map extrapolated over the
	// full input image

	// Calculate intensity min/max range so that saved image can be signed
	if (filename && *filename) {
		GreyValue maxg, ming;

		// AGG: This is wrong, since dist images are shorts
		calc_intensity_range(*im, maxg, ming);

		if (globalVerbosity > -1)
			cout << "Saving distance map " << filename
			<< " with range " << ming << " .. " << maxg <<  endl;
		saveImage(filename, im, ming, maxg);
	}
	// And load the image into the display
	//control->LoadImage(im);
	delete im;	// This was created only to write the file

	//cout << "Distance map created!" << endl;
	return true;
}

#endif

// Outputs a histogram of the current image (for plotting in gnuplot, etc.)
void P3DControl::outputHistogram(const char * filename)
{
	Image3D * image;
	FILE * fp;
	GreyValue * dataPtr;
	int size;
	int * histogram;
	int i;

	image = imagePlanes.getImagePtr();

	if(image == NULL)
		return;

	fp = fopen(filename, "wb");
	if(fp == NULL)
		return;

	dataPtr = image->getVoxels();
	size = image->getXDim() * image->getYDim() * image->getZDim();

	histogram = new int[MAX_GREY_VALUE - MIN_GREY_VALUE + 1];
	if(histogram == NULL)
		return;

	memset(histogram, 0, sizeof(int) * (MAX_GREY_VALUE - MIN_GREY_VALUE + 1));
	for(i = 0; i < size; i++)
	{
		histogram[*dataPtr]++;
		dataPtr++;
	}

	for(i = MIN_GREY_VALUE; i <= MAX_GREY_VALUE; i++)
		fprintf(fp, "%f\t%d\n", (double)i/(double)(MAX_GREY_VALUE - MIN_GREY_VALUE), histogram[i]);

	fclose(fp);
	if (globalVerbosity >= 0)
		cout << "Wrote historgram to file " << filename << endl;
	delete [] histogram;
}

void P3DControl::recordTransformation()
{
	if(object == NULL || recordingTransformation)
		return;

	recordingTransformation = true;

	if(referenceObject == NULL)
		referenceObject = new M3DObject();

	object->select();
	(*referenceObject) = (*object);
	allowSelection(false);
	setEditingMode(SIMILARITY_ONLY_EDITING_MODE);

	if (recordedTransformation == NULL)
		recordedTransformation = new SimilarityTransform3D;
	else
#ifdef AE2_BUILD
		// do not reset the recorded sim transform for ae2
		;
#else
		recordedTransformation->setToIdentity();
#endif
}

void P3DControl::stopRecordingTransformation()
{
	if(object == NULL || referenceObject == NULL ||
		!recordingTransformation)
		return;

	recordingTransformation = false;
	allowSelection(true);
	setEditingMode(ANY_EDITING_MODE);

	delete objectSimTransform;
	objectSimTransform = new SimilarityTransform3D(*recordedTransformation);
}

void P3DControl::applyRecordedTransformation()
{
	if(object == NULL)
		return;

	if (recordedTransformation == NULL)
		recordedTransformation = new SimilarityTransform3D;
	if (objectSimTransform)
		*recordedTransformation = *objectSimTransform;

#ifdef UNFLIPPED
	object->applySimilarity(*recordedTransformation, Vector3D(0.0, 1.0, 0.0));
#else
	object->applySimilarity(*recordedTransformation, Vector3D(0.0, 0.0, 0.0));
#endif
	modelNeedsSaving(true);
	selectionChanged = true;
}

void P3DControl::resetTransformation()
{
	if(referenceObject == NULL)
		return;

	(*object) = (*referenceObject);

	modelNeedsSaving(true);
}

// -------------------------  Functions for viewing CPNS Deformations  -------------------------

bool P3DControl::haveCPNS() {
	if( cpnsStats != NULL )
		return( true ) ;
	else
		return( false ) ;

}

bool P3DControl::cpnsDeform(std::vector<double> & vals)
{

	int nEigenmodes = cpnsStats->getNEigenModes() ;

	double * scores = new double[ cpnsStats->getNEigenModes() ] ;

	for( int i = 0 ; i < cpnsStats->getNEigenModes() ; i++ )
		scores[i] = 0.0 ;

	for( int i = 0 ; i < vals.size() ; i++ )
		scores[i] = vals[i] ;

	M3DObject * newObject = cpnsStats->eigenmodeDeformMean( scores, object ) ;

	if ( newObject == NULL )
		return false;	

	if (object)
		delete object;
	// object = NULL ;

	object = newObject->assign();

	modelNeedsSaving(true);

	delete newObject ;
	delete [] scores ;

	return true ;
}



void P3DControl::resetCPNSStage() {
	
	if (object) {
		object->restore();
		modelIsSaved();
	}
}


// -------------------------  End of functions for viewing CPNS Deformations  -------------------------

// -------------------------  Functions for viewing PGA Deformations  -------------------------

#ifndef BINARY
void P3DControl::clearPGA()
{
	pgaStats->clear();
}
#endif

bool P3DControl::havePGA()
{
	return (haveFigurePGA() || havePrimitivePGA());

}

bool P3DControl::haveFigurePGA()
{
	if (pgaStats == NULL )
		return false;

	if (pgaStats != NULL){
		if (pgaStats->getNumOfPGSets() > 0 )
			return true;
	}

	return false;
}


bool P3DControl::havePrimitivePGA()
{
	if (atomPgaStats == NULL)
		return false;

	if (atomPgaStats->getPGDataPtr(0,0)->numPGs > 0)
		return true;

	return false;
}

void P3DControl::resetResiduePGAStage()
{
#ifdef DEBUG
	cout << "P3DControl::resetResiduePGAStage()\n";
#endif
	residuePGAStage = 0;
	if (pgaStats != NULL)
		pgaStats->resetObjsToMean();
	if (object) {
		object->restore();
		modelIsSaved();
	}
}

int P3DControl::getResiduePGAStage() 
{
	return residuePGAStage;
}

void P3DControl::setResiduePGAStage(int depth) 
{
#ifdef BINARY
	cout << "Error: P3DControl::setResiduePGAStage called" << endl;
#else
#ifdef DEBUG
	cout << "P3DControl::setResiduePGAStage(" << depth << ")\n";
#endif
	if (depth >= pgaStats->getNumOfPGSets() || depth < 0) {
		if (object) {
			object->restore();
			modelIsSaved();
		}
		resetResiduePGAStage();
	}
	else
		residuePGAStage = depth;
#endif
}

char * P3DControl::getResiduePGAStageName()
{
	PGSet * pgSet = pgaStats->getPGSetPtr(residuePGAStage);

	//show the depth of scale 
	if (pgSet->name == NULL)
	{
		pgSet->name = new char[50];
		strcpy(pgSet->name, "Level : ");
		char buf[10];
		sprintf(buf, "%d", residuePGAStage + 1);
		strcat(pgSet->name, buf);
	}
	return pgSet->name;
}

bool P3DControl::applyMeanResidue()
{
	if (object)
		delete object;
	if (pgaStats != NULL){
		object = pgaStats->getLatestObjPtr()->assign();

#ifdef DEBUG
		cout << "applyMeanResidue : resiudePGAStage " << residuePGAStage << endl;
#endif

		bool retVal = pgaStats->applyMeanResidue(object, residuePGAStage);
		if (! retVal) {
			cout << "check : pgaStats.applyMeanResidue()" << endl;
		}
		return retVal;
	}
	else if (atomPgaStats != NULL){
		object = atomPgaStats->getAlignReferenceObj()->assign();
		modelNeedsSaving(true);
		return true;
	}

	modelNeedsSaving(true);

	return false;
}

#ifndef BINARY
void P3DControl::initPGA()
{
	if (object == NULL)
		return;

	if (pgaStats != NULL) {

		//JJ: What does this do?
		pgaStats->setHandPlacementTransform(*handPlacementTransform);

#ifdef DEBUG
		cout << "\nHandPlacement in P3DControl::initPGA()\n";
		handPlacementTransform->getTranslation().print();
		handPlacementTransform->getRotation().print();
		cout << "scale " << handPlacementTransform->getScale() << endl;
#endif
	}
	return;
}
#endif

bool P3DControl::readPGAFile(const char * filename)
{
	M3DObjectFile objectFile;
	M3DObject * PGA_object;

	if (filename == NULL)
		return false;

	PGA_object = objectFile.read(filename);
	if (PGA_object == NULL)
		return false;

	M3DPGAStats * ps = PGA_object->getPGAStats();
	if (ps == NULL)
		return false;
	pgaStats = new M3DPGAStats(*ps);
	// pgaStats does not contain pointers to PGA_object
	delete PGA_object;

	object->setPGAStats(pgaStats);
	if (pgaStats != NULL) {
		object->restore();	// See read()
		pgaStats->setMeanObj(object->loadedObject());	// object->that
		pgaStats->resetObjsToMean();

		// See M3DObjectFile.cpp
		if (isAdaptive())
			pgaStats->convertMeanDiff();
		else
			pgaStats->rescale();

		modelIsSaved();	// See read()
	}

	return true;
}

/*
1. get the current depth
2. If  depth < set.Size
1. get the value from control
2. show the window/sliders, as the slider changes, pass the value to
control and deform the model
Else go back to the 0 depth and reset deformation.
*/
bool P3DControl::pgaDeform(std::vector<double> & vals)
{
	bool retVal;

#ifdef BINARY
	//Hack: only used in Binary Pablo, when the model contains atom stats.
	if (atomPgaStats != NULL) {
		retVal = pgaDeformAtom(vals);
		return retVal;
	}
	else{
		cout<<"no atom PGA!"<<endl;
	}
#endif

	// Must start with intermedObj
	if (object)
		delete object;

	//JJ temporary
	M3DObject * intermedObject = pgaStats->getIntermedObjPtr();
	if (intermedObject == NULL)
		return false;
	object = intermedObject->assign();

	retVal = pgaStats->doPGADeform(object, vals, residuePGAStage, true);
	if (! retVal)
		return retVal;

	modelNeedsSaving(true);

	return retVal;
}

// Atom PGA slider, selecting one atom each time
bool P3DControl::pgaDeformAtom(std::vector<double> & vals)
{
	int i, j, figId,atomId;
	M3DPrimitive * selectAtom;
	atomId = 0;
	int numPrimitives;

	M3DObject * intermedObject = atomPgaStats->getAlignReferenceObj();

	if (intermedObject == NULL)
		return false;

	//
	int numFigs = object ->getFigureCount();	
	for (j = 0;j < numFigs; j++ ){

		numPrimitives = object->getFigurePtr(j)->getPrimitiveCount();

		for (i = 0; i < numPrimitives; i++) {

			selectAtom = object->getFigurePtr(j)->getPrimitivePtr(i);

			if (selectAtom->isSelected()) {
				atomId = i;
				figId = j;				
				i = numPrimitives;
				j =  numFigs;
			}
		}
	}


	object = intermedObject->assign();
	M3DDeformationOptimizer * deformOp = new M3DDeformationOptimizer;
	deformOp->initializeAtomPGAStage( atomPgaStats, figId, object);

	Vector tempVal;
	tempVal.setSize(vals.size());

	for (i = 0; i < vals.size(); i++){
		tempVal(i) = vals[i];//*((double) tuningWt(AtomPGAFactor));
		//	cout<<"x["<<i<<"] ="<<tempVal(i)<<"\n"<<endl;
	}


	M3DFigure * currFigPtr = object->getFigurePtr(figId);
	M3DPrimitive * readyToDeformAtom = (atomPgaStats->getDeltaMeanPrimitivePtr(figId,atomId))->copyPtr(); 

	M3DPrimitive * newPrim = readyToDeformAtom->copyPtr();
	atomPgaStats->doPGADeform(newPrim, tempVal, figId, atomId, readyToDeformAtom);

	M3DFigure * alignReferenceFigurePtr = intermedObject->getFigurePtr(figId);
	deformOp->addDiffToPrimitive(currFigPtr,atomId, alignReferenceFigurePtr,readyToDeformAtom);

	//SimilarityTransform3D * xform = object->getTransformation();
	//	if (xform!= NULL )
	//	object->getFigurePtr(figId)->applySimilarity( *xform, 	Vector3D (0.0, 0.0, 0.0));


	delete newPrim;

	object->deselect();
	object->getFigurePtr(figId)->getPrimitivePtr(atomId)->select();


	modelNeedsSaving(true);

	return true;
}



// -------------------------  End of functions for viewing PGA Deformations  -------------------------

// -------------------------  Start of functions for viewing PCA Deformations  -------------------------

#ifndef BINARY

#ifdef PCA_VISUALIZER

bool P3DControl::initPCA_Deform()
{
	if (auxTileSet == NULL)
		return false;
	if (pca == NULL)
		return false;
	if (pca->status() == false)
		return false;

	pca->setMeanPCA((BYU *) auxTileSet);
	return true;
}

bool P3DControl::pcaDeform(std::vector<float> & vals)
{
	if (auxTileSet == NULL) {
		cout << "Mean PCA file is not loaded\n";
		return false;
	}

	if (pca == NULL) {
		cout << "PCA residue file is not loaded\n";
		return false;
	}

	if (! pca->pcaDeform((BYU*) auxTileSet, vals)) {
		cout << "PCA deformation failed\n";
		return false;
	}

	return true;
}

bool P3DControl::loadPCA_File(const char * filename)
{
	if (pca != NULL)
		delete pca;
	pca = new PCA;

	if (pca->readPCAfile(filename))
		return true;
	else
		return false;
}

#endif	/* PCA_VISUALIZER */

#endif	/* BINARY */

// -------------------------  End of functions for viewing PCA Deformations  -------------------------


//////////////////// MODEL SLIDE SHOW

void P3DControl::startModelTileSlideShow()
{
	slideShowSavedObject = object;
	object = NULL;

	slideShowSavedTileSet = tileSet;	// save current tiles
	tileSet = NULL;

	slideShowSavedEditingMode = editingMode;
	// GST: allow model editing during the show; assume user knows
	//  that these models/tiles might be lost when the show ends.
	//setEditingMode(NO_EDITING_MODE);	
	inSlideShow = true;
}

void P3DControl::addModelSlideShow(const char * filename)
{
	M3DObjectFile objFile;

	modelSlideShow.push_back(objFile.read(filename));
}

void P3DControl::removeModelSlideShow(int modelNumber)
{
	if(modelNumber >= 0 && modelNumber < modelSlideShow.size())
		modelSlideShow.erase(modelSlideShow.begin() + modelNumber);
}

void P3DControl::moveUpModelSlideShow(int modelNumber)
{
	M3DObject * tempObj;

	if(modelNumber > 0 && modelNumber < modelSlideShow.size())
	{
		tempObj = modelSlideShow[modelNumber];
		modelSlideShow[modelNumber] = modelSlideShow[modelNumber - 1];
		modelSlideShow[modelNumber - 1] = tempObj;
	}
}

void P3DControl::moveDownModelSlideShow(int modelNumber)
{
	M3DObject * tempObj;

	if(modelNumber >= 0 && modelNumber < modelSlideShow.size() - 1)
	{
		tempObj = modelSlideShow[modelNumber];
		modelSlideShow[modelNumber] = modelSlideShow[modelNumber + 1];
		modelSlideShow[modelNumber + 1] = tempObj;
	}
}

void P3DControl::clearModelSlideShow()
{
	int i;

	for(i = 0; i < modelSlideShow.size(); i++)
		delete modelSlideShow[i];

	modelSlideShow.clear();
	object = NULL;
}

void P3DControl::changeModelSlideShow(int number)
{
	if(number >= 0 && number < modelSlideShow.size())
		object = modelSlideShow[number];
}

void P3DControl::quitModelSlideShow()
{
	object = slideShowSavedObject;
	slideShowSavedObject = NULL;
	quitTileSlideShow();

	// see comment in startModelTileSlideShow()
	//setEditingMode(ANY_EDITING_MODE);
}

//////////////////// TILE SLIDE SHOW



void P3DControl::addTileSlideShow(const char * filename)
{
	TileSet *t = new TileSet(filename);
	tileSlideShow.push_back(t);
}

void P3DControl::removeTileSlideShow(int tileNumber)
{
	if(tileNumber >= 0 && tileNumber < tileSlideShow.size())
		tileSlideShow.erase(tileSlideShow.begin() + tileNumber);
}

void P3DControl::moveUpTileSlideShow(int tileNumber)
{
	TileSet* tempObj;

	if(tileNumber > 0 && tileNumber < tileSlideShow.size())
	{
		tempObj = tileSlideShow[tileNumber];
		tileSlideShow[tileNumber] = tileSlideShow[tileNumber - 1];
		tileSlideShow[tileNumber - 1] = tempObj;
	}
}

void P3DControl::moveDownTileSlideShow(int tileNumber)
{
	TileSet* tempObj;

	if(tileNumber >= 0 && tileNumber < tileSlideShow.size() - 1)
	{
		tempObj = tileSlideShow[tileNumber];
		tileSlideShow[tileNumber] = tileSlideShow[tileNumber + 1];
		tileSlideShow[tileNumber + 1] = tempObj;
	}
}

void P3DControl::clearTileSlideShow()
{
	int i;

	for(i = 0; i < tileSlideShow.size(); i++)
		delete tileSlideShow[i];

	tileSlideShow.clear();
	tileSet = NULL;
}

void P3DControl::changeTileSlideShow(int number)
{
	if(number >= 0 && number < tileSlideShow.size())
		tileSet = tileSlideShow[number];
}

void P3DControl::quitTileSlideShow()
{
	tileSet = slideShowSavedTileSet;	// restore saved tileset
	slideShowSavedTileSet = NULL;
	inSlideShow = false;
}

////////////////////

#ifndef PRODUCTION_VERSION

#ifndef BINARY
double P3DControl::testConstraintsPenalty(int figureId)
{
	if(object == NULL)
	{
		cerr << "Object is NULL" << endl;
		return 0.0;
	}
	if(figureId < 0 || figureId >= object->getFigureCount())
	{
		cerr << "Invalid Figure" << endl;
		return 0.0;
	}
	if (match == NULL) {
		cerr << "Match is NULL" << endl;
		return 0.0;
	}
	return match->computeConstraintsPenalty(*object, figureId);
}

void P3DControl::testSeurat(int figureId, int level) {
	if(object == NULL)
	{
		cerr << "Object is NULL" << endl;
		return;
	}

	if(figureId < 0 || figureId >= object->getFigureCount())
	{
		cerr << "Invalid Figure" << endl;
		return;
	}

	M3DFigure * figure = object->getFigurePtr(figureId);
	if(match == NULL)
		match = new Match;
	match->testSeurat(figure, level);
}
#endif	/* BINARY */

#endif	/* PRODUCTION_VERSION */

int P3DControl::machineByteOrder() const
{
	return cpuByteOrder();
}

/*	A model is adaptive, when the geometry statistics are computed on the
difference between a day-n model and a day-1 model.  The idea is that
we want to see how a model deforms from a specific reference model.

A model is not adaptive, when the statistics are computed on the models
themselves, not the differences between them.
*/
bool P3DControl::isAdaptive()
{
	if (object == NULL)
		return false;

	if (object->isAdaptive())
		return true;
	else 
		return false;
}

#ifdef AE2_BUILD

void P3DControl::addAnastruct(Anastruct3D * a)
{
	if (anaList == NULL)
		anaList = new AnastructList;
	anaList->addAnastruct(a);
}

bool P3DControl::convertAnastructCoords()
{
	if (anaList == NULL)
		return false;

	Image3D * im = getImagePtr();
	if (im == NULL)
		return false;

	return anaList->convertAnastructCoords(im);
}

void P3DControl::clearAnastructList()
{
	delete anaList;
	anaList = NULL;
}

#endif

bool P3DControl::align(const char * pointsFile, const char * modelFile,
					   const char * imageFile, double closeness, bool pgaOnly,
					   bool scaleOnly, bool scaleTransOnly, const char * outputTransform)
{
	M3DObject * model;
	Image3D * image;
	Image3D * img;

	if (! pointsFile || pointsFile[0] == '\0') {
		cout << "Unable to align; no points file was provided" << endl;
		return false;
	}

	if (modelFile == NULL) {
		if (object == NULL) {
			cout << "Unable to align without a model" << endl;
			return false;
		}
		else model = object;
	}
	else {
		if (! read(modelFile)) {
			cout << "Error: Unable to load model: " << modelFile << '\n';
			return false;
		}
		else {
			model = object;
			if (globalVerbosity > -1)
				cout << "Loaded model " << modelFile << '\n';
		}
	}

	image = getImagePtr();
	if (imageFile == NULL) {
		if (image == NULL) {
			cout << "Unable to align without an image" << endl;
			return false;
		}
		else img = image;
	}
	else {
		if (! loadImage(imageFile)) {
			cout << "Error: Unable to load image: " << imageFile << '\n';
			return false;
		}
		else {
			image = getImagePtr();
			img = image;
			if (globalVerbosity > -1)
				cout << "Loaded image " << imageFile << '\n';
		}
	}

	int order = 0;
	if (scaleOnly || scaleTransOnly)
		order = -1;
#ifdef FAST_ALIGNMENT
	FitUnlabeledPoints align(2, 0, order, closeness, 1.0e-1, 1.0e-1);
#else
	FitUnlabeledPoints align(2, 0, order, closeness, 1.0e-4, 1.0e-5);
#endif
	/*  for (i = 0; i < m; i++) {
	Vector3D v(x[i], y[i], z[i]);
	image->image_3d->worldToModelCoordinates(v);
	cout << v.getX() << ' ' << v.getY() << ' ' << v.getZ() << '\n';
	}
	return false;*/

	FitUnlabeledPoints::PointFitType type = FitUnlabeledPoints::SimPGA;
#ifdef FAST_ALIGNMENT
	type = FitUnlabeledPoints::Trans;
#else
	if (pgaOnly)
		type = FitUnlabeledPoints::PGA;
	else if (scaleOnly)
		type = FitUnlabeledPoints::Scale;
	else if (scaleTransOnly)
		type = FitUnlabeledPoints::ScaleTrans;
#endif
	if (! align.fit(model, pointsFile, img, type)) {
		cout << "Alignment of model " << modelFileName << " failed" << endl;
		return false;
	}
	else
		cout << "Aligned model " << modelFileName << endl;

	model = align.result();
	newModel(model);

	if (outputTransform != NULL) {
		SimilarityTransform3D placementSim;
		align.result(placementSim);	// AGG: Return using objectSimTransform?
		placementSim.writeSimilarity(outputTransform);
		cout << "Saved transformation " << outputTransform << endl;
	}

	return true;
}

void P3DControl::optAll(const char *filename)
{
    if(object == NULL || filename == "")
    {
        return;
    }
    AllImageIO imageFile;
    Image3D* image3D = imageFile.read(filename);
    ImageDistanceMap* binaryDistanceMap = new ImageDistanceMap(image3D);
    binaryDistanceMap->fromImage3D(image3D);
    M3DNEWUOAOptimizer optimizer;

    optimizer.setImage(binaryDistanceMap);
    optimizer.setObject(object);

    optimizer.perform(object);
    object->setModified(true);

    delete binaryDistanceMap;
}

void P3DControl::optAngle(const char* filename)
{
    if(object == NULL || filename == "")
    {
        return;
    }

    AllImageIO imageFile;
    Image3D* image3D = imageFile.read(filename);
    ImageDistanceMap* binaryDistanceMap = new ImageDistanceMap(image3D);
    binaryDistanceMap->fromImage3D(image3D);
    //M3DSpokeLengthOptimizer optimizer;
    M3DSpokeAngleOptimizer optimizer;

    optimizer.setImage(binaryDistanceMap);
    optimizer.setObject(object);

    optimizer.perform(object);
    object->setModified(true);

    delete binaryDistanceMap;
}
void P3DControl::optLength(const char* filename)
{
    if(object == NULL || filename == "")
    {
        return;
    }

    AllImageIO imageFile;
    Image3D* image3D = imageFile.read(filename);
    ImageDistanceMap* binaryDistanceMap = new ImageDistanceMap(image3D);
    binaryDistanceMap->fromImage3D(image3D);
    M3DSpokeLengthOptimizer optimizer;
    //M3DSpokeAngleOptimizer optimizer;

    optimizer.setImage(binaryDistanceMap);
    optimizer.setObject(object);

    optimizer.perform(object);
    object->setModified(true);

    delete binaryDistanceMap;
}
