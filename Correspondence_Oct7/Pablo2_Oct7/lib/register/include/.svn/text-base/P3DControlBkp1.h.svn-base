#ifndef P3D_CONTROL_H
#define P3D_CONTROL_H

// MSVC6: disable compiler warning about debugging symbol names too long
#pragma warning (disable:4786)

// Applies to BINARY:
//  Landmarks are used to initialize the model and for generating a penalty
//  during optimization by comparing the IMAGE landmarks to the MODEL landmarks
//  in order to minimize some distance function.
//  - Original method: image landmarks are entered at runtime and are not saved
//    into files.  Model landmarks may be saved in the model as a special construct.
//    Compiled when LM_METHOD_OBJ is not #define'd
//  - New method: image landmarks are provided as MODELS
//    to facilite moving, visualizing, scaling and creating them. If the user names
//    a figure specially (starting with "LM-"), then Edit->Add Quad Figure will
//    make all 2x2 atoms in the same exact position; other names results in
//    normal models. This special model is kept in control->landmarkObject &
//    is compiled when LM_METHOD_OBJ is #define'd.  Model landmarks are associated
//    with spoke ends and there can be up to 3 per atom.  These are stored in the
//    regular model files.
// Applies to non-BINARY:
//  LM_METHOD_OBJ MODEL landmarks (spoke tips) may be read, written and displayed
//  but do not yet contribute to any optimization penalty.

#include <vector>

// globalBuildParms.h is included here both for the local
// needs of this class, but also because P3DControl.h is
// the top-level API for the Pablo libraries.
#include "globalBuildParms.h"

#include "Vector2D.h"
#include "M3DObject.h"
#include "ImagePlanes.h"
#include "TileSet.h"
#include "P3DUndoList.h"
#include "Match.h"
#include "M3DRegistrationOptimizer.h"
#ifndef BINARY
#include "M3DRegistrationPGAOptimizer.h"
#include "M3DRegistrationSimilarityOptimizer.h"
#endif
#include "M3DMainFigureOptimizer.h"
#include "M3DSimilarityElongationOptimizer.h"
#include "M3DSimilarityPGAOptimizer.h"
#include "M3DSubfigureOptimizer.h"
#include "M3DDeformationOptimizer.h"
#include "M3DSpokeOptimizer.h"
#include "M3DSubfigureTransformation.h"
#include "M3DFigureElongater.h"
#include "M3DVoxelOptimizer.h"
#include "Registry.h"
#ifdef BINARY
#include "M3DPGA.h"
#endif
#ifdef AE2_BUILD
#include "Anastruct.h"
#endif
#include "Tuning.h"


enum OptimizerStage
{
    OPTIMIZER_SETUP_STAGE,
    OPTIMIZER_ENSEMBLE_STAGE,
    OPTIMIZER_OBJECT_STAGE,
    OPTIMIZER_FIGURE_STAGE,
    OPTIMIZER_MAIN_FIGURE_STAGE,
    OPTIMIZER_SUBFIGURE_STAGE,
    OPTIMIZER_PRIMITIVE_STAGE,
	OPTIMIZER_SPOKE_STAGE,
	OPTIMIZER_BOUNDARY_DISPLACEMENT_STAGE,
    OPTIMIZER_END_STAGE,
	NO_OPTIMIZER_STAGE
};

enum EditingMode
{
    ANY_EDITING_MODE,
    SIMILARITY_ONLY_EDITING_MODE,
    NO_EDITING_MODE
};

// List of linkable applications other than Pablo, for which
// special #define's are needed at build time.  AE2 was the
// pre-release name of the program later called ConStruct.
enum application_t { Pablo = 0, AE2, AE2Binary  };


struct imageLandmark_t {
	double coords[3];
	short figureId;
	short number;
    const char * name;
};

// Clip planes limit the view to show (x y z) that satify:
//   (A B C D) * inverse(MODEL) * transpose(x y z w) >= 0
// where MODEL is the current openGL modelview matrix.
// All openGL implementations support at least 6 clip planes.
struct clipPlane {
	int mask;		// enable/disable this plane by matching mask bits.
		// by convention, bit 0x8000 selects all planes
	GLdouble equation[4];	// A, B, C, D
};

#define CLIP_MASK_ALL     0x7E	// all possible masks
#define CLIP_MASK_ENABLED 0x01	// enables all other clips
#define CLIP_MASK_TILE    0x02	// clip tiles
#define CLIP_MASK_MODEL   0x04	// clip model surface
#define CLIP_MASK_AXES    0x08	// clip unit cube axes
#define CLIP_MASK_IMAGE   0x10	// clip image planes
#define CLIP_MASK_ILM     0x20	// clip image landmarks
#define CLIP_MASK_CAMERA  0x40	// clip ALL objects in scene
#define CLIP_MASK_VIEW    0x80	// same as camera but tracks view

class SimilarityTransform3D;
class M3DPGAStats;
class M3DPGAPrimitiveStats;
class M3DPGASpokeStats;
class PCA;
#ifndef BINARY
class SurfacePatchEnsemble;
class MultiPatch;
#endif
class DQFImage;


#if defined(SWIG) && defined(BINARY)
// Don't create accessors for optimRegistry as it doesn't have a copy constructor.
%ignore P3DControl::optimRegistry;
#endif

class P3DControl
{
public:

    P3DControl(int undoLength);
    ~P3DControl();
    bool verifyApplication(application_t id) const;
    bool verifyOrientation() const;   // True if UNFLIPPED is defined

    M3DObject * getObjectPtr() { return object; }
#ifdef OPTIMIZATION_VISUALIZER
	M3DObject * replaceObjectPtr(M3DObject* obj);
#endif

    ImagePlanes * getImagePlanesPtr() { return &imagePlanes; }
    Image3D * getImagePtr() { return imagePlanes.getImagePtr(); }
    bool imagePlanesStatus() { return imagePlanesChanged; }
    void toggleImagePlanesStatus() { imagePlanesChanged = ! imagePlanesChanged; }

	bool evaluateTileSetTypes(int numFiles, const char ** filenames,
		int & min, int & max);
    bool loadTileSet(const char * filename, bool split = false, bool useWorldCoords = false);
	bool appendTileSet(const char * filenames, bool split = false, bool useWorldCoords = false);
	bool saveTileSet(const char * filename, bool useWorldCoords = false);
    TileSet * getTileSetPtr() { return tileSet; }
	void clearTileSet();
	bool generateTileSet();
	bool haveTileSet() const;
	void setTileSetMobile(bool val) {tileSetMobile = val;}
	bool getTileSetMobile() { return tileSetMobile; }

	void toggleTileSetMobileFigure(int figId) {
		tileSetMobileFigure[figId] = !tileSetMobileFigure[figId];
		tileSet->calcCOG(tileSetMobileFigure);
	}

	bool loadAuxTileSet(const char * filename);
	void * getAuxTileSetPtr() { return auxTileSet; }

    void setEditingMode(EditingMode mode) { editingMode = mode; }
    void allowSelection(bool yesNo) { selection_allowed = yesNo; }
	bool selectionIsAllowed() const { return selection_allowed; }

    // Ownership of all loaded images passes to ImagePlanes
    bool loadImage(const char * filename, bool stacked = false, bool mapToWorld = false,
        bool headerOnly = false, bool external = false);
    // This is called to register a loaded image with this class and class ImagePlanes
    bool finishImageLoading(Image3D * image, bool stacked = false,
        bool mapToWorld = false, bool external = false);

#ifndef BINARY
	bool doCommandHistogramTraining(Registry & script);
	bool doCommandMultiPatchTraining(Registry & script);
	bool doCommandDistanceQuantileFunctions(Registry & script);
	bool doCommandDistanceQuantileFunctionImageMaker(Registry & script);
	bool doCommandDistanceQuantileFunctionModelMaker(Registry & script);
	static const char * calcDQFFilename(char * imageFilename, int figID);
#endif

	// Utility function to input an image without this class owning it
	Image3D * readImage(const char * filename, bool stacked = false,
		bool headerOnly = false) const;

	void setInputImageScaling(bool scaleInput);

	bool getIsImageStacked();

	// If external is true, this function doesn't free the Image3D object
	// (for use by programs other than Pablo).
	void discardImage(bool external = false);

	// Write the image provided to the specified file.  If no image is
	// provided, the one owned by this class is saved.  Normally, an
	// intensity range should not be specified, in which the range
	// [MIN_GREY_VALUE, MAX_GREY_VALUE] will be used.
	bool saveImage(const char * filename, Image3D * image = NULL,
		GreyValue lowIntens = 1, GreyValue highIntens = 0);

	bool imageIsBlurred() const { return blurred; }

	// Load and store DQF image files
	DQFImage * readDQF(const char * filename) const;
	bool saveDQF(const char * filename, DQFImage * dqf) const;

    int addQuadFigure(int numRows, int numColumns, char * figureName,
        float * color);

	// Regularize the model
	void regularize(double stepsize, int iterations);

#ifndef BINARY
	bool loadSPEObject(const char * speFilename, Registry &script);
	bool loadMPObject(const char *mpFilename, Registry &script);
#endif

    // Removes all selected figures
    void removeSelectedFigures();

    // Absolute translate and rotate    
    void translate(const Vector3D &trans);
    void rotate(Quat &q);

	void rotateTileSet(Quat &q);
	void pirouetteTileSet(Quat &q);

    void scale(double val);
    void scaleWidth(double val);

    // Elongation affects only selected figures
    void elongateRow(double val);
    void elongateColumn(double val);

    // Operations affecting only the marked primitive
    void setTheta(double val);
    void setElongation(double val);

    void toggleFigure(int figId);
    void togglePrimitive(int primId);

    void selectAll();
    void deselectAll();
    void toggleAll();

	void selectTree(M3DFigureTreeNode * node);
	void deselectTree(M3DFigureTreeNode * node);
	void toggleTree(M3DFigureTreeNode * node);

    // Copies the selected figures into the copy buffer
    void copy();

    // Pastes the contents of the copy buffer into the model
    void paste();

    // Mirrors the selected figures about the y-z plane
    void mirror();

    // Called when starting a transformation operation (adds to undo list)
    int startOperation();

    // Called when a transformation operation is completed
    int endOperation();

    // Undoes previous operations without possibility of redo
    void resetOperations(int npushes);

    // Undoes previous transformation operation
    bool undo();

    // Redoes the previous undone transformation operation
    bool redo();

    // Set the world coordinate positions of the image planes
    void setXCutPlanePos(double pos) { imagePlanes.setXCutPlaneWorldPos(pos); }
    void setYCutPlanePos(double pos) { imagePlanes.setYCutPlaneWorldPos(pos); }
    void setZCutPlanePos(double pos) { imagePlanes.setZCutPlaneWorldPos(pos); }

    // Set the slice positions of the image planes
    void setXCutPlaneSlice(int sliceIndex) { imagePlanes.setXCutPlaneSlice(sliceIndex); }
    void setYCutPlaneSlice(int sliceIndex) { imagePlanes.setYCutPlaneSlice(sliceIndex); }
    void setZCutPlaneSlice(int sliceIndex) { imagePlanes.setZCutPlaneSlice(sliceIndex); }

	// Get the world coordinate positions of the image planes
    double getXCutPlaneWorldPos() { return imagePlanes.getXCutPlaneWorldPos(); }
    double getYCutPlaneWorldPos() { return imagePlanes.getYCutPlaneWorldPos(); }
    double getZCutPlaneWorldPos() { return imagePlanes.getZCutPlaneWorldPos(); }

	// Get the model coordinate positions of the image planes
    double getXCutPlaneModelPos() { return imagePlanes.getXCutPlaneModelPos(); }
    double getYCutPlaneModelPos() { return imagePlanes.getYCutPlaneModelPos(); }
    double getZCutPlaneModelPos() { return imagePlanes.getZCutPlaneModelPos(); }

    // Turn smoothing of the image planes on or off
    void setPrefsSmoothImages(bool smooth) {
        imagePlanes.setSmoothing(smooth);
		imagePlanesChanged = true;
    }

    // Clear current model and start with a blank one
    void newModel();

    // Clear current model and start with a different one
    void newModel(M3DObject *newModel);  // Ownership of newModel passes to P3DControl

    // I/O routines to read/write a model from/to a .m3d file

    bool read(const char * filename, bool mapToWorld = false,
		const char * xformFilename = NULL,
		const char * pgaFilename = NULL);
    bool readAppend(const char * filename, bool mapToWorld, bool & ignoredPGA);
	bool loadSimTransform(const char * xformFilename);
#ifdef AE2_BUILD
	void applySimTransform(SimilarityTransform3D * transform = NULL);
#else
	void applySimTransform();
#endif
	bool write(const char * filename, const char * xformFilename = NULL);

    bool readOldFile(const char * filename);

	void applyWorld(Image3D * image);

    void exportToImage(const char * filename);
#ifndef PRODUCTION_VERSION
    void exportDistanceMap(const char * filename);
#endif
    void exportToBYUTiles(const char * filename, bool quads = false,
		bool useWorldCoords = false, int level = MATCH_POINT_CLOUD_SUBDIVISIONS,
		unsigned long figureMask = 0L, bool separate = false);

	// If the model is unsaved, then it needs to be saved to capture some change.
	// If it is modified, then it needs to be saved because a coordinate, angle,
	// or similar item changed.  Minor changes such as colors and figure names
	// do not cause the model to be "modified".  The primary reason for this
	// distinction is that PGA data is invalidated if the atoms of an MREP move.
    bool modelIsUnsaved() { return modelIsDirty; }
    bool modelIsModified() { return modelWasMoved; }
    void modelIsReplaced() { modelWasMoved = false; }//xx:in order to get rid of
	// PGA when saving the middle resulting mrep of optimization
    void modelNeedsSaving(bool modified = false) {
		modelIsDirty = true;
		modelWasMoved |= modified;
	}
	void modelIsSaved() {
		modelIsDirty = false;
		modelWasMoved = false;
	}
    char * getModelFileName();

	bool testReorderModel();
	bool reorderModel();

    int getMarkedPrimitiveId() { return markedPrimitiveId; }
    void getRelativeMarkedPrimitiveId(int & figIndex, int & primIndex);
    void toggleMarkedPrimitive(int id);
    M3DPrimitive * getMarkedPrimitive();
	void primitiveFiguralCoordinates(int figIndex, int primIndex, int & u, int & v);

    M3DPrimitive * getInterpolatedPrimitive() { return interpolatedPrimitive; }
    void setInterpolatedPrimitive(double u, double v);
    void clearInterpolatedPrimitive();

	bool haveImage();
	bool haveModel();

	// Landmarks
	void allocateImageLandmarks(int count);
	void insertImageLandmark(int at);
	void appendImageLandmark();
	bool setImageLandmark(int landmarkId, double x, double y, double z,
		int figureId, int number, const char * name = NULL);
	imageLandmark_t getImageLandmark(int landmarkId) const;
	int getImageLandmarksCount() const { return numImageLandmarks; }
	void clearImageLandmarks();
//	void clearImageLandmark(int figureId, int landmark);
	bool clearImageLandmark(int index);
	void deleteImageLandmark(int index);
	void deleteImageLandmarks();
	bool landmarkIsPaired(int figureId, int landmark);
//	bool anyLandmarkIsPaired(int figureId);
//	Vector3D landmarkCoordinates(int figureId, int landmark);
	int figureOfLandmark(int index, int & landmark);
	int landmarksAllocation() const { return imageLandmarks.size(); }
	vector<imageLandmark_t> & landmarkCoordinates() { return imageLandmarks; }
#ifndef BINARY
	void applyLandmarkTransformation();
#endif

    OptimizerStage getCurrentOptimizerStage() { return currentOptimizerStage; }
    OptimizerStage nextOptimizerStage();
	void cancelOptimization();
	int getOptimizationFigureNumber() { return optimizerFigureId; }
	void setOptimizerIntensityWindowing(bool toggle);
#ifndef BINARY
	bool getOptimizerIntensityWindowing() const { return optimizeWithWindowing; }
#endif
	OptimizerBase * getOptimzer(OptimizerBase::optimizer_t type);

    bool initializeMatch(MatchType type, int surfaceLevel, char * templateFilename = NULL, 
		const char * profileFilename = NULL, Image3D * trainingImage = NULL);
	MatchType type() const
    {
        if(match != NULL)
            return match->getType();
		else
			return NO_MATCH_TYPE;
    }
	// Filename of the three templates used in registration.
	char * getTemplatesFilename() const { return templates_filename; }
	void setTemplatesFilename(char * filename) { templates_filename = filename; }

#ifndef BINARY
    void setNormalPenaltyWeight(double w) {
        if(match != NULL)
            match->setNormalPenaltyWeight(w);
    }
    void setSlidingPenaltyWeight(double w) {
        if(match != NULL)
            match->setSlidingPenaltyWeight(w);
    }
    void setPenetrationPenaltyWeight(double w) { 
        if(match != NULL)
            match->setPenetrationPenaltyWeight(w);
    }
    void setCreasingPenaltyWeight(double w) {
        if(match != NULL)
            match->setCreasingPenaltyWeight(w);
    }

    double getNormalPenaltyWeight() {
        if(match != NULL)
            return match->getNormalPenaltyWeight();
        return -1.0;
    }
    double getSlidingPenaltyWeight() {
        if(match != NULL)
            return match->getSlidingPenaltyWeight();
        return -1.0;
    }
    double getPenetrationPenaltyWeight() { 
        if(match != NULL)
            return match->getPenetrationPenaltyWeight();
        return -1.0;
    }
    double getCreasingPenaltyWeight() {
        if(match != NULL)
            return match->getCreasingPenaltyWeight();
        return -1.0;
    }
#endif

	// This function is for use only within stages of optimization
	bool testAndResetModelChanged() {		// Self-clearing function
		bool flag = selectionChanged;
		selectionChanged = false;
		return flag;
	}

#ifndef BINARY
	// Model stage optimization functions
    void initializeRegistrationOptimizer();
	//void initializeRegistrationOptimizer(bool usePGA = false);  //JJ : NOT DONE!!!!
	void doRegistrationIterations(int numIterations);
    void setRegistrationPenaltyWeight(double w);
#endif
    bool registrationOptimizerIsFinished();
    double getLastRegistrationMatchValue();
    double getLastRegistrationPenaltyValue();

	// Main figure stage optimization functions
    void doMainFigureIterations(int numIterations, bool verbose = false);
    double getLastMainFigureGeometricMatchValue();
    double getLastMainFigureObjectiveFunctionValue();
#ifndef BINARY
    void clearPenaltyWeights();	// Discard penalty weights
    void setMainFigurePenaltyWeight(gpTune_t tuningParam, double w);
#endif
    bool mainFigureOptimizerIsFinished();
	bool mainFigureStageSkipped();

#ifdef BINARY
	bool initializeMainFigureOptimizer(int doScale, // bitFlags: 1:scale
		const char* pgaFilename);
#else
	bool initializeMainFigureOptimizer(bool rigid = false, bool usePGA = false);
#endif

#ifndef BINARY
	// Subfigure stage optimization functions
    void initializeSubfigureOptimizer();
//    void resetSubfigureOptimizer();
    void doSubfigureIterations(int numIterations);
    double getLastSubfigureMatchValue();
    bool subfigureOptimizerIsFinished();
#endif

	// Atom stage optimization functions
    bool initializeDeformationOptimizer(int appID = 0);
    M3DDeformationOptimizer * getDeformationOptimizer() const { return deformationOptimizer; }
    void resetDeformationOptimizer();
    bool doDeformationIterations(int figureId, bool preview = false, bool verbosity = false);
    double getLastDeformationObjectiveFunctionValue(bool preview = false);
    double getLastDeformationGeometricMatchValue();
#ifndef BINARY
	void setDeformationPenaltyWeight(gpTune_t tuningParam, double w);
#endif

	// Spoke stage optimization functions
    bool initializeSpokeOptimizer(int appID = 0);
    M3DSpokeOptimizer * getSpokeOptimizer() const { return spokeOptimizer; }
    void resetSpokeOptimizer();
    bool doSpokeIterations(int figureId, bool preview = false, bool verbosity = false);
    double getLastSpokeObjectiveFunctionValue(bool preview = false);
    double getLastSpokeGeometricMatchValue();
#ifndef BINARY
	void setSpokePenaltyWeight(gpTune_t tuningParam, double w);
#endif


	// Boundary stage optimization functions
	bool initializeVoxelOptimizer(int surfaceLevel=MATCH_POINT_CLOUD_SUBDIVISIONS);
    void doVoxelIterations(int numIterations, bool preview = false);
	double getLastVoxelMatchValue();
    double getLastVoxelPenaltyValue();
    double getLastVoxelObjectiveFunctionValue();
    void setVoxelPenaltyWeight(double w);
    void setVoxelConstraintsPenaltyWeight(double w);
    bool voxelStageOptimizerIsFinished();

    void setIntensityWindow(double min, double max);
    void getIntensityWindow(double & min, double & max);

	// Multi-figure construction
    bool attachSubfigure(int parentId, int childId,
		M3DFigureTreeNode::SubfigureAttachment_t mode = M3DFigureTreeNode::PROTRUDE);
    bool detachSubfigure(int figureId);

    // Interpenetration avoidance
    void setPartialDrawingMode(bool mode) { drawing_partial_surface = mode; }
    bool getPartialDrawingMode() { return drawing_partial_surface; }
    void togglePartialRendering(bool mode) { drawing_full_surface = mode; }
    bool getPartialSurfaceMode() { return drawing_full_surface; }

    // Subfigure transformations test functions
    void initializeSubfigureTransformation();
    void translateSubfigure(double u, double v);
    void rotateSubfigure(double angle);
    void scaleSubfigure(double scale);
    void elongateSubfigure(double elongation);
    void hingeSubfigure(double angle);

    // Blurs the current image
    void gaussBlur();

#ifdef BINARY
    	// isotropicly resample image to smallest voxel spacing (or to user-spec'd spacing)...
	void resampleImage();
	// resample as above, then pad up to maximum extent
	void resampleImageUnitCubic(float maxExtent = 0);
	// Both overwrite the image file
#endif

	// load image file and return max extent of any dimension.
	double calcMaxExtent(const char * imageFile);

#ifdef BINARY
	// Creates a distance map off the current image
    bool distMap(const char *filename);
#endif

    // Outputs a histogram of the current image (for plotting in gnuplot, etc.)
    void outputHistogram(const char * filename);

    void recordTransformation();
    SimilarityTransform3D * getRecordedTransformation() { return recordedTransformation; }
    void stopRecordingTransformation();
    void setInitialTransformInMatch(SimilarityTransform3D * transform) {
		match->setInitialTransform(transform);
	}
    void applyRecordedTransformation();
    void resetTransformation();

	// PGA methods

    bool readPGAFile(const char * filename);
#ifndef BINARY
	void initPGA();
#endif
    bool pgaDeform(std::vector<double> & vals);

	bool pgaDeformAtom(std::vector<double> & vals);

	bool havePGA();
	bool haveFigurePGA();
	bool havePrimitivePGA();
	void resetResiduePGAStage();
	int getResiduePGAStage();
	char * getResiduePGAStageName();
	void setResiduePGAStage(int depth);
	bool applyMeanResidue();
#ifndef BINARY
	void clearPGA();

#ifdef PCA_VISUALIZER
	//PCA methods
	bool loadPCA_File(const char *filename);
	bool initPCA_Deform();
	bool pcaDeform(std::vector<float> &vals);
#endif
#endif

    // Model slide show
    void startModelTileSlideShow();
    void addModelSlideShow(const char * filename);
    void removeModelSlideShow(int modelNumber);
    void moveUpModelSlideShow(int modelNumber);
    void moveDownModelSlideShow(int modelNumber);
    void clearModelSlideShow();
    void changeModelSlideShow(int number);
    void quitModelSlideShow();
	int getInSlideShow() { return inSlideShow;}

    // Tile slide show
    void addTileSlideShow(const char * filename);
    void removeTileSlideShow(int TileNumber);
    void moveUpTileSlideShow(int TileNumber);
    void moveDownTileSlideShow(int TileNumber);
    void clearTileSlideShow();
    void changeTileSlideShow(int number);
    void quitTileSlideShow();

#ifndef BINARY
    double testConstraintsPenalty(int figureId);
	void testSeurat(int figureId, int level);
#endif

	// Export undo list parameters
	int undoListSize() { return undoList.size(); }
	int undoListMaxSize() { return undoList.maxSize(); }

	int machineByteOrder() const;

	// To check whether the entire sub-figure tree is selected or not
	bool isFigureSubTreeSelected(M3DFigureTreeNode * node);

	void setMatchReferenceObject(M3DObject * _matchReferenceObject) {
		matchReferenceObject = _matchReferenceObject;
	}

#ifndef BINARY
	void setCurrentPGAOptimizationStep(int step) { currentPGAOptimizationStep = step; }
#endif

#ifdef LM_METHOD_OBJ
	// each figure represents a landmark's position and the first atom's
	// radius is interpreted as the landmark's tolerance in the IMAGE.
	// These are not true figures; see landmark comment at top of file.
    M3DObject * landmarkObject;
	void binaryPabloListLandmarks();
	bool binaryPabloLoadLandmarkModel(const char * filename);

    bool readLandmarkModel(const char * filename);
    bool readAppendLandmarkModel(const char * filename);
	void mergeLandmarkModel(M3DObject * newObject);

	// atom selection via mouse is broken when the landmark model exists.
	// LATER: fix atom selection
	void newLandmarkModel();	// forget the landmark model.

	// search either an IMAGE or MODEL landmark for a named landmark,
	// and set lmPosit to it's model coordinates. return if found
	int findLandmarkPositByName(Vector3D &lmPosit, M3DObject * obj, int figureId,
		const char * landmarkName, int objType);

#endif

#ifdef BINARY

	Registry optimRegistry;
	// thin plate spline addition.
	int runThinPlateSpline(Registry &registry);
	bool runBinaryPablo(Registry & script);
	Image3D * getResampedBinaryScript(Registry & script,
		char *binImageFileKey, char * resBinImageFileKey);
	ImageDistanceMap * getDDMScript(Registry & script,
		Image3D * image, char * filenameIn);
	void binaryPabloNextOptimizationStage(int surfaceLevel);
	bool binaryPabloFigureStage(int maxNumIterations,
		int doScale, // bitFlags 1=scale
		const char * pgaFilename);
	bool binaryPabloAtomStage(int maxNumIterations);
	bool binaryPabloSpokeStage(int maxNumIterations);

	ImageDistanceMap * getImageDistMapPointer() {
		ImageDistanceMap * dm;
		if (match != NULL) {
			dm = match->getImageDistMapPointer();
			if (dm != NULL) return dm;
		}
		return NULL;
	}

	Match * getMatchPtr() {
		return ((match != NULL) ? match : NULL);
	}

	void setDistMapPtr(ImageDistanceMap * im);
	ImageDistanceMap * getDistMapPtr();
	bool loadDistMap(const char * filename);
	bool loadAntialiasedDistMap(const char * filename);

#else	/* BINARY */
	bool runPablo(Registry & script);
	bool initializeOptimizer(Registry & script);
	void figureStage(bool mainFigureStage, Registry & script);
	void atomStage(Registry & script);
	void spokeStage(Registry & script);
	void modelStage(bool modelStage, Registry & script);
	void subfigureStage(Registry & script);
	void voxelStage(Registry & script);
#endif	/* BINARY */

    bool align(const char * pointsFile, const char * modelFile = NULL,
		const char * imageFile = NULL, double closeness = 0.04,
		bool pgaOnly = false, bool scaleOnly = false, bool scaleTransOnly = false,
		const char * outputTransform = NULL);

	M3DPGAStats * getPGAStats() const { return pgaStats; }
	M3DPGAPrimitiveStats * getPrimPGAStats() const { return atomPgaStats; }
	// TODO:
	//M3DPGASpokeStats * getSpokePGAStats() const { return spokePGAStats; }

	bool isAdaptive();

	void setGuiCallbacks(void (* guiCB)(int), void (* guiStageCB)(OptimizerStage));

	// after modifying clip planes, call redraw()
	clipPlane* getClipPlanes() { return clips; }

#ifdef AE2_BUILD
    void addAnastruct(Anastruct3D * a);
    void clearAnastructList();
    bool convertAnastructCoords();  // Call after setAnastruct(); P3DControl must have an image
    AnastructList * getAnastructList() const { return anaList; }
#endif


protected:

    M3DObject * object;
	// The similarity transform conceptually goes with the model
	SimilarityTransform3D * objectSimTransform;

    ImagePlanes imagePlanes;
    bool imagePlanesChanged;
	bool selectionChanged;	// Indicates to the optimizer that the model changed
    TileSet * tileSet;
    void * auxTileSet;	// Pointer to a class BYU object

    // List of undo operations
    P3DUndoList undoList;

    M3DObject * objectCopyBuffer;
	bool shiftOnPaste;

    int markedPrimitiveId;

    EditingMode editingMode;
    bool selection_allowed;

	// Data for the templates filename
	char * templates_filename;

    // Data for transformation recording
    SimilarityTransform3D * recordedTransformation;
    bool recordingTransformation;

#ifndef BINARY
	// Similar to transformation recording, here we record the hand placement
	// similarity.  This is so PGA deformations are applied to a non-sim-transformed
	// object.  To be applied to the unmodified, loaded object after the PGA deformations.
	SimilarityTransform3D * handPlacementTransform;
#endif

    // Flag set if the model has changed since last save operation
    bool modelIsDirty;
	bool modelWasMoved;	// Kind of change made

	bool worldWasApplied;	// True if model was mapped to image's coordinate system

    // File names of the current model and image
    char * modelFileName;
	char * pgaModelFileName;
//#ifdef BINARY
//	char * lmrkFileName;
//#endif
    char * imageFileName;

    M3DPrimitive * interpolatedPrimitive;

#ifdef AE2_BUILD
    AnastructList * anaList;
#endif

    bool inOperation;
	bool blurred;

    // Optimizer variables
    OptimizerStage currentOptimizerStage;
    Match * match;
    M3DRegistrationOptimizer * registrationOptimizer;
    M3DMainFigureOptimizer * mainFigureOptimizer;
    M3DSubfigureOptimizer * subfigureOptimizer;
    M3DDeformationOptimizer * deformationOptimizer;
    M3DSpokeOptimizer * spokeOptimizer;
	M3DVoxelOptimizer * voxelStageOptimizer;

	double registrationPenaltyWeight;

#ifndef BINARY
	// Arrays for storing penalties used by the optimizer classes
	double * penaltyWeights;
	gpTune_t * penaltyNames;
	int numPenaltyWeights;
#endif

	double deformationObjFunc;
	double deformationGeomPenalty;

	double spokeObjFunc;
	double spokeGeomPenalty;

    int optimizerTreeIndex;
    int optimizerFigureId;
    int * treeFigures;
    int treeFiguresCount;

	bool optimizeWithWindowing;

	/* CL: currentPGAOptimizationStep keeps track of the current optimization step,
		   which is used by pgaStats to determine which PGA statistics to use.
	*/
	int currentPGAOptimizationStep;

	// Visualizer
	bool colorSurface;

	// Refresh the GUI
	void (* updateGui)(int iteration);	// Update the GUI after each optimizer iteration
	void (* updateGuiStage)(OptimizerStage);	// Update the GUI when a stage finishes

    // Inter-penetration avoidance
    bool drawing_partial_surface;
    bool drawing_full_surface;

    // Subfigure Transformation test transform
    M3DSubfigureTransformation subfigureTransform;

    // Elongater test
    M3DFigureElongater * elongater;
    bool initializeElongater();

	// Only used by recordTransformation(), resetTransformation(), and related functions
    M3DObject * referenceObject;

	// For use in the bootstrap - reusing the same reference object each iteration
	M3DObject * matchReferenceObject;

    // Variables used when picking/editing landmarks
	// COMMENT OUT the code which uses these runtimeLM values SOON -gst
	int numImageLandmarks;		// Current number of specified image landmarks
	vector<imageLandmark_t> imageLandmarks;

    // List of model objects and tile filenames in slide show
    vector<M3DObject *> modelSlideShow;
    vector<TileSet *> tileSlideShow;
    M3DObject * slideShowSavedObject;
    TileSet * slideShowSavedTileSet;
    EditingMode slideShowSavedEditingMode;
    bool inSlideShow;

	bool tileSetMobile;
	bool tileSetMobileFigure[6];	// allow figure to move with mouse;
		// hardcode max # of figures

	void merge(M3DObject * newObject);

	// Pointers to the PGA statistic objects in object or PGAObject 
	M3DPGAStats * pgaStats;
	M3DPGAPrimitiveStats * atomPgaStats;
	M3DPGASpokeStats * spokePgaStats;
	bool attachedStats;    // Allocation flag used with pgaStats and atomPgaStats

	int residuePGAStage;

#ifndef BINARY
#ifdef PCA_VISUALIZER
	PCA * pca;
#endif
#endif

#ifdef BINARY
	ImageDistanceMap * dist_map;
#else
	SurfacePatchEnsemble * spe;
	MultiPatch * mp;
	DQFMatch * dqf;
#endif

	clipPlane clips[6];	// openGL must support at least 6 planes, so use 6

private:

    void verifyOrientationMatches(M3DObject * model);
	void reorderSubModel(M3DFigureTreeNode * node, int * newFigureNums, int & n);
	bool testReorderSubModel(M3DFigureTreeNode * node, int & n);
    void listTreeFigures(M3DFigureTreeNode * node);
    M3DFigureTreeNode * initializeOptimizerTree();
	int currentFigureIndex;

};


#endif

