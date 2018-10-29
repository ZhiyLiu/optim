#ifndef P3D_USER_INTERFACE_CALLBACK_H
#define P3D_USER_INTERFACE_CALLBACK_H

#include "DistanceMap3D.h"
#include "P3DControl.h"
#include "M3DObjectSurfaceRenderer.h"
#include "ControlParms.h"
#include "P3DView.h"

class P3DUserInterface;
class P3DControl;
class SelectedPartialFigures;

// These are geometric "things" on which we can center the view.
// Use them in combinations to identify
// the thing's type and thing components to center upon.
// FIG_MASK is a set of bits interpreted as an integer.
// Used for P3DUserInterfaceCallback::centerThing

#define CEN_MODEL      0x10000000	// thing = object
#define CEN_TILE       0x20000000	// thing = tile
#define CEN_SPACE      0x40000000	// thing = unit cube
#define CEN_ILM        0x80000000	// thing = image landmark;
									// use FIG_MASK and ILM_MASK

#define CEN_ALL        0x01000000	// use all/whole thing
#define CEN_FIG        0x02000000	// use FIG_MASK to select a piece of thing
#define CEN_MARK       0x04000000	// use mark within thing
#define CEN_SELECT     0x08000000	// use selected part of thing

#define CEN_FIG_MASK   0x000000ff	// use this figure within model/tileset
#define CEN_FIG_SHIFT  0			// (FIG_MASK >> FIG_SHIFT) is the encoded int
#define CEN_ILM_MASK   0x0000ff00	// use this image landmark index within FIG_MASK
#define CEN_ILM_SHIFT  8			// (ILM_MASK >> ILM_SHIFT) is the encoded int


class P3DUserInterfaceCallback
{
public:

	enum cut_plane_mode_t {BOTH_WINS = 0, SEPARATE_WIN = 1, MODEL_WIN = 2};

	enum display_tab_t {GENERAL_TAB, SURFACE_TAB,
		IMAGE_TAB, TILESET_TAB};

	enum slice_t {SAGGITAL = 1, CORONAL = 2, AXIAL = 3};

    P3DUserInterfaceCallback();
    ~P3DUserInterfaceCallback();
    void closeWindowsAndExit();

	// Functions for interactive display of scripted optimization
    void setUserInterface(P3DUserInterface * ui) { userInterface = ui; }
	void startInteractive(int matchMethod);
	void stopInteractive();
	void updateOptimizerWindow(OptimizerStage nextStageId);
	void startNextInteractiveStage();
	void showInteractiveResults(int iteration);
	void continuousOptimization();
	void enableOptimzerButtons(OptimizerStage stageId);
	void disableOptimzerButtons(OptimizerStage stageId);

    void set_mouse_button_order(bool left_hand_mouse);
    bool usingLeftHandedMouse() { return mouse_left_handed; }

    void init();
    void show();
	bool exiting() { return halting; }

	P3DControl * getControl() const;

	void adjustInterfaceForBinary();

	void binaryPabloAlignModel();

	// Binary Pablo script functions
	void loadBinaryPabloScript();

    // File Menu commands
    void newModel();
    void loadModel();
    bool loadModel(const char * filename);
    void loadOldModel();
	void loadSimTransform();
	bool loadSimTransform(const char * filename);
    void newLandmarkModel();
    void loadLandmarkModel();
    void loadLandmarkModel(const char * filename);
    void saveModel();
    void saveModelAs();
	void clearReorderWindowMode() { reorderWindowOpen = false; }
	void cancelSave() { saveCanceled = true; }

    void exportToImage();
    void exportDistanceMap();
    void exportToBYUTiles();

    void loadImage(bool mapToWorld = false, bool stacked = false);
    void loadImage(const char * filename, bool mapToWorld = false,
		bool stacked = false, bool silent = false);
    void discardImage();
	void saveBlur();

    void loadTileSet();
	void generateTileSet();
	void saveTileSetAs();

	void loadAuxTileSet();
	void viewAuxTileSet(int val);

	const char * askSingleFilename(control_keys_t startDir,
		char * message, char * pattern, bool ignore = false);

    // Edit Menu commands
    void undo();
    void redo();
    void addQuadFigure();
    void initializeAddQuadDlg();
    void removeSelectedFigures();
    void enableEditingWindows();
    void disableEditingWindows();

	void editModelProperties();
	void reorderModel();
	void embedModelFigureNameEditor(M3DObject * object);
	void editModelFigureNumber();
	void editFigureColor();
	void editModelFigureBlendExtent();
	void editModelFigureBlendAmount();
	void editModelPolarity();
	void editModelFigureTolerance();
	void editModelReblend();
	void editModelAttachmentMode(int mode);
	void finishEditModel();

    void selectAll();
    void deselectAll();
    void toggleAll();
    void setSelectionType(P3DSelectionType type);

    void copy();
    void paste();
    void mirror();

	void applyLoadedTransformation();

	// Landmarks
	void startEditLandmarks();
	void clearSurfaceLandmarks();
	void clearImageLandmarks();
	void setLandmarksSurfaceNumber(int numLandmarks);
	void setLandmarksImageNumber(int numImageLandmarks);
	void updateLandmarkEditorWindow();
	void setLandmarksEditMode(int mode);
	void editLandmarksNext();
	void editLandmarksPrevious();
	void nextImageLandmark();
	void markImageLandmark(int landmark, bool onOff);
	void editLandmarkDelete();
	void editLandmarkGoto();
	int getLandmarksEditMode() { return landmarksEditMode; }
	int getCurrentLandmarkIndex() { return editingLandmark; }
	void finishEditLandmarks();
	void applyLandmarks();

    void setMainWindowLabel();
    void setFigureColor(void * currentColor, enum Fl_Color color_code);
    void setFigureColor(void * currentColor);

    void update();          // Causes everything to be rendered
    void updateObject();    // Causes everything except the surfaces to be rendered
    void updateSurface();   // Causes everything except the implicit surface to be rendered
    void redraw();          // Refreshes the modelView window without rendering anything

    // Viewing controls
    void toggleUnitCube(int val);
    void toggleStdAxes(int val);
    void toggleRocking(int val);
    void resetView();

    void setAxialView();
    void setSaggitalView();
    void setCoronalView();

    void rotateLeft90();
    void rotateRight90();
    void rotateHorizontal180();
    void rotateVertical180();

	void setCenterOfRotation();
	void setCenterThing();

    void setSurfaceStyle(M3DSurfaceStyle style);
    void setSurfaceLevel(int val);
	void setSurfaceLineWidth(int width);
    void toggleBlendSubfigures(int val);
	void toggleDrawBoundary(int val);
	void toggleDrawFigureNames(int val);

    void viewXImagePlane(int val);
    void viewYImagePlane(int val);
    void viewZImagePlane(int val);
    void setXImagePlanePosition();
    void setYImagePlanePosition();
    void setZImagePlanePosition();

	void setImagePlane(enum slice_t slice, int val);

    void viewTileSet(int val);
    void setTileAlpha(double val);
	void useDefaultColor(int val);
	void tileSetSolid(int val);
	void tileSetBright(int val);
	void tileSetDark(int val);
	void tileSetMobile(int val);

	void tileMobileToggleFigure(int figId);

    void setInterfaceModelPath(const char * path);
    void setInterfaceImagePath(const char * path);
    void setInterfaceTilesetPath(const char * path);

	void raiseDisplayTab(enum display_tab_t tab);

	void clipToggleEnable(int plane);
	void clipAValuatorChanged(int plane);
	void clipBValuatorChanged(int plane);
	void clipCValuatorChanged(int plane);
	void clipDValuatorChanged(int plane);
	void clipToggleMask(int plane, int mask);
	void clipToggleTile      (int plane) { clipToggleMask(plane, CLIP_MASK_TILE  );};
	void clipToggleModel     (int plane) { clipToggleMask(plane, CLIP_MASK_MODEL );};
	void clipToggleImagePlane(int plane) { clipToggleMask(plane, CLIP_MASK_IMAGE );};
	void clipToggleScene     (int plane) { clipToggleMask(plane, CLIP_MASK_CAMERA);};
	void clipToggleView      (int plane) { clipToggleMask(plane, CLIP_MASK_VIEW);};
	void clipCaptureView(int plane);
	void clipCaptureView0() {clipCaptureView(0);};
	void clipCaptureView1() {clipCaptureView(1);};
	void clipCaptureView2() {clipCaptureView(2);};
	void clipCaptureView3() {clipCaptureView(3);};
	void clipCaptureView4() {clipCaptureView(4);};
	void clipCaptureView5() {clipCaptureView(5);};
	void clipFlipView(int plane);
	void clipFlipView0() {clipFlipView(0);};
	void clipFlipView1() {clipFlipView(1);};
	void clipFlipView2() {clipFlipView(2);};
	void clipFlipView3() {clipFlipView(3);};
	void clipFlipView4() {clipFlipView(4);};
	void clipFlipView5() {clipFlipView(5);};
	void setClips(int plane, float A, float B, float C, float D);

    // Visibility controls
    void initializeInterfaceVisibility();
    void visibilityAllOn();
    void visibilityAllOff();
	void visibilityAllInverted();
    void atomVectorsVisible(int val);
    void meshConnectorsVisible(int val);
	void atomsVisible(int val);
    void setFigureVisibility(int buttonID, int val);
    void setInterfaceVisibilityCount(int numFigures);
    void updateVisibilityFiguresDisplay(bool reset);
    void setInterfaceFigureName(int figureID, const char * name,
		M3DObject * object, bool isRoot, bool vis);
	void constraintsVisible(int val);
	void landmarksVisible(int val);
	void setVisibilityFigureRange();

    // Preferences setting
    void saveWindowPositions();
    void saveWindowsOpen();
    void showPreferencesEditor();
    void setPrefsRememberPositions(int toggle);
    void setPrefsRememberOpen(int toggle);
	void setPrefsIconifyWindow(int val);
	void setPrefsStandardAxes(int toggle);
    void setPrefsSmoothImages(int toggle);
	void setPrefsRangeDisplayed(int toggle);
	void setPrefsScaleInputImage(int toggle);
	void setPrefsCoordsDisplayed(int toggle);
    void setPrefsByteOrder(int order);
    void setPrefsCompressImages(int toggle);
	void setPrefsConvertImages(int toggle);
	void setPrefsImageFormat(int val);
	void setPrefsConvertFormat(int toggle);
	void setPrefsDrawBoundary(int toggle);
	void setPrefsSliceDefault(int slice, int val);
    void setPrefsCutPlaneBoundary(int toggle);
	void setPrefsCutPlaneBdryWidth(int w);
	void incrementPrefsCutPlaneMode();
    void setPrefsRockingAngle(double angle);
    void setPrefsRockingIncrement(double angle);
    void setPrefsShowAtomVectors(int toggle);
	void setPrefsShowAtoms(int toggle);
	void setPrefsUseFigureColors(int toggle);
    void setPrefsExtraAtomVectors(int toggle);
    void setPrefsAtomVectorsWidth(int width);
    void setPrefsAtomVectorsType(int toggle);
    void setPrefsAtomBVectorsType(int val);
    void setPrefsShowMesh(int toggle);
    void setPrefsMeshType(int toggle);
    void setPrefsMeshColor(enum Fl_Color color_code);
    void setPrefsMeshColor();
    void setPrefsMeshConnectorsWidth(int width);
    void setPrefsSurfaceStyle(M3DSurfaceStyle style);
    void setPrefsSurfaceLevel(int lvl);
	void setPrefsSurfaceLineWidth(int width);
	void setPrefsUseTwoLights(int toggle);
    void setPrefsPartialSurfaceStyle(M3DSurfaceStyle style);
    void setPrefsPartialSurfaceLevel(int val);
    void setPrefsPartialRendering(int val);
	void setPrefsShowConstraints(int toggle);
	void setPrefsShowLandmarks(int toggle);
	void setPrefsLandmarksNarrowWidth(int width);
	void setPrefsLandmarksWideWidth(int width);
	void setPrefsLandmarksColor(enum Fl_Color color_code);
	void setPrefsLandmarksColor();
    void setPrefsBackgroundColor(enum Fl_Color color_code);
    void setPrefsBackgroundColor();
    void setPrefsPartialSurfaceColor(enum Fl_Color color_code);
    void setPrefsPartialSurfaceColor();
	void setUndoLength(int undoLength);
	void setPrefsMouseButtonOrder(bool left_handed);
	void setPrefsReorderWindow(int alwaysNeverPrompt);
	void setPrefsSurfaceSmoothness(int tolerance);
    void setPrefsTilesColor(enum Fl_Color color_code);
    void setPrefsTilesColor();
	void setPrefsBYUOutputType(int type);
	void setPrefsVerbosity(int level);
	void setPrefsSimTransOutputType();

	// Information Window
	void showAboutPabloWindow();
	void initAboutPabloWindow();
	void setPrefsAboutWindow(bool show);

    // Cut Planes Viewing Control
    void viewBBPerpCutPlane();
    void viewBNCutPlane();
    void viewBPerpNCutPlane();
    void viewBPerpY0CutPlane();
    void viewBPerpY1CutPlane();
    void viewInvolutesCutPlane();
	void hideCutPlanes(cut_plane_t plane);
    P3DCutPlaneView * getView(cut_plane_t planeType) const;
	bool redrawCutPlane(cut_plane_t planeType);
	void incrementCutPlanesMode(cut_plane_t plane);
	cut_plane_mode_t getCutPlaneMode(cut_plane_t plane) const {
		return (cut_plane_mode_t) cutPlaneModes[(int) plane];
	}
	void initCutPlanesMode();

    void showPrimitiveEditor();
    bool updateAtomEditorWindow();
    void setTheta();
    void setElongation();
	void updateElongationWindow();

    void interpolatePrimitive(double u, double v);
    void destroyInterpolatedPrimitive();

    // Regularizer
    void regularizeModel(int iterations, double stepsize);

    // Optimizer Controls
    OptimizerStage nextOptimizationStage(bool interactive = false);
    void cancelOptimization();
#ifndef BINARY
    bool initializeOptimizer();
    void setOptimizationFigureName(int figureId);
	void zeroValuators(OptimizerStage s);
	void enableOptimizerWindow();
	bool disableOptimizerWindow();
#endif
	void previewOptimizationStage();
	void setOptimizerIntensityWindowing();

		// These are used for the ensemble and object stages
    void startModelStage(bool modelStage);
    void stopModelStage();
    void resetModelStage(bool modelStage);
    void setModelStagePenaltyWeight(bool modelStage);
#ifndef BINARY
    void clearModelStage(bool modelStage);
#endif

		// These are used for the Main Figure and Figure stages
    void startMainFigureStage(bool mainFigureStage);
    void stopMainFigureStage();
    void setMainFigureStagePenaltyWeight(bool mainFigureStage, MatchType type);
    void setMainFigureStageConstraintsPenaltyWeight(bool mainFigureStage,
		MatchType type);
#ifndef BINARY
    void clearMainFigureStage(bool mainFigureStage);
#endif

    void startSubfigureStage();
    void stopSubfigureStage();
    void setSubfigureStagePenaltyWeight();
    void setSubfigureStageConstraintsPenaltyWeight();
#ifndef BINARY
    void clearSubfigureStage();
#endif

    void startPrimitiveStage();
    void stopPrimitiveStage();
    void setPrimitiveStagePenaltyWeight();
    void setPrimitiveStageConstraintsPenaltyWeight();
    void setPrimitiveStageNeighborPenaltyWeight();
#ifndef BINARY
    void clearPrimitiveStage();
#endif

	void startVoxelStage();
    void stopVoxelStage();
    void setVoxelStagePenaltyWeight();
    void setVoxelStageConstraintsPenaltyWeight();
#ifndef BINARY
    void clearVoxelStage();
#endif

#ifndef BINARY
	void initializeOptimizerSettings();
#endif
    void setPenaltyWeights();
    void revertPenaltyWeights();
	void loadPGAFile();
	void loadPCAFile();

    void intensityWindowSliderChanged();
    void intensityWindowValuatorChanged();

    void initializeAttachSubfigureDialog();
    void advanceAttachSubfigureDialog();
	void setAttachSubfigureMode(int mode);
	M3DFigureTreeNode::SubfigureAttachment_t getAttachSubfigureMode();
    void cancelAttachSubfigureDialog();

    void detachSubfigure();

    // Elongater Test Interface
    void initializeElongater();
    void elongateRow();
    void elongateColumn();
    void elongateFinish();

    // Subfigure Transformation Test Interface
    void initializeSubfigureTransformation();
    void translateSubfigure();
    void rotateSubfigure();
    void scaleSubfigure();
    void elongateSubfigure();
    void hingeSubfigure();

    // Record Transformation Test Interface
    void startRecording();
    void endRecording();
    void resetObject();
    void testTransformation();

    // PGA Deformation Interface
	void pgaInitDeformationWindow();
    void resetPGADeformation();
    void pgaDeform();
	void pgaNextScale();
	void pgaZero();
	void pgaApplyMeanRes();
	void initPGADeform();
	void pgaDeformationWindowDone();
	void enableDeformationWindow();
	void disableDeformationWindow();

	// dibyendu
	// CPNS Deformation Interface	
	void cpnsInitDeformationWindow() ;
	void resetCPNSDeformation();
	void cpnsDeform() ;
	void cpnsZero();
	void initCPNSDeform() ;
	void cpnsDeformationWindowDone();
	void enableCPNSDeformationWindow();
	void disableCPNSDeformationWindow();

    // Zhiyuan Liu
    void optLength();//optimize srep using NEWUOA
    void optAngle();//optimize srep using NEWUOA
    void optAll(); // optimize both length and angle by NEWUOA

	// PCA Deformation Interface
	void initPCADeform();
	void pcaDeform();

    // Model Slide Show Interface
    void startModelTileSlideShow();
    void addModelsToSlideShow();
    void removeModelSlideShow();
    void moveUpModelSlideShow();
    void moveDownModelSlideShow();
    void clearModelSlideShow();
    void changeModelSlideShow();
    void selectModelSlideShow();
    void quitModelSlideShow();

    // Tile Slide Show Interface
    void addTilesToSlideShow();
    void removeTileSlideShow();
    void moveUpTileSlideShow();
    void moveDownTileSlideShow();
    void clearTileSlideShow();
    void changeTileSlideShow();
    void selectTileSlideShow();
    void quitTileSlideShow();

    void outputPOVRay();

    void gaussBlur();

	void distMapAs();
	void distMap(const char *filename);
	bool loadDistMap(const char * filename, bool silent = false);
	void loadDistMap();

    void outputHistogram();
    void initializeTestConstraints();
    void testConstraintsPenalty();
	void testSeurat();
	void printModel(int worldCoords);
	void printImage();

    void makeMovie();
    void screenDump();
	void toggleBlendingRegionDisplay();
	void toggleTileType();

	// Functions used with the (Match) Surface Visualizer
	void importSurfaceValues();
	void discardMatchSurfaces();
	void importSurfaceColorMap();
	void openMatchSurfacesWindow();
	void toggleMatchSurfaceDisplay();
	void setMatchSurfaceColorRange();
	void setMatchSurfaceStyle(M3DSurfaceStyle style);
	void setMatchSurfaceSerialName(char * name);
	void resetMatchSurfaceSerialNumber();
	void matchSurfaceSerialNumberChanged(int val);
	void chooseMatchSurfaceColormap();
	void tieSurfaceToImageSlider(int val);
	void togglePlainMatchSurfaces(int yesNo);

    bool askToSaveModel();

    // Interpenetration avoidance
    bool constraintsWindowSetup();
    void setInterObjectDistance();
    bool initializeConstraints();
    void showPenetrationArea();
    void setInterObjectConstraint();
    void exitConstraintsMode();
 //   void toggleSelectedAtoms(int button);
    double getInterObjectDistance();
    void setPartialSurfaceStyle(M3DSurfaceStyle style);
    void setPartialSurfaceLevel(int val);
    SelectedPartialFigures * getPartialFigures() { return selectedPartialFigures; }
    void togglePartialRendering(int button);
    void clearPartialConstraints();
    void clearAllConstraints();
    void setInterObjectFigure();
    void setDistanceSliderFromFigureSlider();
    int getInterObjectFigureNumber() { return currentFigure; }
    void updateInterObjectFigureSlider(int figureNumber);

    bool imagePlaneSliderFinishedMoving();

	// Special function for testing
	void secretAction();

	void exportSeuratUVPhi();

#ifdef OPTIMIZATION_VISUALIZER
	void setOptVisualizerUI(void * oF_UI) { ovUI = oF_UI; }
#endif
	void openOptimizationVisualizer();

#ifdef BINARY
	bool readDistanceVectorsFile(const char * filename);
#endif

private:

    P3DUserInterface * userInterface;
    char * mainWindowLabel;
#ifdef OPTIMIZATION_VISUALIZER
	void * ovUI;
#endif

	// Elongation slider state
	struct sliderState {
		void * widget;
		void (* pushAction)(void *, bool);
		void (* releaseAction)(void *, bool);
		void * self;
		bool rowCol;
		bool pushed;
	} rowSliderState, columnSliderState;
	static void releaseDispatcher(void * self, bool rowCol);
	void resetRowElongation();
	void resetColumnElongation();
    char ELONGATION_STEP1[46];

    // Optimizer Control variables
    bool optimizerIsRunning;
	bool optimizerWindowInUse;
    bool optimizerStarted;
    bool optimizerStopRequested;
    bool resettingOptimizer;
    int currentMaskFigureId;
	bool previewMode;
	bool lastModeWasPreview;
	char * profileFilename;
	bool interactiveOpt;
	bool doNextInteractiveOpt;
	bool continueInteractiveOpt;
	OptimizerStage currentInteractiveStage;
	OptimizerStage nextInteractiveStage;
	MatchType optimizationType;

	// PGA
    char * pgaFilename;
	void finishPGALoading(const char * filename);
	void activatePGASliders();
	void deactivatePGASliders();

	// CPNS
	void activateCPNSSliders();
	void deactivateCPNSSliders();

	// Cut planes control
	bool cutPlanesShown[6];
    float cut_plane_intensity;
	int cutPlaneModes[6];
	char * cutPlaneModeLabels[3];
	int cutPlaneModePreference;
	void setCutPlanesMode(int p, int value);

    // Attach subfigure dialog variables
    int subfigureAttachScreenId;
    int parentFigureId;
    int childFigureId;
	M3DFigureTreeNode::SubfigureAttachment_t attachmentMode;
    char ATTACH_STEP2[54];
    char ATTACH_STEP3[53];
    char ATTACH_STEP4[78];

	// Atom editor dialog
    char ATOM_EDITOR_STEP1[33];

    int imagePlaneSliderMotion;

    void setIntensityWindow(int min, double relativeMin,
		int max, double relativeMax);

	bool changedBlend;

    // List of figures and other information for interpenetration avoidance
    SelectedPartialFigures * selectedPartialFigures;
    bool partialSurfaceSliderChanged;
    int currentFigure;
    M3DFigure * markedFigurePtr;
    DistanceMap3D map;
    float partialSurfaceColor[3];

    // Copies of file names needed for Display Control/General Tab
    char * displayModelPath;
    char * displayTilesetPath;
    char * displayImagePath;

	// Slice, position and intensity labels needed for Display Control/Image Tab
	char xImagePlaneSliceStr[4];
	char yImagePlaneSliceStr[4];
	char zImagePlaneSliceStr[4];
	char xImagePlanePosnStr[12];
	char yImagePlanePosnStr[12];
	char zImagePlanePosnStr[12];
	char intensityMinStr[12];
	char intensityMaxStr[12];

	// Preferences
    bool recordWindowPositions;
    bool recordWindowsOpen;
	bool mouse_left_handed;
	bool displayModelCoords;
	bool displayIntensityScale;

	// Model names
	class Flvt_Edit * modelNameEditor;

	// Landmarks
	int landmarksEditMode;
	bool priorLandmarkVisibility;
	int editingLandmark;	// Serial number of landmark currently being edited
	class Flvt_Edit * landmarkNameEditor;
	int landmarkCaret;		// Line number, [-1, 5], of marked landmark in editing window 
	void embedLandmarkNameEditor(M3DObject * object);
	void * generateLandmarkEditorCells(M3DObject * object, int numLandmarks);
	static void editLandmarkName(void *, int);
	Fl_Callback * vscrollbar_cb;
	static void landmarkEditScrollbarCallback(Fl_Widget * widget, void * user_data);
	void setLandmarkSlider();
	void initializeSurfaceLandmarkEditing();
	void initializeImageLandmarkEditing();
	void setLandmarkCaret(int row, bool adjustSlider);

	// Visibility
	char visibilityFigureNumStrs[7][3];
	int visibilityFirstFigureNum;

	bool reorderWindowOpen;
	bool saveCanceled;

	int usingModalWindow;
	bool halting;

	int brightTiles;

#ifdef BINARY
	DistanceVectorList * distanceVectorList;	// Diagnostic distance vector list
#endif
	int centerThing;	// center on this model, tile, image, etc.
						// see CEN_* defines

};



#endif

