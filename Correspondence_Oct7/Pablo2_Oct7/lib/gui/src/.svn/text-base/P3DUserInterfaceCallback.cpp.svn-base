#include <time.h>
#include <FL/Fl_Box.H>
#include <FL/Fl_File_Chooser.H>
#include <FL/Flv_Table.H>
#include <FL/Flvt_Edit.H>
#include <FL/Flvt_Edit_Cell.H>
#include <iostream>
#include <iomanip>
#ifdef BINARY
#include <iterator>
#endif
#include "globalBuildParms.h"
#include "P3DUserInterfaceCallback.h"
#include "P3DUserInterface.h"
#include "M3DFigure.h"
#include "M3DEndPrimitive.h"
#include "M3DObjectFile.h"
#ifdef BINARY
#include <FL/Fl_Pixmap.H>
#include <FL/Fl_Tiled_Image.H>
#include "bp_splash.xpm"
#include "DistanceVectorList.h"
#include "utility.h"
#endif
#include "make_windows.h"
#include "ControlParms.h"
#include "renderDefinitions.h"
#include "LandmarkDeformation.h"
#include "button_labels.h"
#include "pablo_version.h"
#include "M3DtoPovray.h"
#include "SelectedPartialFigures.h"
#include "CutPlaneGlobals.h"
#include "menu.h"

#include <typeinfo>

#include "Matrix4D.h"

#ifdef OPTIMIZATION_VISUALIZER
#include "OptVisualizerUI.h"
#endif

#ifdef SURFACE_VISUALIZER
#include "M3DObjectSurfaceVisualizer.h"
#include "SurfaceColorMap.h"
#endif

//#define DEBUG

#ifdef BINARY

#define BINARY_PABLO_REPORT

// The change history - SHOULD moved to pablo.cpp so that all license and version
// info is in the same file

// MSVC6 has a max string limit of 2048 characters. Workaround:
//   char* msg = "a very" "long string" "can be constructed" "with multiple strings";

char * longChangesString = "Changes:\n\
 [Run \"bin_pablo -pd\" for details on switches]\n\
 [continued on http://titan.radonc.unc.edu/colab/Private/Pablo/changes.html]\n\
20050815/GT: version 1.0.14:\n\
 -for compatibility, world coordinate system is now read/written in M3D files:\n\
    model.world.{origin,bound,spacing}\n\
  in place of model.world.{origin,opposite}\n\
 -upgraded to glutBitmapCharacter()\n\
20050726/GT: version 1.0.14:\n\
 -modelRename disabled: fig names are used as filenames were in profile\n\
  extraction and should not contain extra info.\n\
 -image landmark: prefix of LM- on figure name is now optional\n\
20050722/GT: version 1.0.13:\n\
 -Model Window: landmark names allow numbers, lower-case letters and punctuation\n\
 -Model Window: multiple landmarks per atom are now all drawn\n\
 -Subdivision levels now labeled [0..3] instead of [1..4], to match values\n\
  a user types in a control file\n\
"
"20050720/GT: version 1.0.12:\n\
 -tiles:\n\
  -add \"file->save tile\" and \"generate tiles from (binary) image\"\n\
  -add \"-oit file.byu\" switch to generate & save a tile set\n\
  -Model Window: Center button will center tiles if no model is loaded\n\
  -^t now loads tiles from a file instead of generating tiles (now ^g)\n\
 -when model is loaded and -lv enabled, print model landmark positions in\n\
  <wcs> & <ics> coordinate systems\n\
 -modelRename switch: append outModel's filename to obj/fig names only when\n\
  missing from current names\n\
20050711/GT: version 1.0.11:\n\
 -multi-object (via separate figures) code is complete: new switches\n\
    involve interpenetration of neighboring obj's:\n\
  - figureID spec'd the figure being optimized (change to figureName soon)\n\
  - momMask specifies a bitStack mask used to constraint MOM to spec'd objs\n\
  - atom/figureInterpenetrationMatch is a new objective term for figure to avoid\n\
  - interpenetrationMask selects mask of objects not to be penetrated\n\
 -when loading a model, print landmarks positions in ucs and ics coordinates;\n\
  enabled by -lv switch; used to send landmark positions to other programs.\n\
 -add -pQ switch to guarantee non-interactive mode, eg, with -oit\n\
 -change warning from \"Not a valid registry file\" to \"INFO: optional registry\n\
  file missing\"\n\
 -BUG FIX: landmark as arbitrary spoke tips during optimization (only init stage\n\
  was tested in 1.0.9)\n\
"
"20050707/GT: version 1.0.10:\n\
 -KNOWN BUG: using an old image coordinate system; if exact model positions\n\
  are important to you (eg, within half a voxel), please wait for a version\n\
  in which this bug has been fixed.\n\
 -m3d files have 2 points (origin, opposite) which define the\n\
  world coordinate system of the first and last voxel centers of the\n\
  image that was loaded when the model was written. This allows us to\n\
  fit a model to a processed (eg, cropped or resampled) image and\n\
  position it properly on a corresponding unprocessed image;\n\
  these will change very soon!\n\
 -new switch \"-ph FILE\" prints an image's header plus summary info\n\
 -new switch -coordCubeImage enables -coordMaxExtent, off by default,\n\
  so that the default behavior is NOT to cube the image\n\
 -switches -pm & -ph are fast now; read image header only, not voxels\n\
 -bug fix: reads large binary images correctly from AFS\n\
 -BASE CODE for tubes (1xN atom circular cross-section figures)\n\
20050629/GT: version 1.0.9:\n\
 -prep for JL's multi-object changes:\n\
  -landmarks across multi fig's can match a single image landmark\n\
  -warn if no matching image lm is found instead of the reverse\n\
 -extend model landmarks from b-spoke to any spoke of atomIndex:\n\
  -new m3d file field:  model.figure[#].landmark[#].atomT\n\
   -values are -1=Y0 spoke (cyan), 0=bisector spoke (red), 1=Y1 spoke (magenta)\n\
   -warn if an internal atom's bisector spoke is used as a landmark\n\
  -no spoke indication for image landmarks; users could use landmark's name\n\
   to track each spokes purpose, eg, use landmark names of top-b & top-y0.\n\
20050617/GT: version 1.0.8md:\n\
"
" -bug fix: \"Min Dist\" report is again positive\n\
20050616/RS,GT: version 1.0.8:\n\
 -remove unused switch: figureAverageNeighborMatch\n\
 -image resampling was dilating the image; now maintains size.\n\
 -enable reading of signed distance map files -- users should\n\
  realize that they display incorrectly as unsigned images, but\n\
  the pixel intensities are interpreted properly for optimization.\n\
 -object-to-image: when search exceeds object or image\n\
  boundaries, use DDM distance as best guess.\n\
 -bug fix: -renameModel crashed when no name exists\n\
 -disable object-to-image cache & don't print: \"Cache: lookups\"\n\
20050610/RS,GT: version 1.0.7:\n\
 -Image match has been lightly tested - better tests coming...\n\
 -add switch DDMSearchThreshold: the voxel distance under\n\
  which object-to-image image match should stop searching the normal\n\
  and just trust the distance map.\n\
 -bug fix: transparent tiles are visible in front of image planes\n\
 -bug fix: figureSurfaceLevel was using value of atomSurfaceLevel\n\
 -bug fix: \"Min Dist\" report is now always positive\n\
"
"20050606/RS,GT: version 1.0.6:\n\
 -THIS VERSION'S IMAGE MATCH IS NOT CERTIFIED!\n\
  -TEMPORARILY disable reading/writing of signed distance\n\
   map files until we can verify maps are right\n\
  -new switch ignoreDistFile forces calc of new distance map\n\
 -switch renamings: image-to-object and object-to-image\n\
  distance behavior is identical across all stages:\n\
  -figureDMApproxEnable & atomDMApproxEnable -> DDMEnable\n\
  -figureAngleTolerance & atomAngleTolerance -> DDMAngleTolerance\n\
  -figureCenterDistance -> DDMCenter\n\
 -Read raw3 image containing extended header (EH) format - UNTESTED:\n\
  -never scale or shift pixel intensities during read\n\
  -never write EH when writing raw3 files\n\
 -append filename to Model and (first) Figure names when writing;\n\
  use add renameModel switch to leave names alone\n\
20050531/JL,REB: version 1.0.6:\n\
 -Conjugate Gradient bug fixes:\n\
  -GLIMIT was reset to 100.0 (from an experimental value of .05),\n\
   so that Brent Linear Search stepsize is more appropriate\n\
  -gradient direction is no longer normalized during conjugation,\n\
   so it searches in the proper direction now.\n\
"
"20050531/RS,GT: version 1.0.5-sab4-edgeDist:\n\
 -add figureCenterDistance to specify that the Danielsson-style Distance Map\n\
  use obsolete behavior: calc signed 8*distance to voxel centers instead\n\
  of unsigned 1*distance to voxel edges. Distances inside the object\n\
  are negative; outside are positive. Obsolete code will be deleted soon.\n\
  -add ignoreDistFile: overwrite DDM file, if it exists\n\
 -sab4: calc distance from each surface position along it's normal\n\
  by searching binary voxels for an edge:\n\
  -add figureDMApproxEnable, which enables the new search\n\
  -add figureAngleTolerance, which fine-tunes the trigger\n\
  -surface-based distances are NOT saved to the distance map image\n\
  -print of \"Average Dist\" metric excludes total distance\n\
 -default preference: show B vectors for crest atoms only\n\
 -default preference: do not smooth images\n\
 -hide the untested -id switch\n\
 -license 1.03: convert license calc's to double precision;\n\
  key includes site encoding now\n\
";

/* STASH - this part is only seen on web page to keep the rest short
  20050127/GT: version 1.0.4:\n\
 -add note of default value for figureMahalanobisMatch\n\
 -fix license bug, new license version is 1.02 \n\
20050126/GT: version 1.0.3:\n\
 -parse arguments even if license keyFile is absent, to print version\n\
 -multiple host licensing fully functional\n\
"
"20050125/GT: version 1.0.2:\n\
 -suppress confusing message: \"Preferences File: .\"\n\
 -change \"Click X\" to \"Close window\"\n\
20050119/GT: version 1.0.1:\n\
 -add tuning parameter \"voxelSpacing\"\n\
 -support multiple hosts in license\n\
 -bug fix: -pi: round image dimensions instead of truncate\n\
20050117/GT: version 1.0.0:\n\
 -OFFICIAL RELEASE upon which user guide is based\n\
 -requires license to run a control script or the GUI\n\
 -PGA; reads stats embedded in a M3D model file\n\
 -export tiles in world coordinates of image, if present\n\
 -all tunings may be set from a control script file\n\
 -draw landmark name in model window\n\
 -draw landmark b-spoke-end atoms as white balls\n\
 -Add Quad: use name \"LM-top\", to make a\n\
  LandMark named 'top' w/ tolerance R\n\
 -press ^t to tile unblurred binary image.\n\
*/

#else	/* BINARY */

char * longChangesString = "Changes:\n\
	None\n";

extern const char * release;	// Defined in pablo.cpp

#endif	/* BINARY */

// Allowable length of landmark names.  If increased, the size of the
// landmarkNameEditor and its default column width must be increased.
#define MAX_FLVW_STRING_LENGTH	40

#define COLOR_CUTOFF	1.45f

using namespace std;

const char * NONAME = "Unsaved Model";
#ifdef BINARY
const char * MAIN_WINDOW_TITLE = "Binary Pablo - ";
#else
const char * MAIN_WINDOW_TITLE = "Pablo - ";
#endif
const char * SAVE_ON_EXIT_MESSAGE = "Save changes to \"%s\"?";

extern const int MAX_OPTIMIZER_ITERATIONS;

//#define OUTPUT_DEFORMING_MODELS

const int CUT_PLANE_SIZE = 64;      // Dimension of the cut plane visualization images
const int IMAGE_PLANE_COUNTER_TOP = 20;	// Initial no. of delay cycles for slider motion

const char empty[1] = "";

const int DISPLAY_PATH_LINE_LEN = 50;

// Color slots not used in P3DUserInterface.fl or by FLTK (see FL/Enumeration.H)
const int FIGURAL_COLOR = 9;            // Unused color, used for setting figural colors
const int BACKGROUND_COLOR = 10;        // Unused color, used for setting background color
const int PARTIAL_SURFACE_COLOR = 11;   // Unused color, used for setting partial surface color
const int CONNECTORS_COLOR = 12;        // Unused color, used for setting mesh connectors color
const int LANDMARKS_COLOR = 13;         // Unused color, used for setting default landmarks color
// INDIGO is defined as 14 in P3DView.h and cannot be used here for a color code.
const int TILES_COLOR = 16;				// Unused color, used for setting default tiles color
float color[3];     // Temporary storage used below

static char editLandmarkNamesLabel_0[] = {
	"Figure      Landmark Name                                      Paired" };
static char editLandmarkNamesLabel_1[] = {
	"Figure      Landmark Name                                            Paired" };

static char editModelAttachmentModeLabel_0[] = "Protrusion  ";
static char editModelAttachmentModeLabel_1[] = "Indentation ";

static char editModelPolarityLabel_0[] = "Positive polarity";
static char editModelPolarityLabel_1[] = "Negative polarity";

static char cutPlaneModeLabel_0[] = "Show both cut plane views";
static char cutPlaneModeLabel_1[] = "Show cut planes only in model window";
static char cutPlaneModeLabel_2[] = "Show only cut plane windows";

extern const char * release;	// Defined in pablo.cpp

// Note: FLTK, shows the two braces below in the pull-down.  Not using the
// braces, however, makes fl_file_chooser() not find extensions after the
// first one in the initial set.
static char imageFileExtensionList[] =   // See also fltkUtilities.cxx in ConStruct
    "{*.raw3,*.pim,*.gipl,*.mha,*.mhd,*.hdr,plan_im}|*.raw3|*.pim,plan_im|*.gipl|*.mha|*.mhd|*.hdr";
const int NUM_VISIBILITY_FIG_BTNS = 7;


#ifdef BINARY

// Default pga file read from working directory
static const char * DEFAULT_PGA_FILENAME = "pablo.pga";

#endif


// collect filenames from a Fl browser
vector<string> getFilenamesFromBrowser(Fl_File_Browser *browser)
{
	vector<string> files;	
	for (int f=0; f < browser->size(); f++)
	{
		char *file = static_cast<char*> (browser->data(f));
		files.push_back(string(file));
	}
	return files;
}


const char * baseName(const char * filename)
{
	const char * slash = strrchr(filename, '/');
	const char * rslash = strrchr(filename, '\\');
	if (rslash > slash)
		slash = rslash;
	if (slash == NULL) {
		const char * colon = strrchr(filename, ':');
		if (colon == NULL)
			return filename;
		else
			return colon + 1;
	}
	return slash + 1;
}


// These global functions are call-backs used in P3DControl to
// pass information to the Script Display Code below.  They are
// set in P3DUserInterfaceCallback::startInteractive(), which
// is called near the end of main().
void updateObjectInGui(int iteration)
{
	extern P3DUserInterface * ui;

	ui->callback->showInteractiveResults(iteration);
	Fl::check();
}

void updateOptimizerStageInGui(OptimizerStage stage)
{
	extern P3DUserInterface * ui;

	ui->callback->updateOptimizerWindow(stage);
	Fl::check();
}


P3DUserInterfaceCallback::P3DUserInterfaceCallback()
{
    bool left_hand_mouse;

	halting = false;
	usingModalWindow = 0;
    optimizerStopRequested = false;
    optimizerIsRunning = false;
    optimizerWindowInUse = false;
	interactiveOpt = false;
    resettingOptimizer = false;

    cut_plane_intensity = 0.0f;
	for (int i = 0; i < 6; i++) {
		cutPlanesShown[i] = false;
		cutPlaneModes[i] = 0;
	}
	cutPlaneModePreference = 0;
	cutPlaneModeLabels[0] = "2";
	cutPlaneModeLabels[1] = "C";
	cutPlaneModeLabels[2] = "M";
    selectedPartialFigures = NULL;
    partialSurfaceSliderChanged = false;
    markedFigurePtr = NULL; // Only used in constraints code
    imagePlaneSliderMotion = 0;
    displayModelPath = NULL;
    displayImagePath = NULL;
    displayTilesetPath = NULL;
    recordWindowPositions = false;
	recordWindowsOpen = false;
    partialSurfaceColor[0] = (float) globalControl->readDouble(PartialSurfaceColor_R);
    partialSurfaceColor[1] = (float) globalControl->readDouble(PartialSurfaceColor_G);
    partialSurfaceColor[2] = (float) globalControl->readDouble(PartialSurfaceColor_B);
	modelNameEditor = NULL;
	landmarkNameEditor = NULL;
	mainWindowLabel = NULL;
	landmarksEditMode = -1;
	displayModelCoords = globalControl->read(DisplayModelUnits, false);
	displayIntensityScale = globalControl->read(DisplayRangeUnits, false);

	brightTiles = 0;

	centerThing = CEN_MODEL|CEN_ALL;

    // Set up the interactive output for attaching subfigures
#ifdef _WIN32
    left_hand_mouse = (0 != GetSystemMetrics(SM_SWAPBUTTON));
#else
    left_hand_mouse = globalControl->read(LeftHandedMouse, false);
#endif
	set_mouse_button_order(left_hand_mouse);

#ifdef BINARY
	distanceVectorList = NULL;

    pgaFilename = new char[strlen(DEFAULT_PGA_FILENAME) + 1];
    strcpy(pgaFilename, DEFAULT_PGA_FILENAME);
#else
	pgaFilename = NULL;
#endif
	profileFilename = NULL;
	previewMode = false;
}

P3DUserInterfaceCallback::~P3DUserInterfaceCallback()
{
#ifdef DEBUG
	cout << "P3DUserInterfaceCallback::~P3DUserInterfaceCallback()" << endl;
#endif
	if (profileFilename != NULL)
		delete [] profileFilename;

    if (selectedPartialFigures != NULL)
        delete selectedPartialFigures;

    if (displayModelPath != NULL)
        delete [] displayModelPath;
    if (displayImagePath != NULL)
        delete [] displayImagePath;
    if (displayTilesetPath != NULL)
        delete [] displayTilesetPath;

    delete [] mainWindowLabel;
    // markedFigurePtr should never be deleted

    if (pgaFilename)
        delete [] pgaFilename;

#ifdef BINARY
	if (distanceVectorList != NULL) {
		distanceVectorList->plist.clear();
		delete distanceVectorList;
	}
#endif
}

void P3DUserInterfaceCallback::set_mouse_button_order(bool left_hand_mouse) {
	mouse_left_handed = left_hand_mouse;

    strcpy(ATTACH_STEP2, "Mark (ctrl-");
    if (left_hand_mouse)
        strcat(ATTACH_STEP2, "left");
    else
        strcat(ATTACH_STEP2, "right");
    strcat(ATTACH_STEP2, " click) an atom on\nthe parent figure.");

    strcpy(ATTACH_STEP3, "Mark (ctrl-");
    if (left_hand_mouse)
        strcat(ATTACH_STEP3, "left");
    else
        strcat(ATTACH_STEP3, "right");
    strcat(ATTACH_STEP3, " click) an atom on\nthe child figure.");

    strcpy(ATTACH_STEP4, "Select (ctrl-");
    if (left_hand_mouse)
        strcat(ATTACH_STEP4, "right");
    else
        strcat(ATTACH_STEP4, "left");
    strcat(ATTACH_STEP4,
        " click) a row or column of\nhinge atoms on the child figure.");

    strcpy(ATOM_EDITOR_STEP1, "Mark (ctrl-");
    if (left_hand_mouse)
        strcat(ATOM_EDITOR_STEP1, "left");
    else
        strcat(ATOM_EDITOR_STEP1, "right");
    strcat(ATOM_EDITOR_STEP1, " click) any atom");

    strcpy(ELONGATION_STEP1, "Mark (ctrl-");
    if (left_hand_mouse)
        strcat(ELONGATION_STEP1, "left");
    else
        strcat(ELONGATION_STEP1, "right");
    strcat(ELONGATION_STEP1, " click) an atom of the figure");
}

void P3DUserInterfaceCallback::closeWindowsAndExit() {
    if (userInterface == NULL)
        return;

	halting = true;

    if (! askToSaveModel())
		return;    // Canceled

    if (recordWindowPositions)
        saveWindowPositions();
    if (recordWindowsOpen)
        saveWindowsOpen();

#ifdef OPTIMIZATION_VISUALIZER
	if (ovUI)
		((OptVisualizerUI *) ovUI)->hide();
#endif

    if (userInterface->displayControlWindow)
        userInterface->displayControlWindow->hide();

    if (userInterface->addQuadFigureDlg)
        userInterface->addQuadFigureDlg->hide();

    if (userInterface->editModelWindow)
        userInterface->editModelWindow->hide();

    if (userInterface->visibilityControlWindow)
        userInterface->visibilityControlWindow->hide();

    if (userInterface->constraintsWindow)
        userInterface->constraintsWindow->hide();

    if (userInterface->editLandmarksWindow)
        userInterface->editLandmarksWindow->hide();

    if (userInterface->primitiveEditorWindow)
        userInterface->primitiveEditorWindow->hide();

    if (userInterface->bNCutPlaneWindow)
        userInterface->bNCutPlaneWindow->hide();

    if (userInterface->bBPerpCutPlaneWindow)
        userInterface->bBPerpCutPlaneWindow->hide();

    if (userInterface->bPerpNCutPlaneWindow)
        userInterface->bPerpNCutPlaneWindow->hide();

    if (userInterface->bPerpY0CutPlaneWindow)
        userInterface->bPerpY0CutPlaneWindow->hide();

    if (userInterface->bPerpY1CutPlaneWindow)
        userInterface->bPerpY1CutPlaneWindow->hide();

    if (userInterface->InvoluteCutPlaneWindow)
        userInterface->InvoluteCutPlaneWindow->hide();

    if (userInterface->cutPlanesControlWindow)
        userInterface->cutPlanesControlWindow->hide();

    if (userInterface->optimizerControlDialog)
        userInterface->optimizerControlDialog->hide();

    if (userInterface->optimizerSettingsDialog)
        userInterface->optimizerSettingsDialog->hide();

    if (userInterface->attachSubfigureDialog)
        userInterface->attachSubfigureDialog->hide();

    if (userInterface->elongationDialog)
        userInterface->elongationDialog->hide();

    if(userInterface->pgaDialog)
        userInterface->pgaDialog->hide();

    if(userInterface->cpnsDialog)
        userInterface->cpnsDialog->hide();

#ifdef PCA_VISUALIZER
    if(userInterface->pcaDialog)
        userInterface->pcaDialog->hide();
#endif

    if(userInterface->transformationRecordingTestDialog)
        userInterface->transformationRecordingTestDialog->hide();

// The test windows need to go here, with any other windows that have been omitted from this list.
#ifndef PRODUCTION_VERSION
    if (userInterface->interpolatedPrimitiveWindow)
        userInterface->interpolatedPrimitiveWindow->hide();

    if(userInterface->modelSlideShowDialog)
        userInterface->modelSlideShowDialog->hide();

    if(userInterface->regularizerWindow)
        userInterface->regularizerWindow->hide();

    if(userInterface->matchSurfacesDialog)
        userInterface->matchSurfacesDialog->hide();
#endif	/* PRODUCTION_VERSION */

    if (userInterface->prefsImageFilesWindow)
        userInterface->prefsImageFilesWindow->hide();

    if (userInterface->preferencesEditorWindow)
        userInterface->preferencesEditorWindow->hide();

    userInterface->mainWindow->hide();
    userInterface->modelWindow->hide();

#ifdef BINARY
    if (userInterface->aboutBinaryPabloWindow)
        userInterface->aboutBinaryPabloWindow->hide();
#else
    if (userInterface->aboutPabloWindow)
        userInterface->aboutPabloWindow->hide();
#endif

    // Destruction of all windows results in a return from the
    // call to Fl::run() in pablo.cpp.  This is necessary for
    // changes in the control file to be saved.
}

P3DControl * P3DUserInterfaceCallback::getControl() const
{
    if(userInterface == NULL)
        return NULL;

    return userInterface->control;
}

// ----------------------------- Start of Script Display Code --------------------------------

void P3DUserInterfaceCallback::startInteractive(int matchMethod)
{
	getControl()->setGuiCallbacks(updateObjectInGui, updateOptimizerStageInGui);
	interactiveOpt = true;

	// Translate matchMethod numbers into button numbers.
	// See the OptimizerMatchMethod entry in gpTuning.cpp.
	if (matchMethod > 5)
		matchMethod -= 4;
	userInterface->setOptimizerMatchBtns(matchMethod);
	userInterface->nextOptimizationBtn->hide();
	userInterface->previewOptimizationBtn->hide();
	userInterface->continuousOptimizationBtn->show();
	userInterface->startNextInteractiveStageBtn->show();

	Fl::check();
	// Now the Optimizer Window's Setup page will indicate the match
	// method to be used, while P3DControl is loading the files to
	// be used in the optimization.

	continueInteractiveOpt = false;
	currentInteractiveStage = OPTIMIZER_SETUP_STAGE;
}

void P3DUserInterfaceCallback::stopInteractive()
{
	getControl()->setGuiCallbacks(NULL, NULL);
	interactiveOpt = false;
	updateOptimizerWindow(OPTIMIZER_END_STAGE);
	userInterface->startNextInteractiveStageBtn->hide();
	userInterface->continuousOptimizationBtn->hide();
	userInterface->previewOptimizationBtn->show();
	userInterface->nextOptimizationBtn->show();
}

void P3DUserInterfaceCallback::updateOptimizerWindow(OptimizerStage nextStageId)
{
#ifndef BINARY

	OptimizerStage stage;

    nextInteractiveStage = nextStageId;

	if (! interactiveOpt)
		return;

	// Update the Optimizer Control Window, even if it is hidden, since it
	// may be opened at any time.
	stage = nextOptimizationStage(true);
	if (stage != nextStageId)
		cout << "Warning: Mismatch of interactive and script stages"
			<< endl;	// This should never be printed
	enableOptimzerButtons(currentInteractiveStage);
	disableOptimzerButtons(stage);
	currentInteractiveStage = stage;

	if (nextStageId == OPTIMIZER_END_STAGE)
		return;		// stopInteractive() was called
	if (nextStageId == OPTIMIZER_SUBFIGURE_STAGE)
		return;		// TEMPORARY: Because this stage is not presently working
	if (nextStageId == OPTIMIZER_BOUNDARY_DISPLACEMENT_STAGE)
		return;		// TEMPORARY: Because this stage is not presently in use

	if (! userInterface->optimizerControlDialog->shown())
		return;		// Without the window, there are no buttons to be pressed

	Fl::check();	// This makes the interface somewhat more responsive
	if (continueInteractiveOpt)
		return; 	// Continuous operation requested
	doNextInteractiveOpt = false;
	userInterface->startNextInteractiveStageBtn->activate();
	while (! doNextInteractiveOpt)
		Fl::check();	//  Could use a timer here

#endif	/* BINARY */
}

void P3DUserInterfaceCallback::startNextInteractiveStage()
{
	doNextInteractiveOpt = true;
	userInterface->startNextInteractiveStageBtn->deactivate();
}

void P3DUserInterfaceCallback::continuousOptimization()
{
	continueInteractiveOpt = true;
	doNextInteractiveOpt = true;
	userInterface->startNextInteractiveStageBtn->deactivate();
	userInterface->continuousOptimizationBtn->deactivate();
}

void P3DUserInterfaceCallback::enableOptimzerButtons(OptimizerStage stageId)
{
    switch(stageId)
    {
        case OPTIMIZER_OBJECT_STAGE:
			userInterface->objectStartBtn->activate();
			userInterface->objectStopBtn->activate();
			userInterface->objectResetBtn->activate();
			userInterface->objectPenaltyWeightSlider->activate();
            break;

        case OPTIMIZER_ENSEMBLE_STAGE:
			userInterface->modelStartBtn->activate();
			userInterface->modelStopBtn->activate();
			userInterface->modelResetBtn->activate();
			userInterface->modelPenaltyWeightSlider->activate();
            break;

        case OPTIMIZER_FIGURE_STAGE:
			userInterface->figureStageStartBtn->activate();
			userInterface->figureStageStopBtn->activate();
			userInterface->figureStagePenaltyWeightSlider->activate();
			userInterface->figureStageConstraintsPenaltyWeightSlider->activate();
            break;

        case OPTIMIZER_MAIN_FIGURE_STAGE:
			userInterface->mainFigureStartBtn->activate();
			userInterface->mainFigureStopBtn->activate();
			userInterface->mainFigurePenaltyWeightSlider->activate();
			userInterface->mainFigureConstraintsPenaltyWeightSlider->activate();
            break;

        case OPTIMIZER_SUBFIGURE_STAGE:
			userInterface->subfigureStartBtn->activate();
			userInterface->subfigureStopBtn->activate();
			userInterface->subfigurePenaltyWeightSlider->activate();
			userInterface->subfigureConstraintsPenaltyWeightSlider->activate();
            break;

        case OPTIMIZER_PRIMITIVE_STAGE:
			userInterface->primitiveStartBtn->activate();
			userInterface->primitiveStopBtn->activate();
			userInterface->primitivePenaltyWeightSlider->activate();
			userInterface->primitiveConstraintsPenaltyWeightSlider->activate();
			userInterface->primitiveNeighborPenaltyWeightSlider->activate();
            break;

        case OPTIMIZER_BOUNDARY_DISPLACEMENT_STAGE:
			userInterface->voxelStageStartBtn->activate();
			userInterface->voxelStageStopBtn->activate();
			userInterface->voxelStagePenaltyWeightSlider->activate();
			userInterface->voxelStageConstraintsPenaltyWeightSlider->activate();
            break;

        case OPTIMIZER_END_STAGE:
        case OPTIMIZER_SETUP_STAGE:
        case NO_OPTIMIZER_STAGE:
        default:
            break;
    };
}

void P3DUserInterfaceCallback::disableOptimzerButtons(OptimizerStage stageId)
{
    switch(stageId)
    {
        case OPTIMIZER_OBJECT_STAGE:
			userInterface->objectStartBtn->deactivate();
			userInterface->objectStopBtn->deactivate();
			userInterface->objectResetBtn->deactivate();
			userInterface->objectPenaltyWeightSlider->deactivate();
            break;

        case OPTIMIZER_ENSEMBLE_STAGE:
			userInterface->modelStartBtn->deactivate();
			userInterface->modelStopBtn->deactivate();
			userInterface->modelResetBtn->deactivate();
			userInterface->modelPenaltyWeightSlider->deactivate();
            break;

        case OPTIMIZER_FIGURE_STAGE:
			userInterface->figureStageStartBtn->deactivate();
			userInterface->figureStageStopBtn->deactivate();
			userInterface->figureStagePenaltyWeightSlider->deactivate();
			userInterface->figureStageConstraintsPenaltyWeightSlider->deactivate();
            break;

        case OPTIMIZER_MAIN_FIGURE_STAGE:
			userInterface->mainFigureStartBtn->deactivate();
			userInterface->mainFigureStopBtn->deactivate();
			userInterface->mainFigurePenaltyWeightSlider->deactivate();
			userInterface->mainFigureConstraintsPenaltyWeightSlider->deactivate();
            break;

        case OPTIMIZER_SUBFIGURE_STAGE:
			userInterface->subfigureStartBtn->deactivate();
			userInterface->subfigureStopBtn->deactivate();
			userInterface->subfigurePenaltyWeightSlider->deactivate();
			userInterface->subfigureConstraintsPenaltyWeightSlider->deactivate();
            break;

        case OPTIMIZER_PRIMITIVE_STAGE:
			userInterface->primitiveStartBtn->deactivate();
			userInterface->primitiveStopBtn->deactivate();
			userInterface->primitivePenaltyWeightSlider->deactivate();
			userInterface->primitiveConstraintsPenaltyWeightSlider->deactivate();
			userInterface->primitiveNeighborPenaltyWeightSlider->deactivate();
            break;

        case OPTIMIZER_BOUNDARY_DISPLACEMENT_STAGE:
			userInterface->voxelStageStartBtn->deactivate();
			userInterface->voxelStageStopBtn->deactivate();
			userInterface->voxelStagePenaltyWeightSlider->deactivate();
			userInterface->voxelStageConstraintsPenaltyWeightSlider->deactivate();
            break;

        case OPTIMIZER_END_STAGE:
        case OPTIMIZER_SETUP_STAGE:
        case NO_OPTIMIZER_STAGE:
        default:
            break;
    };
}

void P3DUserInterfaceCallback::showInteractiveResults(int iteration)
{
    P3DControl * control;

    control = userInterface->control;
	if (control == NULL)
		return;

    switch(currentInteractiveStage)
    {
        case OPTIMIZER_OBJECT_STAGE:
			userInterface->objectIterationOutput->value(iteration);
			userInterface->objectMatchValueOutput->value(
				control->getLastRegistrationMatchValue());
			userInterface->objectPenaltyValueOutput->value(
				control->getLastRegistrationPenaltyValue());
            break;

        case OPTIMIZER_ENSEMBLE_STAGE:
			userInterface->modelIterationOutput->value(iteration);
			userInterface->modelMatchValueOutput->value(
				control->getLastRegistrationMatchValue());
			userInterface->modelPenaltyValueOutput->value(
				control->getLastRegistrationPenaltyValue());
            break;

        case OPTIMIZER_FIGURE_STAGE:
			userInterface->figureStageIterationOutput->value(iteration);
 			userInterface->figureStagePenaltyValueOutput->value(
				control->getLastMainFigureGeometricMatchValue());
			userInterface->figureStageObjectiveFunctionValueOutput->value(
				control->getLastMainFigureObjectiveFunctionValue());
			userInterface->figureStageImageMatchValueOutput->value(
				control->getLastMainFigureObjectiveFunctionValue()
				- control->getLastMainFigureGeometricMatchValue());
            break;

       case OPTIMIZER_MAIN_FIGURE_STAGE:
			userInterface->mainFigureIterationOutput->value(iteration);
			userInterface->mainFigurePenaltyValueOutput->value(
				control->getLastMainFigureGeometricMatchValue());
			userInterface->mainFigureObjectiveFunctionValueOutput->value(
				control->getLastMainFigureObjectiveFunctionValue());
			userInterface->mainFigureImageMatchValueOutput->value(
				control->getLastMainFigureObjectiveFunctionValue()
				- control->getLastMainFigureGeometricMatchValue());
            break;

        case OPTIMIZER_SUBFIGURE_STAGE:
#ifndef BINARY
			userInterface->subfigureIterationOutput->value(iteration);
			userInterface->subfigureMatchValueOutput->value(
				control->getLastSubfigureMatchValue());
#endif
            break;

        case OPTIMIZER_PRIMITIVE_STAGE:
			userInterface->primitiveIterationOutput->value(iteration);
			userInterface->primitiveObjectiveFunctionValueOutput->value(
				control->getLastDeformationObjectiveFunctionValue());
			userInterface->primativePenaltyValueOutput->value(
				control->getLastDeformationGeometricMatchValue());
			userInterface->primitiveImageMatchValueOutput->value(
				control->getLastDeformationObjectiveFunctionValue()
				- control->getLastDeformationGeometricMatchValue());
            break;

        case OPTIMIZER_BOUNDARY_DISPLACEMENT_STAGE:
			userInterface->voxelStageIterationOutput->value(iteration);
			userInterface->voxelStageMatchValueOutput->value(
				control->getLastVoxelMatchValue());
			userInterface->voxelStageObjectiveFunctionValueOutput->
				value(-control->getLastVoxelObjectiveFunctionValue());
			userInterface->voxelStagePenaltyValueOutput->value(
				control->getLastVoxelPenaltyValue());
            break;

        case OPTIMIZER_END_STAGE:
        case OPTIMIZER_SETUP_STAGE:
        case NO_OPTIMIZER_STAGE:
        default:
            break;
    };
	update();
}
// ----------------------------- End of Script Display Code --------------------------------

void P3DUserInterfaceCallback::saveWindowPositions() {
    int x, y;
	int oldX, oldY;
	int dX, dY;
	bool moved;

    if (userInterface == NULL)
        return;

    if (userInterface->displayControlWindow) {
        x = userInterface->displayControlWindow->x();
        y = userInterface->displayControlWindow->y();
        globalControl->write(DisplayControlWindow_X, x);
        globalControl->write(DisplayControlWindow_Y, y);
    }
    if (userInterface->addQuadFigureDlg) {
        x = userInterface->addQuadFigureDlg->x();
        y = userInterface->addQuadFigureDlg->y();
        globalControl->write(AddQuadFigureWindow_X, x);
        globalControl->write(AddQuadFigureWindow_Y, y);
    }
    if (userInterface->visibilityControlWindow) {
        x = userInterface->visibilityControlWindow->x();
        y = userInterface->visibilityControlWindow->y();
        globalControl->write(VisibilityControlWindow_X, x);
        globalControl->write(VisibilityControlWindow_Y, y);
    }
    if (userInterface->constraintsWindow) {
        x = userInterface->constraintsWindow->x();
        y = userInterface->constraintsWindow->y();
        globalControl->write(ConstraintsWindow_X, x);
        globalControl->write(ConstraintsWindow_Y, y);
    }
    if (userInterface->primitiveEditorWindow) {
        x = userInterface->primitiveEditorWindow->x();
        y = userInterface->primitiveEditorWindow->y();
        globalControl->write(AtomEditorWindow_X, x);
        globalControl->write(AtomEditorWindow_Y, y);
    }
    if (userInterface->bNCutPlaneWindow) {
        x = userInterface->bNCutPlaneWindow->x();
        y = userInterface->bNCutPlaneWindow->y();
        globalControl->write(AtomPlaneWindow_X, x);
        globalControl->write(AtomPlaneWindow_Y, y);
    }
    if (userInterface->bBPerpCutPlaneWindow) {
        x = userInterface->bBPerpCutPlaneWindow->x();
        y = userInterface->bBPerpCutPlaneWindow->y();
        globalControl->write(CrestPlaneWindow_X, x);
        globalControl->write(CrestPlaneWindow_Y, y);
    }
    if (userInterface->bPerpNCutPlaneWindow) {
        x = userInterface->bPerpNCutPlaneWindow->x();
        y = userInterface->bPerpNCutPlaneWindow->y();
        globalControl->write(BPerpNPlaneWindow_X, x);
        globalControl->write(BPerpNPlaneWindow_Y, y);
    }
    if (userInterface->bPerpY0CutPlaneWindow) {
        x = userInterface->bPerpY0CutPlaneWindow->x();
        y = userInterface->bPerpY0CutPlaneWindow->y();
        globalControl->write(StarboardPlaneWindow_X, x);
        globalControl->write(StarboardPlaneWindow_Y, y);
    }
    if (userInterface->bPerpY1CutPlaneWindow) {
        x = userInterface->bPerpY1CutPlaneWindow->x();
        y = userInterface->bPerpY1CutPlaneWindow->y();
        globalControl->write(PortPlaneWindow_X, x);
        globalControl->write(PortPlaneWindow_Y, y);
    }
    if (userInterface->InvoluteCutPlaneWindow) {
        x = userInterface->InvoluteCutPlaneWindow->x();
        y = userInterface->InvoluteCutPlaneWindow->y();
        globalControl->write(InvolutesPlaneWindow_X, x);
        globalControl->write(InvolutesPlaneWindow_Y, y);
    }
    if (userInterface->cutPlanesControlWindow) {
        x = userInterface->cutPlanesControlWindow->x();
        y = userInterface->cutPlanesControlWindow->y();
        globalControl->write(CutPlanesControlWindow_X, x);
        globalControl->write(CutPlanesControlWindow_Y, y);
    }
    if (userInterface->optimizerControlDialog) {
        x = userInterface->optimizerControlDialog->x();
        y = userInterface->optimizerControlDialog->y();
        globalControl->write(OptimizerControlWindow_X, x);
        globalControl->write(OptimizerControlWindow_Y, y);
    }
    if (userInterface->optimizerSettingsDialog) {
        x = userInterface->optimizerSettingsDialog->x();
        y = userInterface->optimizerSettingsDialog->y();
        globalControl->write(PenaltyWeightsWindow_X, x);
        globalControl->write(PenaltyWeightsWindow_Y, y);
    }
    if (userInterface->attachSubfigureDialog) {
        x = userInterface->attachSubfigureDialog->x();
        y = userInterface->attachSubfigureDialog->y();
        globalControl->write(AttachSubfigureWindow_X, x);
        globalControl->write(AttachSubfigureWindow_Y, y);
    }
    if (userInterface->elongationDialog) {
        x = userInterface->elongationDialog->x();
        y = userInterface->elongationDialog->y();
        globalControl->write(ElongationWindow_X, x);
        globalControl->write(ElongationWindow_Y, y);
    }
    if (userInterface->pgaDialog) {
        x = userInterface->pgaDialog->x();
        y = userInterface->pgaDialog->y();
        globalControl->write(PGADeformWindow_X, x);
        globalControl->write(PGADeformWindow_Y, y);
    }

    if (userInterface->cpnsDialog) {
        x = userInterface->cpnsDialog->x();
        y = userInterface->cpnsDialog->y();
        globalControl->write(CPNSDeformWindow_X, x);
        globalControl->write(CPNSDeformWindow_Y, y);
    }

#ifdef PCA_VISUALIZER
    if (userInterface->pcaDialog) {
        x = userInterface->pcaDialog->x();
        y = userInterface->pcaDialog->y();
        globalControl->write(PCADeformWindow_X, x);
        globalControl->write(PCADeformWindow_Y, y);
    }
#endif

#ifdef OPTIMIZATION_VISUALIZER
	if (ovUI) {
		if (((OptVisualizerUI *) ovUI)->optVisualizerUIWindow) {
			x = ((OptVisualizerUI *) ovUI)->optVisualizerUIWindow->x();
			y = ((OptVisualizerUI *) ovUI)->optVisualizerUIWindow->y();
			globalControl->write(OptVisualizerWindow_X, x);
			globalControl->write(OptVisualizerWindow_Y, y);
		}
	}
#endif

// The test windows need to go here, with any other windows that have been omitted from this list.
#ifndef PRODUCTION_VERSION
    if (userInterface->interpolatedPrimitiveWindow) {
        x = userInterface->interpolatedPrimitiveWindow->x();
        y = userInterface->interpolatedPrimitiveWindow->y();
        globalControl->write(InterpolatedPrimitiveWindow_X, x);
        globalControl->write(InterpolatedPrimitiveWindow_Y, y);
    }
    if (userInterface->regularizerWindow) {
        x = userInterface->regularizerWindow->x();
        y = userInterface->regularizerWindow->y();
        globalControl->write(RegularizerWindow_X, x);
        globalControl->write(RegularizerWindow_Y, y);
    }
    if (userInterface->matchSurfacesDialog) {
        x = userInterface->matchSurfacesDialog->x();
        y = userInterface->matchSurfacesDialog->y();
        globalControl->write(MatchSurfaces_X, x);
        globalControl->write(MatchSurfaces_Y, y);
    }
#endif	/* PRODUCTION_VERSION */

	moved = false;
    if (userInterface->preferencesEditorWindow) {
		oldX = globalControl->readInt(PreferencesEditorWindow_X);
		oldY = globalControl->readInt(PreferencesEditorWindow_Y);
        x = userInterface->preferencesEditorWindow->x();
        y = userInterface->preferencesEditorWindow->y();
        globalControl->write(PreferencesEditorWindow_X, x);
        globalControl->write(PreferencesEditorWindow_Y, y);
		if (x != oldX || y != oldY) {
			dX = x - oldX;
			dY = y - oldY;
			moved = true;
		}
    }
    if (userInterface->prefsImageFilesWindow) {
		oldX = globalControl->readInt(ImagePrefsSubEditorWindow_X);
		oldY = globalControl->readInt(ImagePrefsSubEditorWindow_Y);
        x = userInterface->prefsImageFilesWindow->x();
        y = userInterface->prefsImageFilesWindow->y();
		if (moved && x == oldX && y == oldY) {
			// User moved Preferences Editor Window, but not the Image
			// Preferences Subwindow; drag the later along to keep the
			// same relationship between them.
			x = oldX + dX;
			y = oldY + dY;
		}
        globalControl->write(ImagePrefsSubEditorWindow_X, x);
        globalControl->write(ImagePrefsSubEditorWindow_Y, y);
    }
    if (userInterface->mainWindow) {
        x = userInterface->mainWindow->x();
        y = userInterface->mainWindow->y();
        globalControl->write(MainWindow_X, x);
        globalControl->write(MainWindow_Y, y);
    }
    if (userInterface->modelWindow) {
        x = userInterface->modelWindow->x();
        y = userInterface->modelWindow->y();
        globalControl->write(ModelWindow_X, x);
        globalControl->write(ModelWindow_Y, y);
    }
}

void P3DUserInterfaceCallback::saveWindowsOpen() {
    if (userInterface == NULL)
        return;

    if (userInterface->displayControlWindow->shown())
        globalControl->write(DisplayControlWindow_Open, true);
	else
        globalControl->write(DisplayControlWindow_Open, false);

    if (userInterface->addQuadFigureDlg->shown())
        globalControl->write(AddQuadFigureWindow_Open, true);
	else
        globalControl->write(AddQuadFigureWindow_Open, false);

    if (userInterface->visibilityControlWindow->shown())
        globalControl->write(VisibilityControlWindow_Open, true);
	else
        globalControl->write(VisibilityControlWindow_Open, false);

    if (userInterface->constraintsWindow->shown())
        globalControl->write(ConstraintsWindow_Open, true);
	else
        globalControl->write(ConstraintsWindow_Open, false);

    if (userInterface->primitiveEditorWindow->shown())
        globalControl->write(AtomEditorWindow_Open, true);
	else
        globalControl->write(AtomEditorWindow_Open, false);

    if (userInterface->bNCutPlaneWindow->shown())
        globalControl->write(AtomPlaneWindow_Open, true);
	else
        globalControl->write(AtomPlaneWindow_Open, false);

    if (userInterface->bBPerpCutPlaneWindow->shown())
        globalControl->write(CrestPlaneWindow_Open, true);
	else
        globalControl->write(CrestPlaneWindow_Open, false);

    if (userInterface->bPerpNCutPlaneWindow->shown())
        globalControl->write(BPerpNPlaneWindow_Open, true);
	else
        globalControl->write(BPerpNPlaneWindow_Open, false);

    if (userInterface->bPerpY0CutPlaneWindow->shown())
        globalControl->write(StarboardPlaneWindow_Open, true);
	else
        globalControl->write(StarboardPlaneWindow_Open, false);

    if (userInterface->bPerpY1CutPlaneWindow->shown())
        globalControl->write(PortPlaneWindow_Open, true);
	else
        globalControl->write(PortPlaneWindow_Open, false);

    if (userInterface->InvoluteCutPlaneWindow->shown())
        globalControl->write(InvolutesPlaneWindow_Open, true);
	else
        globalControl->write(InvolutesPlaneWindow_Open, false);

    if (userInterface->cutPlanesControlWindow->shown())
        globalControl->write(CutPlanesControlWindow_Open, true);
	else
        globalControl->write(CutPlanesControlWindow_Open, false);

    if (userInterface->optimizerControlDialog->shown())
        globalControl->write(OptimizerControlWindow_Open, true);
	else
        globalControl->write(OptimizerControlWindow_Open, false);

    if (userInterface->optimizerSettingsDialog->shown())
        globalControl->write(PenaltyWeightsWindow_Open, true);
	else
        globalControl->write(PenaltyWeightsWindow_Open, false);

    if (userInterface->attachSubfigureDialog->shown())
        globalControl->write(AttachSubfigureWindow_Open, true);
	else
        globalControl->write(AttachSubfigureWindow_Open, false);

    if (userInterface->elongationDialog->shown())
        globalControl->write(ElongationWindow_Open, true);
	else
        globalControl->write(ElongationWindow_Open, false);

    if (userInterface->pgaDialog->shown())
        globalControl->write(PGADeformWindow_Open, true);
	else
        globalControl->write(PGADeformWindow_Open, false);

    if (userInterface->cpnsDialog->shown())
        globalControl->write(CPNSDeformWindow_Open, true);
	else
        globalControl->write(CPNSDeformWindow_Open, false);

// The test windows need to go here, with any other windows that have been omitted from this list.
#ifndef PRODUCTION_VERSION
    if (userInterface->interpolatedPrimitiveWindow->shown())
        globalControl->write(InterpolatedPrimitiveWindow_Open, true);
	else
        globalControl->write(InterpolatedPrimitiveWindow_Open, false);

    if (userInterface->regularizerWindow->shown())
        globalControl->write(RegularizerWindow_Open, true);
	else
        globalControl->write(RegularizerWindow_Open, false);

    if (userInterface->matchSurfacesDialog->shown())
        globalControl->write(MatchSurfaces_Open, true);
	else
        globalControl->write(MatchSurfaces_Open, false);
#endif	/* PRODUCTION_VERSION */

    if (userInterface->preferencesEditorWindow->shown())
        globalControl->write(PreferencesEditorWindow_Open, true);
	else
        globalControl->write(PreferencesEditorWindow_Open, false);

    if (userInterface->mainWindow->shown())
        globalControl->write(MainWindow_Open, true);
	else
        globalControl->write(MainWindow_Open, false);

    if (userInterface->modelWindow->shown())
        globalControl->write(ModelWindow_Open, true);
	else
        globalControl->write(ModelWindow_Open, false);

#ifdef OPTIMIZATION_VISUALIZER
	if (ovUI) {
		if (((OptVisualizerUI *) ovUI)->optVisualizerUIWindow) {
			if (((OptVisualizerUI *) ovUI)->optVisualizerUIWindow->shown())
				globalControl->write(OptVisualizerWindow_Open, true);
			else
				globalControl->write(OptVisualizerWindow_Open, false);
		}
	}
#endif
}

void P3DUserInterfaceCallback::init()
{
    P3DView * modelView;

    if(userInterface == NULL)
        return;

	int undoLength = globalControl->readInt(UndoListLength);
    userInterface->control = new P3DControl(undoLength);
    if (! userInterface->control->verifyApplication(Pablo)) {
        cout << "Error: Bad configuration; you need to recompile Pablo using";
		cout << " the correct\n    configuration." << endl;
        exit(1);    // OK to exit, since nothing has been done yet.
    }

    modelView = userInterface->modelView;
    if(modelView != NULL)
        modelView->setControl(userInterface->control);

    setMainWindowLabel();
    initializeInterfaceVisibility();

#ifndef BINARY
    initializeOptimizerSettings();
#endif

    bool flag = globalControl->read(PartialSurfaceRendering, false);
    userInterface->control->togglePartialRendering(flag);
    userInterface->setPrefsPartialSurfaceRendering((int) flag);
    userInterface->setPartialSurfaceRendering((int) flag);

	// Read and set image slice settings
    if (globalControl->readBool(AxialSliceDefault)) {
		Fl_Check_Button * checkBox;
		checkBox = userInterface->zImagePlaneCheckBox;
		checkBox->value(1);
		viewZImagePlane(1);
    }
    if (globalControl->readBool(CoronalSliceDefault)) {
		Fl_Check_Button * checkBox;
		checkBox = userInterface->yImagePlaneCheckBox;
		checkBox->value(1);
		viewYImagePlane(1);
    }
    if (globalControl->readBool(SagittalSliceDefault)) {
		Fl_Check_Button * checkBox;
		checkBox = userInterface->xImagePlaneCheckBox;
		checkBox->value(1);
		viewXImagePlane(1);
    }

    if (globalControl->readBool(DrawBoundary))
		toggleDrawBoundary(1);
	else
		toggleDrawBoundary(0);

	// Set display of standard axes
    if (globalControl->readBool(DisplayStdAxes)) {
		toggleStdAxes(true);
        userInterface->stdAxesCheckBox->value(1);
    }
    else {
		toggleStdAxes(false);
        userInterface->stdAxesCheckBox->value(0);
    }

    // Indicate current machine's byte order in the preferences window
	userInterface->byteOrder_Caret0->hide();
	userInterface->byteOrder_Caret1->hide();
	int i =  userInterface->control->machineByteOrder();
	if (i == 1)
		userInterface->byteOrder_Caret0->show();
	else if (i == 0)
		userInterface->byteOrder_Caret1->show();

#ifndef PRODUCTION_VERSION
#ifdef SURFACE_VISUALIZER
	char * mapNames[] = {
		"Loaded", "Grey", " Inverse Grey", "Magenta Green"
	};

	userInterface->matchSurfaceColormapChoice->add(mapNames[0]);
	userInterface->matchSurfaceColormapChoice->add(mapNames[1]);
	userInterface->matchSurfaceColormapChoice->add(mapNames[2]);
	userInterface->matchSurfaceColormapChoice->add(mapNames[3]);

#else	/* SURFACE_VISUALIZER */

	userInterface->openSurfaceVisualizationWindow->hide();
	userInterface->openLoadSurfaceColormapWindow->hide();
	userInterface->openLoadSurfaceValuesWindow->hide();

#endif	/* SURFACE_VISUALIZER */
#endif	/* PRODUCTION_VERSION */
}

#ifdef BINARY

void P3DUserInterfaceCallback::adjustInterfaceForBinary()
{
	// Edit Menu
	userInterface->editLandmarksMenuItem->hide();
	userInterface->constraintsPulldown->hide();

	// Windows Menu
	userInterface->windowsOptControlWin->hide();
	userInterface->windowsOptSettingsWin->hide();

	// Test/Debug Menu
	userInterface->debugTestConstrsWin->hide();
	userInterface->debugSubfigureTransWin->hide();

	// Preferences Editor Window
	userInterface->prefsEditorConstraintsTab->clear();
	userInterface->prefsEditorTabsGrp->remove(userInterface->prefsEditorConstraintsTab);

	// Display Window
	userInterface->blendFiguresBtn->hide();
}

#else

void P3DUserInterfaceCallback::adjustInterfaceForBinary()
{
	// File Menu
	userInterface->loadDistMap->hide();
	userInterface->imageToDistanceMap->hide();
#ifndef LM_METHOD_OBJ
	// LATER: never hide these
	userInterface->loadLandMarkModel->hide();
	userInterface->clearLandMarkModel->hide();
#endif

	// Tools Menu
	userInterface->runScriptMenuItem->hide();

	// Test/Debug Menu
	userInterface->alignModel0Moment->hide();
}

#endif

#define rotateLabel	FL_FREE_LABELTYPE

// This function is called from P3DUserInterface::show().
// See button_labels.cpp for important details.
void P3DUserInterfaceCallback::show()
{
    if(userInterface == NULL)
        return;

	Fl::set_labeltype(rotateLabel, (Fl_Label_Draw_F *) drawRotateLabel, NULL);
	userInterface->rotateLeft90->labeltype(rotateLabel);
	userInterface->rotateRight90->labeltype(rotateLabel);
	userInterface->rotateHorizontal180->labeltype(rotateLabel);
	userInterface->rotateVertical180->labeltype(rotateLabel);

    userInterface->mainWindow->show();
    userInterface->modelWindow->show();
}

void P3DUserInterfaceCallback::raiseDisplayTab(enum display_tab_t tab)
{
    if(userInterface == NULL)
        return;

    Fl_Window * displayControl = userInterface->displayControlWindow;
    if(displayControl == NULL)
        return;

	if (! userInterface->displayControlWindow->shown())
		return;

	switch(tab) {
	    case GENERAL_TAB:
			userInterface->displayControlTabs->value(userInterface->displayControlGeneralTab);	
			break;

	    case SURFACE_TAB:	
			userInterface->displayControlTabs->value(userInterface->displayControlSurfaceTab);	
			break;

	    case IMAGE_TAB:
			userInterface->displayControlTabs->value(userInterface->displayControlImageTab);			
			break;

	    case TILESET_TAB:
			userInterface->displayControlTabs->value(userInterface->displayControlTilesetTab);
			break;
	}
}

void P3DUserInterfaceCallback::newModel()
{
    P3DControl * control;
    P3DView * modelView;

    if(userInterface == NULL)
        return;

    control = userInterface->control;
    if(control == NULL)
        return;

//	finishEditModel();
	finishEditLandmarks();

    if (! askToSaveModel())
		return;    // Canceled

	modelView = userInterface->modelView;
	if(modelView != NULL) {
		control->deleteImageLandmarks();
		setLandmarksImageNumber(0);
	}
	editingLandmark = -1;

    control->newModel();

    resetPGADeformation();

#ifndef BINARY
	finishPGALoading(NULL);
#endif

    setMainWindowLabel();

	modelView->clearBoundaries();
    updateVisibilityFiguresDisplay(true);
    update();
}

/*  This function prompts the user to load a single file.  It either returns the
	filename or NULL if the user cancels the load.  The directory used will be
	selected based on the user's preferences for the type of file requested, and
	unless ignore is true, the preferences will be adjusted to reflect any change
	by the user.  The prompt string should be provided in message.  The final
	argument is the file extension pattern for fl_file_chooser().
*/
const char * P3DUserInterfaceCallback::askSingleFilename(control_keys_t startDir,
	char * message, char * pattern, bool ignore)
{
    char * filename;
	char * str;
	char * dir;

	switch (startDir) {
		case ModelDirectory:
			dir = globalControl->readString(ModelDirectory);
			if (dir[0] == '\0')
				dir = globalControl->readString(ImageDirectory);
			break;

		case ImageDirectory:
			dir = globalControl->readString(ImageDirectory);
			if (dir[0] == '\0')
				dir = globalControl->readString(ModelDirectory);
			break;

		case TileSetDirectory:
			dir = globalControl->readString(TileSetDirectory);
			if (dir[0] == '\0')
				dir = globalControl->readString(ModelDirectory);
			break;

		default:
			return NULL;	// Improper first argument
	}

	str = new char[strlen(dir) + 2];
	strcpy(str, dir);
	strcat(str, "/");

	filename = fl_file_chooser(message, pattern, str);
	delete [] str;

	if (filename == NULL || filename[0] == '\0')
		return NULL;	// User canceled the request

	if (ignore)
		return (const char *) filename;

	str = strrchr(filename, '/');
	str[0] = '\0';	// Temporarily change the fl_file_chooser() buffer

	dir = globalControl->readString(startDir);
	switch (startDir) {
		case ModelDirectory:
			if (0 != strcmp(dir, filename))
				globalControl->write(ModelDirectory, filename);
			break;

		case ImageDirectory:
			if (0 != strcmp(dir, filename)) {
				char * m_dir = globalControl->readString(ModelDirectory);
				if (0 == strcmp(m_dir, filename)) {
					if (dir[0] != '\0')
						globalControl->write(ImageDirectory, "");	// Same directory for images and models
				}
				else
					globalControl->write(ImageDirectory, filename);
			}
			break;

		case TileSetDirectory:
			if (0 != strcmp(dir, filename)) {
				char * m_dir = globalControl->readString(ModelDirectory);
				if (0 == strcmp(m_dir, filename)) {
					if (dir[0] != '\0')
						globalControl->write(TileSetDirectory, "");	// Same directory for images and models
				}
				else
					globalControl->write(TileSetDirectory, filename);
			}
			break;
	}

	str[0] = '/';	// Restore the fl_file_chooser() buffer (and filename)

	return (const char *) filename;
}

// This function is only called from P3DUserInterface.cpp
void P3DUserInterfaceCallback::loadModel()
{
	cout << "loadModel" << endl;
    P3DControl * control;
    P3DView * modelView;
	bool success;
	char * str;
	char * dir;
	int i;
	string textOutput;

    if (userInterface == NULL)
        return;

    control = userInterface->control;
    if (control == NULL)
        return;

	if (usingModalWindow > 0) {
		fl_beep(FL_BEEP_ERROR);
		return;
	}

	modelView = userInterface->modelView;

//	finishEditModel();
	finishEditLandmarks();

	if (! askToSaveModel())
		return;    // Canceled

	dir = globalControl->readString(ModelDirectory);
	if (dir[0] == '\0')
		dir = globalControl->readString(ImageDirectory);
	str = new char[strlen(dir) + 2];
	strcpy(str, dir);
	strcat(str, "/");

	Fl_File_Chooser fc(str, "*.m3d", Fl_File_Chooser::MULTI, "Load Model(s)");
	fc.callback(NULL);

    fc.show();
    while (fc.shown())
        Fl::wait();
	delete [] str;

	if (fc.count() == 0)
		return;		// User canceled request

	// Adjust ModelDirectory preference as needed
	str = (char *) fc.directory();
	dir = globalControl->readString(ModelDirectory);
	if (0 != strcmp(dir, str))
		globalControl->write(ModelDirectory, str);

	// Prepare to load the model(s)
	if (modelView != NULL) {
		control->deleteImageLandmarks();
		setLandmarksImageNumber(0);
		// Turn off blending, which can be slow
		userInterface->modelView->toggleBlendSubfigures(false);
		userInterface->blendFiguresBtn->value(0);
	}
	editingLandmark = -1;

	// Load the first model
	const char * filename = fc.value(1);
	cout << filename << endl;
#ifdef BINARY
	if (! control->read(filename, false)) {
#else
	if (! control->read(filename, true)) {
#endif
		cout << "Unable to load model: " << filename << endl;
		return;
	}

#ifdef BINARY
	char * dvec_fname = new char[strlen(filename) + strlen(".list") + 1];
	sprintf(dvec_fname, "%s.list", filename);
	// Load the distance vectors file if any
	if (readDistanceVectorsFile(dvec_fname)) {
		cout << "Loaded distance vector list file: "
			<< dvec_fname << endl;
		modelView->setDistanceVectorList(distanceVectorList);
	}
	else {
		modelView->setDistanceVectorList(NULL);
	}
#endif

	char s[3] = "s:";
	if (fc.count() == 1)
		s[0] = '\0';
	textOutput = "Loaded model";
	textOutput += s;
	textOutput += ' ';
	if (fc.count() == 1)
		cout << textOutput << filename << endl;
	else
		textOutput += baseName(filename);

	// Load any remaining models
	bool * ignore = new bool[fc.count()];
	success = true;
	for (i = 2; i <= fc.count(); i++) {
		filename = fc.value(i);
		if (filename != NULL && 0 != strlen(filename)) {

#ifdef BINARY
			if (control->readAppend(filename, false, ignore[i - 1]))
#else
			if (control->readAppend(filename, true, ignore[i - 1]))
#endif
			{
				textOutput += ", ";
				textOutput += baseName(filename);
			}
			else
				success = false;
		}
		else {
			success = false;
			ignore[i - 1] = false;
		}
	}
	if (fc.count() > 1)
		cout << textOutput << "\n    from directory " << fc.directory() << endl;

	if (! success)
		cout << "One or more of the specified files did not load" << endl;

#ifndef BINARY
	finishPGALoading(filename);
#endif
	for (i = 2; i <= fc.count(); i++)
		if (ignore[i])
			cout << "Warning: ignored PGA statistics in " << fc.value(i) << endl;
	delete [] ignore;

	setMainWindowLabel();
	if (modelView != NULL)
		modelView->clearBoundaries();
	update();
	updateVisibilityFiguresDisplay(true);
	setInterfaceModelPath(fc.value(fc.count()));	// Use last model's name
	visibilityAllOn();
	raiseDisplayTab(SURFACE_TAB);

    resetPGADeformation();
}

// This function is used for command-line parsing in pablo.cpp
bool P3DUserInterfaceCallback::loadModel(const char * filename)
{
	cout << "LoadModel2" << endl;
    P3DControl * control;
    P3DView * modelView;

    if (userInterface == NULL)
        return false;

    if (filename == NULL || 0 == strlen(filename))
        return false;

    control = userInterface->control;
    if (control == NULL)
        return false;

	modelView = userInterface->modelView;

	finishEditModel();
	if (modelView != NULL) {
		control->deleteImageLandmarks();
		setLandmarksImageNumber(0);
		// Turn off blending, which can be slow
		userInterface->modelView->toggleBlendSubfigures(false);
		userInterface->blendFiguresBtn->value(0);
	}
	finishEditLandmarks();
	editingLandmark = -1;

	cout << "Modified" << endl;
	if (! control->read(filename, false)) {
		cout << "Unable to load model: " << filename << endl;
        return false;
	}

	if (globalVerbosity > -1)
		cout << "Loaded model " << filename << endl;

#ifdef LM_METHOD_OBJ
	control->binaryPabloListLandmarks();
#endif

#ifndef BINARY
	finishPGALoading(filename);
#endif

    cout << "Loaded model " << filename << endl;
    setMainWindowLabel();
	if (modelView != NULL)
		modelView->clearBoundaries();
    update();
    updateVisibilityFiguresDisplay(true);
    setInterfaceModelPath(filename);
    visibilityAllOn();
    raiseDisplayTab(SURFACE_TAB);

    resetPGADeformation();
    return true;
}

#ifdef BINARY
bool P3DUserInterfaceCallback::readDistanceVectorsFile(const char * filename)
{
	if (distanceVectorList) {
		distanceVectorList->plist.clear();
		delete distanceVectorList;
		distanceVectorList = NULL;
	}

	if (! fileExists(filename))
		return false;

	distanceVectorList = new DistanceVectorList;
	DistanceVectorList::DistVectorPList & plist = distanceVectorList->plist;

	ifstream inf(filename);
	if (! inf.bad()) {
		std::copy(istream_iterator<DistanceVector>(inf),
			istream_iterator<DistanceVector>(),
			back_insert_iterator< vector<DistanceVector> >(plist));
		inf.close();
		return true;
	}
	else
		return false;
}
#endif

bool P3DUserInterfaceCallback::loadSimTransform(const char * filename)
{
    P3DControl * control;

    control = userInterface->control;
    if(control == NULL)
        return false;

	if (! control->haveModel()) {
		cout <<
			"A model must have been loaded before a similarity transformation can be loaded\n";
		// fl_beep() is not used when running a script
		return false;
	}

    if (filename != NULL && filename[0] != '\0')
		return control->loadSimTransform(filename);
	else
		return false;
}

void P3DUserInterfaceCallback::loadSimTransform()
{
    P3DControl * control;
	const char * filename;

    if (userInterface == NULL)
        return;

    control = userInterface->control;
    if (control == NULL)
        return;

	if (! control->haveModel()) {
		cout <<
			"A model must have been loaded before a similarity transformation can be loaded\n";
		fl_beep(FL_BEEP_ERROR);
		return;
	}

    filename = askSingleFilename(ModelDirectory, "Load Similarity Transformation", "{*.sim,*.m3d}|*.sim|*.m3d");
    if (filename != NULL && filename[0] != '\0') {
		if (! control->loadSimTransform(filename))
			fl_beep(FL_BEEP_ERROR);
		else {
			userInterface->applyTransMenuItem->activate();
			updateObject();
		}
	}
}

void P3DUserInterfaceCallback::applyLoadedTransformation()
{
    P3DControl * control;

    if (userInterface == NULL)
        return;

    control = userInterface->control;
    if (control == NULL)
        return;

	control->applySimTransform();
	userInterface->applyTransMenuItem->deactivate();
}

void P3DUserInterfaceCallback::deactivateCPNSSliders()
{
	userInterface->cpnsSlider1->deactivate();
	userInterface->cpnsSlider2->deactivate();
	userInterface->cpnsSlider3->deactivate();
	userInterface->cpnsSlider4->deactivate();
	userInterface->cpnsSlider5->deactivate();
}

void P3DUserInterfaceCallback::activateCPNSSliders()
{
	userInterface->cpnsSlider1->activate();
	userInterface->cpnsSlider2->activate();
	userInterface->cpnsSlider3->activate();
	userInterface->cpnsSlider4->activate();
	userInterface->cpnsSlider5->activate();
}


void P3DUserInterfaceCallback::deactivatePGASliders()
{
	userInterface->pgaSlider1->deactivate();
	userInterface->pgaSlider2->deactivate();
	userInterface->pgaSlider3->deactivate();
	userInterface->pgaSlider4->deactivate();
	userInterface->pgaSlider5->deactivate();
}

void P3DUserInterfaceCallback::activatePGASliders()
{
	userInterface->pgaSlider1->activate();
	userInterface->pgaSlider2->activate();
	userInterface->pgaSlider3->activate();
	userInterface->pgaSlider4->activate();
	userInterface->pgaSlider5->activate();
}

#ifndef BINARY
void P3DUserInterfaceCallback::finishPGALoading(const char * filename)
{
	deactivatePGASliders();

    P3DControl * control = userInterface->control;

	if (pgaFilename != NULL) {
		delete [] pgaFilename;
		pgaFilename = NULL;
//		cout << "Discarded geometric statistics" << endl;
	}
	userInterface->usePGABtn->value(0);
	userInterface->loadPGAFileBtn->deactivate();

	userInterface->pgaStageOutput->label(NULL);
	userInterface->pgaStageOutput->clear_visible_focus();
	userInterface->pgaNextBtn->deactivate();

	control->initPGA();
	if (control->havePGA() && filename != NULL)
	{
		cout << "Loaded geometric statistics" << endl;

		// When the statistics are in the model file, gray-out the PGA buttons
		pgaFilename = new char[strlen(filename) + 1];
		strcpy(pgaFilename, filename);

		userInterface->usePGABtn->value(1);
		userInterface->usePGABtn->activate();
		userInterface->pgaMeanResBtn->activate();
		activatePGASliders();

		initPGADeform();
	}
	else {
		userInterface->loadPGAFileBtn->activate();

		userInterface->usePGABtn->deactivate();
		userInterface->pgaMeanResBtn->deactivate();
	}
	userInterface->pgaFileOutput->value(pgaFilename);
}
#endif

void P3DUserInterfaceCallback::setInterfaceModelPath(const char * path) {
    int len;

    if (displayModelPath != NULL)
        delete [] displayModelPath;

	if (path == NULL) {
		displayModelPath = NULL;
		return;
	}
    len = strlen(path);

    displayModelPath = new char[2 + len];

    if (len > DISPLAY_PATH_LINE_LEN) {
        strncpy(displayModelPath, path, DISPLAY_PATH_LINE_LEN);
        displayModelPath[DISPLAY_PATH_LINE_LEN] = '\0';
        strcat(displayModelPath, "\n");
        strcat(displayModelPath, path + DISPLAY_PATH_LINE_LEN);
    }
    else
        strcpy(displayModelPath, path);
    userInterface->modelPathOutput->value(displayModelPath);
}

void P3DUserInterfaceCallback::setInterfaceImagePath(const char * path) {
    int len;

    if (displayImagePath != NULL)
        delete [] displayImagePath;

	if (path != NULL) {
		len = strlen(path);
		displayImagePath = new char[2 + len];
	}
	else {
		len = 0;
		displayImagePath = NULL;
	}

    if (len > DISPLAY_PATH_LINE_LEN) {
        strncpy(displayImagePath, path, DISPLAY_PATH_LINE_LEN);
        displayImagePath[DISPLAY_PATH_LINE_LEN] = '\0';
        strcat(displayImagePath, "\n");
        strcat(displayImagePath, path + DISPLAY_PATH_LINE_LEN);
    }
    else {
		if (path != NULL)
			strcpy(displayImagePath, path);
	}
    userInterface->imagePathOutput->value(displayImagePath);
}

void P3DUserInterfaceCallback::setInterfaceTilesetPath(const char * path) {
    int len;

    len = strlen(path);

    if (displayTilesetPath != NULL)
        delete [] displayTilesetPath;
    displayTilesetPath = new char[2 + len];

    if (len > DISPLAY_PATH_LINE_LEN) {
        strncpy(displayTilesetPath, path, DISPLAY_PATH_LINE_LEN);
        displayTilesetPath[DISPLAY_PATH_LINE_LEN] = '\0';
        strcat(displayTilesetPath, "\n");
        strcat(displayTilesetPath, path + DISPLAY_PATH_LINE_LEN);
    }
    else
        strcpy(displayTilesetPath, path);
    userInterface->tilesetPathOutput->value(displayTilesetPath);
}

#ifdef LM_METHOD_OBJ

void P3DUserInterfaceCallback::newLandmarkModel()
{
    P3DControl * control;

    if(userInterface == NULL)
        return;

    control = userInterface->control;
    if(control == NULL)
        return;

    control->newLandmarkModel();

    update();
}

void P3DUserInterfaceCallback::loadLandmarkModel()
{
    P3DControl * control;
    P3DView * modelView;
	bool success;
	char * str;
	char * dir;

    if(userInterface == NULL)
        return;

    control = userInterface->control;
    if(control == NULL)
        return;

	if (usingModalWindow > 0) {
		fl_beep(FL_BEEP_ERROR);
		return;
	}

	modelView = userInterface->modelView;

//	finishEditModel();
	finishEditLandmarks();

	// NOT NEEDED: landmarkModel is read-only
	//if (! askToSaveModel())
	//	return;    // Canceled

	dir = globalControl->readString(ModelDirectory);
	if (dir[0] == '\0')
		dir = globalControl->readString(ImageDirectory);
	str = new char[strlen(dir) + 2];
	strcpy(str, dir);
	strcat(str, "/");

	Fl_File_Chooser fc(str, "*.m3d", Fl_File_Chooser::MULTI, "Load LANDMARK Model(s)");
	fc.callback(NULL);

    fc.show();
    while (fc.shown())
        Fl::wait();
	delete [] str;

	const char * filename = fc.value(1);
	if (fc.count() == 0)
		return;		// User canceled request

	// Adjust ModelDirectory preference as needed
	//str = (char *) fc.directory();
	//dir = globalControl->readString(ModelDirectory);
	//if (0 != strcmp(dir, filename))
		//globalControl->write(ModelDirectory, fc.directory());

	// Prepare to load the model(s)
	if (modelView != NULL) {
		control->deleteImageLandmarks();
		setLandmarksImageNumber(0);
		// Turn off blending, which can be slow
		userInterface->modelView->toggleBlendSubfigures(false);
		userInterface->blendFiguresBtn->value(0);
	}
	editingLandmark = -1;

	// Load the first model
	if (! control->readLandmarkModel(filename)) {
		cout << "Unable to load landmark model: " << filename << endl;
		return;
	}

	char s[3] = "s:";
	if (fc.count() == 1)
		s[0] = '\0';
	cout << "Loaded landmark model" << s << ' ';
	if (fc.count() == 1)
		cout << filename << endl;
	else
		cout << baseName(filename);

	// Load any remaining models
	success = true;
	for (int i = 2; i <= fc.count(); i++) {
		filename = fc.value(i);
		if (filename != NULL && 0 != strlen(filename)) {

			if (control->readAppendLandmarkModel(filename))
				cout << ", " << baseName(filename);
			else
				success = false;
		}
		else
			success = false;
	}
	if (fc.count() > 1)
		cout << "\n    from directory " << fc.directory() << endl;

	if (! success)
		cout << "One or more of the specified files did not load" << endl;

    //resetPGADeformation();

	setMainWindowLabel();
	if (modelView != NULL)
		modelView->clearBoundaries();
	update();
	updateVisibilityFiguresDisplay(true);
	setInterfaceModelPath(fc.value(fc.count()));	// Use last model's name
	visibilityAllOn();
	raiseDisplayTab(SURFACE_TAB);
}

// This function will be used for command-line parsing in pablo.cpp
void P3DUserInterfaceCallback::loadLandmarkModel(const char * filename)
{
    P3DControl * control;
    P3DView * modelView;

    if(userInterface == NULL)
        return;

    if(filename == NULL || 0 == strlen(filename))
        return;

    control = userInterface->control;
    if(control == NULL)
        return;

	modelView = userInterface->modelView;

	finishEditModel();
	if(modelView != NULL) {
		control->deleteImageLandmarks();
		setLandmarksImageNumber(0);
		// Turn off blending, which can be slow
		userInterface->modelView->toggleBlendSubfigures(false);
		userInterface->blendFiguresBtn->value(0);
	}
	finishEditLandmarks();
	editingLandmark = -1;

    if (control->readLandmarkModel(filename)) {
        cout << "Loaded Landmark model " << filename << endl;
        setMainWindowLabel();
		if(modelView != NULL)
			modelView->clearBoundaries();
        update();
        updateVisibilityFiguresDisplay(true);
        setInterfaceModelPath(filename);
        visibilityAllOn();
        raiseDisplayTab(SURFACE_TAB);

		// print useful info about model just loaded
		if (globalVerbosity > -1)
			cout << "Loaded landmark model " << filename << endl;
		if (globalVerbosity > 0)
		{
			cout << " Found " << control->landmarkObject->getFigureCount() << " landmarks" << endl;

			// print the position of each landmark
			int nli = control->landmarkObject->getFigureCount();
			for (int li = 0; li < nli; li++)
			{
				const char *lmName = control->landmarkObject->getFigurePtr(li)->getName();

				Vector3D lmPosit;
				control->findLandmarkPositByName(lmPosit, control->landmarkObject, 0, (char*) lmName+3, 2);

				cout << "  lm #" << li << ": name=\"" << lmName << "\" position=<ucs>";
				lmPosit.print();
			}
		}

	}
    //resetPGADeformation();
}

#else	/* LM_METHOD_OBJ */

void P3DUserInterfaceCallback::newLandmarkModel()
{
}

void P3DUserInterfaceCallback::loadLandmarkModel()
{
}

#endif	/* LM_METHOD_OBJ */

#ifndef PRODUCTION_VERSION

void P3DUserInterfaceCallback::loadOldModel()
{
    const char * filename;
    P3DControl * control;
    P3DView * modelView;

    if(userInterface == NULL)
        return;

    control = userInterface->control;
    if(control == NULL)
        return;

	if (usingModalWindow > 0) {
		fl_beep(FL_BEEP_ERROR);
		return;
	}

	modelView = userInterface->modelView;

//	finishEditModel();
	finishEditLandmarks();

    if (! askToSaveModel())
		return;    // Canceled

    filename = askSingleFilename(ModelDirectory, "Load Old Model", "*.mod");

    if (filename != NULL && 0 != strlen(filename)) {
		if(modelView != NULL) {
			control->deleteImageLandmarks();
			setLandmarksImageNumber(0);
			// Turn off blending, which can be slow
			userInterface->modelView->toggleBlendSubfigures(false);
			userInterface->blendFiguresBtn->value(0);
		}
		editingLandmark = -1;

        if (control->readOldFile(filename)) {
            cout << "Loaded old-format model " << filename << endl;
            setMainWindowLabel();
			if(modelView != NULL)
				modelView->clearBoundaries();
            update();
            updateVisibilityFiguresDisplay(true);
            setInterfaceModelPath(filename);
            visibilityAllOn();
        }
    }
}

#else

void P3DUserInterfaceCallback::loadOldModel() {
}

#endif	/* PRODUCTION_VERSION */

void P3DUserInterfaceCallback::saveModel()
{
    char * filename;
    P3DControl * control;

    if(userInterface == NULL)
        return;

    control = userInterface->control;
    if(control == NULL)
        return;

	if (control->haveModel() == false) {
		cout << "A model must have been loaded before it can be saved\n";
		fl_beep(FL_BEEP_ERROR);
		return;
	}

    filename = control->getModelFileName();

    if (filename == NULL)
        saveModelAs();
    else {
		if (control->testReorderModel()) {
			switch (globalControl->readInt(ReorderModels)) {
				case -1:
					reorderWindowOpen = true;
					saveCanceled = false;
					userInterface->reorderPopupWindow->show();
					usingModalWindow++;
					// Must wait for a response to the modal popup window
					while (reorderWindowOpen)
						Fl::wait();
					usingModalWindow--;
					if (saveCanceled) {
						delete [] filename;
						return;
					}
					break;
				case 0:
					break;
				case 1:
					// Always reorder
					reorderModel();
					break;
			}
		}
        (void) control->write(NULL);
        cout << "Saved model " << filename << endl;
	}
    delete [] filename;
}

void P3DUserInterfaceCallback::saveModelAs()
{
    const char * filename;
    P3DControl * control;

    if(userInterface == NULL)
        return;

    control = userInterface->control;
    if(control == NULL)
        return;

	if (control->haveModel() == false) {
		cout << "A model must have been loaded before it can be saved\n";
		fl_beep(FL_BEEP_ERROR);
		return;
	}

    filename = askSingleFilename(ModelDirectory, "Save Model As", "*.m3d");
    if(filename != NULL && 0 != strlen(filename)) {
		if (control->testReorderModel()) {
			switch (globalControl->readInt(ReorderModels)) {
				case -1:
					reorderWindowOpen = true;
					saveCanceled = false;
					userInterface->reorderPopupWindow->show();
					// Must wait for a response to the modal popup window
					while (reorderWindowOpen)
						Fl::wait();
					if (saveCanceled)
						return;
					break;
				case 0:
					break;
				case 1:
					// Always reorder
					reorderModel();
					break;
			}
		}
        (void) control->write(filename);
        cout << "Saved model " << filename << endl;
    }

	// Adjust displayed path and file name to agree with P3DControl
    setInterfaceModelPath(filename);
    setMainWindowLabel();
}

void P3DUserInterfaceCallback::exportToImage()
{
    const char * filename;
    P3DControl * control;

    if(userInterface == NULL)
        return;

    control = userInterface->control;
    if(control == NULL)
        return;

	if (control->haveModel() == false) {
		cout << "A model must have been loaded before it can be exported as an image\n";
		fl_beep(FL_BEEP_ERROR);
		return;
	}

    filename = askSingleFilename(ImageDirectory, "Export To Image",
		imageFileExtensionList);
    if(filename != NULL)
        control->exportToImage(filename);
}

#ifndef PRODUCTION_VERSION

void P3DUserInterfaceCallback::exportDistanceMap()
{
    const char * filename;
    P3DControl * control;

    if(userInterface == NULL)
        return;

    control = userInterface->control;
    if(control == NULL)
        return;

	if (control->haveModel() == false) {
		cout << "A model must have been loaded before a distance map can be exported\n";
		fl_beep(FL_BEEP_ERROR);
		return;
	}

    filename = askSingleFilename(ImageDirectory, "Export Distance Map",
		imageFileExtensionList);
    if(filename != NULL)
        control->exportDistanceMap(filename);
}

#else	/* PRODUCTION_VERSION */

void P3DUserInterfaceCallback::exportDistanceMap() {
}

#endif	/* PRODUCTION_VERSION */

void P3DUserInterfaceCallback::exportToBYUTiles()
{
    const char * filename;
    P3DControl * control;
    P3DView * modelView;

    if(userInterface == NULL)
        return;

    control = userInterface->control;
    if(control == NULL)
        return;

	if (control->haveModel() == false) {
		cout << "A model must have been loaded before BYU tiles can be exported\n";
		fl_beep(FL_BEEP_ERROR);
		return;
	}

    filename = askSingleFilename(TileSetDirectory, "Export To BYU Tiles", "*.byu");
    if(filename != NULL) {
		modelView = userInterface->modelView;
		bool byuType = (globalControl->readInt(BYUOutputType) == 1 ? true : false);
		bool byuCoords = (globalControl->readInt(BYUOutputCoords) == 1 ? true : false);
		if(modelView != NULL)
			control->exportToBYUTiles(filename, byuType, byuCoords,
			    modelView->getSurfaceLevel());
		else
			control->exportToBYUTiles(filename, byuType, byuCoords);
	}
}

void P3DUserInterfaceCallback::loadImage(bool mapToWorld, bool stacked)
{
    const char * filename;

    if (stacked)
	    filename = askSingleFilename(ImageDirectory, "Load Stacked Image",
			imageFileExtensionList);
    else
	    filename = askSingleFilename(ImageDirectory, "Load Image",
			imageFileExtensionList);
    if (filename == NULL)
		return;

    loadImage(filename, mapToWorld, stacked);
}

// This function is used for command-line parsing in pablo.cpp.
// It is also a called from the loadImage() above.
void P3DUserInterfaceCallback::loadImage(const char * filename, bool mapToWorld,
	bool stacked, bool silent)
{
    P3DControl * control;
    Image3D * image;
    P3DView * modelView;
    double maxIndex, pos;
    double min, max;
    int oldX, oldY, oldZ;
	int slice;
	bool initialized;


    if(filename == NULL || 0 == strlen(filename))
        return;

    if(userInterface == NULL)
        return;

    control = userInterface->control;
    if(control == NULL)
        return;

	initialized = control->haveImage();
    oldX = (int) userInterface->xImagePlanePositionRoller->value();
    oldY = (int) userInterface->yImagePlanePositionRoller->value();
    oldZ = (int) userInterface->zImagePlanePositionRoller->value();

    // See Image.h for an explanation of stacked images
    if (control->loadImage(filename, stacked, mapToWorld)) {
		if (! silent) {
			cout << "Loaded";
			if (stacked)
				cout << " stacked";
			cout << " image " << filename << endl;
		}
        setInterfaceImagePath(filename);
    }
	else {
		fl_beep(FL_BEEP_ERROR);
		return;
	}

	finishEditLandmarks();
    modelView = userInterface->modelView;
	if(modelView != NULL) {
		control->deleteImageLandmarks();
		setLandmarksImageNumber(0);
	}
	editingLandmark = -1;

    image = control->getImagePtr();
    if(image == NULL)
        return;

	/*  The 3 image slice rollers operate by slice number, so that all slices
		can be displayed, regardless of the size of the image.  The code below
		sets the necessary parameters of the rollers for a new image.  It puts
		the minimum roller positions, slices 0, at the world coordinate origins
		specified by the input image file, which is the center of the first
		voxel.  The highest roller positions will correspond to the center of
		the last voxel and have slice indexes of the dimensions minus one.  The
		rollers are initially set to the middle of their ranges.  If a new image
		is loaded, with the same number of slices in any of the cardinal
		directions, the corresponding roller will not be reset.
	*/
    maxIndex = image->getXDim() - 1;
    if (initialized &&
		(userInterface->xImagePlanePositionRoller->maximum() == maxIndex))
			slice = oldX;
    else
    {
		// Set roller's slice range
        userInterface->xImagePlanePositionRoller->range(0, maxIndex);

		// Set initial roller value
		slice = image->getXDim()/2;
        userInterface->xImagePlanePositionRoller->value(slice);
		sprintf(xImagePlaneSliceStr, "%d", slice);
		userInterface->xImagePlaneSlice->label(xImagePlaneSliceStr);
    }

	control->setXCutPlaneSlice(slice);

	if (displayModelCoords)
		pos = control->getXCutPlaneModelPos();	// Model coordinate of slice
	else
		pos = image->imageXToWorld(slice);	// World coordinate of slice
	sprintf(xImagePlanePosnStr, "%f", pos);
	userInterface->xImagePlanePosition->label(xImagePlanePosnStr);

	// Set bar graph showing relative slice value
	userInterface->xImagePlanePositionGraph->value((float) control->getXCutPlaneModelPos());

#ifdef DEBUG
    cout << "x-roller set to " << userInterface->xImagePlanePositionRoller->value() << endl;
    cout << "x-roller min set to " << userInterface->xImagePlanePositionRoller->minimum() << endl;
    cout << "x-roller max set to " << userInterface->xImagePlanePositionRoller->maximum() << endl;
    cout << "x-roller step = " << userInterface->xImagePlanePositionRoller->step() << endl;
	cout << "getXCutPlaneModelPos() = " << control->getXCutPlaneModelPos() << endl;
	cout << "getXCutPlaneWorldPos() = " << control->getXCutPlaneWorldPos() << endl;
#endif

    maxIndex = image->getYDim() - 1;
    if (initialized &&
		(userInterface->yImagePlanePositionRoller->maximum() == maxIndex))
			slice = oldY;
    else
    {
		// Set roller's slice range
        userInterface->yImagePlanePositionRoller->range(0, maxIndex);

		// Set initial roller value
		slice = image->getYDim()/2;
        userInterface->yImagePlanePositionRoller->value(slice);
		sprintf(yImagePlaneSliceStr, "%d", slice);
		userInterface->yImagePlaneSlice->label(yImagePlaneSliceStr);
    }

	control->setYCutPlaneSlice(slice);

	if (displayModelCoords)
		pos = control->getYCutPlaneModelPos();	// Model coordinate of slice
	else
		pos = image->imageYToWorld(slice);	// World coordinate of slice
	sprintf(yImagePlanePosnStr, "%f", pos);
	userInterface->yImagePlanePosition->label(yImagePlanePosnStr);

	// Set bar graph showing relative slice value
	userInterface->yImagePlanePositionGraph->value((float) control->getYCutPlaneModelPos());

#ifdef DEBUG
    cout << "y-roller set to " << userInterface->yImagePlanePositionRoller->value() << endl;
    cout << "y-roller min set to " << userInterface->yImagePlanePositionRoller->minimum() << endl;
    cout << "y-roller max set to " << userInterface->yImagePlanePositionRoller->maximum() << endl;
    cout << "y-roller step = " << userInterface->yImagePlanePositionRoller->step() << endl;
	cout << "getYCutPlaneModelPos() = " << control->getYCutPlaneModelPos() << endl;
	cout << "getYCutPlaneWorldPos() = " << control->getYCutPlaneWorldPos() << endl;
#endif

    maxIndex = image->getZDim() - 1;
    if (initialized &&
		(userInterface->zImagePlanePositionRoller->maximum() == maxIndex))
			slice = oldZ;
    else
    {
		// Set roller's slice range
        userInterface->zImagePlanePositionRoller->range(0, maxIndex);

		// Set initial roller value
		slice = image->getZDim()/2;
        userInterface->zImagePlanePositionRoller->value(slice);
		sprintf(zImagePlaneSliceStr, "%d", slice);
		userInterface->zImagePlaneSlice->label(zImagePlaneSliceStr);
    }

	control->setZCutPlaneSlice(slice);

	if (displayModelCoords)
		pos = control->getZCutPlaneModelPos();	// Model coordinate of slice
	else
		pos = image->imageZToWorld(slice);	// World coordinate of slice
	sprintf(zImagePlanePosnStr, "%f", pos);
	userInterface->zImagePlanePosition->label(zImagePlanePosnStr);

	// Set bar graph showing relative slice value
	userInterface->zImagePlanePositionGraph->value((float) control->getZCutPlaneModelPos());

#ifdef DEBUG
    cout << "z-roller set to " << userInterface->zImagePlanePositionRoller->value() << endl;
    cout << "z-roller min set to " << userInterface->zImagePlanePositionRoller->minimum() << endl;
    cout << "z-roller max set to " << userInterface->zImagePlanePositionRoller->maximum() << endl;
    cout << "z-roller step = " << userInterface->zImagePlanePositionRoller->step() << endl;
	cout << "getZCutPlaneModelPos() = " << control->getZCutPlaneModelPos() << endl;
	cout << "getZCutPlaneWorldPos() = " << control->getZCutPlaneWorldPos() << endl;
#endif

	/*  The 2 image intensity sliders operate by intensity number, so that all
		intensities can be displayed, regardless of the range of the image.  The
		code below sets the necessary parameters of the sliders for a new image.
		The minimum intensity slider is set to the minimum actual intensity and
		the maximum intensity slider is set to the maximum actual intensity.  The
		right and left arrow keys can be used to make fine adjustments to these
		sliders, once they have been roughly set.  Loading a new image that has
		the same range will cause the slider settings to be retained.
	*/
    GreyValue lo, hi;
    image->range(lo, hi);
    min = (double) image->mapDisplayToActual(lo);
    max = (double) image->mapDisplayToActual(hi);
#ifdef DEBUG
	cout << "Actual intensity range = (" << min << ", " << max << ')' << endl;
#endif
    if (! initialized || (userInterface->intensityMinSlider->minimum() != min) ||
		(userInterface->intensityMaxSlider->maximum() != max))
    {
		// Set initial slider value and output valuator
		userInterface->intensityMinSlider->range(min, max);
		userInterface->intensityMinValuator->range(min, max);
		userInterface->intensityMinSlider->value(min);
		userInterface->intensityMinValuator->value(min);
		userInterface->intensityMaxSlider->range(min, max);
		userInterface->intensityMaxValuator->range(min, max);
		userInterface->intensityMaxSlider->value(max);
		userInterface->intensityMaxValuator->value(max);
    }

	// Trigger updating the image display.
	// This must always be done, even if the min and max don't change,
	// because it is necessary when going from a regular to a stacked
	// image, or vice versa.
	intensityWindowSliderChanged();

	// Display bar graphs
	userInterface->xImagePlanePositionGraph->show();
	userInterface->yImagePlanePositionGraph->show();
	userInterface->zImagePlanePositionGraph->show();

	raiseDisplayTab(IMAGE_TAB);

    userInterface->modelView->updateImagePlanes(true);
}

void P3DUserInterfaceCallback::saveBlur()
{
    P3DControl * control;
    const char * filename;

    if(userInterface == NULL)
        return;

    control = userInterface->control;
    if(control == NULL)
        return;

	if (control->haveImage() == false) {
		cout << "No image is loaded\n";
		fl_beep(FL_BEEP_ERROR);
		return;
	}

	if (control->imageIsBlurred() == false) {
		cout << "The image has not been blurred\n";
		fl_beep(FL_BEEP_ERROR);
		return;
	}

	filename = askSingleFilename(ImageDirectory, "Gaussian Blur Image",
		imageFileExtensionList);
	if (filename == NULL)
		return;

    if (control->saveImage(filename)) {
		if (globalVerbosity >= 0)
			cout << "Saved blurred image file: " << filename << endl;
	}
	else {
		cout << "Could not save blurred image" << endl;
		fl_beep(FL_BEEP_ERROR);
	}
}

void P3DUserInterfaceCallback::discardImage()
{
    P3DControl * control;
    P3DView * modelView;

    if(userInterface == NULL)
        return;

    control = userInterface->control;
    if(control == NULL)
        return;

	control->discardImage();
    setInterfaceImagePath(NULL);

	finishEditLandmarks();
    modelView = userInterface->modelView;
	if(modelView != NULL) {
		control->deleteImageLandmarks();
		setLandmarksImageNumber(0);
	}
	editingLandmark = -1;

    userInterface->xImagePlanePositionRoller->range(0, 100);
    userInterface->xImagePlanePositionRoller->value(50);
	xImagePlaneSliceStr[0] = '\0';
	userInterface->xImagePlaneSlice->label(xImagePlaneSliceStr);
	xImagePlanePosnStr[0] = '\0';
	userInterface->xImagePlanePosition->label(xImagePlanePosnStr);

    userInterface->yImagePlanePositionRoller->range(0, 100);
    userInterface->yImagePlanePositionRoller->value(50);
	yImagePlaneSliceStr[0] = '\0';
	userInterface->yImagePlaneSlice->label(yImagePlaneSliceStr);
	yImagePlanePosnStr[0] = '\0';
	userInterface->yImagePlanePosition->label(yImagePlanePosnStr);

    userInterface->zImagePlanePositionRoller->range(0, 100);
    userInterface->zImagePlanePositionRoller->value(50);
	zImagePlaneSliceStr[0] = '\0';
	userInterface->zImagePlaneSlice->label(zImagePlaneSliceStr);
	zImagePlanePosnStr[0] = '\0';
	userInterface->zImagePlanePosition->label(zImagePlanePosnStr);

    userInterface->intensityMinSlider->range(0, 1000);
    userInterface->intensityMinValuator->range(0, 1000);
    userInterface->intensityMinSlider->value(0);
    userInterface->intensityMinValuator->value(0);
	intensityMinStr[0] = '\0';
	userInterface->intensityMinOutput->label(intensityMinStr);

    userInterface->intensityMaxSlider->range(0, 1000);
    userInterface->intensityMaxValuator->range(0, 1000);
    userInterface->intensityMaxSlider->value(1000);
    userInterface->intensityMaxValuator->value(1000);
	intensityMaxStr[0] = '\0';
	userInterface->intensityMaxOutput->label(intensityMaxStr);

	userInterface->xImagePlanePositionGraph->hide();
	userInterface->yImagePlanePositionGraph->hide();
	userInterface->zImagePlanePositionGraph->hide();

    userInterface->displayControlImageTab->redraw();
    userInterface->modelView->updateImagePlanes(true);
	cout << "Discarded the previous image" << endl;
}

void P3DUserInterfaceCallback::loadTileSet()
{
    P3DControl * control;
    P3DView * modelView;
	bool success;
	char * str;
	char * dir;
	int i;
	bool split;
	char **filenames;
	const char * filename;
	int min, max, count;

    if (userInterface == NULL)
        return;

    control = userInterface->control;
    if (control == NULL)
        return;

    modelView = userInterface->modelView;
	if (modelView == NULL)
		return;

	dir = globalControl->readString(TileSetDirectory);
	if (dir[0] == '\0')
		dir = globalControl->readString(ModelDirectory);
	str = new char[strlen(dir) + 2];
	strcpy(str, dir);
	strcat(str, "/");

	Fl_File_Chooser fc(str, "*.{byu,til}|*.byu|*.til", Fl_File_Chooser::MULTI, "Load Tile Set(s)");
	fc.callback(NULL);

    fc.show();
    while (fc.shown())
        Fl::wait();
	delete [] str;

	count = fc.count();
	if (count == 0)
		return;		// User canceled request

	// Adjust TileSetDirectory preference as needed
	str = (char *) fc.directory();
	dir = globalControl->readString(TileSetDirectory);
	if (0 != strcmp(dir, str))
		globalControl->write(TileSetDirectory, fc.directory());

	split = false;
	if (count > 1) {
		filenames = new char *[count];
		for (int i = 1; i <= count; i++) {
			filenames[i - 1] = new char[1 + strlen(fc.value(i))];
			strcpy(filenames[i - 1], fc.value(i));	// Fl_File_Chooser requires copying
		}

		if (! control->evaluateTileSetTypes(count, (const char **) filenames, min, max))
		{
			fl_beep(FL_BEEP_ERROR);
			return;
		}
		if (min != max)
			split = true;
	}
	else
		filenames = NULL;

	// Load the first tile set
	if (count == 1)
		filename = fc.value();
	else
		filename = filenames[0];
	bool byuCoords = (globalControl->readInt(BYUOutputCoords) == 1 ? true : false);
	if (! control->loadTileSet(filename, split, byuCoords)) {
		cout << "Unable to load tile set: " << filename << endl;
		return;
	}

	char s[3] = "s:";
	if (count == 1)
		s[0] = '\0';
	cout << "Loaded tile set" << s << ' ';
	if (count == 1)
		cout << filename << endl;
	else
		cout << baseName(filename);

	// Load any remaining tile sets
	success = true;
	for (i = 2; i <= count; i++)
	{
		filename = fc.value(i);
		if (filename != NULL && 0 != strlen(filename)) {

			if (control->appendTileSet(filename, split, byuCoords))
				cout << ", " << baseName(filename);
			else
				success = false;
		}
		else
			success = false;
	}
	if (count > 1)
		cout << "\n    from directory " << fc.directory() << endl;

	if (! success) {
		cout << "One or more of the specified files did not load" << endl;
		fl_beep(FL_BEEP_ERROR);
	}

	if (count > 1)
		filename = filenames[0];
	setInterfaceTilesetPath(filename);
	modelView->setTileFileChanged();
	raiseDisplayTab(TILESET_TAB);
	redraw();

	if (count > 1) {
		for (i = 0; i < count; i++)
			delete [] filenames[i];
		delete [] filenames;
	}
}

#ifdef BINARY

void P3DUserInterfaceCallback::loadDistMap()
{
    const char * filename;

    filename = askSingleFilename(ImageDirectory, "Load Distance Map Image",
		imageFileExtensionList);
    if (filename == NULL)
		return;

    loadDistMap(filename);
}

// This function is used for command-line parsing in pablo.cpp.
// It is also a called from the loadDistMap() above.
bool P3DUserInterfaceCallback::loadDistMap(const char * filename, bool silent)
{
    P3DControl * control;
    ImagePlanes * imagePlanes;

    if(filename == NULL || 0 == strlen(filename))
        return false;

    if(userInterface == NULL)
        return false;

    control = userInterface->control;
    if(control == NULL)
        return false;

    // See Image.h for an explanation of stacked images
    if (control->loadDistMap(filename)) {
		if (! silent)
			cout << "Loaded distance map " << filename << endl;
        setInterfaceImagePath(filename);
    }
	else {
		//fl_beep(FL_BEEP_ERROR);
		return false;
	}

    imagePlanes = control->getImagePlanesPtr();
	ImageDistanceMap * distMap = control->getDistMapPtr();  // Should be set by now

    if(distMap == NULL){
		fl_beep(FL_BEEP_ERROR);
		return false;
	}

	return true;
}

#else	/* BINARY */

void P3DUserInterfaceCallback::loadDistMap()
{
}

bool P3DUserInterfaceCallback::loadDistMap(const char * filename, bool silent)
{
	return false;
}

#endif	/* BINARY */

// This is presently used for loading PCA data
void P3DUserInterfaceCallback::loadAuxTileSet()
{
    const char * filename;
    P3DControl * control;

    if(userInterface == NULL)
        return;

    control = userInterface->control;
    if(control == NULL)
        return;

    filename = askSingleFilename(TileSetDirectory, "Load Auxiliary Tile Set", "*.byu");

    if (filename == NULL)
		return;

    if (control->loadAuxTileSet(filename)) {
		P3DView * modelView;

        cout << "Loaded auxiliary tileset " << filename << '\n';
		modelView = userInterface->modelView;
		if (modelView != NULL)
			modelView->setAuxFileChanged();
		raiseDisplayTab(SURFACE_TAB);
    }
	else
		fl_beep(FL_BEEP_ERROR);
}

void P3DUserInterfaceCallback::generateTileSet()
{
    P3DControl * control;
    P3DView * modelView;

    if(userInterface == NULL)
        return;

    control = userInterface->control;
    if(control == NULL)
        return;

    modelView = userInterface->modelView;
	if (modelView == NULL)
		return;
	
	// if no image, erase tiles, else generate new tiles; update displays
	if (! control->haveImage()) {
		cout << "No image present; tile set erased\n";
		fl_beep(FL_BEEP_ERROR);
		control->clearTileSet();
		const char * filename = "";
		setInterfaceTilesetPath(filename);
	} else {
		if (! control->generateTileSet())
			return;
		const char * filename = "- generated -";
		setInterfaceTilesetPath(filename);
	}
	
	modelView->setTileFileChanged();
    raiseDisplayTab(TILESET_TAB);
    redraw();
}

void P3DUserInterfaceCallback::saveTileSetAs()
{
    const char * filename;
    P3DControl * control;

    if(userInterface == NULL)
        return;

    control = userInterface->control;
    if(control == NULL)
        return;

	if (control->haveTileSet() == false) {
		cout << "Warning: A tile set must exist before it can be saved\n";
		fl_beep(FL_BEEP_ERROR);
		return;
	}

    filename = askSingleFilename(TileSetDirectory, "Save Tile Set As", "*.byu");
    if (filename != NULL && 0 != strlen(filename)) {
		bool byuCoords = (globalControl->readInt(BYUOutputCoords) == 1 ? true : false);
        if (control->saveTileSet(filename, byuCoords)) {
			cout << "Saved tile set " << filename << endl;

			// Adjust displayed path and file name to agree with P3DControl
			setInterfaceTilesetPath(filename);
		}
		else
			fl_beep(FL_BEEP_ERROR);
    }
}

// This function is only called from P3DUserInterface.cpp
void P3DUserInterfaceCallback::undo()
{
    P3DControl * control;

    control = userInterface->control;
    if (control == NULL)
        return;

	if (! control->haveModel())
		return;

    if (control->undo()) {
		updateVisibilityFiguresDisplay(false);
		update();
	}
}

// This function is only called from P3DUserInterface.cpp
void P3DUserInterfaceCallback::redo()
{
    P3DControl * control;

    control = userInterface->control;
    if(control == NULL)
        return;

	if (! control->haveModel())
		return;

    if (control->redo()) {
		updateVisibilityFiguresDisplay(false);
		update();
	}
}

void P3DUserInterfaceCallback::enableEditingWindows()
{
    if(userInterface == NULL)
        return;

	userInterface->addQuadFigurePulldown->activate();
	userInterface->addQuadFigureRowInput->activate();
	userInterface->addQuadFigureColumnInput->activate();
	userInterface->addQuadFigureNameInput->activate();
	userInterface->addQuadFigureColorBtn->activate();
	userInterface->addQuadFigureBtn->activate();

	userInterface->attachSubfigurePulldown->activate();
	userInterface->attachProtrusionButton->activate();
	userInterface->attachIndentationButton->activate();
	userInterface->attachNextButton->activate();
	userInterface->attachCancelButton->activate();

	userInterface->detachSubfigurePulldown->activate();

	userInterface->elongationPulldown->activate();
	userInterface->elongateRowSlider->activate();
	userInterface->elongateColumnSlider->activate();

	userInterface->atomEditorPulldown->activate();
	userInterface->atomEditorThetaSlider->activate();
	userInterface->atomEditorElongationSlider->activate();

	if (userInterface->menuBar->visible())
		userInterface->menuBar->redraw();

	userInterface->addQuadFigureDlg->activate();
	if (userInterface->addQuadFigureDlg->visible())
		userInterface->addQuadFigureDlg->redraw();

	userInterface->attachSubfigureDialog->activate();
	if (userInterface->attachSubfigureDialog->visible())
		userInterface->attachSubfigureDialog->redraw();

	userInterface->elongationDialog->activate();
	if (userInterface->elongationDialog->visible())
		userInterface->elongationDialog->redraw();

	userInterface->primitiveEditorWindow->activate();
	if (userInterface->primitiveEditorWindow->visible())
		userInterface->primitiveEditorWindow->redraw();
}

void P3DUserInterfaceCallback::disableEditingWindows()
{
    if(userInterface == NULL)
        return;

	userInterface->addQuadFigurePulldown->deactivate();
	userInterface->addQuadFigureRowInput->deactivate();
	userInterface->addQuadFigureColumnInput->deactivate();
	userInterface->addQuadFigureNameInput->deactivate();
	userInterface->addQuadFigureColorBtn->deactivate();
	userInterface->addQuadFigureBtn->deactivate();

	userInterface->attachSubfigurePulldown->deactivate();
	userInterface->attachProtrusionButton->deactivate();
	userInterface->attachIndentationButton->deactivate();
	userInterface->attachNextButton->deactivate();
	userInterface->attachCancelButton->deactivate();

	userInterface->detachSubfigurePulldown->deactivate();

	userInterface->elongationPulldown->deactivate();
	userInterface->elongateRowSlider->deactivate();
	userInterface->elongateColumnSlider->deactivate();

	userInterface->atomEditorPulldown->deactivate();
	userInterface->atomEditorThetaSlider->deactivate();
	userInterface->atomEditorElongationSlider->deactivate();

	if (userInterface->menuBar->visible())
		userInterface->menuBar->redraw();

	userInterface->addQuadFigureDlg->deactivate();
	if (userInterface->addQuadFigureDlg->visible())
		userInterface->addQuadFigureDlg->redraw();

	userInterface->attachSubfigureDialog->deactivate();
	if (userInterface->attachSubfigureDialog->visible())
		userInterface->attachSubfigureDialog->redraw();

	userInterface->elongationDialog->deactivate();
	if (userInterface->elongationDialog->visible())
		userInterface->elongationDialog->redraw();

	userInterface->primitiveEditorWindow->deactivate();
	if (userInterface->primitiveEditorWindow->visible())
		userInterface->primitiveEditorWindow->redraw();
}

#ifndef BINARY

void P3DUserInterfaceCallback::enableOptimizerWindow()
{
#ifdef DEBUG
	cout << "P3DUserInterfaceCallback::enableOptimizerWindow()" << endl;
#endif
	if (optimizerWindowInUse)
		return;

	if (! userInterface->optimizerControlDialog->shown())
		return;

	userInterface->nextOptimizationBtn->activate();
	userInterface->optimizationCancelBtn->activate();
	userInterface->optimizerSetupStage->activate();
}

bool P3DUserInterfaceCallback::disableOptimizerWindow()
{
    P3DControl * control;
    OptimizerStage currentOptimizerStage;

#ifdef DEBUG
	cout << "P3DUserInterfaceCallback::disableOptimizerWindow()" << endl;
#endif
    control = userInterface->control;
    if (control == NULL)
	    return true;

    currentOptimizerStage = control->getCurrentOptimizerStage();
	if (currentOptimizerStage == NO_OPTIMIZER_STAGE)
		return true;
	if (currentOptimizerStage == OPTIMIZER_END_STAGE)
		return true;
	if (optimizerWindowInUse)
		return false;

	userInterface->nextOptimizationBtn->deactivate();
	userInterface->optimizationCancelBtn->deactivate();
	userInterface->optimizerSetupStage->deactivate();
	return true;
}

void P3DUserInterfaceCallback::enableDeformationWindow()
{
#ifdef DEBUG
	cout << "P3DUserInterfaceCallback::enableDeformationWindow()" << endl;
#endif
//	userInterface->pgaNextBtn->activate();
//	userInterface->pgaStageOutput->activate();

	if (! userInterface->pgaDialog->shown())
		return;

	initPGADeform();
}

#endif

void P3DUserInterfaceCallback::disableDeformationWindow()
{
#ifdef DEBUG
	cout << "P3DUserInterfaceCallback::disableDeformationWindow()" << endl;
#endif
	userInterface->pgaNextBtn->deactivate();
//	userInterface->pgaStageOutput->deactivate();
	userInterface->pgaMeanResBtn->deactivate();
	deactivatePGASliders();
	userInterface->pgaZeroBtn->deactivate();
}

void P3DUserInterfaceCallback::disableCPNSDeformationWindow() {

	deactivateCPNSSliders();
	userInterface->cpnsZeroBtn->deactivate();
}

void P3DUserInterfaceCallback::addQuadFigure()
{
    P3DControl * control;
    int numRows,
        numColumns;
    char * name;

    if(userInterface == NULL)
        return;

    if (! userInterface->addQuadFigureRowInput->active_r())
	    return;

    control = userInterface->control;
    if(control == NULL)
        return;

	finishEditLandmarks();
	editingLandmark = -1;

    numRows = (int) userInterface->addQuadFigureRowInput->value();
    numColumns = (int) userInterface->addQuadFigureColumnInput->value();
    name = (char *) userInterface->addQuadFigureNameInput->value();

    int figureId = control->addQuadFigure(numRows, numColumns, name, color);
    if (figureId >= 0) {
        M3DObject * object;
        M3DFigure * figurePtr;

        object = control->getObjectPtr();
        if (object != NULL) {
			P3DView * modelView;

            figurePtr = object->getFigurePtr(figureId);
            const char * figureName = figurePtr->getName();
            // Note: name can't be used here because it is changed upon reuse.
            setInterfaceFigureName(figureId, figureName, object, true, true);

			modelView = userInterface->modelView;
			if (modelView != NULL)
				modelView->clearBoundaries();
        }
	    updateVisibilityFiguresDisplay(false);
    }
    updateSurface();
}

void P3DUserInterfaceCallback::initializeAddQuadDlg() {
    color[0] = 0.0f;    // Black is not used as a figure color, but is reserved to 
    color[1] = 0.0f;    // indicate that no color has been set.  If no color is
    color[2] = 0.0f;    // specified, the figure comes out white (see function
                        // addQuadFigure() in P3DControl.cpp.
}

void P3DUserInterfaceCallback::setFigureColor(void * currentColorBox, enum Fl_Color color_code)
{
	Fl_Box * currentColor = (Fl_Box *) currentColorBox;

    // The values stored in array color here are used in addQuadFigure() above
    currentColor->labelcolor(FL_BLACK);
    switch (color_code) {
        case FL_BLACK:   color[0] = 0.0f;   // See comments in initializeAddQuadDlg() above
                         color[1] = 0.0f;
                         color[2] = 0.0f;
                         currentColor->labelcolor(FL_WHITE);
                         break;
        case FL_RED:     color[0] = 1.0f;
                         color[1] = 0.0f;
                         color[2] = 0.0f;
                         break;
        case FL_GREEN:   color[0] = 0.0f;
                         color[1] = 1.0f;
                         color[2] = 0.0f;
                         break;
        case FL_BLUE:    color[0] = 0.0f;
                         color[1] = 0.0f;
                         color[2] = 1.0f;
                         currentColor->labelcolor(FL_WHITE);
                         break;
        case FL_YELLOW:  color[0] = 1.0f;
                         color[1] = 1.0f;
                         color[2] = 0.0f;
                         break;
        case FL_MAGENTA: color[0] = 1.0f;
                         color[1] = 0.0f;
                         color[2] = 1.0f;
                         currentColor->labelcolor(FL_WHITE);
                         break;
        case FL_CYAN:    color[0] = 0.0f;
                         color[1] = 1.0f;
                         color[2] = 1.0f;
                         break;
        case FL_GRAY:    color[0] = 0.6666667f;
                         color[1] = 0.6666667f;
                         color[2] = 0.6666667f;
                         break;
    }
    currentColor->color((unsigned char) color_code);
    currentColor->redraw();
}

void editModelPropertiesCallback(void *) {
#ifdef DEBUG
	cout << "P3DUserInterfaceCallback::editModelPropertiesCallback()" << endl;
#endif
	Fl::flush();
}

void P3DUserInterfaceCallback::setFigureColor(void * currentColorBox)
{
    uchar r, g, b;

#ifdef DEBUG
	cout << "P3DUserInterfaceCallback::setFigureColor()" << endl;
#endif
	Fl_Box * currentColor = (Fl_Box *) currentColorBox;

    const enum Fl_Color color_code = (Fl_Color) FIGURAL_COLOR;
    const Fl_Color col = currentColor->color();
    Fl::get_color(col, r, g, b);
    fl_color_chooser("Figure colors", r, g, b);
    // The values stored in array color here are used in addQuadFigure() above
    color[0] = ((float) r)/255.0f;
    color[1] = ((float) g)/255.0f;
    color[2] = ((float) b)/255.0f;
    Fl::set_color(color_code, r, g, b);
    currentColor->color((unsigned char) color_code);
    if (color[0] + color[1] + color[2] < COLOR_CUTOFF)
        currentColor->labelcolor(FL_WHITE);
    else
        currentColor->labelcolor(FL_BLACK);
    currentColor->redraw();
}

void P3DUserInterfaceCallback::editModelProperties() {
    P3DControl * control;
    M3DObject * object;
    const char * name;
	const char blank = '\0';

#ifdef DEBUG
	cout << "P3DUserInterfaceCallback::editModelProperties()" << endl;
#endif
    if(userInterface == NULL)
        return;

    control = userInterface->control;
    if(control == NULL)
        return;

    object = control->getObjectPtr();
    if (object == NULL)
		return;

	embedModelFigureNameEditor(object);

    name = object->getName();
	if (name == NULL)
		name = (char *) &blank;
	userInterface->editNameInput->value(name);

	userInterface->editModelFigureNumber->value(0);
	userInterface->editModelFigureNumber->maximum(object->getFigureCount() - 1);

	changedBlend = false;

	editModelFigureNumber();

	userInterface->editModelDimensionsOutput->clear_visible_focus();

	userInterface->editModelWindow->show();
	usingModalWindow++;
}

void P3DUserInterfaceCallback::reorderModel()
{
    P3DControl * control;

    if(userInterface == NULL)
        return;

    control = userInterface->control;
    if(control == NULL)
        return;

	if (control->haveModel() == false)
		return;

	control->reorderModel();
	control->modelNeedsSaving(true);	// TEMPORARY: PGA statistics don't reorder, so the model cannot be saved
	updateVisibilityFiguresDisplay(true);
}

void P3DUserInterfaceCallback::editFigureColor() {
    P3DControl * control;
    M3DObject * object;

#ifdef DEBUG
	cout << "P3DUserInterfaceCallback::editFigureColor()" << endl;
#endif
	int figureNumber = (int) userInterface->editModelFigureNumber->value();

    control = userInterface->control;
    object = control->getObjectPtr();
	M3DFigure * figure = object->getFigurePtr(figureNumber);
	figure->setColor(color);
	control->modelNeedsSaving();

	P3DView * modelView = userInterface->modelView;
	modelView->setModelColorsChanged();

	update();
}

void P3DUserInterfaceCallback::embedModelFigureNameEditor(M3DObject * object)
{
	char str[4];
	int i;

#ifdef DEBUG
	cout << "P3DUserInterfaceCallback::embedModelFigureNameEditor()" << endl;
#endif
	int numFigures = object->getFigureCount();
	userInterface->modelNameCells = new Flvt_Edit_Cell[numFigures + numFigures];
	Flvt_Edit_Cell * cells = (Flvt_Edit_Cell *) userInterface->modelNameCells;

	int cols = 2;
	for (i = 0; i < numFigures; i++) {
		sprintf(str, "%d.", i);
		cells[i].set(str);
	}
	// The column width setting effectively limits the string length.
	cells[numFigures].setCommonWidth(MAX_FLVW_STRING_LENGTH);
	for (i = 0; i < numFigures; i++) {
		M3DFigure * figure = object->getFigurePtr(i);
		const char * name = figure->getName();
		if (name == NULL)
			cells[i + numFigures].set(" ");
		else
			cells[i + numFigures].set(name);
	}

	// This function call adds the Flvt_Edit widget to the editModelFigureNames group.
	// The column width below is the maximum that fits the group.
	modelNameEditor = flvt_cell_editor(userInterface->editModelFigureNames,
		cells, numFigures, cols, 248, 30);
}

void P3DUserInterfaceCallback::editModelFigureNumber() {
    P3DControl * control;
    M3DObject * object;
    enum Fl_Color color_code;
	float r, g, b;
	M3DFigureTreeNode * node;
	int mode;
	char str[4];
	static char editModelDimensionsStr[8];

#ifdef DEBUG
	cout << "P3DUserInterfaceCallback::editModelFigureNumber()" << endl;
#endif
	int figureNumber = (int) userInterface->editModelFigureNumber->value();

    control = userInterface->control;
    object = control->getObjectPtr();

	M3DFigure * figure = object->getFigurePtr(figureNumber);
	r = figure->getColor()[0];
	g = figure->getColor()[1];
	b = figure->getColor()[2];
    color_code = (Fl_Color) FIGURAL_COLOR;
    Fl::set_color(color_code, (uchar) (r*255.0f), (uchar) (g*255.0f), (uchar) (b*255.0f));
    userInterface->editCurrentColorOutput->color((unsigned char) color_code);
    if (r + g + b < COLOR_CUTOFF)
        userInterface->editCurrentColorOutput->labelcolor(FL_WHITE);
    else
        userInterface->editCurrentColorOutput->labelcolor(FL_BLACK);
    userInterface->editCurrentColorOutput->redraw();

	if (figureNumber > 0)
		node = object->getFigureTreeNode(figureNumber);
	else
		node = NULL;

	if (figureNumber > 0 && node->getParent() != NULL) {
		userInterface->editModelFigureBlendExtent->value(node->getBlendExtent());
		userInterface->editModelFigureBlendAmount->value(node->getBlendAmount());
		userInterface->editModelFigureBlendExtent->activate();
		userInterface->editModelFigureBlendAmount->activate();
	}
	else {
		userInterface->editModelFigureBlendExtent->value(0.5);
		userInterface->editModelFigureBlendAmount->value(0.5);
		userInterface->editModelFigureBlendExtent->deactivate();
		userInterface->editModelFigureBlendAmount->deactivate();
		// userInterface->editModelFigureTolerance must always be active.
	}
	userInterface->editModelFigureBlendGroup->redraw();

	if (figure->hasPositivePolarity()) {
		userInterface->editModelPolarity->label(editModelPolarityLabel_0);
		userInterface->editModelPolarity->value(1);
	}
	else {
		userInterface->editModelPolarity->label(editModelPolarityLabel_1);
		userInterface->editModelPolarity->value(0);
	}
	userInterface->editModelPolarityGroup->redraw();

	if (node == NULL || node->getParent() == NULL)
		userInterface->editModelAttachmentMode->hide();
	else {
		mode = (int) node->getAttachmentMode();
		userInterface->editModelAttachmentMode->value(mode - 1);
		if (mode == 2)
			userInterface->editModelAttachmentMode->label(editModelAttachmentModeLabel_1);
		else
			userInterface->editModelAttachmentMode->label(editModelAttachmentModeLabel_0);
		userInterface->editModelAttachmentMode->show();
		userInterface->editModelAttachmentMode->redraw();
	}

	userInterface->editModelFigureTolerance->value(figure->getTolerance());
	userInterface->editModelFigureTolerance->redraw();

	//
	// FIXME: Change this for tubes. - rrs
	//
	if( typeid(*figure) == typeid(M3DQuadFigure) ) {
		sprintf(editModelDimensionsStr, "%d", dynamic_cast<M3DQuadFigure*>( figure)->getRowCount());
		sprintf(str, "%d", dynamic_cast<M3DQuadFigure*>( figure)->getColumnCount());
		strcat(editModelDimensionsStr, " x ");
		strcat(editModelDimensionsStr, str);
	}
	else if( typeid(*figure) == typeid(M3DTubeFigure) ) {
		sprintf(editModelDimensionsStr, "%d", 1);
		sprintf(str, "%d", dynamic_cast<M3DTubeFigure*>( figure)->getColumnCount());
		strcat(editModelDimensionsStr, " x ");
		strcat(editModelDimensionsStr, str);
	}
	else {
		assert(0);
	}
	userInterface->editModelDimensionsOutput->label(editModelDimensionsStr);
	userInterface->editModelDimensionsOutput->redraw();
}

void P3DUserInterfaceCallback::editModelAttachmentMode(int mode)
{
    P3DControl * control;
    M3DObject * object;
	M3DFigureTreeNode * node;

#ifdef DEBUG
	cout << "P3DUserInterfaceCallback::editModelAttachmentMode()" << endl;
#endif
	if (mode == 1)
		userInterface->editModelAttachmentMode->label(editModelAttachmentModeLabel_1);
	else
		userInterface->editModelAttachmentMode->label(editModelAttachmentModeLabel_0);
	userInterface->editModelAttachmentMode->redraw();

	int figureNumber = (int) userInterface->editModelFigureNumber->value();

    control = userInterface->control;
    object = control->getObjectPtr();
	node = object->getFigureTreeNode(figureNumber);
	node->setAttachmentMode((mode ? M3DFigureTreeNode::INDENT : M3DFigureTreeNode::PROTRUDE));
	control->modelNeedsSaving();
}

void P3DUserInterfaceCallback::editModelFigureBlendExtent() {
    P3DControl * control;
    M3DObject * object;
	M3DFigureTreeNode * node;

#ifdef DEBUG
	cout << "P3DUserInterfaceCallback::editModelFigureBlendExtent()" << endl;
#endif
	int figureNumber = (int) userInterface->editModelFigureNumber->value();

    control = userInterface->control;
    object = control->getObjectPtr();
	node = object->getFigureTreeNode(figureNumber);
	node->setBlendExtent(userInterface->editModelFigureBlendExtent->value());
	control->modelNeedsSaving();
	// For speed, update blended surface only after Done is pressed
	changedBlend = true;
}

void P3DUserInterfaceCallback::editModelFigureBlendAmount() {
    P3DControl * control;
    M3DObject * object;
	M3DFigureTreeNode * node;

#ifdef DEBUG
	cout << "P3DUserInterfaceCallback::editModelFigureBlendAmount()" << endl;
#endif
	int figureNumber = (int) userInterface->editModelFigureNumber->value();

    control = userInterface->control;
    object = control->getObjectPtr();
	node = object->getFigureTreeNode(figureNumber);
	node->setBlendAmount(userInterface->editModelFigureBlendAmount->value());
	control->modelNeedsSaving();
	// For speed, update blended surface only after Done is pressed
	changedBlend = true;
}

void P3DUserInterfaceCallback::editModelFigureTolerance()
{
    P3DControl * control;
    M3DObject * object;
	M3DFigure * figure;

#ifdef DEBUG
	cout << "P3DUserInterfaceCallback::editModelFigureTolerance()" << endl;
#endif
	int figureNumber = (int) userInterface->editModelFigureNumber->value();

    control = userInterface->control;
    object = control->getObjectPtr();
	figure = object->getFigurePtr(figureNumber);
	figure->setTolerance((int) userInterface->editModelFigureTolerance->value());
	control->modelNeedsSaving();

	if (userInterface != NULL && userInterface->modelView != NULL) {
		if (userInterface->modelView->blendingIsOn()) {
			changedBlend = true;
			return;
		}
	}
	updateSurface();
}

void P3DUserInterfaceCallback::editModelReblend() {
	if (changedBlend)
		update();
	changedBlend = false;
}

void P3DUserInterfaceCallback::editModelPolarity()
{
    P3DControl * control;
    M3DObject * object;
	M3DFigure * figure;
	int polarity;

#ifdef DEBUG
	cout << "P3DUserInterfaceCallback::editModelPolarity()" << endl;
#endif
	polarity = userInterface->editModelPolarity->value();
	if (polarity) {
		userInterface->editModelPolarity->label(editModelPolarityLabel_0);
		userInterface->editModelPolarity->value(1);
	}
	else {
		userInterface->editModelPolarity->label(editModelPolarityLabel_1);
		userInterface->editModelPolarity->value(0);
	}
	userInterface->editModelPolarityGroup->redraw();

	int figureNumber = (int) userInterface->editModelFigureNumber->value();

    control = userInterface->control;
    object = control->getObjectPtr();
	figure = object->getFigurePtr(figureNumber);
	figure->setPositivePolarity((polarity ? true : false));
	control->modelNeedsSaving();
}

// This function is called when a new model is loaded or when the Done
// button of the editModelWindow is pushed.

void P3DUserInterfaceCallback::finishEditModel() {
    P3DControl * control;
    M3DObject * object;
    const char * oldName;
    const char * newName;
	const char blank = '\0';

    if(userInterface == NULL)
        return;

	// Needed the first time a model is loaded.
	// Cells is undefined until the Edit Model Window is opened.
	if (! userInterface->editModelWindow->shown())
		return;

	if (userInterface->modelNameCells == NULL)
		return;

    control = userInterface->control;
    if(control == NULL)
        return;

    object = control->getObjectPtr();
    if (object == NULL)
		return;

    oldName = object->getName();
	if (oldName == NULL)
		oldName = &blank;
	newName = userInterface->editNameInput->value();
	if (0 != strcmp(oldName, newName)) {
		object->setName(newName);
		control->modelNeedsSaving();
		setMainWindowLabel();
	}

	Flvt_Edit_Cell * cells = (Flvt_Edit_Cell *) userInterface->modelNameCells;
	int numFigures = object->getFigureCount();
	bool changed = false;
	for (int i = numFigures; i < numFigures + numFigures; i++) {
		if (cells[i].changed()) {
			M3DFigure * figure = object->getFigurePtr(i - numFigures);
			figure->setName(cells[i].textPtr());
			control->modelNeedsSaving();
			changed = true;
		}
	}
	if (changed)
		updateVisibilityFiguresDisplay(false);

	userInterface->editModelWindow->hide();
	usingModalWindow--;

	// Must remove the widget from the group before deleting it
	userInterface->editModelFigureNames->remove(modelNameEditor);
	delete modelNameEditor;
	modelNameEditor = NULL;
	delete [] (Flvt_Edit_Cell *) userInterface->modelNameCells;
	userInterface->modelNameCells = NULL;

	if (changedBlend)
		update();
}

void P3DUserInterfaceCallback::removeSelectedFigures()
{
    P3DControl * control;
//    M3DObject * object;
    P3DView * modelView;
//	M3DFigure * figure;
//	int numFigures;

    if(userInterface == NULL)
        return;

    control = userInterface->control;
    if (control == NULL)
        return;

//    object = control->getObjectPtr();
//    if (object == NULL)
	if (! control->haveModel())
		return;

    modelView = userInterface->modelView;
    if(modelView == NULL)
		return;

	if (control->getImageLandmarksCount() > 0) {
		control->deleteImageLandmarks();
		setLandmarksImageNumber(0);
		finishEditLandmarks();
		editingLandmark = -1;
/*	AGG: won't work because the remaining figureId's in P3DView::imageLandmarks
    may have to be adjusted.  Reconsider if P3DView::imageLandmarks is moved to
	a class Image3D.
		numFigures = object->getFigureCount();
		for (int i = 0; i < numFigures; i++) {
			if (modelView->anyLandmarkIsPaired(i) == true) {
				figure = object->getFigurePtr(i);
				int n = figure->getLandmarkCount();
				for (int l = 0; l < n; l++)
					modelView->deleteImageLandmark(i, l, true);
			}
		}
*/
	}

    control->removeSelectedFigures();
    updateVisibilityFiguresDisplay(true);
    update();
}

void P3DUserInterfaceCallback::selectAll()
{
    P3DControl * control;

    if(userInterface == NULL)
        return;

    control = userInterface->control;
    if(control == NULL)
        return;

	if (userInterface->editLandmarksWindow->shown()) {
		cout << "You cannot select atoms in landmark mode" << endl;
		fl_beep(FL_BEEP_ERROR);
		return;
	}

    control->selectAll();

    updateObject();
}

void P3DUserInterfaceCallback::deselectAll()
{
    P3DControl * control;

    if(userInterface == NULL)
        return;

    control = userInterface->control;
    if(control == NULL)
        return;

    control->deselectAll();

    updateObject();
}

void P3DUserInterfaceCallback::toggleAll()
{
    P3DControl * control;

    if(userInterface == NULL)
        return;

    control = userInterface->control;
    if(control == NULL)
        return;

    control->toggleAll();

    updateObject();
}

void P3DUserInterfaceCallback::setSelectionType(P3DSelectionType type)
{
    P3DControl * control;
    P3DView * modelView;

    if(userInterface == NULL)
        return;

    control = userInterface->control;
    if(control == NULL)
        return;

    modelView = userInterface->modelView;
    if(modelView != NULL)
        modelView->setSelectType(type);
}

void P3DUserInterfaceCallback::copy()
{
    P3DControl * control;

    if(userInterface == NULL)
        return;

    control = userInterface->control;
    if(control == NULL)
        return;

    control->copy();
//    updateObject();
}

void P3DUserInterfaceCallback::paste()
{
    P3DControl * control;

    if(userInterface == NULL)
        return;

    control = userInterface->control;
    if(control == NULL)
        return;

    control->paste();
	updateVisibilityFiguresDisplay(false);
    update();
}

void P3DUserInterfaceCallback::mirror()
{
    P3DControl * control;

    if(userInterface == NULL)
        return;

    control = userInterface->control;
    if(control == NULL)
        return;

    control->mirror();
    update();
}

void P3DUserInterfaceCallback::setMainWindowLabel()
{
    char * filename;
    Fl_Window * mainWindow;
    P3DControl * control;
    M3DObject * object;
    int labelLength;
    const char * name;

    if(userInterface == NULL)
        return;

    mainWindow = userInterface->mainWindow;

    if(mainWindow == NULL)
        return;

    control = userInterface->control;
    if (control == NULL)
		return;

    labelLength = strlen(MAIN_WINDOW_TITLE) + 1;

    filename = control->getModelFileName();

	if (mainWindowLabel)
		delete [] mainWindowLabel;

    object = control->getObjectPtr();
    if (object != NULL)
        name = object->getName();
    else
        name = NULL;

    if (name == NULL || 0 == strcmp(name, "noName")) {
        if (filename != NULL) {
            // Trim leading path from filename before display
            name = strrchr(filename, '/');
#ifdef _MSC_VER
            char * alt = strrchr(filename, '\\');
            if (alt > name)
                name = alt;
#endif
            if (name == NULL)
                name = filename;
            else
                name++;
			labelLength += strlen(name);
			mainWindowLabel = new char[labelLength];
			mainWindowLabel[0] = '\0';
			strcpy(mainWindowLabel, MAIN_WINDOW_TITLE);
            strcat(mainWindowLabel, name);
        }
        else {
			labelLength += strlen(NONAME);
			mainWindowLabel = new char[labelLength];
			mainWindowLabel[0] = '\0';
			strcpy(mainWindowLabel, MAIN_WINDOW_TITLE);
            strcat(mainWindowLabel, NONAME);
		}
    }
    else {
		labelLength += strlen(name);
		mainWindowLabel = new char[labelLength];
		mainWindowLabel[0] = '\0';
		strcpy(mainWindowLabel, MAIN_WINDOW_TITLE);
        strcat(mainWindowLabel, name);
	}

    if (filename != NULL)
       delete [] filename;

    mainWindow->label(mainWindowLabel);
}


void P3DUserInterfaceCallback::update()
{
    if(userInterface == NULL)
        return;

    if(userInterface->modelView != NULL)
        userInterface->modelView->update();
}

void P3DUserInterfaceCallback::updateObject()
{
    if(userInterface == NULL)
        return;

    if(userInterface->modelView != NULL)
        userInterface->modelView->updateObject();
}

void P3DUserInterfaceCallback::updateSurface()
{
    if(userInterface == NULL)
        return;

    if(userInterface->modelView != NULL)
        userInterface->modelView->updateSurface();
}

void P3DUserInterfaceCallback::redraw()
{
    if(userInterface == NULL)
        return;

    gMasterDispListIsValid = false;

    if(userInterface->modelView != NULL)
        userInterface->modelView->redraw();
}

void P3DUserInterfaceCallback::toggleUnitCube(int val)
{
    if(userInterface == NULL)
        return;

    if(userInterface->modelView == NULL)
        return;

    if(val == 0)
        userInterface->modelView->viewUnitCube(false);
    else
        userInterface->modelView->viewUnitCube(true);

    redraw();
}

void P3DUserInterfaceCallback::toggleStdAxes(int val)
{
    if(userInterface == NULL)
        return;

    if(userInterface->modelView == NULL)
        return;

    if(val == 0)
        userInterface->modelView->viewStdAxes(false);
    else
        userInterface->modelView->viewStdAxes(true);

    redraw();
}

void P3DUserInterfaceCallback::toggleRocking(int val)
{
    if(userInterface == NULL)
        return;

    if(userInterface->modelView == NULL)
        return;

    if(val == 0)
        userInterface->modelView->viewRocking(false);
    else
        userInterface->modelView->viewRocking(true);

    redraw();
}

void P3DUserInterfaceCallback::resetView()
{
    if(userInterface == NULL)
        return;

    if(userInterface->modelView == NULL)
        return;

    userInterface->modelView->resetTrackball();

    redraw();
}

void P3DUserInterfaceCallback::setAxialView()
{
    if(userInterface == NULL)
        return;

    if(userInterface->modelView == NULL)
        return;

    userInterface->modelView->setAxialView();
	setImagePlane(AXIAL, 1);

    redraw();
}

void P3DUserInterfaceCallback::setSaggitalView()
{
    if(userInterface == NULL)
        return;

    if(userInterface->modelView == NULL)
        return;

    userInterface->modelView->setSaggitalView();
	setImagePlane(SAGGITAL, 1);

    redraw();
}

void P3DUserInterfaceCallback::setCoronalView()
{
    if(userInterface == NULL)
        return;

    if(userInterface->modelView == NULL)
        return;

    userInterface->modelView->setCoronalView();
	setImagePlane(CORONAL, 1);

    redraw();
}

void P3DUserInterfaceCallback::rotateLeft90()
{
    if(userInterface == NULL)
        return;

    if(userInterface->modelView == NULL)
        return;

    userInterface->modelView->rotateLeft90();

    redraw();
}

void P3DUserInterfaceCallback::rotateRight90()
{
    if(userInterface == NULL)
        return;

    if(userInterface->modelView == NULL)
        return;

    userInterface->modelView->rotateRight90();

    redraw();
}

void P3DUserInterfaceCallback::rotateHorizontal180()
{
    if(userInterface == NULL)
        return;

    if(userInterface->modelView == NULL)
        return;

    userInterface->modelView->rotateHorizontal180();

    redraw();
}

void P3DUserInterfaceCallback::rotateVertical180()
{
    if(userInterface == NULL)
        return;

    if(userInterface->modelView == NULL)
        return;

    userInterface->modelView->rotateVertical180();

    redraw();
}

// draw menu of things we know how to center upon.
// also acts as an inventory of all objects in pablo.
void P3DUserInterfaceCallback::setCenterThing()
{
    if(userInterface == NULL)
        return;
	
    if(userInterface->modelView == NULL)
        return;
	
	P3DControl *control = userInterface->control;
	if (!control)
		return;

	init_menu("Center on:");

	int nFigures;

	///////////////
    // MODEL
	M3DObject *object = control->getObjectPtr();
	// set nFigures to 0 if there's no objects, else limit nFigures
	// to the max number of figures that can be selected via menu
	nFigures = object ? (object->getFigureCount() & CEN_FIG_MASK) : 0;
	add_menu_item("All figures", CEN_MODEL|CEN_ALL, nFigures ? -3 : -1);
	add_menu_item("Selected figure(s)", CEN_MODEL|CEN_SELECT,
		(nFigures && object->isAnySelected()) ? -3 : -1);
	add_menu_item("Marked Atom", CEN_MODEL|CEN_MARK,
		(control->getMarkedPrimitiveId() != -1) ? -3 : -1);

	if (nFigures)
	{
		for (int figNum = 0; figNum < nFigures; figNum++) {
			M3DFigure *fig = object->getFigurePtr(figNum);
			char *item = new char[100];
			sprintf(item, "Fig %d: %.20s", figNum, fig->getName());
			add_menu_item(item, CEN_MODEL | CEN_FIG | (figNum<<CEN_FIG_SHIFT));
			delete[] item;

			// IMAGE LANDMARKS for this figure.
			// stuff into a sub-menu to keep upper level a little cleaner
			// since there can be dozens of landmarks
#ifdef LM_METHOD_OBJ
			const int nIlm = fig->getLandmarkCount();

			for (int i = 0; i < nIlm; i++) {
				const char * name = fig->getLandmarkName(i);
				if (name && *name)
				{
					char *item3 = new char[100];
					sprintf(item3, "  Fig %d Landmarks/%.20s", figNum, name);
					add_menu_item(item3,
						CEN_ILM | (figNum<<CEN_FIG_SHIFT) | (i<<CEN_ILM_SHIFT));
					delete[] item3;
				}
			}
#endif		
		}
	}


	///////////////
    // TILE
	TileSet *til = control->getTileSetPtr();
	add_menu_item("  ", -1, -2);	// white space plus divider line;
		// HACK: each menu item must be uniquely named (this has 2 spaces)
	nFigures = til ? til->figureCount() : 0;
	add_menu_item("All Tiles", CEN_TILE|CEN_ALL, nFigures ? -3 : -1);
	if (til && til->tileCount())
	{
		for (int figNum = 0; figNum < nFigures; figNum++) {
			char item[100];
			sprintf(item, " Tile Figure %d", figNum);
			add_menu_item(item, CEN_TILE|CEN_FIG|(figNum<<CEN_FIG_SHIFT));
		}
	}

	///////////////
    // UNIT CUBE
    add_menu_item("   ", -1, -2);	// white space plus divider line
		// HACK: each menu item must be uniquely named (this has 3 spaces)
	add_menu_item("Unit Cube", CEN_SPACE|CEN_ALL);
	
	///////////////
	centerThing = do_menu();
	//cout << "SELECT=" << showbase << hex << select << dec << noshowbase << endl;
	

	//if (centerThing == -1) centerThing = CEN_MODEL|CEN_ALL;	// disable CANCEL
	setCenterOfRotation();
}

// center on chosen thing
void P3DUserInterfaceCallback::setCenterOfRotation()
{
    if(userInterface == NULL)
        return;
	
    if(userInterface->modelView == NULL)
        return;

	if (centerThing != -1)
	{
		userInterface->modelView->setCenterOfRotation(centerThing);
		redraw();
	}
}

/*  This function sets the label and button value of a figure on the
    visibility control.  Note, that the name is not copied to the
    widget, but rather the pointer is copied.  Argument isRoot should be
    true when the figure is the root figure of a tree.  Argument vis
	should be the current visibility condition of the figure.
*/
void P3DUserInterfaceCallback::setInterfaceFigureName(int figureId,
	const char * name, M3DObject * object, bool isRoot, bool vis)
{
    Fl_Round_Button * tag;
    Fl_Box * number;
    bool useBold;

    switch (figureId - visibilityFirstFigureNum) {
		case 0:	tag = userInterface->Figure0;
				number = userInterface->figureNumber0;
				break;
		case 1:	tag = userInterface->Figure1;
				number = userInterface->figureNumber1;
				break;
		case 2:	tag = userInterface->Figure2;
				number = userInterface->figureNumber2;
				break;
		case 3:	tag = userInterface->Figure3;
				number = userInterface->figureNumber3;
				break;
		case 4:	tag = userInterface->Figure4;
				number = userInterface->figureNumber4;
				break;
		case 5:	tag = userInterface->Figure5;
				number = userInterface->figureNumber5;
				break;
		case 6:	tag = userInterface->Figure6;
				number = userInterface->figureNumber6;
				break;
		default:
				return;
    }

	// Set the figural name
    if (name != NULL && 0 != strlen(name))
        tag->label(name);
    else
        tag->label(empty);

	// Set the figural visibility
	tag->value(vis ? 1 : 0);

	if (object->getFigureTreeCount() < object->getFigureCount())
		useBold = true;
	else
		useBold = false;

	// Highlight the tree roots using bold and a larger font
	if (useBold && isRoot) {
		number->labelfont(FL_HELVETICA_BOLD);
		number->labelsize(14);
	}
	else {
		number->labelfont(FL_HELVETICA);
		number->labelsize(12);
	}

    userInterface->visibilityControlWindow->redraw();
}

void P3DUserInterfaceCallback::setVisibilityFigureRange()
{
    M3DObject * object;
    P3DControl * control;
    int numFigures;

    control = userInterface->control; 
    if (control == NULL)
        return;

    object = control->getObjectPtr();
    if (object == NULL)
        return;

	numFigures = object->getFigureCount();

	visibilityFirstFigureNum = (int) userInterface->visibilityFigureRoller->value();
	if (visibilityFirstFigureNum > numFigures - NUM_VISIBILITY_FIG_BTNS) {
		visibilityFirstFigureNum = numFigures - NUM_VISIBILITY_FIG_BTNS;
		// The widget's minimum value is 0, so don't need to test it here
		userInterface->visibilityFigureRoller->value(visibilityFirstFigureNum);
	}

	updateVisibilityFiguresDisplay(false);
}

/*  This function causes the correct figural names and buttons to be
    displayed on the visibility window.  If reset is true, the roller's
    position will be set to 0; a value of false, may also require this
    to be done.
*/
void P3DUserInterfaceCallback::updateVisibilityFiguresDisplay(bool reset) {
    M3DObject * object;
    M3DFigure * figure;
    P3DControl * control;
    int numFigures, i;
	bool isRoot;

    control = userInterface->control; 
    if (control == NULL)
        return;

    object = control->getObjectPtr();
    if (object == NULL) {
		setInterfaceVisibilityCount(0);
        return;
	}

	// Display the correct number of figure name buttons
	numFigures = object->getFigureCount();
	setInterfaceVisibilityCount(numFigures);

	if (reset ||
		(visibilityFirstFigureNum + NUM_VISIBILITY_FIG_BTNS - 1 >= numFigures))
	{
		// Reset the starting figure number and the roller's position
		visibilityFirstFigureNum = 0;
		userInterface->visibilityFigureRoller->value(0);
	}

	// Specify the numbers beside the figural buttons on the visibility window
	for (i = 0; i < NUM_VISIBILITY_FIG_BTNS; i++)
		sprintf(visibilityFigureNumStrs[i], "%d", visibilityFirstFigureNum + i);
	userInterface->figureNumber0->label(visibilityFigureNumStrs[0]);
	userInterface->figureNumber1->label(visibilityFigureNumStrs[1]);
	userInterface->figureNumber2->label(visibilityFigureNumStrs[2]);
	userInterface->figureNumber3->label(visibilityFigureNumStrs[3]);
	userInterface->figureNumber4->label(visibilityFigureNumStrs[4]);
	userInterface->figureNumber5->label(visibilityFigureNumStrs[5]);
	userInterface->figureNumber6->label(visibilityFigureNumStrs[6]);

	// Set the figural names on the visibility window
    for (i = 0; i < NUM_VISIBILITY_FIG_BTNS; i++) {
		int figureId = i + visibilityFirstFigureNum;
		if (figureId > numFigures)
			break;
        figure = object->getFigurePtr(figureId);
        if (figure == NULL)
            continue;
		isRoot = object->getFigureTreeNode(figureId)->getParent() ? false : true;
        setInterfaceFigureName(figureId, figure->getName(), object, isRoot,
			figure->getVisibility());
    }
}

/*  This function decides which figural buttons to display on the visibility
    window.  It also shows or hides the figural roller.
*/
void P3DUserInterfaceCallback::setInterfaceVisibilityCount(int numFigures)
{
    if (numFigures > 0) {
        userInterface->Figure0->show();
        userInterface->figureNumber0->show();
    }
    else {
        userInterface->Figure0->hide();
        userInterface->figureNumber0->hide();
    }
    if (numFigures > 1) {
        userInterface->Figure1->show();
        userInterface->figureNumber1->show();
    }
    else {
        userInterface->Figure1->hide();
        userInterface->figureNumber1->hide();
    }
    if (numFigures > 2) {
        userInterface->Figure2->show();
        userInterface->figureNumber2->show();
    }
    else {
        userInterface->Figure2->hide();
        userInterface->figureNumber2->hide();
    }
    if (numFigures > 3) {
        userInterface->Figure3->show();
        userInterface->figureNumber3->show();
    }
    else {
        userInterface->Figure3->hide();
        userInterface->figureNumber3->hide();
    }
    if (numFigures > 4) {
        userInterface->Figure4->show();
        userInterface->figureNumber4->show();
    }
    else {
        userInterface->Figure4->hide();
        userInterface->figureNumber4->hide();
    }
    if (numFigures > 5) {
        userInterface->Figure5->show();
        userInterface->figureNumber5->show();
    }
    else {
        userInterface->Figure5->hide();
        userInterface->figureNumber5->hide();
    }
    if (numFigures > 6) {
        userInterface->Figure6->show();
        userInterface->figureNumber6->show();
    }
    else {
        userInterface->Figure6->hide();
        userInterface->figureNumber6->hide();
    }

    if (numFigures > NUM_VISIBILITY_FIG_BTNS)
		userInterface->visibilityFigureRoller->show();
	else
		userInterface->visibilityFigureRoller->hide();
}

void P3DUserInterfaceCallback::initializeInterfaceVisibility() {

    bool val = globalControl->readBool(ShowMeshConnectors);
    userInterface->meshConnectorsVisible->value(val);
    userInterface->modelView->setMedialMeshVisibility(val);

    val = globalControl->readBool(ShowAtomVectors);
    userInterface->atomVectorsVisible->value(val);
    userInterface->modelView->setAtomVectorsVisibility(val);

    val = globalControl->readBool(ShowAtoms);
    userInterface->atomsVisible->value(val);
    userInterface->modelView->setAtomsVisibility(val);

    val = globalControl->readBool(ShowConstraints);
    userInterface->constraintsVisible->value(val);
    constraintsVisible(val);

	// why is globalControl not remembering the value for ShowLandmarks?
	// So forget user's intentions and just hardcode it '1'. -GST 20061108
    val = 1; //globalControl->readBool(ShowLandmarks);
    userInterface->landmarksVisible->value(val);
    landmarksVisible((int) val);

    userInterface->Figure0->show();
    userInterface->figureNumber0->show();
    userInterface->Figure1->show();
    userInterface->figureNumber1->show();
    userInterface->Figure2->show();
    userInterface->figureNumber2->show();
    userInterface->Figure3->hide();
    userInterface->figureNumber3->hide();
    userInterface->Figure4->hide();
    userInterface->figureNumber4->hide();
    userInterface->Figure5->hide();
    userInterface->figureNumber5->hide();
    userInterface->Figure6->hide();
    userInterface->figureNumber6->hide();

	visibilityFirstFigureNum = 0;
	userInterface->visibilityFigureRoller->value(0);

    userInterface->visibilityControlWindow->redraw();
}

void P3DUserInterfaceCallback::visibilityAllOn() {
    userInterface->modelView->setAllFiguresVisibility(true);
    userInterface->Figure0->value(1);
    userInterface->Figure1->value(1);
    userInterface->Figure2->value(1);
    userInterface->Figure3->value(1);
    userInterface->Figure4->value(1);
    userInterface->Figure5->value(1);
    userInterface->Figure6->value(1);
    userInterface->visibilityControlWindow->redraw();
    updateSurface();
}

void P3DUserInterfaceCallback::visibilityAllOff() {
    userInterface->modelView->setAllFiguresVisibility(false);
    userInterface->Figure0->value(0);
    userInterface->Figure1->value(0);
    userInterface->Figure2->value(0);
    userInterface->Figure3->value(0);
    userInterface->Figure4->value(0);
    userInterface->Figure5->value(0);
    userInterface->Figure6->value(0);
    userInterface->visibilityControlWindow->redraw();
    updateSurface();
}

void P3DUserInterfaceCallback::visibilityAllInverted() {
    M3DObject * object;
    P3DControl * control;
    int numFigures, figureId;
	M3DFigure * figure;

    control = userInterface->control;
    if (control == NULL)
        return;

    object = control->getObjectPtr();
    if (object == NULL)
        return;

	numFigures = object->getFigureCount();

	// Invert all figural visibilities
	for (figureId = 0; figureId < numFigures; figureId++) {
		figure = object->getFigurePtr(figureId);
		if (figure->getVisibility())
			figure->setVisibility(false);
		else
			figure->setVisibility(true);
	}

    updateVisibilityFiguresDisplay(false);
    updateSurface();
}

void P3DUserInterfaceCallback::setFigureVisibility(int buttonID, int val) {
    userInterface->modelView->setFiguralVisibility(
		buttonID + visibilityFirstFigureNum,
        (val == 1 ? true : false));
    updateSurface();
}

void P3DUserInterfaceCallback::atomVectorsVisible(int val) {
    userInterface->modelView->setAtomVectorsVisibility(val == 1 ? true : false);
    updateObject();
}

void P3DUserInterfaceCallback::meshConnectorsVisible(int val) {
    userInterface->modelView->setMedialMeshVisibility(val == 1 ? true : false);
    updateObject();
}

void P3DUserInterfaceCallback::atomsVisible(int val) {
    userInterface->modelView->setAtomsVisibility(val == 1 ? true : false);
    updateObject();
}

void P3DUserInterfaceCallback::constraintsVisible(int val) {
    if(userInterface == NULL)
        return;

    if(userInterface->modelView == NULL)
        return;

	if (val)
		userInterface->modelView->setConstraintsVisible(true);
	else
		userInterface->modelView->setConstraintsVisible(false);
}

void P3DUserInterfaceCallback::landmarksVisible(int val) {
    if(userInterface == NULL)
        return;

    if(userInterface->modelView == NULL)
        return;

#ifdef LM_METHOD_OBJ
	globalControl->write(ShowLandmarks, val);
	userInterface->modelView->viewLandmarks(val > 0);	// pass bool
	updateObject();
#else
	if (val) {
		if (! userInterface->editLandmarksWindow->shown())
			userInterface->modelView->viewLandmarks(true);
		priorLandmarkVisibility = true;
	}
	else {
		if (! userInterface->editLandmarksWindow->shown())
			userInterface->modelView->viewLandmarks(false);
		priorLandmarkVisibility = false;
	}
    updateSurface();
#endif

}

void P3DUserInterfaceCallback::setSurfaceStyle(M3DSurfaceStyle style)
{
    if(userInterface == NULL)
        return;

    if(userInterface->modelView == NULL)
        return;

    userInterface->modelView->setSurfaceStyle(style);
    updateObject();
}

#ifndef PRODUCTION_VERSION

// Called when one of the surface style buttons is pressed on the
// Surface Visualizer Control window
void P3DUserInterfaceCallback::setMatchSurfaceStyle(M3DSurfaceStyle style)
{
    if(userInterface == NULL)
        return;

    if(userInterface->modelView == NULL)
        return;

    userInterface->modelView->setMatchSurfaceStyle(style);
    updateObject();
}

#else

void P3DUserInterfaceCallback::setMatchSurfaceStyle(M3DSurfaceStyle style) {
}

#endif	/* PRODUCTION_VERSION */

void P3DUserInterfaceCallback::setSurfaceLevel(int val)
{
    if(userInterface == NULL)
        return;

    if(userInterface->modelView == NULL)
        return;

    userInterface->modelView->setSurfaceLevel(val);
    updateSurface();
}

void P3DUserInterfaceCallback::setSurfaceLineWidth(int width)
{
    if(userInterface == NULL)
        return;

    if(userInterface->modelView == NULL)
        return;

    userInterface->modelView->setSurfaceLineWidth(width);
    updateSurface();
}

void P3DUserInterfaceCallback::toggleBlendSubfigures(int val) {
    if(userInterface == NULL)
        return;

    if(userInterface->modelView == NULL)
        return;

    userInterface->modelView->toggleBlendSubfigures((val == 0 ? false : true));
    updateObject();
}

void P3DUserInterfaceCallback::toggleDrawBoundary(int val) {
    if(userInterface == NULL)
        return;

	userInterface->drawBoundaryBtn->value(val);
	userInterface->drawBoundaryBtn2->value(val);
    if(userInterface->modelView == NULL)
        return;

    userInterface->modelView->toggleDrawBoundary((val == 0 ? false : true));
    updateSurface();
}

void P3DUserInterfaceCallback::toggleDrawFigureNames(int val) {
    if(userInterface == NULL)
        return;

    if(userInterface->modelView == NULL)
        return;

    userInterface->modelView->toggleDrawFigureNames((val == 0 ? false : true));
    updateSurface();
}

void P3DUserInterfaceCallback::viewXImagePlane(int val)
{
    if(userInterface == NULL)
        return;

    if(val == 0)
    {
        if(userInterface->modelView != NULL)
            userInterface->modelView->viewXSlice(false);
    }

    else
    {
        if(userInterface->modelView != NULL)
            userInterface->modelView->viewXSlice(true);
    }

    userInterface->modelView->updateImagePlanes(true);
}

void P3DUserInterfaceCallback::viewYImagePlane(int val)
{
    if(userInterface == NULL)
        return;

    if(val == 0)
    {
        if(userInterface->modelView != NULL)
            userInterface->modelView->viewYSlice(false);
    }

    else
    {
        if(userInterface->modelView != NULL)
            userInterface->modelView->viewYSlice(true);
    }

    userInterface->modelView->updateImagePlanes(true);
}

void P3DUserInterfaceCallback::viewZImagePlane(int val)
{
    if(userInterface == NULL)
        return;

    if(val == 0)
    {
        if(userInterface->modelView != NULL)
            userInterface->modelView->viewZSlice(false);
    }

    else
    {
        if(userInterface->modelView != NULL)
            userInterface->modelView->viewZSlice(true);
    }

    userInterface->modelView->updateImagePlanes(true);
}

void P3DUserInterfaceCallback::setXImagePlanePosition()
{
    P3DControl * control;
    double pos;
    int slice;

    if (userInterface == NULL)
        return;

    control = userInterface->control;
    if (control == NULL)
        return;

    if (! control->haveImage())
        return;

    slice = (int) userInterface->xImagePlanePositionRoller->value();
	sprintf(xImagePlaneSliceStr, "%d", slice);
	userInterface->xImagePlaneSlice->label(xImagePlaneSliceStr);

    control->setXCutPlaneSlice(slice);

    imagePlaneSliderMotion = IMAGE_PLANE_COUNTER_TOP;
	userInterface->xImagePlanePositionRoller->take_focus();	// Not really needed

    userInterface->modelView->updateImagePlanes(false);
	if (displayModelCoords)
		pos = control->getXCutPlaneModelPos();	// Model coordinate of slice
	else
		pos = control->getXCutPlaneWorldPos();	// World coordinate of slice
	sprintf(xImagePlanePosnStr, "%f", pos);
	userInterface->xImagePlanePosition->label(xImagePlanePosnStr);

	// Set bar graph showing relative slice value
	userInterface->xImagePlanePositionGraph->value((float) control->getXCutPlaneModelPos());
}

void P3DUserInterfaceCallback::setYImagePlanePosition()
{
    P3DControl * control;
    double pos;
    int slice;

    if( userInterface == NULL)
        return;

    control = userInterface->control;
    if (control == NULL)
        return;

    if (! control->haveImage())
        return;

    slice = (int) userInterface->yImagePlanePositionRoller->value();
	sprintf(yImagePlaneSliceStr, "%d", slice);
	userInterface->yImagePlaneSlice->label(yImagePlaneSliceStr);

    control->setYCutPlaneSlice(slice);

    imagePlaneSliderMotion = IMAGE_PLANE_COUNTER_TOP;
	userInterface->yImagePlanePositionRoller->take_focus();	// Not really needed

    userInterface->modelView->updateImagePlanes(false);
	if (displayModelCoords)
		pos = control->getYCutPlaneModelPos();	// Model coordinate of slice
	else
		pos = control->getYCutPlaneWorldPos();	// World coordinate of slice
	sprintf(yImagePlanePosnStr, "%f", pos);
	userInterface->yImagePlanePosition->label(yImagePlanePosnStr);

	// Set bar graph showing relative slice value
	userInterface->yImagePlanePositionGraph->value((float) control->getYCutPlaneModelPos());
}

void P3DUserInterfaceCallback::setZImagePlanePosition()
{
    P3DControl * control;
    double pos;
    int slice;

    if (userInterface == NULL)
        return;

    control = userInterface->control;
    if (control == NULL)
        return;

    if (! control->haveImage())
        return;

    slice = (int) userInterface->zImagePlanePositionRoller->value();
	sprintf(zImagePlaneSliceStr, "%d", slice);
	userInterface->zImagePlaneSlice->label(zImagePlaneSliceStr);

    control->setZCutPlaneSlice(slice);

    imagePlaneSliderMotion = IMAGE_PLANE_COUNTER_TOP;
	userInterface->zImagePlanePositionRoller->take_focus();	// Not really needed

    userInterface->modelView->updateImagePlanes(false);
	if (displayModelCoords)
		pos = control->getZCutPlaneModelPos();	// Model coordinate of slice
	else
		pos = control->getZCutPlaneWorldPos();	// World coordinate of slice
	sprintf(zImagePlanePosnStr, "%f", pos);
	userInterface->zImagePlanePosition->label(zImagePlanePosnStr);

	// Set bar graph showing relative slice value
	userInterface->zImagePlanePositionGraph->value((float) control->getZCutPlaneModelPos());
}

#ifdef BadClockFn
static time_t sTime;
#else
static clock_t sTime;
#endif

bool P3DUserInterfaceCallback::imagePlaneSliderFinishedMoving()
{
    if (imagePlaneSliderMotion > 0) {
		if (imagePlaneSliderMotion == IMAGE_PLANE_COUNTER_TOP)
			// Initialize countdown
#ifdef BadClockFn
			sTime = time(NULL);
#else
			sTime = clock();
#endif
        imagePlaneSliderMotion--;
//        cout << '-' << flush;
        if (imagePlaneSliderMotion == 0) {
			// Enforce a 1 second delay
#ifdef BadClockFn
			// The clock() function doesn't seem to work on Solaris
			if ((time(NULL) - sTime) >= 1)
#else
			if (((float) (clock() - sTime))/CLOCKS_PER_SEC > 1.0f)
#endif
				return true;
			else
				imagePlaneSliderMotion = IMAGE_PLANE_COUNTER_TOP - 1;
		}
        else
            return false;
    }
    return false;
}

void P3DUserInterfaceCallback::viewTileSet(int val)
{
    P3DControl * control;

    if(userInterface == NULL)
        return;

    control = userInterface->control;
    if(control == NULL)
        return;

	// if user turns on tiles and there's nothing to view, turn it off again.
	// [just print a message; user might enable tiles before loading them.
	//   -GST 6-nov-2006]
	if (val > 0 && ! control->haveTileSet()) {
		//userInterface->viewTileSetCheckBox->value(0);
		//return;
		cout << "There are currently no tiles to draw: load or generate some." << endl;
	}

    if(userInterface->modelView != NULL)
        userInterface->modelView->viewTileSet(val == 0 ? false : true);

    redraw();
}

void P3DUserInterfaceCallback::useDefaultColor(int val)
{
    if (userInterface == NULL)
        return;

    if (userInterface->modelView != NULL)
        userInterface->modelView->useDefaultColor(val == 0 ? false : true);
}

void P3DUserInterfaceCallback::tileSetSolid(int val)
{
    if (userInterface == NULL)
        return;

    if (userInterface->modelView != NULL)
        userInterface->modelView->tileSetSolid(val == 0 ? false : true);
}

void P3DUserInterfaceCallback::tileSetMobile(int val)
{
    if (userInterface == NULL)
        return;

    P3DControl * control = userInterface->control;
    if (control == NULL)
        return;

    control->setTileSetMobile(val == 0 ? false : true);

    if (userInterface->modelView)
        userInterface->modelView->setTileSetMobile(val == 0 ? false : true);
}

void P3DUserInterfaceCallback::tileMobileToggleFigure(int figId)
{
    if (userInterface == NULL)
        return;

    P3DControl * control = userInterface->control;
    if (control == NULL)
        return;

    control->toggleTileSetMobileFigure(figId);
}

void P3DUserInterfaceCallback::setTileAlpha(double val)
{
    if(userInterface == NULL)
        return;

	// The sqrt makes the display brighter; squaring yields darker
    if(userInterface->modelView != NULL) {
		if (brightTiles > 0)
			userInterface->modelView->setTileAlpha(sqrt(val));
		else if (brightTiles < 0)
			userInterface->modelView->setTileAlpha(val*val);
		else
			userInterface->modelView->setTileAlpha(val);
	}

    redraw();
}

void P3DUserInterfaceCallback::tileSetBright(int val) {
	brightTiles = val;

    if (userInterface == NULL)
        return;

	// Trigger a redrawing of the tileset by P3DView
	setTileAlpha(userInterface->tileAlphaSlider->value());
}

void P3DUserInterfaceCallback::tileSetDark(int val) {
	brightTiles = -val;

    if (userInterface == NULL)
        return;

	// Trigger a redrawing of the tileset by P3DView
	setTileAlpha(userInterface->tileAlphaSlider->value());
}

void P3DUserInterfaceCallback::viewAuxTileSet(int val)
{
    if (userInterface == NULL)
        return;

    if(userInterface->modelView != NULL)
        userInterface->modelView->viewAuxTileSet(val == 0 ? false : true);
    redraw();
}

/////////////////////////////////////////////////////
// CLIP PLANES - limit objects drawn on one side of
// each of 6 half-planes

void P3DUserInterfaceCallback::clipToggleEnable(int plane) {
    if(userInterface == NULL)
        return;

    P3DControl * control = userInterface->control;
    if(control == NULL)
        return;

	clipPlane *clips = control->getClipPlanes();
	clips[plane].mask ^= CLIP_MASK_ENABLED;
	redraw();	
}


void P3DUserInterfaceCallback::clipToggleMask(int plane, int mask){
    if(userInterface == NULL)
        return;

    P3DControl * control = userInterface->control;
    if(control == NULL)
        return;

	clipPlane *clips = control->getClipPlanes();
	clips[plane].mask ^= mask;
	redraw();	
}


// there's prob'ly a less tedious way to implement this... :(

void P3DUserInterfaceCallback::clipAValuatorChanged(int plane)
{
    if(userInterface == NULL)
        return;

    P3DControl * control = userInterface->control;
    if(control == NULL)
        return;

	clipPlane *clips = control->getClipPlanes();
	switch(plane) {
	case 0:
		clips[plane].equation[0] = userInterface->clipA0->value();
		break;
	case 1:
		clips[plane].equation[0] = userInterface->clipA1->value();
		break;
	case 2:
		clips[plane].equation[0] = userInterface->clipA2->value();
		break;
	default:
		cout << "internal err: plane (" << plane << ") out of range" << endl;
		return;
	}
	redraw();
}

void P3DUserInterfaceCallback::clipBValuatorChanged(int plane)
{
    if(userInterface == NULL)
        return;

    P3DControl * control = userInterface->control;
    if(control == NULL)
        return;

	clipPlane *clips = control->getClipPlanes();
	switch(plane) {
	case 0:
		clips[plane].equation[1] = userInterface->clipB0->value();
		break;
	case 1:
		clips[plane].equation[1] = userInterface->clipB1->value();
		break;
	case 2:
		clips[plane].equation[1] = userInterface->clipB2->value();
		break;
	default:
		cout << "internal err: plane (" << plane << ") out of range" << endl;
		return;
	}
	redraw();
}

void P3DUserInterfaceCallback::clipCValuatorChanged(int plane)
{
    if(userInterface == NULL)
        return;

    P3DControl * control = userInterface->control;
    if(control == NULL)
        return;

	clipPlane *clips = control->getClipPlanes();
	switch(plane) {
	case 0:
		clips[plane].equation[2] = userInterface->clipC0->value();
		break;
	case 1:
		clips[plane].equation[2] = userInterface->clipC1->value();
		break;
	case 2:
		clips[plane].equation[2] = userInterface->clipC2->value();
		break;
	default:
		cout << "internal err: plane (" << plane << ") out of range" << endl;
		return;
	}
	redraw();
}

void P3DUserInterfaceCallback::clipDValuatorChanged(int plane)
{
    if(userInterface == NULL)
        return;

    P3DControl * control = userInterface->control;
    if(control == NULL)
        return;

	clipPlane *clips = control->getClipPlanes();
	switch(plane) {
	case 0:
		clips[plane].equation[3] = userInterface->clipD0->value();
		break;
	case 1:
		clips[plane].equation[3] = userInterface->clipD1->value();
		break;
	case 2:
		clips[plane].equation[3] = userInterface->clipD2->value();
		break;
	default:
		cout << "internal err: plane (" << plane << ") out of range" << endl;
		return;
	}
	redraw();
}

// set clips as calc'd elsewhere; updates GUI
void P3DUserInterfaceCallback::setClips(int plane,
										float A, float B, float C, float D)
{
	clipPlane *clips = userInterface->control->getClipPlanes();
	clips[plane].equation[0] = A;
	clips[plane].equation[1] = B;
	clips[plane].equation[2] = C;
	clips[plane].equation[3] = D;

	switch (plane)
	{
	case 0:
		userInterface->clipA0->value(A);
		userInterface->clipB0->value(B);
		userInterface->clipC0->value(C);
		userInterface->clipD0->value(D);
		break;
	case 1:
		userInterface->clipA1->value(A);
		userInterface->clipB1->value(B);
		userInterface->clipC1->value(C);
		userInterface->clipD1->value(D);
		break;
	case 2:
		userInterface->clipA2->value(A);
		userInterface->clipB2->value(B);
		userInterface->clipC2->value(C);
		userInterface->clipD2->value(D);
		break;
	}
	
	redraw();
}

// align one of the clip planes with the current model view;
//  position clip plane so the view's COG is on the plane
//  because there's likely an object there to see the effect
//  of the plane (otherwise you can't see the plane).
void P3DUserInterfaceCallback::clipCaptureView(int plane)
{
	// calc rotation of model view via trackball pose
	float *mat = userInterface->modelView->getTrackball()->getMatrix();
	Matrix4D m4(mat);
	Quat rot;
	double scale;
	Vector3D trans;
	m4.decompose(rot, scale, trans);
	rot.norm();

	// transform view's central beam (from eye to center of view)
	// by view to find normal to clip plane
	Vector3D viewDir(0,0,1);	// initial dir when view is (0,0,0,1)
	rot.rotateVector(viewDir);
	viewDir *= scale;
	viewDir += trans;

	cout << "ROT "; rot.print();
	cout << "SCALE " << scale << endl;
	cout << "TRANS "; trans.print();
	cout << "viewDir "; viewDir.print();

	// project COG onto view direction to find distance along view dir
	// to place plane inside current object; object's COG should be at
	// the center of the unit cube or the clip plane's angle will be off.
	Vector3D cog(.5, .5, .5);
	//M3DObject *obj = userInterface->control->getObjectPtr();
    //if (obj) cog = obj->getCOG();
	setClips(plane, viewDir.getX(), viewDir.getY(), viewDir.getZ(),
		-(cog*viewDir));
}

void P3DUserInterfaceCallback::clipFlipView(int plane)
{
	clipPlane *clips = userInterface->control->getClipPlanes();
	setClips(plane,
		- clips[plane].equation[0],
		- clips[plane].equation[1],
		- clips[plane].equation[2],
		- clips[plane].equation[3]);
}



P3DCutPlaneView * P3DUserInterfaceCallback::getView(cut_plane_t planeType) const
{
    if (userInterface == NULL)
        return NULL;

    P3DCutPlaneView * view;
    switch (planeType) {
        case ATOM_PLANE:
							view = userInterface->bNCutPlaneView;
                            break;
        case CREST_PLANE:
							view = userInterface->bBPerpCutPlaneView;
                            break;
        case BPERP_N_PLANE:
							view = userInterface->bPerpNCutPlaneView;
                            break;
        case STARBOARD_PLANE:
							view = userInterface->bPerpY0CutPlaneView;
                            break;
        case PORT_PLANE:
							view = userInterface->bPerpY1CutPlaneView;
                            break;
        case INVOLUTES_PLANE:
							view = userInterface->InvoluteCutPlaneView;
                            break;
    }
    return view;
}

// Causes the contents of a particular cut plane view to be updated, usually
// because the marked atom has been moved.  The view will be redrawn in the
// separate cut plane window for that view, whether open or not, and on the
// model window, if it is displayed there.
bool P3DUserInterfaceCallback::redrawCutPlane(cut_plane_t planeType)
{
	P3DCutPlaneView * view;
    P3DControl * control;

    control = userInterface->control;
    if (control == NULL)
        return false;

	view = getView(planeType);
    if (! view->rerun(control->getObjectPtr())) {
		hideCutPlanes(planeType);
        return false;
	}
	view->redraw();
	userInterface->modelView->signalRedraw();
	return true;
}

// Called to generate a particular cut plane view and to draw it in the
// separte cut plane window for that view.  The window is not shown.
bool process_cut_plane(P3DControl * control, P3DCutPlaneView * view,
                       M3DPrimitive * markedPrimitive,
                       cut_plane_t view_type, const char * label,
                       float cut_plane_intensity)
{
    // Get image pointer
    Image3D * image;

    image = control->getImagePtr();
    if(image == NULL)
        return false;

    view->init(markedPrimitive, control->getObjectPtr(), image);
    if (view->status() == false) {
        cout << label << " cut-plane view could not be generated" << endl;
        return false;
    }

    bool flag = globalControl->readBool(DrawCutPlaneBoundary);
    view->boundary(flag);

    int w = globalControl->readInt(CutPlaneBoundaryWidth);
    view->boundaryWidth(w);

    view->type((plane_t) view_type);
    view->half_planes(true);
    view->bperp_line(true);
    view->gamma(cut_plane_intensity);
    view->separate_windowing(false);

    if (! view->run(CUT_PLANE_SIZE, CUT_PLANE_SIZE))
        return false;

    view->draw();

    int dim;
    if (view->w() < view->h())
        dim = view->w();
    else
        dim = view->h();

    view->resize(view->x(), view->y(), dim, dim);   // Force the window to be square
#ifdef DEBUG
    cout << "Resized cut-plane window to: " << view->w() << " x "
        << view->h() << endl;
#endif
    return true;
}

// Called from P3DUserInterface.cpp to step the mode of a cut plane to the
// next possible setting: both, separate window only, and view in the model
// window only.
void P3DUserInterfaceCallback::incrementCutPlanesMode(cut_plane_t plane)
{
	int p;
	Fl_Button * btn;
	movable_Fl_Window * win;
	bool viewing;

	p = (int) plane;
	cutPlaneModes[p]++;
	if (cutPlaneModes[p] > 2)
		cutPlaneModes[p] = 0;

	switch (plane) {
		case ATOM_PLANE:	btn = userInterface->cutPlaneAtomModeBtn;
							win = userInterface->bNCutPlaneWindow;
							viewing = userInterface->modelView->viewingBNCutPlane();
							break;

		case CREST_PLANE:	btn = userInterface->cutPlaneCrestModeBtn;
							win = userInterface->bBPerpCutPlaneWindow;
							viewing = userInterface->modelView->viewingBBperpCutPlane();
							break;

		case BPERP_N_PLANE:	btn = userInterface->cutPlaneBPerpNModeBtn;
							win = userInterface->bPerpNCutPlaneWindow;
							viewing = userInterface->modelView->viewingBperpNCutPlane();
							break;

		case STARBOARD_PLANE:
							btn = userInterface->cutPlaneStarboardModeBtn;
							win = userInterface->bPerpY0CutPlaneWindow;
							viewing = userInterface->modelView->viewingBperpY0CutPlane();
							break;

		case PORT_PLANE:	btn = userInterface->cutPlanePortModeBtn;
							win = userInterface->bPerpY1CutPlaneWindow;
							viewing = userInterface->modelView->viewingBperpY1CutPlane();
							break;

		case INVOLUTES_PLANE:
							btn = userInterface->cutPlaneInvolutesModeBtn;
							win = userInterface->InvoluteCutPlaneWindow;
							viewing = userInterface->modelView->viewingInvolutesCutPlane();
							break;
	}

	btn->value(0);
	btn->label(cutPlaneModeLabels[cutPlaneModes[p]]);
	btn->redraw();

	// Update the separate cut plane window
	if (cutPlaneModes[p] != MODEL_WIN) {
		if (viewing)
			win->show();
	}
	else
		win->hide();

	redraw();	// Update the cut plane view in the model window
}

// Called when the cut planes mode preference is set, to reset the buttons in the
// cut plane window to the new setting.
void P3DUserInterfaceCallback::setCutPlanesMode(int p, int value)
{
	cut_plane_t plane;
	Fl_Button * btn;
	movable_Fl_Window * win;
	bool viewing;

	plane = (cut_plane_t) p;
	if (cutPlaneModes[p] == value)
		return;
	cutPlaneModes[p] = value;

	switch (plane) {
		case ATOM_PLANE:	btn = userInterface->cutPlaneAtomModeBtn;
							win = userInterface->bNCutPlaneWindow;
							viewing = userInterface->modelView->viewingBNCutPlane();
							break;

		case CREST_PLANE:	btn = userInterface->cutPlaneCrestModeBtn;
							win = userInterface->bBPerpCutPlaneWindow;
							viewing = userInterface->modelView->viewingBBperpCutPlane();
							break;

		case BPERP_N_PLANE:	btn = userInterface->cutPlaneBPerpNModeBtn;
							win = userInterface->bPerpNCutPlaneWindow;
							viewing = userInterface->modelView->viewingBperpNCutPlane();
							break;

		case STARBOARD_PLANE:
							btn = userInterface->cutPlaneStarboardModeBtn;
							win = userInterface->bPerpY0CutPlaneWindow;
							viewing = userInterface->modelView->viewingBperpY0CutPlane();
							break;

		case PORT_PLANE:	btn = userInterface->cutPlanePortModeBtn;
							win = userInterface->bPerpY1CutPlaneWindow;
							viewing = userInterface->modelView->viewingBperpY1CutPlane();
							break;

		case INVOLUTES_PLANE:
							btn = userInterface->cutPlaneInvolutesModeBtn;
							win = userInterface->InvoluteCutPlaneWindow;
							viewing = userInterface->modelView->viewingInvolutesCutPlane();
							break;
	}

	btn->label(cutPlaneModeLabels[cutPlaneModes[p]]);
	btn->redraw();

	// Update the separate cut plane window
	if (cutPlaneModes[p] != MODEL_WIN) {
		if (viewing)
			win->show();
	}
	else
		win->hide();

	redraw();	// Update the cut plane view in the model window
}

// Called to close a cut plane window
void P3DUserInterfaceCallback::hideCutPlanes(cut_plane_t plane) {
	switch (plane) {
	    case ATOM_PLANE:		userInterface->bNCutPlaneWindow->hide();
								userInterface->cutPlaneBtn0->value(0);
								userInterface->modelView->viewBNCutPlane(false);
								cutPlanesShown[ATOM_PLANE] = false;
								break;
	    case CREST_PLANE:		userInterface->bBPerpCutPlaneWindow->hide();
								userInterface->cutPlaneBtn1->value(0);
								userInterface->modelView->viewBBperpCutPlane(false);
								cutPlanesShown[CREST_PLANE] = false;
								break;
	    case BPERP_N_PLANE:		userInterface->bPerpNCutPlaneWindow->hide();
								userInterface->cutPlaneBtn2->value(0);
								userInterface->modelView->viewBperpNCutPlane(false);
								cutPlanesShown[BPERP_N_PLANE] = false;
								break;
	    case STARBOARD_PLANE:   userInterface->bPerpY0CutPlaneWindow->hide();
								userInterface->cutPlaneBtn3->value(0);
								userInterface->modelView->viewBperpY0CutPlane(false);
								cutPlanesShown[STARBOARD_PLANE] = false;
								break;
	    case PORT_PLANE:		userInterface->bPerpY1CutPlaneWindow->hide();
								userInterface->cutPlaneBtn4->value(0);
								userInterface->modelView->viewBperpY1CutPlane(false);
								cutPlanesShown[PORT_PLANE] = false;
								break;
	    case INVOLUTES_PLANE:	userInterface->InvoluteCutPlaneWindow->hide();
								userInterface->cutPlaneBtn5->value(0);
								userInterface->modelView->viewInvolutesCutPlane(false);
						        cutPlanesShown[INVOLUTES_PLANE] = false;
								break;
	}
    updateObject();
}

// Called when the "Atom" cut plane view is turned on
void P3DUserInterfaceCallback::viewBNCutPlane()
{
    P3DControl * control;
    M3DPrimitive * markedPrimitive;

    if (userInterface == NULL)
        return;

	if (cutPlanesShown[ATOM_PLANE]) {
		hideCutPlanes(ATOM_PLANE);
		return;
	}

    control = userInterface->control;
    if (control == NULL)
        return;

    markedPrimitive = control->getMarkedPrimitive();
    if (markedPrimitive == NULL) {
        cout << "Cut plane visualization requires an atom to be marked" << endl;
		fl_beep(FL_BEEP_ERROR);
        userInterface->cutPlaneBtn0->value(0);
        return;
    }

    P3DCutPlaneView * bNCutPlaneView = userInterface->bNCutPlaneView;

    const char * str = "Atom";
    if (process_cut_plane(control, bNCutPlaneView, markedPrimitive,
        ATOM_PLANE, str, cut_plane_intensity)) {
			if (getCutPlaneMode(ATOM_PLANE) != MODEL_WIN)
				userInterface->bNCutPlaneWindow->show();
			else
				userInterface->bNCutPlaneWindow->hide();
            cutPlanesShown[ATOM_PLANE] = true;
    }
    else
        userInterface->cutPlaneBtn0->value(0);

    userInterface->modelView->viewBNCutPlane(true);
    updateObject();
}

// Called when the "Crest" cut plane view is turned on
void P3DUserInterfaceCallback::viewBBPerpCutPlane()
{
    P3DControl * control;
    M3DPrimitive * markedPrimitive;

    if (userInterface == NULL)
        return;

	if (cutPlanesShown[CREST_PLANE]) {
		hideCutPlanes(CREST_PLANE);
		return;
	}

    control = userInterface->control;
    if (control == NULL)
        return;

    markedPrimitive = control->getMarkedPrimitive();
    if (markedPrimitive == NULL) {
        cout << "Cut plane visualization requires an atom to be marked" << endl;
		fl_beep(FL_BEEP_ERROR);
        userInterface->cutPlaneBtn1->value(0);
        return;
    }

    P3DCutPlaneView * bBPerpCutPlaneView = userInterface->bBPerpCutPlaneView;

    const char * str = "Crest";
    if (process_cut_plane(control, bBPerpCutPlaneView, markedPrimitive,
        CREST_PLANE, str, cut_plane_intensity)) {
			if (getCutPlaneMode(CREST_PLANE) != MODEL_WIN)
				userInterface->bBPerpCutPlaneWindow->show();
			else
				userInterface->bBPerpCutPlaneWindow->hide();
            cutPlanesShown[CREST_PLANE] = true;
    }
    else
        userInterface->cutPlaneBtn1->value(0);

    userInterface->modelView->viewBBperpCutPlane(true);
    updateObject();
}

// Called when the "Bperp-N" cut plane view is turned on
void P3DUserInterfaceCallback::viewBPerpNCutPlane()
{
    P3DControl * control;
    M3DPrimitive * markedPrimitive;

    if (userInterface == NULL)
        return;

	if (cutPlanesShown[BPERP_N_PLANE]) {
		hideCutPlanes(BPERP_N_PLANE);
		return;
	}

    control = userInterface->control;
    if (control == NULL)
        return;

    markedPrimitive = control->getMarkedPrimitive();
    if (markedPrimitive == NULL) {
        cout << "Cut plane visualization requires an atom to be marked" << endl;
		fl_beep(FL_BEEP_ERROR);
        userInterface->cutPlaneBtn2->value(0);
        return;
    }

    P3DCutPlaneView * bPerpNCutPlaneView = userInterface->bPerpNCutPlaneView;

    const char * str = "Bperp-N";
    if (process_cut_plane(control, bPerpNCutPlaneView, markedPrimitive,
        BPERP_N_PLANE, str, cut_plane_intensity)) {
			if (getCutPlaneMode(BPERP_N_PLANE) != MODEL_WIN)
	            userInterface->bPerpNCutPlaneWindow->show();
			else
	            userInterface->bPerpNCutPlaneWindow->hide();
            cutPlanesShown[BPERP_N_PLANE] = true;
    }
    else
        userInterface->cutPlaneBtn2->value(0);

    userInterface->modelView->viewBperpNCutPlane(true);
    updateObject();
}

// Called when the "Starboard" cut plane view is turned on
void P3DUserInterfaceCallback::viewBPerpY0CutPlane()
{
    P3DControl * control;
    M3DPrimitive * markedPrimitive;

    if (userInterface == NULL)
        return;

	if (cutPlanesShown[STARBOARD_PLANE]) {
		hideCutPlanes(STARBOARD_PLANE);
		return;
	}

    control = userInterface->control;
    if (control == NULL)
        return;

    markedPrimitive = control->getMarkedPrimitive();
    if (markedPrimitive == NULL) {
        cout << "Cut plane visualization requires an atom to be marked" << endl;
		fl_beep(FL_BEEP_ERROR);
        userInterface->cutPlaneBtn3->value(0);
        return;
    }

    P3DCutPlaneView * bPerpY0CutPlaneView = userInterface->bPerpY0CutPlaneView;

    const char * str = "Starboard";
    if (process_cut_plane(control, bPerpY0CutPlaneView, markedPrimitive,
        STARBOARD_PLANE, str, cut_plane_intensity)) {
			if (getCutPlaneMode(STARBOARD_PLANE) != MODEL_WIN)
	            userInterface->bPerpY0CutPlaneWindow->show();
			else
	            userInterface->bPerpY0CutPlaneWindow->hide();
            cutPlanesShown[STARBOARD_PLANE] = true;
    }
    else
        userInterface->cutPlaneBtn3->value(0);

    userInterface->modelView->viewBperpY0CutPlane(true);
    updateObject();
}

// Called when the "Port" cut plane view is turned on
void P3DUserInterfaceCallback::viewBPerpY1CutPlane()
{
    P3DControl * control;
    M3DPrimitive * markedPrimitive;

    if (userInterface == NULL)
        return;

	if (cutPlanesShown[PORT_PLANE]) {
		hideCutPlanes(PORT_PLANE);
		return;
	}

    control = userInterface->control;
    if (control == NULL)
        return;

    markedPrimitive = control->getMarkedPrimitive();
    if (markedPrimitive == NULL) {
        cout << "Cut plane visualization requires an atom to be marked" << endl;
		fl_beep(FL_BEEP_ERROR);
        userInterface->cutPlaneBtn4->value(0);
        return;
    }

    P3DCutPlaneView * bPerpY1CutPlaneView = userInterface->bPerpY1CutPlaneView;

    const char * str = "Port";
    if (process_cut_plane(control, bPerpY1CutPlaneView, markedPrimitive,
        PORT_PLANE, str, cut_plane_intensity)) {
			if (getCutPlaneMode(PORT_PLANE) != MODEL_WIN)
	            userInterface->bPerpY1CutPlaneWindow->show();
			else
	            userInterface->bPerpY1CutPlaneWindow->hide();
            cutPlanesShown[PORT_PLANE] = true;
    }
    else
        userInterface->cutPlaneBtn4->value(0);

    userInterface->modelView->viewBperpY1CutPlane(true);
    updateObject();
}

// Called when the "Involutes" cut plane view is turned on
void P3DUserInterfaceCallback::viewInvolutesCutPlane()
{
    P3DControl * control;
    M3DPrimitive * markedPrimitive;

    if (userInterface == NULL)
        return;

	if (cutPlanesShown[INVOLUTES_PLANE]) {
		hideCutPlanes(INVOLUTES_PLANE);
		return;
	}

    control = userInterface->control;
    if (control == NULL)
        return;

    markedPrimitive = control->getMarkedPrimitive();
    if (markedPrimitive == NULL) {
        cout << "Cut plane visualization requires an atom to be marked" << endl;
		fl_beep(FL_BEEP_ERROR);
        userInterface->cutPlaneBtn5->value(0);
        return;
    }

    P3DCutPlaneView * InvoluteCutPlaneView = userInterface->InvoluteCutPlaneView;

    const char * str = "Involutes";
    if (process_cut_plane(control, InvoluteCutPlaneView, markedPrimitive,
        INVOLUTES_PLANE, str, cut_plane_intensity)) {
			if (getCutPlaneMode(INVOLUTES_PLANE) != MODEL_WIN)
	            userInterface->InvoluteCutPlaneWindow->show();
			else
	            userInterface->InvoluteCutPlaneWindow->hide();
            cutPlanesShown[INVOLUTES_PLANE] = true;
    }
    else
        userInterface->cutPlaneBtn5->value(0);

    userInterface->modelView->viewInvolutesCutPlane(true);
    updateObject();
}

// Called when the Edit -> "Atom Editor..." menu item is selected
void P3DUserInterfaceCallback::showPrimitiveEditor()
{
    if(userInterface == NULL)
        return;

    if (updateAtomEditorWindow())
		userInterface->atomEditor1stMessageBtn->label(ATOM_EDITOR_STEP1);
    userInterface->primitiveEditorWindow->show();
}

bool P3DUserInterfaceCallback::updateAtomEditorWindow()
{
    P3DControl * control;
    M3DPrimitive * markedPrimitive;
    Fl_Slider * thetaSlider,
              * elongationSlider;
    Fl_Value_Output * figureIdLabel;
    Fl_Value_Output * primitiveIdLabel;
	Fl_Output * uvLabel;
    int figureId, atomId;
	int u, v;
	static char uvStr[12];

    if(userInterface == NULL)
        return false;
    control = userInterface->control;
    if(control == NULL)
        return false;

    if(userInterface->primitiveEditorWindow == NULL)
        return false;

    thetaSlider = userInterface->atomEditorThetaSlider;
    elongationSlider = userInterface->atomEditorElongationSlider;
    figureIdLabel = userInterface->atomEditorFigureOutput;
    primitiveIdLabel = userInterface->atomEditorAtomIdOutput;
	uvLabel = userInterface->atomEditorUVOutput;

    if(thetaSlider == NULL || elongationSlider == NULL ||
       figureIdLabel == NULL || primitiveIdLabel == NULL)
        return false;

    markedPrimitive = control->getMarkedPrimitive();

    if(markedPrimitive == NULL)
    {
        thetaSlider->deactivate();
        elongationSlider->deactivate();
		userInterface->atomEditorOutputGrp->hide();
		userInterface->atomEditor1stMessageBtn->show();
		userInterface->atomEditor1stMessageBtn->clear_visible_focus();
		userInterface->atomEditor2ndMessageBtn->show();
		userInterface->atomEditor2ndMessageBtn->clear_visible_focus();
    }
    else
    {
        thetaSlider->activate();
        thetaSlider->value(markedPrimitive->getTheta() * R_RADIANS_TO_DEGREES);
        if(markedPrimitive->type() == M3D_END_PRIMITIVE)
        {
            elongationSlider->value((dynamic_cast<M3DEndPrimitive*>(markedPrimitive))->getElongation());
            elongationSlider->activate();
        }
        else
        {
            elongationSlider->value(1.0);
            elongationSlider->deactivate();
        }

        control->getRelativeMarkedPrimitiveId(figureId, atomId);
        figureIdLabel->value(figureId);
        primitiveIdLabel->value(atomId);
		control->primitiveFiguralCoordinates(figureId, atomId, u, v);
		sprintf(uvStr, "(%d, %d)", u, v);
		uvLabel->value(uvStr);

		userInterface->atomEditorOutputGrp->show();
		userInterface->atomEditor1stMessageBtn->hide();
		userInterface->atomEditor2ndMessageBtn->hide();
    }
    return true;
}

void P3DUserInterfaceCallback::setTheta()
{
    P3DControl * control;
    Fl_Slider * thetaSlider;

    if(userInterface == NULL)
        return;

    thetaSlider = userInterface->atomEditorThetaSlider;
    if(thetaSlider == NULL)
        return;

    if (! thetaSlider->active_r())
	    return;

    control = userInterface->control;
    if(control == NULL)
        return;

    control->setTheta(thetaSlider->value() * R_DEGREES_TO_RADIANS);
    updateSurface();
}

void P3DUserInterfaceCallback::setElongation()
{
    P3DControl * control;
    Fl_Slider * elongationSlider;

    if(userInterface == NULL)
        return;

    elongationSlider = userInterface->atomEditorElongationSlider;
    if(elongationSlider == NULL)
        return;

    if (! elongationSlider->active_r())
	    return;

    control = userInterface->control;
    if(control == NULL)
        return;

    control->setElongation(elongationSlider->value());
    updateSurface();
}

void P3DUserInterfaceCallback::regularizeModel(int iterations, double stepsize) {
    P3DControl * control;

    control = userInterface->control;
    if(control == NULL)
        return;

	if (! control->haveModel()) {
		cout << "Regularization requires a model to be loaded" << endl;
		return;
	}

    control->regularize(stepsize, iterations);
	updateObject();
}

#ifndef PRODUCTION_VERSION

void P3DUserInterfaceCallback::interpolatePrimitive(double u, double v)
{
    P3DControl * control;

    control = userInterface->control;
    if(control == NULL)
        return;

    control->setInterpolatedPrimitive(u, v);
    updateObject();
}

void P3DUserInterfaceCallback::destroyInterpolatedPrimitive() {
    P3DControl * control;

    control = userInterface->control;
    if(control == NULL)
        return;

	control->clearInterpolatedPrimitive();
    updateObject();
}

#else	/* PRODUCTION_VERSION */

void P3DUserInterfaceCallback::interpolatePrimitive(double u, double v) {
}

void P3DUserInterfaceCallback::destroyInterpolatedPrimitive() {
}

#endif	/* PRODUCTION_VERSION */

const char unitVectorsStr[] = "unit length";
const char proportionalStr[] = "proportional B + eta";
const char smoothStr[] = "smoothed";
const char roughStr[] = "unsmoothed";
const char yesStr[] = "yes";
const char noStr[] = "no";
const char solidStr[] = "solid";
const char dottedStr[] = "dotted";
const char compressStr[] = "Convert (compress ";
const char uncompressStr[] = "Convert (uncompress ";
const char modelStr[] = "model";
const char worldStr[] = "world";
const char percentStr[] = "percent/100";
const char screenStr[] = "screen";

void P3DUserInterfaceCallback::showPreferencesEditor()
{
    float r, g, b;
    enum Fl_Color color_code;

    if (userInterface == NULL)
        return;

    // Initialize preferences on the "Display" tab
    userInterface->displayPrefsRememberPosnText->label(noStr);
	userInterface->displayPrefsRememberPosnText->clear_visible_focus();
    userInterface->displayPrefsRememberOpenText->label(noStr);
	userInterface->displayPrefsRememberOpenText->clear_visible_focus();

    if (globalControl->readBool(DisplayStdAxes)) {
        userInterface->displayPrefsStandardAxes->value(1);
        userInterface->displayPrefsStandardAxesText->label(yesStr);
    }
    else {
        userInterface->displayPrefsStandardAxes->value(0);
        userInterface->displayPrefsStandardAxesText->label(noStr);
    }
    userInterface->displayPrefsStandardAxesText->clear_visible_focus();

    if (globalControl->readBool(DrawCutPlaneBoundary)) {
        userInterface->imagePrefsCutPlaneBoundary->value(1);
        userInterface->imagePrefsCutPlaneBoundaryText->label(yesStr);
    }
    else {
        userInterface->imagePrefsCutPlaneBoundary->value(0);
        userInterface->imagePrefsCutPlaneBoundaryText->label(noStr);
    }
    userInterface->imagePrefsCutPlaneBoundaryText->clear_visible_focus();

    userInterface->imagePrefsCutPlaneBdryWidth->value(
		globalControl->readInt(CutPlaneBoundaryWidth));

	cutPlaneModePreference = globalControl->readInt(CutPlaneMode);
	userInterface->imagePrefsCutPlaneMode->label(cutPlaneModeLabels[cutPlaneModePreference]);
	if (cutPlaneModePreference == 0)
		userInterface->imagePrefsCutPlaneModeText->label(cutPlaneModeLabel_0);
	else if (cutPlaneModePreference == 2)
		userInterface->imagePrefsCutPlaneModeText->label(cutPlaneModeLabel_1);
	else
		userInterface->imagePrefsCutPlaneModeText->label(cutPlaneModeLabel_2);

    userInterface->displayPrefsRockingAngle->value(
        globalControl->readDouble(RockingAngle));

    userInterface->displayPrefsRockingIncrement->value(
        globalControl->readDouble(RockingIncrement));

    userInterface->displayPrefsRockingIncrement->value(
        globalControl->readDouble(RockingIncrement));

    // Initialize preferences on the "Atom" tab
    if (globalControl->readBool(ShowAtoms)) {
        userInterface->atomPrefsShowAtoms->value(1);
        userInterface->atomPrefsShowAtomsText->label(yesStr);
    }
    else {
        userInterface->atomPrefsShowAtoms->value(0);
        userInterface->atomPrefsShowAtomsText->label(noStr);
    }
    userInterface->atomPrefsShowAtomsText->clear_visible_focus();

    if (globalControl->readBool(FiguresColorAtoms)) {
        userInterface->atomPrefsUseFigureColors->value(1);
        userInterface->atomPrefsUseFigureColorsText->label(yesStr);
	}
    else {
        userInterface->atomPrefsUseFigureColors->value(0);
        userInterface->atomPrefsUseFigureColorsText->label(noStr);
	}
    userInterface->atomPrefsUseFigureColorsText->clear_visible_focus();

    if (globalControl->readBool(ShowAtomVectors)) {
        userInterface->atomPrefsShowAtomVectors->value(1);
        userInterface->atomPrefsShowAtomVectorsText->label(yesStr);
    }
    else {
        userInterface->atomPrefsShowAtomVectors->value(0);
        userInterface->atomPrefsShowAtomVectorsText->label(noStr);
    }
    userInterface->atomPrefsShowAtomVectorsText->clear_visible_focus();

    if (globalControl->readBool(ExtraAtomVectors)) {
        userInterface->atomPrefsExtraAtomVectors->value(1);
        userInterface->atomPrefsExtraAtomVectorsText->label(yesStr);
    }
    else {
        userInterface->atomPrefsExtraAtomVectors->value(0);
        userInterface->atomPrefsExtraAtomVectorsText->label(noStr);
    }
    userInterface->atomPrefsExtraAtomVectorsText->clear_visible_focus();

    if (globalControl->readBool(AtomVectorsType)) {
        userInterface->atomPrefsAtomVectorsType->value(1);
        userInterface->atomPrefsAtomVectorsTypeText->label(proportionalStr);
    }
    else {
        userInterface->atomPrefsAtomVectorsType->value(0);
        userInterface->atomPrefsAtomVectorsTypeText->label(unitVectorsStr);
    }
    userInterface->atomPrefsAtomVectorsTypeText->clear_visible_focus();

    userInterface->setPrefsAtomBVectorsBtns(globalControl->readInt(BVectorsType));

    userInterface->atomPrefsAtomVectorsWidth->value(
        globalControl->readInt(AtomVectorsLineWidth));

    // Initialize preferences on the "Mesh" tab
    if (globalControl->readBool(ShowMeshConnectors)) {
        userInterface->meshPrefsShowMesh->value(1);
        userInterface->meshPrefsShowMeshText->label(yesStr);
    }
    else {
        userInterface->meshPrefsShowMesh->value(0);
        userInterface->meshPrefsShowMeshText->label(noStr);
    }
    userInterface->meshPrefsShowMeshText->clear_visible_focus();

    if (globalControl->readBool(MeshConnectorsType)) {
        userInterface->meshPrefsMeshType->value(1);
        userInterface->meshPrefsMeshTypeText->label(dottedStr);
    }
    else {
        userInterface->meshPrefsMeshType->value(0);
        userInterface->meshPrefsMeshTypeText->label(solidStr);
    }
    userInterface->meshPrefsMeshTypeText->clear_visible_focus();

    userInterface->meshPrefsMeshConnectorsWidth->value(
        globalControl->readInt(MeshConnectorsLineWidth));

    r = (float) globalControl->readDouble(ConnectorsColor_R);
    g = (float) globalControl->readDouble(ConnectorsColor_G);
    b = (float) globalControl->readDouble(ConnectorsColor_B);
    color_code = (Fl_Color) CONNECTORS_COLOR;
    Fl::set_color(color_code, (uchar) (r*255.0f), (uchar) (g*255.0f), (uchar) (b*255.0f));
    userInterface->meshPrefsCurrentMeshColor->color((unsigned char) color_code);
    if (r + g + b < COLOR_CUTOFF)
        userInterface->meshPrefsCurrentMeshColor->labelcolor(FL_WHITE);
    else
        userInterface->meshPrefsCurrentMeshColor->labelcolor(FL_BLACK);

    // Initialize preferences on the "Surface" tab

#ifndef BINARY
    color_code = (Fl_Color) PARTIAL_SURFACE_COLOR;
    Fl::set_color(color_code, (uchar) (partialSurfaceColor[0]*255.0f),
        (uchar) (partialSurfaceColor[1]*255.0f), (uchar) (partialSurfaceColor[2]*255.0f));
    userInterface->constraintPrefsPartialSurfaceColor->color((unsigned char) color_code);
    if (partialSurfaceColor[0] + partialSurfaceColor[1] + partialSurfaceColor[2] < COLOR_CUTOFF)
        userInterface->constraintPrefsPartialSurfaceColor->labelcolor(FL_WHITE);
    else
        userInterface->constraintPrefsPartialSurfaceColor->labelcolor(FL_BLACK);
#endif

    userInterface->surfacePrefsSmoothness->value(globalControl->readInt(SurfaceSmoothnessDefault));

    userInterface->surfacePrefsLineWidth->value(globalControl->readInt(SurfaceLineWidth));

    if (globalControl->readBool(TwoLights)) {
		userInterface->surfacePrefsUseTwoLights->value(1);
		userInterface->surfacePrefsUseTwoLightsText->label(yesStr);
	}
	else {
		userInterface->surfacePrefsUseTwoLights->value(0);
		userInterface->surfacePrefsUseTwoLightsText->label(noStr);
	}
	userInterface->surfacePrefsUseTwoLightsText->clear_visible_focus();

    // Initialize preferences on the "Image" tab

    if (globalControl->readBool(SmoothImages)) {
        userInterface->imagePrefsSmoothImages->value(1);
        userInterface->imagePrefsSmoothImagesText->label(yesStr);
    }
    else {
        userInterface->imagePrefsSmoothImages->value(0);
        userInterface->imagePrefsSmoothImagesText->label(noStr);
    }
    userInterface->imagePrefsSmoothImagesText->clear_visible_focus();

    if (globalControl->readBool(ScaleImages)) {
        userInterface->imagePrefsScaleInput->value(1);
        userInterface->imagePrefsScaleInputText->label(yesStr);
    }
    else {
        userInterface->imagePrefsScaleInput->value(0);
        userInterface->imagePrefsScaleInputText->label(noStr);
    }
    userInterface->imagePrefsScaleInputText->clear_visible_focus();

    if (globalControl->readBool(DrawBoundary))
        userInterface->imagePrefsDrawBoundary->value(1);
    else
        userInterface->imagePrefsDrawBoundary->value(0);

    if (globalControl->readBool(DisplayModelUnits)) {
        userInterface->imagePrefsCoordsDisplayed->value(1);
        userInterface->imagePrefsCoordsDisplayedText->label(modelStr);
    }
    else {
        userInterface->imagePrefsCoordsDisplayed->value(0);
        userInterface->imagePrefsCoordsDisplayedText->label(worldStr);
    }
    userInterface->imagePrefsCoordsDisplayedText->clear_visible_focus();

    if (globalControl->readBool(DisplayRangeUnits)) {
        userInterface->imagePrefsRangeDisplayed->value(1);
        userInterface->imagePrefsRangeDisplayedText->label(screenStr);
    }
    else {
        userInterface->imagePrefsRangeDisplayed->value(0);
        userInterface->imagePrefsRangeDisplayedText->label(percentStr);
    }
    userInterface->imagePrefsRangeDisplayedText->clear_visible_focus();

	int order = globalControl->readInt(ByteOrder);
    if (order < 0 || order > 4) {
		cout << "Invalid ByteOrder preference; assuming native byte order" << endl;
		order = 1;
	}
    userInterface->setPrefsByteOrderBtns(order);

    if (globalControl->readBool(CompressImages)) {
        userInterface->imagePrefsCompressImages->value(1);
        userInterface->imagePrefsCompressImagesText->label(yesStr);
		userInterface->imagePrefsConvertImagesType->label(compressStr);
    }
    else {
        userInterface->imagePrefsCompressImages->value(0);
        userInterface->imagePrefsCompressImagesText->label(noStr);
		userInterface->imagePrefsConvertImagesType->label(uncompressStr);
    }
    userInterface->imagePrefsCompressImagesText->clear_visible_focus();

	if (globalControl->readBool(ConvertImages)) {
		userInterface->imagePrefsConvertImages->value(1);
		userInterface->imagePrefsConvertImagesText->label(yesStr);
	}
	else {
		userInterface->imagePrefsConvertImages->value(0);
		userInterface->imagePrefsConvertImagesText->label(noStr);
	}
	userInterface->imagePrefsConvertImagesText->clear_visible_focus();

	userInterface->imagePrefsImageFormat->value(globalControl->readInt(ImageFormat));

	if (globalControl->readBool(ConvertImageFormat)) {
		userInterface->imagePrefsConvertFormat->value(1);
		userInterface->imagePrefsConvertFormatText->label(yesStr);
	}
	else {
		userInterface->imagePrefsConvertFormat->value(0);
		userInterface->imagePrefsConvertFormatText->label(noStr);
	}
	userInterface->imagePrefsConvertFormatText->clear_visible_focus();

    if (globalControl->readBool(AxialSliceDefault))
        userInterface->imagePrefsZSliceDefault->value(1);
    else
        userInterface->imagePrefsZSliceDefault->value(0);

    if (globalControl->readBool(CoronalSliceDefault))
        userInterface->imagePrefsYSliceDefault->value(1);
    else
        userInterface->imagePrefsYSliceDefault->value(0);

    if (globalControl->readBool(SagittalSliceDefault))
        userInterface->imagePrefsXSliceDefault->value(1);
    else
        userInterface->imagePrefsXSliceDefault->value(0);

    // Initialize preferences on the "Tiles" tab
    r = (float) globalControl->readDouble(TilesColor_R);
    g = (float) globalControl->readDouble(TilesColor_G);
    b = (float) globalControl->readDouble(TilesColor_B);
    color_code = (Fl_Color) TILES_COLOR;
    Fl::set_color(color_code, (uchar) (r*255.0f), (uchar) (g*255.0f), (uchar) (b*255.0f));
    userInterface->tilePrefsCurrentTilesColor->color((unsigned char) color_code);
    if (r + g + b < COLOR_CUTOFF)
        userInterface->tilePrefsCurrentTilesColor->labelcolor(FL_WHITE);
    else
        userInterface->tilePrefsCurrentTilesColor->labelcolor(FL_BLACK);

#ifndef BINARY
    // Initialize preferences on the "Constraints" tab
	if (globalControl->readBool(ShowConstraints)) {
        userInterface->constraintPrefsShowConstraints->value(1);
        userInterface->constraintPrefsShowConstraintsText->label(yesStr);
	}
    else {
        userInterface->constraintPrefsShowConstraints->value(0);
        userInterface->constraintPrefsShowConstraintsText->label(noStr);
	}
    userInterface->constraintPrefsShowConstraintsText->clear_visible_focus();
#endif

    // Initialize preferences on the "Landmarks" tab
	if (globalControl->readBool(ShowLandmarks)) {
        userInterface->landmarkPrefsShowLandmarks->value(1);
        userInterface->landmarkPrefsShowLandmarksText->label(yesStr);
	}
    else {
        userInterface->landmarkPrefsShowLandmarks->value(0);
        userInterface->landmarkPrefsShowLandmarksText->label(noStr);
	}
    userInterface->landmarkPrefsShowLandmarksText->clear_visible_focus();

    userInterface->landmarkPrefsNarrowWidth->value(
        globalControl->readInt(LandmarkNarrowWidth));

    userInterface->landmarkPrefsWideWidth->value(
        globalControl->readInt(LandmarkWideWidth));

    r = (float) globalControl->readDouble(LandmarksColor_R);
    g = (float) globalControl->readDouble(LandmarksColor_G);
    b = (float) globalControl->readDouble(LandmarksColor_B);
    color_code = (Fl_Color) LANDMARKS_COLOR;
    Fl::set_color(color_code, (uchar) (r*255.0f), (uchar) (g*255.0f), (uchar) (b*255.0f));
    userInterface->landmarkPrefsLandmarksColor->color((unsigned char) color_code);
    if (r + g + b < COLOR_CUTOFF)
        userInterface->landmarkPrefsLandmarksColor->labelcolor(FL_WHITE);
    else
        userInterface->landmarkPrefsLandmarksColor->labelcolor(FL_BLACK);

    // Initialize preferences on the "Other" tab
    P3DControl * control;
    control = userInterface->control;
	userInterface->otherPrefsUndoLength->value(control->undoListMaxSize());

#ifdef _WIN32
	char leftHandedMouseComment[] = "(Must be set in the Control Panel on this system).";
	userInterface->otherPrefsLeftComment->value(leftHandedMouseComment);
#endif

	int left_handed;
	left_handed = usingLeftHandedMouse();
#ifdef _WIN32
	left_handed = globalControl->read(LeftHandedMouse, left_handed);
#endif
	userInterface->otherPrefsLeftHandedButton->value(left_handed);
	if (left_handed)
		userInterface->otherPrefsLeftHandedText->label(yesStr);
	else
		userInterface->otherPrefsLeftHandedText->label(noStr);
	userInterface->otherPrefsLeftHandedText->clear_visible_focus();

    userInterface->setPrefsVerbosityBtns(globalControl->readInt(OutputVerbosity));

    userInterface->setPrefsReorderBtns(globalControl->readInt(ReorderModels));

	bool show_about_pablo = globalControl->read(AboutPablo_Open, true);
	userInterface->otherPrefsAboutPabloButton->value(show_about_pablo);
	if (show_about_pablo)
		userInterface->otherPrefsAboutPabloText->label(yesStr);
	else
		userInterface->otherPrefsAboutPabloText->label(noStr);
	userInterface->otherPrefsAboutPabloText->clear_visible_focus();

    r = (float) globalControl->readDouble(BackgroundColor_R);
    g = (float) globalControl->readDouble(BackgroundColor_G);
    b = (float) globalControl->readDouble(BackgroundColor_B);
    color_code = (Fl_Color) BACKGROUND_COLOR;
    Fl::set_color(color_code, (uchar) (r*255.0f), (uchar) (g*255.0f), (uchar) (b*255.0f));
    userInterface->displayPrefsCurrentBkgColor->color((unsigned char) color_code);
    if (r + g + b < COLOR_CUTOFF)
        userInterface->displayPrefsCurrentBkgColor->labelcolor(FL_WHITE);
    else
        userInterface->displayPrefsCurrentBkgColor->labelcolor(FL_BLACK);

	if (globalControl->readBool(SimTransformSeparate)) {
		userInterface->otherPrefsSimTransModelButton->value(0);
		userInterface->otherPrefsSimTransSeparateButton->value(1);
	}
	else {
		userInterface->otherPrefsSimTransModelButton->value(1);
		userInterface->otherPrefsSimTransSeparateButton->value(0);
	}
	if (globalControl->readBool(SimTransformMatrix)) {
		userInterface->otherPrefsSimTransVectorsButton->value(0);
		userInterface->otherPrefsSimTransMatrixButton->value(1);
	}
	else {
		userInterface->otherPrefsSimTransVectorsButton->value(1);
		userInterface->otherPrefsSimTransMatrixButton->value(0);
	}

	// Redraw preference window for display
    userInterface->preferencesEditorWindow->redraw();
    return;
}

void P3DUserInterfaceCallback::setPrefsRememberPositions(int toggle) {
	if (toggle) {
		recordWindowPositions = true;
		userInterface->displayPrefsRememberPosnText->label(yesStr);
	}
	else {
		recordWindowPositions = false;
		userInterface->displayPrefsRememberPosnText->label(noStr);
	}
	userInterface->displayPrefsRememberPosnText->redraw();
}

void P3DUserInterfaceCallback::setPrefsRememberOpen(int toggle) {
	if (toggle) {
		recordWindowsOpen = true;
		userInterface->displayPrefsRememberOpenText->label(yesStr);
	}
	else {
		recordWindowsOpen = false;
		userInterface->displayPrefsRememberOpenText->label(noStr);
	}
	userInterface->displayPrefsRememberOpenText->redraw();
}

void P3DUserInterfaceCallback::setPrefsIconifyWindow(int val)
{
    globalControl->write(IconifyMode, val);

    if(userInterface == NULL)
		return;

	userInterface->setPrefsIconifyBtns(val);
}

void P3DUserInterfaceCallback::setPrefsStandardAxes(int toggle) {
    globalControl->write(DisplayStdAxes, (toggle == 1 ? true : false));

    userInterface->stdAxesCheckBox->value(toggle);
	toggleStdAxes(toggle);

	if (toggle)
		userInterface->displayPrefsStandardAxesText->label(yesStr);
	else
		userInterface->displayPrefsStandardAxesText->label(noStr);
	userInterface->displayPrefsStandardAxesText->redraw();
}

void P3DUserInterfaceCallback::setPrefsSmoothImages(int toggle) {
    P3DControl * control;
    bool smooth;

    smooth = toggle == 1 ? true : false;
    globalControl->write(SmoothImages, smooth);

    if(userInterface == NULL)
		return;

    if (smooth)
        userInterface->imagePrefsSmoothImagesText->label(yesStr);
    else
        userInterface->imagePrefsSmoothImagesText->label(noStr);
    userInterface->imagePrefsSmoothImagesText->redraw();

    control = userInterface->control;
    control->setPrefsSmoothImages(smooth);
    userInterface->modelView->updateImagePlanes(true);	// Force refresh of the image planes
}

void P3DUserInterfaceCallback::setPrefsScaleInputImage(int toggle)
{
    P3DControl * control;
	bool scaleImages;

	scaleImages = toggle == 1 ? true : false;
    globalControl->write(ScaleImages, scaleImages);

    if(userInterface == NULL)
		return;

	if (scaleImages)
		userInterface->imagePrefsScaleInputText->label(yesStr);
	else
		userInterface->imagePrefsScaleInputText->label(noStr);
	userInterface->imagePrefsScaleInputText->redraw();

    control = userInterface->control;
    control->setInputImageScaling(scaleImages);
}

void P3DUserInterfaceCallback::setPrefsCoordsDisplayed(int toggle)
{
	P3DControl * control;

	displayModelCoords = toggle == 1 ? true : false;
    globalControl->write(DisplayModelUnits, displayModelCoords);

    if(userInterface == NULL)
		return;

	if (displayModelCoords)
		userInterface->imagePrefsCoordsDisplayedText->label(modelStr);
	else
		userInterface->imagePrefsCoordsDisplayedText->label(worldStr);
	userInterface->imagePrefsCoordsDisplayedText->redraw();

    control = userInterface->control;
    if(control == NULL)
        return;

	if (control->haveImage()) {
		setXImagePlanePosition();
		setYImagePlanePosition();
		setZImagePlanePosition();
	}
}

void P3DUserInterfaceCallback::setPrefsRangeDisplayed(int toggle)
{
	displayIntensityScale = toggle == 1 ? true : false;
    globalControl->write(DisplayRangeUnits, displayIntensityScale);

    if(userInterface == NULL)
		return;

	if (displayIntensityScale)
		userInterface->imagePrefsRangeDisplayedText->label(screenStr);
	else
		userInterface->imagePrefsRangeDisplayedText->label(percentStr);
	userInterface->imagePrefsRangeDisplayedText->redraw();

	// Trigger updating the image display
	intensityWindowSliderChanged();
}

void P3DUserInterfaceCallback::setPrefsByteOrder(int order)
{
    globalControl->write(ByteOrder, order);
}

void P3DUserInterfaceCallback::setPrefsCompressImages(int toggle) {
    bool compress;

    compress = toggle == 1 ? true : false;
    globalControl->write(CompressImages, compress);

    if (compress) {
        userInterface->imagePrefsCompressImagesText->label(yesStr);
		userInterface->imagePrefsConvertImagesType->label(compressStr);
	}
    else {
        userInterface->imagePrefsCompressImagesText->label(noStr);
		userInterface->imagePrefsConvertImagesType->label(uncompressStr);
	}
	userInterface->prefsEditorImageTab->redraw();
}

void P3DUserInterfaceCallback::setPrefsConvertImages(int toggle)
{
    bool convertCompression;

    convertCompression = toggle == 1 ? true : false;
    globalControl->write(ConvertImages, convertCompression);

    if (convertCompression)
        userInterface->imagePrefsConvertImagesText->label(yesStr);
    else
        userInterface->imagePrefsConvertImagesText->label(noStr);
	userInterface->imagePrefsConvertImagesText->redraw();
}

void P3DUserInterfaceCallback::setPrefsImageFormat(int val)
{
    globalControl->write(ImageFormat, val);
}

void P3DUserInterfaceCallback::setPrefsConvertFormat(int toggle)
{
    bool convert;

    convert = toggle == 1 ? true : false;
    globalControl->write(ConvertImageFormat, convert);

    if (convert)
        userInterface->imagePrefsConvertFormatText->label(yesStr);
    else
        userInterface->imagePrefsConvertFormatText->label(noStr);
	userInterface->prefsEditorImageTab->redraw();
}

void P3DUserInterfaceCallback::setPrefsDrawBoundary(int toggle)
{
	bool drawBoundary;

	drawBoundary = toggle == 1 ? true : false;
    globalControl->write(DrawBoundary, drawBoundary);

    if(userInterface == NULL)
		return;

	toggleDrawBoundary(toggle);
}

void P3DUserInterfaceCallback::setImagePlane(enum slice_t slice, int val)
{
	Fl_Check_Button * checkBox;

	switch (slice) {
		case SAGGITAL:
					checkBox = userInterface->xImagePlaneCheckBox;
					checkBox->value(val);
					viewXImagePlane(val);
					break;
		case CORONAL:
					checkBox = userInterface->yImagePlaneCheckBox;
					checkBox->value(val);
					viewYImagePlane(val);
					break;
		case AXIAL:
					checkBox = userInterface->zImagePlaneCheckBox;
					checkBox->value(val);
					viewZImagePlane(val);
					break;
		default:
					break;
	}
}

void P3DUserInterfaceCallback::setPrefsSliceDefault(int slice, int val)
{
	setImagePlane((enum slice_t) slice, val);
	switch (slice) {
		case SAGGITAL:
					globalControl->write(SagittalSliceDefault, val);
					break;
		case CORONAL:
					globalControl->write(CoronalSliceDefault, val);
					break;
		case AXIAL:
					globalControl->write(AxialSliceDefault, val);
					break;
		default:
					break;
	}
}

void P3DUserInterfaceCallback::setPrefsCutPlaneBoundary(int toggle)
{
    globalControl->write(DrawCutPlaneBoundary, (toggle == 1 ? true : false));

    if (toggle)
        userInterface->imagePrefsCutPlaneBoundaryText->label(yesStr);
    else
        userInterface->imagePrefsCutPlaneBoundaryText->label(noStr);
    userInterface->imagePrefsCutPlaneBoundaryText->redraw();
}

void P3DUserInterfaceCallback::setPrefsCutPlaneBdryWidth(int w)
{
    globalControl->write(CutPlaneBoundaryWidth, w);
	userInterface->modelView->boundaryWidth(w);
    redraw();
}

// Called from P3DUserInterface::make_CutPlaneControl() to make the cut planes
// window reflect the cut planes mode preference, when the window is created.
void P3DUserInterfaceCallback::initCutPlanesMode()
{
    cutPlaneModePreference = globalControl->readInt(CutPlaneMode);
	for (int i = 0; i < 6; i++)
		setCutPlanesMode(i, cutPlaneModePreference);
}

// Called when the cut planes mode preference is changed
void P3DUserInterfaceCallback::incrementPrefsCutPlaneMode()
{
	cutPlaneModePreference++;
	if (cutPlaneModePreference > 2)
		cutPlaneModePreference = 0;

	userInterface->imagePrefsCutPlaneMode->value(0);
	userInterface->imagePrefsCutPlaneMode->label(cutPlaneModeLabels[cutPlaneModePreference]);
	userInterface->imagePrefsCutPlaneMode->redraw();
	if (cutPlaneModePreference == 0)
		userInterface->imagePrefsCutPlaneModeText->label(cutPlaneModeLabel_0);
	else if (cutPlaneModePreference == 2)
		userInterface->imagePrefsCutPlaneModeText->label(cutPlaneModeLabel_1);
	else
		userInterface->imagePrefsCutPlaneModeText->label(cutPlaneModeLabel_2);
	userInterface->imagePrefsCutPlaneModeText->redraw();

    globalControl->write(CutPlaneMode, cutPlaneModePreference);
	for (int i = 0; i < 6; i++)
		setCutPlanesMode(i, cutPlaneModePreference);
}

void P3DUserInterfaceCallback::setPrefsRockingAngle(double angle) {
    globalControl->write(RockingAngle, angle);
    userInterface->modelView->setRockingAngle(angle);
    redraw();
}

void P3DUserInterfaceCallback::setPrefsRockingIncrement(double angle) {
    globalControl->write(RockingIncrement, angle);
    userInterface->modelView->setRockingIncrement(angle);
    redraw();
}

void P3DUserInterfaceCallback::setPrefsAtomVectorsType(int toggle) {
    globalControl->write(AtomVectorsType, (toggle == 1 ? true : false));

    if (toggle)
        userInterface->atomPrefsAtomVectorsTypeText->label(proportionalStr);
    else
        userInterface->atomPrefsAtomVectorsTypeText->label(unitVectorsStr);
    userInterface->atomPrefsAtomVectorsTypeText->redraw();
    userInterface->modelView->setAtomVectorsType((toggle == 1 ? true : false));
    updateObject();
}

void P3DUserInterfaceCallback::setPrefsAtomBVectorsType(int val) {
    globalControl->write(BVectorsType, val);
    userInterface->modelView->setBVectorsType(val);
    updateObject();
}

void P3DUserInterfaceCallback::setPrefsShowAtomVectors(int toggle) {
    globalControl->write(ShowAtomVectors, (toggle == 1 ? true : false));

    if (toggle)
        userInterface->atomPrefsShowAtomVectorsText->label(yesStr);
    else
        userInterface->atomPrefsShowAtomVectorsText->label(noStr);
    userInterface->atomPrefsShowAtomVectorsText->redraw();
    userInterface->atomVectorsVisible->value((toggle == 1 ? true : false));
    userInterface->modelView->setAtomVectorsVisibility((toggle == 1 ? true : false));
    updateObject();
}

void P3DUserInterfaceCallback::setPrefsShowAtoms(int toggle) {
    globalControl->write(ShowAtoms, (toggle == 1 ? true : false));

    if (toggle)
        userInterface->atomPrefsShowAtomsText->label(yesStr);
    else
        userInterface->atomPrefsShowAtomsText->label(noStr);
    userInterface->atomPrefsShowAtomsText->redraw();
    userInterface->atomsVisible->value((toggle == 1 ? true : false));
    userInterface->modelView->setAtomsVisibility((toggle == 1 ? true : false));
    updateObject();
}

void P3DUserInterfaceCallback::setPrefsUseFigureColors(int toggle) {
    globalControl->write(FiguresColorAtoms, (toggle == 1 ? true : false));

    if (toggle)
        userInterface->atomPrefsUseFigureColorsText->label(yesStr);
    else
        userInterface->atomPrefsUseFigureColorsText->label(noStr);
    userInterface->atomPrefsUseFigureColorsText->redraw();
    userInterface->modelView->setFiguresColorAtoms((toggle == 1 ? true : false));
    updateObject();
}

void P3DUserInterfaceCallback::setPrefsExtraAtomVectors(int toggle) {
    globalControl->write(ExtraAtomVectors, (toggle == 1 ? true : false));

    if (toggle)
        userInterface->atomPrefsExtraAtomVectorsText->label(yesStr);
    else
        userInterface->atomPrefsExtraAtomVectorsText->label(noStr);
    userInterface->atomPrefsExtraAtomVectorsText->redraw();
    userInterface->modelView->setShowExtraAtomVectors((toggle == 1 ? true : false));
    updateObject();
}

void P3DUserInterfaceCallback::setPrefsAtomVectorsWidth(int width) {
    globalControl->write(AtomVectorsLineWidth, width);
    userInterface->modelView->setAtomVectorsWidth(width);
    updateObject();
}

void P3DUserInterfaceCallback::setPrefsMeshType(int toggle) {
    globalControl->write(MeshConnectorsType, (toggle == 1 ? true : false));

    if (toggle)
        userInterface->meshPrefsMeshTypeText->label(dottedStr);
    else
        userInterface->meshPrefsMeshTypeText->label(solidStr);
    userInterface->meshPrefsMeshTypeText->redraw();
    userInterface->modelView->setMeshConnectorsType((toggle == 1 ? true : false));
    updateObject();
}

void P3DUserInterfaceCallback::setPrefsMeshConnectorsWidth(int width) {
    globalControl->write(MeshConnectorsLineWidth, width);
    userInterface->modelView->setMeshConnectorsWidth(width);
    updateObject();
}

void P3DUserInterfaceCallback::setPrefsShowMesh(int toggle) {
    globalControl->write(ShowMeshConnectors, (toggle == 1 ? true : false));

    if (toggle)
        userInterface->meshPrefsShowMeshText->label(yesStr);
    else
        userInterface->meshPrefsShowMeshText->label(noStr);
    userInterface->meshPrefsShowMeshText->redraw();
    userInterface->modelView->setShowMeshConnectors((toggle == 1 ? true : false));
    updateObject();
}

void P3DUserInterfaceCallback::setPrefsSurfaceStyle(M3DSurfaceStyle style) {
    globalControl->write(SurfaceStyle, (int) style);
    userInterface->setSurfaceTypeBtns((int) style);
    setSurfaceStyle(style);
}

void P3DUserInterfaceCallback::setPrefsSurfaceLevel(int lvl) {
    globalControl->write(SurfaceLevel, lvl);
    userInterface->setSurfaceLevelBtns(lvl);
    setSurfaceLevel(lvl);
}

void P3DUserInterfaceCallback::setPrefsSurfaceLineWidth(int width)
{
    globalControl->write(SurfaceLineWidth, width);
    setSurfaceLineWidth(width);
}

void P3DUserInterfaceCallback::setPrefsUseTwoLights(int toggle)
{
	bool yesNo = (toggle == 0 ? false : true);
    if (yesNo)
		userInterface->surfacePrefsUseTwoLightsText->label(yesStr);
	else
		userInterface->surfacePrefsUseTwoLightsText->label(noStr);
	userInterface->surfacePrefsUseTwoLightsText->redraw();
    globalControl->write(TwoLights, yesNo);
    userInterface->modelView->useTwoLights(yesNo);
    redraw();
}

void P3DUserInterfaceCallback::setPrefsSurfaceSmoothness(int tolerance) {
    globalControl->write(SurfaceSmoothnessDefault, tolerance);
}

void P3DUserInterfaceCallback::setPrefsPartialRendering(int val) {
    globalControl->write(PartialSurfaceRendering, (val == 1 ? true : false));
    userInterface->setPartialSurfaceRendering((int) val);
    togglePartialRendering(val);
}

void P3DUserInterfaceCallback::setPrefsPartialSurfaceStyle(M3DSurfaceStyle style) {
    globalControl->write(PartialSurfaceStyle, (int) style);
    userInterface->setPartialSurfaceTypeBtns((int) style);
    setPartialSurfaceStyle(style);
}

void P3DUserInterfaceCallback::setPrefsPartialSurfaceLevel(int val) {
    globalControl->write(PartialSurfaceLevel, val);
    userInterface->setPartialSurfaceLevelBtns(val);
    setPartialSurfaceLevel(val);
}

void P3DUserInterfaceCallback::setPrefsShowConstraints(int toggle) {
#ifndef BINARY
    globalControl->write(ShowConstraints, (toggle == 1 ? true : false));

    if (toggle)
        userInterface->constraintPrefsShowConstraintsText->label(yesStr);
    else
        userInterface->constraintPrefsShowConstraintsText->label(noStr);
    userInterface->constraintPrefsShowConstraintsText->redraw();
    userInterface->constraintsVisible->value(toggle);
    userInterface->visibilityControlWindow->redraw();
    constraintsVisible(toggle);
#endif
}

void P3DUserInterfaceCallback::setPrefsShowLandmarks(int toggle)
{
	bool landmarksOn = (toggle == 1 ? true : false);
    globalControl->write(ShowLandmarks, landmarksOn);

    if (landmarksOn)
        userInterface->landmarkPrefsShowLandmarksText->label(yesStr);
    else
        userInterface->landmarkPrefsShowLandmarksText->label(noStr);
    userInterface->landmarkPrefsShowLandmarksText->redraw();
    userInterface->landmarksVisible->value(toggle);
    userInterface->visibilityControlWindow->redraw();
    landmarksVisible(toggle);
	priorLandmarkVisibility = landmarksOn;
}

void P3DUserInterfaceCallback::setPrefsLandmarksNarrowWidth(int width) {
    globalControl->write(LandmarkNarrowWidth, width);
    userInterface->modelView->setLandmarksNarrowWidth(width);
    updateSurface();
}

void P3DUserInterfaceCallback::setPrefsLandmarksWideWidth(int width) {
    globalControl->write(LandmarkWideWidth, width);
    userInterface->modelView->setLandmarksWideWidth(width);
    updateSurface();
}

void P3DUserInterfaceCallback::setPrefsLandmarksColor(enum Fl_Color color_code)
{
    userInterface->landmarkPrefsLandmarksColor->labelcolor(FL_BLACK);
    switch (color_code) {
        case FL_RED:     color[0] = 1.0f;				// Actually orange
                         color[1] = 165.0f/255.0f;
                         color[2] = 0.0f;
						 color_code = (Fl_Color) LANDMARKS_COLOR;
                         Fl::set_color(color_code, 255, 165, 0);
                         break;
        case FL_GREEN:   color[0] = 46.0f/255.0f;		// Actually sea green
                         color[1] = 139.0f/255.0f;
                         color[2] = 87.0f/255.0f;
						 color_code = (Fl_Color) LANDMARKS_COLOR;
                         Fl::set_color(color_code, 46, 139, 87);
                         userInterface->landmarkPrefsLandmarksColor->labelcolor(FL_WHITE);
                         break;
        case FL_BLUE:    color[0] = 100.0f/255.0f;		// Actually cornflower blue
                         color[1] = 149.0f/255.0f;
                         color[2] = 1.0f;
						 color_code = (Fl_Color) LANDMARKS_COLOR;
                         Fl::set_color(color_code, 100, 149, 255);
                         break;
        case FL_YELLOW:  color[0] = 188.0f/255.0f;		// Actually rosy brown
                         color[1] = 143.0f/255.0f;
                         color[2] = 143.0f/255.0f;
						 color_code = (Fl_Color) LANDMARKS_COLOR;
                         Fl::set_color(color_code, 188, 143, 143);
                         break;
        case FL_MAGENTA: color[0] = 160.0f/255.0f;		// Actually purple
                         color[1] = 32.0f/255.0f;
                         color[2] = 240.0f/255.0f;
						 color_code = (Fl_Color) LANDMARKS_COLOR;
                         Fl::set_color(color_code, 160, 32, 240);
                         break;
        case FL_CYAN:    color[0] = 127.0f/255.0f;		// Actually chartreuse
                         color[1] = 1.0f;
                         color[2] = 0.0f;
						 color_code = (Fl_Color) LANDMARKS_COLOR;
                         Fl::set_color(color_code, 127, 255, 0);
                         break;
        case FL_GRAY:    color[0] = 0.6666667f;
                         color[1] = 0.6666667f;
                         color[2] = 0.6666667f;
                         break;
    }
    userInterface->landmarkPrefsLandmarksColor->color((unsigned char) color_code);
    userInterface->landmarkPrefsLandmarksColor->redraw();
    globalControl->write(LandmarksColor_R, color[0]);
    globalControl->write(LandmarksColor_G, color[1]);
    globalControl->write(LandmarksColor_B, color[2]);
    userInterface->modelView->setLandmarksColor(color);
}

void P3DUserInterfaceCallback::setPrefsLandmarksColor()
{
    uchar r, g, b;
    const enum Fl_Color color_code = (Fl_Color) LANDMARKS_COLOR;

    const Fl_Color col = userInterface->landmarkPrefsLandmarksColor->color();
    Fl::get_color(col, r, g, b);

    fl_color_chooser("Landmarks color", r, g, b);
    color[0] = ((float) r)/255.0f;
    color[1] = ((float) g)/255.0f;
    color[2] = ((float) b)/255.0f;
    Fl::set_color(color_code, r, g, b);
    userInterface->landmarkPrefsLandmarksColor->color((unsigned char) color_code);
    if (color[0] + color[1] + color[2] < COLOR_CUTOFF)
        userInterface->landmarkPrefsLandmarksColor->labelcolor(FL_WHITE);
    else
        userInterface->landmarkPrefsLandmarksColor->labelcolor(FL_BLACK);

    userInterface->landmarkPrefsLandmarksColor->redraw();
    globalControl->write(LandmarksColor_R, color[0]);
    globalControl->write(LandmarksColor_G, color[1]);
    globalControl->write(LandmarksColor_B, color[2]);
    userInterface->modelView->setLandmarksColor(color);
}

void P3DUserInterfaceCallback::setPrefsBackgroundColor(enum Fl_Color color_code)
{
	uchar r, g, b;

    userInterface->displayPrefsCurrentBkgColor->labelcolor(FL_WHITE);
    switch (color_code) {
        case FL_BLACK:   color[0] = 0.0f;
                         color[1] = 0.0f;
                         color[2] = 0.0f;
                         break;
        case FL_WHITE:   color[0] = 1.0f;
                         color[1] = 1.0f;
                         color[2] = 1.0f;
                         userInterface->displayPrefsCurrentBkgColor->labelcolor(FL_BLACK);
                         break;
        case FL_DARK3:   color[0] = 0.3333333f;
                         color[1] = 0.3333333f;
                         color[2] = 0.3333333f;
						 r = (uchar) (255.0f*color[0]);
						 g = (uchar) (255.0f*color[1]);
						 b = (uchar) (255.0f*color[2]);
						 Fl::set_color(color_code, r, g, b);
                         break;
        case INDIGO:     color[0] = 0.05f;
                         color[1] = 0.05f;
                         color[2] = 0.1f;
						 r = (uchar) (255.0f*color[0]);
						 g = (uchar) (255.0f*color[1]);
						 b = (uchar) (255.0f*color[2]);
						 Fl::set_color(color_code, r, g, b);
                         break;
    }
    userInterface->displayPrefsCurrentBkgColor->color((unsigned char) color_code);
    userInterface->displayPrefsCurrentBkgColor->redraw();
    userInterface->modelView->setBackgroundColor(color);
    userInterface->modelView->redraw();
    globalControl->write(BackgroundColor_R, color[0]);
    globalControl->write(BackgroundColor_G, color[1]);
    globalControl->write(BackgroundColor_B, color[2]);
}

void P3DUserInterfaceCallback::setPrefsBackgroundColor()
{
    uchar r, g, b;

    const enum Fl_Color color_code = (Fl_Color) BACKGROUND_COLOR;

    const Fl_Color col = userInterface->displayPrefsCurrentBkgColor->color();
    Fl::get_color(col, r, g, b);

    fl_color_chooser("Background color", r, g, b);
    color[0] = ((float) r)/255.0f;
    color[1] = ((float) g)/255.0f;
    color[2] = ((float) b)/255.0f;
    Fl::set_color(color_code, r, g, b);
    userInterface->displayPrefsCurrentBkgColor->color((unsigned char) color_code);
    if (color[0] + color[1] + color[2] < COLOR_CUTOFF)
        userInterface->displayPrefsCurrentBkgColor->labelcolor(FL_WHITE);
    else
        userInterface->displayPrefsCurrentBkgColor->labelcolor(FL_BLACK);
    userInterface->displayPrefsCurrentBkgColor->redraw();
    userInterface->modelView->setBackgroundColor(color);
    userInterface->modelView->redraw();
    globalControl->write(BackgroundColor_R, color[0]);
    globalControl->write(BackgroundColor_G, color[1]);
    globalControl->write(BackgroundColor_B, color[2]);
}

// In this function, the color_code argument implies the color indirectly.
void P3DUserInterfaceCallback::setPrefsPartialSurfaceColor(enum Fl_Color color_code)
{
    userInterface->constraintPrefsPartialSurfaceColor->labelcolor(FL_BLACK);
    switch (color_code) {
        case FL_RED:     partialSurfaceColor[0] = 0.8f;     // Actually pink
                         partialSurfaceColor[1] = 0.545f;
                         partialSurfaceColor[2] = 0.677f;
                         break;
        case FL_GREEN:   partialSurfaceColor[0] = 0.659f;   // Actually tan
                         partialSurfaceColor[1] = 0.547f;
                         partialSurfaceColor[2] = 0.427f;
                         break;
        case FL_BLUE:    partialSurfaceColor[0] = 0.733f;   // Actually light blue
                         partialSurfaceColor[1] = 0.967f;
                         partialSurfaceColor[2] = 0.967f;
                         break;
        case FL_YELLOW:  partialSurfaceColor[0] = 0.768f;   // Actually lime
                         partialSurfaceColor[1] = 0.967f;
                         partialSurfaceColor[2] = 0.287f;
                         break;
        case FL_MAGENTA: partialSurfaceColor[0] = 0.283f;   // Actually olive
                         partialSurfaceColor[1] = 0.396f;
                         partialSurfaceColor[2] = 0.013f;
                         userInterface->constraintPrefsPartialSurfaceColor->labelcolor(FL_WHITE);
                         break;
        case FL_CYAN:    partialSurfaceColor[0] = 0.967f;   // Actually orange
                         partialSurfaceColor[1] = 0.438f;
                         partialSurfaceColor[2] = 0.064f;
                         break;
        case FL_GRAY:    partialSurfaceColor[0] = 0.6666667f;
                         partialSurfaceColor[1] = 0.6666667f;
                         partialSurfaceColor[2] = 0.6666667f;
                         break;
    }
    color_code = (Fl_Color) PARTIAL_SURFACE_COLOR;
    Fl::set_color(color_code, (unsigned char) (partialSurfaceColor[0]*255.0f), (unsigned char)
        (partialSurfaceColor[1]*255.0f), (unsigned char) (partialSurfaceColor[2]*255.0f));
    userInterface->constraintPrefsPartialSurfaceColor->color((unsigned char) color_code);
    userInterface->constraintPrefsPartialSurfaceColor->redraw();
    if (selectedPartialFigures != NULL) {
        selectedPartialFigures->color(partialSurfaceColor);
        P3DControl * control = userInterface->control;
        if (control && control->getPartialDrawingMode())
            updateSurface();
    }
    globalControl->write(PartialSurfaceColor_R, partialSurfaceColor[0]);
    globalControl->write(PartialSurfaceColor_G, partialSurfaceColor[1]);
    globalControl->write(PartialSurfaceColor_B, partialSurfaceColor[2]);
}

void P3DUserInterfaceCallback::setPrefsPartialSurfaceColor()
{
    uchar r, g, b;
    const enum Fl_Color color_code = (Fl_Color) PARTIAL_SURFACE_COLOR;

    const Fl_Color col = userInterface->constraintPrefsPartialSurfaceColor->color();
    Fl::get_color(col, r, g, b);

    fl_color_chooser("Partial surface color", r, g, b);
    partialSurfaceColor[0] = ((float) r)/255.0f;
    partialSurfaceColor[1] = ((float) g)/255.0f;
    partialSurfaceColor[2] = ((float) b)/255.0f;
    Fl::set_color(color_code, r, g, b);
    userInterface->constraintPrefsPartialSurfaceColor->color((unsigned char) color_code);
    if (partialSurfaceColor[0] + partialSurfaceColor[1] + partialSurfaceColor[2] < COLOR_CUTOFF)
        userInterface->constraintPrefsPartialSurfaceColor->labelcolor(FL_WHITE);
    else
        userInterface->constraintPrefsPartialSurfaceColor->labelcolor(FL_BLACK);
    userInterface->constraintPrefsPartialSurfaceColor->redraw();
    if (selectedPartialFigures != NULL) {
        selectedPartialFigures->color(partialSurfaceColor);
        P3DControl * control = userInterface->control;
        if (control && control->getPartialDrawingMode())
            updateSurface();
    }
    globalControl->write(PartialSurfaceColor_R, partialSurfaceColor[0]);
    globalControl->write(PartialSurfaceColor_G, partialSurfaceColor[1]);
    globalControl->write(PartialSurfaceColor_B, partialSurfaceColor[2]);
}

void P3DUserInterfaceCallback::setPrefsMeshColor(enum Fl_Color color_code)
{
    userInterface->meshPrefsCurrentMeshColor->labelcolor(FL_BLACK);
    switch (color_code) {
        case FL_RED:     color[0] = 1.0f;
                         color[1] = 0.0f;
                         color[2] = 0.0f;
                         break;
        case FL_GREEN:   color[0] = 0.0f;
                         color[1] = 1.0f;
                         color[2] = 0.0f;
                         break;
        case FL_BLUE:    color[0] = 0.0f;
                         color[1] = 0.0f;
                         color[2] = 1.0f;
                         userInterface->meshPrefsCurrentMeshColor->labelcolor(FL_WHITE);
                         break;
        case FL_YELLOW:  color[0] = 1.0f;
                         color[1] = 1.0f;
                         color[2] = 0.0f;
                         break;
        case FL_MAGENTA: color[0] = 1.0f;
                         color[1] = 0.0f;
                         color[2] = 1.0f;
                         userInterface->meshPrefsCurrentMeshColor->labelcolor(FL_WHITE);
                         break;
        case FL_CYAN:    color[0] = 0.0f;
                         color[1] = 1.0f;
                         color[2] = 1.0f;
                         break;
        case FL_GRAY:    color[0] = 0.6666667f;
                         color[1] = 0.6666667f;
                         color[2] = 0.6666667f;
                         break;
    }
    userInterface->meshPrefsCurrentMeshColor->color((unsigned char) color_code);
    userInterface->meshPrefsCurrentMeshColor->redraw();
    globalControl->write(ConnectorsColor_R, color[0]);
    globalControl->write(ConnectorsColor_G, color[1]);
    globalControl->write(ConnectorsColor_B, color[2]);
    userInterface->modelView->setMeshConnectorsColor(color);
}

void P3DUserInterfaceCallback::setPrefsMeshColor()
{
    uchar r, g, b;
    const enum Fl_Color color_code = (Fl_Color) CONNECTORS_COLOR;

    const Fl_Color col = userInterface->meshPrefsCurrentMeshColor->color();
    Fl::get_color(col, r, g, b);

    fl_color_chooser("Medial mesh color", r, g, b);
    color[0] = ((float) r)/255.0f;
    color[1] = ((float) g)/255.0f;
    color[2] = ((float) b)/255.0f;
    Fl::set_color(color_code, r, g, b);
    userInterface->meshPrefsCurrentMeshColor->color((unsigned char) color_code);
    if (color[0] + color[1] + color[2] < COLOR_CUTOFF)
        userInterface->meshPrefsCurrentMeshColor->labelcolor(FL_WHITE);
    else
        userInterface->meshPrefsCurrentMeshColor->labelcolor(FL_BLACK);
    userInterface->meshPrefsCurrentMeshColor->redraw();
    globalControl->write(ConnectorsColor_R, color[0]);
    globalControl->write(ConnectorsColor_G, color[1]);
    globalControl->write(ConnectorsColor_B, color[2]);
    userInterface->modelView->setMeshConnectorsColor(color);
}

void P3DUserInterfaceCallback::setUndoLength(int undoLength) {
    globalControl->write(UndoListLength, undoLength);
}

void P3DUserInterfaceCallback::setPrefsMouseButtonOrder(bool left_handed) {
	set_mouse_button_order(left_handed);
    globalControl->write(LeftHandedMouse, left_handed);
	if (left_handed)
		userInterface->otherPrefsLeftHandedText->label(yesStr);
	else
		userInterface->otherPrefsLeftHandedText->label(noStr);
	userInterface->otherPrefsLeftHandedText->redraw();
}

void P3DUserInterfaceCallback::setPrefsReorderWindow(int alwaysNeverPrompt)
{
    globalControl->write(ReorderModels, alwaysNeverPrompt);
    userInterface->setPrefsReorderBtns(alwaysNeverPrompt);
}

void P3DUserInterfaceCallback::setPrefsTilesColor(enum Fl_Color color_code)
{
    userInterface->tilePrefsCurrentTilesColor->labelcolor(FL_BLACK);
    switch (color_code) {
        case FL_RED:     color[0] = 1.0f;			// Actually deep pink
                         color[1] = 0.078431373f;
                         color[2] = 0.576470588f;
 						 color_code = (Fl_Color) TILES_COLOR;
                         Fl::set_color(color_code, 255, 20, 147);
                         break;
        case FL_GREEN:   color[0] = 0.133333333f;	// Actually forest green
                         color[1] = 0.545098039f;
                         color[2] = 0.133333333f;
                         userInterface->tilePrefsCurrentTilesColor->labelcolor(FL_WHITE);
 						 color_code = (Fl_Color) TILES_COLOR;
                         Fl::set_color(color_code, 34, 139, 43);
                         break;
        case FL_BLUE:    color[0] = 0.415686275f;	// Actually slate blue
                         color[1] = 0.352941176f;
                         color[2] = 0.803921569f;
                         userInterface->tilePrefsCurrentTilesColor->labelcolor(FL_WHITE);
 						 color_code = (Fl_Color) TILES_COLOR;
                         Fl::set_color(color_code, 106, 90, 205);
                        break;
        case FL_YELLOW:  color[0] = 1.0f;			// Actually gold
                         color[1] = 0.943137255f;
                         color[2] = 0.0f;
 						 color_code = (Fl_Color) TILES_COLOR;
                         Fl::set_color(color_code, 255, 215, 0);
                         break;
        case FL_MAGENTA: color[0] = 0.647058824f;	// Actually brick
                         color[1] = 0.164705882f;
                         color[2] = 0.164705882f;
                         userInterface->tilePrefsCurrentTilesColor->labelcolor(FL_WHITE);
 						 color_code = (Fl_Color) TILES_COLOR;
                         Fl::set_color(color_code, 165, 42, 42);
                         break;
        case FL_CYAN:    color[0] = 0.498039216f;	// Actually aquamarine
                         color[1] = 1.0f;
                         color[2] = 0.831372549f;
						 color_code = (Fl_Color) TILES_COLOR;
                         Fl::set_color(color_code, 127, 255, 212);
                         break;
        case FL_GRAY:    color[0] = 0.960784314f;	// Actually beige
                         color[1] = 0.960784314f;
                         color[2] = 0.862745098f;
						 color_code = (Fl_Color) TILES_COLOR;
                         Fl::set_color(color_code, 245, 245, 220);
                         break;
    }
    userInterface->tilePrefsCurrentTilesColor->color((unsigned char) color_code);
    userInterface->tilePrefsCurrentTilesColor->redraw();
    globalControl->write(TilesColor_R, color[0]);
    globalControl->write(TilesColor_G, color[1]);
    globalControl->write(TilesColor_B, color[2]);
    userInterface->modelView->setTilesColor(color);
}

void P3DUserInterfaceCallback::setPrefsTilesColor()
{
    uchar r, g, b;
    const enum Fl_Color color_code = (Fl_Color) TILES_COLOR;

    const Fl_Color col = userInterface->tilePrefsCurrentTilesColor->color();
    Fl::get_color(col, r, g, b);

    fl_color_chooser("Tiled surfaces color", r, g, b);
    color[0] = ((float) r)/255.0f;
    color[1] = ((float) g)/255.0f;
    color[2] = ((float) b)/255.0f;
    Fl::set_color(color_code, r, g, b);
    userInterface->tilePrefsCurrentTilesColor->color((unsigned char) color_code);
    if (color[0] + color[1] + color[2] < COLOR_CUTOFF)
        userInterface->tilePrefsCurrentTilesColor->labelcolor(FL_WHITE);
    else
        userInterface->tilePrefsCurrentTilesColor->labelcolor(FL_BLACK);
    userInterface->tilePrefsCurrentTilesColor->redraw();
    globalControl->write(TilesColor_R, color[0]);
    globalControl->write(TilesColor_G, color[1]);
    globalControl->write(TilesColor_B, color[2]);
    userInterface->modelView->setTilesColor(color);
}

void P3DUserInterfaceCallback::setPrefsBYUOutputType(int type)
{
	globalControl->write(BYUOutputType, type);
#ifdef BINARY
    switch(type) {
	    case 0:	userInterface->otherPrefsBYUTrianglesButton->value(1);
			    userInterface->otherPrefsBYUQuadsButton->value(0);
			    break;
	    case 1:	userInterface->otherPrefsBYUTrianglesButton->value(0);
			    userInterface->otherPrefsBYUQuadsButton->value(1);
			    break;
    }
#else
    switch(type) {
	    case 0:	userInterface->otherPrefsBYUTrianglesButton->value(1);
			    userInterface->otherPrefsBYUQuadsButton->value(0);
				globalControl->write(BYUOutputType, type);
			    break;
	    case 1:	userInterface->otherPrefsBYUTrianglesButton->value(0);
			    userInterface->otherPrefsBYUQuadsButton->value(1);
				globalControl->write(BYUOutputType, type);
			    break;
	    case 2:	userInterface->otherPrefsBYUModelCoordButton->value(1);
			    userInterface->otherPrefsBYUWorldCoordButton->value(0);
				globalControl->write(BYUOutputCoords, type - 2);
			    break;
	    case 3:	userInterface->otherPrefsBYUModelCoordButton->value(0);
			    userInterface->otherPrefsBYUWorldCoordButton->value(1);
				globalControl->write(BYUOutputCoords, type - 2);
			    break;
    }
#endif
}

void P3DUserInterfaceCallback::setPrefsSimTransOutputType()
{
	if (userInterface->otherPrefsSimTransModelButton->value() == 1)
		globalControl->write(SimTransformSeparate, false);
	else
		globalControl->write(SimTransformSeparate, true);

	if (userInterface->otherPrefsSimTransVectorsButton->value() == 1)
		globalControl->write(SimTransformMatrix, false);
	else
		globalControl->write(SimTransformMatrix, true);
}

void P3DUserInterfaceCallback::setPrefsVerbosity(int level)
{
    globalControl->write(OutputVerbosity, level);
	globalVerbosity = level;
}

void P3DUserInterfaceCallback::setPrefsAboutWindow(bool show)
{
	globalControl->write(AboutPablo_Open, show);
	if (show) {
	    userInterface->otherPrefsAboutPabloButton->value(1);
		userInterface->otherPrefsAboutPabloText->label(yesStr);
	}
	else {
	    userInterface->otherPrefsAboutPabloButton->value(0);
		userInterface->otherPrefsAboutPabloText->label(noStr);
	}
	userInterface->otherPrefsAboutPabloText->redraw();
}

// ----------------------------- Start of Optimizer Methods --------------------------------

#ifndef BINARY

void P3DUserInterfaceCallback::setOptimizerIntensityWindowing()
{
    P3DControl * control;
    Fl_Button * button;

    control = userInterface->control;
    if(control == NULL)
        return;

	button = userInterface->optimizerWindowedIntensities;
    control->setOptimizerIntensityWindowing(button->value() ? true : false);
}

#ifdef PRODUCTION_VERSION

void P3DUserInterfaceCallback::previewOptimizationStage()
{
}

#else

// Called when the PROBE button in the optimization interface is pressed.
void P3DUserInterfaceCallback::previewOptimizationStage()
{
    P3DControl * control;
    OptimizerStage currentOptimizerStage;
	Fl_Button * prevBtn;

    if(userInterface == NULL)
        return;

    control = userInterface->control;
    if(control == NULL)
        return;

    if(optimizerIsRunning)
        return;

    currentOptimizerStage = control->getCurrentOptimizerStage();

	previewMode = true;
	// Keep the button down, because the atom stage takes some time
	prevBtn = userInterface->previewOptimizationBtn;
	prevBtn->value(1);

    switch(currentOptimizerStage)
    {
        case OPTIMIZER_OBJECT_STAGE:
			startModelStage(false);
            break;

        case OPTIMIZER_ENSEMBLE_STAGE:
			startModelStage(true);
            break;

        case OPTIMIZER_FIGURE_STAGE:
			startMainFigureStage(false);
            break;

        case OPTIMIZER_MAIN_FIGURE_STAGE:
			startMainFigureStage(true);
            break;

        case OPTIMIZER_SUBFIGURE_STAGE:
			startSubfigureStage();
            break;

        case OPTIMIZER_PRIMITIVE_STAGE:
			startPrimitiveStage();
            break;

        case OPTIMIZER_BOUNDARY_DISPLACEMENT_STAGE:
			startVoxelStage();
            break;

        case OPTIMIZER_END_STAGE:
        case OPTIMIZER_SETUP_STAGE:
        case NO_OPTIMIZER_STAGE:
        default:
            break;
    };

	// Done, so release the button
	prevBtn->value(0);
	previewMode = false;
}

#endif	/* PRODUCTION_VERSION */

void P3DUserInterfaceCallback::zeroValuators(OptimizerStage s)
{
    switch(s)
    {
        case OPTIMIZER_OBJECT_STAGE:
			userInterface->objectMatchValueOutput->value(0);
			userInterface->objectPenaltyValueOutput->value(0);
            break;

        case OPTIMIZER_ENSEMBLE_STAGE:
			userInterface->modelMatchValueOutput->value(0);
			userInterface->modelPenaltyValueOutput->value(0);
            break;

        case OPTIMIZER_FIGURE_STAGE:
			userInterface->figureStageImageMatchValueOutput->value(0);
			userInterface->figureStagePenaltyValueOutput->value(0);
			userInterface->figureStageObjectiveFunctionValueOutput->value(0);
            break;

        case OPTIMIZER_MAIN_FIGURE_STAGE:
			userInterface->mainFigureImageMatchValueOutput->value(0);
			userInterface->mainFigurePenaltyValueOutput->value(0);
			userInterface->mainFigureObjectiveFunctionValueOutput->value(0);
            break;

        case OPTIMIZER_SUBFIGURE_STAGE:
			userInterface->subfigureMatchValueOutput->value(0);
            break;

        case OPTIMIZER_PRIMITIVE_STAGE:
			userInterface->primitiveObjectiveFunctionValueOutput->value(0);
			userInterface->primativePenaltyValueOutput->value(0);
			userInterface->primitiveImageMatchValueOutput->value(0);
            break;

        case OPTIMIZER_BOUNDARY_DISPLACEMENT_STAGE:
			userInterface->voxelStageMatchValueOutput->value(0);
			userInterface->voxelStagePenaltyValueOutput->value(0);
			userInterface->voxelStageObjectiveFunctionValueOutput->value(0);
            break;

        case OPTIMIZER_END_STAGE:
        case OPTIMIZER_SETUP_STAGE:
        case NO_OPTIMIZER_STAGE:
        default:
            break;
    };
}

// Called when the NEXT button in the optimization interface is pressed.
// Also called with true from updateOptimizerWindow() under script control.
// The new stage to be run is returned.
OptimizerStage P3DUserInterfaceCallback::nextOptimizationStage(bool interactive)
{
    P3DControl * control;
    OptimizerStage currentOptimizerStage;
    Fl_Group * stage;
    M3DObject * object;
	Fl_Button * nextBtn;
	Fl_Button * prevBtn;
	static bool switchingStages = false;

    if(userInterface == NULL)
        return NO_OPTIMIZER_STAGE;

    control = userInterface->control;
    if(control == NULL)
        return NO_OPTIMIZER_STAGE;

    if(optimizerIsRunning)
        return NO_OPTIMIZER_STAGE;

	// Make certain we have an object
    object = control->getObjectPtr();
    if (object == NULL) {
        cerr << "A model must first be loaded" << endl;
		fl_beep(FL_BEEP_ERROR);
        return NO_OPTIMIZER_STAGE;
    }

	// Make certain we have an image
    if (control->haveImage() == false) {
        cerr << "An image must first be loaded" << endl;
		fl_beep(FL_BEEP_ERROR);
        return NO_OPTIMIZER_STAGE;
    }

	disableEditingWindows();

	// Prohibit reentry; not doing this can result in corruption of
	// the optimizer window
	if (switchingStages)
		return NO_OPTIMIZER_STAGE;
	switchingStages = true;

	if (interactive)
		currentOptimizerStage = currentInteractiveStage;
	else
		currentOptimizerStage = control->getCurrentOptimizerStage();
	//cout << "Leaving stage " << currentOptimizerStage << endl;

	prevBtn = userInterface->previewOptimizationBtn;

    // First handle the stage we are leaving
    switch(currentOptimizerStage)
    {
        case NO_OPTIMIZER_STAGE:
        case OPTIMIZER_SETUP_STAGE:
			// Deactivate the NEXT button because initialization is slow
		    nextBtn = userInterface->nextOptimizationBtn;
			nextBtn->deactivate();

			if (! interactive)
				if (! initializeOptimizer()) {
					// Initialization failed: exit use of the optimization window
					switchingStages = false;
					nextBtn->activate();	// Restore it
					enableEditingWindows();
					return NO_OPTIMIZER_STAGE;
				}
            resettingOptimizer = true;
			optimizerWindowInUse = true;
            stage = userInterface->optimizerSetupStage;
			nextBtn->activate();	// Restore it
#ifndef PRODUCTION_VERSION
			if (! interactive) {
				prevBtn->show();
				prevBtn->activate();
			}
			previewMode = false;
#endif
			lastModeWasPreview = false;
            break;

        case OPTIMIZER_OBJECT_STAGE:
            clearModelStage(false);
//            control->stopRecordingTransformation();
            resettingOptimizer = true;
            stage = userInterface->optimizerObjectStage;
            break;

        case OPTIMIZER_ENSEMBLE_STAGE:
            clearModelStage(true);
//            control->stopRecordingTransformation();
            resettingOptimizer = true;
            stage = userInterface->optimizerEnsembleStage;
            break;

        case OPTIMIZER_FIGURE_STAGE:
            clearMainFigureStage(false);
            resettingOptimizer = true;
            stage = userInterface->optimizerFigureStage;
            break;

        case OPTIMIZER_MAIN_FIGURE_STAGE:
            clearMainFigureStage(true);
            resettingOptimizer = true;
            stage = userInterface->optimizerMainFigureStage;
            break;

        case OPTIMIZER_SUBFIGURE_STAGE:
            clearSubfigureStage();
            resettingOptimizer = true;
            stage = userInterface->optimizerSubfigureStage;
            break;

        case OPTIMIZER_PRIMITIVE_STAGE:
            clearPrimitiveStage();
            resettingOptimizer = true;
            stage = userInterface->optimizerPrimitiveStage;
            break;

        case OPTIMIZER_BOUNDARY_DISPLACEMENT_STAGE:
            clearVoxelStage();
            resettingOptimizer = true;
            stage = userInterface->optimizerVoxelStage;
            break;

        case OPTIMIZER_END_STAGE:
            stage = userInterface->optimizerEndStage;
			break;

        default:
            stage = NULL;
            resettingOptimizer = false;
            break;
    };

    if(stage != NULL)
        stage->hide();

#ifndef PRODUCTION_VERSION
	if (! interactive)
		prevBtn->activate();
#endif

    // Get the next stage we are entering (depends on model)
	if (interactive)
		currentOptimizerStage = nextInteractiveStage;
	else
		currentOptimizerStage = control->nextOptimizerStage();
	//cout << "Entering stage " << currentOptimizerStage << endl;

    // Finally initialize and display user interface for next stage
    switch(currentOptimizerStage)
    {
        case NO_OPTIMIZER_STAGE:
			if (userInterface->pgaDialog->shown())
				disableOptimizerWindow();
			// NEXT button was pressed from OPTIMIZER_END_STAGE, so fall through
        case OPTIMIZER_SETUP_STAGE:
            stage = userInterface->optimizerSetupStage;
            break;

        case OPTIMIZER_OBJECT_STAGE:
            stage = userInterface->optimizerObjectStage;
            break;

        case OPTIMIZER_ENSEMBLE_STAGE:
            stage = userInterface->optimizerEnsembleStage;
            break;

        case OPTIMIZER_FIGURE_STAGE:
			setOptimizationFigureName(control->getOptimizationFigureNumber());
			stage = userInterface->optimizerFigureStage;

			// Deactivate various sliders when doing binary image match
			if (control->type() == BINARY_IMAGE_MATCH) {
				userInterface->figureStagePenaltyWeightSlider->deactivate();
				userInterface->primitivePenaltyWeightSlider->deactivate();
				userInterface->primitiveConstraintsPenaltyWeightSlider->deactivate();
				userInterface->primitiveNeighborPenaltyWeightSlider->deactivate();
			}
			else {
				userInterface->figureStagePenaltyWeightSlider->activate();
				userInterface->primitivePenaltyWeightSlider->activate();
				userInterface->primitiveConstraintsPenaltyWeightSlider->activate();
				userInterface->primitiveNeighborPenaltyWeightSlider->activate();
			}

			// Deactivate constraints slider when there is only one figure tree
			// Also disabled when doing binary image match, which is done on only 1 figure
			if (object != NULL && object->getFigureTreeCount() < 2
				&& optimizationType != HISTO_IMAGE_MATCH)
					userInterface->figureStageConstraintsPenaltyWeightSlider->deactivate();
			else
				userInterface->figureStageConstraintsPenaltyWeightSlider->activate();
            break;

        case OPTIMIZER_MAIN_FIGURE_STAGE:
			setOptimizationFigureName(control->getOptimizationFigureNumber());
			stage = userInterface->optimizerMainFigureStage;
            break;

        case OPTIMIZER_SUBFIGURE_STAGE:
            setOptimizationFigureName(control->getOptimizationFigureNumber());
            stage = userInterface->optimizerSubfigureStage;
            break;

        case OPTIMIZER_PRIMITIVE_STAGE:
            setOptimizationFigureName(control->getOptimizationFigureNumber());
            stage = userInterface->optimizerPrimitiveStage;

			// Deactivate constraints slider when there is only one figure tree
			if (object != NULL && object->getFigureTreeCount() < 2 ||
				control->type() == BINARY_IMAGE_MATCH)
					userInterface->primitiveConstraintsPenaltyWeightSlider->deactivate();
			else
				userInterface->primitiveConstraintsPenaltyWeightSlider->activate();
            break;

        case OPTIMIZER_BOUNDARY_DISPLACEMENT_STAGE:
            setOptimizationFigureName(control->getOptimizationFigureNumber());
            stage = userInterface->optimizerVoxelStage;

			// Deactivate constraints slider when there is only one figure tree
			if (object != NULL && object->getFigureTreeCount() < 2)
				userInterface->voxelStageConstraintsPenaltyWeightSlider->deactivate();
			else
				userInterface->voxelStageConstraintsPenaltyWeightSlider->activate();
            break;

        case OPTIMIZER_END_STAGE:
            stage = userInterface->optimizerEndStage;
#ifndef PRODUCTION_VERSION
			prevBtn->hide();
#endif
			optimizerWindowInUse = false;
			enableDeformationWindow();
			break;

        default:
            stage = NULL;
            break;
    };

	// Zero the valuators of the next stage
	zeroValuators(currentOptimizerStage);

    if(stage != NULL)
        stage->show();

    update();
	switchingStages = false;
	return currentOptimizerStage;
}

// Called when either the Cancel button of the optimizer window
// is pressed, or when that window is dismissed.
void P3DUserInterfaceCallback::cancelOptimization()
{
    P3DControl * control;
    OptimizerStage currentOptimizerStage;
    Fl_Group * stage;
#ifndef PRODUCTION_VERSION
	Fl_Button * prevBtn;
#endif

    if (userInterface == NULL)
        return;

    control = userInterface->control;
    if (control == NULL)
        return;

	if (interactiveOpt) {
		// A script is running.  Either the Cancel Button
		// was pressed or the window was closed.
        userInterface->optimizerControlDialog->hide();
		return;
	}
	// The window may have been dismissed while optimizing
    optimizerStopRequested = true;

#ifdef OPTIMIZATION_VISUALIZER
	if (ovUI)
		((OptVisualizerUI *) ovUI)->hide();
#endif

    if (userInterface->optimizerControlDialog != NULL)
        userInterface->optimizerControlDialog->hide();

    currentOptimizerStage = control->getCurrentOptimizerStage();

    // First handle the stage we are leaving
    switch(currentOptimizerStage)
    {
        case OPTIMIZER_SETUP_STAGE:
            stage = userInterface->optimizerSetupStage;
            break;

        case OPTIMIZER_OBJECT_STAGE:
            clearModelStage(false);
            stage = userInterface->optimizerObjectStage;
            break;

        case OPTIMIZER_ENSEMBLE_STAGE:
            clearModelStage(true);
            stage = userInterface->optimizerEnsembleStage;
            break;

        case OPTIMIZER_FIGURE_STAGE:
            clearMainFigureStage(false);
            stage = userInterface->optimizerFigureStage;
            break;

        case OPTIMIZER_MAIN_FIGURE_STAGE:
            clearMainFigureStage(true);
            stage = userInterface->optimizerMainFigureStage;
            break;

        case OPTIMIZER_SUBFIGURE_STAGE:
            clearSubfigureStage();
            stage = userInterface->optimizerSubfigureStage;
            break;

        case OPTIMIZER_PRIMITIVE_STAGE:
            clearPrimitiveStage();
            stage = userInterface->optimizerPrimitiveStage;
            break;

		case OPTIMIZER_BOUNDARY_DISPLACEMENT_STAGE:
            clearVoxelStage();
            stage = userInterface->optimizerVoxelStage;
            break;

        case OPTIMIZER_END_STAGE:
            stage = userInterface->optimizerEndStage;
            break;

        case NO_OPTIMIZER_STAGE:
        default:
            stage = NULL;
            break;
    };

    if(stage != NULL)
        stage->hide();

    stage = userInterface->optimizerSetupStage;
    if(stage != NULL)
        stage->show();

#ifndef PRODUCTION_VERSION
	prevBtn = userInterface->previewOptimizationBtn;
	prevBtn->hide();
#endif

	enableEditingWindows();

    control->cancelOptimization();
    control->setEditingMode(ANY_EDITING_MODE);
    control->allowSelection(true);

    resettingOptimizer = false;
    optimizerIsRunning = false;
	optimizerWindowInUse = false;
	enableDeformationWindow();
}

bool P3DUserInterfaceCallback::initializeOptimizer()
{
    P3DControl * control;
    const char * filename;
	bool ret = false;
	extern P3DUserInterface * ui;

	optimizationType = NO_MATCH_TYPE;
	profileFilename = NULL;
	int surfaceLevel = userInterface->modelView->getSurfaceLevel();

    control = userInterface->control;
	M3DObject * object = control->getObjectPtr();
	if (! object->verifyInBounds())		// AGG: This is ignoring the possibility of a loaded symmetry transform
		cout << "Warning: the model is larger than or outside the unit cube" << endl;

    if (userInterface->gaussDerivMatchRadioBtn->value() == 1)
        ret = control->initializeMatch(GAUSSIAN_DERIVATIVE_MATCH, surfaceLevel, NULL);
	else if (userInterface->negGaussDerivMatchRadioBtn->value() == 1)
        ret = control->initializeMatch(NEG_GAUSSIAN_DERIVATIVE_MATCH, surfaceLevel, NULL);
	else if (userInterface->simpleMaskMatchRadioBtn->value() == 1) {
		filename = askSingleFilename(ModelDirectory,
			"Load Simple Mask File", "*.smsk");
		if (filename != NULL) {
			profileFilename = new char[1 + strlen(filename)];
			strcpy(profileFilename, filename);
		}
		else {
			if (profileFilename != NULL) {
				delete [] profileFilename;
				profileFilename = NULL;
			}
			cout << "A simple mask file must be loaded" << endl;
			fl_beep(FL_BEEP_ERROR);
			return false;
		}
        ret = control->initializeMatch(SIMPLE_MASK_MATCH, surfaceLevel, NULL,
			profileFilename);
		if (! ret) {
			delete [] profileFilename;
			profileFilename = NULL;
		}
	}
    else if (userInterface->binaryMatchRadioBtn->value() == 1)
        ret = control->initializeMatch(BINARY_IMAGE_MATCH, surfaceLevel);
    else if (userInterface->histogramMatchRadioBtn->value() == 1) {
		/* Histogram match.  This code should work, in principle.  It uses
		   the geometry weight slider to specify the histogram match penalty
		   and the inter-object weight slider to specify Mahalanobis penalty.
		   The D^2, constraint, and profile (mask) penalties are not used. 
		   All other tuning settings are defaults, which may be a problem.
		 */
		filename = askSingleFilename(ModelDirectory, "Load Surface Patch Ensemble File", "*.spe");
		Registry dummy;
		if (! control->loadSPEObject(filename, dummy))
		{
			cout << "A surface patch ensemble file must be loaded" << endl;
			fl_beep(FL_BEEP_ERROR);
			return false;
		}
        ret = control->initializeMatch(HISTO_IMAGE_MATCH, surfaceLevel);
	}
	if (ret) {
		optimizerStarted = true;
		optimizationType = control->type();
	}
	else
		cancelOptimization();
//    control->recordTransformation();
	return ret;
}

void P3DUserInterfaceCallback::setOptimizationFigureName(int figureId) {
    P3DControl * control;
    M3DObject * object;
	char number[3];

    if(userInterface == NULL)
        return;

    control = userInterface->control;
    if(control == NULL)
        return;

    object = control->getObjectPtr();
    if(object == NULL)
        return;

    M3DFigure * figure = object->getFigurePtr(figureId);

	const char * name = figure->getName();
	if (name == NULL || 0 == strlen(name)) {
		(void) sprintf(number, "%d", figureId);
		name = number;
	}

	switch (control->getCurrentOptimizerStage()) {
	    case OPTIMIZER_SETUP_STAGE:
	    case OPTIMIZER_OBJECT_STAGE:
	    case OPTIMIZER_ENSEMBLE_STAGE:
	    case OPTIMIZER_END_STAGE:
				// No action necessary
				break;

        case OPTIMIZER_FIGURE_STAGE:
				userInterface->figureStageFigureName->value(name);
				userInterface->figureStageFigureName->redraw();	// AGG: needed?
				break;

	    case OPTIMIZER_MAIN_FIGURE_STAGE:
				userInterface->mainFigureName->value(name);
				userInterface->mainFigureName->redraw();	// AGG: needed?
				break;

	    case OPTIMIZER_SUBFIGURE_STAGE:
				userInterface->subFigureName->value(name);
				userInterface->subFigureName->redraw();
				break;

	    case OPTIMIZER_PRIMITIVE_STAGE:
				userInterface->primitiveFigureName->value(name);
				userInterface->primitiveFigureName->redraw();
				break;

	    case OPTIMIZER_BOUNDARY_DISPLACEMENT_STAGE:
				userInterface->voxelStageFigureName->value(name);
				userInterface->voxelStageFigureName->redraw();
				break;
	}
	Fl::check();	// AGG: needed?
}

// Used for the initial fit of the full ensemble (model) or a multifigure tree (object)
void P3DUserInterfaceCallback::startModelStage(bool modelStage)
{
    P3DControl * control;
    Fl_Value_Output * iterationOutput;
    Fl_Value_Output * imageMatchValueOutput;
    Fl_Value_Output * geometricPenaltyValueOutput;
    Fl_Button * startBtn;
    Fl_Button * resetBtn;
    Fl_Button * nextBtn;
    Fl_Button * cancelBtn;
	Fl_Button * prevBtn;

    int i;
    int startIteration;

#ifdef OUTPUT_DEFORMING_MODELS
    M3DObjectFile objectFile;
    M3DObject * object;
    char filename[1024];
    char outputFilename[1024];
    char * modelFilename;
#endif


    if (optimizerIsRunning)
        return;

    if (userInterface == NULL)
        return;

    control = userInterface->control;
    if (control == NULL)
        return;

	if (modelStage) {	//JJ : Ensemble stage
		startBtn = userInterface->modelStartBtn;
		resetBtn = userInterface->modelResetBtn;
		iterationOutput = userInterface->modelIterationOutput;
		imageMatchValueOutput = userInterface->modelMatchValueOutput;
		geometricPenaltyValueOutput = userInterface->modelPenaltyValueOutput;
	}
	else {   //JJ : Multi figure stage
		startBtn = userInterface->objectStartBtn;
		resetBtn = userInterface->objectResetBtn;
		iterationOutput = userInterface->objectIterationOutput;
		imageMatchValueOutput = userInterface->objectMatchValueOutput;
		geometricPenaltyValueOutput = userInterface->objectPenaltyValueOutput;
	}
	if (iterationOutput == NULL || imageMatchValueOutput == NULL
		|| geometricPenaltyValueOutput == NULL)
			return;

	nextBtn = userInterface->nextOptimizationBtn;
	cancelBtn = userInterface->optimizationCancelBtn;
	prevBtn = userInterface->previewOptimizationBtn;
    if (startBtn == NULL || resetBtn == NULL || nextBtn == NULL || cancelBtn == NULL)
		return;

    optimizerIsRunning = true;

    // Deactivate (grey) buttons we don't want to get hit
    startBtn->deactivate();
    resetBtn->deactivate();
    nextBtn->deactivate();
    cancelBtn->deactivate();

    Fl::check();

    if (! previewMode)
        startIteration = (int) iterationOutput->value() + 1;
    else
        startIteration = 1;

    control->setEditingMode(NO_EDITING_MODE);
    setModelStagePenaltyWeight(modelStage);
	if (! previewMode)
		if (modelStage)
			cout << "Object ensemble" << flush;
		else
			cout << "Object" << flush;

	if (control->testAndResetModelChanged())
		resettingOptimizer = true;

	if (! previewMode && lastModeWasPreview)
		resettingOptimizer = true;

    if(resettingOptimizer)
    {
		if (! optimizerStarted) {
			int surfaceLevel = userInterface->modelView->getSurfaceLevel();
			(void) control->initializeMatch(GAUSSIAN_DERIVATIVE_MATCH, surfaceLevel, NULL);
		}
		else if (previewMode)
			optimizerStarted = false;

//        control->stopRecordingTransformation();
		control->initializeRegistrationOptimizer();

        resettingOptimizer = false;
    }

    optimizerStopRequested = false;

	int niter = MAX_OPTIMIZER_ITERATIONS;
	if (previewMode)
		niter = 2;
	else
		prevBtn->deactivate();	// Once START has been pressed, previewing is disallowed

    for(i = startIteration; i < niter; i++)
    {
		if (! previewMode) {
			control->doRegistrationIterations(1);
			iterationOutput->value(i);
		}
        imageMatchValueOutput->value(-control->getLastRegistrationMatchValue());
        geometricPenaltyValueOutput->value(control->getLastRegistrationPenaltyValue());

#ifdef OUTPUT_DEFORMING_MODELS
        // Output the current model to file
        modelFilename = control->getModelFileName();
        if(modelFilename == NULL)
            strcpy(filename, "registration");
        else
        {
            strcpy(filename, modelFilename);
            delete [] modelFilename;
        }

        sprintf(outputFilename, "%s/registration_%.4d.m3d", filename, i);

        object = control->getObjectPtr();
        if(object != NULL)
            objectFile.write(outputFilename, *object);
#endif
        update();
        Fl::check();
        if(optimizerStopRequested == true)
            break;

        if(control->registrationOptimizerIsFinished())
            break;
    }

    optimizerIsRunning = false;
	lastModeWasPreview = previewMode;
    updateVisibilityFiguresDisplay(false);
    control->setEditingMode(ANY_EDITING_MODE);

    // Reactivate buttons that we can use now
    startBtn->activate();
    resetBtn->activate();
    nextBtn->activate();
    cancelBtn->activate();

	if (! previewMode)
		cout << " stage" << endl;
}

void P3DUserInterfaceCallback::stopModelStage()
{
    optimizerStopRequested = true;
}

void P3DUserInterfaceCallback::resetModelStage(bool modelStage)
{
    P3DControl * control;
    Fl_Value_Output * iterationOutput;
    Fl_Value_Output * imageMatchValueOutput;
    Fl_Value_Output * geometricPenaltyValueOutput;

    if(optimizerIsRunning)
        return;

    if(userInterface == NULL)
        return;

    control = userInterface->control;
    if(control == NULL)
        return;

	if (modelStage) {
		iterationOutput = userInterface->modelIterationOutput;
		imageMatchValueOutput = userInterface->modelMatchValueOutput;
		geometricPenaltyValueOutput = userInterface->modelPenaltyValueOutput;
	}
	else {
		iterationOutput = userInterface->objectIterationOutput;
		imageMatchValueOutput = userInterface->objectMatchValueOutput;
		geometricPenaltyValueOutput = userInterface->objectPenaltyValueOutput;
	}

//    control->recordTransformation();
    resettingOptimizer = true;

    if(iterationOutput != NULL)
        iterationOutput->value(0);
    if(imageMatchValueOutput != NULL)
        imageMatchValueOutput->value(-control->getLastRegistrationMatchValue());
    if(geometricPenaltyValueOutput != NULL)
        geometricPenaltyValueOutput->value(0);
}

void P3DUserInterfaceCallback::setModelStagePenaltyWeight(bool modelStage)
{
    P3DControl * control;
    Fl_Value_Slider * penaltyWeightSlider;
	double wt;

    control = userInterface->control;
    if (control == NULL)
        return;

	if (modelStage)
		penaltyWeightSlider = userInterface->modelPenaltyWeightSlider;
	else
		penaltyWeightSlider = userInterface->objectPenaltyWeightSlider;

    if (penaltyWeightSlider == NULL)
        return;

    wt = penaltyWeightSlider->value();
    control->setRegistrationPenaltyWeight(wt);
	if (globalVerbosity > 0)
		cout << "Penalty weight = " << wt << '\n';
}

void P3DUserInterfaceCallback::clearModelStage(bool modelStage)
{
    Fl_Value_Output * iterationOutput;
    Fl_Value_Output * imageMatchValueOutput;
    Fl_Value_Output * geometricPenaltyValueOutput;

    Fl_Button * startBtn,
              * resetBtn,
              * nextBtn,
              * cancelBtn;

    if(userInterface == NULL)
        return;

    if(optimizerIsRunning)
        return;

	if (modelStage) {
		startBtn = userInterface->modelStartBtn;
		resetBtn = userInterface->modelResetBtn;
		iterationOutput = userInterface->modelIterationOutput;
		imageMatchValueOutput = userInterface->modelMatchValueOutput;
		geometricPenaltyValueOutput = userInterface->modelPenaltyValueOutput;
	}
	else {
		startBtn = userInterface->objectStartBtn;
		resetBtn = userInterface->objectResetBtn;
		iterationOutput = userInterface->objectIterationOutput;
		imageMatchValueOutput = userInterface->objectMatchValueOutput;
		geometricPenaltyValueOutput = userInterface->objectPenaltyValueOutput;
	}
    nextBtn = userInterface->nextOptimizationBtn;
    cancelBtn = userInterface->optimizationCancelBtn;

    if(iterationOutput != NULL)
        iterationOutput->value(0);
    if(imageMatchValueOutput != NULL)
        imageMatchValueOutput->value(0);
    if(geometricPenaltyValueOutput != NULL)
        geometricPenaltyValueOutput->value(0);

    if(startBtn != NULL)
        startBtn->activate();
    if(resetBtn != NULL)
        resetBtn->activate();
    if(nextBtn != NULL)
        nextBtn->activate();
    if(cancelBtn != NULL)
        cancelBtn->activate();
}

// Used for the mainfigure stage in multifigures and figure stage in simple figures
void P3DUserInterfaceCallback::startMainFigureStage(bool mainFigureStage)
{
    P3DControl * control;
    Fl_Value_Output * iterationOutput;
    Fl_Value_Output * imageMatchValueOutput;
    Fl_Value_Output * geometricPenaltyValueOutput;
    Fl_Value_Output * objectiveFunctionValueOutput;

    Fl_Button * startBtn;
    Fl_Button * nextBtn;
    Fl_Button * cancelBtn;
	Fl_Button * prevBtn;

    int i;
    int startIteration;

#ifdef OUTPUT_DEFORMING_MODELS
    M3DObjectFile objectFile;
    M3DObject * object;
    char filename[1024];
    char outputFilename[1024];
    char * modelFilename;
#endif


    if (optimizerIsRunning)
        return;

    if (userInterface == NULL)
        return;

    control = userInterface->control;
    if (control == NULL)
        return;

	if (mainFigureStage) {
		iterationOutput = userInterface->mainFigureIterationOutput;
		imageMatchValueOutput = userInterface->mainFigureImageMatchValueOutput;
		geometricPenaltyValueOutput = userInterface->mainFigurePenaltyValueOutput;
		objectiveFunctionValueOutput = userInterface->mainFigureObjectiveFunctionValueOutput;

		startBtn = userInterface->mainFigureStartBtn;
	}
	else {
		iterationOutput = userInterface->figureStageIterationOutput;
		imageMatchValueOutput = userInterface->figureStageImageMatchValueOutput;
		geometricPenaltyValueOutput = userInterface->figureStagePenaltyValueOutput;
		objectiveFunctionValueOutput = userInterface->figureStageObjectiveFunctionValueOutput;

		startBtn = userInterface->figureStageStartBtn;
	}

	if (iterationOutput == NULL || imageMatchValueOutput == NULL || geometricPenaltyValueOutput == NULL
		|| objectiveFunctionValueOutput == NULL)
			return;

    nextBtn = userInterface->nextOptimizationBtn;
    cancelBtn = userInterface->optimizationCancelBtn;
	prevBtn = userInterface->previewOptimizationBtn;
    if (startBtn == NULL || nextBtn == NULL || cancelBtn == NULL|| prevBtn == NULL)
		return;

    optimizerIsRunning = true;

    // Deactivate (grey) buttons we don't want to get hit
    startBtn->deactivate();
    nextBtn->deactivate();
    cancelBtn->deactivate();

    Fl::check();

    if (! previewMode)
        startIteration = (int) iterationOutput->value() + 1;
    else
        startIteration = 1;

    control->setEditingMode(NO_EDITING_MODE);
	control->clearPenaltyWeights();

	setMainFigureStagePenaltyWeight(mainFigureStage, optimizationType);
	setMainFigureStageConstraintsPenaltyWeight(mainFigureStage, optimizationType);
	if (! previewMode)
		if (mainFigureStage)
			cout << "Main figure" << flush;
		else
			cout << "Figure" << flush;

	bool binaryImageMatch = (control->type() == BINARY_IMAGE_MATCH);

	if (control->testAndResetModelChanged())
		resettingOptimizer = true;

	if (! previewMode && lastModeWasPreview)
		resettingOptimizer = true;

    if (resettingOptimizer) {
		if (userInterface->usePGABtn->value())
			control->initializeMainFigureOptimizer(binaryImageMatch, true);
        else
            control->initializeMainFigureOptimizer(binaryImageMatch);
        resettingOptimizer = false;
    }
    optimizerStopRequested = false;

	bool negate = true;		// AGG: this variable should be discarded
	if (binaryImageMatch)
		negate = false;
	int niter = MAX_OPTIMIZER_ITERATIONS;
	if (previewMode)
		niter = 2;
	else
		prevBtn->deactivate();	// Once START has been pressed, previewing is disallowed

    for(i = startIteration; i < niter; i++) {
		if (! previewMode) {
			control->doMainFigureIterations(1);
			iterationOutput->value(i);
#ifdef OPTIMIZATION_VISUALIZER
			((OptVisualizerUI *) ovUI)->updateChart();
#endif
		}
		geometricPenaltyValueOutput->value(control->getLastMainFigureGeometricMatchValue());
		if (negate) {
			objectiveFunctionValueOutput->value(-control->getLastMainFigureObjectiveFunctionValue());
			imageMatchValueOutput->value(-(control->getLastMainFigureObjectiveFunctionValue()
				- control->getLastMainFigureGeometricMatchValue()));
		}
		else {
			objectiveFunctionValueOutput->value(control->getLastMainFigureObjectiveFunctionValue());
			imageMatchValueOutput->value(control->getLastMainFigureObjectiveFunctionValue()
				- control->getLastMainFigureGeometricMatchValue());
		}

#ifdef OUTPUT_DEFORMING_MODELS
        // Output the current model to file
        modelFilename = control->getModelFileName();
        if(modelFilename == NULL)
            strcpy(filename, "registration");
        else
        {
            strcpy(filename, modelFilename);
            delete [] modelFilename;
        }

        sprintf(outputFilename, "%s/registration_%.4d.m3d", filename, i);

        object = control->getObjectPtr();
        if(object != NULL)
            objectFile.write(outputFilename, *object);
#endif
        update();
        Fl::check();
        if(optimizerStopRequested == true)
            break;

        if(control->mainFigureOptimizerIsFinished())
            break;
    }

    optimizerIsRunning = false;
	lastModeWasPreview = previewMode;
    updateVisibilityFiguresDisplay(false);
    control->setEditingMode(ANY_EDITING_MODE);

    // Reactivate buttons that we can use now
    startBtn->activate();
    nextBtn->activate();
    cancelBtn->activate();

	if (! previewMode)
		cout << " stage" << endl;
}

void P3DUserInterfaceCallback::stopMainFigureStage()
{
    optimizerStopRequested = true;
}

void P3DUserInterfaceCallback::setMainFigureStagePenaltyWeight(bool mainFigureStage,
	MatchType type)
{
    P3DControl * control;
    Fl_Value_Slider * penaltyWeightSlider;
	double wt;

    control = userInterface->control;
    if (control == NULL)
        return;

	if (mainFigureStage)
		penaltyWeightSlider = userInterface->mainFigurePenaltyWeightSlider;
	else
		penaltyWeightSlider = userInterface->figureStagePenaltyWeightSlider;
    if (penaltyWeightSlider == NULL)
        return;

	if (penaltyWeightSlider->active())
		wt = penaltyWeightSlider->value();
	else
		wt = 0.0;

	if (type == HISTO_IMAGE_MATCH) {
		// Using the geometry weight slider to enter histogram match penalty
		control->setMainFigurePenaltyWeight(GpFigureGeometryWeight, 0.0);
		control->setMainFigurePenaltyWeight(GpFigureHistogramMatch, wt);
		if (globalVerbosity > 0)
			cout << "Histogram match weight = " << wt << '\n';
	}
	else if (type == MULTIPATCH_MATCH) {
		// Using the geometry weight slider to enter histogram match penalty
		control->setMainFigurePenaltyWeight(GpFigureGeometryWeight, 0.0);
		control->setMainFigurePenaltyWeight(GpFigureMultiPatchMatch, wt);
		if (globalVerbosity > 0)
			cout << "MultiPatch match weight = " << wt << '\n';
	}
	else {
		control->setMainFigurePenaltyWeight(GpFigureGeometryWeight, wt);
		if (globalVerbosity > 0)
			cout << "Geometry weight = " << wt << '\n';
	}
}

void P3DUserInterfaceCallback::setMainFigureStageConstraintsPenaltyWeight(bool mainFigureStage,
	MatchType type)
{
    P3DControl * control;
    Fl_Value_Slider * constraintsPenaltyWeightSlider;
	double wt;

    control = userInterface->control;
    if (control == NULL)
        return;

	if (mainFigureStage)
		constraintsPenaltyWeightSlider =
			userInterface->mainFigureConstraintsPenaltyWeightSlider;
	else
		constraintsPenaltyWeightSlider =
			userInterface->figureStageConstraintsPenaltyWeightSlider;
    if (constraintsPenaltyWeightSlider == NULL)
        return;

    if (constraintsPenaltyWeightSlider->active())
		wt = constraintsPenaltyWeightSlider->value();
	else
		wt = 0.0;
	if (type == HISTO_IMAGE_MATCH) {
		// Using the inter-object weight slider to enter Mahalanobis penalty
		control->setMainFigurePenaltyWeight(GpFigureInterObjectWeight, 0.0);
		control->setMainFigurePenaltyWeight(FigureMahalanobisMatch, wt);
		if (globalVerbosity > 0)
			cout << "Mahalanobis weight = " << wt << '\n';
	}
	else {
		control->setMainFigurePenaltyWeight(GpFigureInterObjectWeight, wt);
		if (globalVerbosity > 0)
			cout << "Inter-object weight = " << wt << '\n';
	}
}

void P3DUserInterfaceCallback::clearMainFigureStage(bool mainFigureStage)
{
    Fl_Value_Output * iterationOutput;
    Fl_Value_Output * imageMatchValueOutput;
    Fl_Value_Output * objectiveFunctionValueOutput;

    Fl_Button * startBtn,
              * nextBtn,
              * cancelBtn;

    if (userInterface == NULL)
        return;

    if (optimizerIsRunning)
        return;

	if (mainFigureStage) {
		iterationOutput = userInterface->mainFigureIterationOutput;
		imageMatchValueOutput = userInterface->mainFigureImageMatchValueOutput;
		objectiveFunctionValueOutput = userInterface->mainFigureObjectiveFunctionValueOutput;

		startBtn = userInterface->mainFigureStartBtn;
	}
	else {
		iterationOutput = userInterface->figureStageIterationOutput;
		imageMatchValueOutput = userInterface->figureStageImageMatchValueOutput;
		objectiveFunctionValueOutput = userInterface->figureStageObjectiveFunctionValueOutput;

		startBtn = userInterface->figureStageStartBtn;
	}
	nextBtn = userInterface->nextOptimizationBtn;
	cancelBtn = userInterface->optimizationCancelBtn;

    if(iterationOutput != NULL)
        iterationOutput->value(0);
    if(imageMatchValueOutput != NULL)
        imageMatchValueOutput->value(0);
    if(objectiveFunctionValueOutput != NULL)
        objectiveFunctionValueOutput->value(0);

    if(startBtn != NULL)
        startBtn->activate();
    if(nextBtn != NULL)
        nextBtn->activate();
    if(cancelBtn != NULL)
        cancelBtn->activate();
}

void P3DUserInterfaceCallback::startSubfigureStage()
{
    P3DControl * control;
    Fl_Value_Output * iterationOutput;
    Fl_Value_Output * imageMatchValueOutput;
    Fl_Button * startBtn;
    Fl_Button * nextBtn;
    Fl_Button * cancelBtn;
	Fl_Button * prevBtn;

    int i;
    int startIteration;

#ifdef OUTPUT_DEFORMING_MODELS
    M3DObjectFile objectFile;
    M3DObject * object;
    char filename[1024];
    char outputFilename[1024];
    char * modelFilename;
#endif


    if (optimizerIsRunning)
        return;

    if (userInterface == NULL)
        return;

    control = userInterface->control;
    if (control == NULL)
        return;

    iterationOutput = userInterface->subfigureIterationOutput;
    imageMatchValueOutput = userInterface->subfigureMatchValueOutput;
    if (iterationOutput == NULL || imageMatchValueOutput == NULL)
		return;

    startBtn = userInterface->subfigureStartBtn;
    nextBtn = userInterface->nextOptimizationBtn;
    cancelBtn = userInterface->optimizationCancelBtn;
	prevBtn = userInterface->previewOptimizationBtn;
    if (startBtn == NULL || nextBtn == NULL || cancelBtn == NULL || prevBtn == NULL)
		return;

    optimizerIsRunning = true;

    // Deactivate (grey) buttons we don't want to get hit
    startBtn->deactivate();
    nextBtn->deactivate();
    cancelBtn->deactivate();

    Fl::check();

    if (! previewMode)
        startIteration = (int) iterationOutput->value() + 1;
    else
        startIteration = 1;

    control->setEditingMode(NO_EDITING_MODE);
	control->clearPenaltyWeights();

    setSubfigureStagePenaltyWeight();	// AGG: the penalty weights aren't handled right yet for this stage
	if (! previewMode)
		cout << "Subfigure" << flush;

	if (control->testAndResetModelChanged())
		resettingOptimizer = true;

	if (! previewMode && lastModeWasPreview)
		resettingOptimizer = true;

    if(resettingOptimizer) {
        control->initializeSubfigureOptimizer();
        resettingOptimizer = false;
    }

    optimizerStopRequested = false;

	int niter = MAX_OPTIMIZER_ITERATIONS;
	if (previewMode)
		niter = 2;
	else
		prevBtn->deactivate();	// Once START has been pressed, previewing is disallowed

    for (i = startIteration; i < niter; i++) {
		if (! previewMode) {
			control->doSubfigureIterations(1);
			iterationOutput->value(i);
		}
        imageMatchValueOutput->value(-control->getLastSubfigureMatchValue());
#ifdef OUTPUT_DEFORMING_MODELS
        // Output the current model to file
        modelFilename = control->getModelFileName();
        if (modelFilename == NULL)
            strcpy(filename, "registration");
        else {
            strcpy(filename, modelFilename);
            delete [] modelFilename;
        }

        sprintf(outputFilename, "%s/registration_%.4d.m3d", filename, i);

        object = control->getObjectPtr();
        if(object != NULL)
            objectFile.write(outputFilename, *object);
#endif
        update();
        Fl::check();
        if(optimizerStopRequested == true)
            break;

        if(control->subfigureOptimizerIsFinished())
            break;
    }

    optimizerIsRunning = false;
	lastModeWasPreview = previewMode;
    updateVisibilityFiguresDisplay(false);
    control->setEditingMode(ANY_EDITING_MODE);

    // Reactivate buttons that we can use now
    startBtn->activate();
    nextBtn->activate();
    cancelBtn->activate();

	if (! previewMode)
		cout << " stage" << endl;
}

void P3DUserInterfaceCallback::stopSubfigureStage()
{
    optimizerStopRequested = true;
}

void P3DUserInterfaceCallback::setSubfigureStagePenaltyWeight()
{
    P3DControl * control;
    Fl_Value_Slider * penaltyWeightSlider;

    control = userInterface->control;
    if(control == NULL)
        return;

    penaltyWeightSlider = userInterface->subfigurePenaltyWeightSlider;
    if(penaltyWeightSlider == NULL)
        return;

//    control->setSubfigurePenaltyWeight(penaltyWeightSlider->value());
}

void P3DUserInterfaceCallback::setSubfigureStageConstraintsPenaltyWeight()
{
	// Not yet implemented
}

void P3DUserInterfaceCallback::clearSubfigureStage()
{
    Fl_Value_Output * iterationOutput;
    Fl_Value_Output * imageMatchValueOutput;

    Fl_Button * startBtn,
              * nextBtn,
              * cancelBtn;

    if(userInterface == NULL)
        return;

    if(optimizerIsRunning)
        return;

    iterationOutput = userInterface->subfigureIterationOutput;
    imageMatchValueOutput = userInterface->subfigureMatchValueOutput;

    startBtn = userInterface->subfigureStartBtn;
    nextBtn = userInterface->nextOptimizationBtn;
    cancelBtn = userInterface->optimizationCancelBtn;

    if(iterationOutput != NULL)
        iterationOutput->value(0);
    if(imageMatchValueOutput != NULL)
        imageMatchValueOutput->value(0);

    if(startBtn != NULL)
        startBtn->activate();
    if(nextBtn != NULL)
        nextBtn->activate();
    if(cancelBtn != NULL)
        cancelBtn->activate();
}

void P3DUserInterfaceCallback::startPrimitiveStage()
{
    P3DControl * control;
    M3DObject * object;
	M3DFigure * currFigure;
    Fl_Value_Output * iterationOutput;
    Fl_Value_Output * imageMatchValueOutput;
    Fl_Value_Output * geometricPenaltyValueOutput;
    Fl_Value_Output * objectiveFunctionValueOutput;
    Fl_Button * startBtn;
    Fl_Button * nextBtn;
    Fl_Button * cancelBtn;
	Fl_Button * prevBtn;

    int i, figureId, numFigures;
    int startIteration;

#ifdef OUTPUT_DEFORMING_MODELS
    M3DObjectFile objectFile;
    char filename[1024];
    char outputFilename[1024];
    char * modelFilename;
#endif


    if (optimizerIsRunning)
        return;

    if (userInterface == NULL)
        return;

    control = userInterface->control;
    if (control == NULL)
        return;

    object = control->getObjectPtr();
    if (object == NULL)
        return;

    if (object->isAnySelected() == false) {
        cout << "Something must be selected for deformation" << endl;
		fl_beep(FL_BEEP_ERROR);
        return;
    }

    iterationOutput = userInterface->primitiveIterationOutput;
	objectiveFunctionValueOutput = userInterface->primitiveObjectiveFunctionValueOutput;
	geometricPenaltyValueOutput = userInterface->primativePenaltyValueOutput;
	imageMatchValueOutput = userInterface->primitiveImageMatchValueOutput;
	if (iterationOutput == NULL || objectiveFunctionValueOutput == NULL
		|| geometricPenaltyValueOutput == NULL || imageMatchValueOutput == NULL)
			return;

    startBtn = userInterface->primitiveStartBtn;
    nextBtn = userInterface->nextOptimizationBtn;
    cancelBtn = userInterface->optimizationCancelBtn;
    prevBtn = userInterface->previewOptimizationBtn;
	if (startBtn == NULL || nextBtn == NULL || cancelBtn == NULL || prevBtn == NULL)
		return;

    optimizerIsRunning = true;

    // Deactivate (grey) buttons we don't want to get hit
    startBtn->deactivate();
    nextBtn->deactivate();
    cancelBtn->deactivate();

    Fl::check();

	if (lastModeWasPreview == true)
		iterationOutput->value(0);
	if (! previewMode)
		startIteration = (int) iterationOutput->value() + 1;
	else
		startIteration = 0;

	bool negate = true;
	if (control->type() == BINARY_IMAGE_MATCH)
		negate = false;
	if (negate) {
		objectiveFunctionValueOutput->value(-control->getLastDeformationObjectiveFunctionValue());
		geometricPenaltyValueOutput->value(-control->getLastDeformationGeometricMatchValue());
		imageMatchValueOutput->value(-(control->getLastDeformationObjectiveFunctionValue()
			- control->getLastDeformationGeometricMatchValue()));
	}
	else {
		objectiveFunctionValueOutput->value(control->getLastDeformationObjectiveFunctionValue());
		geometricPenaltyValueOutput->value(control->getLastDeformationGeometricMatchValue());
		imageMatchValueOutput->value(control->getLastDeformationObjectiveFunctionValue()
			- control->getLastDeformationGeometricMatchValue());
	}

    control->setEditingMode(NO_EDITING_MODE);
	control->clearPenaltyWeights();

    setPrimitiveStagePenaltyWeight();
    setPrimitiveStageConstraintsPenaltyWeight();
    setPrimitiveStageNeighborPenaltyWeight();
	if (! previewMode)
		cout << "Atom" << flush;

	if (control->testAndResetModelChanged())
		resettingOptimizer = true;

	if (! previewMode && lastModeWasPreview)
		resettingOptimizer = true;

    if (resettingOptimizer) {
        control->initializeDeformationOptimizer();
        resettingOptimizer = false;
    }

    optimizerStopRequested = false;

	// Find the (selected) figure to be optimized
	M3DDeformationOptimizer * deformationOptimizer = control->getDeformationOptimizer();
	M3DObject * candidateObject = deformationOptimizer->getCandidateObject();
	numFigures = candidateObject->getFigureCount();
    for (figureId = 0; figureId < numFigures; figureId++) {
        currFigure = candidateObject->getFigurePtr(figureId);
        if (currFigure != NULL && currFigure->isAnySelected())
			break;
    }

	int niter = MAX_OPTIMIZER_ITERATIONS;
	if (previewMode)
		niter = currFigure->getPrimitiveCount();
	else
		prevBtn->deactivate();	// Once START has been pressed, previewing is disallowed

	// Optimize the figure
	if (figureId < numFigures)
		for(i = startIteration; i < niter; i++)		// Essentially a while(1) loop
		{
			if (control->doDeformationIterations(figureId, previewMode, true) == true)	// Optimize selected atoms
				optimizerStopRequested = true;

			iterationOutput->value(i);
			if (negate) {
				objectiveFunctionValueOutput->value(-control->getLastDeformationObjectiveFunctionValue(previewMode));
				imageMatchValueOutput->value(-(control->getLastDeformationObjectiveFunctionValue(previewMode)
					- control->getLastDeformationGeometricMatchValue()));
				geometricPenaltyValueOutput->value(-control->getLastDeformationGeometricMatchValue());
			}
			else {
				objectiveFunctionValueOutput->value(control->getLastDeformationObjectiveFunctionValue(previewMode));
				imageMatchValueOutput->value(control->getLastDeformationObjectiveFunctionValue(previewMode)
					- control->getLastDeformationGeometricMatchValue());
				geometricPenaltyValueOutput->value(control->getLastDeformationGeometricMatchValue());
			}

#ifdef OUTPUT_DEFORMING_MODELS
			// Output the current model to file
			modelFilename = control->getModelFileName();
			if(modelFilename == NULL)
				strcpy(filename, "deformation");
			else
			{
				strcpy(filename, modelFilename);
				delete [] modelFilename;
			}

			sprintf(outputFilename, "%s/deformation_%.4d.m3d", filename, i);

			objectFile.write(outputFilename, *object);
#endif
			update();
			Fl::check();
            if(optimizerStopRequested == true) 
				break;
		}

    optimizerIsRunning = false;
	lastModeWasPreview = previewMode;
    updateVisibilityFiguresDisplay(false);
    control->setEditingMode(ANY_EDITING_MODE);

    // Reactivate buttons that we can use now
    startBtn->activate();
    nextBtn->activate();
    cancelBtn->activate();

	if (! previewMode)
		cout << " stage";
	cout << endl;	// Needed because the deformation optimizer prints dots as it runs.
}

void P3DUserInterfaceCallback::stopPrimitiveStage()
{
    optimizerStopRequested = true;
}

void P3DUserInterfaceCallback::setPrimitiveStagePenaltyWeight()
{
    P3DControl * control;
    Fl_Value_Slider * penaltyWeightSlider;
	double wt;

    control = userInterface->control;
    if (control == NULL)
        return;

    penaltyWeightSlider = userInterface->primitivePenaltyWeightSlider;
    if (penaltyWeightSlider == NULL)
        return;

	if (penaltyWeightSlider->active())
		wt = penaltyWeightSlider->value();
	else
		wt = 0.0;
	control->setDeformationPenaltyWeight(GpAtomGeometryWeight, wt);
	if (globalVerbosity > 0)
		cout << "Geometry weight = " << wt << '\n';
}

void P3DUserInterfaceCallback::setPrimitiveStageConstraintsPenaltyWeight()
{
    P3DControl * control;
    Fl_Value_Slider * constraintsPenaltyWeightSlider;
	double wt;

    control = userInterface->control;
    if (control == NULL)
        return;

    constraintsPenaltyWeightSlider =
        userInterface->primitiveConstraintsPenaltyWeightSlider;
    if (constraintsPenaltyWeightSlider == NULL)
        return;

    if (constraintsPenaltyWeightSlider->active())
		wt = constraintsPenaltyWeightSlider->value();
	else
		wt = 0.0;
	control->setDeformationPenaltyWeight(GpAtomInterWeight, wt);
	if (globalVerbosity > 0)
		cout << "Inter-atom weight = " << wt << '\n';
}

void P3DUserInterfaceCallback::setPrimitiveStageNeighborPenaltyWeight()
{
    P3DControl * control;
    Fl_Value_Slider * neighborPenaltyWeightSlider;
	double wt;

    control = userInterface->control;
    if (control == NULL)
        return;

    neighborPenaltyWeightSlider = userInterface->primitiveNeighborPenaltyWeightSlider;
    if (neighborPenaltyWeightSlider == NULL)
        return;

	if (neighborPenaltyWeightSlider->active())
		wt = neighborPenaltyWeightSlider->value();
	else
		wt = 0.0;
	control->setDeformationPenaltyWeight(GpAtomNeighborWeight, wt);
	if (globalVerbosity > 0)
		cout << "Neighbor penalty weight = " << wt << '\n';
}

void P3DUserInterfaceCallback::clearPrimitiveStage()
{
    Fl_Value_Output * iterationOutput;

    Fl_Button * startBtn,
              * nextBtn,
              * cancelBtn;

    if(userInterface == NULL)
        return;

    if(optimizerIsRunning)
        return;

    iterationOutput = userInterface->primitiveIterationOutput;

    startBtn = userInterface->primitiveStartBtn;
    nextBtn = userInterface->nextOptimizationBtn;
    cancelBtn = userInterface->optimizationCancelBtn;

    if(iterationOutput != NULL)
        iterationOutput->value(0);

    if(startBtn != NULL)
        startBtn->activate();
    if(nextBtn != NULL)
        nextBtn->activate();
    if(cancelBtn != NULL)
        cancelBtn->activate();
}

void P3DUserInterfaceCallback::startVoxelStage()
{
    P3DControl * control;
    Fl_Value_Output * iterationOutput;
    Fl_Value_Output * matchValueOutput;
    Fl_Value_Output * penaltyValueOutput;
    Fl_Value_Output * objectiveFunctionValueOutput;
    Fl_Button * startBtn;
    Fl_Button * nextBtn;
    Fl_Button * cancelBtn;
	Fl_Button * prevBtn;

    int i;
    int startIteration;


    if (optimizerIsRunning)
        return;

    if (userInterface == NULL)
        return;

    control = userInterface->control;
    if( control == NULL)
        return;

    if (optimizerIsRunning)
        return;

    iterationOutput = userInterface->voxelStageIterationOutput;
    matchValueOutput = userInterface->voxelStageMatchValueOutput;
    penaltyValueOutput = userInterface->voxelStagePenaltyValueOutput;
    objectiveFunctionValueOutput = userInterface->voxelStageObjectiveFunctionValueOutput;
	if (iterationOutput == NULL || matchValueOutput == NULL
		|| penaltyValueOutput == NULL || objectiveFunctionValueOutput == NULL)
			return;

    startBtn = userInterface->voxelStageStartBtn;
    nextBtn = userInterface->nextOptimizationBtn;
    cancelBtn = userInterface->optimizationCancelBtn;
	prevBtn = userInterface->previewOptimizationBtn;
	if (startBtn == NULL || nextBtn == NULL || cancelBtn == NULL || prevBtn == NULL)
		return;

	optimizerIsRunning = true;
	if (! previewMode)
		cout << "Boundary" << flush;

	if (control->testAndResetModelChanged())
		resettingOptimizer = true;

	if (! previewMode && lastModeWasPreview)
		resettingOptimizer = true;

    if (resettingOptimizer) {
		int surfaceLevel = userInterface->modelView->getSurfaceLevel();
        bool initFlag = control->initializeVoxelOptimizer(surfaceLevel);
		if (! initFlag) {
			optimizerIsRunning = false;
			return;
		}

        resettingOptimizer = false;
    }

    control->setEditingMode(NO_EDITING_MODE);

    // Deactivate (grey) buttons we don't want to get hit
    startBtn->deactivate();
    nextBtn->deactivate();
    cancelBtn->deactivate();

    Fl::check();

	if (lastModeWasPreview == true)
		iterationOutput->value(0);
    if (! previewMode && iterationOutput != NULL)
        startIteration = (int) iterationOutput->value() + 1;
    else
        startIteration = 1;

    setVoxelStagePenaltyWeight();
    setVoxelStageConstraintsPenaltyWeight();

    optimizerStopRequested = false;

	int niter = MAX_OPTIMIZER_ITERATIONS;
	if (previewMode)
		niter = 2;
	else
		prevBtn->deactivate();	// Once START has been pressed, previewing is disallowed

    for (i = startIteration; i < niter; i++) {
        control->doVoxelIterations(100, previewMode);

        iterationOutput->value(i);
		matchValueOutput->value(-control->getLastVoxelMatchValue());
		objectiveFunctionValueOutput->
			value(-control->getLastVoxelObjectiveFunctionValue());
		penaltyValueOutput->value(control->getLastVoxelPenaltyValue());

        update();
        Fl::check();
        if(optimizerStopRequested == true)
            break;

        if(control->voxelStageOptimizerIsFinished())
            break;
    }

    optimizerIsRunning = false;
	lastModeWasPreview = previewMode;
    updateVisibilityFiguresDisplay(false);
    control->setEditingMode(ANY_EDITING_MODE);

    // Reactivate buttons that we can use now
    startBtn->activate();
    nextBtn->activate();
    cancelBtn->activate();

	if (! previewMode)
		cout << " stage" << endl;
}

void P3DUserInterfaceCallback::stopVoxelStage()
{
    optimizerStopRequested = true;
}

void P3DUserInterfaceCallback::setVoxelStagePenaltyWeight()
{
    P3DControl * control;
    Fl_Value_Slider * penaltyWeightSlider;

    control = userInterface->control;
    if(control == NULL)
        return;

    penaltyWeightSlider = userInterface->voxelStagePenaltyWeightSlider;
    if(penaltyWeightSlider == NULL)
        return;

	control->setVoxelPenaltyWeight(penaltyWeightSlider->value());
}

void P3DUserInterfaceCallback::setVoxelStageConstraintsPenaltyWeight()
{
	P3DControl * control;
    Fl_Value_Slider * penaltyWeightSlider;

    control = userInterface->control;
    if(control == NULL)
        return;

    penaltyWeightSlider = userInterface->voxelStageConstraintsPenaltyWeightSlider;
    if(penaltyWeightSlider == NULL)
        return;

	if (penaltyWeightSlider->active())
		control->setVoxelConstraintsPenaltyWeight(penaltyWeightSlider->value());
	else
		control->setVoxelConstraintsPenaltyWeight(0.0);
}

void P3DUserInterfaceCallback::clearVoxelStage()
{
    Fl_Value_Output * iterationOutput;
    Fl_Value_Output * matchValueOutput;
    Fl_Value_Output * penaltyValueOutput;
    Fl_Value_Output * objectiveFunctionValueOutput;

    Fl_Button * startBtn,
              * nextBtn,
              * cancelBtn;

    if(userInterface == NULL)
        return;

    if(optimizerIsRunning)
        return;

    iterationOutput = userInterface->voxelStageIterationOutput;
    matchValueOutput = userInterface->voxelStageMatchValueOutput;
    penaltyValueOutput = userInterface->voxelStagePenaltyValueOutput;
	objectiveFunctionValueOutput = userInterface->voxelStageObjectiveFunctionValueOutput;

    startBtn = userInterface->voxelStageStartBtn;
    nextBtn = userInterface->nextOptimizationBtn;
    cancelBtn = userInterface->optimizationCancelBtn;

    if(iterationOutput != NULL)
        iterationOutput->value(0);
    if(matchValueOutput != NULL)
        matchValueOutput->value(0);
    if(penaltyValueOutput != NULL)
        penaltyValueOutput->value(0);
    if(objectiveFunctionValueOutput != NULL)
        objectiveFunctionValueOutput->value(0);

    if(startBtn != NULL)
        startBtn->activate();
    if(nextBtn != NULL)
        nextBtn->activate();
    if(cancelBtn != NULL)
        cancelBtn->activate();
}

void P3DUserInterfaceCallback::initializeOptimizerSettings()
{
    P3DControl * control;

    control = userInterface->control;
    if(control == NULL)
        return;

	// Set default penalty weights
    userInterface->normalPenaltyInput->value(DEFAULT_NORMAL_PENALTY_WEIGHT);
    userInterface->slidingPenaltyInput->value(DEFAULT_SLIDING_PENALTY_WEIGHT);
    userInterface->penetrationPenaltyInput->value(DEFAULT_PENETRATION_PENALTY_WEIGHT);
    userInterface->creasingPenaltyInput->value(DEFAULT_CREASING_PENALTY_WEIGHT);

    userInterface->skipMultiplierInput->value(DEFAULT_SCHEDULE_MULTIPLIER);
    userInterface->skipThresholdInput->value(DEFAULT_THRESHOLD);
    userInterface->conjGradIterInput->value(DEFAULT_NUM_CONJGRAD_ITERATIONS);

    if(pgaFilename != NULL)
    {
        if(control->readPGAFile(pgaFilename))
            userInterface->pgaFileOutput->value(pgaFilename);
        else
        {
            delete [] pgaFilename;
            pgaFilename = NULL;
        }
    }
}

#else	/* BINARY */

void P3DUserInterfaceCallback::setOptimizerIntensityWindowing() {}
void P3DUserInterfaceCallback::previewOptimizationStage() {}

OptimizerStage P3DUserInterfaceCallback::nextOptimizationStage(bool)
{
	return NO_OPTIMIZER_STAGE;
}

void P3DUserInterfaceCallback::cancelOptimization() {}
void P3DUserInterfaceCallback::startModelStage(bool) {}
void P3DUserInterfaceCallback::stopModelStage() {}
void P3DUserInterfaceCallback::resetModelStage(bool) {}
void P3DUserInterfaceCallback::setModelStagePenaltyWeight(bool) {}
void P3DUserInterfaceCallback::startMainFigureStage(bool) {}
void P3DUserInterfaceCallback::stopMainFigureStage() {}
void P3DUserInterfaceCallback::setMainFigureStagePenaltyWeight(bool,
	MatchType) {}
void P3DUserInterfaceCallback::setMainFigureStageConstraintsPenaltyWeight(
	bool, MatchType) {}
void P3DUserInterfaceCallback::startSubfigureStage() {}
void P3DUserInterfaceCallback::stopSubfigureStage() {}
void P3DUserInterfaceCallback::setSubfigureStagePenaltyWeight() {}
void P3DUserInterfaceCallback::setSubfigureStageConstraintsPenaltyWeight() {}
void P3DUserInterfaceCallback::startPrimitiveStage() {}
void P3DUserInterfaceCallback::stopPrimitiveStage() {}
void P3DUserInterfaceCallback::setPrimitiveStagePenaltyWeight() {}
void P3DUserInterfaceCallback::setPrimitiveStageConstraintsPenaltyWeight() {}
void P3DUserInterfaceCallback::setPrimitiveStageNeighborPenaltyWeight() {}
void P3DUserInterfaceCallback::startVoxelStage() {}
void P3DUserInterfaceCallback::stopVoxelStage() {}
void P3DUserInterfaceCallback::setVoxelStagePenaltyWeight() {}
void P3DUserInterfaceCallback::setVoxelStageConstraintsPenaltyWeight() {}

#endif	/* BINARY */

// ----------------------------- End of Optimizer Methods --------------------------------

#ifdef BINARY

void P3DUserInterfaceCallback::setPenaltyWeights() {
}

void P3DUserInterfaceCallback::revertPenaltyWeights() {
}

#else	/* BINARY */

void P3DUserInterfaceCallback::setPenaltyWeights()
{
    P3DControl * control;

    if(userInterface == NULL)
        return;

    control = userInterface->control;
    if(control == NULL)
        return;

    control->setNormalPenaltyWeight(userInterface->normalPenaltyInput->value());
    control->setSlidingPenaltyWeight(userInterface->slidingPenaltyInput->value());
    control->setPenetrationPenaltyWeight(userInterface->penetrationPenaltyInput->value());
    control->setCreasingPenaltyWeight(userInterface->creasingPenaltyInput->value());

    M3DDeformationOptimizer * defOpt = control->getDeformationOptimizer();
    if (defOpt != NULL) {
        defOpt->setScheduleMultiplier((int) userInterface->skipMultiplierInput->value());
        defOpt->setThreshold(userInterface->skipThresholdInput->value());
        defOpt->setConjGradientIterations((int) userInterface->conjGradIterInput->value());
    }
}

void P3DUserInterfaceCallback::revertPenaltyWeights()
{
    P3DControl * control;

    if(userInterface == NULL)
        return;

    control = userInterface->control;
    if(control == NULL)
        return;

    userInterface->normalPenaltyInput->value(control->getNormalPenaltyWeight());
    userInterface->slidingPenaltyInput->value(control->getSlidingPenaltyWeight());
    userInterface->penetrationPenaltyInput->value(control->getPenetrationPenaltyWeight());
    userInterface->creasingPenaltyInput->value(control->getCreasingPenaltyWeight());

    M3DDeformationOptimizer * defOpt = control->getDeformationOptimizer();
    if(defOpt != NULL)
    {
        userInterface->skipMultiplierInput->value(defOpt->getScheduleMultiplier());
        userInterface->skipThresholdInput->value(defOpt->getThreshold());
        userInterface->conjGradIterInput->value(defOpt->getConjGradientIterations());
    }
    else
    {
        userInterface->skipMultiplierInput->value(-1);
        userInterface->skipThresholdInput->value(-1);
        userInterface->conjGradIterInput->value(-1);
    }

    Fl::check();
}

#endif	/* BINARY */

void P3DUserInterfaceCallback::loadPGAFile()
{
    P3DControl * control;
    const char * filename;

    control = userInterface->control;
    if (control == NULL)
        return;

    filename = askSingleFilename(ModelDirectory, "Load PGA File", "*.m3d", false);

	if (filename != NULL)
    {
        if (! control->readPGAFile(filename))
            fl_alert("File \'%s\' does not contain PGA statistics.", filename);
        else
        {
#ifndef BINARY
			// Loaded PGA data
			finishPGALoading(filename);
#endif
            if (pgaFilename != NULL)
                delete [] pgaFilename;
            pgaFilename = new char[strlen(filename) + 1];
            strcpy(pgaFilename, filename);
		    userInterface->pgaFileOutput->value(pgaFilename);

#ifdef BINARY
			initPGADeform();
#else
			activatePGASliders();
#endif
        }
    }
}

// ----------------------------- Start of Constraints Code --------------------------------

// Called when the constraint distance slider is changed
void P3DUserInterfaceCallback::setInterObjectDistance()
{
    P3DControl * control;

    control = userInterface->control;
    if(control == NULL)
        return;

    partialSurfaceSliderChanged = true;

    // Trigger draw() in P3DView.cpp to be called.  The inter-object
    // distance is used there for the partial rendering.
    updateSurface();
}

void P3DUserInterfaceCallback::updateInterObjectFigureSlider(int figureNumber) {
    P3DControl * control;
	const char blank[2] = " ";
	const char all[] = "(all figures)";

    userInterface->constraintsFigureSlider->value(figureNumber);

	if (figureNumber < 0) {
        userInterface->constraintsFigureSliderText->value(all);
		return;
	}

    control = userInterface->control;
    if(control == NULL)
        return;

    M3DObject * object = control->getObjectPtr();
    if (object == NULL)
        return;

    M3DFigure * figure = object->getFigurePtr(figureNumber);
    const char * name = figure->getName();
	if (name)
        userInterface->constraintsFigureSliderText->value(name);
	else
        userInterface->constraintsFigureSliderText->value(blank);
}

// Called when the figure slider is changed
void P3DUserInterfaceCallback::setInterObjectFigure()
{
    P3DControl * control;
    int markedPrimitiveId;
    int markedFigureId;

    control = userInterface->control;
    if(control == NULL)
        return;

    // The slider is never allowed to rest on the marked figure's number
    int val = (int) userInterface->constraintsFigureSlider->value();
    control->getRelativeMarkedPrimitiveId(markedFigureId, markedPrimitiveId);
    if (val == markedFigureId) {
        if (markedFigureId == userInterface->constraintsFigureSlider->maximum()) {
			updateInterObjectFigureSlider(currentFigure);
            return;
        }
        if (val > currentFigure) val++;
        else val--;
    }
    currentFigure = val;
	updateInterObjectFigureSlider(currentFigure);

    selectedPartialFigures->displayFigure(currentFigure);   // For rendering
    setDistanceSliderFromFigureSlider();
    updateSurface();
}

// Called from the user interface when the constraints panel is activated
bool P3DUserInterfaceCallback::constraintsWindowSetup()
{
    P3DControl * control;
    int markedPrimitiveId;

    control = userInterface->control;
    if (control == NULL)
        return false;

    M3DObject * object = control->getObjectPtr();
    if (object == NULL)
        return false;

    if (object->getFigureCount() < 2) {
        cout << "There must be two or more figures to use constraints" << endl;
		fl_beep(FL_BEEP_ERROR);
        return false;
    }

    // Test that there is a marked primitive
    markedPrimitiveId = control->getMarkedPrimitiveId();
    if (markedPrimitiveId == INVALID_PRIMITIVE_ID) {
        cout <<
            "Constraint specification requires an atom in the central object to be marked" 
            << endl;
		fl_beep(FL_BEEP_ERROR);
        return false;
    }

    // The marked figure cannot be omitted from the figure slider
    userInterface->constraintsFigureSlider->maximum(object->getFigureCount() - 1);

    // Allocate data structure for storing interfigural data
    if (selectedPartialFigures == NULL) {
        selectedPartialFigures = new SelectedPartialFigures;
        selectedPartialFigures->color(partialSurfaceColor);
    }

    return initializeConstraints();
}

// Called whenever the "Set governor" button is pressed and from constraintsWindowSetup()
bool P3DUserInterfaceCallback::initializeConstraints()
{
    P3DControl * control;
    int markedPrimitiveId;
    int markedFigureId;

#ifdef DEBUG
	cout << "P3DUserInterfaceCallback::initializeConstraints()\n";
#endif
    control = userInterface->control;
    if (control == NULL)
        return false;

    M3DObject * object = control->getObjectPtr();
    if (object == NULL)
        return false;

    // Get the ID of the marked primitive and the figure containing it
    control->getRelativeMarkedPrimitiveId(markedFigureId, markedPrimitiveId);
    if (markedPrimitiveId == INVALID_PRIMITIVE_ID) {
        cout <<
            "Constraint specification requires an atom in the central object to be marked" 
            << endl;
		fl_beep(FL_BEEP_ERROR);
        return false;
    }

    int num_figures = object->getFigureCount();
#ifdef DEBUG
    cout << "Object has " << num_figures << " figures." << endl;
#endif
    selectedPartialFigures->initialize(num_figures - 1);

    // Obtain and set the default slider value
    currentFigure = (int) userInterface->constraintsFigureSlider->value();
    // The figure slider value is never the marked figure number
    if (currentFigure == markedFigureId) {
        currentFigure++;
        if (currentFigure > userInterface->constraintsFigureSlider->maximum())
            currentFigure = currentFigure - 2;
    }
    updateInterObjectFigureSlider(currentFigure);
    selectedPartialFigures->displayFigure(currentFigure);
    selectedPartialFigures->markedFigure(markedFigureId);

    // Get a pointer to the marked figure
#ifdef DEBUG
    cout << "Marked figure is number " << markedFigureId << endl;
#endif
    markedFigurePtr = object->getFigurePtr(markedFigureId);
    if(markedFigurePtr == NULL)
        return false;           // Should never happen (see above)

    // Create an object containing a copy of the marked figure
    M3DObject * tmp_object = new M3DObject;
    M3DFigure * tmp_figure = markedFigurePtr->assign();
    tmp_object->addFigure(tmp_figure);

    // Get the rendering level
    if (userInterface->modelView == NULL)
        return false;
    int level = userInterface->modelView->getPartialSurfaceLevel();

    // Create a distance map from that object's boundary surface
    BoundingBox bbox;
    int mapSize[3];

    map.defineBBox(object, bbox, mapSize, CONSTRAINTS_RESOLUTION_BIAS + level);
    map.initializeWithBBox(tmp_object, &bbox, mapSize, CONSTRAINTS_RESOLUTION_BIAS + level);
    map.createMap(true);    // Produce an r-normalized map
    const int * map_size_ptr = map.getMapSize();
    if (map_size_ptr[0] <= 0) {
        delete tmp_object;
        cout << "map_size of 0!" << endl;
        return false;
    }
    if (map.status() == false) {
        delete tmp_object;
        cout << "No distance map!" << endl;
        return false;
    }
    delete tmp_object;

    // Output constraint count to Constraints window
    InterfiguralConstraints & ifc = markedFigurePtr->constraints();
    userInterface->constraintCount->value((double) ifc.size());
	if (ifc.size() == 0)
		userInterface->constraintCount->textcolor(FL_RED);
	else
		userInterface->constraintCount->textcolor(FL_BLACK);

    showPenetrationArea();
    control->setPartialDrawingMode(true);
    updateSurface();

    return true;

}

double P3DUserInterfaceCallback::getInterObjectDistance()
{
    P3DControl * control;
    Fl_Slider * interObjectDistanceSlider;
    Fl_Value_Output * interObjectDistanceValue;

    control = userInterface->control;
    if (control == NULL)
        return -1.0;

    interObjectDistanceSlider = userInterface->interObjectDistanceSlider;
    if (interObjectDistanceSlider == NULL)
        return -1.0;

    interObjectDistanceValue = userInterface->interObjectDistanceValue;
    if (interObjectDistanceValue == NULL)
        return -1.0;

    interObjectDistanceValue->value(interObjectDistanceSlider->value());
    return interObjectDistanceSlider->value();
}

void P3DUserInterfaceCallback::setPartialSurfaceStyle(M3DSurfaceStyle style) {
    if(userInterface == NULL)
        return;

    if(userInterface->modelView == NULL)
        return;

    userInterface->modelView->setPartialSurfaceStyle(style);
    updateSurface();
}

// Called when any of the Level buttons on the constraints window are pushed
void P3DUserInterfaceCallback::setPartialSurfaceLevel(int val)
{
    if(userInterface == NULL)
        return;

    if(userInterface->modelView == NULL)
        return;

    userInterface->modelView->setPartialSurfaceLevel(val);

    showPenetrationArea();  // AGG: this is buggy: changing the level can cause the slider to move!
    updateSurface();
}

// Called when the Set button is pressed
void P3DUserInterfaceCallback::setInterObjectConstraint() {
    P3DControl * control;

    control = userInterface->control;
    if (control == NULL)
        return;

    // Get the ID of the marked primitive and the figure containing it
    int markedFigureId, markedPrimitiveId;
    control->getRelativeMarkedPrimitiveId(markedFigureId, markedPrimitiveId);

    // Save specified constraint
    double slider_dist = getInterObjectDistance();
    InterfiguralConstraints & ifc = markedFigurePtr->constraints();
    cout << "For governor " << markedFigureId << ", ";
    if (currentFigure < 0) {
        M3DObject * object = control->getObjectPtr();
        int num_figures = object->getFigureCount();
        if (ifc.size() == 0) {
            for (int figureId = 0; figureId < num_figures; figureId++) {
                if (figureId == markedFigureId) continue;
                ifc.addFigure(figureId, (float) slider_dist);
            }
        }
        else {
            for (int figureId = 0; figureId < num_figures; figureId++) {
                if (figureId == markedFigureId) continue;
                if (! ifc.updateFigure(figureId, (float) slider_dist))
                    ifc.addFigure(figureId, (float) slider_dist);
            }
        }
        cout << "constraints of all figures set to " << slider_dist << endl;
    }
    else {
        if (ifc.size() == 0) {
            ifc.addFigure(currentFigure, (float) slider_dist);
        }
        else {
            if (! ifc.updateFigure(currentFigure, (float) slider_dist))
                ifc.addFigure(currentFigure, (float) slider_dist);
        }
        cout << "figure " << currentFigure << " constraint set to " << slider_dist << endl;
    }
    userInterface->constraintCount->value((double) ifc.size());
	if (ifc.size() == 0)
		userInterface->constraintCount->textcolor(FL_RED);
	else
		userInterface->constraintCount->textcolor(FL_BLACK);
    control->modelNeedsSaving();
    if (userInterface->modelView == NULL)
        return;
	userInterface->modelView->updateObject();	// Needed to render the new constraint
}

// Called when the constraints window is deactivated
void P3DUserInterfaceCallback::exitConstraintsMode() {
    P3DControl * control;

    control = userInterface->control;
    if (control == NULL)
        return;

    control->setPartialDrawingMode(false);
    updateSurface();    // AGG: this assumes the model is not edited during the
                        // constraints operation, so that any implicit surface
                        // remains valid.  It would be a good idea to assure this.
}

// Called when the "Clear all" button is pressed
void P3DUserInterfaceCallback::clearAllConstraints() {
    P3DControl * control;
	M3DObject * object;
    M3DFigure * figurePtr;

    control = userInterface->control;
    if (control == NULL)
        return;

    object = control->getObjectPtr();
    if(object == NULL)
        return;

	bool found = false;
    for (int i = 0; i < object->getFigureCount(); i++) {
		figurePtr = object->getFigurePtr(i);
		InterfiguralConstraints & ifc = figurePtr->constraints();
		if (ifc.size() > 0) {
			ifc.clear();
			found = true;
		}
	}
	if (! found)
		return;

    userInterface->constraintCount->value(0.0);
    userInterface->constraintCount->textcolor(FL_RED);
	setDistanceSliderFromFigureSlider();
    updateSurface();
    control->modelNeedsSaving();
}

// Called when the "Clear" button is pressed
void P3DUserInterfaceCallback::clearPartialConstraints() {
    P3DControl * control;

    control = userInterface->control;
    if (control == NULL)
        return;

    M3DObject * object = control->getObjectPtr();
    if(object == NULL)
        return;

    int markedPrimitiveId = control->getMarkedPrimitiveId();
    if (markedPrimitiveId == INVALID_PRIMITIVE_ID) {
        cout <<
            "Clearing a constraint specification requires an atom in the central object to be marked" 
            << endl;
		fl_beep(FL_BEEP_ERROR);
        return;
    }

    int markedFigureId;

    // Get the ID of the marked primitive and the figure containing it
    control->getRelativeMarkedPrimitiveId(markedFigureId, markedPrimitiveId);
    cout << "Clearing constraints for figure " << markedFigureId << endl;
    markedFigurePtr = object->getFigurePtr(markedFigureId);

    InterfiguralConstraints & ifc = markedFigurePtr->constraints();
    ifc.clear();
    userInterface->constraintCount->value(0.0);
    userInterface->constraintCount->textcolor(FL_RED);
	setDistanceSliderFromFigureSlider();
    updateSurface();
    control->modelNeedsSaving();
}

// Called from initializeConstraints() and setPartialSurfaceLevel().
// This function recomputes the distance lists stored in selectedPartialFigures. 
void P3DUserInterfaceCallback::showPenetrationArea() {
    int i, j;
    P3DControl * control;

#ifdef DEBUG
    cout << "P3DUserInterfaceCallback::showPenetrationArea() called\n";
#endif
    if (userInterface == NULL)
        return;

    control = userInterface->control;
    if (control == NULL)
        return;

    M3DObject * object = control->getObjectPtr();
    if(object == NULL)
        return;

    // Get the ID of the marked primitive and the figure containing it
    int markedPrimitiveId, markedFigureId;
    control->getRelativeMarkedPrimitiveId(markedFigureId, markedPrimitiveId);
    // The user accidentally may have no marked atom
    if (markedPrimitiveId == INVALID_PRIMITIVE_ID) {
        cout <<
            "Constraint specification requires an atom in the central object to be marked" 
            << endl;
		fl_beep(FL_BEEP_ERROR);
        return;                 // AGG: Perhaps should call exitConstraintsMode() here
    }
    userInterface->Governor->value(markedFigureId);     // Output to constraints window
#ifdef DEBUG
    cout << "P3DUserInterfaceCallback::showPenetrationArea(): markedPrimitiveId = "
        << markedPrimitiveId << endl;
#endif

    using namespace ThallCode;

    int numBPoints, numQuads;
    Bpoint * bpoints;
    Bpoint * quads;
    Pointlist_server2 ** pLists =
        (Pointlist_server2 **) selectedPartialFigures->pointLists();

    // Get the rendering level
    if (userInterface->modelView == NULL)
        return;
    int level = userInterface->modelView->getPartialSurfaceLevel();

    // Compute vectors of distances from boundary points to marked figure
    // for all figures, except the marked figure.
    int figureId;
    double min = 1.0;
    double max = 0.0;
    int num_figures = object->getFigureCount();
    i = 0;
    // AGG: This puts all figures, except the marked figure, on the selectedPartialFigures
    // list.  It would be faster to only put the currently selected figure on the list,
    // adding others as they are selected by the figure slider.
    for (figureId = 0; figureId < num_figures; figureId++)
    {
        M3DFigure * f = object->getFigurePtr(figureId);
        if (f == NULL) continue;
        if (figureId == markedFigureId) continue;   // Ignore the governor
        selectedPartialFigures->figure(i, figureId);

        Xferlist * xferList = convertM3DtoXfer(f);

        // Compute boundary quadrilaterals for each figure selected.  The
        // number of levels must be 'level + CONSTRAINTS_RESOLUTION_BIAS'
		// to correspond with these same calls in
		// M3DObjectSurfaceRenderer::partial_render(), called
        // from P3DView.cpp.
        pLists[i]->init(xferList);
#ifdef DEBUG
        cout << "Initializing point list for figure " << figureId << " at level "
            << level + CONSTRAINTS_RESOLUTION_BIAS << endl;
#endif
        pLists[i]->ComputeSubdivBoundaryTiles(level + CONSTRAINTS_RESOLUTION_BIAS);
        pLists[i]->subdivtileinfo(&numQuads, &quads);
        pLists[i]->ComputeSubdivPointCloud(level + CONSTRAINTS_RESOLUTION_BIAS);
        pLists[i]->subdivboundaryinfo(&numBPoints, &bpoints);
#ifdef DEBUG
        cout << "Figure " << figureId << " has " << numBPoints
            << " boundary points and " << numQuads << " tiles." << endl;
#endif

        // Set up distances vector for this figure
        selectedPartialFigures->figureSize(i, numBPoints, numQuads);
//        int *map_size_ptr = map.getMapSize();
        for (j = 0; j < 4*numQuads; j++)
        {
            double x = quads[j].pnt[0];
            double y = quads[j].pnt[1];
            double z = quads[j].pnt[2];
            double d = (double) map.getDistance(x, y, z);
            selectedPartialFigures->quad_distance(i, j, d);
            if (d > max) max = d;
            else if (d < min) min = d;
        }
        for (j = 0; j < numBPoints; j++)
        {
            double x = bpoints[j].pnt[0];
            double y = bpoints[j].pnt[1];
            double z = bpoints[j].pnt[2];
            double d = (double) map.getDistance(x, y, z);
            selectedPartialFigures->point_distance(i, j, d);
        }
        delete [] (xferList->atomlist);
        delete xferList;
        i++;
    }

    // Set the distance slider's maximum value
#ifdef DEBUG
    cout << "Setting slider max = " << max << "   min = " << min << endl;
#endif
    userInterface->interObjectDistanceSlider->minimum(min);
    userInterface->interObjectDistanceSlider->maximum(max);
    userInterface->interObjectDistanceSlider->step(max - min, 100);
    userInterface->interObjectDistanceSlider->redraw();
#ifdef DEBUG
    cout << "Slider step = " << userInterface->interObjectDistanceSlider->step() << endl;
#endif

    // Set the distance slider's value
    setDistanceSliderFromFigureSlider();

#ifdef DEBUG
    cout << "Drawing figures ";
    for (i = 0; i < selectedPartialFigures->number(); i++)
        cout << selectedPartialFigures->figureID(i) << " ";
    cout << endl;
#endif
}

/*
void P3DUserInterfaceCallback::toggleSelectedAtoms(int button) {
    if (button == 1)
        cout << "Toggled 1" << endl;
    if (button == 2)
        cout << "Toggled 2" << endl;
}
*/

void P3DUserInterfaceCallback::togglePartialRendering(int button) {
    P3DControl * control;

    if (userInterface == NULL)
        return;

    control = userInterface->control;
    if (control == NULL)
        return;

    control->togglePartialRendering((button ? true : false));
    updateSurface();
}

// Called from showPenetrationArea() and setInterObjectFigure()
void P3DUserInterfaceCallback::setDistanceSliderFromFigureSlider() {
    double d, min, max;
    int i, fig;

    // This function assumes the user interface is present and that
    // selectedPartialFigures contains at least one figure.

    // Determine an appropriate value for the distance slider
    fig = (int) userInterface->constraintsFigureSlider->value();
    if (fig < 0) {
        // The slider is set to -1 (selecting all governed figures)
        max = userInterface->interObjectDistanceSlider->maximum();
        min = userInterface->interObjectDistanceSlider->minimum();
        d = (max + min)/2.0;
    }
    else {
        InterfiguralConstraints & ifc = markedFigurePtr->constraints();
        d = -1.0;
        for (i = 0; i < ifc.size(); i++) {
            if (ifc.figure(i) == fig) {
                d = ifc.distance(i);
                break;
            }
        }
        if (d < 0.0) {
            // This figure has no constraint
			Fl_Value_Output * interObjectDistanceValue;

			interObjectDistanceValue = userInterface->interObjectDistanceValue;
			if (interObjectDistanceValue != NULL)
				interObjectDistanceValue->value(-1.0);
			d = -1.0;	// May (intentionally) be less than the minimum
        }
    }

    // Set the distance slider's value
#ifdef DEBUG
    cout << "Reset slider value = " << d << endl;
#endif
    userInterface->interObjectDistanceSlider->value(d);
    userInterface->interObjectDistanceSlider->redraw();
}

// ----------------------------- End of Constraints Code --------------------------------

void P3DUserInterfaceCallback::intensityWindowValuatorChanged()
{
    P3DControl * control;
    Fl_Valuator * minValuator;
    Fl_Valuator * maxValuator;
	Image3D * image;
	int min, max;
	double relativeMin, relativeMax;
	int m;

    if (userInterface == NULL)
        return;

    control = userInterface->control;
    if (control == NULL)
        return;

    image = control->getImagePtr();
    if (image == NULL)
        return;

    minValuator = userInterface->intensityMinValuator;
    maxValuator = userInterface->intensityMaxValuator;

	min = (int) minValuator->value();
	max = (int) maxValuator->value();
	userInterface->intensityMinSlider->value((double) min);
	userInterface->intensityMaxSlider->value((double) max);

	// Clamp the valuator to its range
	m = (int) minValuator->minimum();
	if (min < m) {
		min = m;
		minValuator->value(min);
	}
	else {
		m = (int) minValuator->maximum();
		if (max > m) {
			max = m;
			maxValuator->value(max);
		}
	}

	relativeMin = image->mapActualToRelative(min);
	relativeMax = image->mapActualToRelative(max);

	setIntensityWindow((int) minValuator->value(), relativeMin,
		(int) maxValuator->value(), relativeMax);
}

void P3DUserInterfaceCallback::intensityWindowSliderChanged()
{
    P3DControl * control;
    Fl_Slider * minSlider;
    Fl_Slider * maxSlider;
	Image3D * image;
	int min, max;
	double relativeMin, relativeMax;
	double oldMin, oldMax;

    if (userInterface == NULL)
        return;

    control = userInterface->control;
    if (control == NULL)
        return;

    image = control->getImagePtr();
    if (image == NULL)
        return;

    minSlider = userInterface->intensityMinSlider;
    maxSlider = userInterface->intensityMaxSlider;

	min = (int) minSlider->value();
	max = (int) maxSlider->value();
	userInterface->intensityMinValuator->value((double) min);
	userInterface->intensityMaxValuator->value((double) max);

	relativeMin = image->mapActualToRelative(min);
	relativeMax = image->mapActualToRelative(max);

	// This is needed because FLTK sliders seem to sometimes not take
	// the focus automatically when used
	control->getIntensityWindow(oldMin, oldMax);
	if (oldMin != relativeMin)
		minSlider->take_focus();
	else
		maxSlider->take_focus();

	setIntensityWindow((int) minSlider->value(), relativeMin,
		(int) maxSlider->value(), relativeMax);
}

void P3DUserInterfaceCallback::setIntensityWindow(int min, double relativeMin,
												  int max, double relativeMax)
{
    P3DControl * control;
	Image3D * image;

    control = userInterface->control;
    control->setIntensityWindow(relativeMin, relativeMax);
#ifdef DEBUG
    cout << "Intensity window set to [" << relativeMin << ", " << relativeMax
        << ')' << endl;
#endif
    image = control->getImagePtr();

	if (displayIntensityScale) {
		GreyValue intens = image->mapActualToDisplay(min);
		sprintf(intensityMinStr, "%d", intens);
		intens = image->mapActualToDisplay(max);
		sprintf(intensityMaxStr, "%d", intens);
	}
	else {
		sprintf(intensityMinStr, "%.4lf", relativeMin);
		sprintf(intensityMaxStr, "%.4lf", relativeMax);
	}
	userInterface->intensityMinOutput->label(intensityMinStr);
	userInterface->intensityMaxOutput->label(intensityMaxStr);

    userInterface->modelView->updateImagePlanes(false);
}

// Called when the Edit -> "Attach a Subfigure" menu item is selected
void P3DUserInterfaceCallback::initializeAttachSubfigureDialog()
{
    if(userInterface == NULL)
        return;

    P3DControl * control;
    control = userInterface->control;
    if(control == NULL)
        return;

    Fl_Window * dialog = userInterface->attachSubfigureDialog;
    if(dialog == NULL)
        return;

    Fl_Button * step = userInterface->attachStepTextOutput;
    if(step == NULL)
        return;

    step->clear_visible_focus();
    step->label("Step 1");
	attachmentMode = getAttachSubfigureMode();
	userInterface->attachInstructionsTextOutput->hide();
	userInterface->attachModeGroup->show();

    deselectAll();

	userInterface->attachNextButton->label("Next");
    dialog->redraw();
    subfigureAttachScreenId = 1;	// The step number
}

// This function is called when a button of the first step of the
// Attach Subfigure Dialog is pressed.
void P3DUserInterfaceCallback::setAttachSubfigureMode(int mode)
{
    if(userInterface == NULL)
        return;

    if (! userInterface->attachProtrusionButton->active_r())
	    return;

	// mode of 1 ==> PROTRUDE
	// mode of 2 ==> INDENT
	attachmentMode = (M3DFigureTreeNode::SubfigureAttachment_t) mode;
}

M3DFigureTreeNode::SubfigureAttachment_t P3DUserInterfaceCallback::getAttachSubfigureMode()
{
    if(userInterface == NULL)
        return M3DFigureTreeNode::UNATTACHED;	// Actually an error report

    if (! userInterface->attachProtrusionButton->active_r())
	    return M3DFigureTreeNode::UNATTACHED;

	if (1 == (int) userInterface->attachProtrusionButton->value())
		return M3DFigureTreeNode::PROTRUDE;
	else
		return M3DFigureTreeNode::INDENT;
}

/*  Thus function does not permit moving a figure from one place in
    a tree to another place in the same tree (because of the test
    in case 2, below).  Note that control->attachSubfigure() does
    permit such an action.
*/
void P3DUserInterfaceCallback::advanceAttachSubfigureDialog()
{
    if(userInterface == NULL)
        return;

    if (! userInterface->attachNextButton->active_r())
	    return;

    P3DControl * control = userInterface->control;
    if(control == NULL)
        return;

    M3DObject * object = control->getObjectPtr();
    Fl_Window * dialog = userInterface->attachSubfigureDialog;
    Fl_Button * step = userInterface->attachStepTextOutput;
    Fl_Output * instr = userInterface->attachInstructionsTextOutput;
    if(object == NULL || dialog == NULL || step == NULL || instr == NULL)
        return;

    int primId, num;
	M3DQuadFigure * figure;
    P3DView * modelView;
	bool error;

    switch(subfigureAttachScreenId)	// This is the step being left
    {
        case 1:
			instr->show();
			userInterface->attachModeGroup->hide();

            step->label("Step 2");
            instr->value(ATTACH_STEP2);

			control->setEditingMode(NO_EDITING_MODE);
			control->allowSelection(false);
            break;

        case 2:
            parentFigureId = -1;
            control->getRelativeMarkedPrimitiveId(parentFigureId, primId);
            if(parentFigureId == -1)
            {
                fl_alert("Please mark a figure for the parent.");
                subfigureAttachScreenId--;
                break;
            }

            step->label("Step 3");
            instr->value(ATTACH_STEP3);

			control->setEditingMode(NO_EDITING_MODE);
			control->allowSelection(false);
            break;

        case 3:
            childFigureId = -1;
            control->getRelativeMarkedPrimitiveId(childFigureId, primId);
            if(childFigureId == -1 || childFigureId == parentFigureId)
            {
                fl_alert("Please mark a different figure for the child.");
                subfigureAttachScreenId--;
                break;
            }
            if (object->testForPotentialLoop(parentFigureId, childFigureId) == true) {
				// Note: the loop test will fail if the loop only contains the parent
				// and child.  The later test in object->attachFigureTreeNode() reached
				// through the call to control->attachSubfigure() below will catch this case.
                fl_alert("This would generate a loop; please mark a different child figure.");
                subfigureAttachScreenId--;
                break;
            }

            step->label("Step 4");
            instr->value(ATTACH_STEP4);
			userInterface->attachNextButton->label("Finish");

            control->setEditingMode(ANY_EDITING_MODE);	// Needed for row/column selection in stage 4
			control->allowSelection(true);	// Needed for row/column selection in stage 4 and next line
			deselectAll();
            break;

        case 4:

            figure = dynamic_cast<M3DQuadFigure*>(object->getFigurePtr(childFigureId));
			num = figure->numberSelected();
            if (num == 0) {
                fl_alert("Please select atoms on the child figure.");
				deselectAll();
                subfigureAttachScreenId--;
                break;
            }

			error = false;
            if (figure->rowSelected() < 0) {
				if (figure->columnSelected() < 0)
					error = true;
				else {
					if (num != figure->getRowCount())
						error = true;
				}
			}
			else {
				if (num != figure->getColumnCount())
					error = true;
			}

            if (error) {
                fl_alert("Please select a full row or column on the child figure.");
				deselectAll();
                subfigureAttachScreenId--;
                break;
            }

			// Delete all image landmarks.  AGG: This is not optimal!
			modelView = userInterface->modelView;
			finishEditModel();
			if(modelView != NULL) {
				control->deleteImageLandmarks();
				setLandmarksImageNumber(0);
			}
			finishEditLandmarks();
			editingLandmark = -1;

            if (control->attachSubfigure(parentFigureId, childFigureId, attachmentMode)) {
				updateSurface();
				updateVisibilityFiguresDisplay(false);
			}

			initializeAttachSubfigureDialog();	// Reset the attachment procedure
			subfigureAttachScreenId--;
            break;

        default:
            break;
    };

    dialog->redraw();
    subfigureAttachScreenId++;
}

void P3DUserInterfaceCallback::cancelAttachSubfigureDialog()
{
    if(userInterface == NULL)
        return;

    if (! userInterface->attachCancelButton->active_r())
	    return;

    P3DControl * control = userInterface->control;
    if(control == NULL)
        return;

    control->setEditingMode(ANY_EDITING_MODE);
    control->allowSelection(true);

    initializeAttachSubfigureDialog();	// Reset the attachment procedure
}


void P3DUserInterfaceCallback::detachSubfigure()
{
    int primId;
	int detachFigureId;

    if(userInterface == NULL)
        return;

    P3DControl * control = userInterface->control;
    if(control == NULL)
        return;

    M3DObject * object = control->getObjectPtr();
    if(object == NULL)
        return;

	detachFigureId = -1;
	control->getRelativeMarkedPrimitiveId(detachFigureId, primId);
	if (detachFigureId == -1)
	{
		fl_alert("Please mark a primitive on the subfigure to be detached.");
		return;
	}
	if (object->isRootNode(detachFigureId))
	{
		fl_alert("Please mark a subfigure.");
		return;
	}

	if (control->detachSubfigure(detachFigureId)) {
		updateObject();
		updateVisibilityFiguresDisplay(false);
		deselectAll();
		userInterface->modelView->updateObject();
		M3DFigure * figure = object->getFigurePtr(detachFigureId);
		figure->select();
	}
}

void P3DUserInterfaceCallback::updateElongationWindow()
{
    P3DControl * control;
    M3DPrimitive * markedPrimitive;
	int figureId, atomId;

    if(userInterface == NULL)
        return;

    control = userInterface->control;
    if(control == NULL)
        return;

    markedPrimitive = control->getMarkedPrimitive();

    if(markedPrimitive == NULL)
    {
		userInterface->elongationFigureBtn->hide();
		userInterface->elongationFigureNumberOutput->hide();
		userInterface->elongation1stMessageBtn->show();
		userInterface->elongation1stMessageBtn->clear_visible_focus();
		userInterface->elongation2ndMessageBtn->show();
		userInterface->elongation2ndMessageBtn->clear_visible_focus();
		userInterface->elongateRowSlider->deactivate();
		userInterface->elongateColumnSlider->deactivate();
    }
    else
    {
		userInterface->elongationFigureBtn->show();
		userInterface->elongationFigureBtn->clear_visible_focus();
        control->getRelativeMarkedPrimitiveId(figureId, atomId);
		userInterface->elongationFigureNumberOutput->value((double) figureId);
		userInterface->elongationFigureNumberOutput->show();
		userInterface->elongation1stMessageBtn->hide();
		userInterface->elongation2ndMessageBtn->hide();
		userInterface->elongateRowSlider->activate();
		userInterface->elongateColumnSlider->activate();
    }

    userInterface->elongationDialog->redraw();
}

void P3DUserInterfaceCallback::elongateRow()
{
    if(userInterface == NULL)
        return;

    Fl_Slider * slider = userInterface->elongateRowSlider;
    if (! slider->active_r())
	    return;

    P3DControl * control = userInterface->control;
    if(control == NULL)
        return;

    control->elongateRow(slider->value());
    update();
}

void P3DUserInterfaceCallback::elongateColumn()
{
    if(userInterface == NULL)
        return;

    Fl_Slider * slider = userInterface->elongateColumnSlider;
    if (! slider->active_r())
	    return;

    P3DControl * control = userInterface->control;
    if(control == NULL)
        return;

    control->elongateColumn(slider->value());
    update();
}

typedef void (*voidFn)(void*);

void buttonReleased(void * arg) {		// Non-class function used below
	struct sliderState {
		void * widget;
		void (* pushAction)(void *, bool);
		void (* releaseAction)(void *, bool);
		void * self;
		bool rowCol;
		bool pushed;
	} * s;

    s = (sliderState *) arg;
	if (Fl::pushed() == s->widget) {
		if (s->pushed)
			return;
		// The slider button was pushed
		s->pushed = true;
		if (s->pushAction)
			s->pushAction(s->self, s->rowCol);
		return;
	}
	if (s->pushed) {
		// The slider button was released
		s->pushed = false;
		if (s->releaseAction)
			s->releaseAction(s->self, s->rowCol);
		return;
	}
}

void P3DUserInterfaceCallback::elongateFinish()
{
    if(userInterface == NULL)
        return;

	Fl::remove_check((voidFn) buttonReleased, (void *) &rowSliderState);
	Fl::remove_check((voidFn) buttonReleased, (void *) &columnSliderState);

	userInterface->elongationDialog->hide();
}

void P3DUserInterfaceCallback::releaseDispatcher(void * self, bool rowCol)
{
	if (rowCol)
		((P3DUserInterfaceCallback *) self)->resetRowElongation();
	else
		((P3DUserInterfaceCallback *) self)->resetColumnElongation();
}

void P3DUserInterfaceCallback::resetRowElongation()
{
	// Released row slider button
    Fl_Slider * slider = userInterface->elongateRowSlider;
	slider->value(1.0);
}

void P3DUserInterfaceCallback::resetColumnElongation()
{
	// Released column slider button
    Fl_Slider * slider = userInterface->elongateColumnSlider;
	slider->value(1.0);
}

// Called when the Edit -> "Elongation..." menu item is selected
void P3DUserInterfaceCallback::initializeElongater()
{
    if(userInterface == NULL)
        return;

    P3DControl * control = userInterface->control;
    if(control == NULL)
        return;

	// Set up to get callbacks to resetRowElongation() and resetColumnElongation(),
	// when the slider buttons are released.  This allows the slider positions to
	// be reset as soon as the buttons are released.
	rowSliderState.pushed = false;
	rowSliderState.widget = (void *) userInterface->elongateRowSlider;
	rowSliderState.pushAction = NULL;
	rowSliderState.self = this;
	rowSliderState.rowCol = true;
	rowSliderState.releaseAction = releaseDispatcher;
	Fl::add_check((voidFn) buttonReleased, (void *) &rowSliderState);

	columnSliderState.pushed = false;
	columnSliderState.widget = (void *) userInterface->elongateColumnSlider;
	columnSliderState.pushAction = NULL;
	columnSliderState.self = this;
	columnSliderState.rowCol = false;
	columnSliderState.releaseAction = releaseDispatcher;
	Fl::add_check((voidFn) buttonReleased, (void *) &columnSliderState);

	userInterface->elongation1stMessageBtn->label(ELONGATION_STEP1);
	updateElongationWindow();
	userInterface->elongationDialog->show();
}

void P3DUserInterfaceCallback::gaussBlur()
{
    P3DControl * control;

    if(userInterface == NULL)
        return;

    control = userInterface->control;
    if(control == NULL)
        return;

	if (control->haveImage() == false) {
		cout << "An image must have been loaded before it can be blurred\n";
		fl_beep(FL_BEEP_ERROR);
		return;
	}

    control->gaussBlur();
    userInterface->modelView->updateImagePlanes(true);
}

#ifdef BINARY

void P3DUserInterfaceCallback::binaryPabloAlignModel()
{
	P3DControl* control = userInterface->control;
    if (control == NULL)
        return;

	// Use Method of Moments to initialize object to image;
	//   aligns centroid (position) ONLY
	MatchUtility::initializeObjectToImage(control->getObjectPtr(), control->getImagePtr(), 1, 1);
	userInterface->modelView->updateObject();	// Show object in new spot!
}

void P3DUserInterfaceCallback::distMapAs()
{
    const char * filename;

    filename = askSingleFilename(ImageDirectory, "Save Distance Map As", "*.raw3");
	if (filename != NULL)
		distMap(filename);
	else
		distMap("ddm.raw3");
}

// Asks control to create a distance map
void P3DUserInterfaceCallback::distMap(const char *filename)
{
    P3DControl * control;

    if(userInterface == NULL)
        return;

    control = userInterface->control;
    if(control == NULL)
        return;

	if (control->haveImage() == false) {
		cout << "An image must have been loaded before a distance map can be calculated\n";
		fl_beep(FL_BEEP_ERROR);
		return;
	}

    control->distMap(filename);
    userInterface->modelView->updateImagePlanes(true);
}

#else	 /* BINARY */

void P3DUserInterfaceCallback::binaryPabloAlignModel()
{
}

void P3DUserInterfaceCallback::distMapAs()
{
}

void P3DUserInterfaceCallback::distMap(const char *filename)
{
}

#endif	 /* BINARY */

void P3DUserInterfaceCallback::outputHistogram()
{
    P3DControl * control;
    const char * filename;

    if(userInterface == NULL)
        return;

    control = userInterface->control;
    if(control == NULL)
        return;

    filename = askSingleFilename(ModelDirectory, "Output Histogram", "*.hist");
	if (filename != NULL)
		control->outputHistogram(filename);
}


void P3DUserInterfaceCallback::makeMovie()
{
    Fl_Menu_Item * screenDumpItem;

    if(userInterface == NULL)
        return;

    screenDumpItem = userInterface->makeMovieToggle;

    if(screenDumpItem->value() == 0)
        userInterface->modelView->stopScreenDump();
    else
        userInterface->modelView->startScreenDump();
}

void P3DUserInterfaceCallback::screenDump()
{
    if(userInterface == NULL)
        return;

    userInterface->modelView->scheduleScreenDump();
	redraw();	// Necessary to avoid having the Tools Menu visible in the dump
}

void P3DUserInterfaceCallback::startModelTileSlideShow()
{
    userInterface->control->startModelTileSlideShow();
    changeModelSlideShow();
    update();
}


/////////////////  MODEL SLIDE SHOW


void P3DUserInterfaceCallback::addModelsToSlideShow()
{
    Fl_File_Browser * files;
    int count, i;
    char relative[1024];
    P3DControl * control;
	char * str;
	char * dir;


    control = userInterface->control;
    if(control == NULL)
        return;

	dir = globalControl->readString(ModelDirectory);
	if (dir[0] == '\0')
		dir = globalControl->readString(ImageDirectory);
	str = new char[strlen(dir) + 2];
	strcpy(str, dir);
	strcat(str, "/");

    Fl_File_Chooser fc(str, "*.m3d", Fl_File_Chooser::MULTI, "Add Models To Slide Show");
	fc.callback(NULL);

    fc.show();
    while (fc.shown())
        Fl::wait();
	delete [] str;

	const char * filename = fc.value(1);
	if (fc.count() == 0)
		return;		// User canceled request

	// Adjust ModelDirectory preference as needed
	str = (char *) fc.directory();
	dir = globalControl->readString(ModelDirectory);
	if (0 != strcmp(dir, str))
		globalControl->write(ModelDirectory, fc.directory());

    files = userInterface->modelSlideShowFileBrowser;

    count = fc.count();
    if(count > 0)
    {
        for(i = 1; i <= count; i ++)
        {
            fl_filename_relative(relative, sizeof(relative), fc.value(i));

            files->add(relative,
                Fl_File_Icon::find(fc.value(i), Fl_File_Icon::PLAIN));

            control->addModelSlideShow(relative);
        }

        files->redraw();
    }

    userInterface->modelSlideShowCounter->bounds(1.0, 
        userInterface->modelSlideShowCounter->maximum() + (double) count);
    userInterface->modelSlideShowCounter->value(1.0);
    changeModelSlideShow();
}

void P3DUserInterfaceCallback::removeModelSlideShow()
{
    P3DControl * control;
    Fl_Counter * counter;
    Fl_File_Browser * files;
    int modelNumber;


    control = userInterface->control;
    counter = userInterface->modelSlideShowCounter;
    files = userInterface->modelSlideShowFileBrowser;

    // Nothing to remove
    if(counter->value() == 0.0)
        return;

    modelNumber = (int) counter->value() - 1;
    control->removeModelSlideShow(modelNumber);
    files->remove(modelNumber + 1);

    if(modelNumber == counter->maximum() - 1)
        counter->value(modelNumber);

    if(counter->maximum() == 1.0)
    {
        counter->bounds(0.0, 0.0);
        counter->value(0.0);
    }
    else
    {
        counter->bounds(1.0, counter->maximum() - 1.0);
        changeModelSlideShow();
    }
}

void P3DUserInterfaceCallback::moveUpModelSlideShow()
{
    P3DControl * control;
    Fl_Counter * counter;
    Fl_File_Browser * files;
    int modelNumber;


    control = userInterface->control;
    counter = userInterface->modelSlideShowCounter;
    files = userInterface->modelSlideShowFileBrowser;

    // Can't move
    if(counter->value() <= 1.0)
        return;

    modelNumber = (int) counter->value() - 1;
    files->move(modelNumber + 1, modelNumber);
    files->select(modelNumber);
    control->moveUpModelSlideShow(modelNumber);
    counter->value(counter->value() - 1.0);

}

void P3DUserInterfaceCallback::moveDownModelSlideShow()
{
    P3DControl * control;
    Fl_Counter * counter;
    Fl_File_Browser * files;
    int modelNumber;


    control = userInterface->control;
    counter = userInterface->modelSlideShowCounter;
    files = userInterface->modelSlideShowFileBrowser;

    // Can't move
    if(counter->value() == 0.0 || counter->value() == counter->maximum())
        return;

    modelNumber = (int) counter->value() - 1;
    files->move(modelNumber + 1, modelNumber + 2);
    files->select(modelNumber + 2);
    control->moveDownModelSlideShow(modelNumber);
    counter->value(counter->value() + 1.0);

}

void P3DUserInterfaceCallback::clearModelSlideShow()
{
    userInterface->control->clearModelSlideShow();
    userInterface->modelSlideShowFileBrowser->clear();
    userInterface->modelSlideShowCounter->bounds(0.0, 0.0);
    userInterface->modelSlideShowCounter->value(0.0);
    update();
}

void P3DUserInterfaceCallback::changeModelSlideShow()
{
    P3DControl * control;
    Fl_Counter * counter;
    Fl_File_Browser * files;
    int modelNumber;


    control = userInterface->control;
    counter = userInterface->modelSlideShowCounter;
    files = userInterface->modelSlideShowFileBrowser;

    modelNumber = (int) counter->value() - 1;
    control->changeModelSlideShow(modelNumber);
    files->select(modelNumber + 1);

	if (int(userInterface->slideShowCenter->value()))
	{
	    userInterface->modelView->setCenterOfRotation();
	}

	if (int(userInterface->slideShowSync->value()))
	{
        userInterface->tileSlideShowCounter->value(counter->value());
		changeTileSlideShow();	// side-effect: calls update()
	}
	update();
}

void P3DUserInterfaceCallback::selectModelSlideShow()
{
    P3DControl * control;
    Fl_Counter * counter;
    Fl_File_Browser * files;
    int modelNumber;


    control = userInterface->control;
    counter = userInterface->modelSlideShowCounter;
    files = userInterface->modelSlideShowFileBrowser;

    modelNumber = (int) files->value() - 1;
    if(modelNumber == -1)
    {
        files->select((int) counter->value());
        return;
    }

    counter->value(modelNumber + 1.0);
    control->changeModelSlideShow(modelNumber);
    update();
}

void P3DUserInterfaceCallback::quitModelSlideShow()
{
    userInterface->control->quitModelSlideShow();
    update();
}


/////////////////  TILE SLIDE SHOW


void P3DUserInterfaceCallback::addTilesToSlideShow()
{
    Fl_File_Browser * files;
    int count, i;
    char relative[1024];
    P3DControl * control;
	char * str;
	char * dir;


    control = userInterface->control;
    if(control == NULL)
        return;

	dir = globalControl->readString(TileSetDirectory);
	if (dir[0] == '\0')
		dir = globalControl->readString(ImageDirectory);
	str = new char[strlen(dir) + 2];
	strcpy(str, dir);
	strcat(str, "/");

    Fl_File_Chooser fc(str, "*.byu", Fl_File_Chooser::MULTI, "Add Tiles To Slide Show");
	fc.callback(NULL);

    fc.show();
    while (fc.shown())
        Fl::wait();
	delete [] str;

	const char * filename = fc.value(1);
	if (fc.count() == 0)
		return;		// User canceled request

	// Adjust TileDirectory preference as needed
	str = (char *) fc.directory();
	dir = globalControl->readString(TileSetDirectory);
	if (0 != strcmp(dir, str))
		globalControl->write(TileSetDirectory, fc.directory());

    files = userInterface->tileSlideShowFileBrowser;

    count = fc.count();
    if(count > 0)
    {
        for(i = 1; i <= count; i ++)
        {
            fl_filename_relative(relative, sizeof(relative), fc.value(i));

            files->add(relative,
                Fl_File_Icon::find(fc.value(i), Fl_File_Icon::PLAIN));

            control->addTileSlideShow(relative);
			cout << "Loading tiles: " << relative << endl;
        }

        files->redraw();
    }

    userInterface->tileSlideShowCounter->bounds(1.0, 
        userInterface->tileSlideShowCounter->maximum() + (double) count);
    userInterface->tileSlideShowCounter->value(1.0);
    changeTileSlideShow();
}

void P3DUserInterfaceCallback::removeTileSlideShow()
{
    P3DControl * control;
    Fl_Counter * counter;
    Fl_File_Browser * files;
    int TileNumber;


    control = userInterface->control;
    counter = userInterface->tileSlideShowCounter;
    files = userInterface->tileSlideShowFileBrowser;

    // Nothing to remove
    if(counter->value() == 0.0)
        return;

    TileNumber = (int) counter->value() - 1;
    control->removeTileSlideShow(TileNumber);
    files->remove(TileNumber + 1);

    if(TileNumber == counter->maximum() - 1)
        counter->value(TileNumber);

    if(counter->maximum() == 1.0)
    {
        counter->bounds(0.0, 0.0);
        counter->value(0.0);
    }
    else
    {
        counter->bounds(1.0, counter->maximum() - 1.0);
        changeTileSlideShow();
    }
}

void P3DUserInterfaceCallback::moveUpTileSlideShow()
{
    P3DControl * control;
    Fl_Counter * counter;
    Fl_File_Browser * files;
    int TileNumber;


    control = userInterface->control;
    counter = userInterface->tileSlideShowCounter;
    files = userInterface->tileSlideShowFileBrowser;

    // Can't move
    if(counter->value() <= 1.0)
        return;

    TileNumber = (int) counter->value() - 1;
    files->move(TileNumber + 1, TileNumber);
    files->select(TileNumber);
    control->moveUpTileSlideShow(TileNumber);
    counter->value(counter->value() - 1.0);
}

void P3DUserInterfaceCallback::moveDownTileSlideShow()
{
    P3DControl * control;
    Fl_Counter * counter;
    Fl_File_Browser * files;
    int TileNumber;


    control = userInterface->control;
    counter = userInterface->tileSlideShowCounter;
    files = userInterface->tileSlideShowFileBrowser;

    // Can't move
    if(counter->value() == 0.0 || counter->value() == counter->maximum())
        return;

    TileNumber = (int) counter->value() - 1;
    files->move(TileNumber + 1, TileNumber + 2);
    files->select(TileNumber + 2);
    control->moveDownTileSlideShow(TileNumber);
    counter->value(counter->value() + 1.0);

}

void P3DUserInterfaceCallback::clearTileSlideShow()
{
    userInterface->control->clearTileSlideShow();
    userInterface->tileSlideShowFileBrowser->clear();
    userInterface->tileSlideShowCounter->bounds(0.0, 0.0);
    userInterface->tileSlideShowCounter->value(0.0);
	userInterface->modelView->setTileFileChanged();
    update();
}

void P3DUserInterfaceCallback::changeTileSlideShow()
{
    P3DControl * control;
    Fl_Counter * counter;
    Fl_File_Browser * files;
    int TileNumber;


    control = userInterface->control;
    counter = userInterface->tileSlideShowCounter;
    files = userInterface->tileSlideShowFileBrowser;


    TileNumber = (int) counter->value() - 1;
    control->changeTileSlideShow(TileNumber);
    files->select(TileNumber + 1);
	userInterface->modelView->setTileFileChanged();
	if (int(userInterface->slideShowCenter->value()))
		userInterface->modelView->setCenterOfRotation();
    update();
}

void P3DUserInterfaceCallback::selectTileSlideShow()
{
    P3DControl * control;
    Fl_Counter * counter;
    Fl_File_Browser * files;
    int TileNumber;


    control = userInterface->control;
    counter = userInterface->tileSlideShowCounter;
    files = userInterface->tileSlideShowFileBrowser;

    TileNumber = (int) files->value() - 1;
    if(TileNumber == -1)
    {
        files->select((int) counter->value());
        return;
    }

    counter->value(TileNumber + 1.0);
    control->changeTileSlideShow(TileNumber);
	userInterface->modelView->setTileFileChanged();
    update();
}

void P3DUserInterfaceCallback::quitTileSlideShow()
{
    userInterface->control->quitTileSlideShow();
	userInterface->modelView->setTileFileChanged();
    update();
}


////////////////////


void P3DUserInterfaceCallback::outputPOVRay()
{
    P3DControl * control;
    M3DObject * object;
    M3DtoPovray povrayOut;
    char * filename;

    control = userInterface->control;
    if(control == NULL)
        return;

    object = control->getObjectPtr();
    if(object == NULL)
        return;

    filename = fl_file_chooser("Export to POV-Ray", "*.pov", ".");

    povrayOut.setObject(object);
    povrayOut.output(filename);
}

// ---------------------------  Start of CPNS Code  ---------------------------

void P3DUserInterfaceCallback::cpnsZero()
{
	resetCPNSDeformation() ;
	// userInterface->cpnsZeroBtn->activate() ;
	cpnsDeform() ;
}

// For Adaptive Pablo, this zeros the PGA sliders and deactivates the Zero
// Button.  For Regular Pablo, the mean model is also initialized.
void P3DUserInterfaceCallback::resetCPNSDeformation()
{
    userInterface->cpnsSlider1->value(0);
    userInterface->cpnsSlider2->value(0);
    userInterface->cpnsSlider3->value(0);
    userInterface->cpnsSlider4->value(0);
    userInterface->cpnsSlider5->value(0);
	// userInterface->cpnsZeroBtn->deactivate();
}

void P3DUserInterfaceCallback::initCPNSDeform()
{
	resetCPNSDeformation();

	P3DControl * control = userInterface->control;

	if (control == NULL)
		return;
	if (! control->haveCPNS() )
		return;

	control->resetCPNSStage();
}

void P3DUserInterfaceCallback::cpnsDeform() {

    bool worked;

    if (userInterface == NULL)
        return;

    P3DControl* control = userInterface->control;
    if (control == NULL)
        return;

	if ( ! control->haveCPNS() )
		return;

    std::vector<double> vals;

    vals.push_back(userInterface->cpnsSlider1->value());
    vals.push_back(userInterface->cpnsSlider2->value());
    vals.push_back(userInterface->cpnsSlider3->value());
    vals.push_back(userInterface->cpnsSlider4->value());
    vals.push_back(userInterface->cpnsSlider5->value());

	cout << "Coeffs : "
		 << setw(4) << setprecision(2) << userInterface->cpnsSlider1->value()
		 << setw(4) << setprecision(2) << userInterface->cpnsSlider2->value()
		 << setw(4) << setprecision(2) << userInterface->cpnsSlider3->value()
		 << setw(4) << setprecision(2) << userInterface->cpnsSlider4->value()
		 << setw(4) << setprecision(2) << userInterface->cpnsSlider5->value()
		 << endl;

    worked = control->cpnsDeform(vals);
    if (! worked) {
        fl_alert("Cannot deform the current model.\nSee \'Windows->Optimizer Settings\'");
        resetCPNSDeformation();
    }

	vals.clear() ;

    update();
}

// Called when the CPNS Deformation Window is opened
void P3DUserInterfaceCallback::cpnsInitDeformationWindow()
{
}

void P3DUserInterfaceCallback::cpnsDeformationWindowDone()
{

    if(userInterface == NULL)
        return;

	disableCPNSDeformationWindow();

}

// ---------------------------  End of CPNS Code  ---------------------------

// ---------------------------  Start of PGA Code  ---------------------------

// For Adaptive Pablo, this zeros the PGA sliders and deactivates the Zero
// Button.  For Regular Pablo, the mean model is also initialized.
void P3DUserInterfaceCallback::resetPGADeformation()
{
#ifdef DEBUG
	cout << "P3DUserInterfaceCallback::resetPGADeformation()" << endl;
#endif

    userInterface->pgaSlider1->value(0);
    userInterface->pgaSlider2->value(0);
    userInterface->pgaSlider3->value(0);
    userInterface->pgaSlider4->value(0);
    userInterface->pgaSlider5->value(0);
	userInterface->pgaZeroBtn->deactivate();
}

#ifndef BINARY

// Called when the Zero Button in the PGA Deformation Window is pressed
void P3DUserInterfaceCallback::pgaZero()
{
#ifdef DEBUG
	cout << "P3DUserInterfaceCallback::pgaZero()" << endl;
#endif
	resetPGADeformation();
	userInterface->pgaZeroBtn->activate();
	(void) pgaDeform();
}

#else	/* BINARY */

void P3DUserInterfaceCallback::pgaZero()
{

}

#endif	/* BINARY */

void P3DUserInterfaceCallback::initPGADeform()
{
#ifdef DEBUG
	cout << "P3DUserInterfaceCallback::initPGADeform()" << endl;
#endif
	resetPGADeformation();

	P3DControl * control = userInterface->control;
#ifndef BINARY
	Fl_Button * pgaText = userInterface->pgaStageOutput;
	Fl_Button * pgaMeanResBtn = userInterface->pgaMeanResBtn;
#endif

	if (control == NULL)
		return;
	if (! control->havePGA())
		return;
	control->resetResiduePGAStage();

#ifndef BINARY
	pgaMeanResBtn->activate();
	deactivatePGASliders();

	pgaText->label(control->getResiduePGAStageName());
	pgaText->redraw();
#endif
}

#ifndef BINARY

// This function is called when the "Apply Mean of Residue" button
// on the PGA Deformation Window is pressed.
void P3DUserInterfaceCallback::pgaApplyMeanRes()
{
	bool worked;

#ifdef DEBUG
	cout << "P3DUserInterfaceCallback::pgaApplyMeanRes()" << endl;
#endif
	Fl_Button * pgaNextBtn = userInterface->pgaNextBtn;
	Fl_Button * pgaMeanResBtn = userInterface->pgaMeanResBtn;

	//Apply mean of the residue to day 1 model fit
	//At this point, object and intermedObj, latestObj in M3DPGAStats class are out of sync!
	P3DControl * control = userInterface->control;

	if (control == NULL)
		return;

	if (! control->havePGA()){
		fl_alert("PGA is not set.\nSee \'Windows->Optimizer Settings\'");	
		return;
	}

	worked = control->applyMeanResidue();
    if (! worked)
    {
        fl_alert("PGA file does not match the current model.\nSee \'Windows->Optimizer Settings\'");
        resetPGADeformation();
    }

	pgaMeanResBtn->deactivate();
	pgaNextBtn->activate();
	activatePGASliders();
	userInterface->pgaZeroBtn->activate();

	update();
}

#else	/* BINARY */

void P3DUserInterfaceCallback::pgaApplyMeanRes()
{
}

#endif	/* BINARY */

void P3DUserInterfaceCallback::pgaDeform()
{
    bool worked;

#ifdef DEBUG
	cout << "P3DUserInterfaceCallback::pgaDeform()" << endl;
#endif
    if (userInterface == NULL)
        return;

    P3DControl* control = userInterface->control;
    if (control == NULL)
        return;

	if (! control->havePGA())
		return;

    std::vector<double> vals;
    vals.push_back(userInterface->pgaSlider1->value());
    vals.push_back(userInterface->pgaSlider2->value());
    vals.push_back(userInterface->pgaSlider3->value());
    vals.push_back(userInterface->pgaSlider4->value());
    vals.push_back(userInterface->pgaSlider5->value());

	cout << "Coeffs : "
		 << setw(4) << setprecision(2) << userInterface->pgaSlider1->value()
		 << setw(4) << setprecision(2) << userInterface->pgaSlider2->value()
		 << setw(4) << setprecision(2) << userInterface->pgaSlider3->value()
		 << setw(4) << setprecision(2) << userInterface->pgaSlider4->value()
		 << setw(4) << setprecision(2) << userInterface->pgaSlider5->value()
		 << endl;

    worked = control->pgaDeform(vals);
    if (! worked)
    {
        fl_alert("Cannot deform the current model.\nSee \'Windows->Optimizer Settings\'");
        resetPGADeformation();
#ifndef BINARY
		userInterface->pgaZeroBtn->activate();
#endif
    }

    update();
}

/*
	1. get the current depth
	2. If  depth < set.Size
		1. get the value from control
		2. show the window/sliders, as the slider changes, pass the value to
		   control and deform the model 
	   Else go back to the 0 depth and reset deformation.
 */
void P3DUserInterfaceCallback::pgaNextScale()
{
#ifndef BINARY
	P3DControl * control;
	Fl_Button * pgaText;

	int depth;

    if (userInterface == NULL)
        return;

    control = userInterface->control;
    if (control == NULL)
        return;

	if (! control->havePGA())
		return;

	// Print out the current stage
    depth = control->getResiduePGAStage();

#ifdef DEBUG
	cout << "Leaving stage " << depth << endl;
#endif

    // Get the next stage we are entering
	control->setResiduePGAStage(depth+1);
	pgaText = userInterface->pgaStageOutput;
	pgaText->label(control->getResiduePGAStageName());
	pgaText->redraw();

#ifdef DEBUG
	cout << "Entering stage " << control->getResiduePGAStage() << endl;
#endif

    // Finally initialize and display user interface for next stage
	resetPGADeformation();

	Fl_Button * pgaNextBtn = userInterface->pgaNextBtn;
	Fl_Button * pgaMeanResBtn = userInterface->pgaMeanResBtn;

	pgaMeanResBtn->activate();
	deactivatePGASliders();

	update();
#endif	/* BINARY */
}

// Called when the PGA Deformation Window is opened
void P3DUserInterfaceCallback::pgaInitDeformationWindow()
{
#ifdef DEBUG
	cout << "P3DUserInterfaceCallback::pgaInitDeformationWindow()" << endl;
#endif
#ifndef BINARY
    if(userInterface == NULL)
        return;

	userInterface->pgaStageOutput->label(NULL);
	userInterface->pgaStageOutput->clear_visible_focus();
	userInterface->pgaStageNameLabel->clear_visible_focus();

	disableDeformationWindow();
	if (disableOptimizerWindow())
		initPGADeform();
#endif
}

void P3DUserInterfaceCallback::pgaDeformationWindowDone()
{
#ifdef DEBUG
	cout << "P3DUserInterfaceCallback::pgaCloseDeformationWindow()" << endl;
#endif
    if(userInterface == NULL)
        return;

	disableDeformationWindow();
#ifndef BINARY
	enableOptimizerWindow();
#endif
}

// ---------------------------  End of PGA Code  ---------------------------

// ---------------------------  Start PCA Code  ----------------------------

#ifdef PCA_VISUALIZER

void P3DUserInterfaceCallback::initPCADeform()
{
	P3DControl * control;

#ifdef DEBUG
	cout << "P3DUserInterfaceCallback::initPCADeform()" << endl;
#endif
	userInterface->pcaSlider1->value(0);
    userInterface->pcaSlider2->value(0);
    userInterface->pcaSlider3->value(0);
    userInterface->pcaSlider4->value(0);
    userInterface->pcaSlider5->value(0);

	control = userInterface->control;
    if (control == NULL)
        return;

	if (control->initPCA_Deform())
		userInterface->pcaDialog->show();
	else
		fl_alert("A PCA mean object and a PCA residue file must be loaded before deforming");
}

void P3DUserInterfaceCallback::loadPCAFile()
{
	P3DControl * control;
    const char * filename;

    control = userInterface->control;
    if (control == NULL)
        return;

    filename = askSingleFilename(ModelDirectory, "Load PCA File", "*.pca", false);
	if (filename != NULL)
    {
        if (! control->loadPCA_File(filename))
            fl_alert("File \'%s\' does not contain PCA statistics.", filename);
    }
}

void P3DUserInterfaceCallback::pcaDeform()
{
	P3DControl * control;

#ifdef DEBUG
	cout << "P3DUserInterfaceCallback::pcaDeform()" << endl;
#endif
    if (userInterface == NULL)
        return;

    control = userInterface->control;
    if (control == NULL)
        return;

    std::vector<float> vals;
    vals.push_back(userInterface->pcaSlider1->value());
    vals.push_back(userInterface->pcaSlider2->value());
    vals.push_back(userInterface->pcaSlider3->value());
    vals.push_back(userInterface->pcaSlider4->value());
    vals.push_back(userInterface->pcaSlider5->value());

    if (! control->pcaDeform(vals))
    {
        fl_alert("Cannot deform the current model.\nSee \'Windows->Optimizer Settings\'");
        initPCADeform();
    }

    update();
}

#else

void P3DUserInterfaceCallback::initPCADeform() {
}

void P3DUserInterfaceCallback::loadPCAFile() {
}

void P3DUserInterfaceCallback::pcaDeform() {
}

#endif	/* PCA_VISUALIZER */

// ---------------------------  End of PCA Code  ---------------------------

void P3DUserInterfaceCallback::startRecording()
{
    P3DControl * control;

    if(userInterface == NULL)
        return;

    control = userInterface->control;
    if(control == NULL)
        return;

    control->recordTransformation();

    update();
}

void P3DUserInterfaceCallback::endRecording()
{
    P3DControl * control;

    if(userInterface == NULL)
        return;

    control = userInterface->control;
    if(control == NULL)
        return;

    control->stopRecordingTransformation();
}

void P3DUserInterfaceCallback::resetObject()
{
    P3DControl * control;

    if(userInterface == NULL)
        return;

    control = userInterface->control;
    if(control == NULL)
        return;

    control->resetTransformation();

    update();
}

void P3DUserInterfaceCallback::testTransformation()
{
    P3DControl * control;

    if(userInterface == NULL)
        return;

    control = userInterface->control;
    if(control == NULL)
        return;

    control->resetTransformation();
    control->applyRecordedTransformation();

    update();
}

// Called from discardMatchSurfaces()
void P3DUserInterfaceCallback::resetMatchSurfaceSerialNumber()
{
	userInterface->matchSurfaceSerialNumber->value(0);
    matchSurfaceSerialNumberChanged(0);
}

#ifndef PRODUCTION_VERSION

void P3DUserInterfaceCallback::initializeSubfigureTransformation()
{
    P3DControl * control;

    control = userInterface->control;
    if(control == NULL)
        return;

    control->initializeSubfigureTransformation();
}

void P3DUserInterfaceCallback::translateSubfigure()
{
    P3DControl * control;
    Fl_Counter * uCounter;
    Fl_Counter * vCounter;

    static double lastUVal = 0.0;
    static double lastVVal = 0.0;

    double uTrans,
           vTrans;

    control = userInterface->control;
    if(control == NULL)
        return;

    uCounter = userInterface->uTranslationCounter;
    vCounter = userInterface->vTranslationCounter;

    if(uCounter == NULL || vCounter == NULL)
        return;

    uTrans = uCounter->value() - lastUVal;
    vTrans = vCounter->value() - lastVVal;

    control->translateSubfigure(uTrans, vTrans);

    lastUVal = uCounter->value();
    lastVVal = vCounter->value();

    update();
}

void P3DUserInterfaceCallback::rotateSubfigure()
{
    P3DControl * control;
    Fl_Counter * rotationCounter;

    static double lastRotationVal = 0.0;

    double rotationAngle;

    control = userInterface->control;
    if(control == NULL)
        return;

    rotationCounter = userInterface->surfaceRotationCounter;

    if(rotationCounter == NULL)
        return;

    rotationAngle = (rotationCounter->value() - lastRotationVal) *
                    R_DEGREES_TO_RADIANS;

    control->rotateSubfigure(rotationAngle);

    lastRotationVal = rotationCounter->value();

    update();
}

void P3DUserInterfaceCallback::scaleSubfigure()
{
    P3DControl * control;
    Fl_Counter * scaleCounter;

    static double lastScaleVal = 0.0;

    double scale;

    control = userInterface->control;
    if(control == NULL)
        return;

    scaleCounter = userInterface->subfigureScaleCounter;

    if(scaleCounter == NULL)
        return;

    scale = scaleCounter->value() - lastScaleVal;
    if(scale >= 0)
        scale = scale + 1.0;
    else
        scale = 1.0 / (1.0 - scale);

    control->scaleSubfigure(scale);

    lastScaleVal = scaleCounter->value();

    update();
}

void P3DUserInterfaceCallback::elongateSubfigure()
{
    P3DControl * control;
    Fl_Counter * elongateCounter;
    double elongate;

    static double lastElongateVal = 0.0;

    if(userInterface == NULL)
        return;

    control = userInterface->control;
    if(control == NULL)
        return;

    elongateCounter = userInterface->subfigureElongationCounter;

    if(elongateCounter == NULL)
        return;

    elongate = elongateCounter->value() - lastElongateVal;
    if(elongate >= 0)
        elongate = elongate + 1.0;
    else
        elongate = 1.0 / (1.0 - elongate);

    control->elongateSubfigure(elongate);

    lastElongateVal = elongateCounter->value();

    update();
}

void P3DUserInterfaceCallback::hingeSubfigure()
{
    P3DControl * control;
    Fl_Counter * hingeCounter;
    double hingeAngle;

    if(userInterface == NULL)
        return;

    static double lastHingeVal = 0.0;

    control = userInterface->control;
    if(control == NULL)
        return;

    hingeCounter = userInterface->hingeRotationCounter;
    if(hingeCounter == NULL)
        return;

    hingeAngle = (hingeCounter->value() - lastHingeVal) *
		R_DEGREES_TO_RADIANS;
    control->hingeSubfigure(hingeAngle);

    lastHingeVal = hingeCounter->value();

    update();
}

void P3DUserInterfaceCallback::initializeTestConstraints()
{
    P3DControl * control;

    if(userInterface == NULL)
        return;

    control = userInterface->control;
    if(control == NULL)
        return;

    if(control->haveModel() == false) {
        cerr << "A model must first be loaded" << endl;
		fl_beep(FL_BEEP_ERROR);
        return;
    }

    if(control->haveImage() == false) {
        cerr << "An image must first be loaded" << endl;
		fl_beep(FL_BEEP_ERROR);
        return;
    }

	int surfaceLevel = userInterface->modelView->getSurfaceLevel();

    (void) control->initializeMatch(GAUSSIAN_DERIVATIVE_MATCH, surfaceLevel);
}

void P3DUserInterfaceCallback::testConstraintsPenalty()
{
#ifndef BINARY
    Fl_Value_Input * figureIdInput;
    Fl_Value_Output * constraintsOutput;
    P3DControl * control;
    double value;

    if(userInterface == NULL)
        return;

    control = userInterface->control;
    if(control == NULL)
        return;

    figureIdInput = userInterface->testConstraintsFigureNumberInput;
    if(figureIdInput == NULL)
        return;

    value = control->testConstraintsPenalty((int) figureIdInput->value());

    constraintsOutput = userInterface->testConstraintsOutput;
    if(constraintsOutput == NULL)
        return;

    constraintsOutput->value(value);
#endif
}

void P3DUserInterfaceCallback::testSeurat()
{
#ifndef BINARY
    Fl_Value_Input * figureIdInput;
    Fl_Value_Input * levelInput;
    P3DControl * control;

    if(userInterface == NULL)
        return;

    control = userInterface->control;
    if(control == NULL)
        return;

    figureIdInput = userInterface->testSeuratFigureNumberInput;
    if(figureIdInput == NULL)
        return;

    levelInput = userInterface->testSeuratLevelInput;
    if(levelInput == NULL)
        return;

    control->testSeurat((int) figureIdInput->value(), (int) levelInput->value());
#endif
}

void P3DUserInterfaceCallback::printModel(int worldCoords)
{
    P3DControl * control;
    M3DObject * object;
	Image3D * image;
	int figIndex;
	int primIndex;

    if(userInterface == NULL)
        return;

    control = userInterface->control;
    if(control == NULL)
        return;

    object = control->getObjectPtr();
	if (object== NULL) {
		cout << "An object must first be loaded" << endl;
		fl_beep(FL_BEEP_ERROR);
		return;
	}
	if (worldCoords) {
		image = control->getImagePtr();
		if (image == NULL) {
			cout << "An image must be loaded to compute world coordinates" << endl;
			fl_beep(FL_BEEP_ERROR);
			return;
		}
	}
	else
		image = NULL;
	if (image != NULL)
		M3DPrimitive::setWorld(image);
	else
		M3DPrimitive::setWorld(NULL);

	cout << "\n==========  Current Model ==========\n";
	control->getRelativeMarkedPrimitiveId(figIndex, primIndex);
	object->print(figIndex, primIndex);
}

void P3DUserInterfaceCallback::printImage()
{
    P3DControl * control;
    Image3D * image;

    if(userInterface == NULL)
        return;

    control = userInterface->control;
    if(control == NULL)
        return;

	image = userInterface->control->getImagePtr();
	if(image == NULL)
		return;

	image->print(displayImagePath);
}

void P3DUserInterfaceCallback::toggleBlendingRegionDisplay()
{
    Fl_Menu_Item * blendRegionItem;

    if(userInterface == NULL)
        return;

    blendRegionItem = userInterface->showBlendingRegion;
    if (blendRegionItem->value() == 0)
        userInterface->modelView->setBlendingRegionDisplay(false);
    else
        userInterface->modelView->setBlendingRegionDisplay(true);
}

// Called when the Display Match Surface button is pressed on the
// Surface Visualizer Control window.  Also called from
// discardMatchSurfaces().
void P3DUserInterfaceCallback::toggleMatchSurfaceDisplay()
{
	if (userInterface->showMatchSurfacesButton->value() == 0)
		userInterface->modelView->toggleSurfaceColoring(false);
	else
		userInterface->modelView->toggleSurfaceColoring(true);

	updateSurface();
}

// Called when either intensity slider on the Surface Visualizer
// Control window is used
void P3DUserInterfaceCallback::setMatchSurfaceColorRange()
{
	double min, max;

	if (0 != (int) userInterface->tieToImageSlidersButton->value()) {
		Image3D * image;
		Fl_Slider * minSlider;
		Fl_Slider * maxSlider;

		image = userInterface->control->getImagePtr();
		if(image == NULL)
			return;

		minSlider = userInterface->intensityMinSlider;
		maxSlider = userInterface->intensityMaxSlider;

		if (minSlider == NULL || maxSlider == NULL)
			return;

		min = image->mapActualToRelative((int) minSlider->value());
		max = image->mapActualToRelative((int) maxSlider->value());
	}
	else {
		Fl_Value_Slider * minSlider;
		Fl_Value_Slider * maxSlider;

		minSlider = userInterface->matchIntensityMinSlider;
		maxSlider = userInterface->matchIntensityMaxSlider;

		if (minSlider == NULL || maxSlider == NULL)
			return;

		min = minSlider->value();
		max = maxSlider->value();
	}

	userInterface->modelView->setSurfaceColorRange((float) min,
		(float) max);
}

// Called when the match surface (visualization) serial number is changed
// on the Surface Visualizer Control
void P3DUserInterfaceCallback::matchSurfaceSerialNumberChanged(int val)
{
#ifdef SURFACE_VISUALIZER
	if (objectSurfaceVisualizer == NULL) return;

	Fl_Counter * counter = userInterface->matchSurfaceSerialNumber;
	int newValue = objectSurfaceVisualizer->selectBySerialNumber(val);
	if (val != newValue)
		counter->value(newValue);

	Fl_Button * textOutput = userInterface->matchSurfaceSerialName;
	textOutput->label(objectSurfaceVisualizer->visualizationName());
	textOutput->redraw();

	updateSurface();
#endif
}

// Available for use with the Surface Visualizer Control window.  Not
// currently used.
void P3DUserInterfaceCallback::setMatchSurfaceSerialName(char * name)
{
/*    if (userInterface == NULL)
        return;

	Fl_Button * textOutput = userInterface->matchSurfaceSerialName;

	if (name == NULL || name[0] == '\0')
		textOutput->label("");
	else
		textOutput->label(name);
	textOutput->redraw();*/
}

// Called when the Surface Visualizer Control is opened by selecting
// the entry on the Test/Debug pulldown 
void P3DUserInterfaceCallback::openMatchSurfacesWindow()
{
#ifdef SURFACE_VISUALIZER
	int i = (int) surfaceColorMap.currentMap();
	userInterface->matchSurfaceColormapChoice->value(i);
	userInterface->matchSurfacesDialog->show();
#endif
}

// Called when a selection is made from the matchSurfaceColormapChoice
// pulldown on the Surface Visualizer Control window
void P3DUserInterfaceCallback::chooseMatchSurfaceColormap()
{
#ifdef SURFACE_VISUALIZER
	int i = userInterface->matchSurfaceColormapChoice->value();
	surfaceColorMap.specifyMap((SurfaceColorMap::color_map_type) i);

	updateSurface();
#endif	/* SURFACE_VISUALIZER */
}

// Called when the tieToImageSlidersButton checkbox is toggled on the
// Surface Visualizer Control window
void P3DUserInterfaceCallback::tieSurfaceToImageSlider(int val)
{
	if (val == 0)
		userInterface->modelView->tieSurfaceToImageSlider(false);
	else
		userInterface->modelView->tieSurfaceToImageSlider(true);
	setMatchSurfaceColorRange();
}

void P3DUserInterfaceCallback::togglePlainMatchSurfaces(int yesNo)
{
	userInterface->modelView->togglePlainMatchSurfaces((bool) (yesNo != 0));
}

void P3DUserInterfaceCallback::toggleTileType()
{
    Fl_Menu_Item * flatTilesItem;

    if(userInterface == NULL)
        return;

    flatTilesItem = userInterface->useFlatTiles;
    if (flatTilesItem->value() == 0)
        userInterface->modelView->toggleTileType(false);
    else
        userInterface->modelView->toggleTileType(true);
}

#else	/* PRODUCTION_VERSION */

void P3DUserInterfaceCallback::initializeSubfigureTransformation() {
}

void P3DUserInterfaceCallback::translateSubfigure() {
}

void P3DUserInterfaceCallback::rotateSubfigure() {
}

void P3DUserInterfaceCallback::scaleSubfigure() {
}

void P3DUserInterfaceCallback::elongateSubfigure() {
}

void P3DUserInterfaceCallback::hingeSubfigure() {
}

void P3DUserInterfaceCallback::initializeTestConstraints() {
}

void P3DUserInterfaceCallback::testConstraintsPenalty() {
}

void P3DUserInterfaceCallback::testSeurat() {
}

void P3DUserInterfaceCallback::printModel(int worldCoords) {
}

void P3DUserInterfaceCallback::printImage() {
}

void P3DUserInterfaceCallback::toggleBlendingRegionDisplay() {
}

void P3DUserInterfaceCallback::toggleMatchSurfaceDisplay() {
}

void P3DUserInterfaceCallback::setMatchSurfaceColorRange() {
}

void P3DUserInterfaceCallback::matchSurfaceSerialNumberChanged(int val) {
}

void P3DUserInterfaceCallback::setMatchSurfaceSerialName(char * name) {
}

void P3DUserInterfaceCallback::openMatchSurfacesWindow() {
}

void P3DUserInterfaceCallback::chooseMatchSurfaceColormap() {
}

void P3DUserInterfaceCallback::tieSurfaceToImageSlider(int val) {
}

void P3DUserInterfaceCallback::togglePlainMatchSurfaces(int yesNo) {
}

void P3DUserInterfaceCallback::toggleTileType() {
}

#endif	/* PRODUCTION_VERSION */

// This function returns false if the action was canceled
bool P3DUserInterfaceCallback::askToSaveModel()
{
    P3DControl * control;
    char * filename;
    char * msg;

    int size;
    int ret;

    if(userInterface == NULL)
        return true;

    control = userInterface->control;
    if(control == NULL)
        return true;

    if (control->modelIsUnsaved())
    {
		// This test is needed because FLTK can be forced to forget
		// that a modal window is open by making a change of focus
		// using the mouse.
		if (usingModalWindow > 0) {
			fl_beep(FL_BEEP_ERROR);
			return false;
		}
        filename = control->getModelFileName();
        if (filename != NULL) {
            size = strlen(SAVE_ON_EXIT_MESSAGE) + strlen(filename) + 1;
            msg = new char[size];
            msg[0] = '\0';
            sprintf(msg, SAVE_ON_EXIT_MESSAGE, filename);
            delete [] filename;
        }
        else {
            size = strlen(SAVE_ON_EXIT_MESSAGE) + strlen(NONAME) + 1;
            msg = new char[size];
            msg[0] = '\0';
            sprintf(msg, SAVE_ON_EXIT_MESSAGE, NONAME);
        }

        ret = fl_choice(msg, "No", "Yes", "Cancel");
        delete [] msg;
        if (ret == 1)		// Yes
            saveModel();
		else if (ret == 2)	// Cancel
			return false;
    }
	return true;
}

// Show the long version history to stdout and the
// short version history in the "About Pablo" window
void P3DUserInterfaceCallback::initAboutPabloWindow() {
#ifdef BINARY
	userInterface->aboutBinaryPabloVersion->label(revision);
	userInterface->aboutBinaryPabloRelease->label(release);

	userInterface->aboutBinaryPabloSplashGrp->image(new Fl_Tiled_Image(
		new Fl_Pixmap((const char * const *)bp_splash)));

	userInterface->aboutBinaryPabloVersionLabel->clear_visible_focus();
	userInterface->aboutBinaryPabloReleaseLabel->clear_visible_focus();
#else
	userInterface->aboutPabloHeader->label(
		"                           Pablo\n  is a product of the Medical Image\n   Display &&\
 Analysis Group of the\n    University of North Carolina at\nChapel Hill.  \251 2002 - 2009 UNC-CH.");
	userInterface->aboutPabloAuthors->label(
		"Authors:  Eli Broadhurst, Ed Chaney, Brad Davis,\n\
Tom Fletcher, Mark Foskey, Graham Gash,\nQiong Han, Ja Yeon Jeong, \
Sarang Joshi, Josh\nLevy, Shawn Liu, Xiaoxiao Liu, Conglin Lu,\n\
Derek Merck, Stephen Pizer, Rohit Saboo,\nJoshua Stough, Martin Styner, Andrew Thall,\n\
Gregg Tracton, and Paul Yushkevich.");
	userInterface->aboutPabloRelease->label(release);
	userInterface->aboutPabloVersion->label(revision);
	userInterface->aboutPabloHeader->clear_visible_focus();
	userInterface->aboutPabloAuthors->clear_visible_focus();
	userInterface->aboutPabloVersion->clear_visible_focus();
#endif
}

void P3DUserInterfaceCallback::showAboutPabloWindow() {
    if (userInterface == NULL)
        return;
#ifdef BINARY
    userInterface->aboutBinaryPabloWindow->show();
#else
    userInterface->aboutPabloWindow->show();
#endif
}

// ----------------------------- Start of Landmarks Code --------------------------------

#ifndef LM_METHOD_OBJ

// Called when Edit -> Edit Landmarks is selected
void P3DUserInterfaceCallback::startEditLandmarks() {
    P3DControl * control;
    M3DObject * object;
    P3DView * modelView;
	int numLandmarks, numImageLandmarks;

    if (userInterface == NULL)
        return;

    modelView = userInterface->modelView;
    if (modelView == NULL)
		return;

	if (userInterface->editLandmarksWindow->shown())
		return;

    modelView->setLandmarkMode(true);

    control = userInterface->control;
    if (control == NULL)
        return;

    object = control->getObjectPtr();
    if (object == NULL)
		return;
	if (object->getFigureCount() <= 0)
		return;

	deselectAll();	// Prohibits removal of figures
	userInterface->modelView->updateObject();

	if (modelView->getSurfaceStyle() == M3D_NONE) {
		setSurfaceStyle(M3D_SOLID);
		userInterface->setSurfaceTypeBtns(M3D_SOLID);
	}

	// Turn off blending so landmarks can be displayed
	userInterface->modelView->toggleBlendSubfigures(false);
	userInterface->blendFiguresBtn->value(0);

	modelView->initializeEditLandmarks();

	embedLandmarkNameEditor(object);

	// Preset the editing buttons
	numLandmarks = object->getLandmarkCount();
	numImageLandmarks = control->getImageLandmarksCount();

	if (numLandmarks > 0 && numLandmarks != numImageLandmarks) {
		// Need to add image landmarks
		setLandmarksEditMode(1);
		if (editingLandmark < 0)	// New model was loaded
			editingLandmark = 0;
	}
	else						// Need to add surface landmarks
		setLandmarksEditMode(0);

	userInterface->editLandmarks_Caret0->hide();
	userInterface->editLandmarks_Caret1->hide();
	userInterface->editLandmarks_Caret2->hide();
	userInterface->editLandmarks_Caret3->hide();
	userInterface->editLandmarks_Caret4->hide();
	userInterface->editLandmarks_Caret5->hide();
	landmarkCaret = -1;
	setLandmarkCaret(editingLandmark, false);	// Assumes editingLandmark is correct

	setLandmarksSurfaceNumber(numLandmarks);
	setLandmarksImageNumber(numImageLandmarks);

	int x = userInterface->editLandmarkNames->x();
	int y = userInterface->editLandmarkNames->y();
	int h = userInterface->editLandmarkNames->h();
	if (numLandmarks <= 6) {
		userInterface->editLandmarkNames->resize(x, y, 300, h);
		userInterface->editLandmarkNames->label(editLandmarkNamesLabel_0);
	}
	else {
		userInterface->editLandmarkNames->resize(x, y, 317, h);
		userInterface->editLandmarkNames->label(editLandmarkNamesLabel_1);
	}
	priorLandmarkVisibility = (userInterface->landmarksVisible->value() ? true : false);
	userInterface->modelView->viewLandmarks(true);
	userInterface->editLandmarksWindow->show();
}

// Called when the landmarks surface and image buttons are pressed
void P3DUserInterfaceCallback::setLandmarksEditMode(int mode)
{
#ifdef DEBUG
	cout << "P3DUserInterfaceCallback::setLandmarksEditMode(" << mode << ')' << endl;
#endif

/*    P3DControl * control;
    M3DObject * object;
	M3DFigure * figure;

	if (landmarksEditMode == mode) {
		// Turn landmark editing off
		userInterface->landmarksEditMode_0Btn->value(0);
		userInterface->landmarksEditMode_1Btn->value(0);
		landmarksEditMode = -1;

		// Unmark the marked landmark
		control = userInterface->control;
		object = control->getObjectPtr();

		int numFigures = object->getFigureCount();
		for (int figureId = 0; figureId < numFigures; figureId++) {
			figure = object->getFigurePtr(figureId);
			if (figure->getMarkedLandmark() >= 0) {
				figure->unmarkLandmark();
				break;
			}
		}
		updateSurface();
		return;
	}
*/
	landmarksEditMode = mode;
	if (mode) {
		userInterface->landmarksEditMode_0Btn->value(0);
		userInterface->landmarksEditMode_1Btn->value(1);
		initializeImageLandmarkEditing();
	}
	else {
		userInterface->landmarksEditMode_0Btn->value(1);
		userInterface->landmarksEditMode_1Btn->value(0);
		initializeSurfaceLandmarkEditing();
	}
}

void P3DUserInterfaceCallback::initializeSurfaceLandmarkEditing()
{
    P3DControl * control;
    M3DObject * object;
	M3DFigure * figure;
	int numFigures, figureId;
	int n;

#ifdef DEBUG
	cout << "P3DUserInterfaceCallback::initializeSurfaceLandmarkEditing()" << endl;
#endif
    control = userInterface->control;
    object = control->getObjectPtr();

	n = object->getMarkedLandmark(figureId);
//	cout << "Landmark " << l << " on figure " << figureId << " was marked" << endl;
	if (figureId >= 0) {
		figure = object->getFigurePtr(figureId);
		figure->unmarkLandmark();
	}

	// Determine which landmark to designate, by cycling through the
	// list of surface landmarks.
	editingLandmark = -1;
	numFigures = object->getFigureCount();
	for (figureId = 0; figureId < numFigures; figureId++) {
        figure = object->getFigurePtr(figureId);
		n = figure->getLandmarkCount();
		if (n > 0) {
			figure->markLandmark(0);	// Mark the first one found
			editingLandmark = 0;
			break;
		}
	}
	setLandmarkCaret(editingLandmark, true);

//	n = object->getMarkedLandmark(figureId);
//	cout << "Landmark " << n << " on figure " << figureId << " is marked" << endl;
	if (editingLandmark == -1)
		return;		// There are no landmarks, since there is no first landmark

    updateSurface();
}

void P3DUserInterfaceCallback::initializeImageLandmarkEditing()
{
    P3DControl * control;
    P3DView * modelView;
    M3DObject * object;
	M3DFigure * figure;
	int numFigures, figureId;
	int n, l;

#ifdef DEBUG
	cout << "P3DUserInterfaceCallback::initializeImageLandmarkEditing()" << endl;
#endif
    control = userInterface->control;
    object = control->getObjectPtr();
    modelView = userInterface->modelView;

	l = object->getMarkedLandmark(figureId);
//	cout << "Landmark " << l << " on figure " << figureId << " was marked" << endl;
	if (figureId >= 0) {
		figure = object->getFigurePtr(figureId);
		figure->unmarkLandmark();
	}

	// Determine which landmark to designate, by cycling through the surface
	// landmarks to find the first one without a corresponding image landmark.
	editingLandmark = -1;
	bool found = false;
	numFigures = object->getFigureCount();
	for (figureId = 0; figureId < numFigures; figureId++) {
        figure = object->getFigurePtr(figureId);
		n = figure->getLandmarkCount();
		for (l = 0; l < n; l++) {
			editingLandmark++;
			if (! control->landmarkIsPaired(figureId, l)) {
				found = true;
//				cout << "Landmark " << l << " on figure " << figureId << " is not paired" << endl;
				figure->markLandmark(l);
				break;
			}
		}
		if (found)
			break;
	}
	if (! found) {	// There are no unpaired landmarks
		if (editingLandmark >= 0) {
			// There are landmarks; find and mark the first one
			editingLandmark = 0;
			for (figureId = 0; figureId < numFigures; figureId++) {
				figure = object->getFigurePtr(figureId);
				n = figure->getLandmarkCount();
				if (n > 0) {
					figure->markLandmark(0);
					break;
				}
			}
		}
	}

	setLandmarkCaret(editingLandmark, true);

//	n = object->getMarkedLandmark(figureId);
//	cout << "Landmark " << n << " on figure " << figureId << " is marked" << endl;
	if (editingLandmark == -1)
		return;		// There are no landmarks, since there is no first landmark

    updateSurface();
}

// Called from P3DView after an image landmark is added to move the caret
// to the next unpaired landmark
void P3DUserInterfaceCallback::nextImageLandmark()
{
    P3DControl * control;
    P3DView * modelView;
    M3DObject * object;
	M3DFigure * figure;
	int figureId, numLandmarks;
	int n, l, oldIndex;

#ifdef DEBUG
	cout << "P3DUserInterfaceCallback::nextImageLandmark()" << endl;
#endif
	if (editingLandmark < 0) {
		initializeImageLandmarkEditing();
		return;
	}

    control = userInterface->control;
    object = control->getObjectPtr();

	numLandmarks = object->getLandmarkCount();
	if (numLandmarks == 1)
		return;

	// Unmark the current landmark
    modelView = userInterface->modelView;
	figureId = control->figureOfLandmark(editingLandmark, l);
	figure = object->getFigurePtr(figureId);
	figure->unmarkLandmark();

	oldIndex = editingLandmark;	// Remember the previous serial number

	// Cycle through the surface landmarks to find the next one that
	// is not paired to an image landmark.
	editingLandmark++;
	if (editingLandmark >= numLandmarks)
		editingLandmark = 0;
	bool found = false;
	while (editingLandmark != oldIndex) {
		figureId = control->figureOfLandmark(editingLandmark, l);
		if (figureId < 0) {
			// This landmark is unpaired
			found = true;
			int count = 0;
			// Find its corresponding surface landmark
			for (figureId = 0; figureId < object->getFigureCount(); figureId++) {
				figure = object->getFigurePtr(figureId);
				n = figure->getLandmarkCount();
				if (editingLandmark < count + n) {
					l = editingLandmark - count;
//					cout << "Marking landmark number " << l << " on figure " << figureId << endl;
					figure->markLandmark(l);
					break;
				}
				count += n;
			}
		}
		if (found)
			break;

		editingLandmark++;
		if (editingLandmark >= numLandmarks)
			editingLandmark = 0;
	};

	if (! found) {	// There are no unpaired landmarks
		editingLandmark++;
		if (editingLandmark >= numLandmarks)
			editingLandmark = 0;
		figureId = control->figureOfLandmark(editingLandmark, l);
		figure = object->getFigurePtr(figureId);
		figure->markLandmark(l);
//		cout << "Marking landmark number " << l << " on figure " << figureId << endl;
	}

	setLandmarkCaret(editingLandmark, true);

//	n = object->getMarkedLandmark(figureId);
//	cout << "Landmark " << n << " on figure " << figureId << " is marked" << endl;
    updateSurface();
}

// This function forces the landmarkNameEditor scrollbar to adjust so the caret line is visible
void P3DUserInterfaceCallback::setLandmarkSlider()
{
	int top;
	Fl_Scrollbar * s;
	Fl_Valuator * v;
	Flv_List * l;

	top = landmarkNameEditor->scrollbar.value();
	s = &(landmarkNameEditor->scrollbar);
	v = (Fl_Valuator *) s;
	l = (Flv_List *) landmarkNameEditor;
	if (editingLandmark < top)
		v->value(editingLandmark);
	else	// This assumes editingLandmark > bottom and an landmarkNameEditor length of 6 lines
		v->value(editingLandmark - 5);
	l->top_row(s->value());
	l->damage(FL_DAMAGE_CHILD);
}

void P3DUserInterfaceCallback::setLandmarkCaret(int row, bool adjustSlider)
{
	int top;

//	cout << "setLandmarkCaret(" << row << ", " << adjustSlider << ')' << endl;
	switch (landmarkCaret) {
	  case 0:	userInterface->editLandmarks_Caret0->hide();
				break;
	  case 1:	userInterface->editLandmarks_Caret1->hide();
				break;
	  case 2:	userInterface->editLandmarks_Caret2->hide();
				break;
	  case 3:	userInterface->editLandmarks_Caret3->hide();
				break;
	  case 4:	userInterface->editLandmarks_Caret4->hide();
				break;
	  case 5:	userInterface->editLandmarks_Caret5->hide();
				break;
	  default:	break;
	}

	landmarkCaret = -1;

	if (row < 0)
		return;

	if (! landmarkNameEditor)
		return;

	top = landmarkNameEditor->scrollbar.value();
//	cout << "top = " << top << endl;
	row -= top;
	switch (row) {
	  case 0:	userInterface->editLandmarks_Caret0->show();
				break;
	  case 1:	userInterface->editLandmarks_Caret1->show();
				break;
	  case 2:	userInterface->editLandmarks_Caret2->show();
				break;
	  case 3:	userInterface->editLandmarks_Caret3->show();
				break;
	  case 4:	userInterface->editLandmarks_Caret4->show();
				break;
	  case 5:	userInterface->editLandmarks_Caret5->show();
				break;
	  default:	// Slider has been adjusted so the current landmark name is not showing
				row = -1;
				break;
	}
	landmarkCaret = row;

	if (adjustSlider) {
		// If the slider is positioned so the entry for editingLandmark is
		// not visible, the previous code will set landmarkCaret to -1.
		if (landmarkCaret < 0) {
			// Adjust the landmarks slider and set the caret
			setLandmarkSlider();
			setLandmarkCaret(editingLandmark, false);
		}
	}
}

// Called when the Next button on the landmarks window is pressed
void P3DUserInterfaceCallback::editLandmarksNext()
{
    P3DControl * control;
    M3DObject * object;
	M3DFigure * figure;
	int index, n;
	int numFigures, figureId;
	int numLandmarks;
	bool found;

#ifdef DEBUG
	cout << "P3DUserInterfaceCallback::editLandmarksNext()" << endl;
#endif
	if (editingLandmark < 0)
		return;

//	if (landmarksEditMode < 0)	// Both edit-type buttons off
//		return;

    control = userInterface->control;
    object = control->getObjectPtr();

	editingLandmark++;
	if (editingLandmark >= object->getLandmarkCount())
		editingLandmark = 0;

	// Cycle through the surface landmarks, clearing the old mark and setting
	// the new one to correspond with the new value of editingLandmark.  Also
	// update the caret to match.
	numFigures = object->getFigureCount();
	numLandmarks = object->getLandmarkCount();
	index = 0;
	found = false;
	for (figureId = 0; figureId < numFigures; figureId++) {
        figure = object->getFigurePtr(figureId);
		n = figure->getLandmarkCount();
		if (n > 0)
			figure->unmarkLandmark();
		if (! found && editingLandmark < index + n) {
			figure->markLandmark(editingLandmark - index);	// Mark the current landmark
			setLandmarkCaret(editingLandmark, true);
			found = true;
		}
		index += n;
	}

//	n = object->getMarkedLandmark(figureId);
//	cout << "Landmark " << n << " on figure " << figureId << " is marked" << endl;
	updateSurface();
}

// Called when the Previous button on the landmarks window is pressed
void P3DUserInterfaceCallback::editLandmarksPrevious()
{
    P3DControl * control;
    M3DObject * object;
	M3DFigure * figure;
	int index, n;
	int numLandmarks;
	int numFigures, figureId;
	bool found;

#ifdef DEBUG
	cout << "P3DUserInterfaceCallback::editLandmarksPrevious()" << endl;
#endif
	if (editingLandmark < 0)
		return;

//	if (landmarksEditMode < 0)	// Both edit-type buttons off
//		return;

    control = userInterface->control;
    object = control->getObjectPtr();

	editingLandmark--;
	if (editingLandmark < 0)
		editingLandmark = object->getLandmarkCount() - 1;

	// Cycle through the surface landmarks, clearing the old mark and setting
	// the new one to correspond with the new value of editingLandmark.  Also
	// update the caret to match.
	numFigures = object->getFigureCount();
	numLandmarks = object->getLandmarkCount();
	index = 0;
	found = false;
	for (figureId = 0; figureId < numFigures; figureId++) {
        figure = object->getFigurePtr(figureId);
		n = figure->getLandmarkCount();
		if (n > 0)
			figure->unmarkLandmark();
		if (! found && editingLandmark < index + n) {
			figure->markLandmark(editingLandmark - index);	// Mark the current landmark
			setLandmarkCaret(editingLandmark, true);
			found = true;
		}
		index += n;
	}

//	n = object->getMarkedLandmark(figureId);
//	cout << "Landmark " << n << " on figure " << figureId << " is marked" << endl;
	updateSurface();
}

// Called when the Delete Landmark button is pressed
void P3DUserInterfaceCallback::editLandmarkDelete()
{
    P3DControl * control;
    M3DObject * object;
	M3DFigure * figure;
    P3DView * modelView;
	int figureId;

#ifdef DEBUG
	cout << "P3DUserInterfaceCallback::editLandmarkDelete()" << endl;
#endif
	if (editingLandmark < 0)
		return;

//	if (landmarksEditMode < 0)	// Both edit-type buttons off
//		return;

    control = userInterface->control;
    object = control->getObjectPtr();
    modelView = userInterface->modelView;

#ifdef DEBUG
	cout << "Deleting landmark " << editingLandmark << endl;
#endif

	// Delete the marked landmark
	int n = object->getMarkedLandmark(figureId);
	if (n >= 0) {
	    figure = object->getFigurePtr(figureId);
		if (landmarksEditMode == 0) {	// Surface Mode
			// Delete surface landmark and any corresponding image landmark
			figure->clearLandmark(n);	// Unmark and delete the landmark
			control->modelNeedsSaving();
			control->deleteImageLandmark(editingLandmark);
			setLandmarksImageNumber(control->getImageLandmarksCount());	// May have changed
			updateLandmarkEditorWindow();

			// Determine which landmark should next be designated for editing
			editingLandmark--;
			if (editingLandmark < 0)
				editingLandmark = object->getLandmarkCount() - 1;
			editLandmarksNext();		// Calls updateSurface(), unless
			if (editingLandmark < 0)	//   there are no more landmarks.
				updateSurface();
		}
		else {	// Image Mode
			markImageLandmark(editingLandmark, false);	// Clear the landmark's pairing display
			// Delete image landmark
			if (control->clearImageLandmark(editingLandmark))
				setLandmarksImageNumber(control->getImageLandmarksCount());
			// The mark should not move, since only the image landmark was deleted
			redraw();
		}
	}
//	n = object->getMarkedLandmark(figureId);
//	cout << "Landmark " << n << " on figure " << figureId << " is marked" << endl;
}

void P3DUserInterfaceCallback::editLandmarkGoto() {
    P3DControl * control;
    M3DObject * object;
	M3DFigure * figure;
    P3DView * modelView;
	int figureId;

	if (editingLandmark < 0)
		return;

    control = userInterface->control;
    object = control->getObjectPtr();
    modelView = userInterface->modelView;

	int n = object->getMarkedLandmark(figureId);
	if (n >= 0) {
	    figure = object->getFigurePtr(figureId);
		Vector3D v = figure->getLandmark(n);
		modelView->translateToCenter(v);
	}
}

void * P3DUserInterfaceCallback::generateLandmarkEditorCells(M3DObject * object, int numLandmarks)
{
    P3DControl * control;
	char str[3];
	int i, index;
	M3DFigure * figure;
	int figureId;

#ifdef DEBUG
	cout << "generateLandmarkEditorCells() called for " << numLandmarks << " landmarks" << endl;
#endif
    control = userInterface->control;

	// Assumes numLandmarks > 0
	userInterface->landmarkCells = new Flvt_Edit_Cell[4*numLandmarks];
	Flvt_Edit_Cell * cells = (Flvt_Edit_Cell *) userInterface->landmarkCells;

	// The column width setting effectively limits the string length.
	cells[numLandmarks].setCommonWidth(MAX_FLVW_STRING_LENGTH);
	index = 0;
	int numFigures = object->getFigureCount();
	for (figureId = 0; figureId < numFigures; figureId++) {
		figure = object->getFigurePtr(figureId);
		int lcount = figure->getLandmarkCount();
		for (i = 0; i < lcount; i++) {
			sprintf(str, "%d", figureId);			// First column is figureID
			cells[index].set(str);
			const char * name = figure->getLandmarkName(i);	// Second column is landmark name
			if (name == NULL)
				cells[index + numLandmarks].set(" ");
			else
				cells[index + numLandmarks].set(name);
			if (control->landmarkIsPaired(figureId, i))
				str[0] = '+';	// Paired with image landmark
			else
				str[0] = ' ';	// No image landmark
			str[1] = '\0';
			cells[index + numLandmarks + numLandmarks].set(str);    // Third column indicates pairing
			sprintf(str, "%d.", i);		// Fourth column is landmark number on figure and not visible
			cells[index + 3*numLandmarks].set(str);
			index++;
		}
	}
	return (void *) cells;	// Avoids having to define Flvt_Edit_Cell in the .h file
}

// This function is called when the landmarkNameEditor vertical scrollbar is changed
void P3DUserInterfaceCallback::landmarkEditScrollbarCallback(Fl_Widget * widget, void * user_data)
{
	P3DUserInterfaceCallback * self;

	self = (P3DUserInterfaceCallback *) user_data;
	Fl_Scrollbar * s;
	s = (Fl_Scrollbar *) widget;
	self->vscrollbar_cb(widget, user_data);
//	cout << "Scrollbar = " << s->value() << endl;

	// The whole point of this function is to update the caret
	self->setLandmarkCaret(self->editingLandmark, false);
}

// This function fills in the 3rd column entry in the landmark name editor table.
// A '+' will be displayed if onOff is true, otherwise a blank.
void P3DUserInterfaceCallback::markImageLandmark(int landmark, bool onOff)
{
	char str[2];
    P3DControl * control;
    M3DObject * object;

#ifdef DEBUG
	cout << "markImageLandmark(" << landmark << ", " << onOff << ')' << endl;
#endif
	str[1] = '\0';
	Flvt_Edit_Cell * cells = (Flvt_Edit_Cell *) userInterface->landmarkCells;
    control = userInterface->control;
    object = control->getObjectPtr();
	int numLandmarks = object->getLandmarkCount();

	if (onOff)
		str[0] = '+';
	else
		str[0] = ' ';
	cells[landmark + numLandmarks + numLandmarks].set(str);    // Third column indicates pairing
	landmarkNameEditor->redraw();
}

// Called from P3DView and when the number of surface landmarks changes
void P3DUserInterfaceCallback::updateLandmarkEditorWindow()
{
    P3DControl * control;
    M3DObject * object;
	Flvt_Edit_Cell * cells;

#ifdef DEBUG
	cout << "updateLandmarkEditorWindow()" << endl;
#endif
    control = userInterface->control;
    object = control->getObjectPtr();

	int numLandmarks = object->getLandmarkCount();
	setLandmarksSurfaceNumber(numLandmarks);

	delete [] (Flvt_Edit_Cell *) userInterface->landmarkCells;
	if (numLandmarks == 0) {
		userInterface->landmarkCells = NULL;
		if (landmarkNameEditor != NULL) {
			setLandmarkCaret(-1, false);
			userInterface->editLandmarkNames->remove(landmarkNameEditor);
			delete landmarkNameEditor;
			landmarkNameEditor = NULL;
			userInterface->editLandmarkNames->redraw();
		}
		return;
	}

	int x = userInterface->editLandmarkNames->x();
	int y = userInterface->editLandmarkNames->y();
//	int w = userInterface->editLandmarkNames->w();
	int h = userInterface->editLandmarkNames->h();
	if (numLandmarks <= 6) {
		userInterface->editLandmarkNames->resize(x, y, 300, h);
		userInterface->editLandmarkNames->label(editLandmarkNamesLabel_0);
	}
	else {
		userInterface->editLandmarkNames->resize(x, y, 317, h);
		userInterface->editLandmarkNames->label(editLandmarkNamesLabel_1);
	}

	cells = (Flvt_Edit_Cell * ) generateLandmarkEditorCells(object, numLandmarks);
	if (landmarkNameEditor == NULL) {
		int cols = 4;
		// This function call adds the Flvt_Edit widget to the editLandmarkNames group.
		// The column width below is the maximum that fits the group.
		landmarkNameEditor = flvt_cell_editor(userInterface->editLandmarkNames,
			cells, numLandmarks, cols, 248, 24, false, NULL, &editLandmarkName,
			(void *) this);
		landmarkNameEditor->col_style[3].width(0);
		vscrollbar_cb = landmarkNameEditor->scrollbar.callback();
		landmarkNameEditor->scrollbar.when(FL_WHEN_CHANGED);	// Probably not needed
		landmarkNameEditor->scrollbar.callback(&landmarkEditScrollbarCallback, this);
	}
	else {
		landmarkNameEditor->rows(numLandmarks);
		landmarkNameEditor->initialize(cells, numLandmarks);
		for (int i = numLandmarks; i < numLandmarks + numLandmarks; i++)	// Name column only
			cells[i].setCellChangeCallback(editLandmarkName, i);
		landmarkNameEditor->redraw();
		setLandmarkCaret(0, false);
	}
}

// This function is called when the user changes any landmark name.
// Tag is the index into the 2nd column (landmark name), treating the
// column-major table as a vector.
void P3DUserInterfaceCallback::editLandmarkName(void * classObject, int tag)
{
    P3DControl * control;
    M3DObject * object;
	P3DUserInterfaceCallback * self;
	P3DUserInterface * ui;
	Flvt_Edit_Cell * cells;
	const char * str;

//	cout << "Cell " << tag << " changed" << endl;
	self = (P3DUserInterfaceCallback *) classObject;
	ui = self->userInterface;
	control = ui->control;
    object = control->getObjectPtr();
	int numLandmarks = object->getLandmarkCount();
	cells = (Flvt_Edit_Cell *) ui->landmarkCells;

	// Record changes in the model
	str = cells[tag + numLandmarks + numLandmarks].textPtr();	// 4th column entry
	int n = atoi(str);	// Number of landmark within figure
//	cout << "n = " << n << endl;
	str = cells[tag - numLandmarks].textPtr();	// 1st column entry
	int figureId = atoi(str);	// Figure number
//	cout << "figureId = " << figureId << endl;
	M3DFigure * figure = object->getFigurePtr(figureId);
	figure->setLandmarkName(n, cells[tag].textPtr());
	control->modelNeedsSaving();
}

void P3DUserInterfaceCallback::embedLandmarkNameEditor(M3DObject * object)
{
	Flvt_Edit_Cell * cells;

#ifdef DEBUG
	cout << "embedLandmarkNameEditor()" << endl;
#endif
	int numLandmarks = object->getLandmarkCount();
	if (numLandmarks == 0) {
		userInterface->landmarkCells = NULL;
		landmarkNameEditor = NULL;
	}
	else {
		int cols = 4;
		cells = (Flvt_Edit_Cell * ) generateLandmarkEditorCells(object, numLandmarks);
		// This function call adds the Flvt_Edit widget to the editModelFigureNames group.
		// The column width below is the maximum that fits the group.
		landmarkNameEditor = flvt_cell_editor(userInterface->editLandmarkNames,
			cells, numLandmarks, cols, 248, 24, false, NULL, &editLandmarkName,
			(void *) this);
		landmarkNameEditor->col_style[3].width(0);
		vscrollbar_cb = landmarkNameEditor->scrollbar.callback();
		landmarkNameEditor->scrollbar.when(FL_WHEN_CHANGED);	// Not actually needed
		landmarkNameEditor->scrollbar.callback(&landmarkEditScrollbarCallback, this);
	}
}

// Called when the Clear Surface Landmarks button is pressed
void P3DUserInterfaceCallback::clearSurfaceLandmarks() {
    P3DControl * control;
    P3DView * modelView;
    M3DObject * object;
    M3DFigure * figure;
	int numFigures;

    control = userInterface->control;
    object = control->getObjectPtr();

	int numLandmarks = object->getLandmarkCount();
	if (numLandmarks == 0)
		return;

	if (landmarksEditMode != 0) {
		cout << "You must be in Surface Mode to delete all landmarks" << endl;
		fl_beep(FL_BEEP_ERROR);
		return;
	}

    modelView = userInterface->modelView;
	control->deleteImageLandmarks();
	setLandmarksImageNumber(0);

	numFigures = object->getFigureCount();
	for (int figureId = 0; figureId < numFigures; figureId++) {
        figure = object->getFigurePtr(figureId);
		figure->clearLandmarks();
	}
	control->modelNeedsSaving();

	updateLandmarkEditorWindow();
    updateSurface();
}

// Called when the Clear Image Landmarks button is pressed
void P3DUserInterfaceCallback::clearImageLandmarks() {
    P3DControl * control;
    P3DView * modelView;
    M3DObject * object;

    control = userInterface->control;
    object = control->getObjectPtr();

	int numLandmarks = object->getLandmarkCount();
	if (numLandmarks == 0)
		return;

    modelView = userInterface->modelView;
	if (control->getImageLandmarksCount() == 0)
		return;

	control->clearImageLandmarks();
	setLandmarksImageNumber(0);
	for (int i = 0; i < numLandmarks; i++)
		markImageLandmark(i, false);

	redraw();
}

void P3DUserInterfaceCallback::setLandmarksSurfaceNumber(int numLandmarks) {
		userInterface->landmarksSurfaceNumber->value(numLandmarks);
}

void P3DUserInterfaceCallback::setLandmarksImageNumber(int numImageLandmarks) {
		userInterface->landmarksImageNumber->value(numImageLandmarks);
}

void P3DUserInterfaceCallback::finishEditLandmarks() {
    P3DView * modelView;
    P3DControl * control;
    M3DObject * object;

	if (! userInterface->editLandmarksWindow->shown())
		return;

    modelView = userInterface->modelView;
    modelView->setLandmarkMode(false);

	landmarksEditMode = -1;		// Disables mouse selection of landmarks in P3DView.cpp

    control = userInterface->control;
    object = control->getObjectPtr();
	if (object == NULL)
		return;

	int numLandmarks = object->getLandmarkCount();
	Flvt_Edit_Cell * cells = (Flvt_Edit_Cell *) userInterface->landmarkCells;

	if (cells != NULL) {
		// Scan 2nd column entries, which contain landmark names 
		for (int i = numLandmarks; i < numLandmarks + numLandmarks; i++) {
			if (cells[i].changed())
				editLandmarkName((void *) this, i);
		}
	}

	userInterface->editLandmarksWindow->hide();
	userInterface->modelView->viewLandmarks(priorLandmarkVisibility);
	updateSurface();

	// Must remove the widget from the group before deleting it
	if (landmarkNameEditor != NULL) {
		userInterface->editLandmarkNames->remove(landmarkNameEditor);
		delete landmarkNameEditor;
		landmarkNameEditor = NULL;
	}
	if (cells != NULL) {
		delete [] (Flvt_Edit_Cell *) userInterface->landmarkCells;
		userInterface->landmarkCells = NULL;
	}
}

void P3DUserInterfaceCallback::applyLandmarks()
{
    P3DControl * control;

    control = userInterface->control;
	control->applyLandmarkTransformation();
	updateObject();
}

#else	/* BINARY */

void P3DUserInterfaceCallback::startEditLandmarks() {
}

void P3DUserInterfaceCallback::setLandmarksEditMode(int) {
}

void P3DUserInterfaceCallback::nextImageLandmark() {
}

void P3DUserInterfaceCallback::editLandmarksNext() {
}

void P3DUserInterfaceCallback::editLandmarksPrevious() {
}

void P3DUserInterfaceCallback::editLandmarkDelete() {
}

void P3DUserInterfaceCallback::editLandmarkGoto() {
}

//void * P3DUserInterfaceCallback::generateLandmarkEditorCells(M3DObject * object, int numLandmarks) {
//	return NULL;
//}

//void P3DUserInterfaceCallback::landmarkEditScrollbarCallback(Fl_Widget * widget, void * user_data) {
//}

void P3DUserInterfaceCallback::markImageLandmark(int,bool) {
}

void P3DUserInterfaceCallback::updateLandmarkEditorWindow() {
}

//void P3DUserInterfaceCallback::editLandmarkName(void * classObject, int tag) {
//}

//void P3DUserInterfaceCallback::embedLandmarkNameEditor(M3DObject * object) {
//}

void P3DUserInterfaceCallback::clearSurfaceLandmarks() {
}

void P3DUserInterfaceCallback::clearImageLandmarks() {
}

void P3DUserInterfaceCallback::setLandmarksSurfaceNumber(int numLandmarks) {
}

void P3DUserInterfaceCallback::setLandmarksImageNumber(int numImageLandmarks) {
}

void P3DUserInterfaceCallback::finishEditLandmarks() {
}

void P3DUserInterfaceCallback::applyLandmarks() {
}

#endif	/* BINARY */

// ----------------------------- End of Landmarks Code --------------------------------

// Special function for testing
// This function is called by typing CTRL-SHIFT-T in any Pablo window
void P3DUserInterfaceCallback::secretAction()
{
	// Put your test code here
//	cout << "Test\n";
}

void P3DUserInterfaceCallback::openOptimizationVisualizer()
{
#ifdef OPTIMIZATION_VISUALIZER
	if (ovUI) {
		((OptVisualizerUI *) ovUI)->callback->updateLabels();
		((OptVisualizerUI *) ovUI)->show();
		}
#endif
}

void P3DUserInterfaceCallback::loadBinaryPabloScript()
{
	cout << "loadBinaryPabloScript() is not implemented\n";
}

// Called to import surface color indexes for display using the
// Surface Visualizer Control window
void P3DUserInterfaceCallback::importSurfaceValues()
{
#ifdef SURFACE_VISUALIZER
	if (objectSurfaceVisualizer == NULL) return;

	int serNum = objectSurfaceVisualizer->getFreeSerialNumber();

	if (serNum < 0) {
		cerr << "Warning: Could not load a match surface; reached maximum\n"
			<< "   number of visualizations allowed" << endl;
		return;
	}

	Fl_File_Chooser fc(globalControl->readString(ModelDirectory), "*.txt",
	Fl_File_Chooser::MULTI, "Load Surface Values File");
	fc.callback(NULL);

	fc.show();
	while (fc.shown())
		Fl::wait();

	if (fc.count() == 0)
		return;		// User canceled request

	P3DControl * control = getControl();
	if (! control) return;

	// AGG: Bug: this is limited to a single figure
	if (control->getObjectPtr() && control->getObjectPtr()->getFigurePtr(0))
	{
		M3DFigure * fig = control->getObjectPtr()->getFigurePtr(0);
		Xferlist * xfer = convertM3DtoXfer(fig);

		int subdivLevel = userInterface->modelView->getSurfaceLevel();

		if (xfer) {							
			// AGG: Bug: an array is not needed
			ThallCode::Pointlist_server2** pList = new ThallCode::Pointlist_server2*[2];

			pList[0] = new ThallCode::Pointlist_server2;
			pList[0]->init(xfer);
			pList[0]->ComputeSubdivPointCloud(subdivLevel);

//			pList[1] = new ThallCode::Pointlist_server2;
//			pList[1]->init(xfer);
//			pList[1]->ComputeSubdivPointCloud(subdivLevel);

			delete [] (xfer->atomlist);
			delete xfer;

			Bpoint* bpList = 0;
			int numPoints;
			pList[0]->subdivboundaryinfo(&numPoints, &bpList);

			for (int idx = 1; idx <= fc.count(); idx++) {
				const char* filename = fc.value(idx);
				serNum = objectSurfaceVisualizer->getFreeSerialNumber();
				if (serNum < 0) {
					// we don't have anywhere to load this file into
					break;
				} 
				objectSurfaceVisualizer->initializeSurfaceMap(serNum, (char *) filename,
					userInterface->control->getObjectPtr(), 0 /* figure id */, 
					subdivLevel, numPoints, pList, false/*auto intensity window */);

				ifstream in(filename);
				double val;
				for (int countRead = 0; countRead < numPoints; countRead++) {
					in >> val;
					if (in) {
						objectSurfaceVisualizer->setSurfaceColor(0 /* figure */, countRead, val, serNum);
					}
					else {	// AGG: Why not just issue an error report here?
						objectSurfaceVisualizer->setSurfaceColor(0 /* figure */, countRead,
							(countRead&1) /* checkerboard indicates not enough data */, serNum);
					}
				}
				cout << "Loaded " << countRead << " of " << numPoints
					<< " values at surface level " << subdivLevel << ",\n   from file "
					<< filename << endl;
				userInterface->matchSurfaceSerialNumber->maximum(serNum);
			}
		}
	}
#endif	/* SURFACE_VISUALIZER */
}

// Called when the Discard button is pressed on the
// Surface Visualizer Control window
void P3DUserInterfaceCallback::discardMatchSurfaces()
{
#ifdef SURFACE_VISUALIZER
	if (objectSurfaceVisualizer == NULL) return;
	userInterface->showMatchSurfacesButton->value(0);
	toggleMatchSurfaceDisplay();
	userInterface->matchSurfaceSerialNumber->maximum(0);
	objectSurfaceVisualizer->clear();
	resetMatchSurfaceSerialNumber();
#endif	/* SURFACE_VISUALIZER */
}

void P3DUserInterfaceCallback::importSurfaceColorMap() 
{
#ifdef SURFACE_VISUALIZER
	const char * filename = askSingleFilename(ModelDirectory, "Load Surface Colormap File", "*.txt");

	if (filename)
		surfaceColorMap.loadFromFile(filename);
#endif /* SURFACE_VISUALIZER */
}

#ifndef PRODUCTION_VERSION

void P3DUserInterfaceCallback::exportSeuratUVPhi()
{
	const char * filename = askSingleFilename(ModelDirectory, "Export U V Phi coords", "{*.txt}|*.txt");

	P3DControl * control = getControl();
	if (control) {
		if (control->getObjectPtr() && filename) {
		    std::ofstream out(filename);
		    for (int figId = 0; figId < control->getObjectPtr()->getFigureCount(); figId++) {
			M3DFigure * fig = control->getObjectPtr()->getFigurePtr(figId);
			if (fig) {

				Xferlist * xfer = convertM3DtoXfer(fig);

				if (xfer) {
					ThallCode::Pointlist_server2* pList = new ThallCode::Pointlist_server2();
					pList->init(xfer);
					pList->ComputeSubdivPointCloud(userInterface->modelView->getSurfaceLevel()); // subdiv level
					delete [] (xfer->atomlist);
					delete xfer;

					Bpoint* bpList = 0;
					int numPoints;
					pList->subdivboundaryinfo(&numPoints, &bpList);

					for (int i = 0; i < numPoints; i++) {
						out << figId << " " << bpList[i].u << " " << bpList[i].v << " " << bpList[i].t << " ";

						M3DQuadFigure* quad_fig = dynamic_cast<M3DQuadFigure*>(fig);
                        // print the atom ID
						if (quad_fig != 0) {
						 out << (quad_fig->indexOfAtom((int) bpList[i].v, (int) bpList[i].u));
						}
						else {
						 // we have a tube
						 out << (bpList[i].u);
						}

                        // print the spoke ID
                        out << " " <<  i << endl; 
					}
				}
			}
		}
		out.close();
		}
	}
}

#else	/* PRODUCTION_VERSION */

void P3DUserInterfaceCallback::exportSeuratUVPhi()
{
}

#endif	/* PRODUCTION_VERSION */

