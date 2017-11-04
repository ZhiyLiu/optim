#include <iostream>
#include <vector>
#include <FL/Fl_Tooltip.H>
#include "M3DObjectRenderer.h"
#define D_POINTLIST_SERVER2
#define D_DIATOMGRID
#define D_XFERLIST
#define D_POINTLIST_SERVERB
#include "Shapedepend.h"
#include "M3DObjectSurfaceRenderer.h"
#include "M3DBlendedRenderer.h"
#include "P3DUserInterface.h"
#include "P3DUserInterfaceCallback.h"
#include "utility.h"
#include "BYU.h"
#include "TileSet.h"
#include "TileSetRenderer.h"
#include "ControlParms.h"
#include "CutPlanes.h"
#include "P3DView.h"
#include "M3DObjectSurfaceVisualizer.h"
#ifdef BINARY
#include "DistanceVectorList.h"
#endif
#ifdef COLOR_CONTOURS
#include "SurfaceColorMap.h"
#endif

//#define DEBUG
//#define DRAWING_DEBUG


#define MAX_NUM_LANDMARKS	12	// Maximum no. of image landmarks; implied limit for surface landmarks 
#define MAX_NUM_LANDMARK_HITS	16	// Maximum no. of surface tile intersections per ray
#define LANDMARK_SIZE	0.015	// Model-space line length used to draw landmarks
#define LANDMARK_SURFACE_LEVEL	4	// Feedback-rendering level for landmark computation

#define FRUSTRUM_DEPTH	1000		// Half-depth of the viewing frustrum

// Orange is the default color of landmarks.   See P3DUserInterfaceCallback.cpp.
const float DEFAULT_LANDMARK_COLOR[3] = { 1.0f, 165.0f/255.0f, 0.0f };


using namespace std;


// Non-class function
void checkBlendingFunction()
{
	extern P3DUserInterface * ui;
	ui->modelView->checkToStopBlending();
}


P3DView::P3DView(int x, int y, int w, int h, ViewType _viewType, char *l)
            : Fl_Gl_Window(x, y, w, h, l)
{
    control = NULL;
    callback = NULL;

	twoLights = globalControl->readBool(TwoLights);
    primitiveVectorsOn = false;
    medialSheetOn = false;
	atomsOn = true;
    stdAxesOn = false;
    unitCubeOn = false;
    xSliceOn = false;
    ySliceOn = false;
    zSliceOn = false;
    landmarksOn = false;
	drawBoundary = false;
	cutPlaneBoundaryWidth = globalControl->readInt(CutPlaneBoundaryWidth);
//    arbSliceOn = false;

	labelFigures = true;	// Sync with default value of drawFigureNamesBtn
    tileSetOn = true;	// Sync with default value of checkButton
	defaultTilesColor[0] = globalControl->readDouble(TilesColor_R);
	defaultTilesColor[1] = globalControl->readDouble(TilesColor_G);
	defaultTilesColor[2] = globalControl->readDouble(TilesColor_B);
	defaultTilesColorChanged = true;
	tileColorsChanged = false;
	currentFigureCount = 0;
	tileColors = NULL;
	useDfltTileSetColor = false;
	renderTileSetSolid = false;
	tileRenderer = NULL;

	rockingOn = false;
	constraintArrows = false;	// Equals initial value of userInterface->constraintsVisible

	bNCutPlaneOn = false;
	bBperpCutPlaneOn = false;
	bperpNCutPlaneOn = false;
	bperpY0CutPlaneOn = false;
	bperpY1CutPlaneOn = false;
	involutesCutPlaneOn = false;
	texture_done = false;

	tileType = GL_SMOOTH;
	blendedRendererRunning = false;	// Needed for setSurfaceLevel()

    surfaceStyle = (M3DSurfaceStyle)
		globalControl->read(SurfaceStyle, (int) M3D_NONE);
	if (surfaceStyle < M3D_NONE || surfaceStyle > M3D_SOLID)
		surfaceStyle = M3D_NONE;
	setSurfaceLevel(globalControl->read(SurfaceLevel, 1));
    surfaceLineWidth = (M3DSurfaceStyle)
		globalControl->read(SurfaceLineWidth, 1);

    partialSurfaceStyle = (M3DSurfaceStyle)
		globalControl->read(PartialSurfaceStyle, (int) M3D_NONE);
	if (partialSurfaceStyle < M3D_NONE || partialSurfaceStyle > M3D_SOLID)
		partialSurfaceStyle = M3D_NONE;
	setPartialSurfaceLevel(globalControl->read(PartialSurfaceLevel, 1));

	auxTileSetOn = false;
	auxTileSetDone = false;
	blendSurface = false;
	surfaceDone = false;
	blendingLevel = 0;
	surfaceLevelChanged = false;
	surfaceBlendedLevelChanged = false;
	blendedSurfaceRenderer.setInterruptCallback(checkBlendingFunction);
	interruptOccurred = false;

#ifndef PRODUCTION_VERSION
	colorSurface = false;
	outlineSurfaceTiles = true;
	surfaceTiedToImageSlider = false;
	surfaceColorMin = 0.0;
	surfaceColorMax = 1.0;
	colorSurfaceStyle = M3D_COLORED_SOLID;
#endif

    tileAlpha = 1.0f;
	tileSetDone = true;

	currentAngle = 0.0;	
	deflectionAngle = globalControl->read(RockingAngle, 20.0);
	incrementalAngle = globalControl->read(RockingIncrement, 0.5);
	direction = false;

    mouseDown = false;

    selectionType = SELECT_PRIMITIVE;
    objectDone = true;
	modelIsBeingEdited = false;
	modelWasEdited = false;

    isDrawing = false;
	damaged = Redraw;
	must_redraw = false;

    trackball.setViewDims(VIEW_WIDTH, VIEW_HEIGHT);
    resetTrackball();

    viewType = _viewType;
/*  if(viewType == XY_ORTHO || viewType == YZ_ORTHO || viewType == XZ_ORTHO)
    {
        trackball.setRotationType(ROTATE_2D);
        trackball.disableRotation();
    }*/

    init = false;

	cutPlaneTextureNames[0] = 0;

	// Read model window background color.  Default is INDIGO, as
	// used in file P3DUserInterfaceCallback.cpp.
	backgroundColor[0] = (float) globalControl->read(BackgroundColor_R, 0.05);
	backgroundColor[1] = (float) globalControl->read(BackgroundColor_G, 0.05);
	backgroundColor[2] = (float) globalControl->read(BackgroundColor_B, 0.1);

	// Read the landmarks color.  The default is rather ORANGE.
	defaultLandmarkColor[0] = (float) globalControl->read(LandmarksColor_R, DEFAULT_LANDMARK_COLOR[0]);
	defaultLandmarkColor[1] = (float) globalControl->read(LandmarksColor_G, DEFAULT_LANDMARK_COLOR[1]);
	defaultLandmarkColor[2] = (float) globalControl->read(LandmarksColor_B, DEFAULT_LANDMARK_COLOR[2]);

	// Initialize the landmark widths.  In the control file they are
	// integers; in this class they are floats for OpenGL.
	int j = globalControl->read(LandmarkNarrowWidth, 1);
	defaultLandmarkWidth = (float) j;
	j = globalControl->read(LandmarkWideWidth, 3);
	markedLandmarkWidth = (float) j;

	// Landmarking variables
	landmarkSurface = NULL;
	distThreshold = NULL;
	landmarks = NULL;
	landmarkMode = false;

    dumpScreen = false;
	screenDumpScheduled = 0;
#ifndef PRODUCTION_VERSION
	displayBlendRegion = false;
#endif

#ifdef BINARY
	distanceVectorList = NULL;
	distanceVectorsDone = true;
#endif
	tileSetMobile = false;
	
	centerMode = CEN_MODEL | CEN_ALL;
}

P3DView::~P3DView()
{
#ifdef DEBUG
	cout << "P3DView::~P3DView()" << endl;
#endif
	if (landmarkSurface != NULL)
		delete landmarkSurface;
	if (distThreshold != NULL)
		delete distThreshold;
	delete [] landmarks;
	delete [] tileColors;
	delete tileRenderer;
	// distanceVectorList does not belong to this class
}

void P3DView::setPartialSurfaceLevel(int level) {
	// M3DSurfaceObjectRenderer requires 1 <= partialSurfaceLevel <= 4
	if (level < 0)
		partialSurfaceLevel = 0;
	else if (level > 3)
		partialSurfaceLevel = 3;
	else partialSurfaceLevel = level;
	surfaceDone = false;
}

void P3DView::setSurfaceLevel(int level) {
	// M3DSurfaceObjectRenderer requires 1 <= surfaceLevel <= 4
	if (level < 0)
		surfaceLevel = 0;
	else if (level > 3)
		surfaceLevel = 3;
	else surfaceLevel = level;
	surfaceDone = false;
	surfaceLevelChanged = true;	// Used for drawing surface boundaries in slice planes
	if (blendedRendererRunning)
		surfaceBlendedLevelChanged = true;		// Used to interrupt blending
	auxTileSetDone = false;
}

void P3DView::toggleTileType(bool val)
{
	tileType = (val ? GL_FLAT : GL_SMOOTH);
	updateSurface();
}

void P3DView::toggleBlendSubfigures(bool val) {
	blendSurface = val;
	surfaceDone = false;
}

#ifndef PRODUCTION_VERSION

void P3DView::toggleSurfaceColoring(bool onOff)
{
	colorSurface = onOff;
	surfaceDone = false;
	updateSurface();
}

void P3DView::tieSurfaceToImageSlider(bool yesno)
{
	surfaceTiedToImageSlider = yesno;
	surfaceDone = false;
	updateSurface();
}

void P3DView::setSurfaceColorRange(float min, float max) {
	surfaceColorMin = min;
	surfaceColorMax = max;
	if (! surfaceTiedToImageSlider) {
		surfaceDone = false;
		updateSurface();
	}
}

void P3DView::getSurfaceIntensityRange(float & min, float & max)
{
	if (surfaceTiedToImageSlider) {
		if (control != NULL) {
			double lo, hi;

			control->getIntensityWindow(lo, hi);
			if (lo > hi) {
				lo = 0.0;
				hi = 1.0;
			}
			min = (float) lo;
			max = (float) hi;
		}
		else {
			min = 0.0f;
			max = 1.0f;
		}
	}
	else {
		min = surfaceColorMin;
		max = surfaceColorMax;
	}
}

#endif	/* PRODUCTION_VERSION */

void P3DView::setControl(P3DControl * cntrl)
{
    control = cntrl;
    trackball.setControl(cntrl);
    objectDone = false;
}

void P3DView::setViewType(ViewType type)
{
/*  if(type == XY_ORTHO || type == YZ_ORTHO || type == XZ_ORTHO)
    {
        trackball.setRotationType(ROTATE_2D);
        trackball.disableRotation();
    }
*/
    viewType = type;
}

void P3DView::setAtomVectorsType(bool toggle) {
	M3DObjectRenderer::setAtomVectorsType(toggle);
}

void P3DView::setBVectorsType(int val) {
	M3DObjectRenderer::setBVectorsType(val);
}

void P3DView::setShowExtraAtomVectors(bool toggle) {
	M3DObjectRenderer::setExtraAtomVectors(toggle);
}

void P3DView::setFiguresColorAtoms(bool toggle) {
	M3DObjectRenderer::setUseFigureColor((float) toggle);
}

void P3DView::setAtomVectorsWidth(int width) {
	M3DObjectRenderer::setAtomVectorsLineWidth((float) width);
}

void P3DView::setShowMeshConnectors(bool toggle) {
	M3DObjectRenderer::setMeshVisibility(toggle);
}

void P3DView::setMeshConnectorsType(bool toggle) {
	M3DObjectRenderer::setMeshConnectorsType(toggle);
}

void P3DView::setMeshConnectorsWidth(int width) {
	M3DObjectRenderer::setMeshConnectorsLineWidth(width);
}

void P3DView::setMeshConnectorsColor(float color[3]) {
	M3DObjectRenderer::setMeshConnectorsColor(color);
	updateObject();
}

void P3DView::setBackgroundColor(float color[3]) {
	backgroundColor[0] = color[0];
	backgroundColor[1] = color[1];
	backgroundColor[2] = color[2];
}

void P3DView::setLandmarksColor(const float color[3]) {
	defaultLandmarkColor[0] = color[0];
	defaultLandmarkColor[1] = color[1];
	defaultLandmarkColor[2] = color[2];
	updateSurface();
}

void P3DView::setLandmarksNarrowWidth(int width) {
	defaultLandmarkWidth = (float) width;
}

void P3DView::setLandmarksWideWidth(int width) {
	markedLandmarkWidth = (float) width;
}

// Reset the trackball to zero rotation
void P3DView::resetTrackball()
{
    if(trackball.getCurrScale() != 1.0f)
        objectDone = false;

    transformMatrix[0][0] = 1.0f;
    transformMatrix[0][1] = 0.0f;
    transformMatrix[0][2] = 0.0f;
    transformMatrix[0][3] = 0.0f;
    transformMatrix[1][0] = 0.0f;
    transformMatrix[1][1] = 1.0f;
    transformMatrix[1][2] = 0.0f;
    transformMatrix[1][3] = 0.0f;
    transformMatrix[2][0] = 0.0f;
    transformMatrix[2][1] = 0.0f;
    transformMatrix[2][2] = 1.0f;
    transformMatrix[2][3] = 0.0f;
    transformMatrix[3][0] = (float) -(STD_AXES_LENGTH / 2.0);
    transformMatrix[3][1] = (float) -(STD_AXES_LENGTH / 2.0);
    transformMatrix[3][2] = (float) -(STD_AXES_LENGTH / 2.0);
    transformMatrix[3][3] = 1.0f;

    trackball.reset(&(transformMatrix[0][0]));
}

void P3DView::setAxialView()
{
    if(trackball.getCurrScale() != 1.0f)
        objectDone = false;

    transformMatrix[0][0] = -2.0f;
    transformMatrix[0][1] = 0.0f;
    transformMatrix[0][2] = 0.0f;
    transformMatrix[0][3] = 0.0f;
    transformMatrix[1][0] = 0.0f;
    transformMatrix[1][1] = 2.0f;
    transformMatrix[1][2] = 0.0f;
    transformMatrix[1][3] = 0.0f;
    transformMatrix[2][0] = 0.0f;
    transformMatrix[2][1] = 0.0f;
    transformMatrix[2][2] = -2.0f;
    transformMatrix[2][3] = 0.0f;
    transformMatrix[3][0] = (float) STD_AXES_LENGTH;
    transformMatrix[3][1] = (float) -STD_AXES_LENGTH;
    transformMatrix[3][2] = (float) STD_AXES_LENGTH;
    transformMatrix[3][3] = 1.0f;

    trackball.reset(&(transformMatrix[0][0]));
}

void P3DView::setSaggitalView()
{
    if(trackball.getCurrScale() != 1.0f)
        objectDone = false;

    transformMatrix[0][0] = 0.0f;
    transformMatrix[0][1] = 0.0f;
    transformMatrix[0][2] = -2.0f;
    transformMatrix[0][3] = 0.0f;
    transformMatrix[1][0] = 2.0f;
    transformMatrix[1][1] = 0.0f;
    transformMatrix[1][2] = 0.0f;
    transformMatrix[1][3] = 0.0f;
    transformMatrix[2][0] = 0.0f;
    transformMatrix[2][1] = -2.0f;
    transformMatrix[2][2] = 0.0f;
    transformMatrix[2][3] = 0.0f;
    transformMatrix[3][0] = (float) -STD_AXES_LENGTH;
    transformMatrix[3][1] = (float) STD_AXES_LENGTH;
    transformMatrix[3][2] = (float) STD_AXES_LENGTH;
    transformMatrix[3][3] = 1.0f;

    trackball.reset(&(transformMatrix[0][0]));
}

void P3DView::setCoronalView()
{
    if(trackball.getCurrScale() != 1.0f)
        objectDone = false;

    transformMatrix[0][0] = 2.0f;
    transformMatrix[0][1] = 0.0f;
    transformMatrix[0][2] = 0.0f;
    transformMatrix[0][3] = 0.0f;
    transformMatrix[1][0] = 0.0f;
    transformMatrix[1][1] = 0.0f;
    transformMatrix[1][2] = 2.0f;
    transformMatrix[1][3] = 0.0f;
    transformMatrix[2][0] = 0.0f;
    transformMatrix[2][1] = -2.0f;
    transformMatrix[2][2] = 0.0f;
    transformMatrix[2][3] = 0.0f;
    transformMatrix[3][0] = (float) -STD_AXES_LENGTH;
    transformMatrix[3][1] = (float) STD_AXES_LENGTH;
    transformMatrix[3][2] = (float) -STD_AXES_LENGTH;
    transformMatrix[3][3] = 1.0f;

    trackball.reset(&(transformMatrix[0][0]));
}

void P3DView::rotateLeft90()
{
	Vector3D center;

	getCenter(center);
	Vector3D negCenter = -center;
	translateToCenter(negCenter);
    rotateAboutCenter(90.0);
	translateToCenter(center);
}

void P3DView::rotateRight90()
{
	Vector3D center;

	getCenter(center);
	Vector3D negCenter = -center;
	translateToCenter(negCenter);
    rotateAboutCenter(-90.0);
	translateToCenter(center);
}

void P3DView::rotateHorizontal180()
{
	Vector3D center;

	getCenter(center);
	Vector3D negCenter = -center;
	translateToCenter(negCenter);
    rotateAboutVertical(180.0);
	translateToCenter(center);
}

void P3DView::rotateVertical180()
{
	Vector3D center;

	getCenter(center);
	Vector3D negCenter = -center;
	translateToCenter(negCenter);
    rotateAboutHorizontal(180.0);
	translateToCenter(center);
}

// find something upon which to center the view.
// arg 'which' determines what to center upon:
//  0: use last 'which'
//  MODEL & ALL: entire model
//  MODEL & FIG: a specified figure FIG
//  MODEL & MARK: marked atom
//  MODEL & SELECT: all selected atoms
//  TILE & ALL: entire tileset
//  TILE & FIG: a specified tile figure FIG
//  SPACE & ALL: unit cube
//  ILM: a spec'd image LANDMARK in a spec'd FIG

void P3DView::setCenterOfRotation(int which)
{
	if (!control)
		return;
	if (!which)		// use last centerMode if called has no preference
	{
		if (!centerMode) return;
		which = centerMode;
	}
	else
		centerMode = which;

	bool foundCenter = false;
	Vector3D posn(0,0,0);
	M3DObject * object;
	M3DFigure * figure;

	if (which & CEN_ILM)
	{
		object = control->getObjectPtr();
		//if (which & CEN_FIG)
		{
#ifdef LM_METHOD_OBJ
			int figNum = (which & CEN_FIG_MASK) >> CEN_FIG_SHIFT;
			int ilmNum = (which & CEN_ILM_MASK) >> CEN_ILM_SHIFT;
			figure = object->getFigurePtr(figNum);
			int atomNum = figure->getLandmarkAtomIndex(ilmNum);
			M3DPrimitive *atom = figure->getPrimitivePtr(atomNum);
			
			int tipT = (int) figure->getLandmarkAtomT(ilmNum);	// AGG: Should this round?
			
			Vector3D y0 = atom->getY0();
			Vector3D bi = atom->getExtendedB();
			Vector3D y1 = atom->getY1();
			if (tipT ==  1.0) posn = y1;
			else if (tipT == -1.0) posn = y0;
			else if (tipT ==  0.0) posn = bi;
			else if (tipT < 0)    posn = bi + -tipT * (y0 - bi); 
			else                  posn = bi +  tipT * (y1 - bi);
			
			posn = atom->getX() + posn;
			foundCenter = true;
			cout << "Centering on image landmark (" << figure->getLandmarkName(ilmNum)
				<< ") of figure #" << ilmNum << ": ";
#endif
		}
		// LATER: all ilm's
	}

	if (which & CEN_MODEL)
	{
		object = control->getObjectPtr();
		if (!object) return;
		int markedFigureId, primId;
		int nFigures;
		
		if (which & CEN_ALL)
		{
			posn = object->getCOG(true);
			foundCenter = true;
			cout << "Centering on COG of whole model: ";
		}
		
		else if (which & CEN_SELECT)
		{
			nFigures = object->getFigureCount();
			int n = 0;
			for (int j = 0; j < nFigures; j++) {
				figure = object->getFigurePtr(j);
				if (figure->isSelected()) {
					posn += figure->getCOG(true);
					n++;
				}
			}
			if (n > 0)
			{
				posn /= n;
				foundCenter = true;
				cout << "Centering on COG of " << n << " selected figure(s): ";
			}
		}

		else if (which & CEN_MARK)
		{
			markedFigureId = -1;
			control->getRelativeMarkedPrimitiveId(markedFigureId, primId);
			if (markedFigureId >= 0) {
				M3DPrimitive * atom = object->getPrimitivePtr(markedFigureId, primId);
				posn = atom->getX();
				foundCenter = true;
				cout << "CENTERING ON marked atom: ";
			}
		}
		
		else if (which & CEN_FIG)
		{
			int figNum = which & CEN_FIG_MASK;
			figure = object->getFigurePtr(figNum);
			posn = figure->getCOG(true);
			foundCenter = true;
			cout << "Centering on COG of figure #" << figNum << ": ";
		}
	}
	else if (which & CEN_TILE)
	{
		TileSet * tileSet = control->getTileSetPtr();
		if (!tileSet || !tileSet->tileCount())
			return;
		if (which & CEN_ALL)
		{
			tileSet->calcCOG();
			posn = tileSet->getCOG();
			foundCenter = true;
			cout << "Centering on tile set: ";
		}
		else if (which & CEN_FIG)
		{
			int figNum = which & CEN_FIG_MASK;
			bool figs[6] = {0,0,0,0,0,0};
			figs[figNum] = true;
			tileSet->calcCOG(figs);
			
			posn = tileSet->getCOG();
			foundCenter = true;
			cout << "Centering on COG of tile figure " << figNum << ": ";
		}
	}
	else if (which & CEN_SPACE)
	{
		if (which & CEN_ALL)
		{
			posn.set(.5, .5, .5);
			foundCenter = true;
			cout << "CENTERING ON unit cube: ";
		}
		// LATER: corners/faces of unit go here (prob'ly not useful)
	}
	
	if (foundCenter)
	{
		posn.print();
		translateToCenter(posn);
	}
}

// Handle mouse events
int P3DView::handle(int eventId)
{
    int button;
    MouseEventStruct mouse;
    char * dropFile,
         * ext;

    switch(eventId)
    {
        // Allow Drag-and-Drop of m3d and raw3 files
        case FL_DND_ENTER:
        case FL_DND_DRAG:
        case FL_DND_RELEASE:
            return 1;

        case FL_PASTE:
            dropFile = (char *) Fl::event_text();
            ext = strrchr(dropFile, '.');
            if(ext != NULL)
            {
                if(strcmp(ext, ".m3d") == 0)
                {
                    callback->loadModel(dropFile);
                    return 1;
                }
                else if(strcmp(ext, ".raw3") == 0)
                {
                    callback->loadImage(dropFile);
                    return 1;
                }
            }
            break;

        case FL_PUSH:
            button = Fl::event_button();

            mouseDown = false;

            // Left mouse click w/ Ctrl key down
            // -- selects a primitive/slice/figure
            if(button == 1 && Fl::event_state(FL_CTRL) && (! Fl::event_state(FL_SHIFT))
				&& (! getLandmarkMode()))
            {
                mouseDown = true;
                lastMouseX = currMouseX = Fl::event_x();
                lastMouseY = currMouseY = Fl::event_y();
                return 1;
            }

            // Right mouse click w/ Ctrl key down -- marks an origin primitive
            else if(button == 3 && Fl::event_state(FL_CTRL) && (! Fl::event_state(FL_SHIFT))
				&& (! getLandmarkMode()))
            {
                lastMouseX = currMouseX = Fl::event_x();
                lastMouseY = currMouseY = Fl::event_y();

                mouseSelect(true);

                redraw();
            }

            else
            {
                mouse.button = button;
                mouse.width = w();
                mouse.height = h();
                mouse.type = MOUSE_PUSH;
                mouse.x = Fl::event_x();
                mouse.y = Fl::event_y();
                mouse.viewType = viewType;

                mouse.modifiers = NO_KEY;
                if(Fl::event_state(FL_CTRL))
                    mouse.modifiers |= CTRL_KEY;
                if(Fl::event_state(FL_SHIFT))
                    mouse.modifiers |= SHIFT_KEY;
                if(Fl::event_state(FL_ALT))
                    mouse.modifiers |= ALT_KEY;

                trackball.mouseEvent(mouse);

				if (trackball.getModeClass() == true)
					modelIsBeingEdited = true;

                return 1;
            }

            break;

        // Used to draw a box to select a set of primitives
        case FL_DRAG:
            if(isDrawing)
                return 0;

            currMouseX = Fl::event_x();
            currMouseY = Fl::event_y();

            if(trackball.getMode() != NONE)
            {
                mouse.button = Fl::event_button();
                mouse.width = w();
                mouse.height = h();
                mouse.type = MOUSE_DRAG;
                mouse.x = currMouseX;
                mouse.y = currMouseY;
                mouse.viewType = viewType;

                mouse.modifiers = NO_KEY;
                if(Fl::event_state(FL_CTRL))
                    mouse.modifiers |= CTRL_KEY;
                if(Fl::event_state(FL_SHIFT))
                    mouse.modifiers |= SHIFT_KEY;
                if(Fl::event_state(FL_ALT))
                    mouse.modifiers |= ALT_KEY;

                trackball.mouseEvent(mouse);
            }

			if (tileSetMobile)
			{
				setTileFileChanged();
				update();
			}
			else
			{
				updateObject();	// does not move tiles
			}

            break;

        case FL_RELEASE:
            button = Fl::event_button();

            // User released the left mouse -- do selection
            if(button == 1 && trackball.getMode() == NONE && mouseDown)
                mouseSelect(false);

            else if(button == 1 && getLandmarkMode()
				&& Fl::event_state(FL_CTRL) && (! Fl::event_state(FL_SHIFT)))
			{
                currMouseX = Fl::event_x();
                currMouseY = Fl::event_y();
                mouseSelectTile();	// Place landmark on figural surface
			}

			else if(button == 3 && getLandmarkMode()
				&& Fl::event_state(FL_CTRL) && (! Fl::event_state(FL_SHIFT)))
			{
                currMouseX = Fl::event_x();
                currMouseY = Fl::event_y();
                mouseSelectImage();	// Place landmark on image slice
			}

            else if(trackball.getMode() != NONE)
            {
                mouse.button = button;
                mouse.width = w();
                mouse.height = h();
                mouse.type = MOUSE_RELEASE;
                mouse.x = Fl::event_x();
                mouse.y = Fl::event_y();
                mouse.viewType = viewType;

                mouse.modifiers = NO_KEY;
                if(Fl::event_state(FL_CTRL))
                    mouse.modifiers |= CTRL_KEY;
                if(Fl::event_state(FL_SHIFT))
                    mouse.modifiers |= SHIFT_KEY;
                if(Fl::event_state(FL_ALT))
                    mouse.modifiers |= ALT_KEY;
                if(Fl::event_state(FL_META))
                    mouse.modifiers |= META_KEY;

                trackball.mouseEvent(mouse);

				modelWasEdited = modelIsBeingEdited;
				modelIsBeingEdited = false;
            }

            mouseDown = false;

			if (modelWasEdited)
			{
				if (tileSetMobile)
				{
					setTileFileChanged();
					update();
				}
				else
				{
					updateObject();	// does not move tiles
				}
			}
			else
				redraw();

            break;

		case FL_SHORTCUT:	// Key pressed
#ifdef sun
			if ((unsigned short) Fl::event_key() == FL_F - 172) {		// F12 key
#else
			if (Fl::event_key() == FL_F + 12) {		// F12 key
#endif
				if (insideModelWindow())
					installBubblePopUp();
				return 1;
			}
			break;

		case FL_KEYUP:		// Key released
#ifdef sun
			if ((unsigned short) Fl::event_key() == FL_F - 172) {	// F12 key
#else
			if (Fl::event_key() == FL_F + 12) {		// F12 key
#endif
				removeBubblePopUp();
				return 1;
			}
			break;

        default:
            break;
    };
    return Fl_Gl_Window::handle(eventId);
}

// Local debugging routine
void glReport() {
	float p[4];

	cout << "glReport:" << '\n';
	glGetFloatv(GL_CURRENT_COLOR, p);
	cout << "\tcolor = (" << p[0] << ", " << p[1] << ", " << p[2] << ')' << endl;
	glGetFloatv(GL_LINE_WIDTH, p);
	cout << "\tline width = " << p[0] << endl;
}

// Draw a box when the mouse is dragged
// Used to select a set of primitives
void P3DView::drawMouseBox(int x1, int y1, int x2, int y2)
{
    double dx1, dx2, dy1, dy2;

#ifdef DEBUG
	cout << "P3DView::drawMouseBox()" << endl;
#endif
    dx1 = x1 * ((double) VIEW_WIDTH / (double) w()) - (VIEW_WIDTH / 2.0);
    dx2 = x2 * ((double) VIEW_WIDTH / (double) w()) - (VIEW_WIDTH / 2.0);
    dy1 = -y1 * ((double) VIEW_HEIGHT / (double) h()) + (VIEW_HEIGHT / 2.0);
    dy2 = -y2 * ((double) VIEW_HEIGHT / (double) h()) + (VIEW_HEIGHT / 2.0);

    glEnable(GL_LINE_STIPPLE);
    glLineStipple(3,0x5555);

    glColor3f(1.0, 1.0, 0.0);

    glBegin(GL_LINE_LOOP);

    glVertex3d(dx1, dy1, FRUSTRUM_DEPTH/2.0);
    glVertex3d(dx2, dy1, FRUSTRUM_DEPTH/2.0);
    glVertex3d(dx2, dy2, FRUSTRUM_DEPTH/2.0);
    glVertex3d(dx1, dy2, FRUSTRUM_DEPTH/2.0);

    glEnd();

    glLineStipple(3,0xAAAA);

    glColor3f(0.0, 0.0, 0.0);

    glBegin(GL_LINE_LOOP);

    glVertex3d(dx1, dy1, FRUSTRUM_DEPTH/2.0);
    glVertex3d(dx2, dy1, FRUSTRUM_DEPTH/2.0);
    glVertex3d(dx2, dy2, FRUSTRUM_DEPTH/2.0);
    glVertex3d(dx1, dy2, FRUSTRUM_DEPTH/2.0);

    glEnd();

    glDisable(GL_LINE_STIPPLE);
}

// This function determines which if two primitives is front-most
// with respect to the user's viewpoint.
bool P3DView::in_front(int testId, int primitiveId)
{
	double x, y;
	double z1, z2;	// Screen coordinates
	M3DPrimitive * atom;
	Vector3D coords;
	GLint viewport[4];
	GLdouble modelview[16];
	GLdouble projection[16];

    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glLoadIdentity();

    trackball.doViewTransform();

    // Convert from pixels to model units
	glGetIntegerv(GL_VIEWPORT, viewport);
	glGetDoublev(GL_MODELVIEW_MATRIX, modelview);
	glGetDoublev(GL_PROJECTION_MATRIX, projection);
    glPopMatrix();

	// Determine which primitive selection is front-most
    M3DObject * object = control->getObjectPtr();
	atom = object->getPrimitivePtr(testId);
	coords = atom->getX();

	gluProject(coords.getX(), coords.getY(), coords.getZ(),
		modelview, projection, viewport, &x, &y, &z1);

	atom = object->getPrimitivePtr(primitiveId);
	coords = atom->getX();
	gluProject(coords.getX(), coords.getY(), coords.getZ(),
		modelview, projection, viewport, &x, &y, &z2);

	if (z1 < z2)
		return true;
	else
		return false;
}

// Manually select or mark primitives indicated by clicking the mouse buttons.
// This function is also used to determine the figure and atom numbers
// reported in the model view's tooltip when F12 is pressed.
void P3DView::mouseSelect(bool chooseOrigin)
{
    double xConvFactor, yConvFactor;
    double dCurrX, dCurrY, dLastX, dLastY;
    double delta;
	extern P3DUserInterface * ui;


	glRenderMode(GL_RENDER);
    if(control == NULL)
        return;
    M3DObject * object = control->getObjectPtr();
    if(object == NULL)
        return;

    // Convert from pixels to model units
    GLuint * selectBuf = new GLuint[BUFSIZE];
    memset(selectBuf, -1, BUFSIZE * sizeof(GLuint));

    glSelectBuffer(BUFSIZE, selectBuf);

#ifdef DEBUG
	cout << "pixel box: (" << lastMouseX << ", " << lastMouseY << ") -> ("
		<< currMouseX << ", " << currMouseY << ")\n";
#endif
    xConvFactor = (double) VIEW_WIDTH / (double) w();
    yConvFactor = (double) VIEW_HEIGHT / (double) h();

    dCurrX = ((double) currMouseX * xConvFactor) - (VIEW_WIDTH / 2.0);
    dCurrY = - ((double) currMouseY * yConvFactor) + (VIEW_HEIGHT / 2.0);

    if(currMouseX == lastMouseX)
    {
        delta = (LARGE_POINT_SIZE / 2.0);
        dLastX = dCurrX - delta;
        dCurrX += delta;
    }

    else
        dLastX = ((double) lastMouseX * xConvFactor) - (VIEW_WIDTH / 2.0);

    if(currMouseY == lastMouseY)
    {
        delta = (LARGE_POINT_SIZE / 2.0);
        dLastY = dCurrY - delta;
        dCurrY += delta;
    }

    else
        dLastY = - ((double) lastMouseY * yConvFactor) + (VIEW_HEIGHT / 2.0);
#ifdef DEBUG
	cout << "world box: (" << dLastX << ", " << dLastY << ") -> (" << dCurrX
		<< ", " << dCurrY << ")\n";
#endif

    glRenderMode(GL_SELECT);

    glInitNames();
    glPushName(0);

    glPushMatrix();

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(dLastX, dCurrX, dLastY, dCurrY, -FRUSTRUM_DEPTH, FRUSTRUM_DEPTH);

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

    trackball.doViewTransform();

    drawObject(false);

    glPopMatrix();
    glFlush();

    int figIndex;
    int testFigIndex;
    int primitiveId;

    M3DFigure * fig;

    vector<unsigned int> toggledFigures;
    vector<unsigned int>::iterator it;
    bool toggleFigure;

    int j;
	// Specify the first relevant location in selectBuf
	if (primitiveVectorsOn || medialSheetOn)
		j = 7;
	else
		j = 3;
	if (selectionType == SHOW_ID || chooseOrigin) {
    	idTooltip[0] = '\0';
    	primitiveId = -1;
    	while(selectBuf[j] != -1 && j < BUFSIZE)  // -1 indicates an empty selectBuf slot
    	{
    		if(selectBuf[j] == -2)  // -2 indicates an erroneously filled slot
    		{
    			j += 4;
    			continue;
    		}

    		// Find the front-most atom, when atoms overlap
    		if (primitiveId < 0)
    			primitiveId = selectBuf[j];
    		else {
    			if (in_front(selectBuf[j], primitiveId))
    				primitiveId = selectBuf[j];
    		}

    		// Every fourth slot in the selectBuf contains a primitive number
    		// The other slots contain unused information
    		j += 4;
    	}

		if(chooseOrigin)
		{
			control->toggleMarkedPrimitive(primitiveId);
			if (ui->primitiveEditorWindow->shown())
				callback->updateAtomEditorWindow();
			if (ui->elongationDialog->shown())
				callback->updateElongationWindow();
		}
    	else if (primitiveId >= 0)	// SHOW_ID
		{
    		// Display the figure and atom Id's in the tooltip 
    		figIndex = 0;
    		fig = object->getFigurePtr(figIndex);
    		while(fig != NULL && primitiveId >= fig->getPrimitiveCount())
    		{
    			primitiveId -= fig->getPrimitiveCount();
    			figIndex++;
    			fig = object->getFigurePtr(figIndex);
    		}

    		sprintf(idTooltip, "fig%d, hub%d", figIndex, primitiveId);
    	}
    }
    else {
		// Using the mouse to select figure(s) or atom(s)

		// toggle each prim once, even if listed multiple times in selectBuf
		// (which happens due to model landmark tip balls)
		int primCount = object->getPrimitiveCount();
		bool *primList = new bool[primCount];
		memset((void*)primList, 0, sizeof(bool)*primCount);

		if (! control->selectionIsAllowed())
			fl_beep(FL_BEEP_ERROR);
    	else while(selectBuf[j] != -1 && j < BUFSIZE)  // -1 indicates an empty selectBuf slot
    	{
    		if(selectBuf[j] == -2)  // -2 indicates an erroneously filled slot
    		{
    			j += 4;
    			continue;

			}

			primitiveId = selectBuf[j];
			// Note that for atom selection, hidden atoms are not excluded, because
			// selection can involve all atoms inside a drawn box.

			// Select/deselect the primitive that was clicked on
			if(selectionType == SELECT_PRIMITIVE) {
				if (!primList[primitiveId]) {
					control->togglePrimitive(primitiveId);
					primList[primitiveId] = true;
				}
			}

			// Select/deselect the whole figure
			else if(selectionType == SELECT_FIGURE)
			{
				figIndex = 0;

				fig = object->getFigurePtr(figIndex);
				while(fig != NULL && primitiveId >= fig->getPrimitiveCount())
				{
					primitiveId -= fig->getPrimitiveCount();
					figIndex++;
					fig = object->getFigurePtr(figIndex);
				}

				it = toggledFigures.begin();
				toggleFigure = true;
				while(it != toggledFigures.end())
				{
					testFigIndex = (*it);
					if(testFigIndex == figIndex)
					{
						toggleFigure = false;
						break;
					}

					it++;
				}

				if(toggleFigure)
				{
					fig = object->getFigurePtr(figIndex);
					if(fig != NULL)
						control->toggleFigure(figIndex);

					toggledFigures.insert(toggledFigures.end(), figIndex);
				}
			}

			// Every fourth slot in the selectBuf contains a primitive number
			// The other slots contain unused information
			j += 4;
		}
		delete [] primList;
    }
    delete [] selectBuf;

    glRenderMode(GL_RENDER);

    updateObject();
}

// Local debugging routine
void P3DView::printMatrices()
{
	GLint viewport[4];
	GLdouble modelview[16];
	GLdouble projection[16];
	int j;

	glGetIntegerv(GL_VIEWPORT, viewport);
	glGetDoublev(GL_MODELVIEW_MATRIX, modelview);
	glGetDoublev(GL_PROJECTION_MATRIX, projection);
	// Note: OpenGL store matrices in column-major order.  See glLoadMatrix
	// and glMultMatrix in the "OpenGL Reference Manual".
	cout << "Modelview matrix:\n";
	for (j = 0; j < 4; j++) {
		cout << '\t';
		for (int k = 0; k < 3; k++)
			cout << modelview[k+4*j] << ", ";
		cout << modelview[3+4*j] << '\n';
	}
	cout << '\n';
	cout << "Projection matrix:\n";
	for (j = 0; j < 4; j++) {
		cout << '\t';
		for (int k = 0; k < 3; k++)
			cout << projection[k+4*j] << ", ";
		cout << projection[3+4*j] << '\n';
	}
	cout << '\n';
	cout << "Viewport matrix:\n";
	cout << '\t';
	for (j = 0; j < 3; j++)
		cout << viewport[j] << ", ";
	cout << viewport[3] << '\n';
}

bool P3DView::intersectQuad(Point & pt0, Point & pt1, Point & pt2, Point & pt3, Point & a, Point & b)
{
	Point q = intersect_line_with_plane(a, b, pt0, pt1, pt2);
	if (inside_triangle(q, pt0, pt1, pt2)) {
		landmarks[numLandmarkHits] = q;
		if (numLandmarkHits < MAX_NUM_LANDMARK_HITS - 1)
			numLandmarkHits++;
		else {
			cout << "Too many possible landmarks at the selected coordinates" << endl;
			fl_beep(FL_BEEP_ERROR);
		}
		return true;
	}
	else {
		q = intersect_line_with_plane(a, b, pt0, pt2, pt3);
		if (inside_triangle(q, pt0, pt2, pt3)) {
			landmarks[numLandmarkHits] = q;
			if (numLandmarkHits < MAX_NUM_LANDMARK_HITS - 1)
				numLandmarkHits++;
			else {
				cout << "Too many possible landmarks at the selected coordinates" << endl;
				fl_beep(FL_BEEP_ERROR);
			}
			return true;
		}
	}
	return false;
}

// Manually select model landmarks
void P3DView::mouseSelectTile()
{
	GLint viewport[4];
	GLdouble modelview[16];
	GLdouble projection[16];
	double a[3], b[3];
	int figures[MAX_NUM_LANDMARK_HITS];

#ifdef DEBUG
	cout << "mouseSelectTile()\n";
#endif
    if (control == NULL)
        return;

	if (callback->getLandmarksEditMode() != 0)
		return;

    M3DObject * object = control->getObjectPtr();
    if (object == NULL)
        return;

    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glLoadIdentity();

    trackball.doViewTransform();

    // Convert from pixels to model units
	glGetIntegerv(GL_VIEWPORT, viewport);
	glGetDoublev(GL_MODELVIEW_MATRIX, modelview);
	glGetDoublev(GL_PROJECTION_MATRIX, projection);

#ifdef DEBUG
	printMatrices();
	cout << "Picked " << currMouseX << ", " << currMouseY << endl;
#endif
	int y = viewport[3] - currMouseY - 1;
	gluUnProject((GLdouble) currMouseX, (GLdouble) y, 0.0,
		modelview, projection, viewport, &a[0], &a[1], &a[2]);
	gluUnProject((GLdouble) currMouseX, (GLdouble) y, 1.0,
		modelview, projection, viewport, &b[0], &b[1], &b[2]);
#ifdef DEBUG
	cout << "Intersection line: (" << a[0] << ", " << a[1] << ", " << a[2] << ") - (";
	cout << b[0] << ", " << b[1] << ", " << b[2] << ')' << endl;
#endif

//	float depth;
//	glReadPixels(currMouseX, y, 1, 1, GL_DEPTH_COMPONENT, GL_FLOAT, &depth);
//	cout << "depth = " << depth << endl;

	// Test to find minimum distance to surface
	int ntiles;
	Bpoint * tile;
	bool hit;
	int count;
	int nFigure, vrtx0;
	M3DFigure * figure;

	Point va(a[0], a[1], a[2]);
	Point vb(b[0], b[1], b[2]);

#ifdef DEBUG
	int index = 0;
#endif
	count = 0;
	numLandmarkHits = 0;
	for (nFigure = 0; nFigure < object->getFigureCount(); nFigure++) {
		figure = object->getFigurePtr(nFigure);
		if (figure->getVisibility() == false)	// Skip invisible figures
			continue;

		ThallCode::Pointlist_server2 * pList = landmarkSurface->getPList(nFigure);
		if (pList == NULL)
			continue;
		pList->subdivtileinfo(&ntiles, &tile);
		hit = false;
		for (vrtx0 = 0; vrtx0 < 4*ntiles; vrtx0 += 4) {	// Loop through figural tiles
			Point pt0(tile[vrtx0].pnt[0], tile[vrtx0].pnt[1], tile[vrtx0].pnt[2]);
			Point pt1(tile[vrtx0 + 1].pnt[0], tile[vrtx0 + 1].pnt[1], tile[vrtx0 + 1].pnt[2]);
			Point pt2(tile[vrtx0 + 2].pnt[0], tile[vrtx0 + 2].pnt[1], tile[vrtx0 + 2].pnt[2]);
			Point pt3(tile[vrtx0 + 3].pnt[0], tile[vrtx0 + 3].pnt[1], tile[vrtx0 + 3].pnt[2]);
			double d = pt_to_line_distance(pt0, va, vb);
			if (d < distThreshold[nFigure]) {
				hit = intersectQuad(pt0, pt1, pt2, pt3, va, vb);
				if (hit) {
					figures[count] = nFigure;
					count++;
#ifdef DEBUG
					index = vrtx0/4;
					cout << "Intersecting tile is number " << index << " on figure " << nFigure
						<< " at (" << landmarks[numLandmarkHits].getX() << ", "
						<< landmarks[numLandmarkHits].getY() << ", " << landmarks[numLandmarkHits].getZ()
						<< ')' << endl;
#endif
					continue;
				}
			}

			d = pt_to_line_distance(pt1, va, vb);
			if (d < distThreshold[nFigure]) {
				hit = intersectQuad(pt0, pt1, pt2, pt3, va, vb);
				if (hit) {
					figures[count] = nFigure;
					count++;
#ifdef DEBUG
					index = vrtx0/4;
					cout << "Intersecting tile is number " << index << " on figure " << nFigure
						<< " at (" << landmarks[numLandmarkHits].getX() << ", "
						<< landmarks[numLandmarkHits].getY() << ", " << landmarks[numLandmarkHits].getZ()
						<< ')' << endl;
#endif
					continue;
				}
			}

			d = pt_to_line_distance(pt2, va, vb);
			if (d < distThreshold[nFigure]) {
				hit = intersectQuad(pt0, pt1, pt2, pt3, va, vb);
				if (hit) {
					figures[count] = nFigure;
					count++;
#ifdef DEBUG
					index = vrtx0/4;
					cout << "Intersecting tile is number " << index << " on figure " << nFigure
						<< " at (" << landmarks[numLandmarkHits].getX() << ", "
						<< landmarks[numLandmarkHits].getY() << ", " << landmarks[numLandmarkHits].getZ()
						<< ')' << endl;
#endif
					continue;
				}
			}

			d = pt_to_line_distance(pt3, va, vb);
			if (d < distThreshold[nFigure]) {
				hit = intersectQuad(pt0, pt1, pt2, pt3, va, vb);
				if (hit) {
					figures[count] = nFigure;
					count++;
#ifdef DEBUG
					index = vrtx0/4;
					cout << "Intersecting tile is number " << index << " on figure " << nFigure
						<< " at (" << landmarks[numLandmarkHits].getX() << ", "
						<< landmarks[numLandmarkHits].getY() << ", " << landmarks[numLandmarkHits].getZ()
						<< ')' << endl;
#endif
					continue;
				}
			}
		}
	}

	// Determine which figural intersection is front-most
	if (count > 0) {
		double x, y, z;
		double zmin;	// Screen coordinates
		int infront, j;

		gluProject(landmarks[numLandmarkHits - count].getX(),
			landmarks[numLandmarkHits - count].getY(),
			landmarks[numLandmarkHits - count].getZ(),
			modelview, projection, viewport, &x, &y, &z);
#ifdef DEBUG
		cout << '(' << landmarks[numLandmarkHits - count].getX() << ", "
			<< landmarks[numLandmarkHits - count].getY() << ", "
			<< landmarks[numLandmarkHits - count].getZ() << ") -> ( "
			<< x << ", " << y << ", " << z << ")\n";
#endif
		zmin = z;
		infront = 0;
		for (j = 1; j < count; j++) {
			gluProject(landmarks[numLandmarkHits - count + j].getX(),
				landmarks[numLandmarkHits - count + j].getY(),
				landmarks[numLandmarkHits - count + j].getZ(),
				modelview, projection, viewport, &x, &y, &z);
#ifdef DEBUG
			cout << '(' << landmarks[numLandmarkHits - count + j].getX() << ", "
				<< landmarks[numLandmarkHits - count + j].getY() << ", "
				<< landmarks[numLandmarkHits - count + j].getZ() << ") -> ( "
				<< x << ", " << y << ", " << z << ")\n";
#endif
			if (z < zmin) {
				zmin = z;
				infront = j;
			}
		}

		// Keep only the front-most point (zmin)
#ifdef DEBUG
		cout << "Selected landmark (" << landmarks[infront].getX() << ", "
			<< landmarks[infront].getY() << ", "
			<< landmarks[infront].getZ() << ") on figure "
			<< figures[infront] << endl;
#endif
		figure = object->getFigurePtr(figures[infront]);
		figure->addLandmark(landmarks[infront]);
		control->modelNeedsSaving();

		// Insert the slot for an image landmark in the correct position in the list
		int num = figure->getLandmarkCount() - 1;	// The one just added is always last
		for (j = 0; j < figures[infront]; j++) {
			figure = object->getFigurePtr(j);
			num += figure->getLandmarkCount();
		}
		control->insertImageLandmark(num);
		callback->updateLandmarkEditorWindow();
	}
	else {
		cout << "Not close to any tile" << endl;
		fl_beep(FL_BEEP_ERROR);
	}

    glPopMatrix();
	if (count > 0)
		updateSurface();
}

// Decide whether or not to install tooltip
bool P3DView::insideModelWindow()
{
	extern P3DUserInterface * ui;
	int x, y;
	Fl_Widget * wid;

	wid = Fl::belowmouse();
	if (wid == NULL)
		return false;

    if (control == NULL)
		return false;
    M3DObject * object = control->getObjectPtr();
    if (object == NULL)
		return false;

	Fl::get_mouse(x, y);
	x -= ui->modelWindow->x() + ui->modelView->x();
	y -= ui->modelWindow->y() + ui->modelView->y();

	// All these tests are needed because FLTK sometimes drops events
	if ((wid != ui->modelWindow && wid != ui->modelView) ||
		x < 0 || y < 0 || x >= ui->modelView->w() || y >= ui->modelView->h()) {
		return false;
	}
    return true;
}

void P3DView::installBubblePopUp()
{
	extern P3DUserInterface * ui;

	Fl::get_mouse(mouseTipX, mouseTipY);
	currMouseX = mouseTipX;
	currMouseY = mouseTipY;
	currMouseX -= ui->modelWindow->x() + ui->modelView->x();
	currMouseY -= ui->modelWindow->y() + ui->modelView->y();
	lastMouseX = currMouseX;
	lastMouseY = currMouseY;

	P3DSelectionType selectType;
	selectType = getSelectType();
	setSelectType(SHOW_ID);
	mouseSelect(false);
	setSelectType(selectType);

	int y;
	if (mouseTipY < ui->modelWindow->h()/2)
		y = mouseTipY + 15;
	else
		y = mouseTipY - 65;
	// The Fl_Tooltip class ignores x
	Fl_Tooltip::enter_area(ui->modelWindow, 0, y, 0, 0, idTooltip);
	return;
}

void P3DView::removeBubblePopUp()
{
	extern P3DUserInterface * ui;

	Fl_Tooltip::enter_area(ui->modelWindow, 0, 0, 0, 0, NULL);
	return;
}

// Manually select image landmarks
void P3DView::mouseSelectImage()
{
#ifndef LM_METHOD_OBJ
	
	GLint viewport[4];
	GLdouble modelview[16];
	GLdouble projection[16];
	double a[3], b[3];
    ImagePlanes * imagePlanes;
	int count, infront, j;

#ifdef DEBUG
	cout << "mouseSelectImage()\n";
#endif
    if(control == NULL)
        return;

	if (callback->getLandmarksEditMode() != 1)
		return;

    M3DObject * object = control->getObjectPtr();
    if(object == NULL)
        return;

	int num_landmarks = object->getLandmarkCount();
	if (control->getImageLandmarksCount() >= num_landmarks) {
		if (num_landmarks <= 0)
			cout << "Surface landmarks must be added before image landmarks" << endl;
		else
			cout << "All surface landmarks have matching image landmarks" << endl;
		fl_beep(FL_BEEP_ERROR);
		return;
	}

    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glLoadIdentity();

    trackball.doViewTransform();

    // Convert from pixels to model units
	glGetIntegerv(GL_VIEWPORT, viewport);
	glGetDoublev(GL_MODELVIEW_MATRIX, modelview);
	glGetDoublev(GL_PROJECTION_MATRIX, projection);

	// The mouse pick is considered to cast a ray into the viewport
#ifdef DEBUG
	printMatrices();
	cout << "Picked " << currMouseX << ", " << currMouseY << endl;
#endif

//	glGetDoublev(GL_DEPTH_RANGE, a);	// Should be the default
//	cout << "GL Depth range: " << a[0] << ", " << a[1] << endl;

	int iy = viewport[3] - currMouseY - 1;
	gluUnProject((GLdouble) currMouseX, (GLdouble) iy, 0.0,
		modelview, projection, viewport, &a[0], &a[1], &a[2]);
	gluUnProject((GLdouble) currMouseX, (GLdouble) iy, 1.0,
		modelview, projection, viewport, &b[0], &b[1], &b[2]);
#ifdef DEBUG
	cout << "viewport[3] = " << viewport[3] << ", iy = " << iy << endl;
	cout << "Intersection line: (" << a[0] << ", " << a[1] << ", " << a[2] << ") - (";
	cout << b[0] << ", " << b[1] << ", " << b[2] << ')' << endl;
#endif

    glPopMatrix();

	Point va(a[0], a[1], a[2]);	// Point on near clipping plane of viewport
	Point vb(b[0], b[1], b[2]);	// Point on far clipping plane of viewport

#ifdef DEBUG
	cout << "va = "; va.print();
	cout << "vb = "; vb.print();
	double w[3];
	gluProject(a[0], a[1], a[2], modelview, projection, viewport, &w[0], &w[1], &w[2]);
	cout << "Reprojection of va = (" << w[0] << ", " << w[1] << ", " << w[2] << ")\n";
	gluProject(b[0], b[1], b[2], modelview, projection, viewport, &w[0], &w[1], &w[2]);
	cout << "Reprojection of vb = (" << w[0] << ", " << w[1] << ", " << w[2] << ")\n";
#endif

	// Determine which slice is displaying
    imagePlanes = control->getImagePlanesPtr();
    if (imagePlanes == NULL) {
		cout << "An image must be loaded to set image landmarks" << endl;
		fl_beep(FL_BEEP_ERROR);
        return;
	}

    Image3D * image = imagePlanes->getImagePtr();
    if (image == NULL) {
		cout << "An image must be loaded to set image landmarks" << endl;
		fl_beep(FL_BEEP_ERROR);
        return;
	}

	Point q;
	double xpos, ypos, zpos;
	double image_landmarks[9];

	count = 0;	// Number of intersections of the ray with displayed slices
	numLandmarkHits = 0;
    if (xSliceOn)
    {
		xpos = imagePlanes->getXCutPlaneModelPos();
		Point p0(xpos, 0.0, 0.0);
		Point p1(xpos, 0.0, 1.0);
		Point p2(xpos, 1.0, 1.0);
#ifdef DEBUG
		cout << "X slice:\n";
		gluProject(xpos, 0.0, 0.0, modelview, projection, viewport, &w[0], &w[1], &w[2]);
		cout << "Projection of P0 = (" << w[0] << ", " << w[1] << ", " << w[2] << ")\n";
		gluProject(xpos, 0.0, 1.0, modelview, projection, viewport, &w[0], &w[1], &w[2]);
		cout << "Projection of P1 = (" << w[0] << ", " << w[1] << ", " << w[2] << ")\n";
		gluProject(xpos, 1.0, 1.0, modelview, projection, viewport, &w[0], &w[1], &w[2]);
		cout << "Projection of P2 = (" << w[0] << ", " << w[1] << ", " << w[2] << ")\n";
#endif
		q = intersect_line_with_plane(va, vb, p0, p1, p2);
//		cout << "q = "; q.print();
		if (q.getY() >= 0.0 && q.getY() <= 1.0 && q.getZ() >= 0.0 && q.getZ() <= 1.0) {
			image_landmarks[3*numLandmarkHits] = q.getX();
			image_landmarks[3*numLandmarkHits + 1] = q.getY();
			image_landmarks[3*numLandmarkHits + 2] = q.getZ();
			if (numLandmarkHits < MAX_NUM_LANDMARKS - 1) {
				numLandmarkHits++;
				count++;
			}
			else
				cout << "Too many landmarks specified" << endl;
		}
		// else clip
	}
    if (ySliceOn)
    {
		ypos = imagePlanes->getYCutPlaneModelPos();
		Point p0(0.0, ypos, 0.0);
		Point p1(0.0, ypos, 1.0);
		Point p2(1.0, ypos, 1.0);
#ifdef DEBUG
		cout << "Y slice:\n";
		gluProject(0.0, ypos, 0.0, modelview, projection, viewport, &w[0], &w[1], &w[2]);
		cout << "Projection of P0 = (" << w[0] << ", " << w[1] << ", " << w[2] << ")\n";
		gluProject(0.0, ypos, 1.0, modelview, projection, viewport, &w[0], &w[1], &w[2]);
		cout << "Projection of P1 = (" << w[0] << ", " << w[1] << ", " << w[2] << ")\n";
		gluProject(1.0, ypos, 1.0, modelview, projection, viewport, &w[0], &w[1], &w[2]);
		cout << "Projection of P2 = (" << w[0] << ", " << w[1] << ", " << w[2] << ")\n";
#endif
		q = intersect_line_with_plane(va, vb, p0, p1, p2);
//		cout << "q = "; q.print();
		if (q.getX() >= 0.0 && q.getX() <= 1.0 && q.getZ() >= 0.0 && q.getZ() <= 1.0) {
			image_landmarks[3*numLandmarkHits] = q.getX();
			image_landmarks[3*numLandmarkHits + 1] = q.getY();
			image_landmarks[3*numLandmarkHits + 2] = q.getZ();
			if (numLandmarkHits < MAX_NUM_LANDMARKS - 1) {
				numLandmarkHits++;
				count++;
			}
			else
				cout << "Too many landmarks specified" << endl;
		}
		// else clip
	}
    if (zSliceOn)
    {
		zpos = imagePlanes->getZCutPlaneModelPos();
		Point p0(0.0, 0.0, zpos);
		Point p1(0.0, 1.0, zpos);
		Point p2(1.0, 1.0, zpos);
#ifdef DEBUG
		cout << "Z slice:\n";
		gluProject(0.0, 0.0, zpos, modelview, projection, viewport, &w[0], &w[1], &w[2]);
		cout << "Projection of P0 = (" << w[0] << ", " << w[1] << ", " << w[2] << ")\n";
		gluProject(0.0, 1.0, zpos, modelview, projection, viewport, &w[0], &w[1], &w[2]);
		cout << "Projection of P1 = (" << w[0] << ", " << w[1] << ", " << w[2] << ")\n";
		gluProject(1.0, 1.0, zpos, modelview, projection, viewport, &w[0], &w[1], &w[2]);
		cout << "Projection of P2 = (" << w[0] << ", " << w[1] << ", " << w[2] << ")\n";
#endif
		q = intersect_line_with_plane(va, vb, p0, p1, p2);
//		cout << "q = "; q.print();
		if (q.getX() >= 0.0 && q.getX() <= 1.0 && q.getY() >= 0.0 && q.getY() <= 1.0) {
			image_landmarks[3*numLandmarkHits] = q.getX();
			image_landmarks[3*numLandmarkHits + 1] = q.getY();
			image_landmarks[3*numLandmarkHits + 2] = q.getZ();
			if (numLandmarkHits < MAX_NUM_LANDMARKS - 1) {
				numLandmarkHits++;
				count++;
			}
			else
				cout << "Too many landmarks specified" << endl;
		}
		// else clip
	}

	if (count == 0) {
		cout << "Selection not close to any image plane" << endl;
		fl_beep(FL_BEEP_ERROR);
		return;
	}

	double x, y, z;
	double zmin;	// Screen coordinates

	// Determine which ray-slice intersection is front-most
	gluProject(image_landmarks[3*(numLandmarkHits - count)],
		image_landmarks[3*(numLandmarkHits - count) + 1],
		image_landmarks[3*(numLandmarkHits - count) + 2],
		modelview, projection, viewport, &x, &y, &z);
#ifdef DEBUG
	cout << '(' << image_landmarks[3*(numLandmarkHits - count)] << ", "
		<< image_landmarks[3*(numLandmarkHits - count) + 1] << ", "
		<< image_landmarks[3*(numLandmarkHits - count) + 2] << ") -> ( "
		<< x << ", " << y << ", " << z << ")\n";
#endif
	zmin = z;
	infront = 0;
	for (j = 1; j < count; j++) {
		gluProject(image_landmarks[3*(numLandmarkHits - count + j)],
			image_landmarks[3*(numLandmarkHits - count + j) + 1],
			image_landmarks[3*(numLandmarkHits - count + j) + 2],
			modelview, projection, viewport, &x, &y, &z);
#ifdef DEBUG
		cout << '(' << image_landmarks[3*(numLandmarkHits - count + j)] << ", "
			<< image_landmarks[3*(numLandmarkHits - count + j) + 1] << ", "
			<< image_landmarks[3*(numLandmarkHits - count + j) + 2] << ") -> ( "
			<< x << ", " << y << ", " << z << ")\n";
#endif
		if (z < zmin) {
			zmin = z;
			infront = j;
		}
	}
#ifdef DEBUG
	cout << "Selected landmark " << infront << " at: "
		<< image_landmarks[3*(numLandmarkHits - count + infront)] << ", "
		<< image_landmarks[3*(numLandmarkHits - count + infront) + 1] << ", "
		<< image_landmarks[3*(numLandmarkHits - count + infront) + 2] << ")\n";
#endif

	int figureId;

	updateImagePlanes(false);
	int l = object->getMarkedLandmark(figureId);
	if (l < 0) {
		// No landmark is marked - this should never happen
		numLandmarkHits--;
		cout << "Warning: not in image mode; no landmark was set\n";
		fl_beep(FL_BEEP_ERROR);
		return;
	}
	else {
		if (control->landmarkIsPaired(figureId, l)) {
			cout << "The current landmark already has an image landmark" << endl;
			fl_beep(FL_BEEP_ERROR);
			return;
		}
		int landmarkId = callback->getCurrentLandmarkIndex();
		// Keep only the front-most point (zmin)
		if (control->setImageLandmark(landmarkId,
			image_landmarks[3*(numLandmarkHits - count + infront)],
			image_landmarks[3*(numLandmarkHits - count + infront) + 1],
			image_landmarks[3*(numLandmarkHits - count + infront) + 2],
			figureId, l) == true)
		{
#ifdef DEBUG
			cout << "Selected image landmark " << landmarkId << " is number " << l
				<< " for figure " << figureId << endl;
#endif
			callback->markImageLandmark(landmarkId, true);	// Set the landmark window Paired column
			callback->nextImageLandmark();	// Update the landmark window caret
		}
		else
			cout << "There is no such landmark!\n";

		// Update displayed counter
		callback->setLandmarksImageNumber(control->getImageLandmarksCount());
	}
#endif
}

// Called from P3DUserInterfaceCallback::startEditLandmarks(),
// when Edit -> Edit Landmarks is selected.
void P3DView::initializeEditLandmarks() {
#ifndef LM_METHOD_OBJ
	int ntiles;
	Bpoint * tile;

#ifdef DEBUG
	cout << "initializeEditLandmarks()" << endl;
#endif
    M3DObject * object = control->getObjectPtr();
    if(object == NULL)
        return;

	// Temporary storage for holding possible landmarks per mouse pick
	if (landmarks == NULL)
		landmarks = new Vector3D[MAX_NUM_LANDMARK_HITS];

	if (landmarkSurface != NULL)
		delete landmarkSurface;
	landmarkSurface = new M3DObjectSurfaceRenderer(object, LANDMARK_SURFACE_LEVEL);
	landmarkSurface->setSurfaceStyle(M3D_SOLID);
	landmarkSurface->render(true);

	// Generate per-figure distance thresholds for fast determination of
	// tiles to be tested for mouse-pick intersection.  The threshold is
	// the maximum length of sides of all rendered quads of a figure.
	if (distThreshold != NULL)
		delete distThreshold;
	distThreshold = new double[object->getFigureCount()];
	for (int nFigure = 0; nFigure < object->getFigureCount(); nFigure++) {

		ThallCode::Pointlist_server2 * pList = landmarkSurface->getPList(nFigure);
		if (pList == NULL) {
			distThreshold[nFigure] = 0.0;
			cout << "Error: missing pList\n";
			continue;
		}
		pList->subdivtileinfo(&ntiles, &tile);
		distThreshold[nFigure] = 0.0;
		for (int vertex = 0; vertex < 4*ntiles; vertex += 4) {
			Vector3D pt0(tile[vertex].pnt[0], tile[vertex].pnt[1], tile[vertex].pnt[2]);
			Vector3D pt1(tile[vertex + 1].pnt[0], tile[vertex + 1].pnt[1], tile[vertex + 1].pnt[2]);
			Vector3D pt2(tile[vertex + 2].pnt[0], tile[vertex + 2].pnt[1], tile[vertex + 2].pnt[2]);
			Vector3D pt3(tile[vertex + 3].pnt[0], tile[vertex + 3].pnt[1], tile[vertex + 3].pnt[2]);
			double d0 = inter_pt_distance(pt0, pt1);
			double d1 = inter_pt_distance(pt1, pt2);
			double d2 = inter_pt_distance(pt2, pt3);
			double d3 = inter_pt_distance(pt3, pt0);
			if (d0 > distThreshold[nFigure])
				distThreshold[nFigure] = d0;
			if (d1 > distThreshold[nFigure])
				distThreshold[nFigure] = d1;
			if (d2 > distThreshold[nFigure])
				distThreshold[nFigure] = d2;
			if (d3 > distThreshold[nFigure])
				distThreshold[nFigure] = d3;
		}
	}

	// Set up the image landmarks vector, if necessary.
	int num_landmarks = object->getLandmarkCount();
	if (num_landmarks > 0) {
		if (control->landmarksAllocation() > num_landmarks) {
			// A figure was lost by deletion or merging.
			// No effort is made to save the image landmarks in this case.
			control->allocateImageLandmarks(num_landmarks);
		}
		else if (control->landmarksAllocation() < num_landmarks) {
			// A new figure was added
			for (int j = control->landmarksAllocation(); j < num_landmarks; j++)
				control->appendImageLandmark();
		}
	}
#endif
}

void P3DView::clearBoundaries()
{
	xBoundary.clear();  // See drawBoundaryCurves
	// Don't need to clear yBoundary and zBoundary
	boundaryCurvesDone = false;
	// This clears the display list without deleting it
	glNewList(gBoundaryCurveDisplayList, GL_COMPILE);
	glEndList();
}

void P3DView::setBoundariesLevel()
{
	if (control == NULL)
		return;

    M3DObject * object = control->getObjectPtr();
    if(object == NULL)
        return;

	for (int figureNum = 0; figureNum < object->getFigureCount(); figureNum++) {
		xBoundary[figureNum].setLevel(surfaceLevel);
		yBoundary[figureNum].setLevel(surfaceLevel);
		zBoundary[figureNum].setLevel(surfaceLevel);
	}
}

// Set up generation of the 3 cut plane views for the Model Window
void P3DView::initTextures()
{
    ImagePlanes * imagePlanes;

    if(control == NULL)
        return;

    imagePlanes = control->getImagePlanesPtr();

    if(imagePlanes == NULL)
        return;

    int size = imagePlanes->getImagePlaneSize();
    if(size <= 0)
        return;

	// The calls below don't draw the textures, they just set up the texture
	// objects for later drawing.  See draw().
    glBindTexture(GL_TEXTURE_2D, textureNames[0]);

    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);

    glBindTexture(GL_TEXTURE_2D, textureNames[1]);

    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);

    glBindTexture(GL_TEXTURE_2D, textureNames[2]);

    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);

/*    glBindTexture(GL_TEXTURE_2D, textureNames[3]);

    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
*/
	imagePlanesDone = false;

    if (control->imagePlanesStatus())
		control->toggleImagePlanesStatus();
}

void P3DView::rockingMotion() {
	if (direction) {
		currentAngle -= incrementalAngle;
		if (currentAngle < -deflectionAngle) {
			currentAngle += incrementalAngle + incrementalAngle;
			direction = false;
		}
	}
	else {
		currentAngle += incrementalAngle;
		if (currentAngle > deflectionAngle) {
			currentAngle -= (incrementalAngle + incrementalAngle);
			direction = true;
		}
	}
	glRotated(currentAngle, 0.0, 1.0, 0.0);
}


// Drawing routines


void P3DView::checkToStopBlending()
{
	Fl::check();
	if (blendSurface == false) {
		blendedSurfaceRenderer.interrupt();    // The blend button was released
		interruptOccurred = true;
	}
	else if (surfaceBlendedLevelChanged) {
		blendedSurfaceRenderer.interrupt();    // The blending level changed
		surfaceBlendedLevelChanged = false;
		surfaceDone = false;
		interruptOccurred = true;
	}
}

void P3DView::drawAuxSurfaces()
{
	int numVerts, numEdgePoints, numFaces;
	int i, n, numTiles;
	const float * vertices;
	const int * edgePoints;
	const Vector3D * tileList;
	Bpoint * pntlist;
	ThallCode::PindexList * faces;
	Vector3D * normalsPerFace;
	Vector3D * normals;
	unsigned char * cnt;
	Vector3D u, v, w;

	BYU * auxTileSet = (BYU *) control->getAuxTileSetPtr();
	if (auxTileSet == NULL)
		return;

	numVerts = auxTileSet->vertexCount();
	vertices = auxTileSet->vertexList();
	numEdgePoints = auxTileSet->edgeCount();
	edgePoints = auxTileSet->edgeList();
#ifdef DEBUG
	cout << "numVerts = " << numVerts << "  numEdgePoints = " << numEdgePoints
		<< '\n';
#endif
	// Specify faces
	n = 0;
	int  edges_per_tile = auxTileSet->edgesPerTile();
	numFaces = numEdgePoints/edges_per_tile;
	faces = new ThallCode::PindexList[numFaces];
	for (i = 0; i < numEdgePoints; i += edges_per_tile) {
		faces[n].Create(edges_per_tile);
		faces[n][0] = edgePoints[i];
		faces[n][1] = edgePoints[i + 1];
		faces[n][2] = edgePoints[i + 2];
		if (edges_per_tile == 4)
			faces[n][3] = edgePoints[i + 3];
		n++;
	}

	tileList = auxTileSet->tileList();
	numTiles = auxTileSet->tileListLength();
#ifdef DEBUG
	cout << "numTiles = " << numTiles << '\n';
#endif

	// Calculate normals at the vertex coordinates.
	// The input data is assumed to use the same ordering of points
	// for each tile, with respect to the direction of the tile normal.
	// Quads are assumed to be flat, so only 3 of the 4 points are
	// used to compute the normal.
	normalsPerFace = new Vector3D[numFaces];
	n = 0;
	for (i = 0; i < numTiles; i++) {
		const Vector3D & p0 = tileList[n];
		const Vector3D & p1 = tileList[n + 1];
		const Vector3D & p2 = tileList[n + 2];
		u = p2 - p0;
		v = p1 - p0;
		w = v.cross(u);
		w.normalize();
		normalsPerFace[i] = w;
		n += edges_per_tile;
	}

	normals = new Vector3D[numVerts];
	cnt = new unsigned char[numVerts];
	for (i = 0; i < numVerts; i++) {
		normals[i].set(0.0, 0.0, 0.0);
		cnt[i] = 0;
	}
	n = 0;
	if (edges_per_tile == 3) {	// Tris
		for (i = 0; i < numTiles; i++) {
			int i0 = edgePoints[n++];
			int i1 = edgePoints[n++];
			int i2 = edgePoints[n++];
			w = normalsPerFace[i];
			normals[i0] += w;
			normals[i1] += w;
			normals[i2] += w;
			cnt[i0]++;
			cnt[i1]++;
			cnt[i2]++;
		}
	}
	else {	// Quads
		for (i = 0; i < numTiles; i++) {
			int i0 = edgePoints[n++];
			int i1 = edgePoints[n++];
			int i2 = edgePoints[n++];
			int i3 = edgePoints[n++];
			w = normalsPerFace[i];
			normals[i0] += w;
			normals[i1] += w;
			normals[i2] += w;
			normals[i3] += w;
			cnt[i0]++;
			cnt[i1]++;
			cnt[i2]++;
			cnt[i3]++;
		}
	}
	for (i = 0; i < numVerts; i++) {
		if (cnt[i] < 3) {
			cout << "Invalid topology detected: vertex " << i << " is used "
				<< (int) cnt[i] << " times\n";
			delete [] normals;
			delete [] normalsPerFace;
			delete [] cnt;
			delete [] faces;
			return;
		}
		normals[i] /= cnt[i];
		double l = normals[i].norm();
		if (l == 0.0) {
			cout << "Bad normal\n";
			normals[i].set(1.0, 0.0, 0.0);
		}
		else
			normals[i] /= l;
	}
	delete [] normalsPerFace;
	delete [] cnt;

#ifdef DEBUG
	cout << "n = " << n << '\n';
#endif

	// Specify vertices
	pntlist = new Bpoint[numVerts];
	n = 0;
	for (i = 0; i < numVerts; i++) {
		pntlist[i].setvals(vertices[n], vertices[n + 1], vertices[n + 2],
			normals[i].getX(), normals[i].getY(), normals[i].getZ(),
			0.0, 0.0, 0.0, 0.0);
		n += 3;
	}
	delete [] normals;

	// Draw the auxiliary tile surface
	M3DObjectSurfaceRenderer surfRenderer;
	surfRenderer.setSurfaceLevel(surfaceLevel);
	surfRenderer.setSurfaceStyle(surfaceStyle);
	surfRenderer.setLineWidth(surfaceLineWidth);
	surfRenderer.renderFigureFromControlPoints(numVerts, numFaces, pntlist, faces);

	delete [] pntlist;
	delete [] faces;
}

#ifdef BINARY

void P3DView::drawDistanceVectors()
{
	if (distanceVectorList == NULL)
		return;
	
	if (surfaceStyle == M3D_NONE || surfaceStyle == M3D_POINTCLOUD)
		return;

	DistanceVectorList::DistVectorPList & plist = distanceVectorList->plist;
	if (plist.size() == 0)
		return;
	else {
		glPushAttrib(GL_LINE_BIT);
		glLineWidth(1.0);
		glBegin(GL_LINES);
		for (int i = 0; i != plist.size(); i++) {
			DistanceVector & p = plist[i];
			if (p.method == 1)
				glColor3d(0.8, 0.5, 0.0);	// Orange
			else
				glColor3d(0.0, 1.0, 0.0);	// Green
			glVertex3d(p.p.getX(), p.p.getY(), p.p.getZ());
			glVertex3d(p.p.getX() + p.grad.getX(), p.p.getY() + p.grad.getY(), p.p.getZ() + p.grad.getZ());
			if (p.method == 2)
				glColor3d(0.8, 0.0, 0.5);	// Brick red
			else
				glColor3d(0.5, 0.5, 1.0);	// Medium blue
			glVertex3d(p.p.getX(), p.p.getY(), p.p.getZ());
			glVertex3d(p.p.getX() + p.n.getX(), p.p.getY() + p.n.getY(), p.p.getZ() + p.n.getZ());
		}
		glEnd();
		glPopAttrib();
	}
}

#endif

void P3DView::drawSurfaces()
{
    M3DObject * object;

#if defined(DEBUG) || defined(DRAWING_DEBUG)
	cout << "drawSurfaces() called;";
#ifdef DRAWING_DEBUG
	cout << " surfaceStyle = " << surfaceStyle;
#endif
	cout << endl;
#endif
    object = control->getObjectPtr();

	// Draw the surface
	if (control->getPartialDrawingMode() == true)
	{
#ifdef DRAWING_DEBUG
		cout << "PartialDrawingMode = true" << endl;
#endif
		double slider_dist = callback->getInterObjectDistance();

		M3DObjectSurfaceRenderer surfRenderer(object, surfaceLevel);
		surfRenderer.setSurfaceStyle(surfaceStyle);
		surfRenderer.setLineWidth(surfaceLineWidth);
		surfRenderer.setPartialSurfaceStyle(partialSurfaceStyle);
		surfRenderer.partial_render(callback->getPartialFigures(),
			slider_dist, control->getPartialSurfaceMode(), partialSurfaceLevel);
	}
	else {
		// Except when constraints are being specified (above), surfaceStyle
		// determines whether or not the surface is to be rendered.
		if (surfaceStyle != M3D_NONE) {
			if (! blendSurface ||
				object->getFigureCount() == object->getFigureTreeCount() ) {
				surfaceRenderer.setObject(object);
				surfaceRenderer.setSurfaceLevel(surfaceLevel);
				surfaceRenderer.setLineWidth(surfaceLineWidth);
#ifndef PRODUCTION_VERSION
				if (colorSurface) {
					float min, max;
					surfaceRenderer.setSurfaceStyle(colorSurfaceStyle);
					surfaceRenderer.outlineTiles(outlineSurfaceTiles);
					getSurfaceIntensityRange(min, max);
					surfaceRenderer.surfaceColorRange(min, max);
				}
				else
#endif
					surfaceRenderer.setSurfaceStyle(surfaceStyle);
				surfaceRenderer.render();
#ifndef LM_METHOD_OBJ
				// Maximum-distance landmark
				if (control->getMatchPtr() != NULL)
					surfaceRenderer.drawDots(control->getMatchPtr()->getMaxDistance(), 1);
#endif
				// Landmark drawing requires an M3DObjectSurfaceRenderer object, so
				// cannot be done when the rendering style is M3D_NONE.  This could
				// be changed by always having an M3DObjectSurfaceRenderer object.
				if (landmarksOn) {
					surfaceRenderer.drawDots(LANDMARK_SIZE, defaultLandmarkWidth,
						landmarkMode, markedLandmarkWidth,
						(surfaceStyle == M3D_SOLID ? NULL : defaultLandmarkColor));
				}
			}
			else {
				if (! blendedRendererRunning) {
					blendedSurfaceRenderer.setSurfaceStyle(surfaceStyle);
					blendedSurfaceRenderer.setLineWidth(surfaceLineWidth);
#ifndef PRODUCTION_VERSION
					blendedSurfaceRenderer.setCurrentScale(trackball.getCurrScale());
					blendedSurfaceRenderer.showBlendingRegion(displayBlendRegion);
#endif
					if (object->isModified() || blendingLevel != surfaceLevel) {
						blendedSurfaceRenderer.init(object, surfaceLevel);
						blendedRendererRunning = true;
						if (blendedSurfaceRenderer.render())
							blendingLevel = surfaceLevel;
						blendedRendererRunning = false;
					}
					else {
						blendedRendererRunning = true;
						(void) blendedSurfaceRenderer.render();
						blendedRendererRunning = false;
					}
					if (landmarksOn) {
						M3DObjectSurfaceRenderer surfRenderer(object, surfaceLevel);
						surfRenderer.drawDots(LANDMARK_SIZE, defaultLandmarkWidth,
							landmarkMode, markedLandmarkWidth,
							(surfaceStyle == M3D_SOLID ? NULL : defaultLandmarkColor));
					}
				}
			}
		}
	}

#ifdef LM_METHOD_OBJ
	// Now draw the image landmark (new method) object figures' surfaces
	if (surfaceStyle != M3D_NONE) {
	    object = control->landmarkObject;
		M3DObjectSurfaceRenderer surfRenderer(object, surfaceLevel);
		surfRenderer.setSurfaceStyle(surfaceStyle);
		surfRenderer.setLineWidth(surfaceLineWidth);
		surfRenderer.render();
	}
#endif
}

typedef void (*voidFn)(void*);

void signalInterruptRedraw(void * arg) {		// Non-class function used below
	extern P3DUserInterface * ui;

	Fl::remove_check((voidFn) signalInterruptRedraw, arg);
	ui->modelView->update();
}

void signalDamageRedraw(void * arg) {		// Non-class function used below
	extern P3DUserInterface * ui;
	redraw_t damaged;

    damaged = (redraw_t) (long) arg;
	Fl::remove_check((voidFn) signalDamageRedraw, arg);

	switch (damaged) {
		case Redraw:
						ui->modelView->redraw();
						break;
		case UpdateObject:
						ui->modelView->updateObject();
						break;
		case IgnoreRedraw:
		default:
						break;
	}
}


// selectively enable/disable groups of clip planes
void P3DView::enableClips(int mask,		// select planes
				 int enable)	// what to do to selected planes
{
	if (!control) { return; }

	struct clipPlane *clips = control->getClipPlanes();
	for (int c=0; c < 6; c++) {
		int plane = GL_CLIP_PLANE0 + c;
		if (clips[c].mask & mask) {

			// masks match, so enable or disable the clip plane
			if (enable && (clips[c].mask & CLIP_MASK_ENABLED)) {
				glClipPlane(plane, clips[c].equation);
				glEnable(plane);
#ifdef DEBUG
				cout << "clip plane " << c << ":"
					<< " mask=" << clips[c].mask
					<< " A=" << clips[c].equation[0]
					<< " B=" << clips[c].equation[1]
					<< " C=" << clips[c].equation[2]
					<< " D=" << clips[c].equation[3]
					<< endl;
#endif
			}
			else
				glDisable(plane);

		}
	}
}

void P3DView::draw()
{
    M3DObject * object;

#ifdef DEBUG
	cout << "P3DView::draw()\n";
#endif
	if (isDrawing) {
		damaged = Redraw;
		return;
	}
    isDrawing = true;

    glShadeModel(tileType);

    // Just do initialization code once
    if(! init)
    {
        glMatrixMode(GL_PROJECTION);
        glLoadIdentity();

        glMatrixMode(GL_MODELVIEW);
        glLoadIdentity();

        float position0[] = {1.0f, 1.0f, 1.5f, 0.0f};	// On the right
        float position1[] = {-1.0f, 1.0f, 1.5f, 0.0f};	// On the left
        int shininess = 100;

        float mat[] = {1.0f, 1.0f, 1.0f, 1.0f};
        float mat_spec[] = {0.6f, 1.0f, 0.6f, 1.0f};


		glEnable(GL_LIGHTING);
		if (twoLights) {
			float ambient[] = {0.1f, 0.1f, 0.1f, 1.0f};
			float diffuse[] = {0.3f, 0.3f, 0.3f, 1.0f};
			float specular[] = {0.2f, 0.2f, 0.2f, 1.0f};

			glLightfv(GL_LIGHT0, GL_AMBIENT, ambient);
			glLightfv(GL_LIGHT0, GL_DIFFUSE, diffuse);
			glLightfv(GL_LIGHT0, GL_SPECULAR, specular);
			glLightfv(GL_LIGHT0, GL_POSITION, position0);
			glEnable(GL_LIGHT0);

			glLightfv(GL_LIGHT1, GL_AMBIENT, ambient);
			glLightfv(GL_LIGHT1, GL_DIFFUSE, diffuse);
			glLightfv(GL_LIGHT1, GL_SPECULAR, specular);
			glLightfv(GL_LIGHT1, GL_POSITION, position1);
			glEnable(GL_LIGHT1);
		}
		else {
			float ambient[] = {0.2f, 0.2f, 0.2f, 1.0f};
			float diffuse[] = {0.7f, 0.7f, 0.7f, 1.0f};
			float specular[] = {0.3f, 0.3f, 0.3f, 1.0f};

			glLightfv(GL_LIGHT0, GL_AMBIENT, ambient);
			glLightfv(GL_LIGHT0, GL_DIFFUSE, diffuse);
			glLightfv(GL_LIGHT0, GL_SPECULAR, specular);
			glLightfv(GL_LIGHT0, GL_POSITION, position0);
			glEnable(GL_LIGHT0);

			glDisable(GL_LIGHT1);
		}

		glMaterialfv(GL_FRONT, GL_AMBIENT, mat);			// AGG: This appears to have been new (not in BP).  What's the effect?
        glMaterialfv(GL_FRONT, GL_DIFFUSE, mat);
        glMaterialfv(GL_FRONT, GL_SPECULAR, mat_spec);
        glMateriali(GL_FRONT, GL_SHININESS, shininess);
        glDisable(GL_LIGHTING);

        glEnable(GL_NORMALIZE);
        glEnable(GL_DEPTH_TEST);
        glDepthFunc(GL_LEQUAL);

        glEnable(GL_COLOR_MATERIAL);
        glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);

        if(displayListsInit == false)
        {
            // Print OpenGL information
#ifdef DRAWING_DEBUG
            printf("OpenGL Version %s\n", glGetString(GL_VERSION));
            printf("OpenGL Vendor: %s\n", glGetString(GL_VENDOR));
            printf("OpenGL Renderer: %s\n\n", glGetString(GL_RENDERER));
#endif

			gSurfaceDisplayList = glGenLists(1);
#ifdef DRAWING_DEBUG
            printf("Surface Display List Value: %d\n", gSurfaceDisplayList);
#endif

			gAuxSurfaceDisplayList = glGenLists(1);
#ifdef DRAWING_DEBUG
            printf("Auxiliary Surface Display List Value: %d\n", gAuxSurfaceDisplayList);
#endif

	        gModelDisplayList = glGenLists(1);
#ifdef DRAWING_DEBUG
            printf("Model Display List Value: %d\n", gModelDisplayList);
#endif

            gImagePlanesDisplayList = glGenLists(1);
#ifdef DRAWING_DEBUG
            cout << "Image Planes Display List Value: " << gImagePlanesDisplayList << endl;
#endif

            gBoundaryCurveDisplayList = glGenLists(1);
#ifdef DRAWING_DEBUG
            cout << "Boundary Curve Display List Value: " << gBoundaryCurveDisplayList << endl;
#endif

#ifdef BINARY
			gDistanceVectorDisplayList	= glGenLists(1);
#ifdef DRAWING_DEBUG
            cout << "Distance Vector Display List Value: " << gDistanceVectorDisplayList << endl;
#endif
#endif

//          glGenTextures(4, (GLuint *)&(textureNames[0]));
            glGenTextures(3, (GLuint *)&(textureNames[0]));
#ifdef DRAWING_DEBUG
            cout << "Generating textures: val = " << textureNames[0] << endl;
#endif

            gMasterDisplayList = glGenLists(1);
#ifdef DRAWING_DEBUG
            printf("Master Display List Value: %d\n", gMasterDisplayList);
#endif

            displayListsInit = true;
        }

		glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
		glPixelStorei(GL_UNPACK_ROW_LENGTH, 0);
		glPixelStorei(GL_UNPACK_SKIP_ROWS, 0);
		glPixelStorei(GL_UNPACK_SKIP_PIXELS, 0);

		// This prevents an OpenGL Error from occurring in later
		// calls to this function.  It appears to be a bug in FLTK,
		// or OpenGL which only happens the very first time function
		// drawString() is used.  The effect of the error was that a
		// rectangle was always drawn before the string.  To see it,
		// comment this out, compile, and display a model having a
		// a figure name or a spoke-end landmark.
		M3DObjectRenderer::drawString(" ");

        init = true;
    }

	enableClips(CLIP_MASK_CAMERA, 1);	// enable all-objects clip
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glViewport(0, 0, w(), h());
    glOrtho(-VIEW_WIDTH / 2.0, VIEW_WIDTH / 2.0, -VIEW_HEIGHT / 2.0,
            VIEW_HEIGHT / 2.0, -FRUSTRUM_DEPTH, FRUSTRUM_DEPTH);
	//enableClips(CLIP_MASK_ALL, 0);	// disable all planes
	//enableClips(CLIP_MASK_CAMERA, 1);	// enable all-objects clip

    glClearColor(backgroundColor[0], backgroundColor[1], backgroundColor[2], 1.0);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    if (control == NULL) {
        isDrawing = false;
        return;
    }

	// Prepare to generate image plane views for the Model Window.  The
	// actual drawing is done by the later call to drawImageSlices().
    if (control->imagePlanesStatus())
        initTextures();

    object = control->getObjectPtr();

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

	// copy tile colors from object (if object present and
	//  per-figure colors needed) then compile the tiles' GL list(s) if needed.
    TileSet * tileSet = control->getTileSetPtr();
    if (tileSet) {
		if (tileRenderer == NULL)
			tileRenderer = new TileSetRenderer;

		if (defaultTilesColorChanged) {
			tileRenderer->setDefaultColor(defaultTilesColor);
			defaultTilesColorChanged = false;
		}

		if (object != NULL) {
			int nfigs;
			bool redrawTiles;

			redrawTiles = ! tileSetDone;
			nfigs = object->getFigureCount();
			if (currentFigureCount != nfigs) {
				currentFigureCount = nfigs;
				delete [] tileColors;
				tileColors = new const float *[nfigs];
				tileColorsChanged = true;
				redrawTiles = true;
			}
			if (tileColorsChanged) {
				for (int i = 0; i < currentFigureCount; i++) {	
					M3DFigure * fig = object->getFigurePtr(i);
					tileColors[i] = fig->getColor();
				}
				tileColorsChanged = false;
				redrawTiles = true;
			}
			if (redrawTiles) {
				tileRenderer->compileDispList(tileSet, object->getFigureCount(),
					tileColors);
			}
		}
		else {
			if (! tileSetDone)
			{
				tileRenderer->compileDispList(tileSet);
			}
		}
		tileSetDone = true;
	}

    // If the view scale is changing, then we have to redraw the object
    // (point-size is scale relative)
    if (trackball.getCurrScale() != gLastScale) {
        objectDone = false;
        gMasterDispListIsValid = false;
        gLastScale = trackball.getCurrScale();
    }

    if (rockingMode())
	    rockingMotion();

    trackball.doViewTransform();

/*  if(viewType == YZ_ORTHO)
        trackball.rotate(90.0, 0.0, 1.0, 0.0);
    else if(viewType == XZ_ORTHO)
    {
        trackball.rotate(90.0, 1.0, 0.0, 0.0);
        trackball.translate(0.0, 0.0, -1.0);
    }*/

    M3DPrimitive * interpolatedPrimitive = control->getInterpolatedPrimitive();
    if (interpolatedPrimitive != NULL) {

		M3DPrimitiveRenderer primitiveRenderer;

        primitiveRenderer.setPrimitive(interpolatedPrimitive);
        primitiveRenderer.setPointScaleFactor(1.0 / gLastScale);
        primitiveRenderer.render();
    }

#ifdef DRAWING_DEBUG
	cout << "gMasterDispListIsValid=" << gMasterDispListIsValid;
	cout << " objectDone=" << objectDone;
	if (object)
		cout << " object->isModified()=" << object->isModified();
	cout << " surfaceDone=" << surfaceDone << endl;
#endif

    if (gMasterDispListIsValid == true) {
        glCallList(gMasterDisplayList);
	}
    else
    {
#ifdef DRAWING_DEBUG
		printMatrices();
#endif
        if (object != NULL)
        {
			updateAllCutPlanes();

            if (object->isModified() || (! objectDone)) {
				glNewList(gModelDisplayList, GL_COMPILE);
				drawObject(constraintArrows);
				// Not calling drawFigureNames() and drawObjectLandmarks() here,
				// although they could be called with some changes.
				glEndList();
				objectDone = true;
			}

#ifndef PRODUCTION_VERSION
			if (colorSurface && ! imagePlanesDone)
				surfaceDone = false;
#endif
            if (object->isModified() || (! surfaceDone)) {
				glNewList(gSurfaceDisplayList, GL_COMPILE);
				drawSurfaces();
				glEndList();
				if (! interruptOccurred)
					surfaceDone = true;
			}

#ifdef BINARY
			if (object->isModified() || (! distanceVectorsDone)) {
				glNewList(gDistanceVectorDisplayList, GL_COMPILE);
				drawDistanceVectors();
				glEndList();
				distanceVectorsDone = true;
			}
#endif

            if (object->isModified() || surfaceLevelChanged || ! boundaryCurvesDone) {
                glNewList(gBoundaryCurveDisplayList, GL_COMPILE);
                if (drawBoundary && (xSliceOn || ySliceOn || zSliceOn))
                    drawBoundaryCurves(object->isModified() || surfaceLevelChanged);
                glEndList();
                boundaryCurvesDone = true;
            }

			if (object->isModified())
				object->setModified(false);
        }
#ifdef BINARY
		else {
			glNewList(gDistanceVectorDisplayList, GL_COMPILE);
			glEndList();
        }
#endif

		if (! auxTileSetDone) {
			glNewList(gAuxSurfaceDisplayList, GL_COMPILE);
			drawAuxSurfaces();	// Does nothing if no BYU file is loaded
			glEndList();
			auxTileSetDone = true;
		}

        // If image slices are on, scale image to unit cube and draw it
        if (! imagePlanesDone) {
            glNewList(gImagePlanesDisplayList, GL_COMPILE);

            glPushMatrix();
            glScaled(STD_AXES_LENGTH, STD_AXES_LENGTH, STD_AXES_LENGTH);

//            glColor4f(1.0f, 1.0f, 1.0f, imageAlpha);
            glColor3f(1.0f, 1.0f, 1.0f);
            if(xSliceOn || ySliceOn || zSliceOn)
                drawImageSlices();	// Generate image plane views for the Model Window

            glPopMatrix();

            glEndList();

            imagePlanesDone = true;
        }

		// The call lists created above are nested inside the master list
        glNewList(gMasterDisplayList, GL_COMPILE_AND_EXECUTE);

		glEnable(GL_BLEND);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
		glEnable(GL_POINT_SMOOTH);
		glEnable(GL_LINE_SMOOTH);
		//glEnable(GL_POLYGON_SMOOTH);
		glHint(GL_POINT_SMOOTH_HINT, GL_NICEST);
		glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
		//glHint(GL_POLYGON_SMOOTH_HINT, GL_NICEST);

		enableClips(CLIP_MASK_AXES,1);
        if (unitCubeOn)
        {
            glColor3f(1.0, 1.0, 1.0);
            glBegin(GL_LINES);

			if (! stdAxesOn)
			{
				glVertex3d(0.0, 0.0, 0.0);
				glVertex3d(STD_AXES_LENGTH, 0.0, 0.0);
				
				glVertex3d(0.0, 0.0, 0.0);
				glVertex3d(0.0, STD_AXES_LENGTH, 0.0);
				
				glVertex3d(0.0, 0.0, 0.0);
				glVertex3d(0.0, 0.0, STD_AXES_LENGTH);
			}

            glVertex3d(0.0, STD_AXES_LENGTH, 0.0);
            glVertex3d(STD_AXES_LENGTH, STD_AXES_LENGTH, 0.0);

            glVertex3d(0.0, STD_AXES_LENGTH, 0.0);
            glVertex3d(0.0, STD_AXES_LENGTH, STD_AXES_LENGTH);

            glVertex3d(STD_AXES_LENGTH, 0.0, 0.0);
            glVertex3d(STD_AXES_LENGTH, STD_AXES_LENGTH, 0.0);

            glVertex3d(STD_AXES_LENGTH, 0.0, 0.0);
            glVertex3d(STD_AXES_LENGTH, 0.0, STD_AXES_LENGTH);

            glVertex3d(0.0, STD_AXES_LENGTH, STD_AXES_LENGTH);
            glVertex3d(STD_AXES_LENGTH, STD_AXES_LENGTH, STD_AXES_LENGTH);

            glVertex3d(0.0, STD_AXES_LENGTH, STD_AXES_LENGTH);
            glVertex3d(0.0, 0.0, STD_AXES_LENGTH);

            glVertex3d(STD_AXES_LENGTH, STD_AXES_LENGTH, STD_AXES_LENGTH);
            glVertex3d(STD_AXES_LENGTH, STD_AXES_LENGTH, 0.0);

            glVertex3d(STD_AXES_LENGTH, STD_AXES_LENGTH, STD_AXES_LENGTH);
            glVertex3d(STD_AXES_LENGTH, 0.0, STD_AXES_LENGTH);

            glVertex3d(STD_AXES_LENGTH, 0.0, STD_AXES_LENGTH);
            glVertex3d(0.0, 0.0, STD_AXES_LENGTH);

            glEnd();

        }
        if (stdAxesOn)
        {
			Image3D *image = control->getImagePtr();
			float tickStep, tickMax = STD_AXES_LENGTH;
			if (image) {
				// length of 1cm in model coords (actually, in standard axes len coords)
				tickStep = STD_AXES_LENGTH/image->maxExtent();
			}

			glColor3f(1.0, 0.0, 0.0);
			if (image) {
				tickMax = STD_AXES_LENGTH*image->getXExtent()/image->maxExtent();
				
				// draw a tick mark every 5cm along X axis until image data ends
				glLineWidth(1.0);
				for (float tick=5*tickStep; tick<tickMax; tick += 5*tickStep)
				{
					glBegin(GL_LINES);
					glVertex3d(tick, tickStep, 0.0);
					glVertex3d(tick, 0.0, 0.0);
					glVertex3d(tick, 0.0, 0.0);
					glVertex3d(tick, 0.0, tickStep);
					glEnd();
				}
			}
			
			// draw whole/remaining axis in white to indicate no image data is there
			glLineWidth(3.0);
			glBegin(GL_LINES);
			glVertex3d(0.0, 0.0, 0.0);
			glVertex3d(tickMax, 0.0, 0.0);
			glColor3f(1.0, 1.0, 1.0);
			glVertex3d(tickMax, 0.0, 0.0);
			glVertex3d(STD_AXES_LENGTH, 0.0, 0.0);
			glEnd();
			
			
			glColor3f(1.0, 0.0, 1.0);
			if (image) {
				tickMax = STD_AXES_LENGTH*image->getYExtent()/image->maxExtent();
				
				// draw a tick mark every 5cm along Y axis until image data ends
				glLineWidth(1.0);
				for (float tick=5*tickStep; tick<tickMax; tick += 5*tickStep)
				{
					glBegin(GL_LINES);
					glVertex3d(tickStep, tick, 0.0);
					glVertex3d(0.0 , tick, 0.0);
					glVertex3d(0.0 , tick, 0.0);
					glVertex3d(0.0, tick, tickStep);
					glEnd();
				}
			}
			
			// draw whole/remaining axis in white to indicate no image data is there
			glLineWidth(3.0);
			glBegin(GL_LINES);
			glVertex3d(0.0, 0.0, 0.0);
			glVertex3d(0.0, tickMax, 0.0);
			glColor3f(1.0, 1.0, 1.0);
			glVertex3d(0.0, tickMax, 0.0);
			glVertex3d(0.0, STD_AXES_LENGTH, 0.0);
			glEnd();
			
			
			glColor3f(0.0, 0.0, 1.0);
			if (image) {
				tickMax = STD_AXES_LENGTH*image->getZExtent()/image->maxExtent();
				
				// draw a tick mark every 5cm along Z axis until image data ends
				glLineWidth(1.0);
				for (float tick=5*tickStep; tick<tickMax; tick += 5*tickStep)
				{
					glBegin(GL_LINES);
					glVertex3d(tickStep, 0.0       , tick);
					glVertex3d(0.0       , 0.0       , tick);
					glVertex3d(0.0       , 0.0       , tick);
					glVertex3d(0.0       , tickStep, tick);
					glEnd();
				}
			}
			
			// draw whole/remaining axis in white to indicate no image data is there
			glLineWidth(3.0);
			glBegin(GL_LINES);
			glVertex3d(0.0, 0.0, 0.0);
			glVertex3d(0.0, 0.0, tickMax);
			glColor3f(1.0, 1.0, 1.0);
			glVertex3d(0.0, 0.0, tickMax);
			glVertex3d(0.0, 0.0, STD_AXES_LENGTH);
			glEnd();

			glLineWidth(1.0);

        }
		enableClips(CLIP_MASK_AXES,0);

		enableClips(CLIP_MASK_MODEL,1);
        if (surfaceDone)
			glCallList(gSurfaceDisplayList);	// Draw the model surface
		enableClips(CLIP_MASK_MODEL,0);

		if (auxTileSetOn)
			glCallList(gAuxSurfaceDisplayList);	// Draw the auxiliary surface

#ifdef BINARY
		glCallList(gDistanceVectorDisplayList);	// Now draw the distance vectors
#endif

		enableClips(CLIP_MASK_IMAGE,1);
        if (imagePlanesDone)
			glCallList(gImagePlanesDisplayList);
		enableClips(CLIP_MASK_IMAGE,0);

        if (objectDone) {
			glCallList(gModelDisplayList);	// Draw the model
#ifdef LM_METHOD_OBJ
			drawObjectLandmarks();
#endif
			if (labelFigures)
				drawFigureNames();
		}

        if (boundaryCurvesDone)
			glCallList(gBoundaryCurveDisplayList);

		if (landmarksOn)
			drawImageLandmarks();

		drawAllCutPlanes();

        if (tileSetOn && tileRenderer)
		{
			enableClips(CLIP_MASK_TILE,1);
            tileRenderer->render(tileAlpha, useDfltTileSetColor, renderTileSetSolid);
			enableClips(CLIP_MASK_TILE,0);
		}

        glEndList();	// gMasterDisplayList
		if (! interruptOccurred)
			gMasterDispListIsValid = true;
		else
			gMasterDispListIsValid = false;

        glDisable(GL_BLEND);
    }

    if (mouseDown)
    {
        glPushMatrix();
        glLoadIdentity();
        glDisable(GL_DEPTH_TEST);

        drawMouseBox(lastMouseX, lastMouseY, currMouseX, currMouseY);

        glEnable(GL_DEPTH_TEST);

        glPopMatrix();
    }

    isDrawing = false;

    if (dumpScreen)
        screenDump();

	if (interruptOccurred)
	{
		Fl::add_check((voidFn) signalInterruptRedraw, NULL);
		interruptOccurred = false;
	}

	if (must_redraw) {
		// A cut plane window was updated, so must redraw it in the model window
		must_redraw = false;
		damaged = UpdateObject;
	}

	if (damaged > IgnoreRedraw)
	{
		// Note, damaged is passed as an address, and so, by value
		Fl::add_check((voidFn) signalDamageRedraw, (void *) damaged);
		damaged = IgnoreRedraw;
	}
}

void P3DView::drawObject(bool showConstraintArrows)
{
    M3DObject * object;
    M3DObjectRenderer objectRenderer;
    M3DPrimitiveRenderer primRenderer;

    if (control == NULL)
        return;

    object = control->getObjectPtr();
    if (object) {
		objectRenderer.setObjectPtr(object);
		objectRenderer.setCurrentScale(trackball.getCurrScale());

		objectRenderer.setAtomVectorVisibility(primitiveVectorsOn);
		objectRenderer.setMeshVisibility(medialSheetOn);
		objectRenderer.setAtomVisibility(atomsOn);

		primRenderer.setPointColor(1.0f, 0.0f, 1.0f);
		primRenderer.setSelectedPointColor(1.0f, 0.0f, 1.0f);
		objectRenderer.setMarkedPrimitiveId(control->getMarkedPrimitiveId());
		objectRenderer.setMarkedPrimitiveRenderer(&primRenderer);
		objectRenderer.setConstraintArrows(showConstraintArrows);

		objectRenderer.render();
	}
}

void P3DView::drawImageSlices()
{
//    M3DFigure * figure;
	double xpos, ypos, zpos;
    double xExt, yExt, zExt;
    unsigned char * plane;
	Image3D * image;

#ifdef DEBUG
	cout << "P3DView::drawImageSlices()" << endl;
#endif
    if(control == NULL)
        return;

    ImagePlanes * imagePlanes = control->getImagePlanesPtr();

    if(imagePlanes == NULL)
        return;

    int size = imagePlanes->getImagePlaneSize();
    if(size <= 0)
        return;

	// Get model-space coordinates of current slices
	if (imagePlanes->getCutPlaneModelPositions(xpos, ypos, zpos) == false)
		return;
	image = imagePlanes->getImagePtr();
    xExt = image->getXModelExtent();	// Model-space axis lengths
    yExt = image->getYModelExtent();
    zExt = image->getZModelExtent();
/*	cout << "drawImageSlices(): xpos = " << xpos << "  ypos = " << ypos
		<< "  zpos = " << zpos << endl;
	cout << "drawImageSlices(): xExt = " << xExt << "  yExt = " << yExt
		<< "  zExt = " << zExt << endl;
*/
/*
	initBoundaries();

	if (surfaceLevelChanged)
		setBoundariesLevel();

    M3DObject * object = control->getObjectPtr();
	if (object != NULL)
	{
		for (int figNum = 0; figNum < object->getFigureCount(); figNum++)
		{
			xBoundary[figNum].init2(object, surfaceLevel);
			yBoundary[figNum].init2(object, surfaceLevel);
			zBoundary[figNum].init2(object, surfaceLevel);
		}
	}
*/
    if(xSliceOn && !trackball.isTracking())
    {
        // Draw image plane
        plane = imagePlanes->getXCutPlane(imageStable);
        glBindTexture(GL_TEXTURE_2D, textureNames[0]);
        glTexImage2D(GL_TEXTURE_2D, 0, 1, size, size, 0, GL_LUMINANCE,
            GL_UNSIGNED_BYTE, plane);
/*
        // Draw object boundary
		if (drawBoundary && (surfaceStyle == M3D_NONE || surfaceStyle == M3D_POINTCLOUD)
			&& object != NULL)
		{
//			if (surfaceLevelChanged || imagePlanes->xRedrawStatus()) {
				Point origin(xpos, 0.0, 0.0);
				Point corner(xpos, 1.0, 0.0);
				Point end(xpos, 1.0, 1.0);
				for (int figureNum = 0; figureNum < object->getFigureCount(); figureNum++)
                {
                    figure = object->getFigurePtr(figureNum);
                    if(figure == NULL || !(figure->getVisibility()))
                        continue;

					xBoundary[figureNum].run(figureNum, origin, corner, end, size, size);
                }
				xBoundaryDrawingMode = false;
//			}
//			else {
//				// Assure CutPlanes::reset() is only called once per change of the slider
//				if (! xBoundaryDrawingMode) {
//					for (int figureNum = 0; figureNum < object->getFigureCount(); figureNum++)
//						xBoundary[figureNum].reset();
//					xBoundaryDrawingMode = true;
//				}
//			}
		}
*/
    }

    if(ySliceOn && !trackball.isTracking())
    {
        // Draw image plane
        plane = imagePlanes->getYCutPlane(imageStable);
        glBindTexture(GL_TEXTURE_2D, textureNames[1]);
        glTexImage2D(GL_TEXTURE_2D, 0, 1, size, size, 0, GL_LUMINANCE,
            GL_UNSIGNED_BYTE, plane);
/*
        // Draw object boundary
		if (drawBoundary && (surfaceStyle == M3D_NONE || surfaceStyle == M3D_POINTCLOUD)
			&& object != NULL)
		{
//			if (surfaceLevelChanged || imagePlanes->yRedrawStatus()) {
				Point origin(0.0, ypos, 0.0);
				Point corner(1.0, ypos, 0.0);
				Point end(1.0, ypos, 1.0);
				for (int figureNum = 0; figureNum < object->getFigureCount(); figureNum++)
                {
                    figure = object->getFigurePtr(figureNum);
                    if(figure == NULL || !(figure->getVisibility()))
                        continue;

                    yBoundary[figureNum].run(figureNum, origin, corner, end, size, size);
                }
				yBoundaryDrawingMode = false;
//			}
//			else {
//				// Assure CutPlanes::reset() is only called once per change of the slider
//				if (! yBoundaryDrawingMode) {
//					for (int figureNum = 0; figureNum < object->getFigureCount(); figureNum++)
//						yBoundary[figureNum].reset();
//					yBoundaryDrawingMode = true;
//				}
//			}
		}
*/
    }

    if(zSliceOn && !trackball.isTracking())
    {
        // Draw image plane
        plane = imagePlanes->getZCutPlane(imageStable);
        glBindTexture(GL_TEXTURE_2D, textureNames[2]);
        glTexImage2D(GL_TEXTURE_2D, 0, 1, size, size, 0, GL_LUMINANCE,
            GL_UNSIGNED_BYTE, plane);
/*
        // Draw object boundary
		if (drawBoundary && (surfaceStyle == M3D_NONE || surfaceStyle == M3D_POINTCLOUD)
			&& object != NULL)
		{
//			if (surfaceLevelChanged || imagePlanes->zRedrawStatus()) {
				Point origin(0.0, 0.0, zpos);
				Point corner(0.0, 1.0, zpos);
				Point end(1.0, 1.0, zpos);
				for (int figureNum = 0; figureNum < object->getFigureCount(); figureNum++)
                {
                    figure = object->getFigurePtr(figureNum);
                    if(figure == NULL || !(figure->getVisibility()))
                        continue;

					zBoundary[figureNum].run(figureNum, origin, corner, end, size, size);
                }
				zBoundaryDrawingMode = false;
//			}
//			else {
//				// Assure CutPlanes::reset() is only called once per change of the slider
//				if (! zBoundaryDrawingMode) {
//					for (int figureNum = 0; figureNum < object->getFigureCount(); figureNum++)
//						zBoundary[figureNum].reset();
//					zBoundaryDrawingMode = true;
//				}
//			}
		}
*/
    }

/*  if(arbSliceOn && !trackball.isTracking())
    {
        Vector3D ll, lr, ul, ur;
        texPlanes->getArbCutPlaneCoord(ll, lr, ul, ur);

        plane = texPlanes->getArbCutPlane();
        glBindTexture(GL_TEXTURE_2D, texnames[3]);
        glTexImage2D(GL_TEXTURE_2D, 0, 1, size, size, 0, GL_LUMINANCE,
            GL_UNSIGNED_BYTE, plane);
    }
*/
    glColor3f(1.0, 1.0, 1.0);
    glEnable(GL_TEXTURE_2D);
    glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);

    glEnable(GL_POLYGON_OFFSET_FILL);
    glPolygonOffset(1.0, 1.0);
	// Note: the following 3 blocks determine the handedness of
	// the display of images slices.
    if(xSliceOn)
    {
        glBindTexture(GL_TEXTURE_2D, textureNames[0]);

		// Generate tile onto which image is texture mapped
        glBegin(GL_QUADS);
        glTexCoord2d(0.0, 0.0); glVertex3d(xpos, 0.0, 0.0);
        glTexCoord2d(1.0, 0.0); glVertex3d(xpos, yExt, 0.0);
        glTexCoord2d(1.0, 1.0); glVertex3d(xpos, yExt, zExt);
        glTexCoord2d(0.0, 1.0); glVertex3d(xpos, 0.0, zExt);
        glEnd();
    }

    if(ySliceOn)
    {
        glBindTexture(GL_TEXTURE_2D, textureNames[1]);

		// Generate tile onto which image is texture mapped
        glBegin(GL_QUADS);
        glTexCoord2d(0.0, 0.0); glVertex3d(0.0, ypos, 0.0);
        glTexCoord2d(1.0, 0.0); glVertex3d(xExt, ypos, 0.0);
        glTexCoord2d(1.0, 1.0); glVertex3d(xExt, ypos, zExt);
        glTexCoord2d(0.0, 1.0); glVertex3d(0.0, ypos, zExt);
        glEnd();
    }

    if(zSliceOn)
    {
        glBindTexture(GL_TEXTURE_2D, textureNames[2]);

		// Generate tile onto which image is texture mapped
        glBegin(GL_QUADS);
        glTexCoord2d(0.0, 0.0); glVertex3d(0.0, 0.0, zpos);
        glTexCoord2d(1.0, 0.0); glVertex3d(xExt, 0.0, zpos);
        glTexCoord2d(1.0, 1.0); glVertex3d(xExt, yExt, zpos);
        glTexCoord2d(0.0, 1.0); glVertex3d(0.0, yExt, zpos);
        glEnd();
	}

			// AGG: Add optional drawing of image borders here			**********************

/*  if(arbSliceOn && object->getOriginPrimitivePtr() != NULL)
    {
        glBindTexture(GL_TEXTURE_2D, texnames[3]);

        Vector3D ll, lr, ul, ur;
        texPlanes->getArbCutPlaneCoord(ll, lr, ul, ur);

        glBegin(GL_QUADS);
        glTexCoord2f(0.0, 0.0); glVertex3f(ll.getX(), ll.getY(), ll.getZ());
        glTexCoord2f(1.0, 0.0); glVertex3f(lr.getX(), lr.getY(), lr.getZ());
        glTexCoord2f(1.0, 1.0); glVertex3f(ur.getX(), ur.getY(), ur.getZ());
        glTexCoord2f(0.0, 1.0); glVertex3f(ul.getX(), ul.getY(), ul.getZ());
        glEnd();

        glDisable(GL_TEXTURE_2D);

        glLineWidth(2.0f);
        glColor3f(1.0, 1.0, 1.0);

        glBegin(GL_LINE_STRIP);
        glVertex3f(ll.getX(), ll.getY(), ll.getZ());
        glVertex3f(lr.getX(), lr.getY(), lr.getZ());
        glVertex3f(ur.getX(), ur.getY(), ur.getZ());
        glVertex3f(ul.getX(), ul.getY(), ul.getZ());
        glVertex3f(ll.getX(), ll.getY(), ll.getZ());
        glEnd();

        glLineWidth(1.0f);
    }
*/
    glDisable(GL_POLYGON_OFFSET_FILL);
    glDisable(GL_TEXTURE_2D);

	// Draw lines to show the intersection of image planes
    if(xSliceOn && ySliceOn)
    {
        glBegin(GL_LINES);
	    glColor3d(0.0, 0.0, 1.0);	// Blue for Z
        glVertex3d(xpos, ypos, 0.0);
        glVertex3d(xpos, ypos, zExt);
        glEnd();

/*        glBegin(GL_QUADS);
        glVertex3d(xpos, 0.0, 0.0);
        glVertex3d(xpos, yExt, 0.0);
        glVertex3d(xpos, yExt, zExt);
        glVertex3d(xpos, 0.0, zExt);
        glEnd();
*/
    }

    if(xSliceOn && zSliceOn)
    {
        glBegin(GL_LINES);
	    glColor3d(1.0, 0.0, 1.0);	// Magenta for Y
        glVertex3d(xpos, 0.0, zpos);
        glVertex3d(xpos, yExt, zpos);
        glEnd();
    }

    if(ySliceOn && zSliceOn)
    {
        glBegin(GL_LINES);
	    glColor3d(1.0, 0.0, 0.0);	// Red for X
        glVertex3d(0.0, ypos, zpos);
        glVertex3d(xExt, ypos, zpos);
        glEnd();
    }

    glColor3d(1.0, 1.0, 1.0);
}

void P3DView::updateAllCutPlanes()
{
	if (modelIsBeingEdited)
		return;

	if (modelWasEdited)
		modelWasEdited = false;
	else
		return;

	// Draw the cut plane images
	if (bNCutPlaneOn)
		(void) callback->redrawCutPlane(ATOM_PLANE);
	if (bBperpCutPlaneOn)
		(void) callback->redrawCutPlane(CREST_PLANE);
	if (bperpNCutPlaneOn)
		(void) callback->redrawCutPlane(BPERP_N_PLANE);
	if (bperpY0CutPlaneOn)
		(void) callback->redrawCutPlane(STARBOARD_PLANE);
	if (bperpY1CutPlaneOn)
		(void) callback->redrawCutPlane(PORT_PLANE);
	if (involutesCutPlaneOn)
		(void) callback->redrawCutPlane(INVOLUTES_PLANE);
}

void P3DView::drawAllCutPlanes()
{
	// Draw the cut plane images
	if (bNCutPlaneOn && callback->getCutPlaneMode(ATOM_PLANE)
		!= P3DUserInterfaceCallback::SEPARATE_WIN)
			drawCutPlane(ATOM_PLANE);
	if (bBperpCutPlaneOn && callback->getCutPlaneMode(CREST_PLANE)
		!= P3DUserInterfaceCallback::SEPARATE_WIN)
			drawCutPlane(CREST_PLANE);
	if (bperpNCutPlaneOn && callback->getCutPlaneMode(BPERP_N_PLANE)
		!= P3DUserInterfaceCallback::SEPARATE_WIN)
			drawCutPlane(BPERP_N_PLANE);
	if (bperpY0CutPlaneOn && callback->getCutPlaneMode(STARBOARD_PLANE)
		!= P3DUserInterfaceCallback::SEPARATE_WIN)
			drawCutPlane(STARBOARD_PLANE);
	if (bperpY1CutPlaneOn && callback->getCutPlaneMode(PORT_PLANE)
		!= P3DUserInterfaceCallback::SEPARATE_WIN)
			drawCutPlane(PORT_PLANE);
	if (involutesCutPlaneOn && callback->getCutPlaneMode(INVOLUTES_PLANE)
		!= P3DUserInterfaceCallback::SEPARATE_WIN)
			drawCutPlane(INVOLUTES_PLANE);
}

// This function draws a Open GL display list
void P3DView::drawBoundaryCurves(bool changed)
{
	M3DObject * object;
    M3DFigure * figure;
	float * color;
	double xpos, ypos, zpos;
    int figNum;


#ifdef DEBUG
	cout << "P3DView::drawBoundaryCurves()" << endl;
#endif

	if (! drawBoundary)
		return;

    if(control == NULL)
        return;

    object = control->getObjectPtr();
	if (object == NULL)
		return;

//	cout << "cutPlaneBoundaryWidth = " << cutPlaneBoundaryWidth << endl;
	int numFigures = object->getFigureCount();

    ImagePlanes * imagePlanes = control->getImagePlanesPtr();

    if(imagePlanes == NULL)
        return;

    int size = imagePlanes->getImagePlaneSize();
    if(size <= 0)
        return;

	// Get model-space coordinates of current slices
	if (imagePlanes->getCutPlaneModelPositions(xpos, ypos, zpos) == false)
		return;

	if (changed || xBoundary.size() == 0) {
	    xBoundary.resize(numFigures);
	    yBoundary.resize(numFigures);
	    zBoundary.resize(numFigures);

	    for(figNum = 0; figNum < object->getFigureCount(); figNum++)
	    {
		xBoundary[figNum].init2(object, surfaceLevel);
		yBoundary[figNum].init2(object, surfaceLevel);
		zBoundary[figNum].init2(object, surfaceLevel);
	    }
	}

    if (xSliceOn)
    {
        // Draw object boundary
		if (drawBoundary && (surfaceStyle == M3D_NONE || surfaceStyle == M3D_POINTCLOUD)
			&& object != NULL)
		{
			Point origin(xpos, 0.0, 0.0);
			Point corner(xpos, 1.0, 0.0);
			Point end(xpos, 1.0, 1.0);
			for (int figureNum = 0; figureNum < object->getFigureCount(); figureNum++)
            {
                figure = object->getFigurePtr(figureNum);
                if(figure == NULL || !(figure->getVisibility()))
                    continue;

				xBoundary[figureNum].run(figureNum, origin, corner, end, size, size);
            }
//			xBoundaryDrawingMode = false;
		}
    }
    if (ySliceOn)
    {
        // Draw object boundary
	    if (drawBoundary && (surfaceStyle == M3D_NONE || surfaceStyle == M3D_POINTCLOUD)
		    && object != NULL)
	    {
		    Point origin(0.0, ypos, 0.0);
		    Point corner(1.0, ypos, 0.0);
		    Point end(1.0, ypos, 1.0);
		    for (int figureNum = 0; figureNum < object->getFigureCount(); figureNum++)
            {
                figure = object->getFigurePtr(figureNum);
                if(figure == NULL || !(figure->getVisibility()))
                    continue;

                yBoundary[figureNum].run(figureNum, origin, corner, end, size, size);
            }
//		    yBoundaryDrawingMode = false;
	    }
    }
    if (zSliceOn)
    {
        // Draw object boundary
		if (drawBoundary && (surfaceStyle == M3D_NONE || surfaceStyle == M3D_POINTCLOUD)
			&& object != NULL)
		{
			Point origin(0.0, 0.0, zpos);
			Point corner(0.0, 1.0, zpos);
			Point end(1.0, 1.0, zpos);
			for (int figureNum = 0; figureNum < object->getFigureCount(); figureNum++)
            {
                figure = object->getFigurePtr(figureNum);
                if(figure == NULL || !(figure->getVisibility()))
                    continue;

				zBoundary[figureNum].run(figureNum, origin, corner, end, size, size);
            }
//			zBoundaryDrawingMode = false;
		}
    }

    surfaceLevelChanged = false;


#ifdef COLOR_CONTOURS
	// When COLOR_CONTOURS is defined, I want to use the colors from the
	// object surface visualizer to shade the contours

	// First I will draw a thick line with a fixed color
	// then I'll draw a thinner line in the right color
	if (cutPlaneBoundaryWidth > 1)
		  glLineWidth((float) (1.5 * cutPlaneBoundaryWidth)); 
	glColor3f(0.25, 0.25, 0.25); // From SPIE_PAPER
	glBegin(GL_LINES);
	if (xSliceOn) {
		for (int figureNum = 0; figureNum < numFigures; figureNum++) {
            figure = object->getFigurePtr(figureNum);
            if(figure == NULL || !(figure->getVisibility()))
                continue;
			if (xBoundary[figureNum].boundary_point_count() > 0) {
				for (int j = 1; j < xBoundary[figureNum].boundary_point_count(); j += 2) {
					Point p = xBoundary[figureNum].boundary_point(j - 1);
					Point q = xBoundary[figureNum].boundary_point(j);
					glVertex3d(p.getX(), p.getY(), p.getZ());
					glVertex3d(q.getX(), q.getY(), q.getZ());
				}
			}
		}
	}

	if (ySliceOn) {
		for (int figureNum = 0; figureNum < numFigures; figureNum++) {
            figure = object->getFigurePtr(figureNum);
            if(figure == NULL || !(figure->getVisibility()))
                continue;
			if (yBoundary[figureNum].boundary_point_count() > 0) {
				for (int j = 1; j < yBoundary[figureNum].boundary_point_count(); j += 2) {
					Point p = yBoundary[figureNum].boundary_point(j - 1);
					Point q = yBoundary[figureNum].boundary_point(j);
					glVertex3d(p.getX(), p.getY(), p.getZ());
					glVertex3d(q.getX(), q.getY(), q.getZ());
				}
			}
		}
	}

	if (zSliceOn) {
		for (int figureNum = 0; figureNum < numFigures; figureNum++) {
            figure = object->getFigurePtr(figureNum);
            if(figure == NULL || !(figure->getVisibility()))
                continue;
			if (zBoundary[figureNum].boundary_point_count() > 0) {
				for (int j = 1; j < zBoundary[figureNum].boundary_point_count(); j += 2) {
					Point p = zBoundary[figureNum].boundary_point(j - 1);
					Point q = zBoundary[figureNum].boundary_point(j);
					glVertex3d(p.getX(), p.getY(), p.getZ());
					glVertex3d(q.getX(), q.getY(), q.getZ());
				}
			}
		}
	}

	glEnd();


	if (cutPlaneBoundaryWidth > 1)
		glLineWidth((float) (cutPlaneBoundaryWidth));

	int numColors = surfaceColorMap.size();

	std::vector<float> & r = *(surfaceColorMap.red());
	std::vector<float> & g = *(surfaceColorMap.green());
	std::vector<float> & b = *(surfaceColorMap.blue());

#endif /* COLOR_CONTOURS */

	glBegin(GL_LINES);
	if (xSliceOn) {
		for (int figureNum = 0; figureNum < numFigures; figureNum++) {
            figure = object->getFigurePtr(figureNum);
            if(figure == NULL || !(figure->getVisibility()))
                continue;
			if (xBoundary[figureNum].boundary_point_count() > 0) {
				color = xBoundary[figureNum].boundary_color();
				glColor3fv(color);
				for (int j = 1; j < xBoundary[figureNum].boundary_point_count(); j += 2) {
					Point p = xBoundary[figureNum].boundary_point(j - 1);
					Point q = xBoundary[figureNum].boundary_point(j);
#ifdef COLOR_CONTOURS
					float surfVal = objectSurfaceVisualizer->getColorValue(
						figureNum, xBoundary[figureNum].bpoint_index(j - 1), 0.5);
					double dColorIdx = (numColors - 1) * surfVal;
					dColorIdx = (dColorIdx < 0) ? 0 : ((dColorIdx >= (numColors-1)) ? (numColors-1) : dColorIdx);	// AGG: needed?
					int iColorIdx = int(dColorIdx);
					double colorIdxOffset = dColorIdx - iColorIdx;
					glColor3f((1 - colorIdxOffset)*r[iColorIdx] + (colorIdxOffset*r[iColorIdx + 1]),
						(1 - colorIdxOffset)*g[iColorIdx] + (colorIdxOffset*g[iColorIdx + 1]),
						(1 - colorIdxOffset)*b[iColorIdx] + (colorIdxOffset*b[iColorIdx + 1]));
					
#endif /* COLOR_CONTOURS */
					glVertex3d(p.getX(), p.getY(), p.getZ());
#ifdef COLOR_CONTOURS
					surfVal = objectSurfaceVisualizer->getColorValue(figureNum,
						xBoundary[figureNum].bpoint_index(j), 0.5);
					dColorIdx = (numColors - 1) * surfVal;
					dColorIdx = (dColorIdx < 0) ? 0
						: ((dColorIdx >= (numColors-1)) ? (numColors-1) : dColorIdx);
					iColorIdx = int(dColorIdx);
					colorIdxOffset = dColorIdx - iColorIdx;
					glColor3f((1 - colorIdxOffset)*r[iColorIdx] + (colorIdxOffset*r[iColorIdx + 1]),
						(1 - colorIdxOffset)*g[iColorIdx] + (colorIdxOffset*g[iColorIdx + 1]),
						(1 - colorIdxOffset)*b[iColorIdx] + (colorIdxOffset*b[iColorIdx + 1]));
					
#endif /* COLOR_CONTOURS */
					glVertex3d(q.getX(), q.getY(), q.getZ());
				}
			}
		}
	}

	if (ySliceOn) {
		for (int figureNum = 0; figureNum < numFigures; figureNum++) {
            figure = object->getFigurePtr(figureNum);
            if(figure == NULL || !(figure->getVisibility()))
                continue;
			if (yBoundary[figureNum].boundary_point_count() > 0) {
				color = yBoundary[figureNum].boundary_color();
				glColor3fv(color);
				for (int j = 1; j < yBoundary[figureNum].boundary_point_count(); j += 2) {
					Point p = yBoundary[figureNum].boundary_point(j - 1);
					Point q = yBoundary[figureNum].boundary_point(j);
#ifdef COLOR_CONTOURS
					float surfVal = objectSurfaceVisualizer->getColorValue(
						figureNum, yBoundary[figureNum].bpoint_index(j - 1), 0.5);
					double dColorIdx = (numColors - 1) * surfVal;
					dColorIdx = (dColorIdx < 0) ? 0
						: ((dColorIdx >= (numColors-1)) ? (numColors-1) : dColorIdx);
					int iColorIdx = int(dColorIdx);
					double colorIdxOffset = dColorIdx - iColorIdx;
					glColor3f((1 - colorIdxOffset)*r[iColorIdx] + (colorIdxOffset*r[iColorIdx + 1]),
						(1 - colorIdxOffset)*g[iColorIdx] + (colorIdxOffset*g[iColorIdx + 1]),
						(1 - colorIdxOffset)*b[iColorIdx] + (colorIdxOffset*b[iColorIdx + 1]));
					
#endif /* COLOR_CONTOURS */
					glVertex3d(p.getX(), p.getY(), p.getZ());
#ifdef COLOR_CONTOURS
					surfVal = objectSurfaceVisualizer->getColorValue(
						figureNum, yBoundary[figureNum].bpoint_index(j), 0.5);
					dColorIdx = (numColors - 1) * surfVal;
					dColorIdx = (dColorIdx < 0) ? 0
						: ((dColorIdx >= (numColors-1)) ? (numColors-1) : dColorIdx);
					iColorIdx = int(dColorIdx);
					colorIdxOffset = dColorIdx - iColorIdx;
					glColor3f((1 - colorIdxOffset)*r[iColorIdx] + (colorIdxOffset*r[iColorIdx + 1]),
						(1 - colorIdxOffset)*g[iColorIdx] + (colorIdxOffset*g[iColorIdx + 1]),
						(1 - colorIdxOffset)*b[iColorIdx] + (colorIdxOffset*b[iColorIdx + 1]));
					
#endif /* COLOR_CONTOURS */
					glVertex3d(q.getX(), q.getY(), q.getZ());
				}
			}
		}
	}

	if (zSliceOn) {
		for (int figureNum = 0; figureNum < numFigures; figureNum++) {
            figure = object->getFigurePtr(figureNum);
            if(figure == NULL || !(figure->getVisibility()))
                continue;
			if (zBoundary[figureNum].boundary_point_count() > 0) {
				color = zBoundary[figureNum].boundary_color();
				glColor3fv(color);
				for (int j = 1; j < zBoundary[figureNum].boundary_point_count(); j += 2) {
					Point p = zBoundary[figureNum].boundary_point(j - 1);
					Point q = zBoundary[figureNum].boundary_point(j);
#ifdef COLOR_CONTOURS
					float surfVal = objectSurfaceVisualizer->getColorValue(
						figureNum, zBoundary[figureNum].bpoint_index(j - 1), 0.5);
					double dColorIdx = (numColors - 1) * surfVal;
					dColorIdx = (dColorIdx < 0) ? 0
						: ((dColorIdx >= (numColors-1)) ? (numColors-1) : dColorIdx);
					int iColorIdx = int(dColorIdx);
					double colorIdxOffset = dColorIdx - iColorIdx;
					glColor3f((1 - colorIdxOffset)*r[iColorIdx] + (colorIdxOffset*r[iColorIdx + 1]),
						(1 - colorIdxOffset)*g[iColorIdx] + (colorIdxOffset*g[iColorIdx + 1]),
						(1 - colorIdxOffset)*b[iColorIdx] + (colorIdxOffset*b[iColorIdx + 1]));
					
#endif /* COLOR_CONTOURS */
					glVertex3d(p.getX(), p.getY(), p.getZ());
#ifdef COLOR_CONTOURS
					surfVal = objectSurfaceVisualizer->getColorValue(figureNum,
						zBoundary[figureNum].bpoint_index(j), 0.5);
					dColorIdx = (numColors - 1) * surfVal;
					dColorIdx = (dColorIdx < 0) ? 0
						: ((dColorIdx >= (numColors-1)) ? (numColors-1) : dColorIdx);
					iColorIdx = int(dColorIdx);
					colorIdxOffset = dColorIdx - iColorIdx;
					glColor3f((1 - colorIdxOffset)*r[iColorIdx] + (colorIdxOffset*r[iColorIdx + 1]),
						(1 - colorIdxOffset)*g[iColorIdx] + (colorIdxOffset*g[iColorIdx + 1]),
						(1 - colorIdxOffset)*b[iColorIdx] + (colorIdxOffset*b[iColorIdx + 1]));
					
#endif /* COLOR_CONTOURS */
					glVertex3d(q.getX(), q.getY(), q.getZ());
				}
			}
		}
	}

	glEnd();
	if (cutPlaneBoundaryWidth > 1)
		glLineWidth(1.0f);
}

void P3DView::drawFigureNames()
{
    M3DObject * object;
    M3DObjectRenderer objectRenderer;

    if (control == NULL)
        return;

    object = control->getObjectPtr();
    if (object == NULL)
        return;

    objectRenderer.setObjectPtr(object);
    objectRenderer.setCurrentScale(trackball.getCurrScale());

    objectRenderer.renderFigureNames();
}

// Draw spoken-end landmarks
void P3DView::drawObjectLandmarks()
{
    M3DObject * object;
    M3DObjectRenderer objectRenderer;
    M3DPrimitiveRenderer primRenderer;

    if (control == NULL)
        return;

    object = control->getObjectPtr();
    if (object == NULL)
        return;

    objectRenderer.setObjectPtr(object);
    objectRenderer.setCurrentScale(trackball.getCurrScale());
    objectRenderer.setAtomVisibility(true);
    objectRenderer.renderLandmarks();
}

void P3DView::drawImageLandmarks()
{
	double xpos, ypos, zpos;
	double d;
	M3DObject * object;
	int figureId, j;
	M3DFigure * figure;
	bool flag;

#ifdef DEBUG
	cout << "P3DView::drawImageLandmarks()" << endl;
#endif
    if(control == NULL)
        return;
	if (control->getImageLandmarksCount() == 0)
		return;

	if ((! xSliceOn) && (! ySliceOn) && (! zSliceOn))
		return;

    ImagePlanes * imagePlanes = control->getImagePlanesPtr();
    if(imagePlanes == NULL)
        return;
    Image3D * image = imagePlanes->getImagePtr();
//    if(image == NULL)
//        return;

    object = control->getObjectPtr();
	if (object == NULL)
		return;

	// Get model-space coordinates of current slices
	if (imagePlanes->getCutPlaneModelPositions(xpos, ypos, zpos) == false)
		return;

//	cout << "drawImageLandmarks(): pos = (" << xpos << ", " << ypos
//	<< ", " << zpos << ')' << endl;

	vector<imageLandmark_t> & imageLandmarks = control->landmarkCoordinates();

	glBegin(GL_LINES);
	glColor3fv(defaultLandmarkColor);
	flag = false;

	for (j = 0; j < imageLandmarks.size(); j++)
	{
		int slice[3];

		// Ignore unpaired landmarks, regardless of the figure number
		if (! control->landmarkIsPaired(imageLandmarks[j].figureId, j))
			continue;

		Point p(imageLandmarks[j].coords);
		image->modelToImageCoordinates(p, slice);

		if (xSliceOn) {
			double x = imagePlanes->getXCutPlaneIndex();
			d = x - slice[0];
			if (d < 0.5 && d >= -0.5) {
				figureId = imageLandmarks[j].figureId;
				if (figureId >= 0) {
					figure = object->getFigurePtr(figureId);
					if (figure->getMarkedLandmark() == imageLandmarks[j].number) {
						glEnd();
						glLineWidth(markedLandmarkWidth);
						glBegin(GL_LINES);
						flag = true;
					}
					glVertex3d(xpos, imageLandmarks[j].coords[1] - LANDMARK_SIZE,
						imageLandmarks[j].coords[2] - LANDMARK_SIZE);
					glVertex3d(xpos, imageLandmarks[j].coords[1] + LANDMARK_SIZE,
						imageLandmarks[j].coords[2] + LANDMARK_SIZE);
					glVertex3d(xpos, imageLandmarks[j].coords[1] + LANDMARK_SIZE,
						imageLandmarks[j].coords[2] - LANDMARK_SIZE);
					glVertex3d(xpos, imageLandmarks[j].coords[1] - LANDMARK_SIZE,
						imageLandmarks[j].coords[2] + LANDMARK_SIZE);
#ifdef DEBUG
					cout << "Drawing image landmark on X slice at (" << xpos << ", "
						<< imageLandmarks[j].coords[1] << ", " << imageLandmarks[j].coords[2]
						<< ")\n";
#endif
					if (flag) {
						glEnd();
						glLineWidth(defaultLandmarkWidth);
						glBegin(GL_LINES);
						flag = false;
					}
				}
			}
		}

		if (ySliceOn) {
			double y = imagePlanes->getYCutPlaneIndex();
			d = y - slice[1];
			if (d < 0.5 && d >= -0.5) {
				figureId = imageLandmarks[j].figureId;
				if (figureId >= 0) {
					figure = object->getFigurePtr(figureId);
					if (figure->getMarkedLandmark() == imageLandmarks[j].number) {
						glEnd();
						glLineWidth(markedLandmarkWidth);
						glBegin(GL_LINES);
						flag = true;
					}
					glVertex3d(imageLandmarks[j].coords[0] - LANDMARK_SIZE, ypos,
						imageLandmarks[j].coords[2] - LANDMARK_SIZE);
					glVertex3d(imageLandmarks[j].coords[0] + LANDMARK_SIZE, ypos,
						imageLandmarks[j].coords[2] + LANDMARK_SIZE);
					glVertex3d(imageLandmarks[j].coords[0] + LANDMARK_SIZE, ypos,
						imageLandmarks[j].coords[2] - LANDMARK_SIZE);
					glVertex3d(imageLandmarks[j].coords[0] - LANDMARK_SIZE, ypos,
						imageLandmarks[j].coords[2] + LANDMARK_SIZE);
#ifdef DEBUG
					cout << "Drawing image landmark on Y slice at (" 
						<< imageLandmarks[j].coords[0] << ", " << ypos << ", "
						<< imageLandmarks[j].coords[2] << ")\n";
#endif
					if (flag) {
						glEnd();
						glLineWidth(defaultLandmarkWidth);
						glBegin(GL_LINES);
						flag = false;
					}
				}
			}
		}

		if (zSliceOn) {
			double z = imagePlanes->getZCutPlaneIndex();	// Rounded slice number
			d = z - slice[2];
			if (d < 0.5 && d >= -0.5) {
				figureId = imageLandmarks[j].figureId;
				if (figureId >= 0) {
					figure = object->getFigurePtr(figureId);
					if (figure->getMarkedLandmark() == imageLandmarks[j].number) {
						glEnd();
						glLineWidth(markedLandmarkWidth);
						glBegin(GL_LINES);
						flag = true;
					}
					glVertex3d(imageLandmarks[j].coords[0] - LANDMARK_SIZE,
						imageLandmarks[j].coords[1] - LANDMARK_SIZE, zpos);
					glVertex3d(imageLandmarks[j].coords[0] + LANDMARK_SIZE,
						imageLandmarks[j].coords[1] + LANDMARK_SIZE, zpos);
					glVertex3d(imageLandmarks[j].coords[0] + LANDMARK_SIZE,
						imageLandmarks[j].coords[1] - LANDMARK_SIZE, zpos);
					glVertex3d(imageLandmarks[j].coords[0] - LANDMARK_SIZE,
						imageLandmarks[j].coords[1] + LANDMARK_SIZE, zpos);
#ifdef DEBUG
					cout << "Drawing image landmark on Z slice at ("
						<< imageLandmarks[j].coords[0] << ", "
						<< imageLandmarks[j].coords[1] << ", " << zpos << ")\n";
#endif
					if (flag) {
						glEnd();
						glLineWidth(defaultLandmarkWidth);
						glBegin(GL_LINES);
						flag = false;
					}
				}
			}
		}
	}
    glEnd();
}

void P3DView::drawCutPlane(cut_plane_t planeType)
{
	P3DCutPlaneView * view;

#ifdef DEBUG
	cout << "P3DView::drawCutPlane()\n";
#endif
	view = callback->getView(planeType);

	Image3D * vis = view->visualization();
	if (vis == NULL)
		return;

	// Determine corners of the plane to be displayed
	Vector3D p1 = view->origin();
	Vector3D p2 = view->corner();
	Vector3D p3 = view->end();
	Vector3D p4 = p1 + (p3 - p2);

	bool half_plane = false;
	if (planeType == (int) bperp_y0 || planeType == (int) bperp_y1)
	{
		p3 += (p2 - p3) * 0.5;
		p4 += (p1 - p4) * 0.5;
		half_plane = true;
	}

#ifdef DRAWING_DEBUG
	cout << "p1 = "; p1.print();
	cout << "p2 = "; p2.print();
	cout << "p3 = "; p3.print();
	cout << std::flush;
#endif

    glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);

    glEnable(GL_POLYGON_OFFSET_FILL);
    glPolygonOffset(1.0, 1.0);

	texture_done = false;		// AGG: fix this!

	bool smooth = globalControl->readBool(SmoothImages);
	view->paint_image((void *)&(cutPlaneTextureNames[0]), &texture_done,
		&p4, &p3, &p2, &p1, smooth, half_plane);
#ifdef DRAWING_DEBUG
	cout << "Texture # " << cutPlaneTextureNames[0] << endl;
#endif

    glDisable(GL_POLYGON_OFFSET_FILL);

	// Draw a box around the visualization
	glBegin(GL_LINES);
	glVertex3d(p1.getX(), p1.getY(), p1.getZ());
	glVertex3d(p2.getX(), p2.getY(), p2.getZ());
	glVertex3d(p2.getX(), p2.getY(), p2.getZ());
	glVertex3d(p3.getX(), p3.getY(), p3.getZ());
	glVertex3d(p3.getX(), p3.getY(), p3.getZ());
	glVertex3d(p4.getX(), p4.getY(), p4.getZ());
	glVertex3d(p4.getX(), p4.getY(), p4.getZ());
	glVertex3d(p1.getX(), p1.getY(), p1.getZ());
	glEnd();
}

void P3DView::setFiguralVisibility(int figureID, bool yesno) {
    M3DObject * object;
    M3DFigure * fig;

	if (control == NULL)
        return;

    object = control->getObjectPtr();
	if (object == NULL)
		return;

	if (figureID < object->getFigureCount()) {
        fig = object->getFigurePtr(figureID);
		fig->setVisibility(yesno);
	}
}

void P3DView::setAllFiguresVisibility(bool yesno) {
    M3DObject * object;
    M3DFigure * fig;

	if (control == NULL)
        return;

    object = control->getObjectPtr();
	if (object == NULL)
		return;

	int numFigures = object->getFigureCount();
	for (int figureID = 0; figureID < numFigures; figureID++) {
        fig = object->getFigurePtr(figureID);
		fig->setVisibility(yesno);
	}
}

void P3DView::setConstraintsVisible(bool yesno) {
	constraintArrows = yesno;
	updateObject();
}

void P3DView::getCenter(Vector3D & v)
{
	GLint viewport[4];
	GLdouble modelview[16];
	GLdouble projection[16];
	double c[3];

	glGetIntegerv(GL_VIEWPORT, viewport);
	glGetDoublev(GL_MODELVIEW_MATRIX, modelview);
	glGetDoublev(GL_PROJECTION_MATRIX, projection);

	double centerX = viewport[2]/2.0;
	double centerY = viewport[3]/2.0;
	centerY = viewport[3] - centerY - 1.0;
	gluUnProject(centerX, centerY, 0.5, modelview, projection, viewport,
		&c[0], &c[1], &c[2]);
	v.set(c);
}

void P3DView::translateToCenter(Vector3D & v)
{
	Vector3D c;

	getCenter(c);
	v = c - v;
	trackball.translateView(v.getX(), v.getY(), v.getZ());

	redraw();
}

void P3DView::rotateAboutCenter(double angle)
{
	GLint viewport[4];
	GLdouble modelview[16];
	GLdouble projection[16];
	double a[3];
	double b[3];

	glGetIntegerv(GL_VIEWPORT, viewport);
	glGetDoublev(GL_MODELVIEW_MATRIX, modelview);
	glGetDoublev(GL_PROJECTION_MATRIX, projection);

	double centerX = viewport[2]/2.0;
	double centerY = viewport[3]/2.0;
	centerY = viewport[3] - centerY - 1.0;
	gluUnProject(centerX, centerY, 1.0, modelview, projection, viewport,
		&a[0], &a[1], &a[2]);
	gluUnProject(centerX, centerY, 0.0, modelview, projection, viewport,
		&b[0], &b[1], &b[2]);

	b[0] -= a[0];
	b[1] -= a[1];
	b[2] -= a[2];
	trackball.rotateView(angle, b[0], b[1], b[2]);

	redraw();
}

void P3DView::rotateAboutHorizontal(double angle)
{
	GLint viewport[4];
	GLdouble modelview[16];
	GLdouble projection[16];
	double a[3];
	double b[3];

	glGetIntegerv(GL_VIEWPORT, viewport);
	glGetDoublev(GL_MODELVIEW_MATRIX, modelview);
	glGetDoublev(GL_PROJECTION_MATRIX, projection);

	double centerY = viewport[3]/2.0;
	centerY = viewport[3] - centerY - 1.0;
	gluUnProject(0.0, centerY, 0.5, modelview, projection, viewport,
		&a[0], &a[1], &a[2]);
	gluUnProject(1.0, centerY, 0.5, modelview, projection, viewport,
		&b[0], &b[1], &b[2]);

	b[0] -= a[0];
	b[1] -= a[1];
	b[2] -= a[2];
	trackball.rotateView(angle, b[0], b[1], b[2]);

	redraw();
}

void P3DView::rotateAboutVertical(double angle)
{
	GLint viewport[4];
	GLdouble modelview[16];
	GLdouble projection[16];
	double a[3];
	double b[3];

	glGetIntegerv(GL_VIEWPORT, viewport);
	glGetDoublev(GL_MODELVIEW_MATRIX, modelview);
	glGetDoublev(GL_PROJECTION_MATRIX, projection);

	double centerX = viewport[2]/2.0;
	gluUnProject(centerX, 0.0, 0.5, modelview, projection, viewport,
		&a[0], &a[1], &a[2]);
	gluUnProject(centerX, 1.0, 0.5, modelview, projection, viewport,
		&b[0], &b[1], &b[2]);

	b[0] -= a[0];
	b[1] -= a[1];
	b[2] -= a[2];
	trackball.rotateView(angle, b[0], b[1], b[2]);

	redraw();
}

void P3DView::screenDump()
{
    static int num = 0;
    char filename[16];
    unsigned char * data,
                  * data2;
    FILE * fp;
    int width, height, size;

    width = w() - 1;
    height = h() - 1;

    size = width * height * 3;

    sprintf(filename, "out%.3d.ppm", num);
    data = new unsigned char[size];
    data2 = new unsigned char[size];
    glPixelStorei(GL_PACK_ALIGNMENT, 1);
    glReadPixels(0, 0, width, height, GL_RGB, GL_UNSIGNED_BYTE,
		(GLvoid*)(&(data[0])));

    int j, k, index = 0;
    for(j = width * (height - 1) * 3; j >= 0; j -= width * 3)
    {
        for(k = j; k < j + width * 3; k++)
        {
            data2[k] = data[index];
            index++;
        }
    }

    fp = fopen(filename, "wb");
    if(fp != NULL)
    {
        fprintf(fp, "P6\n%d %d\n255\n", width, height);
        fwrite((void *) data2, 1, size, fp);
        fclose(fp);
    }
	else {
		cout << "Could not open dump file" << endl;
		fl_beep(FL_BEEP_ERROR);
	}
	if (globalVerbosity >= 0)
		cout << "Model window dumped to file " << filename << endl;

    num++;
    delete data;
    delete data2;
}

#ifndef PRODUCTION_VERSION

void P3DView::setBlendingRegionDisplay(bool onOff) {
	displayBlendRegion = onOff;
	updateSurface();
}

#endif	/* PRODUCTION_VERSION */

