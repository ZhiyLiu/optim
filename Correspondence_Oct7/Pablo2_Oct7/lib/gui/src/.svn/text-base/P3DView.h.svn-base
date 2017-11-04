#ifndef P3D_VIEW_H
#define P3D_VIEW_H

#include <FL/Fl.H>
#include <FL/Fl_Window.H>
#include <FL/Fl_Gl_Window.H>
#include <FL/fl_ask.H>

#include "P3DDisplayGlobals.h"
#include "M3DObject.h"
#include "Trackball.h"
#include "M3DBlendedRenderer.h"
#include "M3DObjectSurfaceRenderer.h"
#include "CutPlaneGlobals.h"

#define INDIGO 14        /* FLTK definition of default background color.
						    This definition must not conflict with those
							in P3DUserInterfaceCallback.cpp.
						  */
#define MAX_TOOLTIP_SIZE	16


class P3DUserInterfaceCallback;
class M3DBlendedRenderer;
class TileSetRenderer;

enum P3DSelectionType
{
    SELECT_PRIMITIVE,
    SELECT_FIGURE,
	SHOW_ID
};

enum redraw_t
{
	IgnoreRedraw,
    Redraw,
    UpdateObject
};


class P3DView : public Fl_Gl_Window
{
	friend void timeoutCallback(void *);

public:
    P3DView(int x, int y, int w, int h,
            ViewType _viewType = PERSPECTIVE, char *l = 0);

    ~P3DView();

    void setControl(P3DControl * cntrl);
    void setUserInterfaceCallback(P3DUserInterfaceCallback * cb)
                                 { callback = cb; }

    // Initialize the boundary intersection curves
	void clearBoundaries();

	void useTwoLights(bool yesno) { init = false; twoLights = yesno; }

    void setViewType(ViewType type);
	void setAtomVectorsType(bool toggle);
	void setBVectorsType(int val);
	void setShowExtraAtomVectors(bool toggle);
	void setAtomVectorsWidth(int width);
	void setShowMeshConnectors(bool toggle);
	void setMeshConnectorsType(bool toggle);
	void setMeshConnectorsWidth(int width);
	void setBackgroundColor(float color[3]);
	void setMeshConnectorsColor(float color[3]);
	void setFiguresColorAtoms(bool toggle);

	// Figural visibility controls
	void setFiguralVisibility(int figureID, bool yesno);
	void setAllFiguresVisibility(bool yesno);

	// Constraint cones visibility control
	void setConstraintsVisible(bool yesno);

    // Set whether to display the axes of primitives
    void setAtomVectorsVisibility(bool yesno)
    {
        primitiveVectorsOn = yesno;
    }

    // Set whether to display atoms
    void setAtomsVisibility(bool yesno)
    {
        atomsOn = yesno;
    }

    // Set whether to display medial sheet
    void setMedialMeshVisibility(bool yesno)
    {
        medialSheetOn = yesno;
    }

    // Set rocking angle
    void setRockingAngle(double angle)
    {
        deflectionAngle = angle;
    }

    // Set rocking angle increment (controls speed)
    void setRockingIncrement(double angle)
    {
        incrementalAngle = angle;
    }

	// Full surface rendering
    void setSurfaceStyle(M3DSurfaceStyle style) {
		surfaceStyle = style;
		surfaceDone = false;
        boundaryCurvesDone = false;
		auxTileSetDone = false;
	}
	void setSurfaceLevel(int level);
	void setSurfaceLineWidth(int width) {
		surfaceLineWidth = width;
		surfaceDone = false;
		auxTileSetDone = false;
		gMasterDispListIsValid = false;
	}
	void toggleBlendSubfigures(bool val);
	void toggleDrawBoundary(bool val) {
        drawBoundary = val;
        boundaryCurvesDone = false;
    }
	void boundaryWidth(int w) {
        cutPlaneBoundaryWidth = w;
        boundaryCurvesDone = false;
    }
	void toggleDrawFigureNames(bool val) {
        labelFigures = val;
		surfaceDone = false;
    }
    M3DSurfaceStyle getSurfaceStyle() { return surfaceStyle; }
	int getSurfaceLevel() { return surfaceLevel; }
	bool blendingIsOn() { return blendSurface; }
	void toggleTileType(bool val);

	// Partial surface rendering (for interpenetration constraint setting)
    void setPartialSurfaceStyle(M3DSurfaceStyle style) { partialSurfaceStyle = style; }
	void setPartialSurfaceLevel(int level);
    M3DSurfaceStyle getPartialSurfaceStyle() { return partialSurfaceStyle; }
	int getPartialSurfaceLevel() { return partialSurfaceLevel; }

#ifndef PRODUCTION_VERSION
	// Surface coloring for debugging
	void toggleSurfaceColoring(bool onOff);
	void setSurfaceColorRange(float min, float max);
	void tieSurfaceToImageSlider(bool yesno);
    void setMatchSurfaceStyle(M3DSurfaceStyle style) {
		colorSurfaceStyle = style;
		surfaceDone = false;
        boundaryCurvesDone = false;
	}

	void togglePlainMatchSurfaces(bool onOff) {
		outlineSurfaceTiles = !onOff;
		updateSurface();
	}
#endif

    // Set whether to display the standard axes
    void viewStdAxes(bool yesno) { stdAxesOn = yesno; }

    // Set whether to use motion cueing
    void viewRocking(bool yesno) { rockingOn = yesno; }
    bool rockingMode() { return rockingOn; }

    // Set whether to display the unit cube
    void viewUnitCube(bool yesno) { unitCubeOn = yesno; }

    // Set whether to display the image slices
    void viewXSlice(bool yesno) { xSliceOn = yesno; }
    void viewYSlice(bool yesno) { ySliceOn = yesno; }
    void viewZSlice(bool yesno) { zSliceOn = yesno; }
//    void viewArbSlice(bool yesno) { arbSliceOn = yesno; }

    // Set whether to display the landmarks
    void viewLandmarks(bool yesno) { landmarksOn = yesno; }

	void viewBNCutPlane(bool yesno) { bNCutPlaneOn = yesno; }
	void viewBBperpCutPlane(bool yesno) { bBperpCutPlaneOn = yesno; }
	void viewBperpNCutPlane(bool yesno) { bperpNCutPlaneOn = yesno; }
	void viewBperpY0CutPlane(bool yesno) { bperpY0CutPlaneOn = yesno; }
	void viewBperpY1CutPlane(bool yesno) { bperpY1CutPlaneOn = yesno; }
	void viewInvolutesCutPlane(bool yesno) { involutesCutPlaneOn = yesno; }

	bool viewingBNCutPlane() const { return bNCutPlaneOn; }
	bool viewingBBperpCutPlane() const { return bBperpCutPlaneOn; }
	bool viewingBperpNCutPlane() const { return bperpNCutPlaneOn; }
	bool viewingBperpY0CutPlane() const { return bperpY0CutPlaneOn; }
	bool viewingBperpY1CutPlane() const { return bperpY1CutPlaneOn; }
	bool viewingInvolutesCutPlane() const { return involutesCutPlaneOn; }

    // Tile set viewing controls
    void viewTileSet(bool yesno) { tileSetOn = yesno; }
    void setTileAlpha(float val) { tileAlpha = val; }
	void setTilesColor(float * color) {
		defaultTilesColor[0] = color[0];
		defaultTilesColor[1] = color[1];
		defaultTilesColor[2] = color[2];
		defaultTilesColorChanged = true;
		gMasterDispListIsValid = false;
        redraw();
	}
	void setModelColorsChanged() { tileColorsChanged = true; }
	void useDefaultColor(bool yesno) {
		useDfltTileSetColor = yesno;
		gMasterDispListIsValid = false;
        redraw();
	}
	void tileSetSolid(bool yesno) {
		renderTileSetSolid = yesno;
		gMasterDispListIsValid = false;
        redraw();
	}
	void setTileFileChanged() { tileSetDone = false; }

    // Auxiliary tile set viewing controls
    void viewAuxTileSet(bool yesno) {
		auxTileSetOn = yesno;
		auxTileSetDone = false;
	}
	void setAuxFileChanged() { auxTileSetDone = false; }

//    void setCutPlaneIntensity(float val) { cutPlaneIntensity = val; }

    // Set whether mouse click should select primitive, slice or figure
    void setSelectType(P3DSelectionType sType) { selectionType = sType; }
    P3DSelectionType getSelectType() const { return selectionType; }

    // Reset the trackball to zero rotation
    void resetTrackball();
	Trackball* getTrackball() { return &trackball; }

    // Set the trackball to a orthogonal view
    void setAxialView();
    void setSaggitalView();
    void setCoronalView();

    void rotateLeft90();
    void rotateRight90();
    void rotateHorizontal180();
    void rotateVertical180();

	void setCenterOfRotation(int which=0);

	// masks
#define CLIP_MASK_ALL     0x7E	// all possible masks
#define CLIP_MASK_ENABLED 0x01	// enables all other clips
#define CLIP_MASK_TILE    0x02	// clip tiles
#define CLIP_MASK_MODEL   0x04	// clip model surface
#define CLIP_MASK_AXES    0x08	// clip unit cube axes
#define CLIP_MASK_IMAGE   0x10	// clip image planes
#define CLIP_MASK_ILM     0x20	// clip image landmarks
#define CLIP_MASK_CAMERA  0x40	// clip ALL objects in scene
#define CLIP_MASK_VIEW    0x80	// same as camera but tracks view

	void enableClips(int mask,	// see masks, above
		int enable);

    // Handle mouse events
    int handle(int eventId);

    void draw();
	void checkToStopBlending();

    // Called when the object and surface need updating and redrawing
    void update() {
		surfaceDone = false;
		auxTileSetDone = false;
		blendingLevel = 0;
        gMasterDispListIsValid = false;
        objectDone = false;
        boundaryCurvesDone = false;
#ifdef BINARY
		distanceVectorsDone = false;
#endif
        redraw();
    }

    // Called when the object and standard surface need updating and redrawing
	// Does not re-render the implicit surface
    void updateSurface() {
		surfaceDone = false;
        gMasterDispListIsValid = false;
        objectDone = false;
        boundaryCurvesDone = false;
        redraw();
    }

    // Called when only the object needs updating and redrawing
    void updateObject() {
        gMasterDispListIsValid = false;
        objectDone = false;
        boundaryCurvesDone = false;
        redraw();
    }

    // Called when only the image plane slider was moved.  An extra call
    // is made after a delay to cause the plane to be interpolated.  (See
    // ImagePlanes.cpp).
    void updateImagePlanes(bool final) {
        gMasterDispListIsValid = false;
        imagePlanesDone = false;
        imageStable = final;
        boundaryCurvesDone = false;
        redraw();
    }

    // Make sure window is always square
    void resize(int x, int y, int w, int h) {
        if (w > h)
            w = h;
        else if (h > w)
            h = w;

//        cout << "View window: " << w << " x " << h << endl;
        Fl_Gl_Window::resize(x, y, w, h);
    }

    void show() {
        gMasterDispListIsValid = false;
        init = false;
        objectDone = false;
        boundaryCurvesDone = false;
        surfaceDone = false;
        if (control != NULL) {
            M3DObject * object = control->getObjectPtr();
            if (object != NULL)
                object->setModified(true);
        }
        Fl_Gl_Window::show();
        redraw();
    }

	// Cause a re-computation of the display list to happen after the current
	// redraw finishes.  This is used to cause the cut plane views to be
	// redrawn after a cut plane window has been recomputed.
    void signalRedraw() {
		must_redraw = true;
	}

    void startScreenDump() { dumpScreen = true; }
    void stopScreenDump() { dumpScreen = false; }
	void scheduleScreenDump() { screenDumpScheduled = 2; }
    void screenDump();

	// Landmarks
    void initializeEditLandmarks();
    void setLandmarkMode(bool yesNo) { landmarkMode = yesNo; }
    bool getLandmarkMode() { return landmarkMode; }
	void setLandmarksColor(const float color[3]);
	void setLandmarksNarrowWidth(int width);
	void setLandmarksWideWidth(int width);

#ifndef PRODUCTION_VERSION
	void setBlendingRegionDisplay(bool onOff);
#endif

	void getCenter(Vector3D & v);	   // Set v to object coordinates of view center
	void rotateAboutCenter(double angle);		// Rotate around viewport center
	void rotateAboutHorizontal(double angle);	// Rotate up/down in viewport
	void rotateAboutVertical(double angle);		// Rotate right/left in viewport
	void translateToCenter(Vector3D & v);		// Translate to v - current_position

	void printMatrices();

	void installBubblePopUp();
	void removeBubblePopUp();
	bool insideModelWindow();

	const char * getMouseMode() const {		// Debugging function
		return trackball.getMouseMode();
	}

#ifdef BINARY
	void setDistanceVectorList(DistanceVectorList * dvlist) {
		distanceVectorList = dvlist;
		distanceVectorsDone = false;
	}
#endif

	void setTileSetMobile(bool val) { tileSetMobile = val; }

private:

	int mouseTipX, mouseTipY;

    // Draw a box when the mouse is dragged
    void drawMouseBox(int x1, int y1, int x2, int y2);

    // Manually select primitives
    void mouseSelect(bool chooseOrigin);
	void mouseSelectImage();

	// Determine landmark boundary intersections
	void mouseSelectTile();
	bool intersectQuad(Point & pt0, Point & pt1, Point & pt2, Point & pt3,
		Point & a, Point & b);

    void drawObject(bool showConstraintArrows = false);
    void drawFigureNames();
    void drawImageSlices();
	void drawBoundaryCurves(bool changed);
	void drawAllCutPlanes();
#ifdef BINARY
	void drawDistanceVectors();
#endif
	void updateAllCutPlanes();
	void drawCutPlane(cut_plane_t planeType);
	void rockingMotion();
	void drawConstraintArrows(M3DObject * object);

    // Flag to determine whether GL stuff is initialized
    bool init;

    // Flags indicating what to display
    bool primitiveVectorsOn;
    bool medialSheetOn;
	bool atomsOn;
    bool stdAxesOn;
	bool rockingOn;
    bool unitCubeOn;
    bool xSliceOn;
    bool ySliceOn;
    bool zSliceOn;
//    bool arbSliceOn;

    bool tileSetOn;
	bool tileSetDone;
	TileSetRenderer * tileRenderer;
	int currentFigureCount;
	const float ** tileColors;
	bool useDfltTileSetColor;
	bool tileColorsChanged;
	bool renderTileSetSolid;
    bool auxTileSetOn;
    bool auxTileSetDone;

	bool landmarksOn;

	// Flags for the cut planes display
	bool bNCutPlaneOn;
	bool bBperpCutPlaneOn;
	bool bperpNCutPlaneOn;
	bool bperpY0CutPlaneOn;
	bool bperpY1CutPlaneOn;
	bool involutesCutPlaneOn;
	bool texture_done;
//	float cutPlaneIntensity;

	// Normal rendering flags
	int surfaceLevel;
    M3DSurfaceStyle surfaceStyle;
	int surfaceLineWidth;
	bool blendSurface;
	bool surfaceBlendedLevelChanged;
	bool drawBoundary;
	int cutPlaneBoundaryWidth;
	bool surfaceDone;
	bool constraintArrows;
	bool labelFigures;

	M3DObjectSurfaceRenderer surfaceRenderer;
	bool twoLights;

#ifndef PRODUCTION_VERSION
	// Colored surface rendering for debugging
	bool colorSurface;
	bool outlineSurfaceTiles;
	float surfaceColorMin, surfaceColorMax;
	bool surfaceTiedToImageSlider;
    M3DSurfaceStyle colorSurfaceStyle;
	void getSurfaceIntensityRange(float & min, float & max);
#endif

	// Objects for drawing figural boundaries in image plane windows
	vector<CutPlanes> xBoundary;
	vector<CutPlanes> yBoundary;
	vector<CutPlanes> zBoundary;
//	bool xBoundaryDrawingMode;
//	bool yBoundaryDrawingMode;
//	bool zBoundaryDrawingMode;
	bool surfaceLevelChanged;

	// Interpenetration avoidance rendering flags
	int partialSurfaceLevel;
    M3DSurfaceStyle partialSurfaceStyle;

    // Mouse Handling data
    bool mouseDown;
    int  lastMouseX, lastMouseY;
    int  currMouseX, currMouseY;

	// Centering Mode
	int centerMode;		// last thing centered upon

    // Flag to indicate the mouse is actively editing the model
	bool modelIsBeingEdited;
	bool modelWasEdited;

	char idTooltip[MAX_TOOLTIP_SIZE];
    P3DSelectionType selectionType;

	// Rocking motion parameters
	double currentAngle;
	double deflectionAngle;
	double incrementalAngle;
	bool direction;

    Trackball trackball;

    // Saves the scale value of the view at the last draw
    double lastScale;

    // Alpha value for the tile set and the image
    float tileAlpha;
//    float imageAlpha;

    P3DControl * control;
    P3DUserInterfaceCallback * callback;

    // Set to false if the object view properties were modified since last draw
    bool objectDone;

    // Set to false if image planes have changed since last draw
    bool imagePlanesDone;

	// Indicates that the image has stopped changing (e.g. by slider action)
	bool imageStable;

    // Set to false if boundary curves have changed since last draw
    bool boundaryCurvesDone;

#ifdef BINARY
	// Set to false if the distance vectors have changed since last draw
	bool distanceVectorsDone;
	DistanceVectorList * distanceVectorList;	// Diagnostic distance vector list
#endif

	// set true whenever tiles are allowed to move (which rebuilds tile display list)
	bool tileSetMobile;

    // OpenGL texture objects
    unsigned int textureNames[4];
	unsigned int cutPlaneTextureNames[6];

    float transformMatrix[4][4];

    // Flag to indicate we are in the middle of the draw() function
    bool isDrawing;

	redraw_t damaged;
	bool must_redraw;

    // This is set for the orthogonal views
    ViewType viewType;

	// Model window background color
	GLclampf backgroundColor[3];
	// Tileset color
	GLclampf defaultTilesColor[3];
	bool defaultTilesColorChanged;

	// Variables used when picking/editing landmarks
    bool landmarkMode;
	M3DObjectSurfaceRenderer * landmarkSurface;
	double * distThreshold;
	int numLandmarkHits;		// Number of surface or image hits
	Vector3D * landmarks;		// Surface or image hits
	float defaultLandmarkColor[3];
	float defaultLandmarkWidth;
	float markedLandmarkWidth;

	// Blending renderer
	M3DBlendedRenderer blendedSurfaceRenderer;
	int blendingLevel;
	bool blendedRendererRunning;
	bool interruptOccurred;

	void drawSurfaces();
	void drawAuxSurfaces();
	void setBoundariesLevel();
	void drawObjectLandmarks();
    void drawImageLandmarks();

	GLenum tileType;	// Interpolated normals or flat

    bool dumpScreen;
	int screenDumpScheduled;

#ifndef PRODUCTION_VERSION
	bool displayBlendRegion;
#endif

    bool in_front(int testId, int primitiveId);

    void initTextures();
};

#endif

