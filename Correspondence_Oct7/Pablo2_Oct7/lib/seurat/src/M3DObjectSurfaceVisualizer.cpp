#include "globalBuildParms.h"
#include "utility.h"
#include "M3DQuadFigure.h"
#include "M3DObject.h"

#define D_XFERLIST
#define D_POINTLIST_SERVER2
#include "Shapedepend.h"

#include "M3DObjectSurfaceVisualizer.h"


/*
	This class adds surface visualization capabilities to Match and to
	P3DUserInterface.

	It is only available when SURFACE_VISUALIZER is not defined, which is
	true for many build configurations.  In particular, it is never
	available when PRODUCTION_VERSION is defined.  See the class header
	file, M3DObjectSurfaceVisualizer.h.  
*/


//#define DEBUG

// If this is changed, corresponding changes must also be made in P3DUserCallback
#define MAX_MAP_COUNT	256	/* The length of vector surfaceColorsMap */

M3DObjectSurfaceVisualizer * objectSurfaceVisualizer = NULL;

#ifdef SURFACE_VISUALIZER
M3DObjectSurfaceVisualizer::SurfaceColorMappings
	** M3DObjectSurfaceVisualizer::surfaceColorsMap = NULL;
int  M3DObjectSurfaceVisualizer::currentMap = 0;
bool M3DObjectSurfaceVisualizer::excessColors = false;
bool M3DObjectSurfaceVisualizer::unequalSizes = false;
int M3DObjectSurfaceVisualizer::w = 1;   // Line/dot width
#endif	/* SURFACE_VISUALIZER */

int InitM3DObjectSurfaceVisualizer::counter = -1;

using namespace std;
using namespace ThallCode;


InitM3DObjectSurfaceVisualizer::InitM3DObjectSurfaceVisualizer()
{
	if (++counter == 0)
		objectSurfaceVisualizer = new M3DObjectSurfaceVisualizer;
}

InitM3DObjectSurfaceVisualizer::~InitM3DObjectSurfaceVisualizer()
{
	if (--counter < 0)
		delete objectSurfaceVisualizer;
}

M3DObjectSurfaceVisualizer::M3DObjectSurfaceVisualizer()
{
#ifdef SURFACE_VISUALIZER
	if (surfaceColorsMap == NULL) {
		surfaceColorsMap = new SurfaceColorMappings*[MAX_MAP_COUNT];
		for (int i = 0; i < MAX_MAP_COUNT; i++)
			surfaceColorsMap[i] = NULL;
	}
#endif	/* SURFACE_VISUALIZER */
}

M3DObjectSurfaceVisualizer::~M3DObjectSurfaceVisualizer()
{
#ifdef SURFACE_VISUALIZER
	for (int i = 0; i < MAX_MAP_COUNT; i++)
		deleteSurfaceColorsMap(i);
	delete [] surfaceColorsMap;
#endif	/* SURFACE_VISUALIZER */
}

#ifdef SURFACE_VISUALIZER

void M3DObjectSurfaceVisualizer::clear()
{
	for (int i = 0; i < MAX_MAP_COUNT; i++) {
		deleteSurfaceColorsMap(i);
		surfaceColorsMap[i] = NULL;
	}
	currentMap = 0;
}

void M3DObjectSurfaceVisualizer::deleteSurfaceColorsMap(int map)
{
    int i;

	if (surfaceColorsMap[map] != NULL) {
		for (i = 0; i < surfaceColorsMap[map]->numFigures; i++)
			if (surfaceColorsMap[map]->colors[i] != NULL)
				delete [] surfaceColorsMap[map]->colors[i];
		delete [] surfaceColorsMap[map]->colors;
		delete [] surfaceColorsMap[map]->numColors;
		delete [] surfaceColorsMap[map]->allocation;
		delete [] surfaceColorsMap[map]->ranged;
		delete [] surfaceColorsMap[map]->name;
		if (surfaceColorsMap[map]->pListAllocated)
			delete [] surfaceColorsMap[map]->pList;

		delete surfaceColorsMap[map];
	}
}

int M3DObjectSurfaceVisualizer::getFreeSerialNumber() {
	for (int i = 0; i < MAX_MAP_COUNT; i++) {
		if (surfaceColorsMap[i] == NULL) {
			return i;
		}
	}
	return -1;
}

bool M3DObjectSurfaceVisualizer::allocatePList(int serialNumber,
    M3DObject * object, int surfaceLevel, Pointlist_server2 ** p_lists)
{
	int fig;

	if (object == NULL) {
		cout << "An M3DObject is required to initialize surface coloring\n";
		return false;
	}
	surfaceColorsMap[serialNumber] = new SurfaceColorMappings;
	surfaceColorsMap[serialNumber]->numFigures = object->getFigureCount();
	surfaceColorsMap[serialNumber]->colors = new float *[surfaceColorsMap[serialNumber]->numFigures];
	surfaceColorsMap[serialNumber]->numColors = new int[surfaceColorsMap[serialNumber]->numFigures];
	surfaceColorsMap[serialNumber]->allocation = new int[surfaceColorsMap[serialNumber]->numFigures];
	surfaceColorsMap[serialNumber]->ranged = new bool[surfaceColorsMap[serialNumber]->numFigures];
	surfaceColorsMap[serialNumber]->name = NULL;

	if (p_lists == NULL) {
		surfaceColorsMap[serialNumber]->pList =
			new Pointlist_server2 *[surfaceColorsMap[serialNumber]->numFigures];
		for (fig = 0; fig < surfaceColorsMap[serialNumber]->numFigures; fig++)
			surfaceColorsMap[serialNumber]->pList[fig] = NULL;
		surfaceColorsMap[serialNumber]->pListAllocated = true;
	}
	else
		surfaceColorsMap[serialNumber]->pListAllocated = false;

	for (fig = 0; fig < surfaceColorsMap[serialNumber]->numFigures; fig++) {
		surfaceColorsMap[serialNumber]->colors[fig] = NULL;
		surfaceColorsMap[serialNumber]->numColors[fig] = 0;
		surfaceColorsMap[serialNumber]->allocation[fig] = 0;
		surfaceColorsMap[serialNumber]->ranged[fig] = true;
	}
	surfaceColorsMap[serialNumber]->surfaceLevel = surfaceLevel;

	return true;
}

// Evaluate the object to update cached surface
void M3DObjectSurfaceVisualizer::cacheSurface(M3DObject * object, int figureId,
    int surfaceLevel, bool computeSurface)
{
	// Update p_lists of existing surfaceColorsMap's for all visualizations
	for (int serialNumber = 0; serialNumber < MAX_MAP_COUNT; serialNumber++) {
		if (surfaceColorsMap[serialNumber]) {

			// Determine the number of surface points for which colors must be stored
			M3DFigure * figure = object->getFigurePtr(figureId);
			Xferlist * xferList = convertM3DtoXfer(figure);

			if (computeSurface) {
				surfaceColorsMap[serialNumber]->pList[figureId]->init(xferList);
				surfaceColorsMap[serialNumber]->pList[figureId]->
					ComputeSubdivPointCloud(surfaceLevel);
			}
			surfaceColorsMap[serialNumber]->pList[figureId]->UpdateSubdivPointCloud(
				surfaceLevel, xferList);	// AGG: Needed?

			delete [] (xferList->atomlist);
			delete xferList;
		}
	}
}

bool M3DObjectSurfaceVisualizer::initializeSurfaceMap(int serialNumber, char * name,
	M3DObject * object, int figureId, int surfaceLevel, int size,
	Pointlist_server2 ** p_lists, bool rescale)
{
	int i;

#ifdef DEBUG
		cout << "M3DObjectSurfaceVisualizer::initializeSurfaceMap()\n";
#endif

	if (object == NULL) {
		cout << "Error: no object provided when initializing a surface visualization"
			<< endl;
		return false;
	}
	if (figureId >= object->getFigureCount()) {
		cout << "Error: invalid figure number when initializing a surface visualization"
			<< endl;
		return false;
	}

	if (surfaceColorsMap[serialNumber] != NULL) {

		// If the surface level has changed or figures were added, discard all surface data
		if (surfaceColorsMap[serialNumber] != NULL
			&& (surfaceLevel != surfaceColorsMap[serialNumber]->surfaceLevel
			|| object->getFigureCount() > surfaceColorsMap[serialNumber]->numFigures)) {
				deleteSurfaceColorsMap(serialNumber);
				surfaceColorsMap[serialNumber] = NULL;
		}
	}

	if (surfaceColorsMap[serialNumber] == NULL) {
		// Initialize the map for a particular serial number

		if (! allocatePList(serialNumber, object, surfaceLevel, p_lists))
			return false;

		if (name != NULL && name [0] != '\0') {
			int len = strlen(name);
			surfaceColorsMap[serialNumber]->name = new char[1 + len];
			for (i = 0; i < len; i++)
				surfaceColorsMap[serialNumber]->name[i] = name[i];
			surfaceColorsMap[serialNumber]->name[len] = '\0';
		}

		excessColors = false;
		unequalSizes = false;
	}

	// Determine which point list to use
	bool computeSurface = false;
	if (! surfaceColorsMap[serialNumber]->pListAllocated)
		surfaceColorsMap[serialNumber]->pList = p_lists;
	else {
		if (surfaceColorsMap[serialNumber]->pList[figureId] == NULL) {
			surfaceColorsMap[serialNumber]->pList[figureId] = new Pointlist_server2;
			computeSurface = true;
		}
	}

	if (size == 0 || computeSurface) {
		cacheSurface(object, figureId, surfaceLevel, computeSurface);
		int s =  boundarySize(figureId);  // Determine the size using the current map
		if (size != 0 && s != size) {
			if (! unequalSizes) {
				cout << "Warning: wrong number of surface colors provided;"
					<< " surface coloring may be incorrect\n";
				unequalSizes = true;
			}
		}
		else
			size = s;
#ifdef DEBUG
		cout << "Pointlist size = " << size << '\n';
#endif
	}

	// Verify that enough space for colors is present in the array
	if (surfaceColorsMap[serialNumber]->colors[figureId] != NULL) {
		if (surfaceColorsMap[serialNumber]->allocation[figureId] < size) {
			delete [] surfaceColorsMap[serialNumber]->colors[figureId];
			surfaceColorsMap[serialNumber]->colors[figureId] = NULL;
			surfaceColorsMap[serialNumber]->allocation[figureId] = 0;
		}
	}

	// If there is no array for colors, create one
	if (surfaceColorsMap[serialNumber]->colors[figureId] == NULL) {
#ifdef DEBUG
		cout << "Allocating " << size << " surface colors\n";
#endif
		surfaceColorsMap[serialNumber]->colors[figureId] = new float[size];
		surfaceColorsMap[serialNumber]->allocation[figureId] = size;
	}

	surfaceColorsMap[serialNumber]->numColors[figureId] = 0;
	surfaceColorsMap[serialNumber]->ranged[figureId] = ! rescale;

	return true;
}

// Note: the serial number is not needed here, since the number of boundary
// points is the same for each visualization.
int M3DObjectSurfaceVisualizer::boundarySize(int figureId)
{
		Bpoint *plist;
		int s;

		if (! surfaceColorsMap[currentMap] ||
			! surfaceColorsMap[currentMap]->pList[figureId])
				return 0;

		surfaceColorsMap[currentMap]->pList[figureId]->subdivboundaryinfo(&s, &plist);
		return s;
}

void M3DObjectSurfaceVisualizer::setSurfaceColor(int figureId, int surfaceColorIndex,
												 float color, int serialNumber)
{
	// For speed this function doesn't bother to verify that the serialNumber,
	// figureId, and surfaceColorIndex are valid
	if (surfaceColorIndex >= surfaceColorsMap[serialNumber]->allocation[figureId]) {
		if (! excessColors) {
			cout << "Error: too many surface colors provided;"
				<< " surface coloring will be incorrect\n";
			excessColors = true;
		}
		return;
	}
	surfaceColorsMap[serialNumber]->colors[figureId][surfaceColorIndex++] = color;
	surfaceColorsMap[serialNumber]->numColors[figureId] = surfaceColorIndex;
}

double M3DObjectSurfaceVisualizer::getColorValue(int serialNumber, int figureId,
	int index, double defValue = 0.5)
{
	if (! surfaceColorsMap[serialNumber])
		return defValue;
	if (index >= surfaceColorsMap[serialNumber]->allocation[figureId])
		return defValue;
	return surfaceColorsMap[serialNumber]->colors[figureId][index];
}

const char * M3DObjectSurfaceVisualizer::visualizationName() const
{
		if (! surfaceColorsMap[currentMap])
			return NULL;

		return surfaceColorsMap[currentMap]->name;
}

int M3DObjectSurfaceVisualizer::selectBySerialNumber(int serialNumber)
{
	if (serialNumber < 0 || serialNumber >= MAX_MAP_COUNT) {
		cout << "Invalid visualization map number specified\n";
		return currentMap;
	}

	if (surfaceColorsMap && surfaceColorsMap[serialNumber])
		currentMap = serialNumber;
	return currentMap;
}

bool M3DObjectSurfaceVisualizer::initializeRendering(int figureId, int level,
	float min, float max)
{
	if (surfaceColorsMap[currentMap] == NULL) {
			std::cout << "Error: no visualizations have been initialized\n";
			return false;
	}

	if (figureId >= surfaceColorsMap[currentMap]->numFigures ||
		surfaceColorsMap[currentMap]->colors[figureId] == NULL) {
			std::cout << "Error: surface colors not available for figure "
				<< figureId << '\n';
			return false;
	}

	if (level != surfaceColorsMap[currentMap]->surfaceLevel) {
			std::cout << "Error: surface level changed from "
				<< surfaceColorsMap[currentMap]->surfaceLevel << " to "
				<< level << "; try reinitializing the optimizer\n";
			return false;
	}

	for (int nfig = 0; nfig < surfaceColorsMap[currentMap]->numFigures; nfig++)
	{
		if (! surfaceColorsMap[currentMap]->ranged[nfig]) {
			calc_intensity_range(surfaceColorsMap[currentMap]->colors[nfig],
				surfaceColorsMap[currentMap]->numColors[nfig], min, max);
#ifdef DEBUG
			cout << "Surface intensity range rescaled to [" << min << ", "
				<< max << "]\n";
#endif
			float scale = 1.0/(max - min);	// Map to [0, 1]
			for (int i = 0; i < surfaceColorsMap[currentMap]->numColors[nfig]; i++)
				surfaceColorsMap[currentMap]->colors[nfig][i] =
					(surfaceColorsMap[currentMap]->colors[nfig][i] - min)*scale;
			surfaceColorsMap[currentMap]->ranged[nfig] = true;
		}
	}

	return true;
}

#endif	/* SURFACE_VISUALIZER */

void M3DObjectSurfaceVisualizer::renderColoredSubdivPoints(int figureId, int level,		// AGG: This renders incorrectly.  Why?
	float min, float max)
{
#ifdef SURFACE_VISUALIZER
	if (! initializeRendering(figureId, level, min, max)) {
		if (surfaceColorsMap[currentMap] && surfaceColorsMap[currentMap]->pList[figureId])
			surfaceColorsMap[currentMap]->pList[figureId]->glRenderSubdivPoints();
		return;
	}
	surfaceColorsMap[currentMap]->pList[figureId]->lineWidth(w);
	surfaceColorsMap[currentMap]->pList[figureId]->glRenderColoredSubdivPoints(
		surfaceColorsMap[currentMap]->colors[figureId],
		surfaceColorsMap[currentMap]->numColors[figureId], min, max);
#ifdef DEBUG
		cout << "Visualization surface rendered with intensity range of [" << min << ", "
			<< max << "]\n";
#endif
#endif	/* SURFACE_VISUALIZER */
}

void M3DObjectSurfaceVisualizer::renderColoredSubdivWires(int figureId, int level,
	float min, float max)
{
#ifdef SURFACE_VISUALIZER
	if (! initializeRendering(figureId, level, min, max)) {
		if (surfaceColorsMap[currentMap] && surfaceColorsMap[currentMap]->pList[figureId])
			surfaceColorsMap[currentMap]->pList[figureId]->glRenderSubdivWires();
		return;
	}
	surfaceColorsMap[currentMap]->pList[figureId]->ComputeSubdivBoundaryTiles(
		surfaceColorsMap[currentMap]->surfaceLevel);
	surfaceColorsMap[currentMap]->pList[figureId]->lineWidth(w);
	surfaceColorsMap[currentMap]->pList[figureId]->glRenderColoredSubdivWires(
		surfaceColorsMap[currentMap]->colors[figureId],
		surfaceColorsMap[currentMap]->numColors[figureId], min, max);
#ifdef DEBUG
		cout << "Visualization surface rendered with intensity range of [" << min << ", "
			<< max << "]\n";
#endif
#endif	/* SURFACE_VISUALIZER */
}

void M3DObjectSurfaceVisualizer::renderColoredSubdivTiles(int figureId, int level,
	float min, float max, bool outlineTiles, double displacement)
{
#ifdef SURFACE_VISUALIZER
	if (! initializeRendering(figureId, level, min, max)) {
		if (surfaceColorsMap[currentMap] && surfaceColorsMap[currentMap]->pList[figureId])
			surfaceColorsMap[currentMap]->pList[figureId]->glRenderSubdivTiles(displacement);
		return;
	}
	surfaceColorsMap[currentMap]->pList[figureId]->ComputeSubdivBoundaryTiles(
		surfaceColorsMap[currentMap]->surfaceLevel);
	surfaceColorsMap[currentMap]->pList[figureId]->glRenderColoredSubdivTiles(
		surfaceColorsMap[currentMap]->colors[figureId],
		surfaceColorsMap[currentMap]->numColors[figureId], min, max,
		outlineTiles, displacement);
#ifdef DEBUG
		cout << "Visualization surface rendered with intensity range of [" << min << ", "
			<< max << "]\n";
#endif
#endif	/* SURFACE_VISUALIZER */
}

void M3DObjectSurfaceVisualizer::lineWidth(int width)
{
#ifdef SURFACE_VISUALIZER
	w = width;
#endif	/* SURFACE_VISUALIZER */
}

