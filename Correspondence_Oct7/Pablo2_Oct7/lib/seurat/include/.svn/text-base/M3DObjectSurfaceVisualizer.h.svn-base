#ifndef M3DOBJECT_SURFACE_VISUALIZER_H
#define M3DOBJECT_SURFACE_VISUALIZER_H




/*  Class M3DObjectSurfaceVisualizer

	This class provides a capability for displaying various values computed
	in Match.cpp on the object surface.  Capability is also provided in
	P3DUserInterface for loading color value files for display on the object
	surface.

	Pablo contains a single, global instance of this class.  Defined by this
	include file, it is named "objectSurfaceVisualizer".  This is the only
	interface to this object, i.e., it is not accessible through P3DControl.
	The reason for this is both to simplify the code, and because it is
	intended to be a local debugging tool, not often used.  It is not
	compiled in Release configurations of Pablo.

	To use the capability provided here, the initializeSurfaceMap() function
	must be called once per figure to be displayed, and the surface colors
	(actually grey values) must be provided by a series of calls to
	setSurfaceColor().  These are the only two member functions of this class
	ordinarily to be used.  Their usage is described below.  Three surface
	visualizations are built into Match, and can be activated by uncommenting
	the definition of SURFACE_VISUALIZER in globalBuildParms.h.


	bool initializeSurfaceMap(int serialNumber, char * name,
			M3DObject * object, int figureId, int surfaceLevel, int size = 0,
			ThallCode::Pointlist_server2 ** p_lists = NULL,
			bool rescale = true)

		This function acts much as a constructor does in other classes.  It
		assures one-time allocation of the necessary class variables.  Thus,
		it may be used within a loop; successive calls are handled gracefully.
		The returned value will be false, only when the object is NULL.  If
		the number of figures in the model or the surface level change, the
		data in this class will be discarded automatically, so that successive
		calls to initializeSurfaceMap() will cause new data to be stored for
		the renderings.

		int serialNumber
		  Multiple visualizations may be defined, and each one must have a
		  unique, serial number.  Thus, with only one visualization, the
		  value 0 would be used.

		char * name
		  This is a short name for the visualization to be displayed in the
		  user's interface.  It may be NULL or an empty string, if no name
		  is desired.

		M3DObject * object
		  The object representation to be rendered.

		int figureId
		  The number of the particular figure in the object being initialized.

		int surfaceLevel
		  The surface level-of-detail to be used in the rendering.  In Match,
		  this should probably always equal Match::surface_level_for_match.
		  Note that this will not always correspond to the level selected by
		  the user on the Display Window's Surface Tab.  All figures are
		  expected to have the same level.

		int size
		  This indicates the number of surface points to be supplied.  It may
		  be zero, in which case the number will be calculated based on the
		  indicated surface level.  Normally, it should be assigned the
		  number of surface points corresponding to the point lists used.
		  If it is not provided, the value used can be obtained by calling
		  boundarySize().

		ThallCode::Pointlist_server2 ** p_lists
		  The set of point lists to be used for the rendering.  Typically, the
		  one in Match will be used.  However, if this is NULL, then point
		  lists will be allocated for the rendering.

		  There is one possible side effect of providing a point list array.
		  This class will call UpdateSubdivPointCloud() and may call
		  ComputeSubdivBoundaryTiles() on the supplied point lists (See
		  ThallCode::Pointlist_server2).  However, this will normally not be a
		  problem and should not affect computations done in Match.

		bool rescale
		  Most data computed in match, are not image intensities but things
		  such as objective function values, for which the range is not
		  predictable.  To support the visualization of such quantities, the
		  default action (rescale = true) is to map the intensities to [0, 1]
		  before rendering.  The application of slider settings occurs after
		  this initial mapping.

		  However, if rescale is false, the input data are assumed to be image
		  intensities.  In this case, they will not be scaled before display,
		  other than according to slider values supplied during the rendering.
		  Image intensities rendered when a value of true for rescale will
		  produce an incorrect visualization.

		  To get a correct rendering of image intensities, it often will be
		  necessary that the slider values on the Image Tab of the Display
		  Window be in their default position before the data is computed
		  for visualization.  This is because the usual means of getting image
		  intensities is to call Image3D::getWindowedInterpolatedVoxelValue(),
		  which provides values based on the slider values.  Instead, if
		  getInterpolatedVoxelValue() is used, the image intensities will be
		  correct regardless of the initial slider setting.  Note that in
		  either case, division by MAX_GREY_VALUE is necessary before the
		  intensities are provided to setSurfaceColor().


	void setSurfaceColor(int figureId, int surfaceColorIndex, float color,
			int serialNumber = 0)

		This function should be used to supply the surface data values after
		initialization.  It will generally be called repeatedly inside a loop
		over all surface points.  If the surface level in Match increases
		after this function is used, it will be necessary to reinitialize
		before calling this function with the new data.  If this is not done,
		the proper surface will not be displayed.

		int figureId
		  The number of the particular figure in the object being initialized.

		int surfaceColorIndex
		  The number of the particular point (Pointlist index) on the surface
		  of the object for which the data value is being supplied.

		float color
		  The surface data value corresponding to the values of the previous
		  two variables.

		int serialNumber
		  The serial number for a particular rendeing corresponding to the call
		  to initializeSurfaceMap().


	If calls to these two functions are installed permanently in Pablo, they
	should be surrounded by:
		#ifndef PRODUCTION_VERSION
		#endif
	There will be no obvious effect if this is not done, but the program will
	run slightly slower if it is not, since the visualizer is not compiled
	into the production version of Pablo.


	Examples:

	  Here is a sample of the typical use of this function.  In this case, the
	  objective function computed in Match::computeMainFigureMatch() is
	  displayed.

		  // These lines added near the start of computeMainFigureMatch()
		  objectSurfaceVisualizer->initializeSurfaceMap(0, "Objective Fn",
			  targetObject, figureId, surface_level_for_match, numPoints,
			  pList);

		  // This line is added at the end of the second stage for-loop in
		  //   computeMainFigureMatch().  The loop index is "j".
		  objectSurfaceVisualizer->setSurfaceColor(figureId, j,
			  val1 * ((val2 - mu) / rms) * maskElement->volelem);

	  In the next example, the image intensities are to be displayed.  Code
	  insertions are in the same locations in Match::computeMainFigureMatch()
	  as those above.  As many things as possible have been changed for this
	  example relative to the previous one.  Note that this code is intended
	  to be added to the example above, since it specifies the rendering
	  serial number to be 1.

		  // Because the surface level is explicitly 3, the highest level (3)
		  // must be selected on the Display Window's Surface Tab before the
		  // optimizer is initialized.
		  objectSurfaceVisualizer->initializeSurfaceMap(1, "Intensities",
			  targetObject, figureId, 3, 0, NULL, false);

		  // Variable val2 contains the image intensity at point j
		  objectSurfaceVisualizer->setSurfaceColor(figureId, j, val2, 1);


	Known problem:

		The surface does not move properly with the object.  This is because this
		class does not use the point lists in class M3DObjectSurfaceRenderer used
		for normal surface rendering.
*/


#ifdef PRODUCTION_VERSION
#define SURFACE_VISUALIZER
#endif
#ifdef AE2_BUILD
#undef SURFACE_VISUALIZER
#endif

class M3DObject;

class M3DObjectSurfaceVisualizer
{

#ifdef SURFACE_VISUALIZER

	public:

		// Initialize for generation of diagnostic surface colors (greyscale
		// values) for debugging
		bool initializeSurfaceMap(int serialNumber, char * name,
			M3DObject * object, int figureId, int surfaceLevel, int size = 0,
			ThallCode::Pointlist_server2 ** p_lists = NULL,
			bool rescale = true);

		// Get an unused serial number.  This is useful so that the GUI can load
		// a set of surface color values from a file into an empty slot.
		int getFreeSerialNumber();

		// Specify diagnostic surface colors (greyscale values at surface points)
		// for debugging.  The values provided will be treated as indexes into
		// a colormap.
		void setSurfaceColor(int figureId, int surfaceColorIndex, float color,
			int serialNumber = 0);

		// Returns the point list size of the indicated figure
		int boundarySize(int figureId);

		// Evaluate the object to update the surface.  This gives a way to
		// change the surface level used in rendering all visualizations of
		// any figure.  It forces a recomputing of the Pointlists for this figure.
		// If the Pointlist was provided in initializeSurfaceMap(), then there
		// may be unintended side-effects outside this class.  It is also
		// possible to change the contents of a figure or model and then call
		// this function, but such behavior is risky and discouraged.  The
		// difficulty with this function is that if the level is changed, then
		// the current visualization will not longer match the new subdivision.
		// So, this should only be used when that can also be changed to the
		// new surface level number.
		void cacheSurface(M3DObject * object, int figureId, int surfaceLevel) {
			cacheSurface(object, figureId, surfaceLevel, true);
		}

		void clear();	// Discard all surface color value sets

	private:

		struct SurfaceColorMappings		// Storage for one visualization (serial #)
		{
			int numFigures;		// Number of figures to be rendered
			int surfaceLevel;	   // Surface level for rendering the visualization
			bool pListAllocated;   // False if user provided the pLists; true if generated
			float ** colors;	// Color list per figure; one entry per point in each list
			int * numColors;	// Number of entries in 'colors'; one per figure
			int * allocation;	// Length of point list; one per figure
			bool * ranged;		// True, unless a rescaling is needed; one per figure
			ThallCode::Pointlist_server2 ** pList;	// Computed or provided points
			char * name;		// Visualization name
		};

		static SurfaceColorMappings ** surfaceColorsMap;	// Array; one entry per serial #
		static int currentMap;
		static bool excessColors;
		static bool unequalSizes;
		static int w;

		void deleteSurfaceColorsMap(int map);
		bool initializeRendering(int figureId, int level, float min, float max);

		// Make room for a cached copy of the Pointlists (all figures)
		bool allocatePList(int serialNumber, M3DObject *object, int surfaceLevel,
			ThallCode::Pointlist_server2 ** p_lists = NULL);

		// Evaluate the object to update cached surface.
		// Only internal routines should use computeSurface of false.
		void cacheSurface(M3DObject *object, int figureId, int surfaceLevel,
			bool computeSurface);

#endif	/* SURFACE_VISUALIZER */

	public:

		// The remaining functions are called from within Pablo and may
		// generally be ignored by users
		M3DObjectSurfaceVisualizer();
		~M3DObjectSurfaceVisualizer();

#ifdef SURFACE_VISUALIZER
		// Returns the name of the current visualization
		const char * visualizationName() const;

		// Specify which visualization is to be rendered.  The resulting
		// serial number is returned.  If the number provided is not returned,
		// then the function failed.
		int selectBySerialNumber(int serialNumber);

		// Get the color in the range [0, 1] associated with a boundary point.  If the
		// maps aren't set up properly, this returns the default value, 0.5.
		double getColorValue(int serialNumber, int figureId, int index, double defValue);
		double getColorValue(int figureId, int index, double defValue = 0.5) {
			return getColorValue(currentMap, figureId, index, defValue);
		}

#endif	/* SURFACE_VISUALIZER */

		// Set the width of lines and dots used for rendering
		void lineWidth(int width);

		void renderColoredSubdivPoints(int figureId, int level,
			float min, float max);
		void renderColoredSubdivWires(int figureId, int level,
			float min, float max);
		void renderColoredSubdivTiles(int figureId, int level,
			float min, float max, bool outlineTiles, double displacement);
};


extern M3DObjectSurfaceVisualizer * objectSurfaceVisualizer;


// This class is used for automatic initialization of objectSurfaceVisualizer.
// It should never be used directly. 
class InitM3DObjectSurfaceVisualizer {
	public:
		InitM3DObjectSurfaceVisualizer();
		~InitM3DObjectSurfaceVisualizer();
	private:
		static int counter;
};

static InitM3DObjectSurfaceVisualizer visualizerInit;


#endif

