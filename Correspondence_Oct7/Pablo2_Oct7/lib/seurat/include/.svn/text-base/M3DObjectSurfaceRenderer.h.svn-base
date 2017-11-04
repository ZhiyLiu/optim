#ifndef M3D_OBJECT_SURFACE_RENDERER_H
#define M3D_OBJECT_SURFACE_RENDERER_H


/*	Class M3DObjectSurfaceRenderer.

	This class uses class Pointlist_server2 to produce OpenGL renderings of
	MREP models.  Blending is not supported; each figure is separately
	rendered.
*/

enum M3DSurfaceStyle
{
    M3D_NONE = 0,
    M3D_POINTCLOUD,
    M3D_WIREFRAME,
    M3D_SOLID,
	// The following are used for object surface visualizations
    M3D_COLORED_POINTCLOUD,
    M3D_COLORED_WIREFRAME,
    M3D_COLORED_SOLID
};

class M3DObject;
class ThallCode::Pointlist_server2;
class SelectedPartialFigures;

class M3DObjectSurfaceRenderer
{
	friend class M3DBlendedRenderer;

public:

	// The constructor accepts a pointer to the multifigural object to
	// be rendered.  The surface level number must be specified in the
	// range [1, 4].
	M3DObjectSurfaceRenderer(M3DObject * obj = NULL, int lvl = 1);
	virtual ~M3DObjectSurfaceRenderer();

    ThallCode::Diatomgrid * convertFigure(M3DFigure * fig);

	void setObject(M3DObject * obj);
    void setSurfaceLevel(int level) { surfaceLevel = level; }
    void setSurfaceStyle(M3DSurfaceStyle style) { surfaceStyle = style; }
    void setPartialSurfaceStyle(M3DSurfaceStyle style) {
		partialSurfaceStyle = style;
	}

	// Set threshold range for object surface visualizations
	void surfaceColorRange(float min, float max) {
		surfaceColorsMin = min;
		surfaceColorsMax = max;
	}

    void setLineWidth(int width) { lineWidth = width; }

	// Outline tiles in object surface visualizations (default is yes)
	void outlineTiles(bool yesNo) { tile_outlines = yesNo; }

	// Render all surfaces.  If feedback is true, nothing is drawn and the
	// getPList() function can be used to access the generated coordinate
	// lists.  The last two arguments are only used when feedback is true.
	// If ignoreVisibility is false, only visibile figures will be used.
	// If ignoreBoundary is false, any existing boundary displacements will
	// be added to the resulting coordinates.

	// An instance of this class can only be used one way:  either feedback
	// is true or it is false.  Mixing these two kinds of rendering in the
	// same renderer may give unpredictable results.  Note that in either
	// case, the getPList() function may be used.
	void render(bool feedback = false, bool ignoreVisibility = true,
		bool ignoreBoundary = true);

	// Render surfaces to the specifying of constraints.  The marked figure
	// will be rendered normally.  The (current) figure, for which the
	// constraints are being set, will be partially or fully rendered,
	// depending on the value of complete_figures.  When partially rendered,
	// the cull_distance will determine the boundary at which the rendering
	// will stop.
	void partial_render(SelectedPartialFigures * figure_list,
		double cull_distance, bool complete_figures, int lvl);

	// Functions to access the point lists generated in feedback mode
	// by render(true), or the point lists used for regular rendering
	// produced by render(false).  The returned point lists belong to this
	// class, and so should not be deleted or modified.
	ThallCode::Pointlist_server2 * getPList(int figureId);

	// Function to place landmarks as dots on the rendered surface.  They are
	// taken from the figures to be rendered.  If size is greater than zero,
	// a 6-pointed jack is drawn with diameter 2*size, instead of coloring a
	// single spot.  If distinguishMarked is true, any marked landmarks will
	// be drawn as octahedrons to distinguish them, with line widths of
	// markedWidth.  The pixel width of these objects is the current line
	// width.  If size is 0.0, the pixel width of the dot produced will be 1.
	// If the OpenGL color is not given, it will be set to contrast with the
	// figural color.  
	void drawDots(double size = 0.0, float defaultWidth = 1.0f,
		bool distinguishMarked = false, float markedWidth = 3.0f,
		const float * color = NULL);

	// Function to place dots on the rendered surface.  The spot_list must be
	// an array of x,y,z coordinates of length 3*numSpots.  If size is greater
	// than zero, a 6-pointed jack is drawn with diameter 2*size, instead of
	// coloring a single pixel.  If the OpenGL color is not given, it must
	// be set externally before calling this function.  NumBigSpot may be
	// used to specify a dot to be rendered differently from the others.
	// It's value must refer to the X-coordinate of the dot.
	void drawDots(double * spot_list, int numSpots, double size = 0.0,
		float defaultWidth = 1.0f, float markedWidth = 3.0f,
		const float * color = NULL, int numBigSpot = -1);

	void renderFigureFromControlPoints(int numverts, int numfaces,
		Bpoint *pntlist, ThallCode::PindexList *faces);

private:

	void clean();
	void clean(int index);
	void createList(int length);

	void renderFigure(int index);
	ThallCode::Pointlist_server2 * computeFigure(int index,
		bool ignoreVisibility, bool ignoreBoundary);

	M3DObject * object;
	int * rows;
	int * cols;
	int num_figures;

	int surfaceLevel;
	int lineWidth;
    M3DSurfaceStyle surfaceStyle;
    M3DSurfaceStyle partialSurfaceStyle;
	bool tile_outlines;

	ThallCode::Pointlist_server2 ** list;
	int list_length;

	float surfaceColorsMin, surfaceColorsMax;
};

#endif

