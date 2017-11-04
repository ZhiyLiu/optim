#ifndef TILESET_RENDERER_H
#define TILESET_RENDERER_H



/*	Class TileSetRenderer is used to generate Open GL renderings
	of TileSet objects.

    The renderer must be initialized by calling compileDispList()
	and optionally setDefaultColor(), before the actual drawing
	is performed by calls to render().  The colors used for rendering
	will be those provided through compileDispList(), unless the
	useDefault flag of render() is true, in which case the default
	color will be used for all figures.  Unless a default color is
	provided, blue will be used.  The default color will also be
	used for any figures having a figure number greater than the
	number of colors provided.  Unless solid is true, the rendering
	will use transparency, controlled by the value, [0, 1], of alpha.

    Function clear() discards any existing display lists and resets
	the default color to blue.
*/

class Vector3D;
class TileSet;

class TileSetRenderer
{
public:

    TileSetRenderer();
    ~TileSetRenderer();
	void clear();

	void setDefaultColor(const float * color);

	void compileDispList(TileSet * ts, int numFigColors = 0,
		const float ** figureColors = NULL);	// Set up for rendering
    void render(float alpha = 0.0f, bool useDefault = false,
		bool solid = false);

private:

    void compileFigures(int numFigures); 
	void compileFigure(int figureId);
    void compilePart(int lo, int hi);

    const Vector3D * coords;
	const int * tileCounts;
    int numTiles;
	bool quads;
	float dflt_color[3];
	float * pDfltColor[1];
	const float ** colors;
	int ncolors;

    int initialDispList;
	int nlists;

	void render(int index, float alpha, bool useDefault, bool solid);

	TileSetRenderer(const TileSetRenderer & ts);		// Undefined
	TileSetRenderer & operator=(const TileSetRenderer & ts);	// Undefined
};



#endif

