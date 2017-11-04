#ifndef TILESET_H
#define TILESET_H

#include "Vector3D.h"
#include "Quat.h"


/*	Class TileSet can load a tile set from a file (.byu or .til)
	or generate one using marching cubes.  A tile set can contain
	multiple tiles figures.  When a tile set is generated or loaded
	from a file, any tile set previous contained in the class will
	be discarded.  (Note that .byu files may be created directly by
	exporting from a model.  This class is bypassed in that case.)

    Class TileSetRenderer may be used to generate displays of
	tile sets.
*/

class Vector3D;
class Image3D;

class TileSet
{
public:

    TileSet();
    TileSet(const char * filename, bool split = false, const Image3D* image = NULL);
    ~TileSet();

    bool loadFile(const char * filename, bool split = false, const Image3D* image = NULL);
    bool loadAppendFile(const char * filename, bool split = false, const Image3D* image = NULL);
    bool saveFile(const char * filename, bool asQuads, Image3D * image = NULL);
	// For a closed file; returns edges/tile
	int type(const char * filename, bool & isBYU) const;

	int figureCount() { return numFigures; }
	int tileCount() { return numTiles; }
	const Vector3D * getTiles() const { return coords; }
	const int * getTileCounts() const { return tileCounts; }
	int pointsPerTile() const { return quads ? 4 : 3; }
	bool haveQuads() const { return quads; }
	bool byuFormat() const { return have_byu; }

	// The range of coordinates of the union of the tiled figures
	void getBoundingBox(Vector3D & origin, Vector3D & opposite) {
		origin = leftBottomBack;
		opposite = rightTopFront;
	}

	Vector3D getCOG();
	void calcCOG();
	void calcCOG(bool tileSetMobileFigure[6]);

	// A (psuedo) marching cubes to aid in visualization of binary images
	bool simpleMarchingCubes(Image3D * image, bool asQuads = true);

	// modify the tiles: should these be in class BYU instead?
    void translateBy(const Vector3D & vTrans, bool tileSetMobileFigure[6]);
    void scaleBy(double scaleFactor, bool tileSetMobileFigure[6]);
    void scaleBy(double scaleFactor, const Vector3D & center, bool tileSetMobileFigure[6]);
    void rotateBy(const Quat &q, bool tileSetMobileFigure[6]);
    void rotateBy(const Quat &q, const Vector3D & center, bool tileSetMobileFigure[6]);
    void pirouetteBy(const Quat &q, bool tileSetMobileFigure[6]);

	// tell BYU that coord[] and numTiles changed
	void setGenerated();

private:

    Vector3D * coords;
	const int * tileCounts;
	int numFigures;
    int numTiles;
    bool local;
	bool quads;
	void * byu;
	int tc;
	bool have_byu;

    // Bounding box
    Vector3D leftBottomBack, rightTopFront;

	TileSet(const TileSet & ts);		// Undefined
	TileSet & operator=(const TileSet & ts);	// Undefined
};

#endif

