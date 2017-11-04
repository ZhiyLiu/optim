#ifndef BYU_H
#define BYU_H


#include "Vector3D.h"

class Vector3D;
class Image3D;


/*	This class supports the I/O of BYU files.  It also provides
	a simple marching cubes algorithm for the generation of
	BYU files from Image3D objects.

	Functions vertexList() and edgeList() return the actual data
	read from a BYU input file, a list of points and a list of
	polygons implemented as indexes into the list of points.

	Function tileList returns the merger of these two arrays,
	namely a list of tiles, where each tile is a set of points.

	The tiles will be made of either triangles or quadrangles,
	as indicated by functions type() and edgesPerTile().

    The operations of reading a BYU file and generating figures
	by marching cubes cannot be combined.  Each of them will
	cause data from the other to be discarded.

    Using the tie() function, an external file pointer may be
	used for all I/O.  If this is done, the caller is responsible
	for opening and closing the file.
*/

class M3DObject;
class Image3D;

class BYU
{
public:

	enum BYUPoly_t { Tris, Quads, NoTiles };

    BYU();
    BYU(const char * filename, bool split = false);

    virtual ~BYU();
	void clear();	// Discard all data

    bool read(const char * filename, bool split = false, const Image3D* image = NULL);	// Closed or tied
	void tie(void * file_ptr);		// Use an external FILE*
	void untie();					// Release an external FILE*
	bool merge(const BYU & byu);

	BYUPoly_t type() const;					// For the read file
	BYUPoly_t type(const char * filename);	// For a closed or tied file
	int edgesPerTile() { return edges_per_tile; }

	int figureCount() const { return nfigs; }	// Number of parts
	// Numbers of tiles per part (figure)
	int tileCount(int figureId) const { return tileCounts[figureId]; }
	const int * tilesPerFigure() const { return tileCounts; }

	int tileListLength() const { return numTiles; }	// Tiles in tileList() return
	const Vector3D * tileList() const { return coords; }	// Array of tiles

	int vertexCount() const { return numPoints; }	// Points vertexList() return
	const float * vertexList() const { return vertices; }	// Array of XYZ-triples

	int edgeCount() const { return nedges; }		// Length of edgeList() return
	const int * edgeList() const { return polygons; }	// Indexes into vertex list

	// call this if you change coords[] or numTiles
	void setGenerated() {generated = true;}

	// A (psuedo) marching cubes to aid in visualization of binary images
	bool simpleMarchingCubes(Image3D * image, bool asQuads = true);
	//void print(int tileCount);

	// The minimum and maximum coordinates of the union of the tiled figures
	void getBoundingBox(Vector3D & origin, Vector3D & opposite) const {
		origin = leftBottomBack;
		opposite = rightTopFront;
	}

//	bool exportModelToBYU(const char * filename, BYUPoly_t type, int level,
//		M3DObject * object, const bool discardEnds,
//		unsigned long figureMask = 0L, bool separate = false,
//		Image3D * image = NULL);
	enum point_t { XYZ, UVT, NORMALS };
	bool exportModelToBYU(const char * filename, BYUPoly_t type, int level,
		M3DObject * object, const bool discardEnds,
		unsigned long figureMask, bool separate,
		Image3D * image, point_t pt_type = XYZ );

	// Names used by exportModelToBYU() for writing separate output files
	const char * baseName() { return basename; }
	const char * extension() { return extn; }

	bool write(const char * filename, BYUPoly_t outputType, int numFigures,
		int * numPolygons, Vector3D * polygonList, BYUPoly_t type,
		Vector3D * origin = NULL, Vector3D * bound = NULL,
		Image3D * image = NULL);

private:

	BYU(const BYU & byu);				// Not implemented
	BYU & operator=(const BYU & byu);	// Not implemented

	void dump(bool verbose = false);	// For debugging

	void allocTiles(int wantNumTiles, int edgesPerTile);
	void addTile(Vector3D * newTile, int edgesPerTile);
	Vector3D * calcQuadTile(float x, float y, float z, int dir, Vector3D * m);

	bool write1(const char * filename, BYUPoly_t typeIn, bool split,
		int numParts, int * numVertices, int * numPolygons, float * vertexList,
		int * polygonList, Vector3D * origin = NULL, Vector3D * bound = NULL);

	bool write2(const char * filename, BYUPoly_t typeIn, bool split,
		int numParts, int * numVertices, int * numPolygons, float ** vertexLists,
		int ** polygonLists, Vector3D * origin = NULL, Vector3D * bound = NULL);

	bool construct();
	int uniqueVerts(int nTiles, float * points, int * tiles);
	void setBaseAndExt(const char * base, const char * ext);

	int nfigs, numPoints, nedges, edges_per_tile;
	float * vertices;
	int * polygons;
	int * tileCounts;
    Vector3D * coords;
    int numTiles;
    int tileAllocation;
	bool quads;
	void * fp;
	bool fptied;
	bool generated;
	char * basename;
	char * extn;

    // Bounding box
    Vector3D leftBottomBack, rightTopFront;
};


#endif

