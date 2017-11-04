#ifdef WIN32
#include <windows.h>
#endif
#include <stdio.h>
#include <iostream>
#include <sys/stat.h>
#include "Image3D.h"
#include "SubdivBoundary.h"
#include "M3DFigure.h"
#include "M3DObject.h"
#include "BYU.h"
#define D_POINTLIST_SERVER2
#define D_XFERLIST
#include "Shapedepend.h"
#include <assert.h>

#ifndef S_ISREG
#define S_ISREG(m)        (((S_IFREG & m) != 0) ? true : false)
#endif

#define MIN_TILESET_SIZE	1000

//#define DEBUG


using namespace std;

static char msg[] = "Not a valid BYU file -- see file:" __FILE__ " on line:";
static char msg1[] = "BYU file contained no tiles\n";

BYU::BYU()
{
	vertices = NULL;
	polygons = NULL;
    coords = NULL;
	tileAllocation = 0;
    numTiles = 0;
    numPoints = 0;
    nedges = 0;
	tileCounts = NULL;
	nfigs = 0;
	fp = NULL;
	fptied = false;
	edges_per_tile = 0;
	generated = false;
	basename = NULL;
	extn = NULL;
}

BYU::BYU(const char * filename, bool split)
{
	vertices = NULL;
	polygons = NULL;
    coords = NULL;
	tileAllocation = 0;
    numTiles = 0;
    numPoints = 0;
    nedges = 0;
	tileCounts = NULL;
	nfigs = 0;
	fp = NULL;
	fptied = false;
	edges_per_tile = 0;
	generated = false;
	basename = NULL;
	extn = NULL;

    (void) read(filename, split);
}

BYU::~BYU()
{
    if (vertices != NULL)
        delete [] vertices;
    if (polygons != NULL)
        delete [] polygons;
    if (coords != NULL)
        delete [] coords;
    if (tileCounts != NULL)
        delete [] tileCounts;
	if (basename != NULL) {
		delete [] basename;
		delete [] extn;
	}
}

void BYU::tie(void * file_ptr)
{
	fp = file_ptr;
	if (fp != NULL)
		fptied = true;
	else
		fptied = false;
}

void BYU::untie()
{
	if (! fptied)
		return;
	fp = NULL;
	fptied = false;
}

void BYU::setBaseAndExt(const char * base, const char * ext)
{
	int len;

	if (basename != NULL) {
		delete [] basename;
		delete [] extn;
	}
	len = strlen(base);
	basename = new char[len + 1];
	strcpy(basename, base);
	len = strlen(ext);
	extn = new char[len + 1];
	strcpy(extn, ext);
}

/*	Determine the type of a file.  This function only reads the header of
	a file to figure out the type.  This is needed for deciding whether to
	split quads into triangles when the file is read.

	The filename argument will only be used if an external file pointer has
	not been provided using tie(), in which case the filename will typically
	be NULL.
*/
BYU::BYUPoly_t BYU::type(const char * filename)
{
    struct stat buf;
    int nf, np, nt, ne;
	int tileEdges;

	if (fp == NULL) {
		if (filename == NULL)
			return NoTiles;
		else {
			if (stat(filename, &buf) != 0 || S_ISREG(buf.st_mode) == false) {
				cout << msg << __LINE__ << endl;
				return NoTiles;
			}
		}
	}

	if (fp == NULL)
		fp = fopen(filename, "r");
    if (fp == NULL)
        return NoTiles;

    // Read as a BYU format tile file
	fseek((FILE *) fp, 0L, SEEK_SET);	// rewind
    fscanf((FILE *) fp, "%d %d %d %d\n", &nf, &np, &nt, &ne);
	if (! fptied) {
		fclose((FILE *) fp);
		fp = NULL;
	}
	else
		fseek((FILE *) fp, 0L, SEEK_SET);	// rewind

    if (nt <= 0)
		tileEdges = 0;
	else
		tileEdges = ne/nt;

	if (tileEdges == 3)
		return Tris;
	else if (tileEdges == 4)
		return Quads;
	else {
		cout << msg1;
		return NoTiles;
	}
}

/*  Input a BYU file.  This function primarily sets values into variables nfigs,
	numPoints, and numTiles.  It also creates arrays tileCounts, vertices,
	polygons, and coords.  The data from the file will be available two ways:
	as a list of tiles made of sets of vertices in array coords, and as separate
	lists of vertices and sets of indexes forming tiles, in arrays vertices
	and polygons, respectively.

	If the BYU file was written in world coordinates, the data will be converted
	back to model coordinates.  The bounding box will also be calculated and
	made available in model coordinates.

	The input may contain quad- or tri-vertex tiles.  Variable edges_per_tile
	will be set accordingly.  If split is true and the file contains quads,
	the quads will be split, so that the tiles effectively returned by this
	function will be triangles.

	The filename argument will only be used if an external file pointer has
	not been provided using tie().


    Implementation Details:

	The memory allocations are as follows.  For the basic BYU format:
		tileCounts = new int[nfigs],
		vertices = new float[3*numPoints], and
		polygons = new int[nedges], where nedges is numTiles*edges_per_tile.
	For the TileSet form of the data:
		coords = new Vector3D[numTiles*edges_per_tile].
*/
bool BYU::read(const char * filename, bool split, const Image3D* image)
{
    struct stat buf;
    int v1, v2, v3, v4;
    float x, y, z;
	int i;
	int type;

	if (fp == NULL) {
		if (filename == NULL)
			return false;
		else {
			if (stat(filename, &buf) != 0 || S_ISREG(buf.st_mode) == false) {
				cout << msg << __LINE__ << endl;
				return false;
			}
		}
	}

	if (fp == NULL)
		fp = fopen(filename, "r");
    if (fp == NULL)
        return false;

	if (numTiles > 0)
		clear();

#ifdef DEBUG
	cout << "Loading BYU tiles\n";
#endif

    // Read as a BYU format tile file
    fscanf((FILE *) fp, "%d %d %d %d\n", &nfigs, &numPoints, &numTiles, &nedges);
    if (numTiles <= 0) {
		if (! fptied)
			fclose((FILE *) fp);
        cout << msg1;
        return false;
	}
#ifdef DEBUG
	cout << "nfigs = " << nfigs << "    numPoints = " << numPoints << endl;
#endif
	tileCounts = new int[nfigs];
	edges_per_tile = nedges/numTiles;
	if (edges_per_tile < 3 || edges_per_tile > 4) {
		if (! fptied)
			fclose((FILE *) fp);
        cout << msg << __LINE__ << endl;
        return false;
	}

    quads = false;
	if (edges_per_tile == 4)	// We use either triangles or quads only
		quads = true;
	if (split && ! quads)
		split = false;		// Ignore splitting of triangles

	for (i = 0; i < nfigs; i++) {
		fscanf((FILE *) fp, "%d %d", &v1, &v2);
		tileCounts[i] = v2 - v1 + 1;
	}

	if (split)
		numTiles += numTiles;
    coords = new Vector3D[numTiles*edges_per_tile];
#ifdef DEBUG
	cout << "Length of coords = " << edges_per_tile << "x" << numTiles << " = "
		<< numTiles*edges_per_tile << endl;
#endif

    // Read all vertices
#ifdef DEBUG
	cout << "Allocating " << 3*numPoints << " vertices" << endl;
#endif
	vertices = new float[3*numPoints];

    for (i = 0; i < 3*numPoints; i += 3)
    {
        if (3 != fscanf((FILE *) fp, "%f %f %f",
			&vertices[i], &vertices[i + 1], &vertices[i + 2]))
		{
			clear();
			if (! fptied)
				fclose((FILE *) fp);
			cout << msg << __LINE__ << endl;
            return false;
		}
    }

    // Read each tile and poke it's vertices into the coords array
    int n, pt;
	int polygonIndex;
	bool first;

	if (split)
		nedges = (nedges*6)/4;
	polygons = new int[nedges];

	first = true;
	polygonIndex = 0;
	pt = 0;
    for (n = 0; n < nfigs; n++) {
		for (i = 0; i < tileCounts[n]; i++)
		{
			fscanf((FILE *) fp, "%d %d %d", &v1, &v2, &v3);
			v1--;
			v2--;
			if (quads) {
				fscanf((FILE *) fp, " %d", &v4);
				v3--;
				v4 = -v4 - 1;
			}
			else
				v3 = -v3 - 1;
			polygons[polygonIndex++] = v1;
			polygons[polygonIndex++] = v2;
			polygons[polygonIndex++] = v3;
			if (quads) {
				if (split) {
					polygons[polygonIndex++] = v1;
					polygons[polygonIndex++] = v3;
					polygons[polygonIndex++] = v4;
				}
				else
					polygons[polygonIndex++] = v4;
			}

			// Update bounding box as necessary and insert tiles into coords
			x = vertices[3*v1];
			y = vertices[3*v1 + 1];
			z = vertices[3*v1 + 2];
			if (first) {
				leftBottomBack.set(x, y, z);
				rightTopFront.set(x, y, z);
				if (quads) {	// Check first vertex line for final minus sign
					if (-v4 - 1 >= 0) {
						clear();
						if (! fptied)
							fclose((FILE *) fp);
						cout << msg << __LINE__ << endl;
						return false;
					}
				}
				else {
					if (-v3 - 1 >= 0) {
						clear();
						if (! fptied)
							fclose((FILE *) fp);
						cout << msg << __LINE__ << endl;
						return false;
					}
				}
				first = false;
			}
			else {
				leftBottomBack.minPoint(Vector3D(x,y,z));
				rightTopFront.maxPoint(Vector3D(x,y,z));
			}

			coords[pt++].set(x, y, z);
			
			x = vertices[3*v2];
			y = vertices[3*v2 + 1];
			z = vertices[3*v2 + 2];
			leftBottomBack.minPoint(Vector3D(x,y,z));
			rightTopFront.maxPoint(Vector3D(x,y,z));
			coords[pt++].set(x, y, z);

			x = vertices[3*v3];
			y = vertices[3*v3 + 1];
			z = vertices[3*v3 + 2];
			leftBottomBack.minPoint(Vector3D(x,y,z));
			rightTopFront.maxPoint(Vector3D(x,y,z));
			coords[pt++].set(x, y, z);
			
			if (! quads) 
				continue;

			if (split) {
				coords[pt].set(coords[pt - 3]);
				pt++;
				coords[pt++].set(x, y, z);
			}

			x = vertices[3*v4];
			y = vertices[3*v4 + 1];
			z = vertices[3*v4 + 2];
			leftBottomBack.minPoint(Vector3D(x,y,z));
			rightTopFront.maxPoint(Vector3D(x,y,z));
			coords[pt++].set(x, y, z);
		}
	}

	// Extension to original BYU format
	if (1 != fscanf((FILE *) fp, "%d", &type))
		type = 0;

	if (type == 1) {
		double a[3];
		Vector3D origin, bound;
		double modelToWorldScale;

		// Read extra line and convert tiles to model coordinates
		fscanf((FILE *) fp, "%lf %lf %lf", &a[0], &a[1], &a[2]);
		origin = Vector3D(a[0], a[1], a[2]);
		fscanf((FILE *) fp, "%lf %lf %lf", &a[0], &a[1], &a[2]);
		bound = Vector3D(a[0], a[1], a[2]);
		Vector3D extents = bound - origin;

		// Find the longest extent.  Note: The Y and Z axes may have
		// negative inter-voxel spacings; X cannot.
		modelToWorldScale = extents.getX();
		double m = fabs(extents.getY());
		if (m > modelToWorldScale)
			modelToWorldScale = m;
		m = fabs(extents.getZ());
		if (m > modelToWorldScale)
			modelToWorldScale = m;

		if (modelToWorldScale <= 0.0)	// AGG: Why is this needed?
			modelToWorldScale = 1.0;	// Incorrect value? assume "no scaling"

#ifdef DEBUG
		cout << "Converting BYU tile vertices from world coordinates\n";
		cout << "modelToWorldScale = " << modelToWorldScale << '\n';
		cout << "origin = "; origin.print();
		cout << "bound = "; bound.print();
#endif

		// This code must conform to that in Image3D::worldToModelCoordinates()

		// Compute the signed scale factors
		extents.setX(1.0);	// X cannot be flipped	
		if (extents.getY() < 0.0)	
			extents.setY(-1.0);
		else
			extents.setY(1.0);
		if (extents.getZ() < 0.0)	
			extents.setZ(-1.0);
		else
			extents.setZ(1.0);
		extents /= modelToWorldScale;

		// Convert the data to model coordinates
		for (i = 0; i < pt; i++) {
			coords[i] -= origin;
			coords[i] = coords[i].vprod(extents);
		}

		leftBottomBack -= origin;
		leftBottomBack = leftBottomBack.vprod(extents);

		rightTopFront -= origin;
		rightTopFront = rightTopFront.vprod(extents);
	}
	else if(image) {
		// Scale by world->model tx of loaded image.
		for (i = 0; i < pt; i++) {
			image->worldToModelCoordinates(coords[i]);
		}
	}


	if (split) {
		edges_per_tile = 3;
		quads = false;
		for (i = 0; i < nfigs; i++)
			tileCounts[i] += tileCounts[i];
	}

/*	polygonIndex = 0;
	cout << "Array vertices:\n";
	for (i = 0; i < 3*numPoints; i += 3)
		cout << vertices[i] << ", " << vertices[i + 1] << ", " << vertices[i + 2] << '\n';
	cout << "Array polygons:\n";
    for (n = 0; n < nfigs; n++) {
		cout << "  Figure " << n << '\n';
		for (i = 0; i < tileCounts[n]*edges_per_tile; i += edges_per_tile) {
			cout << polygons[polygonIndex++] << ", "  << polygons[polygonIndex++]
				<< ", " << polygons[polygonIndex++];
			if (quads)
				cout << ", " << polygons[polygonIndex++];
			cout << '\n';
		}
	}
*/

	if (! fptied) {
		fclose((FILE *) fp);
		fp = NULL;
	}
//	dump(true);
    return true;
}

void BYU::dump(bool verbose)
{
	int i;

	cout << "BYU contents:\n";
	cout << "  nfigs = " << nfigs << '\n';
	cout << "  numPoints = " << numPoints << '\n';
	cout << "  nedges = " << nedges << '\n';
	cout << "  numTiles = " << numTiles << '\n';
	cout << "  edges_per_tile = " << edges_per_tile << "  ";
	if (quads)
		cout << "Quads   ";
	else
		cout << "Tries   ";
	if (fptied)
		cout << "Tied   ";
	if (generated)
		cout << "Generated";
	cout << '\n';
	if (verbose) {
		cout << "  Vertices\n";
		for (i = 0; i < 3*numPoints; i += 3)
			cout << '(' << vertices[i] << ", " << vertices[i + 1] << ", " << vertices[i + 2] << ")\n";
		cout << "  Polygons\n";
		if (quads)
			for (i = 0; i < nedges; i += 4)
				cout << polygons[i] << ' ' << polygons[i + 1] << ' ' << polygons[i + 2] << ' ' << polygons[i + 3] << '\n';
		else
			for (i = 0; i < nedges; i += 3)
				cout << polygons[i] << ' ' << polygons[i + 1] << ' ' << polygons[i + 2] << '\n';
		cout << "  TileSet tiles\n";
		for (i = 0; i < nedges; i++)
			coords[i].print();
	}
	cout << '\n';
}

bool BYU::merge(const BYU & byu)
{
	int i, j, k;
	int start;
	float * new_vertices;
	int * new_polygons;
	int * new_tileCounts;
    Vector3D * new_coords;

	if (numTiles <= 0) {
		cout << "Target of BYU merge contains no tiles" << endl;
		return false;	// Since there is no assignment operator
	}

	if (fptied) {
		cout << "Cannot merge to a tied BYU object" << endl;
		return false;	// Since this class cannot merge files
	}

	if (edges_per_tile != byu.edges_per_tile) {
		cout << "Merging requires BYU objects to be the same edges_per_tile" << endl;
		return false;
	}

//	dump(true);
//	((BYU &) byu).dump(true);
	start = numPoints;

	// KEEP FOR REFERENCE: original code did not multiply by 3
	//	new_vertices = new float[numPoints + byu.numPoints];
	//	for (i = 0; i < numPoints; i++)
	//		new_vertices[i] = vertices[i];
	//	j = numPoints;
	//	for (i = 0; i < byu.numPoints; i++)
	//		new_vertices[j++] = byu.vertices[i];

	new_vertices = new float[3*(numPoints + byu.numPoints)];
	for (i = 0; i < 3*numPoints; )
	{
		new_vertices[i] = vertices[i++];
		new_vertices[i] = vertices[i++];
		new_vertices[i] = vertices[i++];
	}
	j = 3*numPoints;
	for (i = 0; i < 3*byu.numPoints; )
	{
		new_vertices[j++] = byu.vertices[i++];
		new_vertices[j++] = byu.vertices[i++];
		new_vertices[j++] = byu.vertices[i++];
	}

	
	delete [] vertices;
	vertices = new_vertices;
	numPoints += byu.numPoints;

	new_polygons = new int[nedges + byu.nedges];
	for (i = 0; i < nedges; i++)
		new_polygons[i] = polygons[i];
	j = nedges;
	k = 0;
	for (i = 0; i < byu.nedges; i += edges_per_tile) {
		new_polygons[j++] = start + byu.polygons[k++];
		new_polygons[j++] = start + byu.polygons[k++];
		new_polygons[j++] = start + byu.polygons[k++];
		if (quads)
			new_polygons[j++] = start + byu.polygons[k++];
	}
	delete [] polygons;
	polygons = new_polygons;

	new_coords = new Vector3D[nedges + byu.nedges];
	for (i = 0; i < nedges; i++)
		new_coords[i] = coords[i];
	j = nedges;
	for (i = 0; i < byu.nedges; i++)
		new_coords[j++] = byu.coords[i];
	delete [] coords;
	coords = new_coords;
	numTiles += byu.numTiles;
	nedges += byu.nedges;

	new_tileCounts = new int[nfigs + byu.nfigs];
	for (i = 0; i < nfigs; i++)
		new_tileCounts[i] = tileCounts[i];
	j = nfigs;
	for (i = 0; i < byu.nfigs; i++)
		new_tileCounts[j++] = byu.tileCounts[i];
	delete [] tileCounts;
	tileCounts = new_tileCounts;
	nfigs += byu.nfigs;

	leftBottomBack.minPoint(byu.leftBottomBack);
	rightTopFront.maxPoint(byu.rightTopFront);

	return true;
}

/*	This function converts from the data being in array coords and
	numTiles, to it being in vertices, polygons, numTiles, numPoints,
	coords, and tileCounts, just as after loading a BYU file.  Only
	nedges is not set.  Arrays polygons, vertices and tileCounts are
	assumed to have not previously been allocated.

	Note: the data is considered to be a single part, even though
	the data may have multiple connected components.
*/
bool BYU::construct()
{
	if (numTiles <= 0)
		return false;

	// Compute the length of array coord
	numPoints = numTiles*edges_per_tile;

	polygons = new int[numPoints];
	// This may make vertices larger than necessary
	vertices = new float[3*numPoints];

	tileCounts = new int[1];
	tileCounts[0] = numTiles;

	// Populate the arrays using the data in coords
	numPoints = uniqueVerts(numTiles, vertices, polygons);
	return true;
}

/*	This function converts coords into lists of discrete points, three
	coordinates per point, and tiles, with each tile represented as a
	set of indexes into the vertex list.  It returns the number of
	vertices in array points.  The first argument should contain the
	number of tiles in the coords array and thus returned in the last
	argument.  Variable edges_per_tile must contain 3 or 4.  This is
	theoretically an n^2 algorithm, but because of the reverse order
	of the inner loop, it is normally not so bad in practice for most
	medical objects.
*/
int BYU::uniqueVerts(int nTiles, float * points, int * tiles)
{
	int i, j, k, index;
	double x, y, z;
	int nPoints, nt;
	bool vmatch;

	nPoints = 0;
	nt = 0;
	for (i = 0; i < nTiles; i ++) {
		if (nTiles == 30000)
			cout << "calc'ing uniqueVerts on " << nTiles << " tiles..." << flush;
		if (nTiles > 30000 && i % 5000 == 0)
			cout << ' ' << i << flush;

		index = edges_per_tile*i;
		for (j = 0; j < edges_per_tile; j++) {
			x = coords[index + j].getX();
			y = coords[index + j].getY();
			z = coords[index + j].getZ();
			k = nPoints - 1;
			vmatch = false;
			while ((k > 0) && (! vmatch)) {
				if ((fabs(points[k] - z) < 0.000001) &&
					(fabs(points[k - 1] - y) < 0.000001) &&
					(fabs(points[k - 2] - x) < 0.000001))
						vmatch = true;
				else
				   k -= 3;
			}
			if (vmatch) {
				tiles[nt++] = k/3;
			}
			else {
				points[nPoints] = x;
				points[nPoints + 1] = y;
				points[nPoints + 2] = z;
				tiles[nt++] = nPoints/3;
				nPoints = nPoints + 3;
			}
		}
	}
	if (nTiles > 30000)
		cout << endl;
	return nPoints/3;
}

/*	Output a BYU file.  This function provides the primary way to write the
	contents of this class into a BYU file specified by filename.  The data
	may have been produced by calling the marching cubes function or may have
	been passed in from an outside source, typically after loading either
	another BYU file or a .til file.

	The number of parts should be provided in numFigures.  Array numPolygons
	must contain one entry per figure, each entry being the number of tiles
	the figure contains.

	The outputType argument should be set to indicate whether triangular or
	quadrilateral tiles should be written.  Argument type indicates the form
	of the data in polygonList.  If quads are to be written as triangles,
	the tiles will be split in half.  Triangular tiles cannot be written as
	Quads; the required merge capability is not provided.

	If the data is a result of calling simpleMarchingCubes(), then argument
	polygonList should be array coords, obtained by calling tileList().

	If origin and bound specify a model-space bounding box for the data
	and an image is provided, then the data will be written in world
	coordinates.
*/
bool BYU::write(const char * filename, BYUPoly_t outputType, int numFigures,
	int * numPolygons, Vector3D * polygonList, BYUPoly_t type,
	Vector3D * origin, Vector3D * bound, Image3D * image)
{
	float ** points;
	int ** tiles;
	bool ret;
	int i;
	bool haveCoords;


#ifdef DEBUG
	cout << "BYU::write()" << '\n';
#endif
	if (filename == NULL)
		return false;

	if (type == NoTiles)
		return false;

	if (numFigures <= 0)
		return false;

	if (outputType == Quads && type == Tris) {
		cout << "Cannot merge triangular tiles to get quad-tiles\n";
		return false;
	}

	if (numPolygons[0] <= 0)
		return false;

	if (generated || polygonList != coords) {
		/*	Tile data was generated by marching cubes (as quads) or read from
			a .til file (as triangles).  In the case of .til files, class
			TileSet owns the polygonList array; for generated tile sets, it
			is the coords array of this class.  In both cases, the separate
			vertex and tile index arrays must be constructed.
		*/
#ifdef DEBUG
		cout << "    case 1\n";
#endif
		if (! generated) {
			clear();
			haveCoords = false;	// A .til file was read; coords will be NULL
		}
		else	// simpleMarchingCubes() called clear()
			haveCoords = true;	// Tiles were generated by this BYU object

		coords = polygonList;
		if (type == Quads)
			edges_per_tile = 4;
		else
			edges_per_tile = 3;
		numTiles = numPolygons[0];
		nfigs = numFigures;

		if (! construct())
			return false;
		points = &vertices;
		tiles = &polygons;

		if (origin && bound && image) {
			// Convert to world coordinates
#ifdef DEBUG
			cout << "Converting BYU tile vertices to world coordinates\n";
#endif
			for (i = 0; i < 3*numPoints; i += 3) {
				Vector3D pnt(vertices[i], vertices[i + 1], vertices[i + 2]);
				image->modelToWorldCoordinates(pnt);
				vertices[i] = (float) pnt.getX();
				vertices[i + 1] = (float) pnt.getY();
				vertices[i + 2] = (float) pnt.getZ();
			}
		}
		else {
			origin = NULL;
			bound = NULL;
		}

		ret = write2(filename, type, (type == outputType ? false: true), numFigures,
			&numPoints, numPolygons, points, tiles, origin, bound);
		if (! haveCoords) {
			coords = NULL;	// A .til file was read; the class does not own coords
			clear();		// Delete arrays created by construct()
		}
		// else don't clear because the generated tileset may be being displayed
	}
	else {	// A .byu file was read
		bool split;

		// The tile data was read from a file by this class object.  No
		// array construction is necessary.
#ifdef DEBUG
		cout << "    case 2\n";
#endif
		if (type == Quads && outputType == Tris)
			split = true;
		else
			split = false;

		tiles = new int *[numFigures];
		tiles[0] = polygons;
		for (i = 1; i < numFigures; i++)
			tiles[i] = polygons + tileCounts[i - 1]*edges_per_tile;

		// Build array num_points containing the number of vertices in each part
		int * num_points = new int[numFigures];
		int max = 0;
		int low = 0;
		int n = 0;
		for (i = 0; i < numFigures - 1; i++) {
			for (int j = tileCounts[i]*edges_per_tile - 1; j >= low; j--) {
				int v = tiles[i][j];	// Vertex index
				if (v > max)
					max = v;
			}
			num_points[i] = max + 1;
			n += num_points[i];
			low = tileCounts[i]*edges_per_tile;
		}
		num_points[numFigures - 1] = numPoints - n;

		if (origin && bound && image) {
#ifdef DEBUG
			cout << "Converting BYU tile vertices to world coordinates\n";
#endif
			// Convert to world coordinates
			for (i = 0; i < 3*numPoints; i += 3) {
				Vector3D pnt(vertices[i], vertices[i + 1], vertices[i + 2]);
				image->modelToWorldCoordinates(pnt);
				vertices[i] = (float) pnt.getX();
				vertices[i + 1] = (float) pnt.getY();
				vertices[i + 2] = (float) pnt.getZ();
			}
		}
		else {
			origin = NULL;
			bound = NULL;
		}

/*		cout << "Array num_points:\n";
		for (n = 0; n < numFigures; n++)
			cout << num_points[n] <<' ';
		cout << '\n';
		cout << "Array tileCounts:\n";
		for (n = 0; n < numFigures; n++)
			cout << tileCounts[n] <<' ';
		cout << '\n';
*/
		ret = write1(filename, type, split, numFigures, num_points, tileCounts,
			vertices, polygons, origin, bound);
		delete [] tiles;
		delete [] num_points;
	}

	return ret;
}

/*  Produce a BYU file containing tiles of the specified type of the model
    boundary rendered at the specified level.  If an image is provided, its
	world coordinate system will be recorded in the output file.  The tiles
	will be generated as quads and then written to the specified file using
	the function above.  If type is Tris, then the tiles will be converted
	to triangles at the time of output.

	Argument figureMask is a bit mask for the selection of figures to be
	exported.  If bit 0 is 1, then figure 0 will be exported, and so on,
	unless all bits are zero, in which case all figures are exported.  (This
	function should only be called when intending to output at least one
	figure).  Use of the mask requires that there be no more figures in the
	object than the size of the argument.  No test is performed to verify
	that this is true.

	If separate is true, the figures will be put into separate files.  The
	name of each file will be constructed by inserting the figure number
	(always 2 digits) just before the extension, or at the end of the
	filename, if no extension is present.
*/
bool BYU::exportModelToBYU(const char * filename, BYUPoly_t type, int level,
	M3DObject * object, const bool discardEnds,
	unsigned long figureMask, bool separate, Image3D * image,
	point_t pt_type )
{
    ThallCode::Pointlist_server2 * pList;
    M3DFigure * figure;
	int n, npoints, numtiles;
    Bpoint * pointset;
	int * tile_counts;
	int * point_counts;
	bool useWorldCoords;
    int numFigures, i, figureNo, nFig;
	int * figureNumbers;
	Vector3D o, b;
	Vector3D * origin;
	Vector3D * bound;
	float ** points;
	int ** tiles;
	int * index;
	unsigned long bitSelector;

	if (type == NoTiles)
		return false;

    if(object == NULL)
        return false;

	if (fp == NULL && filename == NULL)
		return false;

	if (image != NULL)
		useWorldCoords = true;
	else
		useWorldCoords = false;

    numFigures = object->getFigureCount();
    pList = new ThallCode::Pointlist_server2[numFigures];

	//
	// An array that stores the new index for the vertex after some vertices have been dropped.
	//
	int** renamedVertex	= new int*[numFigures];

    point_counts = new int[numFigures];	  // Counts of vertices per figure
    points = new float *[numFigures];	  // One list of triples per figure
	figureNumbers = new int[numFigures];  // Saves actual figure Ids to be written
	bitSelector = 1L;
	nFig = 0;
	for (figureNo = 0; figureNo < numFigures; bitSelector <<= 1, figureNo++)
	{
		if (figureMask > 0L && !(figureMask & bitSelector))
			continue;

		figure = object->getFigurePtr(figureNo);

		// Generate the surfaces
		Xferlist * xferList = convertM3DtoXfer(figure);
		pList[nFig].init(xferList);
		pList[nFig].ComputeSubdivPointCloud(level);
		pList[nFig].ComputeSubdivBoundaryTiles(level);
		delete [] (xferList->atomlist);
		delete xferList;

		// Get boundary displacement information
		SubdivBoundary * boundary = figure->getBoundaryPtr();
		if (boundary != NULL)
		{
			Displacements * disp = boundary->getDisplacements(level);
			if (disp != NULL)
				pList[nFig].applyDisplacements(disp->getNumPts(), disp->getVals());
		}

		// Generate the vertex lines (list of points)
		pList[nFig].subdivboundaryinfo(&npoints, &pointset);

		int skipCount	= 0;
		int storedPoints	= npoints;
		renamedVertex[figureNo]	= new int[npoints];
		if(discardEnds) {
			for( i = 0; i < npoints; i++ ) {
				if( pList[nFig].thisMesh->isEndMedcoord( pointset[i]) ) {
					skipCount++;
					renamedVertex[figureNo][i]	= -1;
				}
				else {
					renamedVertex[figureNo][i]	= i - skipCount;
				}
			}
			storedPoints	-= skipCount;
		}
		else {
			for( i = 0; i < npoints; i++ ) {
				renamedVertex[figureNo][i]	= i;	// self-reference.
			}
		}

		point_counts[nFig] = storedPoints;
		points[nFig] = new float[3*storedPoints];
		cout << "Writing " << storedPoints << "(" << npoints << ") vertices for figure " << figureNo << '\n';
		if (pt_type == XYZ) {
			if (useWorldCoords) {
	#ifdef DEBUG
				cout << "Converting BYU tile vertices to world coordinates\n";
	#endif
				for (i = 0, n = 0; i < npoints; i++) {
					if( renamedVertex[nFig][i] != -1 ) {
						Vector3D pnt(pointset[i].pnt);
						image->modelToWorldCoordinates(pnt);
						points[nFig][n++] = (float) pnt.getX();
						points[nFig][n++] = (float) pnt.getY();
						points[nFig][n++] = (float) pnt.getZ();
					}
				}
			}
			else {
				for (i = 0, n = 0; i < npoints; i++) {
					if( renamedVertex[nFig][i] != -1 ) {
						points[nFig][n++] = (float) pointset[i].pnt[0];
						points[nFig][n++] = (float) pointset[i].pnt[1];
						points[nFig][n++] = (float) pointset[i].pnt[2];
					}
				}
			}
		}
		else if (pt_type == UVT) {
			for (i = 0, n = 0; i < npoints; i++) {
				if( renamedVertex[nFig][i] != -1 ) {
					points[nFig][n++] = (float) pointset[i].u;
					points[nFig][n++] = (float) pointset[i].v;
					points[nFig][n++] = (float) pointset[i].t;
				}
			}
		}
		else if (pt_type == NORMALS) {
			for (i = 0, n = 0; i < npoints; i++) {
				if( renamedVertex[nFig][i] != -1 ) {
					points[nFig][n++] = (float) pointset[i].norm[0];
					points[nFig][n++] = (float) pointset[i].norm[1];
					points[nFig][n++] = (float) pointset[i].norm[2];
				}
			}
		}
		else {
			cout << "Unknown point type : " << pt_type << endl;
			return false;
		}
		assert( n == 3*storedPoints );

		figureNumbers[nFig] = figureNo;
		nFig++;
	}
	numFigures = nFig;	// Fewer may be written than in the object
	if (numFigures == 0)	// Invalid bit mask
		return false;

	// Generate the polygon lines (quad tiles) as indexes into list of points
    tile_counts = new int[numFigures];    // Counts of tiles per figure
    tiles = new int *[numFigures];    // One set of 4 indexes per figure
	for (nFig = 0; nFig < numFigures; nFig++) {
		figure = object->getFigurePtr(figureNumbers[nFig]);
		pList[nFig].subdivtileindexinfo(&numtiles, &index);
		// Do not export ends ?
		if( discardEnds ) {
			int* vindex	= new int[4*numtiles];
			for( i = 0, n = 0; i < 4*numtiles; i += 4 ) {
				if( renamedVertex[nFig][index[i]]   == -1 || renamedVertex[nFig][index[i+1]] == -1 ||
					renamedVertex[nFig][index[i+2]] == -1 || renamedVertex[nFig][index[i+3]] == -1 ) {
					// Skip this tile, it contains a vertex that has been skipped.
				}
				else {
					// include this tile.
					vindex[4*n]		= renamedVertex[nFig][index[i]];
					vindex[4*n+1]	= renamedVertex[nFig][index[i+1]];
					vindex[4*n+2]	= renamedVertex[nFig][index[i+2]];
					vindex[4*n+3]	= renamedVertex[nFig][index[i+3]];
					n++;
				}
			}
			tile_counts[nFig] = n;
			tiles[nFig] = vindex;
		}
		else {
			tile_counts[nFig] = numtiles;
			tiles[nFig] = index;
		}
	}

	if (useWorldCoords) {
		o = image->getModelOrigin();
		b = image->getModelBound();
		origin = &o;
		bound = &b;
#ifdef DEBUG
		cout << "origin = "; origin->print();
		cout << "bound = "; bound->print();
#endif
	}
	else {
		origin = NULL;
		bound = NULL;
	}

	// Output the data
	bool ret = false;
	if (separate) {
		float * pts;
		int * tls;
		char * name;
		const char * ext;
		char * next;
		int len;
		char fmt1[] = "%d%s";
		char fmt2[] = "%.2d%s";
		char * fmt;

		name = new char[strlen(filename) + 3];
        ext = strrchr(filename, '.');
        if (ext == NULL)
			ext = filename + strlen(filename);
		len = ext - filename;
		strncpy(name, filename, len);
		name[len] = '\0';
		next = name + len;

		setBaseAndExt(name, ext);

		if (numFigures < 10)
			fmt = fmt1;
		else
			fmt = fmt2;

		for (nFig = 0; nFig < numFigures; nFig++) {
			// Add figure number and extension to the file name
			sprintf(next, fmt, figureNumbers[nFig], ext);

			pts = points[nFig];
			tls = tiles[nFig];
			ret = ret && write2(name, Quads, (type == Tris ? true : false), 1,
				point_counts + nFig, tile_counts + nFig, &pts, &tls, origin, bound);
		}
	}
	else
		ret = write2(filename, Quads, (type == Tris ? true : false), numFigures,
			point_counts, tile_counts, points, tiles, origin, bound);

	for (figureNo = 0; figureNo < numFigures; figureNo++) {
		delete [] points[figureNo];
		delete [] renamedVertex[figureNo];
		// Do not export ends ?
		if(discardEnds) {
			delete[] tiles[figureNo];
		}
	}
	delete [] points;
	delete [] tiles;
	delete [] point_counts;
	delete [] tile_counts;
	delete [] pList;
	delete [] figureNumbers;
	delete [] renamedVertex;

	return ret;
}

/*  Low-level output function for BYU files.  A BYU file can describe multiple figures,
    called parts.  Each part must have a list of vertices and a list of polygons; each
	polygon is represented as a set of indexes into the list of vertices.  These items
	should be provided to this function as follows:
		numParts		The number of parts to be output.
		numVertices		An array of length numParts containing the numbers of
							vertices of each part.  There are 3 values per vertex.
		numPolygons		An array of length numParts containing the numbers of
							polygons of each part.  The BYU format can support any
							number of points per polygon, but this function accepts
							only 3 or 4 vertices per polygon.  Quad-tiles may also
							be split in half to produce triangles (see below).
		vertexList		An array containing the vertices of all parts.  The length of
							the array should be three times the sum of the elements in
							numVertices, since there are 3 coordinates per vertex.
		polygonList		An array containing the polygons of all parts.  Polygons are
							represented as indexes into the vertexList array.  This
							array should have a length 3 or 4 times the sum of entries
							in array numPolygons, corresponding to the 3 or 4 vertices
							per input polygon.

	This function can write its output as 4-tiles, 3-tiles, or can split 4-tiles to
	produce pairs of 3-vertex tiles.  The typeIn argument should be used to indicate
	the current form of the polygon data.  If typeIn is Quad, then split may be used
	to indicate that the 4-tiles are to be split.

    If an origin and bound are provided, they will be used to specify an extension to
	the regular BYU format.  The tile vertices will be assumed to already be in world
	coordinates; they will not be converted.

	This function is much the same as write2(), differing only in  the structure of
	arguments vertexLists and polygonLists.
*/
bool BYU::write1(const char * filename, BYUPoly_t typeIn, bool split, int numParts,
	int * numVertices, int * numPolygons, float * vertexList, int * polygonList,
	Vector3D * origin, Vector3D * bound)
{
	int i;
	int partNo;
	int nverts, npolys;

#ifdef DEBUG
	cout << "BYU::write1()" << '\n';
#endif
	if (typeIn == NoTiles)
		return false;
	if (typeIn == Tris && split)
		return false;

	if (fp == NULL && filename == NULL)
		return false;

    if (fp == NULL)
		fp = fopen(filename, "wb");
    if (fp == NULL)
        return false;

	nverts = 0;    // Length of array vertexList
	npolys = 0;    // Length of array polygonList
	for (i = 0; i < numParts; i++) {
		nverts += numVertices[i];
		npolys += numPolygons[i];
	}

	// Output the header lines
	// First, numbers of: figures, vertices, quads/triangles, and edges.
	if (typeIn == Quads) {
		if (split)
			fprintf((FILE *) fp, "%7d %7d %7d %7d\n", numParts, nverts, 2*npolys, 6*npolys);
		else
			fprintf((FILE *) fp, "%7d %7d %7d %7d\n", numParts, nverts, npolys, 4*npolys);
	}
	else
		fprintf((FILE *) fp, "%7d %7d %7d %7d\n", numParts, nverts, npolys, 3*npolys);
	// Second, the sequence of beginning and ending polygon numbers for each part
	i = 1;
	for (partNo = 0; partNo < numParts; partNo++) {
		// The BYU format requires 7 bytes per value
		if (typeIn == Quads && split) {
			// Splitting quads
			fprintf((FILE *) fp, "%7d %7d ", i, 2*numPolygons[partNo] + i - 1);
			i += 2*numPolygons[partNo];
		}
		else {	// Tris or quads
			fprintf((FILE *) fp, "%7d %7d ", i, numPolygons[partNo] + i - 1);
			i += numPolygons[partNo];
		}
	}
	fprintf((FILE *) fp, "\n");

	// Output the vertex lines (list of points)
	for (i = 0; i < 3*nverts; i += 3)
		fprintf((FILE *) fp, "%f %f %f\n", vertexList[i], vertexList[i + 1],
			vertexList[i + 2]);

	// Output the polygon lines (indexes into list of vertices)
	// Note that the BYU format requires indexing counted from 1, instead of 0
	int start = 1;	// The BYU format requires 1-origin indexing
	if (typeIn == Quads) {	// quads
		if (split) {
#ifdef DEBUG
			cout << "Writing " << 2*npolys << " 3-tiles\n";
#endif
			for (i = 0; i < 4*npolys; i += 4) {
				fprintf((FILE *) fp, "%d %d -%d\n", start + polygonList[i + 3],
					start + polygonList[i + 0], start + polygonList[i + 1]);
				fprintf((FILE *) fp, "%d %d -%d\n", start + polygonList[i + 1],
					start + polygonList[i + 2], start + polygonList[i + 3]);
			}
		}
		else {
#ifdef DEBUG
			cout << "Writing " << npolys << " 4-tiles\n";
#endif
//			cout << "Polygons for part " << partNo << '\n';
			for (i = 0; i < 4*npolys; i += 4) {
//				cout << polygonList[i + 0] << ", " << polygonList[i + 1] << ", "
//					<< polygonList[i + 2] << ", " << polygonList[i + 3] << '\n';
				fprintf((FILE *) fp, "%d %d %d -%d\n", start + polygonList[i + 0],
					start + polygonList[i + 1], start + polygonList[i + 2],
					start + polygonList[i + 3]);
			}
		}
	}
	else {	// triangles
#ifdef DEBUG
		cout << "Writing " << npolys << " 3-tiles\n";
#endif
		for (i = 0; i < 3*npolys; i += 3)
			fprintf((FILE *) fp, "%d %d -%d\n", start + polygonList[i + 0],
				start + polygonList[i + 1], start + polygonList[i + 2]);
	}

	if (origin != NULL && bound != NULL) {
		// An extra line is added at the end, when using world coordinates.
		// This is an extension of the original BYU format necessary to
		// make it possible for Pablo to read the BYU file and place the
		// resulting tiles correctly in model space.  This writes the
		// positions of the extreme corners of model space in world
		// coordinates.
		fprintf((FILE *) fp, "1 %lf %lf %lf %lf %lf %lf\n", origin->getX(), origin->getY(),
			origin->getZ(), bound->getX(), bound->getY(), bound->getZ());
	}
	else
		fprintf((FILE *) fp, "0\n");

	if (! fptied) {
		fclose((FILE *) fp);
		fp = NULL;
	}
	return true;
}

/*  Low-level output function for BYU files.  A BYU file can describe multiple figures,
    called parts.  Each part must have a list of vertices and a list of polygons; each
	polygon is represented as a set of indexes into the list of vertices.  These items
	should be provided to this function as follows:
		numParts		The number of parts to be output.
		numVertices		An array of length numParts containing the numbers of
							vertices of each part.  There are 3 values per vertex.
		numPolygons		An array of length numParts containing the numbers of
							polygons of each part.  The BYU format can support any
							number of points per polygon, but this function accepts
							only 3 or 4 vertices per polygon.  Quad-tiles may also
							be split in half to produce triangles (see below).
		vertexLists		An array of length numParts containing arrays of vertices
							of each part.  Each array should have a length that is
							three times that given in the corresponding entry in
							numVertices, since there are 3 coordinates per vertex.
		polygonLists	An array of length numParts containing arrays of polygons
							of each part.  Polygons are represented as indexes into
							the vertexLists array.  Each array should have a length
							3 or 4 times that given in the corresponding entry of
							numPolygons, as there are 3 or 4 vertices per input
							polygon.

	This function can write its output as 4-tiles, 3-tiles, or can split 4-tiles to
	produce pairs of 3-vertex tiles.  The typeIn argument should be used to indicate
	the current form of the polygon data.  If typeIn is Quad, then split may be used
	to indicate that the 4-tiles are to be split.

    If an origin and bound are provided, they will be used to specify an extension to
	the regular BYU format.  The tile vertices will be assumed to already be in world
	coordinates; they will not be converted.

	This function is much the same as write1(), differing only in  the structure of
	arguments vertexLists and polygonLists.
*/
bool BYU::write2(const char * filename, BYUPoly_t typeIn, bool split, int numParts,
	int * numVertices, int * numPolygons, float ** vertexLists, int ** polygonLists,
	Vector3D * origin, Vector3D * bound)
{
	int i;
	int partNo;
	int nverts, npolys;

#ifdef DEBUG
	cout << "BYU::write2()" << '\n';
#endif
	if (typeIn == NoTiles)
		return false;
	if (typeIn == Tris && split)
		return false;

	if (fp == NULL && filename == NULL)
		return false;

    if (fp == NULL)
		fp = fopen(filename, "wb");
    if (fp == NULL)
        return false;

	nverts = 0;
	npolys = 0;
	for (i = 0; i < numParts; i++) {
		nverts += numVertices[i];
		npolys += numPolygons[i];
	}

	// Output the header lines
	// First, numbers of: figures, vertices, quads/triangles, and edges.
	if (typeIn == Quads) {
		if (split)
			fprintf((FILE *) fp, "%7d %7d %7d %7d\n", numParts, nverts, 2*npolys, 6*npolys);
		else
			fprintf((FILE *) fp, "%7d %7d %7d %7d\n", numParts, nverts, npolys, 4*npolys);
	}
	else
		fprintf((FILE *) fp, "%7d %7d %7d %7d\n", numParts, nverts, npolys, 3*npolys);
	// Second, the sequence of beginning and ending polygon numbers for each part
	i = 1;
	for (partNo = 0; partNo < numParts; partNo++) {
		// The BYU format requires 7 bytes per value
		if (typeIn == Quads && split) {
			// Splitting quads
			fprintf((FILE *) fp, "%7d %7d ", i, 2*numPolygons[partNo] + i - 1);
			i += 2*numPolygons[partNo];
		}
		else {	// Tris or quads
			fprintf((FILE *) fp, "%7d %7d ", i, numPolygons[partNo] + i - 1);
			i += numPolygons[partNo];
		}
	}
	fprintf((FILE *) fp, "\n");

	// Output the vertex lines (list of points)
	for (partNo = 0; partNo < numParts; partNo++) {
		for(i = 0; i < 3*numVertices[partNo]; i += 3) {
			fprintf((FILE *) fp, "%f %f %f\n", vertexLists[partNo][i],
				vertexLists[partNo][i + 1], vertexLists[partNo][i + 2]);
		}
	}

	// Output the polygon lines (indexes into list of vertices)
	// Note that the BYU format requires indexing counted from 1, instead of 0
	int start = 1;	// The BYU format requires 1-origin indexing
	if (typeIn == Quads) {	// quads
		if (split) {
			for (partNo = 0; partNo < numParts; partNo++) {
#ifdef DEBUG
				cout << "Writing " << 2*numPolygons[partNo] << " 3-tiles for figure "
					<< partNo << '\n';
#endif
				for (i = 0; i < 4*numPolygons[partNo]; i += 4) {
					fprintf((FILE *) fp, "%d %d -%d\n",
						start + polygonLists[partNo][i + 3],
						start + polygonLists[partNo][i + 0],
						start + polygonLists[partNo][i + 1]);
					fprintf((FILE *) fp, "%d %d -%d\n",
						start + polygonLists[partNo][i + 1],
						start + polygonLists[partNo][i + 2],
						start + polygonLists[partNo][i + 3]);
				}
				start += numVertices[partNo];
			}
		}
		else {
			for (partNo = 0; partNo < numParts; partNo++) {
#ifdef DEBUG
				cout << "Writing " << numPolygons[partNo] << " 4-tiles for figure "
					<< partNo << '\n';
#endif
//				cout << "Polygons for part " << partNo << '\n';
				for (i = 0; i < 4*numPolygons[partNo]; i += 4) {
//					cout << polygonLists[partNo][i + 0] << ", " << polygonLists[partNo][i + 1] << ", "
//						<< polygonLists[partNo][i + 2] << ", " << polygonLists[partNo][i + 3] << '\n';
					fprintf((FILE *) fp, "%d %d %d -%d\n", start + polygonLists[partNo][i + 0],
						start + polygonLists[partNo][i + 1], start + polygonLists[partNo][i + 2],
						start + polygonLists[partNo][i + 3]);
				}
				start += numVertices[partNo];
			}
		}
	}
	else {	// triangles
		for (partNo = 0; partNo < numParts; partNo++) {
#ifdef DEBUG
			cout << "Writing " << numPolygons[partNo] << " 3-tiles for figure "
				<< partNo << '\n';
#endif
			for (i = 0; i < 3*numPolygons[partNo]; i += 3) {
				fprintf((FILE *) fp, "%d %d -%d\n",
					start + polygonLists[partNo][i + 0], start + polygonLists[partNo][i + 1],
					start + polygonLists[partNo][i + 2]);
			}
			start += numVertices[partNo];
		}
	}

	if (origin != NULL && bound != NULL) {
		// An extra line is added at the end, when using world coordinates.
		// This is an extension of the original BYU format necessary to
		// make it possible for Pablo to read the BYU file and place the
		// resulting tiles correctly in model space.  This writes the
		// positions of the extreme corners of model space in world
		// coordinates.
		fprintf((FILE *) fp, "1 %lf %lf %lf %lf %lf %lf\n", origin->getX(), origin->getY(),
			origin->getZ(), bound->getX(), bound->getY(), bound->getZ());
	}
	else
		fprintf((FILE *) fp, "0\n");

	if (! fptied) {
		fclose((FILE *) fp);
		fp = NULL;
	}
	return true;
}

BYU::BYUPoly_t BYU::type() const
{
	if (edges_per_tile == 3)
		return Tris;
	if (edges_per_tile == 4)
		return Quads;
	return NoTiles;
}

void BYU::clear()
{
	if (coords) {
		delete [] polygons;
		polygons = NULL;
		delete [] coords;
		coords = NULL;
		tileAllocation = 0;
		delete [] vertices;
		vertices = NULL;
		delete [] tileCounts;
		tileCounts = NULL;
	}
	numTiles = 0;
	nfigs = 0;
	edges_per_tile = 0;
}

// ------------------------ Start of Marching Cubes Code ------------------------

inline GreyValue VOX(int x, int y, int z, Image3D * image)
{
	return image->getVoxelValue(x, y, z);
}

/*	Visualize a binary image (0 outside, 1 inside) as a blocky, quad-tiled surface.

	Running this discards any existing tiles; the result is considered to be a
	single set of tiles stored in coords, even if there are actually multiple
	figures.

	This functions sets variables edges_per_tile, nfigs (to 1), tileCounts, coords,
	numTiles, leftBottomBack, and rightTopFront.

	[LATER] add a mask to indicate which bits to detect on, for bitStack images.
*/
bool BYU::simpleMarchingCubes(Image3D * image, bool asQuads)
{
	float ux, uy, uz; 	// unit-cube coord of (x,y,z)

	if (! image)
		return false;

	if (! asQuads) {
		cout << "\nsimpleMarchingCubes: generation of triangular tiles not implemented\n";
		return false;
	}

	clear();	// Discard any existing tiles
	quads = asQuads;
	if (quads)
		edges_per_tile = 4;
	else
		edges_per_tile = 3;		// This is never reached (see above)
	// The tiles produced will be quads, because 4 is hard coded everywhere below

	// foreach voxel V:
	// - foreach dimension in [x, y, z]
	//   - look "back" (towards the origin) in the dim to neighbor N;
	//     if (intensity(N) != intensity(V), draw a tile between the 2 voxels

	int tiles = 0;

	// m is the mesh size: the size of a voxel in unit cube space
	int xdim = image->getXDim();
	int ydim = image->getYDim();
	int zdim = image->getZDim();
	int maxDim = xdim;
	if (ydim > maxDim)
		maxDim = ydim;
	if (zdim > maxDim)
		maxDim = zdim;

	// Convert world-coordinate voxel spacing to model spacing to
	// handle anisotopic voxels
	Vector3D spacing(image->getXSpacing(), image->getYSpacing(),
		image->getZSpacing());
	spacing += image->getModelOrigin();
	image->worldToModelCoordinates(spacing);
	// Spacing now contains the voxel sizes in model coords

	// If (current intensity != neighbor intensity) draw a tile between current
	// voxel and the neighbor.
	// [LATER] this does not allow for the binary to hit the edge of the image
	// Modified by AGG to threshold 0 (versus nonzero), in case a gray image
	// is used, which would otherwise yield too many tiles.
	for (int z=1; z < zdim; z++)
	{
		uz = float(z)*spacing.getZ();
		for (int y=1; y < ydim; y++)
		{
			uy = float(y)*spacing.getY();
			for (int x=1; x < xdim; x++)
			{
				ux = float(x)*spacing.getX();
				int gotOne = 0;
				GreyValue me = VOX(x, y, z, image);
				if (me)
				{
					// I'm inside - any neighbors outside?
					if (0 == VOX(x - 1, y, z, image))
					{
						tiles++; gotOne = 1;
						addTile(calcQuadTile(ux, uy, uz, 1, &spacing), 4);
					}
					if (0 == VOX(x, y - 1, z, image))
					{
						tiles++; gotOne = 1;
						addTile(calcQuadTile(ux, uy, uz, 2, &spacing), 4);
					}
					if (0 == VOX(x, y, z - 1, image))
					{
						tiles++; gotOne = 1;
						addTile(calcQuadTile(ux, uy, uz, 3, &spacing), 4);
					}
				}

				else
				{
					// I'm outside - any neighbors inside?
					if (0 != VOX(x - 1, y, z, image))
					{
						tiles++; gotOne = 1;
						addTile(calcQuadTile(ux, uy, uz, -1, &spacing), 4);
					}
					if (0 != VOX(x, y - 1, z, image))
					{
						tiles++; gotOne = 1;
						addTile(calcQuadTile(ux, uy, uz, -2, &spacing), 4);
					}
					if (0 != VOX(x, y, z - 1, image))
					{
						tiles++; gotOne = 1;
						addTile(calcQuadTile(ux, uy, uz, -3, &spacing), 4);
					}
				}

				if (gotOne && !(tiles % 1000))
					cout << '.' << flush;
			}
		}
	}

	tileCounts = new int[1];
	tileCounts[0] = tiles;
	cout << " " << tiles << " tiles\n";
	nfigs = 1;
	generated = true;
	return true;
}

/*
void BYU::print(int tileCount)
{
	printf("TILES: %d %s tiles\n", numTiles,
		(quads ? "quad" : "tri"));
	printf(" leftBottomBack: (%f %f %f)\n",
		leftBottomBack.getX(),
		leftBottomBack.getY(),
		leftBottomBack.getZ());
	printf(" rightTopFront: (%f %f %f)\n",
		rightTopFront.getX(),
		rightTopFront.getY(),
		rightTopFront.getZ());

	int ept = quads ? 4 : 3;		// edges_per_tile
	for (int t=0; t < ept*tileCount; t += ept)
	{
		if (quads)
			printf(" tile %d: (%f %f %f) - (%f %f %f) - (%f %f %f)\n",
			t/ept,
			coords[t+0].getX(), coords[t+0].getY(), coords[t+0].getZ(), 
			coords[t+1].getX(), coords[t+1].getY(), coords[t+1].getZ(), 
			coords[t+2].getX(), coords[t+2].getY(), coords[t+2].getZ(), 
			coords[t+3].getX(), coords[t+3].getY(), coords[t+3].getZ());
		else
			printf(" tile %d: (%f %f %f) - (%f %f %f) - (%f %f %f)\n",
			t/ept,
			coords[t+0].getX(), coords[t+0].getY(), coords[t+0].getZ(), 
			coords[t+1].getX(), coords[t+1].getY(), coords[t+1].getZ(), 
			coords[t+2].getX(), coords[t+2].getY(), coords[t+2].getZ());
	}
}
*/

// make sure there's at least room for wantNumTiles tiles;
//  copy any existing old points to this new space.
// caller updates numTiles AFTER calling this
void BYU::allocTiles(int wantNumTiles, int edgesPerTile)
{
	Vector3D * newCoords;

	// numTiles is the existing tile count; wantNumTiles is the desired count.
	// tileAllocation is the number allocated in coords[].
	if (wantNumTiles <= tileAllocation)
		return;

	if (tileAllocation > 0) {
		newCoords = new Vector3D[2*tileAllocation*edgesPerTile];
		tileAllocation += tileAllocation;
	}
	else {
		newCoords = new Vector3D[MIN_TILESET_SIZE*edgesPerTile];
		tileAllocation = MIN_TILESET_SIZE;
	}

#ifdef DEBUG
	std::cout << tileAllocation << '\n';
#endif
	if (numTiles)	// anything to copy?
	{
		memcpy(newCoords, coords,
			numTiles * edgesPerTile * sizeof(Vector3D));
		delete [] coords;
	}
	coords = newCoords;
}

// Extend tile set by a 4-vertex tile.
// This assumes the type of type (tri or quad) matches edgesPerTile.
void BYU::addTile(Vector3D *newTile, int edgesPerTile)
{
    if (coords == NULL)
	{
		// First call of addTile()
		leftBottomBack = rightTopFront = newTile[0];
/*
		leftBottomBack.set(newTile[0].getX(), newTile[0].getY(),
			newTile[0].getZ());
		rightTopFront. set(newTile[0].getX(), newTile[0].getY(),
			newTile[0].getZ());
*/
	}
	allocTiles(numTiles + 1, edgesPerTile);

	for (int pt=0; pt < edgesPerTile; pt++)
	{
		coords[numTiles*edgesPerTile + pt] = newTile[pt];
		leftBottomBack.minPoint(newTile[pt]);
		rightTopFront.maxPoint(newTile[pt]);
/*
		float x = newTile[pt].getX();
		float y = newTile[pt].getY();
		float z = newTile[pt].getZ();

		coords[numTiles*edgesPerTile + pt].set(x, y, z);

		// Extend bounding box to cover new vertices
		if(x < leftBottomBack.getX())
			leftBottomBack.setX(x);
		if(y < leftBottomBack.getY())
			leftBottomBack.setY(y);
		if(z < leftBottomBack.getZ())
			leftBottomBack.setZ(z);

		if(x > rightTopFront.getX())
			rightTopFront.setX(x);
		if(y > rightTopFront.getY())
			rightTopFront.setY(y);
		if(z > rightTopFront.getZ())
			rightTopFront.setZ(z);
*/
	}	
	numTiles++;
}

// Figure out a quad tile between point (x,y,z) and it's neighbor in the
// 'dir' direction (1=x, 2=y, 3=z). If dir is negative, invert tile normal.
// m is the size of the mesh (a value for each of x,y,z).
Vector3D *BYU::calcQuadTile(float x, float y, float z, int dir, Vector3D * m)
{
   static Vector3D q[4];    // a quad tile to return

   // These "points" are halfway between (x,y,z) and it's neighbors
   // in the +x, +y, or +z direction.  Combining subsets of this point
   // and the (x,y,z) point results in all possible corners of this voxel
   // in the "positive axis" direction.
   float xm = x + m->getX();
   float ym = y + m->getY();
   float zm = z + m->getZ();

   int absDir = (dir > 0) ? dir : -dir;
   switch (absDir)
   {
   case 1:	// y-z tile
			q[2].setX(x );
			q[2].setY(y );
			q[2].setZ(z );

			q[1].setX(x );
			q[1].setY(ym);
			q[1].setZ(z );
			q[0].setX(x );
			q[0].setY(ym);
			q[0].setZ(zm);
			q[3].setX(x );
			q[3].setY(y );
			q[3].setZ(zm);
			break;
	case 2:	// x-z tile
			q[0].setX(x );
			q[0].setY(y );
			q[0].setZ(z );

			q[1].setX(xm);
			q[1].setY(y );
			q[1].setZ(z );
			q[2].setX(xm);
			q[2].setY(y );
			q[2].setZ(zm);
			q[3].setX(x );
			q[3].setY(y );
			q[3].setZ(zm);
			break;
	case 3:	// x-y tile
			q[2].setX(x );
			q[2].setY(y );
			q[2].setZ(z );

			q[1].setX(xm);
			q[1].setY(y );
			q[1].setZ(z );
			q[0].setX(xm);
			q[0].setY(ym);
			q[0].setZ(z );
			q[3].setX(x );
			q[3].setY(ym);
			q[3].setZ(z );
			break;
   }

   if (dir < 0)		// Invert tile by swapping vertices 0 and 2
   {
       Vector3D temp(q[0]);
       q[0] = q[2];
       q[2] = temp;
   }

   return q;
}

// ------------------------ End of Marching Cubes Code ------------------------



