#ifdef WIN32
# include <windows.h>
#endif
#ifdef __APPLE__
#include <OpenGL/gl.h>
#else
#include <GL/gl.h>
#endif
#include <stdio.h>
#include <sys/stat.h>
#include <iostream>
#include "Image3D.h"
#include "BYU.h"
#include "TileSet.h"


#ifndef S_ISREG
#define S_ISREG(m)        (((S_IFREG & m) != 0) ? true : false)
#endif


#include <string.h>

//#define DEBUG


using namespace std;

TileSet::TileSet()
{
	numFigures = 0;
    numTiles = 0;
	byu = NULL;
	have_byu = false;
	coords = NULL;
	tileCounts = NULL;
	local = false;	// Allocation of array coords
}

TileSet::TileSet(const char * filename, bool split, const Image3D* image)
{
	numFigures = 0;
    numTiles = 0;
	byu = NULL;
	have_byu = false;
	coords = NULL;
	tileCounts = NULL;
	local = false;	// Allocation of array coords

    (void) loadFile(filename, split, image);
}

TileSet::~TileSet()
{
	BYU * byuPtr;

    if (byu != NULL) {
		byuPtr = (BYU *) byu;
        delete byuPtr;
	}
	if (local && coords != NULL)
		delete [] coords;
	// colors should not be deleted
	// tileCounts should not be deleted
}

int TileSet::type(const char * filename, bool & isBYU) const
{
    struct stat buf;

	isBYU = false;
	if (filename == NULL)
		return 0;

	if (stat(filename, &buf) != 0 || S_ISREG(buf.st_mode) == false) {
		cout << "Not a valid file.\n";
		return 0;
	}

    // Determine file format from name or contents
    if (0 == strcmp(filename + strlen(filename) - 4, ".til"))
		return 3;
    if (0 == strcmp(filename + strlen(filename) - 4, ".byu")) {
		BYU b;

		isBYU = true;
		BYU::BYUPoly_t t = b.type(filename);
		if (t == BYU::Quads)
			return 4;
		else if (t == BYU::Tris)
			return 3;
		else
			return 0;
	}
	else {
		FILE * fp;
		int nf, np, nt, ne;

		fp = fopen(filename, "r");
		if (fp == NULL)
			return 0;

        if (4 == fscanf(fp, "%d %d %d %d\n", &nf, &np, &nt, &ne))
			isBYU = true;
		else
			isBYU = false;

	    fclose(fp);
		return ne/nt;
	}
}

// Center of gravity is the average vertex position
Vector3D TileSet::getCOG()
{
	Vector3D cog;
	Vector3D origin;
	Vector3D opposite;

	if (coords == NULL)		// no tiles!?!
		return Vector3D(.5, .5, .5);

	getBoundingBox(origin, opposite);
	return (origin + opposite) / 2.0;
}


// set values used by getCOG() by (re)examining vertices of spec'd figures
void TileSet::calcCOG(bool figs[6])
{
	int lo=0, hi;
	int ptsPerTile = pointsPerTile();
	int firstTime = 1;
	
	for (int fig=0; fig < numFigures; fig++)
	{
		hi = lo + tileCounts[fig] - 1;
		if (figs[fig])
		{
			for (int t = lo; t <= hi; t++)
				for (int pt = 0; pt < ptsPerTile; pt++)
				{
					int v = t*ptsPerTile + pt;
					if (firstTime) {
						firstTime = 0;
						leftBottomBack = rightTopFront = coords[v];
					}
					else {
						leftBottomBack.minPoint(coords[v]);
						rightTopFront.maxPoint(coords[v]);
					}
				}
		}
		lo = hi + 1;
	}
}
	
// set values used by getCOG() by (re)examining ALL vertices
void TileSet::calcCOG()
{
	int maxTile = numTiles * pointsPerTile();
	if (!maxTile)
		return;

	leftBottomBack = rightTopFront = coords[0];
	for (int t=0; t < maxTile; t++)
	{
		leftBottomBack.minPoint(coords[t]);
		rightTopFront.maxPoint(coords[t]);
	}
}


bool TileSet::loadFile(const char * filename, bool split, const Image3D* image)
{
    struct stat buf;
    int edges_per_tile;
    int v1, v2, v3;
    float x, y, z;
	int i;
	BYU * byuPtr;
	FILE * fp;

	if (filename == NULL)
		return false;

	if (stat(filename, &buf) != 0 || S_ISREG(buf.st_mode) == false) {
		cout << "Not a valid file.\n";
		return false;
	}

	fp = fopen(filename, "r");
	if (fp == NULL)
		return false;

    // Determine file format from name or contents
    if (0 == strcmp(filename + strlen(filename) - 4, ".byu"))
        have_byu = true;
    else if (0 == strcmp(filename + strlen(filename) - 4, ".til"))
        have_byu = false;
	else {
        if (4 == fscanf(fp, "%d %d %d %d\n", &i, &v1, &v2, &v3))
			have_byu = true;
		else
			have_byu = false;
		fseek(fp, 0L, SEEK_SET);	// rewind
	}

    if (have_byu) {
#ifdef DEBUG
		cout << "Loading BYU tiles\n";
#endif
		// Use class BYU for file input
		byuPtr = (BYU *) byu;
		if (byuPtr == NULL) {
			byuPtr = new BYU;
			byu = byuPtr;
		}
		byuPtr->tie(fp);

		// First find the number of tiles
		if (! byuPtr->read(NULL, split, image))	// This clears any previous BYU data
			numTiles = 0;
		else
			numTiles = byuPtr->tileListLength();
    }        
    else {
#ifdef DEBUG
		cout << "Loading UNC tiles\n";
#endif
		edges_per_tile = 3;

        // Read as a UNC-format tile file:
        // count # of points to figure # of tiles
        float x, y, z;
        int ret;
        while(fp != NULL)
        {
            ret = fscanf(fp, "%f %f %f", &x, &y, &z);
            if(ret == 0 || ret == EOF)
                    break;

            ret = fscanf(fp, "%f %f %f", &x, &y, &z);
            if(ret == 0 || ret == EOF)
                    break;

            ret = fscanf(fp, "%f %f %f", &x, &y, &z);
            if(ret == 0 || ret == EOF)
                    break;

            numTiles++;
        }
		tc = numTiles;
    }

    if (numTiles == 0) {
		numFigures = 0;
		fclose(fp);
        cout << "File contained no tiles.\n";
        return false;
	}

    if (have_byu) {
		coords = (Vector3D *) byuPtr->tileList();
		local = false;
		quads = (byuPtr->edgesPerTile() == 4);
		byuPtr->getBoundingBox(leftBottomBack, rightTopFront);
		numFigures = byuPtr->figureCount();
		tileCounts = byuPtr->tilesPerFigure();
		byuPtr->untie();
    }
    else {
		coords = new Vector3D[numTiles * edges_per_tile];
		local = true;
		quads = false;

		fseek(fp, 0L, SEEK_SET);	// rewind
        for (i = 0; i < numTiles * 3; i++)
        {
            fscanf(fp, "%f %f %f", &x, &y, &z);

            coords[i].set(x, y, z);
            if (i == 0) {
                leftBottomBack = rightTopFront = coords[0];
            }
            else {
				leftBottomBack.minPoint(coords[i]);
				rightTopFront.maxPoint(coords[i]);
            }

        }
 		numFigures = 1;
		tileCounts = &tc;
   }

    fclose(fp);
    return true;
}

bool TileSet::loadAppendFile(const char * filename, bool split, const Image3D* image)
{
	BYU * byuPtr;
	BYU * byuPtr_2;

	byuPtr = (BYU *) byu;
	if (byuPtr == NULL)
		return false;

	TileSet ts(filename, split, image);
	if (ts.tileCount() <= 0)
		return false;

	byuPtr_2 = (BYU *) ts.byu;
	if (byuPtr_2 == NULL) {
		cout << "Only BYU files can be merged" << endl;
		return false;
	}

	if (! byuPtr->merge(*byuPtr_2))
		return false;

	numTiles = byuPtr->tileListLength();
	coords = (Vector3D *) byuPtr->tileList();
	byuPtr->getBoundingBox(leftBottomBack, rightTopFront);
	numFigures = byuPtr->figureCount();
	tileCounts = byuPtr->tilesPerFigure();

	return true;
}

bool TileSet::saveFile(const char * filename, bool asQuads, Image3D * image)
{
	BYU * byuPtr;
	BYU::BYUPoly_t type;
	BYU::BYUPoly_t typeOut;

	if (filename == NULL)
		return false;

	if (numTiles <= 0) {
		cout << "Warning: no tiles to save" << endl;
		return false;
	}

	byuPtr = (BYU *) byu;
	if (byuPtr == NULL) {
		byuPtr = new BYU;
		byu = byuPtr;
	}

	type = (quads ? BYU::Quads : BYU::Tris);
	typeOut = (asQuads ? BYU::Quads : BYU::Tris);
	if (type == BYU::Tris && typeOut == BYU::Quads) {
		cout << "Warning: cannot convert triangles to quads; writing BYU data as triangles" << endl;
		typeOut = BYU::Tris;
	}
	if (image == NULL)
		return byuPtr->write(filename, typeOut, numFigures, &numTiles, coords,
			type, &leftBottomBack, &rightTopFront);
	else
		return byuPtr->write(filename, typeOut, numFigures, &numTiles, coords,
			type, &leftBottomBack, &rightTopFront, image);
}

bool TileSet::simpleMarchingCubes(Image3D * image, bool asQuads)
{
	BYU * byuPtr;

	byuPtr = (BYU *) byu;
	if (byuPtr == NULL) {
		byuPtr = new BYU;
		byu = byuPtr;
	}

	if (! byuPtr->simpleMarchingCubes(image, asQuads)) {
		coords = NULL;
		numTiles = 0;
		numFigures = 0;
		return false;
	}
	else {
		coords = (Vector3D *) byuPtr->tileList();
		numTiles = byuPtr->tileListLength();
		local = false;
		quads = asQuads;
		tileCounts = byuPtr->tilesPerFigure();
		numFigures = byuPtr->figureCount();
	}

	return true;
}

/////////////////////////////////////////////////////
// translate, scale and rotate tiles, per figure
/////////////////////////////////////////////////////

void TileSet::translateBy(const Vector3D & vTrans, bool figs[6])
{
	int lo=0, hi;
	int ptsPerTile = pointsPerTile();

	for (int fig=0; fig < numFigures; fig++)
	{
		hi = lo + tileCounts[fig] - 1;
		if (figs[fig])
		{
			for (int t = lo; t <= hi; t++)
				for (int pt = 0; pt < ptsPerTile; pt++)
				  coords[t*ptsPerTile + pt] += vTrans;
		}
		lo = hi + 1;
	}

	calcCOG(figs);
}

// NOT in class
static void scaleBy(double scaleFactor, const Vector3D & center, Vector3D &x)
{
	x = ((x - center) * scaleFactor) + center;
}

void TileSet::scaleBy(double scaleFactor, bool figs[6])
{
	scaleBy(scaleFactor, getCOG(), figs);
}


void TileSet::scaleBy(double scaleFactor, const Vector3D & center, bool figs[6])
{
	int lo=0, hi;
	int ptsPerTile = pointsPerTile();

	for (int fig=0; fig < numFigures; fig++)
	{
		hi = lo + tileCounts[fig] - 1;
		if (figs[fig])
		{
			for (int t = lo; t <= hi; t++)
				for (int pt = 0; pt < ptsPerTile; pt++)
					::scaleBy(scaleFactor, center, coords[t*ptsPerTile + pt]);
		}
		lo = hi + 1;
	}

	calcCOG(figs);
}

// NOT in class: rotate vector x by q about center
static void rotateBy(const Quat& q, const Vector3D& center, Vector3D& x)
{
	x -= center;
	q.rotateVector(x);
	x += center;
}

// calc COG of selected figure-tiles and rotate tiles about that
void TileSet::rotateBy(const Quat &q, bool figs[6])
{
	calcCOG(figs);
	rotateBy(q, getCOG(), figs);
}

// rotate selected figure-tiles about center, then update COG
void TileSet::rotateBy(const Quat &q, const Vector3D & center, bool figs[6])
{
	int lo=0, hi;
	int ptsPerTile = pointsPerTile();

	for (int fig=0; fig < numFigures; fig++)
	{
		hi = lo + tileCounts[fig] - 1;
		if (figs[fig])
		{
			for (int t = lo; t <= hi; t++)
				for (int pt = 0; pt < ptsPerTile; pt++)
					::rotateBy(q, center, coords[t*ptsPerTile + pt]);
		}
		lo = hi + 1;
	}

	calcCOG(figs);
}

// rotate each selected tile-figure about it's own COG
void TileSet::pirouetteBy(const Quat &q, bool figs[6])
{
	bool onePerFig[6];		// copy of figs[] with only 1 fig "on"
	for (int fig=0; fig < numFigures; fig++)
	{
		if (figs[fig])
		{
			for (int theFig=0; theFig < numFigures; theFig++)
				onePerFig[theFig] = (fig == theFig);
			rotateBy(q, onePerFig);
		}
	}
}

void TileSet::setGenerated()
{
	if (byu) { 
		BYU *byuPtr = (BYU *)byu;
		byuPtr->setGenerated();
	}
}


