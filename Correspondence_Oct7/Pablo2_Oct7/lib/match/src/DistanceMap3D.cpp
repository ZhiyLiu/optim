//#include <sys/timeb.h>
//#include <unistd.h>

#include <list>
#include <iostream>
#include <math.h>
#include "Danielsson.h"
#include "M3DFigure.h"
#include "DistanceMap3D.h"
#ifdef BINARY
#include "snapshot.h"
#endif
#include "Tuning.h"

#define D_XFERLIST
#define D_POINTLIST_SERVERB
#define D_SUBDIVSURF
#include "Shapedepend.h"


//#define DEBUG


const Vector3D DEFAULT_BBOX_ORIGIN = Vector3D(0.0, 0.0, 0.0);
const Vector3D DEFAULT_BBOX_SIZE = Vector3D(1.0, 1.0, 1.0);
const int DEFAULT_MAP_SIZE[3] = {128, 128, 128};



using namespace std;



std::ostream & operator<< (std::ostream & ostr, BlendedFiguralCoordinate & fc_ref) {
	ostr << "Figural coordinate ("
		<< fc_ref.u << ", "
		<< fc_ref.v << ", "
		<< fc_ref.t << "):    distance = "
		<< fc_ref.distance
		<< "\n  tiles: " << fc_ref.tileId
		<< ' ' << (int) fc_ref.tileId2
		<< "    radii: " << fc_ref.r
		<< ' ' << fc_ref.r2
		;
	if (fc_ref.inBlendRegion)
		ostr << "    in blend region";
	ostr << endl;
	return ostr;
}

std::ostream & operator<< (std::ostream & ostr, FiguralCoordinate & fc_ref) {
	ostr << "Figural coordinate ("
		<< fc_ref.u << ", "
		<< fc_ref.v << ", "
		<< fc_ref.t << "):    distance = "
		<< fc_ref.distance
		<< "\n  tiles: " << fc_ref.tileId
		//<< ' ' << (int) fc_ref.tileId2
		<< "    radii: " << fc_ref.r
		//<< ' ' << fc_ref.r2
		;
	//if (fc_ref.inBlendRegion)
	//	ostr << "    in blend region";
	ostr << endl;
	return ostr;
}


float DistanceMap3D::bilerp_in_z(float x, float y, int z)
{
    unsigned int ix, iy;
    float   dx, dy;

    ix = (unsigned int) x;
    iy = (unsigned int) y;
    float a = getDistance(ix, iy, z);
    float b = getDistance((ix + 1), iy, z);
    float c = getDistance(ix, (iy + 1), z);
    float d = getDistance((ix + 1), (iy + 1), z);
    dx = x - (float) ix;
    dy = y - (float) iy;

    float i1 = (1.0 - dx)*a + dx*b;
    float i2 = (1.0 - dx)*c + dx*d;
    return (1.0 - dy)*i1 + dy*i2;
}


float DistanceMap3D::trilerp(float x, float y, float z)
{
    unsigned int    ix, iy, iz;
    float   a, b;
    float   dz;

    ix = (unsigned int) x;
    iy = (unsigned int) y;
    iz = (unsigned int) z;

    a = bilerp_in_z(x, y, iz);
    b = bilerp_in_z(x, y, iz + 1);
    dz = z - (float) iz;

    return (1.0 - dz)*a + dz*b;
}


DistanceMap3D::DistanceMap3D()
{
    figuralCoordinateMap = NULL;
    object = NULL;
    mapSize[0] = 0;
    size = 0;

    min_mapsize = MIN_MAP_SIZE;
    max_mapsize = MAX_MAP_SIZE;
	r_normalized = false;
	//cout << "FiguralCoord size = " << sizeof(FiguralCoordinate) << endl;
	//cout << "BlendedFiguralCoord size = " << sizeof(BlendedFiguralCoordinate) << endl;
}


DistanceMap3D::DistanceMap3D(int mapSizeMinimum, int mapSizeMaximum)
{
    figuralCoordinateMap = NULL;
    object = NULL;
    mapSize[0] = 0;
    size = 0;

    min_mapsize = mapSizeMinimum;
    max_mapsize = mapSizeMaximum;
	r_normalized = false;
}


DistanceMap3D::DistanceMap3D(M3DObject * _object, int granularity, double cutoff,
                             int mapSizeMinimum, int mapSizeMaximum)
{
    figuralCoordinateMap = NULL;
    mapSize[0] = 0;
    size = 0;

    object = new M3DObject(*_object);
    numMedialSubdivisions = granularity;
    factor = cutoff;
    min_mapsize = mapSizeMinimum;
    max_mapsize = mapSizeMaximum;
    bboxSize[0] = DEFAULT_BBOX_SIZE.getX();
    bboxSize[1] = DEFAULT_BBOX_SIZE.getY();
    bboxSize[2] = DEFAULT_BBOX_SIZE.getZ();
	r_normalized = false;
}


void DistanceMap3D::cleanup()
{
    if (figuralCoordinateMap != NULL) {
        delete [] figuralCoordinateMap;
        figuralCoordinateMap = NULL;
        size = 0;
    }

    if (object != NULL) {
        delete object;
        object = NULL;
    }

	for (int i = 0; i < tileLists.size(); i++) {
		delete tileLists[i];
	    tileLists[i] = NULL;
	}
	tileLists.clear();
	figureIdList.clear();
}


DistanceMap3D::~DistanceMap3D()
{
    cleanup();
}


void DistanceMap3D::initialize(M3DObject * _object, int granularity, double cutoff,
                               int mapSizeMinimum, int mapSizeMaximum)
{
    if (_object == NULL) {
        cout << "DistanceMap3D object not initialized\n";
        return;
    }

    cleanup();

    object = new M3DObject(*_object);
    numMedialSubdivisions = granularity;
    factor = cutoff;
    min_mapsize = mapSizeMinimum;
    max_mapsize = mapSizeMaximum;

    initialize();
}


void DistanceMap3D::initialize()
{
    BoundingBox bbox;

    if (figuralCoordinateMap != NULL) {
        delete [] figuralCoordinateMap;
        figuralCoordinateMap = NULL;
        size = 0;
    }
    (void) defineBBox(object, bbox, mapSize, numMedialSubdivisions, factor);

    bboxOrigin = bbox.origin;
    bboxSize[0] = bbox.size.getX();
    bboxSize[1] = bbox.size.getY();
    bboxSize[2] = bbox.size.getZ();
    bboxVoxelsize = bbox.voxelsize;

#ifdef DEBUG
    cout << "Mapsize: " << mapSize[0] << ", " << mapSize[1] << ", "  << mapSize[2]
        <<  ", voxel size " << bboxVoxelsize << endl;
#endif
}


BoundingBox DistanceMap3D::getBoundingBox() {
    BoundingBox bbox;

    bbox.origin = bboxOrigin;
    bbox.size.setX(bboxSize[0]);
    bbox.size.setY(bboxSize[1]);
	bbox.size.setZ(bboxSize[2]);
    bbox.voxelsize = bboxVoxelsize;
	return bbox;
}


// Determine a bounding box for an object within the limits of the maximum and
// minimum bounding box sizes specified to the constructor.
bool DistanceMap3D::defineBBox(M3DObject * _object, BoundingBox & _bbox,
                               int _mapSize[3], int granularity, double cutoff)
{
    int figureCount, i, k;
    int numBPoints;
    Bpoint * bpoints;
    double x, y, z;
    double minX, minY, minZ, minR;
    double maxX, maxY, maxZ;

    if (_object == NULL) {
        cout << "DistanceMap3D object not initialized\n";
        return false;
    }

    minX = minY = minZ = minR = 1000000.0;
    maxX = maxY = maxZ = 0.0;

    // Determine the size of a bounding box around the object
    M3DFigure * figure;
    figureCount = _object->getFigureCount();
    for (k = 0; k < figureCount; k++)
    {
        figure = _object->getFigurePtr(k);

        if (figure == NULL)
            continue;

        using namespace ThallCode;

        Xferlist * xferList = convertM3DtoXfer(figure);
        Pointlist_server2 pList;

		// It might be just as well to always use a constant instead of granularity below.
        pList.init(xferList);
        pList.ComputeSubdivPointCloud(granularity);
        pList.subdivboundaryinfo(&numBPoints, &bpoints);

        for (i = 0; i < numBPoints; i++)
        {
            x = bpoints[i].pnt[0] + bpoints[i].norm[0] * bpoints[i].rad * cutoff;
            y = bpoints[i].pnt[1] + bpoints[i].norm[1] * bpoints[i].rad * cutoff;
            z = bpoints[i].pnt[2] + bpoints[i].norm[2] * bpoints[i].rad * cutoff;

            if (x < minX) minX = x;
            else if (x > maxX) maxX = x;

            if (y < minY) minY = y;
            else if (y > maxY) maxY = y;

            if (z < minZ) minZ = z;
            else if (z > maxZ) maxZ = z;
        }

        // Now determine minimum radius of primitives in the object
        int numAtoms = figure->getPrimitiveCount();
        for (i = 0; i < numAtoms; i++)
        {
            M3DPrimitive * prim = figure->getPrimitivePtr(i);
            double rad = prim->getR();
            if (rad < minR)
				minR = rad;
        }

        delete [] (xferList->atomlist);
        delete xferList;
    }
    _bbox.size.setX(maxX - minX);
    _bbox.size.setY(maxY - minY);
    _bbox.size.setZ(maxZ - minZ);
    _bbox.origin.set(minX, minY, minZ);

    // Compute the dimensions of the distance map, _mapSize
    double max_length = _bbox.size.getX();
    if (_bbox.size.getY() > max_length)
        max_length = _bbox.size.getY();
    if (_bbox.size.getZ() > max_length)
        max_length = _bbox.size.getZ();

    // To compute mapSize, first compute a tentative map size as if the
    // smallest radius in the object were oriented along the largest
    // axis of the volume.  Try at least 3 voxels for minimal radius.
    // A 1 voxel thickness is the real limit - exceeding it will cause
    // problems when the map is used.
    double max_size = 3 * max_length / minR;
    if (max_size < min_mapsize)
		max_size = min_mapsize;
    if (max_size > max_mapsize) {
        max_size = max_mapsize;
        if (max_size < max_length / minR)
            cerr << "Warning: distance map is too small for this model." << endl;
    }
    _bbox.voxelsize = max_length/max_size; // should be minR/3 if everything is ok

    _mapSize[0] = (int) (_bbox.size.getX()/_bbox.voxelsize);
    _mapSize[1] = (int) (_bbox.size.getY()/_bbox.voxelsize);
    _mapSize[2] = (int) (_bbox.size.getZ()/_bbox.voxelsize);

    return true;
}


// This initialize allows user to set the mapsize and the bounding box size.
// This returns -1 upon failure, 1 if the map size had to be adjusted, and
// 0 for normal success.
int DistanceMap3D::initializeWithBBox(M3DObject * _object, BoundingBox * _bbox,
                            const int _mapSize[3], int granularity, double cutoff,
							bool ignore)
{
    int figureCount, i, k;
    double minR;

    if (_object == NULL) {
        cout << "DistanceMap3D object not initialized" << endl;
        return -1;
    }

    cleanup();

    object = new M3DObject(*_object);
    numMedialSubdivisions = granularity;
    factor = cutoff;

    bboxOrigin = _bbox->origin;
    bboxSize[0] = _bbox->size.getX();
    bboxSize[1] = _bbox->size.getY();
    bboxSize[2] = _bbox->size.getZ();

    mapSize[0] = _mapSize[0];
    mapSize[1] = _mapSize[1];
    mapSize[2] = _mapSize[2];

    // Compute minR for computing the distance map's voxel size
    minR = 1.0;
    M3DFigure * figure;
    figureCount = object->getFigureCount();
    for (k = 0; k < figureCount; k++)
    {
        figure = object->getFigurePtr(k);

        if (figure == NULL)
            continue;

        using namespace ThallCode;

        Xferlist * xferList = convertM3DtoXfer(figure);

        // Now determine minimum radius of primitives in the object
        int numAtoms = figure->getPrimitiveCount();
        for (i = 0; i < numAtoms; i++)
        {
            M3DPrimitive * prim = figure->getPrimitivePtr(i);
            double rad = prim->getR();
            if (rad < minR)
				minR = rad;
        }

        delete [] (xferList->atomlist);
        delete xferList;
    }

    double max_length = bboxSize[0];
    if (bboxSize[1] > max_length)
        max_length = bboxSize[1];
    if (bboxSize[2] > max_length)
        max_length = bboxSize[2];

    // Try at least 3 voxels for minimal radius (1 is HARD LIMIT-> crash otherwise)
    double max_size = 3 * max_length / minR;
    if (max_size < min_mapsize) max_size = min_mapsize;
    if (max_size > max_mapsize) {
        max_size = max_mapsize;
        if (max_size < max_length / minR)
            cerr << "Warning: distance map is too small for this model." << endl;
    }
    bboxVoxelsize = max_length/max_size; // should be minR/3 if everything is ok
#ifdef DEBUG
    cout << "voxel size = " << bboxVoxelsize << endl;
#endif

    i = (int) (bboxSize[0]/bboxVoxelsize);
    bool adjust = false;
    if (mapSize[0] != i) {
		if (! ignore)
			mapSize[0] = i;
        adjust = true;
    }
    i = (int) (bboxSize[1]/bboxVoxelsize);
    if (mapSize[1] != i) {
		if (! ignore)
			mapSize[1] = i;
        adjust = true;
    }
    i = (int) (bboxSize[2]/bboxVoxelsize);
    if (mapSize[2] != i) {
		if (! ignore)
			mapSize[2] = i;
        adjust = true;
    }
    if (adjust) {
		if (ignore) {
			bboxVoxelsize = bboxSize[0]/mapSize[0];
            cout << "Warning: Distance map size not adequate for desired resolution" << endl;
		}
		else
            cout << "Had to adjust distance map size to " << mapSize[0]
                << " x " << mapSize[1] << " x " << mapSize[2] << endl;
	}
	return (adjust ? 1 : 0);
}


// qiong: This function can be kept intact, updates should be added into 'createTreeNodeMap' instead.
//		Compute the distance map for the entire object complex. If there're multiple objects, 
//		the maps for all objects will be merged into one overall map
bool DistanceMap3D::createMap(bool r_normalize)
{
    int i;
    int treeCount;
    M3DFigureTreeNode * treeNode;
    FiguralCoordinate * tempCoordMap;


    size = mapSize[0] * mapSize[1] * mapSize[2];

    if(object == NULL || size == 0 || figuralCoordinateMap != NULL)
    {
        cout << "DistanceMap3D object not properly initialized" << endl;
        return false;
    }

    figuralCoordinateMap = new FiguralCoordinate[size];
	if (figuralCoordinateMap == NULL) {
		size = 0;
        cout << "Insufficient memory for figural coordinates map" << endl;
        return false;
    }

    tempCoordMap = NULL;
    treeCount = object->getFigureTreeCount();

    treeNode = object->getFigureTreeRoot(0);
    if (! createTreeNodeMap(treeNode, figuralCoordinateMap, r_normalize)) {
		delete figuralCoordinateMap;
		figuralCoordinateMap = NULL;
		r_normalized = false;
		return false;
	}

	bool ret = true;
    if(treeCount > 1)
    {
        tempCoordMap = new FiguralCoordinate[size];
		if (tempCoordMap == NULL)
			cout << "Insufficient memory for figural coordinates map" << endl;
		else {
			for(i = 1; i < treeCount; i++)
			{
				treeNode = object->getFigureTreeRoot(i);
				if(treeNode == NULL)
					continue;

				ret = createTreeNodeMap(treeNode, tempCoordMap, r_normalize);
				if (ret == false)
					break;

				distanceMapUnion(figuralCoordinateMap, tempCoordMap);
			}
        }
    }

    if (tempCoordMap != NULL)
        delete [] tempCoordMap;
	r_normalized = r_normalize;
	if (ret == false) {
		delete figuralCoordinateMap;
		figuralCoordinateMap = NULL;
		r_normalized = false;
	}
	return ret;
}

// qiong: function updated.
//		for single-figure object, createFigureMaps is still used
//		otherwise, the function should create the map for the blended multi-figure object
bool DistanceMap3D::createTreeNodeMap(M3DFigureTreeNode * treeNode,
                    FiguralCoordinate * map, bool r_normalize)
{
    if(treeNode == NULL)
        return false;

	if(treeNode->getChildCount()==0)
	{
		// for single-figure object, sam function createFigureMaps is called
		if (! createFigureMaps(treeNode->getFigureId(), map, r_normalize))
			return false;
	}
	else
	{
		// I broke this: FiguralCoordinate included blend features, even
		// if no sub-figure exists. Since there's one of these for *each*
		// voxel in the distance map, RAM was exceeded (one run required 800MB
		// RAM just for the FiguralCoordinate array). So I moved all those blend
		// features to BlendedFiguralCoordinate, striping them out of FiguralCoordinate
		// (reducing it's size from 48 bytes -> 24 bytes) and now we can calc
		// distance maps with a smaller bboxVoxelSize. Someone else (qiong) should
		// polymorphize the code to make a BlendFiguralCoordiate only when
		// needed, or a FiguralCoordinate otherwise.  -GST 20070619
#ifdef _BROKEN
		// otherwise, a new function createBlendedObjectMaps is created
		if(! createBlendedObjectMaps(treeNode->getFigureId(), treeNode, map, r_normalize))
#endif // _BROKEN
			cout << "DistanceMap3D::createTreeNodeMap is BROKEN for subfigures (qiong should fix)"
			<< endl;
			return false;

	}
	return true;
}


// Return distance for a point given in coordinates ranging [0.0, 1.0].  If
// the point is outside the bounding box of the figure, a "city-block" approximation
// of the distance is returned and inBounds will be set to false.

float DistanceMap3D::getDistance(double x, double y, double z, bool & inBounds) {
    float xIndex, yIndex, zIndex;
    float xIndexSave, yIndexSave, zIndexSave;
	bool xover, yover, zover;
	bool xunder, yunder, zunder;
	float d;

	// Convert coordinates to voxels
    xIndex = (float) ((x - bboxOrigin.getX())/bboxVoxelsize);
    yIndex = (float) ((y - bboxOrigin.getY())/bboxVoxelsize);
    zIndex = (float) ((z - bboxOrigin.getZ())/bboxVoxelsize);

	// The trilerp function always uses the cell obtained by
	// truncation of the arguments and the next highest one
	// in all three directions.
	if (xIndex < 0.0f) {
		xunder = true;
		xIndexSave = xIndex;
		xIndex = 0.0f;
	}
	else {
		xunder = false;

		if (((int) xIndex) >= mapSize[0] - 1) {
			xover = true;
			xIndexSave = xIndex;
			xIndex = mapSize[0] - 2;	// See note above
		}
		else
			xover = false;
	}

	if (yIndex < 0.0f) {
		yunder = true;
		yIndexSave = yIndex;
		yIndex = 0.0f;
	}
	else {
		yunder = false;

		if (((int) yIndex) >= mapSize[1] - 1) {
			yover = true;
			yIndexSave = yIndex;
			yIndex = mapSize[1] - 2;	// See note above
		}
		else
			yover = false;
	}

	if (zIndex < 0.0f) {
		zunder = true;
		zIndexSave = zIndex;
		zIndex = 0.0f;
	}
	else {
		zunder = false;

		if (((int) zIndex) >= mapSize[2] - 1) {
			zover = true;
			zIndexSave = zIndex;
			zIndex = mapSize[2] - 2;	// See note above
		}
		else
			zover = false;
	}

    d = trilerp(xIndex, yIndex, zIndex);
    inBounds = !(xunder || xover || yunder || yover || zunder || zover);

	// These increments are, of course, not correct, but since they
	// only occur at the border, they should suffice.  This is also
	// far more accurate than just returning INFINITE_DISTANCE.
	if (! inBounds) {
		if (mapIsNormalized()) {
			int x, y, z;
			double rad;
			FiguralCoordinate * fc;

			x = (int) (0.5 + xIndex);	// Ideally, rad should be trilerped
			y = (int) (0.5 + yIndex);
			z = (int) (0.5 + zIndex);
			fc = &figuralCoordinateMap[x + mapSize[0]*(y + mapSize[1]*z)];
			rad = fc->r;
#ifdef _BLENDED
			if (fc->inBlendRegion) {
				rad += fc->r2;		// Simplest to use the average
				rad *= 0.5;
			}
#endif // _BLENDED
			d *= rad;
			if (xunder)
				d += fabs(xIndexSave)*bboxVoxelsize + bboxOrigin.getX();
			else if (xover)
				d += (xIndexSave - mapSize[0])*bboxVoxelsize + bboxOrigin.getX();
			if (yunder)
				d += fabs(yIndexSave)*bboxVoxelsize + bboxOrigin.getY();
			else if (yover)
				d += (yIndexSave - mapSize[1])*bboxVoxelsize + bboxOrigin.getY();
			if (zunder)
				d += fabs(zIndexSave)*bboxVoxelsize + bboxOrigin.getZ();
			else if (zover)
				d += (zIndexSave - mapSize[2])*bboxVoxelsize + bboxOrigin.getZ();
			d /= rad;
		}
		else {
			if (xunder)
				d += fabs(xIndexSave)*bboxVoxelsize + bboxOrigin.getX();
			else if (xover)
				d += (xIndexSave - mapSize[0])*bboxVoxelsize + bboxOrigin.getX();
			if (yunder)
				d += fabs(yIndexSave)*bboxVoxelsize + bboxOrigin.getY();
			else if (yover)
				d += (yIndexSave - mapSize[1])*bboxVoxelsize + bboxOrigin.getY();
			if (zunder)
				d += fabs(zIndexSave)*bboxVoxelsize + bboxOrigin.getZ();
			else if (zover)
				d += (zIndexSave - mapSize[2])*bboxVoxelsize + bboxOrigin.getZ();
		}
	}

	return d;
}

float DistanceMap3D::getDistance(double x, double y, double z)
{
    bool temp;
    return getDistance(x, y, z, temp);
}

float DistanceMap3D::getApproxDistance(double x, double y, double z) {
    int i, j, k;

    i = (int) ((float) ((x - bboxOrigin.getX())/bboxVoxelsize));
    j = (int) ((float) ((y - bboxOrigin.getY())/bboxVoxelsize));
    k = (int) ((float) ((z - bboxOrigin.getZ())/bboxVoxelsize));

    if (i < 0 || j < 0 || k < 0)
       return INFINITE_DISTANCE;

	if (i >= mapSize[0])
       return INFINITE_DISTANCE;

	if (j >= mapSize[1])
       return INFINITE_DISTANCE;

	if (k >= mapSize[2])
       return INFINITE_DISTANCE;

    return getDistance(i, j, k);
}


// Takes modelCoord and converts it into figural coordinates (u,v,t) in figCoord
bool DistanceMap3D::figuralCoordinates(Vector3D & figCoord,
                                       const Vector3D & modelCoord)
{
	FiguralCoordinate * fcPtr = figuralCoordinates(modelCoord);
	if (fcPtr == NULL)
		return false;

    figCoord.setX(fcPtr->u);
    figCoord.setY(fcPtr->v);
    figCoord.setZ(fcPtr->t);
	return true;
}


bool DistanceMap3D::figuralCoordinates(FiguralCoordinate & figCoord,
                                       const Vector3D & modelCoord)
{
	FiguralCoordinate * fcPtr = figuralCoordinates(modelCoord);
	if (fcPtr == NULL)
		return false;

    figCoord = *fcPtr;
	return true;
}


FiguralCoordinate * DistanceMap3D::figuralCoordinates(const Vector3D & modelCoord)
{
    Vector3D newCoord;
    int xIndex, yIndex, zIndex;
    int index;

    newCoord = (modelCoord - bboxOrigin);
    newCoord /= bboxVoxelsize;

    xIndex = (int) (newCoord.getX());
    yIndex = (int) (newCoord.getY());
    zIndex = (int) (newCoord.getZ());

    if (xIndex < 0)
        return NULL;
	else if (xIndex >= mapSize[0])
        return NULL;
	if (yIndex < 0)
        return NULL;
	else if (yIndex >= mapSize[1])
        return NULL;
    if (zIndex < 0)
        return NULL;
	else if (zIndex >= mapSize[2])
        return NULL;

    index = xIndex + mapSize[0] * (yIndex + mapSize[1] * zIndex);

    return &figuralCoordinateMap[index];
}


FiguralCoordinate * DistanceMap3D::figuralCoordinates(const Vector3D & modelCoord, bool & outside)
{
    Vector3D newCoord;
    int xIndex, yIndex, zIndex;
    int index;

    newCoord = (modelCoord - bboxOrigin);
    newCoord /= bboxVoxelsize;

    xIndex = (int) (newCoord.getX());
    yIndex = (int) (newCoord.getY());
    zIndex = (int) (newCoord.getZ());

	outside = true;
    if (xIndex < 0) {
        xIndex = 0;
		outside = false;
	}
	else if (xIndex >= mapSize[0]) {
        xIndex = mapSize[0] - 1;
		outside = false;
	}
	if (yIndex < 0) {
        yIndex = 0;
		outside = false;
	}
	else if (yIndex >= mapSize[1]) {
        yIndex = mapSize[1] - 1;
		outside = false;
	}
    if (zIndex < 0) {
        zIndex = 0;
		outside = false;
	}
	else if (zIndex >= mapSize[2]) {
        zIndex = mapSize[2] - 1;
		outside = false;
	}

    index = xIndex + mapSize[0] * (yIndex + mapSize[1] * zIndex);

    return &figuralCoordinateMap[index];
}


bool DistanceMap3D::figuralCoordinateIndex(const Vector3D & modelCoord, Vector3D & mapCoord)
{
    Vector3D newCoord;
    int xIndex, yIndex, zIndex;
	bool outside;

    newCoord = (modelCoord - bboxOrigin);
    newCoord /= bboxVoxelsize;

    xIndex = (int) (newCoord.getX());
    yIndex = (int) (newCoord.getY());
    zIndex = (int) (newCoord.getZ());

	outside = true;
    if (xIndex < 0) {
        xIndex = 0;
		outside = false;
	}
	else if (xIndex >= mapSize[0]) {
        xIndex = mapSize[0] - 1;
		outside = false;
	}
	if (yIndex < 0) {
        yIndex = 0;
		outside = false;
	}
	else if (yIndex >= mapSize[1]) {
        yIndex = mapSize[1] - 1;
		outside = false;
	}
    if (zIndex < 0) {
        zIndex = 0;
		outside = false;
	}
	else if (zIndex >= mapSize[2]) {
        zIndex = mapSize[2] - 1;
		outside = false;
	}

	mapCoord.setX(xIndex);
	mapCoord.setY(yIndex);
	mapCoord.setZ(zIndex);

	return outside;
}


void DistanceMap3D::getDistanceMap(float ** mapPtr) {
	*mapPtr = new float[size];
	for (int i = 0; i < size; i++)
		(*mapPtr)[i] = figuralCoordinateMap[i].distance;
}


bool DistanceMap3D::createFigureMaps(unsigned short figureId,
         FiguralCoordinate * map, bool r_normalize)
{
    M3DFigure * figure;
    Xferlist * xferList;
    Bpoint * bpointList;
    int tileCount, pointCount;
    int index, i, ratio;
    int x, y, z;
    char * image;
    unsigned long * bpointIndexMap;
    ThallCode::Pointlist_server2 * p_list;
	const char * msg = "Insufficient memory to create figural map";

    if(object == NULL)
        return false;

    figure = object->getFigurePtr(figureId);
    if(figure == NULL)
        return false;

    using namespace ThallCode;

    image = (char *) calloc(size, sizeof(char));
	if (image == NULL) {
		cout << msg << endl;
		return false;
	}

    bpointIndexMap = new unsigned long[size];
	if (bpointIndexMap == NULL) {
		free(image);
		cout << msg << endl;
		return false;
	}

#ifdef DEBUG
    printf("Creating figure %d map.\n", figureId);
#endif

    xferList = convertM3DtoXfer(figure);

    p_list = new Pointlist_server2;
    p_list->init(xferList);
    p_list->ComputeSubdivBoundaryTiles(numMedialSubdivisions);
    p_list->subdivtileinfo(&tileCount, &bpointList);
    pointCount = 4*tileCount;
    ratio = 4;
#ifdef DEBUG
	cout << "createFigureMaps(): tileCount = " << tileCount << endl;
#endif

    for (i = 0; i < pointCount; i++)
    {
        x = (int) ((bpointList[i].pnt[0] - bboxOrigin.getX()) / bboxVoxelsize);
        y = (int) ((bpointList[i].pnt[1] - bboxOrigin.getY()) / bboxVoxelsize);
        z = (int) ((bpointList[i].pnt[2] - bboxOrigin.getZ()) / bboxVoxelsize);

        // Build binary image of the surface for distance map calculation
        index = x + mapSize[0]*(y + mapSize[1]*z);
        image[index] = 1;
        bpointIndexMap[index] = i;
    }
//if (figureId == 0) {
//snapshot("C:/temp/image.raw3", (unsigned char *) image, mapSize[0], mapSize[1], mapSize[2]);
//exit(0);
//}

//	cout << "Push_back p_list at 0x" << p_list << endl;
	tileLists.push_back(p_list);	// AGG, why can't this go after the function calls?
	figureIdList.push_back(figureId);

	bool ret;
    if (r_normalize)
        ret = createNormalizedSignedDistanceBpointMaps(image, bpointIndexMap, bpointList,
		    map, figureId, ratio);
    else
        ret = createRawSignedDistanceBpointMaps(image, bpointIndexMap, bpointList,
		    map, figureId, ratio);

    free(image);
    delete [] bpointIndexMap;
    delete [] (xferList->atomlist);
    delete xferList;

	return ret;
}


/*  Create the distance map for a multifigural object.
	Instead of using the merge of the individual distance maps of the main 
	figure and subfigures, the map is calculated from the blended surface.
	NOTE: at this point multi-level support is not yet in place.
*/
bool DistanceMap3D::createBlendedObjectMaps(unsigned short figureId,
					M3DFigureTreeNode * treeNode, FiguralCoordinate * map, bool r_normalize)
{
    Bpoint * bpointList;
    int tileCount, pointCount;
    int index, i, ratio;
    int x, y, z;
    char * image;
    unsigned long * bpointIndexMap;
    ThallCode::Pointlist_serverB * p_list;
	const char * msg = "Insufficient memory to create figural map";

    if(object == NULL)
        return false;

	int *unused=new int[object->getFigureCount()];

    using namespace ThallCode;

    image = (char *) calloc(size, sizeof(char));
	if (image == NULL) {
		cout << msg << endl;
		return false;
	}

    bpointIndexMap = new unsigned long[size];
	if (bpointIndexMap == NULL) {
		free(image);
		cout << msg << endl;
		return false;
	}

#ifdef DEBUG
    printf("Creating multi-figure object %d map.\n", figureId);
#endif

	// Create an object containing a copy of the figure tree
	M3DObject *subObject=new M3DObject;
	subObject->addTree(treeNode, object, unused); //subObject->addTree(root, object); changed according to changed in M3DObject 2003/07/15

    p_list = new Pointlist_serverB;
    p_list->init(subObject);

	// level 5 is way too high to use
	// it takes a lot of time to calculate however does not make much difference, so level 4 will be used instead
	// hardcoded '-1'
	// might be adapted to dynamically determine level to use by thresholding the size of the tiles
    p_list->ComputeSubdivBoundaryTiles(numMedialSubdivisions > 4 ? numMedialSubdivisions - 1 : numMedialSubdivisions); 
    p_list->subdivtileinfo(&tileCount, &bpointList);

	delete [] unused;
	delete subObject;
	subObject=NULL;

    pointCount = 4*tileCount;
    ratio = 4;

#ifdef DEBUG
	cout << "createFigureMaps(): tileCount = " << tileCount << endl;
#endif

    for (i = 0; i < pointCount; i++)
    {
        x = (int) ((bpointList[i].pnt[0] - bboxOrigin.getX()) / bboxVoxelsize);
        y = (int) ((bpointList[i].pnt[1] - bboxOrigin.getY()) / bboxVoxelsize);
        z = (int) ((bpointList[i].pnt[2] - bboxOrigin.getZ()) / bboxVoxelsize);

        // Build binary image of the surface for distance map calculation
        index = x + mapSize[0]*(y + mapSize[1]*z);
        image[index] = 1;
        bpointIndexMap[index] = i;
    }

	// qiong: NULL is pushed back instead of the p_list itself
	tileLists.push_back(NULL);	//tileLists.push_back(p_list);	// AGG, why can't this go after the function calls?
	figureIdList.push_back(figureId);

	bool ret;
    if (r_normalize)
        ret = createNormalizedSignedDistanceBpointMaps(image, bpointIndexMap, bpointList,
		    map, figureId, ratio);
    else
        ret = createRawSignedDistanceBpointMaps(image, bpointIndexMap, bpointList,
		    map, figureId, ratio);

    free(image);

    delete [] bpointIndexMap;

	// qiong: ??
	delete p_list;

	return ret;
}


void DistanceMap3D::tiles(int tileListId, int *ntiles, Bpoint **bptlistptr) {
        tileLists[tileListId]->subdivtileinfo(ntiles, bptlistptr);
};


bool DistanceMap3D::createNormalizedSignedDistanceBpointMaps(char * image,
    unsigned long * bpointIndexMap, Bpoint * bpointList, FiguralCoordinate * map,
	unsigned short figureId, int ratio)
{
    short * odx;
    short * ody;
    short * odz;
    int index, i;
    float innerProd;
    Bpoint * bpointPtr;

    float dist;


#ifdef BINARY
	if ((int) tuningWt(BpDebug1) == 1) {
		cout << "Snapshot bdryIm2.raw3" << endl;
		snapshot("bdryIm2.raw3", (unsigned char *) image, mapSize[0], mapSize[1], mapSize[2],
			0, 255, 1.0f, 1.0f, 1.0f);
	}

    if (edt3ddan(image, mapSize[0], mapSize[1], mapSize[2], 1., 1., 1., 0, &odx, &ody, &odz) != 0) {
		cout << "Insufficient memory for distance map" << endl;
		return false;
	}
#else
    if (edt3ddan(image, mapSize[0], mapSize[1], mapSize[2], 0, &odx, &ody, &odz) != 0) {
		cout << "Computation of distance map failed" << endl;
		return false;
	}
#endif
//snapshot("raw_dmap.raw3", image, mapSize[0], mapSize[1], mapSize[2]);

#ifdef DEBUG
    cout << "createNormalizedSignedDistanceBpointMaps(): size = " << size << endl;
#endif
	for (i = 0; i < size; i++)
    {
		// Make every point in the volume point to the surface
        index = i + odx[i] + mapSize[0]*(ody[i] + mapSize[1]*odz[i]);	// Find point on surface
        bpointIndexMap[i] = bpointIndexMap[index];

		// Note: storing dist*dist to avoid this sqrt gives only a relatively small (3%)
		// speedup for large maps.
        dist = (float) (sqrt((double) (odx[i]*odx[i] + ody[i]*ody[i] + odz[i]*odz[i])) *
                        bboxVoxelsize);
		bpointPtr = bpointList + bpointIndexMap[i];

        // Normalize distance map by the local radius
        map[i].r = (float) bpointPtr->rad;
        dist /= (float) bpointPtr->rad;

        innerProd = odx[i]*bpointPtr->norm[0] + ody[i]*bpointPtr->norm[1]
            + odz[i]*bpointPtr->norm[2];

        if(innerProd > 0.0)
            dist = -dist;

        map[i].distance = dist;		// Voxels/R
        map[i].u = bpointPtr->u;
        map[i].v = bpointPtr->v;
        map[i].t = bpointPtr->t;
	    map[i].tileId = bpointIndexMap[i]/ratio;
#ifdef _BLENDED
        map[i].inBlendRegion = false;
#endif // _BLENDED
    }

    free(odx);
    free(ody);
    free(odz);

	return true;
}


bool DistanceMap3D::createRawSignedDistanceBpointMaps(char * image,
    unsigned long * bpointIndexMap, Bpoint * bpointList, FiguralCoordinate * map,
	unsigned short figureId, int ratio)
{
    short * odx;
    short * ody;
    short * odz;
    int index, i;
    float innerProd;
    Bpoint * bpointPtr;

    float dist;


#ifdef BINARY
	// mapSpacing might not be right here... -gst 20050615
    if (edt3ddan(image, mapSize[0], mapSize[1], mapSize[2], 1., 1., 1.,
		0, &odx, &ody, &odz) != 0)
	{
		cout << "Insufficient memory for distance map" << endl;
		return false;
	}
#else
    if (edt3ddan(image, mapSize[0], mapSize[1], mapSize[2], 0, &odx, &ody, &odz) != 0) {
		cout << "Computation of distance map failed" << endl;
		return false;
	}
#endif

#ifdef DEBUG
    cout << "createRawSignedDistanceBpointMaps(): size = " << size << endl;
#endif
    for(i = 0; i < size; i++)
    {
		// Make every point in the volume point to the surface
		short x = odx[i];
		short y = ody[i];
		short z = odz[i];

        index = i + x + mapSize[0]*(y + mapSize[1]*z);	// Find point on surface
		if (index < i) {
			cout << "INDEX < I " << index << " " << i << endl;
		}
        bpointIndexMap[i] = bpointIndexMap[index];

        dist = (float) (sqrt((double) (x*x + y*y + z*z)) * bboxVoxelsize);
		bpointPtr = bpointList + bpointIndexMap[i];

		// Do not normalize
        map[i].r = (float) bpointPtr->rad;

        innerProd =
			  x*bpointPtr->norm[0]
			+ y*bpointPtr->norm[1]
            + z*bpointPtr->norm[2];

        if(innerProd > 0.0)
            dist = -dist;

        map[i].distance = dist;		// Voxels
        map[i].u = bpointPtr->u;
        map[i].v = bpointPtr->v;
        map[i].t = bpointPtr->t;
	    map[i].tileId = bpointIndexMap[i]/ratio;
#ifdef _BLENDED
        map[i].inBlendRegion = false;
#endif // _BLENDED
    }

    free(odx);
    free(ody);
    free(odz);

	return true;
}


void DistanceMap3D::range(float & min, float & max) {
    register int i;
	register float val0, val1;

    if (status() == false) return;

	min = figuralCoordinateMap[0].distance;
	max = figuralCoordinateMap[0].distance;

	int start = size & 0x01;

	// Compute the intensity range of the distance map
    for(i = start; i < size; i += 2)
    {
        val0 = figuralCoordinateMap[i].distance;
        val1 = figuralCoordinateMap[i + 1].distance;

        if(val0 < val1) {
			if(val0 < min)
				min = val0;
			if(val1 > max)
				max = val1;
		}
        else {
			if(val1 < min)
				min = val1;
			if(val0 > max)
				max = val0;
		}
    }
}


void DistanceMap3D::distanceMapUnion(FiguralCoordinate * map1,
                                     FiguralCoordinate * map2)
{
#ifdef DEBUG
    printf("Taking map union.\n");
#endif

    for(int i = 0; i < size; i++)
    {
        if(map1[i].distance > map2[i].distance)
            map1[i] = map2[i];
    }
}


void DistanceMap3D::distanceMapBlendRaw(BlendedFiguralCoordinate * map1,
                                     BlendedFiguralCoordinate * map2,
                                     float blendExtent, float blendAmount)
{
    float minVal,
          correction,
          u, u2, u3;

    register float d1r, d2r;
    register float d1, d2;

#ifdef DEBUG
    printf("Blending maps.\n");
#endif

    for(int i = 0; i < size; i++)
    {
        d1 = map1[i].distance;
        d2 = map2[i].distance;
        d1r = d1/map1[i].r;
        d2r = d2/map2[i].r;

        minVal = (d1 < d2) ? d1 : d2;

        u = fabs(d1r - d2r);

        if (u <= blendExtent)
        {
            u = u / blendExtent;
            u2 = u * u;
            u3 = u2 * u;

            correction = blendExtent *
                (blendAmount * (0.5 - 1.5*u2 + u3) - 0.5 * (u - 2.0*u2 + u3));

			// Since the blending equation requires r-proportional distances, an
			// adjustment must be made.  For now mutliply by average r, which gives
			// only an approximate result.
            map1[i].distance = minVal - correction*(map1[i].r + map2[i].r)/2.0;
            map1[i].inBlendRegion = true;
			if (d1r > d2r) {
				map1[i].u = map2[i].u;
				map1[i].v = map2[i].v;
				map1[i].t = map2[i].t;
			}
        }
		else {
			if (d1r > d2r) {
				map1[i].u = map2[i].u;
				map1[i].v = map2[i].v;
				map1[i].t = map2[i].t;
				map1[i].distance = map2[i].distance;
			    map1[i].inBlendRegion = false;
			}
		}

		// Tileset 2 is always the subfigure
		map1[i].r2 = map2[i].r;
		map1[i].tileId2 = map2[i].tileId;
    }
}


void DistanceMap3D::distanceMapBlendNormalized(BlendedFiguralCoordinate * map1,
                                     BlendedFiguralCoordinate* map2,
                                     float blendExtent, float blendAmount)
{
    float minVal,
          correction,
          u, u2, u3;

    register float d1, d2;

#ifdef DEBUG
    printf("Blending maps.\n");
#endif

    for(int i = 0; i < size; i++)
    {
        d1 = map1[i].distance;
        d2 = map2[i].distance;

        minVal = (d1 < d2) ? d1 : d2;

        u = fabs(d1 - d2);

        if (u <= blendExtent)
        {
            u = u / blendExtent;
            u2 = u * u;
            u3 = u2 * u;

            correction = blendExtent *
                (blendAmount * (0.5 - 1.5*u2 + u3) - 0.5 * (u - 2.0*u2 + u3));

            map1[i].distance = minVal - correction;
            map1[i].inBlendRegion = true;
			if (d1 > d2) {
				map1[i].u = map2[i].u;
				map1[i].v = map2[i].v;
				map1[i].t = map2[i].t;
			}
        }
		else {
			if (d1 > d2) {
				map1[i].u = map2[i].u;
				map1[i].v = map2[i].v;
				map1[i].t = map2[i].t;
				map1[i].distance = map2[i].distance;
			    map1[i].inBlendRegion = false;
			}
		}

		// Tileset 2 is always the subfigure
		map1[i].r2 = map2[i].r;
		map1[i].tileId2 = map2[i].tileId;
    }
}

