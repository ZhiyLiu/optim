#ifndef DISTANCE_MAP_3D_H
#define DISTANCE_MAP_3D_H

#include "M3DObject.h"
#include "renderDefinitions.h"

#define D_XFERLIST
#define D_POINTLIST_SERVER2
#include "Shapedepend.h"


/*    Class DistanceMap3D.

    This class produces 3D distance maps, according to various user-supplied
    parameters.  To use this class, one of the initialize functions must be
    called before createMap() is called.  Thereafter, the various access
    functions can be used to get values from the map.

    To limit the size of the computed map, a bounding box is used.  The
    default size is slightly larger than the object, with medial radius (r)
    twice the value of the cutoff argument.  A separate defineBBox() function
    is provided for determining a bounding box on an object that is a superset
    of the figures used for computing the map.

    The object can contain multiple figures.  If two figures are connected,
    the distance map will be based on a blending of the figures according to
    the blending parameters stored in each hinge figure.

    If createMap() is called with no argument, each voxel of the distance map
    will be normalized by the value of r of the nearest medial atom.  Setting
    r_normalize to false, results in a map in Euclidean space.  The second
    argument is used to indicate whether or not subdivision tiling is to be
    used for boundary point generation.

    Once the map is created, distances can be obtained by the 4 getDistance()
    members.  The two taking integer arguments treat the distance map like a
    simple 3D array, while the other two produce interpolated results.

    Functions are also provided to convert model coordinates (the range [0, 1])
    into both figural coordinates and map coordinates.  Function figureID()
    provides a way of determining which figures a point is near.
*/




// A distance map is an array of FiguralCoordinate objects.  For a figural tree,
// it is produced by blending maps of a the main figure and its subfigure.  This
// class presently only supports a single subfigure.  Once the maps are produced
// for figural trees, they are combined without blending.

// BlendedFiguralCoordinate is an attempt to save RAM when no sub-figures exist

class FiguralCoordinate
{
public:

    float u;	// Figural coordinates of the map point wrt. the nearest figure
    float v;
    float t;

    float distance;	// Distance of a map point to the nearest of the two figures
	float r;		// Medial radius of this point wrt. the main figure

    unsigned short tileId;		// The closest tile of the main figure

	friend std::ostream & operator<< (std::ostream & ostr, FiguralCoordinate & fc_ref);

};

class BlendedFiguralCoordinate
{
public:

    double u;	// Figural coordinates of the map point wrt. the nearest figure
    double v;
    double t;

    float distance;	// Distance of a map point to the nearest of the two figures
	float r;		// Medial radius of this point wrt. the main figure
	float r2;		// Medial radius of this point wrt. the subfigure

    unsigned short tileId;		// The closest tile of the main figure
    unsigned short tileId2;		// The closest tile of the subfigure

    bool inBlendRegion;		// True if inside the figural blending region

	friend std::ostream & operator<< (std::ostream & ostr, BlendedFiguralCoordinate & fc_ref);

};



struct BoundingBox
{
    Vector3D origin;	// Origin of the bounding box in model coordinates
	Vector3D size;		// Lengths of axes of bounding box in model coordinates
	double voxelsize;	// Model coordinate size of a distance map voxel
};



extern const Vector3D DEFAULT_BBOX_ORIGIN;
extern const Vector3D DEFAULT_BBOX_SIZE;
extern const int DEFAULT_MAP_SIZE[3];

const int MAX_MAP_SIZE = 512;
const int MIN_MAP_SIZE = 128;

const double INFINITE_DISTANCE = 1e8;


class DistanceMap3D
{
public:

    DistanceMap3D();
    DistanceMap3D(int mapSizeMinimum, int mapSizeMaximum);
	// Granularity is the level of subdivision of the medial mesh used for
	// computing tiles at the boundary, when generating the image on which
	// the distance map will be computed.
    DistanceMap3D(M3DObject * _object, int granularity = MEDIAL_SUBDIVISIONS,
        double cutoff = 0.6, int mapSizeMinimum = MIN_MAP_SIZE,
        int mapSizeMaximum = MAX_MAP_SIZE);
    virtual ~DistanceMap3D();

    void initialize();
    void initialize(M3DObject * _object, int granularity = MEDIAL_SUBDIVISIONS,
        double cutoff = 0.6, int mapSizeMinimum = MIN_MAP_SIZE,
        int mapSizeMaximum = MAX_MAP_SIZE);

    // This initialize allows user to set the mapsize and the bounding box size.
    // This returns -1 upon failure, 1 if the map size had to be adjusted, and
    // 0 for normal success.  Set ignore to true to avoid the map size adjustment.
    int initializeWithBBox(M3DObject * _object, BoundingBox * _bbox = NULL,
        const int _mapSize[3] = DEFAULT_MAP_SIZE,
		int granularity = MEDIAL_SUBDIVISIONS, double cutoff = 0.6,
		bool ignore = false);

    // Determine a bounding box for an object within the limits of the maximum and
    // minimum bounding box sizes specified to the constructor.  The results of this
    // function are intended to be passed to a later call of initializeWithBBox().
    bool defineBBox(M3DObject * _object, BoundingBox & _bbox, int _mapSize[3],
		int granularity = MEDIAL_SUBDIVISIONS, double cutoff = 0.6);

	// Functions to provide access to the tiles, the verticies of which are
	// used to generate the distance map.
    void tiles(int tileListId, int *ntiles, Bpoint **bptlistptr);
	int tileListCount() { return tileLists.size(); }
	ThallCode::Pointlist_server2 * tileList(int tileListId) {
		return tileLists[tileListId];
	}

	// Distance map creation functions
    bool createMap(bool r_normalize = true);
    bool status() {    	// Returns true if a distance map has been computed
        return figuralCoordinateMap == NULL ? false : true;
    }

    int length() { return size; }	// Returns number of voxels in the distance map
	BoundingBox getBoundingBox();
	bool mapIsNormalized() { return r_normalized; }

    // Access via distance map coordinates
    float getDistance(int index) {             // Uninterpolated
        return figuralCoordinateMap[index].distance;
	}
    float getDistance(int x, int y, int z) {             // Uninterpolated
        return figuralCoordinateMap[x + mapSize[0]*(y + mapSize[1]*z)].distance;
	}

    // Access via model (0.0 - 1.0) coordinates
	// Tri-linearly interpolated inside the bounding box and extrapolated outside
    // inBounds is true if the coordinates are inside the bounding box
    float getDistance(double x, double y, double z, bool & inBounds);
    // Same as above, without the inBounds parameter
    float getDistance(double x, double y, double z);
	// Uninterpolated - can return INFINITE_DISTANCE
    float getApproxDistance(double x, double y, double z);

    const int * getMapSize() { return mapSize; }	// Returns array of 3 axial voxel lengths
    double getMapVoxelsize() { return bboxVoxelsize; }
    void range(float & min, float & max);   // Returns the range of map values

    // Converts model coordinates into figural coordinates (u, v, t) in figCoord
    bool figuralCoordinates(Vector3D & figCoord, const Vector3D & modelCoord);

    // Functions to get the FiguralCoordinate object corresponding to model coordinates
	// The next function returns false if the model coordinates are beyond the bounding box
    bool figuralCoordinates(FiguralCoordinate & figCoord,
        const Vector3D & modelCoord);
	// Model coordinates beyond the bounding box return NULL 
    FiguralCoordinate * figuralCoordinates(const Vector3D & modelCoord);
	// Model coordinates beyond the bounding box are clipped to the bounding box 
    FiguralCoordinate * figuralCoordinates(const Vector3D & modelCoord, bool & outside);

	// Returns a pointer to the figural coordinate map
    FiguralCoordinate * getFiguralCoordinateMap() { return figuralCoordinateMap; }

	// Convert real coordinates to map coordinates
	// The return is false if modelCoord was outside the bounding box
	bool figuralCoordinateIndex(const Vector3D & modelCoord, Vector3D & mapCoord);

	// Given a tile set number from the figural coordinate class, get the figure number
	int figureID(int tileSetID) { return figureIdList[tileSetID]; }

	// Creates array mapPtr using new[] and fills it with the distances produced
	// by calling createMap().
	void getDistanceMap(float ** mapPtr);

private:

	// Create the distance map for a multi-figure object.
	// NOTE: at this point only 2-level support is working.
	bool createBlendedObjectMaps(unsigned short figureId, M3DFigureTreeNode * treeNode,
                    FiguralCoordinate * map, bool r_normalize);

    bool createTreeNodeMap(M3DFigureTreeNode * treeNode,
        FiguralCoordinate * map, bool r_normalize);
    bool createFigureMaps(unsigned short figureId, FiguralCoordinate * map,
        bool r_normalize);

	bool createRawSignedDistanceBpointMaps(char * image,
        unsigned long * bpointIndexMap, Bpoint * bpointList,
        FiguralCoordinate * map, unsigned short figureId, int ratio);
	bool createNormalizedSignedDistanceBpointMaps(char * image,
        unsigned long * bpointIndexMap, Bpoint * bpointList,
        FiguralCoordinate * map, unsigned short figureId, int ratio);

    void distanceMapUnion(FiguralCoordinate * map1, FiguralCoordinate * map2);

    void distanceMapBlendRaw(BlendedFiguralCoordinate * map1, BlendedFiguralCoordinate * map2,
        float blendExtent, float blendAmount);
    void distanceMapBlendNormalized(BlendedFiguralCoordinate * map1, BlendedFiguralCoordinate * map2,
        float blendExtent, float blendAmount);

    void cleanup();

    // Master voxelized map that contains all the information we need
    FiguralCoordinate * figuralCoordinateMap;
    int mapSize[3];
	int size;		// Voxels in the distance map
	bool r_normalized;

    M3DObject * object;

    // Object bounding cube
    double bboxSize[3];
    double bboxVoxelsize;
    Vector3D bboxOrigin;
    double factor;
    int min_mapsize, max_mapsize;

    // Defines the fineness of the interpolation
    int numMedialSubdivisions;

    float bilerp_in_z(float x, float y, int z);
    float trilerp(float x, float y, float z);

	std::vector<ThallCode::Pointlist_server2 *> tileLists;
	std::vector<unsigned short> figureIdList;
};


#endif

