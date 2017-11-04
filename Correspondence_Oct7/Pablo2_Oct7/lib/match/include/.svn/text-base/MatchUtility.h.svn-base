#ifndef MATCH_UTILITY
#define MATCH_UTILITY


class ImageDistanceMap;

class MatchUtility
{
public:

	// BINARY PABLO
	// All this stuff gets used in BP

	// Image Match functions
	// For output ONLY - there is another distance calculation directly in Match
	// that is must faster because it doesn't have to init a new pList every time
	static double calculateAverageDistance(M3DObject * object, M3DObject *
		comparisonObject, int figureId, ImageDistanceMap * dist, int level,
		PrimNeighborhoodDefn neiType = EDGES_HAVE_2_NEIGHBORS,
		DistanceType distType = AVE_GEODESIC_DIST, int verbose = 0);

	// This is the main entry into the Method of Moments initializer routines
	static SimilarityTransform3D initializeObjectToImage(M3DObject * object,
		Image3D * image, int moments = 1|2|4,	// bitFlags: 1=translate 2=scale 4=rotate
		int verbose = 1);

	static unsigned int nPoints;
	static double max_dist;
	static double min_dist;
	static double avg_dist;
	static double avg_dist_squared;
};


#endif

