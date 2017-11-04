#include <math.h>
#include <iostream>
#include "M3DObject.h"
#include "SimilarityTransform3D.h"
#include "matrix.h"
#include "MethodOfMoments.h"
#include "MatchUtility.h"


unsigned int MatchUtility::nPoints;
double MatchUtility::max_dist;
double MatchUtility::min_dist;
double MatchUtility::avg_dist;
double MatchUtility::avg_dist_squared;

using namespace std;


// Compute the method of moments transform to initialize the object to the
// binary image

SimilarityTransform3D MatchUtility::initializeObjectToImage(M3DObject * object,
	Image3D * image, int moments, int verbose)
{
    SimilarityTransform3D momTransform;

	if (moments) {
		Vector3D figureCenter, imageCenter, initialTranslation;
		double figureVolume, imageVolume, initialScale;
		Quat figureRotation, imageRotation, initialRotation;

		MethodOfMoments::objectMoments(object, figureVolume, figureCenter, figureRotation);
		
		MethodOfMoments::imageMoments(image, imageVolume, imageCenter, imageRotation);


		// Translation is just the difference of the two center's
		initialTranslation = imageCenter - figureCenter;

		// Note here that the SCALE is the CUBE ROOT of the RATIO of VOLUMES
		initialScale = pow(imageVolume/figureVolume, 1.0/3.0);

		// Determine the smallest angle rotation we can make to align the figure and image cov's
		MethodOfMoments::leastAngleRotation(imageRotation, figureRotation, initialRotation);

		if (verbose > -1) {
			cout << " MOM: Translation " << (moments&1 ? "(used)":"(not used)") << " : ";
			initialTranslation.print();
			cout << " MOM: Scale       " << (moments&2 ? "(used)":"(not used)") << " : ";
			cout << initialScale << endl;
			cout << " MOM: Rotation    " << (moments&4 ? "(used)":"(not used)") << " : ";
			initialRotation.print();
		}

		// Setup the similarity xform & apply it
		if (moments & 1)
			momTransform.setTranslation(initialTranslation);
		if (moments & 2)
			momTransform.setScale(initialScale);
		if (moments & 4)
			momTransform.setRotation(initialRotation);
	}
	else
		if (verbose > -1)
			cout << " MOM: user requested no moments" << endl;

	return momTransform;
}

// Rohit: This method displays the summary stats in the end.
// Well it's not really the summary but the result of the final
// call to Match::computeBinaryImageMatch(...)

double MatchUtility::calculateAverageDistance(M3DObject * object,
	M3DObject * comparisonObject, int figureId, ImageDistanceMap * dist_map,
	int level, PrimNeighborhoodDefn neiType, DistanceType distType, int verbose)
{
	if (verbose > -1)
	{
		// And the regularity calculation
		cout << " Reg: ";
		object->dist2FromAveOfNeighbors(neiType, figureId, distType, true);
		cout << " Ref: ";

		object->dist2FromObject(comparisonObject, figureId, distType, true);

		cout << " OBJECT-TO-IMAGE MATCH VALUES @ LVL " << level << "---" << endl;
		cout << "  Average Dist:     " << avg_dist << " (over " << nPoints << " points)" << endl;
		cout << "  Average Dist^2:   " << avg_dist_squared << " (match)" << endl;
		cout << "  Max Dist:         " << max_dist << endl;
		cout << "  Min Dist:         " << min_dist << endl;
		cout << endl;
	}

	return -avg_dist_squared;
}


