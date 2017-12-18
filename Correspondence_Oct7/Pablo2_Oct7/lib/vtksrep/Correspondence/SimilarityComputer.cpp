/* The purpose of this class is to compute similarity between two images
 *
 * Zhiyuan Liu
 * 2017.12
 */
#include "SimilarityComputer.h"
#include "M3DSpoke.h"
#include "ImageDistanceMap.h"

SimilarityComputer::SimilarityComputer():
    mDistanceImage(NULL){}

void SimilarityComputer::setTargetImage(ImageDistanceMap* target)
{
    mDistanceImage = target;
}

void SimilarityComputer::setSrepModel(std::vector<M3DSpoke>& spokes)
{
    mSpokes.assign(spokes);
}

bool SimilarityComputer::compute(double *similarityMeasure)
{
    // step 1: validate data
    if(mDistanceImage == NULL || mSpokes.empty())
    {
        std::cout << "[Error]SimilarityComputer requires both target image and reference model" << std::endl;
        return false;
    }
    double match = 0.0;
    double dilationFactor = 0.0;
    toolsfunc tools;

   // step 2: compute the sum of square distances as the match value
    std::vector<Vector3D> surfacePoints;
    for(int i = 0; i < mSpokes.size(); ++i)
    {
        surfacePoints.push_back(mSpokes[i].getB());
    }

    int pointsNum = surfacePoints.size();
    if(pointsNum == 0)
    {
        // ERROR
        std::cout << "[ERROR]Implied boundary after interpolation is empty in method: M3DNewuoaOptimizer::computeSRepImageMatch!" << std::endl;
        return false;
    }

    double max_dist = 0.0, min_dist = +1e20, avg_dist = 0.0;

    for(int i = 0; i < pointsNum; ++i)
    {
        double distance = fabs(mDistanceImage->getDistance(surfacePoints[i]) - dilationFactor);

        double distance = 0.0;
        if(distance > max_dist)
        {
            max_dist = distance;
        }

        if(distance < min_dist)
        {
            min_dist = distance;
        }

        match += distance * distance;
        avg_dist += distance;
    }

    
    avg_dist /= pointsNum;
    match /= pointsNum;

    MatchUtility::max_dist = max_dist;
	MatchUtility::min_dist = min_dist;
	MatchUtility::nPoints = pointsNum;
	MatchUtility::avg_dist = avg_dist;
	MatchUtility::avg_dist_squared = match;
    
    similarityMeasure = match;
    return true;
}