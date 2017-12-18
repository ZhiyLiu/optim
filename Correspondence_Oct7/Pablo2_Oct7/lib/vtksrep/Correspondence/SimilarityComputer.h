/* The purpose of this class is to compute similarity between two objects.
 *
 * Zhiyuan Liu
 * 2017.12
 */
#ifndef _SIMILARITYCOMPUTER_H_
#define _SIMILARITYCOMPUTER_H_

#include <vector>

class M3DSpoke;
class ImageDistanceMap;

class SimilarityComputer {
public:
    SimilarityComputer();

    // Description: Target image want to compare
    // Input: Signed distance image
    void setTargetImage(ImageDistanceMap* distanceImage);

    // Description: srep model want to compare to target image above
    // Input: all spokes after interpolation
    void setSrepModel(std::vector<M3DSpoke>& spokes);

    // Description: start to compute similarity measure
    // Output: similarity measure
    // Returns: true if succeeded, otherwise false
    bool compute(double *similarityMeasure);

private:
    ImageDistanceMap* mDistanceImage;
    std::vector<M3DSpoke> mSpokes; // spokes after interpolation
};

#endif

