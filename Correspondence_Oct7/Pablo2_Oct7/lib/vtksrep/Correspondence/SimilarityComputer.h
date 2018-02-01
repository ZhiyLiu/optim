/* The purpose of this class is to compute similarity between two objects.
 *
 * Zhiyuan Liu
 * 2017.12
 */
#ifndef _SIMILARITYCOMPUTER_H_
#define _SIMILARITYCOMPUTER_H_

#include <vector>
#include <string>
#include "vtksrep.h"
#include <vtkSmartPointer.h>

class M3DObject;
class ImageDistanceMap;
class M3DSpoke;
class M3DQuadFigure;
class M3DQuadPrimitive;
class M3DFigure;

class SimilarityComputer {
public:
    SimilarityComputer();
    ~SimilarityComputer();

    // Description: Target image want to compare
    // Input: Signed distance image
    void setTargetImage(ImageDistanceMap* distanceImage);

    // Description: srep model want to compare to target image above
    // Input: all spokes after interpolation
    void setSrepModel(M3DObject* spokes);

    void setTargetFigureIndex(int figureId);

    void setFigure(M3DQuadFigure* figure);

    // Description: start to compute similarity measure
    // Output: similarity measure
    // Returns: true if succeeded, otherwise false
    bool compute(double *similarityMeasure);

private:
    // Description: interpolate around a spoke
    // Input: figure, quadAtom, atomId, spokeId
    // Output: neighbor spokes
    void getRelevantSpokes(M3DFigure* figure, M3DQuadPrimitive* quadAtom, int atomId, int spokeId, std::vector<M3DSpoke*>* outSpokes);

    // Description: interpolate crest spokes
    // Input: interpolation level
    // Output: crestSpokes
    void GetCrestSpokes(vtkSmartPointer<vtkSRep> srepfig, int level, vector<M3DSpoke*>& spokes, bool istube = false, vtkIdType atomId = -1, int spokeId = -1);

    // Description: generate vtk srep data structure
    // Input: figure which needs interpolation
    // Input: u/v 
    // Output: new spoke in the middle of the quad
    M3DSpoke* interpolateSpoke(M3DFigure *figure, double u, double v, int side);
    
    // Description: generate vtk srep data structure
    // Input: quad figure
    // Output: vtk srep
    // Return: vtk id
    vtkIdType GetVtkSrepFig(M3DFigure* figure, vtkSmartPointer<vtkSRep>& srepfig, int atomId=-1);

    // Description: get sum of squared distance between spokes and image
    double getSSD(std::vector<M3DSpoke*>& spokes);
private:
    ImageDistanceMap* mDistanceImage;
    M3DObject*        mSreps; // spokes after interpolation
    int               mFigureIndex;
    M3DQuadFigure*    mFigure;
};



#endif

