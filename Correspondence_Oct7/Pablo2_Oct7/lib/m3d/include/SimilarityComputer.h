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
    // Output: imageMatch term: image distance
    // Output: image normal match term
    // Returns: true if succeeded, otherwise false
    bool compute(double *imageDist, double *normalMatch);

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

        // Description: get sum of squared distance between spokes and image
    double getSpokeNormalMatch(std::vector<M3DSpoke*>& spokes);

	double h1(double s)
	{
		return 2*(s * s * s) - 3*(s * s) + 1;
	}

	double h1p(double s)
	{
		return 6*(s * s) - 6*(s);
	}

	double h2(double s)
	{
		return -2*(s * s * s) + 3*(s * s);
	}

	double h2p(double s)
	{
		return -6*(s * s) + 6*(s);
	}

	double h3(double s)
	{
		return (s * s * s) - 2*(s * s) + s;		
	}

	double h3p(double s)
	{
		return 3*(s * s) - 4*(s) + 1;
	}

	double h4(double s)
	{
		return (s * s * s) - (s * s);
	}

	double h4p(double s)
	{
		return 3*(s * s) - 2*(s);
	}

private:
    ImageDistanceMap* mDistanceImage;
    M3DObject*        mSreps; // spokes after interpolation
    int               mFigureIndex;
    M3DQuadFigure*    mFigure;
};



#endif

