/* This is an wrapper of optimizer using NEWUOA
 * The main purpose of wrapper is to compute cost according to specific objective function
 * 
 * Zhiyuan Liu
 * Nov. 2017 
*/

#ifndef M3DNEWUOAOPTIMIZER_H
#define M3DNEWUOAOPTIMIZER_H

#include <vector>
#include <string>
#include "OptimizerBase.h"
#include "vtksrep.h"
#include <vtkSmartPointer.h>


class M3DQuadFigure;
class ImageDistanceMap;

class M3DNewuoaOptimizer
{
public:
    M3DNewuoaOptimizer();

    M3DNewuoaOptimizer(M3DObject* sreps);

    // Description: Return a cost value computed according to pre-defined entropy minimization objective function
    // Parameter: @coeff: input coefficients in each iterariton
    // @w1/w2: weight for geometric entropy and regularity
    double getObjectiveFunctionValue(const double *wholeCoeff, double w1, double w2);

    // Description: Override operator (). Required by min_newuoa.
    // Parameter: @coeff: the pointer to coefficients
    double operator () (double *coeff);

    // major entry to invoke this optimizer
    // Each object might contain multiple figures
    // Each figure is an srep. But currently each object has only one srep
    // return: error code
    // output: s-rep model at optimum
    int perform(M3DObject* outputModel);

    // Set object want to optimize
    // Each object might contain multiple figures
    // Each figure is an srep. But currently each object has only one srep
    void setObject(M3DObject* sreps);

    // Description: Set signed distance image related to srep model
    // Input: signed distance image
    void setImage(ImageDistanceMap* imageDistanceMap);
private:
    double computeSradPenalty();
    
    void interpolateSRep(std::vector<M3DSpoke> *outputSpokes);
    
     // update a figure after each optimization 
    void updateFigure(const double *coeff, int figureId);
    
    // Description: interpolate crest spokes
    // Input: interpolation level
    // Output: crestSpokes
    void interpolateCrestSpokes(int interpolationLevel, std::vector<M3DSpoke>* crestSpokes);

    // Description: generate vtk srep data structure
    // Input: quad figure
    // Output: vtk srep
    void generateVtkSrep(M3DQuadFigure* quadfig, vtkSmartPointer<vtkSRep>& srepfig);
private:
    M3DObject*              mSreps;
    ImageDistanceMap*       mSignedDistanceImage;
    std::vector<M3DSpoke>   mSpokesAfterInterp; // all spokes after interpolation
    int                     mFigureIndex = 0; // The figure user want to optimize, currently have only one
};

#endif /* M3DNEWUOAOPTIMIZER_H */