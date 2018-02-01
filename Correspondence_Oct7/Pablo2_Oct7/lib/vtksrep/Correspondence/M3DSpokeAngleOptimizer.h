/* The purpose of this class is to optimize spoke angle by NEWUOA
 *
 * Zhiyuan Liu
 * 2018.01
 */
#ifndef _M3DSPOKEANGLEOPTIMIZER_H_
#define _M3DSPOKEANGLEOPTIMIZER_H_

#include <vector>
#include <string>
#include "OptimizerBase.h"
#include "vtksrep.h"
#include <vtkSmartPointer.h>


class M3DQuadFigure;
class ImageDistanceMap;

class M3DSpokeAngleOptimizer
{
public:
    M3DSpokeAngleOptimizer();

    M3DSpokeAngleOptimizer(M3DObject* sreps);

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

     // update a figure after each optimization 
    void updateFigure(const double *coeff, int figureId);

private:
    M3DObject*              mSreps;
    std::vector<M3DSpoke>   mSpokes;
    ImageDistanceMap*       mSignedDistanceImage;
    int                     mFigureIndex = 0; // The figure user want to optimize, currently have only one
};

#endif
