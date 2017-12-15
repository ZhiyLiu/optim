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

class M3DNewuoaOptimizer : public OptimizerBase
{
public:
    M3DNewuoaOptimizer();

    M3DNewuoaOptimizer(char* rootDir, const std::vector< std::string > &inputSreps, int interpolationLevel, int spokeType);
    void initialize();

    // Description: Return a cost value, with hard-coded weight currently. The cost function was iterated automaticly by the optimizer, each time it change the values of coeff, we use the new coeff
    //Parameter: coeff is a N*46 dimesion arrary. Store the 46 spoke's u or v of all sreps. We will split this array to get each 46 for one srep.
    double getCost(const double *coeff);

    // Description: Return a cost value computed according to pre-defined entropy minimization objective function
    // Parameter: @coeff: input coefficients in each iterariton
    // @w1/w2: weight for geometric entropy and regularity
    double getObjectiveFunctionValue(const double *wholeCoeff, double w1, double w2);

    // Description: Override operator (). Required by min_newuoa.
    // Parameter: @coeff: the pointer to coefficients
    double operator () (double *coeff);

    // major entry to invoke this optimizer
    // input: logFileName- the whole path of  output log file
    // input: initalOpt- whether or not this is first iteration
    // return: error code
    int perform(const std::string& srep);

private:
    double computeSRepImageMatch(double dilationFactor);
    double computeSradPenalty();
    void interpolateSRep(std::vector<M3DSpoke> *outputSpokes);
    void updateSpokes(const double *coeff); // update spokes after each optimization
private:

    char* mVarFileDir;
    int   mInterpolationLevel;
    int   mSpokeType; // 0-upSide 1-downSide 2-crest spoke
    int   mTotalDimensions;
    int   mIterationCounter;
    M3DQuadFigure*    mShiftingQuadFig;
    std::vector<double> mSubdivisions; // get subdivision ( 2D array) according to interpolation level
    std::vector<std::string> mInputSreps;
    std::vector<M3DQuadFigure *> mQuadFigList; //holding the input sreps. Reading in before iteration.
    std::vector<vtkSmartPointer<vtkSRep> > mSrepFigList; // hoding srepfig, used for crest spoke interpolate.
};

#endif /* M3DNEWUOAOPTIMIZER_H */