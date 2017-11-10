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

    double getCost(const double *coeff);

    double getObjectiveFunctionValue(const double *wholeCoeff, double w1, double w2) const;

    double operator () (double *coeff);

    // major entry to invoke this optimizer
    // input: logFileName- the whole path of  output log file
    // input: initalOpt- whether or not this is first iteration
    // return: error code
    int perform(const char* logFileName, bool initialOpt);

private:
    bool isCorrectMove(double *coeff, int length, double moveDis);
private:

    char* mVarFileDir;
    int   mInterpolationLevel;
    int   mSpokeType;
    int   mTotalDimensions;
    int   mIterationCounter;
    M3DQuadFigure*    mShiftingQuadFig;
    std::vector<double> mSubdivisions; // get subdivision ( 2D array) according to interpolation level
    std::vector<std::string> mInputSreps;
    std::vector<M3DQuadFigure *> mQuadFigList; //holding the input sreps. Reading in before iteration.
    std::vector<vtkSmartPointer<vtkSRep> > mSrepFigList; // hoding srepfig, used for crest spoke interpolate.
};

#endif /* M3DNEWUOAOPTIMIZER_H */