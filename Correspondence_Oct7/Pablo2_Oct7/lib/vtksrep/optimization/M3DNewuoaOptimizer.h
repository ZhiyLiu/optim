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

class M3DNewoaOptimizer : public OptimizerBase
{
public:
    M3DNewoaOptimizer(const char* rootDir, const std::vector< std::string > &inputSreps, int interpolationLevel, int spokeType);
    void initialize();

    double getCost(const double coeff);

    double getObjectiveFunctionValue(const double wholeCoeff, double w1, double w2) const;

    double operator () (const double stepSize);

    // major entry to invoke this optimizer
    // input: logFileName- the whole path of  output log file
    // input: initalOpt- whether or not this is first iteration
    // return: error code
    int perform(const char* logFileName, bool initialOpt);

private:
    bool isCorrectMove(const double stepSize);
private:
    M3DNewoaOptimizer();
    char* mRootDir;
    int   mInterpolationLevel;
    int   mSpokeType;
    int   mTotalDimensions;
    int   mIterationCounter;
    std::vector<std::string> inputSreps;
}

#endif /* M3DNEWUOAOPTIMIZER_H */