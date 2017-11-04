#ifndef OPTIMIZATIONUSINGNEWUOA_H
#define OPTIMIZATIONUSINGNEWUOA_H

#include "objectivefunction.h"




class optimizationusingnewuoa
{
public:
    optimizationusingnewuoa();
    optimizationusingnewuoa(const char* rootDir, vector< std::string > inputSreps, int interpolationLevel, int regEntType);

    double cost(double *coeff);

    // a member function used to call cost_function(this) and pass value to min_newuoa.
    int callNEWUOA(const char * logFileName, bool initialOpt);


private:

    int totalVars; //variables number on all of the input sreps.

    const char* rootDir;

    objectivefunction objFunc;

    int srepNum;
    int varsNum;

};

#endif // OPTIMIZATIONUSINGNEWUOA_H



