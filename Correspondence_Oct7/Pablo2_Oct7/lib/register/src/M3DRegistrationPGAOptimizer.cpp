#include <fstream>
#include <stdio.h>
#include <math.h>
#include "OptimizerBase.h"
#include "M3DRegistrationPGAOptimizer.h"

using namespace std;

M3DRegistrationPGAOptimizer::M3DRegistrationPGAOptimizer(M3DPGAStats * pgaPtr, int order)
	: M3DRegistrationOptimizer()
{
    boundsFunction = new M3DRegistrationPGAProblemBoundsFunction;
    problem = new M3DRegistrationPGAProblem();
    ((M3DRegistrationPGAProblem *) problem) -> setPGA(pgaPtr, order);
    numParameters = ((M3DRegistrationPGAProblem *) problem)->parameterCount();
}


