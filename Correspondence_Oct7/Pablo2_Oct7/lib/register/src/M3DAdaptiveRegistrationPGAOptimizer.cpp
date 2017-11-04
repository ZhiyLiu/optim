#include <fstream>
#include <stdio.h>
#include <math.h>
#include "OptimizerBase.h"
#include "M3DAdaptiveRegistrationPGAOptimizer.h"

using namespace std;

M3DAdaptiveRegistrationPGAOptimizer::M3DAdaptiveRegistrationPGAOptimizer(M3DPGAStats * pgaPtr, int order)
	: M3DRegistrationOptimizer(true)
{
    boundsFunction = new M3DAdaptiveRegistrationPGAProblemBoundsFunction;

	problem = new M3DAdaptiveRegistrationPGAProblem(pgaPtr, order);
	numParameters = ((M3DAdaptiveRegistrationPGAProblem *) problem)->parameterCount();
}


