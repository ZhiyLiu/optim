#include <fstream>
#include <math.h>
#include "OptimizerBase.h"
#include "M3DRegistrationSimilarityOptimizer.h"


M3DRegistrationSimilarityOptimizer::M3DRegistrationSimilarityOptimizer()
	: M3DRegistrationOptimizer()
{
    numParameters = 7;		// similarity transform (rotation(3), scale(1), and translation(3))
    problem = new M3DRegistrationSimilarityProblem();
    boundsFunction = new M3DRegistrationSimilarityProblemBoundsFunction;
}


