#include <fstream>
#include <math.h>
#include "OptimizerBase.h"
#include "M3DSimilarityElongationOptimizer.h"


M3DSimilarityElongationOptimizer::M3DSimilarityElongationOptimizer(Match * m,
	int count) : M3DMainFigureOptimizer(m, count)
{
    problem = new M3DSimilarityElongationProblem;
    // Default parameter count is 8; it equals 6 for binary image matching
	problem->setParameterCount(parameterCount());
    boundsFunction = new M3DSimilarityElongationProblemBoundsFunction(static_cast<M3DSimilarityElongationProblem*>(problem));
}


