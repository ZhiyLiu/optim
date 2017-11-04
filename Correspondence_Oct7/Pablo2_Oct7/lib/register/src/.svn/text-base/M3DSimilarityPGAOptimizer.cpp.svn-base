#include <fstream>
#include <stdio.h>
#include <math.h>
#include "OptimizerBase.h"
#include "M3DSimilarityPGAOptimizer.h"
#include "M3DTubeSimilarityPGAProblem.h"


using namespace std;

M3DSimilarityPGAOptimizer::M3DSimilarityPGAOptimizer(Match * m, M3DPGAStats * pga)
	: M3DMainFigureOptimizer(m)
{
//    boundsFunction = new M3DSimilarityPGAProblemBoundsFunction;
//    problem = new M3DSimilarityPGAProblem;
	problem = new M3DTubeSimilarityPGAProblem;
	boundsFunction = new M3DTubeSimilarityPGAProblemBoundsFunction(static_cast<M3DTubeSimilarityPGAProblem*>(problem));
	static_cast<M3DSimilarityPGAProblem *>(problem)->setPGA(pga);
	// problem->parameterCount() cannot be called until after
	// problem->initialize() is called.
}

