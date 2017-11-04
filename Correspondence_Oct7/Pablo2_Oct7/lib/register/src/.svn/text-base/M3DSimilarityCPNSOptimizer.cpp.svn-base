#include <fstream>
#include <stdio.h>
#include <math.h>
#include "OptimizerBase.h"
#include "M3DSimilarityCPNSOptimizer.h"

using namespace std;

M3DSimilarityCPNSOptimizer::M3DSimilarityCPNSOptimizer(Match *m, M3DCPNSStats *cpns)
	: M3DMainFigureOptimizer(m)

{
	boundsFunction = new M3DSimilarityCPNSProblemBoundsFunction();
	problem = new M3DSimilarityCPNSProblem();
	static_cast<M3DSimilarityCPNSProblem *>(problem)->setCPNS(cpns);
}