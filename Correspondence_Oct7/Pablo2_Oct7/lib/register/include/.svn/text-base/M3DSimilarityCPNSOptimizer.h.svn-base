#ifndef M3D_SIMILARITY_CPNS_OPTIMIZER_H
#define M3D_SIMILARITY_CPNS_OPTIMIZER_H

#include "M3DSimilarityCPNSProblem.h"
#include "M3DMainFigureOptimizer.h"
#include <vector>

class M3DSimilarityCPNSOptimizer : public M3DMainFigureOptimizer
{
public:

	M3DSimilarityCPNSOptimizer(Match *m, M3DCPNSStats *cpns);
	~M3DSimilarityCPNSOptimizer() { }

protected:
	const double getNumericalDerivativeEpsilon() { return 1e-7; }
};

#endif