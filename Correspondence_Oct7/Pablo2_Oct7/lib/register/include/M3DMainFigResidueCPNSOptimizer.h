#ifndef M3D_MAIN_FIG_RESIDUE_CPNS_OPTIMIZER_H
#define M3D_MAIN_FIG_RESIDUE_CPNS_OPTIMIZER_H

#include "M3DMainFigResidueCPNSProblem.h"
#include "M3DMainFigureOptimizer.h"
#include <vector>

class M3DMainFigResidueCPNSOptimizer : public M3DMainFigureOptimizer
{
public:
	M3DMainFigResidueCPNSOptimizer(Match *m, M3DCPNSStats *cpnsPtr, int order);
	~M3DMainFigResidueCPNSOptimizer() { };

	bool initialize(M3DObject * referenceObject, int treeIndex);

	virtual M3DObject * createTargetObject(const Vector &x);

protected:
	const double getNumericalDerivativeEpsilon() { return 0.005; }
};

#endif