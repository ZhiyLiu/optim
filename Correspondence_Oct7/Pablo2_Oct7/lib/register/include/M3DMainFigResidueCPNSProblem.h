#ifndef M3D_MAIN_FIG_RESIDUE_CPNS_PROBLEM_H
#define M3D_MAIN_FIG_RESIDUE_CPNS_PROBLEM_H

#include "M3DMainFigureProblem.h"
#include "M3DCPNSStats.h"
#include <vector>

class M3DMainFigResidueCPNSProblemBoundsFunction : public BoundsFunction
{
public:
	M3DMainFigResidueCPNSProblemBoundsFunction();
};

class M3DMainFigResidueCPNSProblem : public M3DMainFigureProblem
{
public:
	M3DMainFigResidueCPNSProblem(M3DCPNSStats * cpnsPtr, int _order);
	~M3DMainFigResidueCPNSProblem() {}

	int parameterCount() const;

	virtual double evaluate(const Vector &x);
	virtual M3DObject * createTargetObject(const Vector &x, int &figureId, bool predict);

	M3DCPNSStats * const getCPNS() { return cpns; }
	int const getOrder() {return order;}

	friend class M3DMainFigResidueCPNSProblemBoundsFunction;

private:
	M3DCPNSStats *cpns;
	int order;
};

#endif