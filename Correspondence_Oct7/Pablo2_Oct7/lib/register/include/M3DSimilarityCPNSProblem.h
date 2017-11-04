#ifndef M3D_SIMILARITY_CPNS_PROBLEM_H
#define M3D_SIMILARITY_CPNS_PROBLEM_H

#include "M3DMainFigureProblem.h"
#include "M3DCPNSStats.h"
#include <vector>

class M3DSimilarityCPNSProblemBoundsFunction : public BoundsFunction
{
public:
	M3DSimilarityCPNSProblemBoundsFunction();
};

class M3DSimilarityCPNSProblem : public M3DMainFigureProblem
{
public:
	M3DSimilarityCPNSProblem();
	~M3DSimilarityCPNSProblem();

	int parameterCount() const;

	void initialize( Match * _match, M3DObject * referenceObject, int treeIndex);

	virtual double evaluate(const Vector &x);
	virtual Components evaluateTerms( const Vector &x );

	M3DObject * createTargetObject(const Vector & x, int & figureId, bool predict = false);
	virtual std::vector<double> cpnsCoefficients(const Vector &x);

	void setCPNS(M3DCPNSStats * cpnsStat);

	friend class M3DSimilarityCPNSProblemBoundsFunction;

protected: 

	void initializeCPNS();

	M3DCPNSStats * cpns;

#ifdef BINARY
	M3DObject * similarityObject;
#endif
};

#endif