#ifndef M3D_TUBE_SIMILARITY_PGA_PROBLEM_H
#define M3D_TUBE_SIMILARITY_PGA_PROBLEM_H

#include "M3DSimilarityPGAProblem.h"
#include "M3DPGAStats.h"
#include <vector>

class M3DTubeSimilarityPGAProblem;

class M3DTubeSimilarityPGAProblemBoundsFunction : public M3DSimilarityPGAProblemBoundsFunction
{
	M3DTubeSimilarityPGAProblem* problem;
public:
	M3DTubeSimilarityPGAProblemBoundsFunction( M3DTubeSimilarityPGAProblem* _problem ) 
		: problem(_problem) {}
	double bound(int i);
	double factor(int i);
};

class M3DTubeSimilarityPGAProblem : public M3DSimilarityPGAProblem
{
public:
    M3DTubeSimilarityPGAProblem();
    ~M3DTubeSimilarityPGAProblem();

    int parameterCount() const;

    M3DObject * createTargetObject(const Vector & x, int & figureId,
		bool predict = false);
	std::vector<double> pgaCoefficients(const Vector &x );

private:
	/**
	 * Whether to optimize over phi or not.
	 * This mode matters only when working
	 * with tubes and using landmark penalty with landmarks along the length
	 * or local regions
	 */
	mutable bool tubePhiMode;
	friend class M3DTubeSimilarityPGAProblemBoundsFunction;
};

#endif

