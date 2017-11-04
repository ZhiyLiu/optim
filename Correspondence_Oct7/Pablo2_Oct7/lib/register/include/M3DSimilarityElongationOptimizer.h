#ifndef M3D_SIMILARITY_ELONGATION_OPTIMIZER_H
#define M3D_SIMILARITY_ELONGATION_OPTIMIZER_H

#include "M3DMainFigureOptimizer.h"
#include "M3DSimilarityElongationProblem.h"


class M3DSimilarityElongationOptimizer: public M3DMainFigureOptimizer
{
public:
    M3DSimilarityElongationOptimizer(Match * m, int count = 8);
	~M3DSimilarityElongationOptimizer() { }

protected:
	const double getNumericalDerivativeEpsilon() { return 1e-7; }
};


#endif

