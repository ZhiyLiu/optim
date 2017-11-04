#ifndef M3D_SIMILARITY_PGA_OPTIMIZER_H
#define M3D_SIMILARITY_PGA_OPTIMIZER_H

#include "M3DSimilarityPGAProblem.h"
#include "M3DMainFigureOptimizer.h"
#include <vector>


class M3DSimilarityPGAOptimizer : public M3DMainFigureOptimizer
{
public:

    M3DSimilarityPGAOptimizer(Match * m, M3DPGAStats * pga);
	~M3DSimilarityPGAOptimizer() { }

protected:
	const double getNumericalDerivativeEpsilon() { return 1e-7; }
};


#endif

