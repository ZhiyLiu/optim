#ifndef M3D_MAIN_FIG_RESIDUE_PGA_OPTIMIZER_H
#define M3D_MAIN_FIG_RESIDUE_PGA_OPTIMIZER_H

#include "M3DMainFigResiduePGAProblem.h"
#include "M3DMainFigureOptimizer.h"
#include <vector>


class M3DMainFigResiduePGAOptimizer : public M3DMainFigureOptimizer
{
	public:
		M3DMainFigResiduePGAOptimizer(Match * m, M3DPGAStats * pgaPtr, int order);
		~M3DMainFigResiduePGAOptimizer() { }

		bool initialize(M3DObject * referenceObject, int treeIndex);

		virtual M3DObject * createTargetObject(const Vector & x);

	protected:
		const double getNumericalDerivativeEpsilon() { return 0.005; }
};


#endif

