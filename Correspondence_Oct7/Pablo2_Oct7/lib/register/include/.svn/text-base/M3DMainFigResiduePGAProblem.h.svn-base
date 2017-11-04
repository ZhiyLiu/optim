#ifndef M3D_MAIN_FIG_RESIDUE_PGA_PROBLEM_H
#define M3D_MAIN_FIG_RESIDUE_PGA_PROBLEM_H

#include "M3DMainFigureProblem.h"
#include "M3DPGAStats.h"
#include <vector>


class M3DMainFigResiduePGAProblemBoundsFunction : public BoundsFunction
{
public:
	M3DMainFigResiduePGAProblemBoundsFunction();
};

class M3DMainFigResiduePGAProblem : public M3DMainFigureProblem
{
public:
    M3DMainFigResiduePGAProblem(M3DPGAStats * pgaPtr, int _order);
    ~M3DMainFigResiduePGAProblem() {}

    int parameterCount() const;

    // The vector parameter is evaluated as a (7+k)-tuple:
    // (translation in x,y,z; rotation around x,y,z; scale;
	//  then the k most significant eigenmodes of deformation)
    virtual double evaluate(const Vector &x);

    virtual M3DObject * createTargetObject(const Vector & x, int & figureId, bool predict);

	M3DPGAStats * const getPGA() { return pga; }
	int const getOrder() { return order; }

    friend class M3DMainFigResiduePGAProblemBoundsFunction;

private:
    M3DPGAStats * pga;
	int order;
};


#endif

