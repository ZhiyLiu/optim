#ifndef M3D_QUAD_INTERPOLATOR_H
#define M3D_QUAD_INTERPOLATOR_H

#include "M3DQuadFigure.h"
//#include "Bezier1D.h"
#include "Bezier2D.h"

class M3DQuadInterpolator
{
public:
	M3DQuadInterpolator() { bez2D = NULL; radiusBez2D = NULL; }

	M3DQuadInterpolator(
		M3DPrimitive * p1, M3DPrimitive * p2,
		M3DPrimitive * p3, M3DPrimitive * p4,
		int xmesh = 0, int ymesh = 0);

	M3DQuadInterpolator(
		const M3DQuadFigure &figure, int startPrimId,
		int xmesh = 0, int ymesh = 0);

	~M3DQuadInterpolator();

	bool init(const M3DQuadFigure &figure, int startPrimId,
			  int xmesh = 0, int ymesh = 0);

	Bezier2D * getBez2D() { return bez2D; }

	// Given a quad of primitives and a u and v, interpolates primitive
	// along the interpolated, parameterized medial sheet (Bezier patch)
	M3DPrimitive * interpolatePrim(double u, double v);

	M3DPrimitive * getPrimitive(int i, int j);

private:
	Bezier2D * bez2D;
    Bezier2D * radiusBez2D;

	// Initializes the bezier patch for the interpolated medial sheet
	void initMedialPatch(M3DPrimitive * prim1, M3DPrimitive * prim2,
                         M3DPrimitive * prim3, M3DPrimitive * prim4,
                         int xmesh, int ymesh);
};

#endif


