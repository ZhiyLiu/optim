#ifndef _FUNCTION_EXPLORER_
#define _FUNCTION_EXPLORER_

#include "optima.h"

#ifdef OPTIMIZATION_VISUALIZER

class FunctionExplorer
{
public:
	// Explore this function in the plane including point v, 
	// and a vector orthogonal to dir
	static void explorePlane(Function * f, const Vector & v, const Vector & dir) {
		explorePlane(f, v, dir, createOrthogonalVector(dir));
	}
	static void explorePlane(Function * f, const Vector & v, const Vector & dir1, const Vector & dir2);

	static const Vector createOrthogonalVector(const Vector & v);
};

#endif


#endif	/* _FUNCTION_EXPLORER_ */
