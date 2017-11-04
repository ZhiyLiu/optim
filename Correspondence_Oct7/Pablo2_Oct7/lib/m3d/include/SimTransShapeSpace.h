#ifndef _SIM_TRANS_SHAPE_SPACE_H_
#define _SIM_TRANS_SHAPE_SPACE_H_


#include "M3DObject.h"
#include "optima.h"


// This [base] class provides a way to generate all possible simtranses.
// A subclass can restrict this.

// Classes derived from this one include RigidShapeSpace, ScaleTransShapeSpace,
// NormalAlignedShapeSpace, and RotationAndScaleShapeSpace.


class SimTransShapeSpace
{
	public:

		// Override these two functions to limit the shape space
		virtual int parameterCount();	// Number of parameters that define the space
		virtual SimilarityTransform3D createSimTrans(const Vector & v); 
};


#endif	/* _SIM_TRANS_SHAPE_SPACE_H_ */

