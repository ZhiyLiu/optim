#ifndef _SCALE_TRANS_SHAPE_SPACE_H_
#define _SCALE_TRANS_SHAPE_SPACE_H_

#include "SimTransShapeSpace.h"

// This class provides a way to generate all possible transformations
// involving scaling and translation only
class ScaleTransShapeSpace : public SimTransShapeSpace
{
	public:

		enum ScaleTransType { Scale, Translation, Both };

		ScaleTransShapeSpace();
		void virtual setType(ScaleTransType t) { type = t; }

		// override these two functions to limit the shape space
		virtual int parameterCount(); // how many parameters define the space 
		virtual SimilarityTransform3D createSimTrans(const Vector &v); 

	private:

		enum ScaleTransType type;

};

#endif /* _SCALE_TRANS_SHAPE_SPACE_H_ */
