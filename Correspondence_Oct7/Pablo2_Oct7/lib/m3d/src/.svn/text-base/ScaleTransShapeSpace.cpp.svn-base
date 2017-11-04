#include "ScaleTransShapeSpace.h"

ScaleTransShapeSpace::ScaleTransShapeSpace()
{
	type = Both;
}

int ScaleTransShapeSpace::parameterCount()
{
	if (type == Scale)
		return 1;
	if (type == Translation)
		return 3;

	return 4;  // Translation (3) + Scale (1) 
}

SimilarityTransform3D ScaleTransShapeSpace::createSimTrans(const Vector &v)
{
	Vector3D trans;
	double scale;
	const Quat rot;

	switch (type) {
		case Scale:
					scale = exp(v(0));
					break;
		case Translation:
					trans = Vector3D(v(0), v(1), v(2));
					scale = 1.0;
					break;
		case Both:
					trans = Vector3D(v(0), v(1), v(2));
					scale = exp(v(3));
					break;
	}

	return SimilarityTransform3D(trans, rot, scale);
}
