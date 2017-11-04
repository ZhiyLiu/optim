#include "SimTransShapeSpace.h"


int SimTransShapeSpace::parameterCount() { 
	return 7;	// Translation (3) + Axis of Rotation in spherical coordinates (2)
				// + Angle(1) + Scale (1) 
}

SimilarityTransform3D SimTransShapeSpace::createSimTrans(const Vector & v) {
	
	Vector3D trans(v(0), v(1), v(2));
	Vector3D axis(cos(v(3)) * sin(v(4)), sin(v(3)) * sin(v(4)), cos(v(4)));
	double angle(v(5) / 64.0);			// AGG: Division by 64 is purely arbitrary
	Quat rot;
	rot.setAxisAngle(axis, angle);
	// rot.normalize();
	double scale(exp(v(6)));
	
	return SimilarityTransform3D(trans, rot, scale);
}

