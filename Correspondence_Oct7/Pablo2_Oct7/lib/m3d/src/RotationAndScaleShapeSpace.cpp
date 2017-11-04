#include "RotationAndScaleShapeSpace.h"

int RotationAndScaleShapeSpace::parameterCount() { 
  return 4;  // Axis of Rotation (2) + Angle(1) + Scale (1) 
}


SimilarityTransform3D RotationAndScaleShapeSpace::createSimTrans(const Vector &v) {

  Vector3D axis(cos(v(0)) * sin(v(1)), sin(v(0)) * sin(v(1)), cos(v(1)));
  double angle(v(2) / 64.0);
  Quat rot;
  rot.setAxisAngle(axis, angle);
  // rot.normalize();
  double scale(exp(v(3)));

  Vector3D trans(0, 0, 0);

  Vector3D newRefCog(refCog);
	  
  // rotate
  rot.rotateVector(newRefCog);

  // scale
  newRefCog *= scale;

  trans = targetCog - newRefCog;

  return SimilarityTransform3D(trans, rot, scale);
}
