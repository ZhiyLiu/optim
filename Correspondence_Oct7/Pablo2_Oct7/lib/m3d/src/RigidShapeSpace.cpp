#include "RigidShapeSpace.h"

int RigidShapeSpace::parameterCount() { 
  return 6;  // Translation (3) + Axis of Rotation (2) + Angle(1) 
}


SimilarityTransform3D RigidShapeSpace::createSimTrans(const Vector &v) {

  Vector3D trans(v(0), v(1), v(2));
  Vector3D axis(cos(v(3)) * sin(v(4)), sin(v(3)) * sin(v(4)), cos(v(4)));
  double angle(exp(v(5)) - 1);
  Quat rot;
  rot.setAxisAngle(axis, angle);
  rot.normalize();
  double scale(1);

  return SimilarityTransform3D(trans, rot, scale);
}
