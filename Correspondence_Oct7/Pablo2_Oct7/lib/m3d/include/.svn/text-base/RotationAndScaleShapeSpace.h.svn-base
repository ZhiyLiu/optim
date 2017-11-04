#ifndef _ROTATION_AND_SCALE_SHAPE_SPACE_H_
#define _ROTATION_AND_SCALE_SHAPE_SPACE_H_

#include "SimTransShapeSpace.h"

// This class provides a way to generate all transformations with
// a fixed translation, by varying rotation and scale
class RotationAndScaleShapeSpace : public SimTransShapeSpace {
public:
  RotationAndScaleShapeSpace() {
    targetCog.set(0,0,0);
	refCog.set(0,0,0);
  }


  void setTargetCog(const Vector3D& v) {
    targetCog = (v);
  }

  void setRefCog(const Vector3D& v) {
    refCog = (v);
  }

  // override these two functions to limit the shape space
  virtual int parameterCount(); // how many parameters define the space 
  virtual SimilarityTransform3D createSimTrans(const Vector &v); 

protected:


  Vector3D refCog;
  Vector3D targetCog;


};

#endif /* _ROTATION_AND_SCALE_SHAPE_SPACE_H_ */
