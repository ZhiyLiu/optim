#ifndef _RIGID_SHAPE_SPACE_H_
#define _RIGID_SHAPE_SPACE_H_

#include "SimTransShapeSpace.h"


// This class provides a way to generate all possible rigid transformations
class RigidShapeSpace : public SimTransShapeSpace
{
public:

  // override these two functions to limit the shape space
  virtual int parameterCount(); // how many parameters define the space 
  virtual SimilarityTransform3D createSimTrans(const Vector &v); 

};


#endif /* _RIGID_SHAPE_SPACE_H_ */

