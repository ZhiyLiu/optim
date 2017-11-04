#ifndef _NORMAL_ALIGNED_SHAPE_SPACE_H_
#define _NORMAL_ALIGNED_SHAPE_SPACE_H_

#include "GeodesicDistanceFunction.h"
#include "SimTransShapeSpace.h"


class NormalAlignedShapeSpace : public SimTransShapeSpace
{

  typedef vector<M3DPrimitive*> primList; 

public:

  NormalAlignedShapeSpace(primList ref, primList target, vector<int> idx);  
  // override these two functions to limit the shape space
  virtual int parameterCount(); // how many parameters define the space 
  virtual SimilarityTransform3D createSimTrans(const Vector &v); 

private:

  Vector3D N;                 //Normal vector
  Vector3D refCOG, targetCOG; //COG of 'ref' and 'target'
  Quat     NQuat;             //Quaternion to represent 2D rotation relative to N
};



#endif /* _NORMAL_SHAPE_SPACE_H_ */

