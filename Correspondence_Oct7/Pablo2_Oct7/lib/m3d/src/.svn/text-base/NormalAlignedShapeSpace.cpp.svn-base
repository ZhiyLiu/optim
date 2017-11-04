#include <fstream>
#include <iostream>

#include "M3DObject.h"
#include "M3DObjectFile.h"
#include "NormalAlignedShapeSpace.h"

using std::cout;
using std::endl;
using std::vector;
using std::ios;


NormalAlignedShapeSpace::NormalAlignedShapeSpace(primList ref, primList target, vector<int> idx) {
	
  //Move atoms in Vector(M3DPrimitive*) to M3DPrimitive** 
  //to use ::atomAverage() member function in M3DQuad/TubePrimitive classes
  int numAtoms = idx.size();
  M3DPrimitive** refList    = new M3DPrimitive*[numAtoms];
  M3DPrimitive** targetList = new M3DPrimitive*[numAtoms];

  primList::iterator refIter = ref.begin();
  int i=0;
  for(; refIter != ref.end(); refIter++){
    if (*refIter){
      refList[i] = (*refIter);
      i++;
    }
  } 

  primList::iterator targetIter = target.begin();
  i=0;
  for(; targetIter != target.end(); targetIter++){
    if (*targetIter){
      targetList[i] = (*targetIter);
      i++;
    }
  }

  M3DPrimitive *refMean, *targetMean;
  /* // Create new copies of same type (quad/tube)
     refMean = refList[0]->copyPtr();
     targetMean = targetList[0]->copyPtr();

     refMean->atomAverage(numAtoms, refList);
     targetMean->atomAverage(numAtoms, targetList);
     */
  //debug purpose, normal and rotation should be all about the current atom, instead of the mean atom
  refMean = refList[numAtoms-1]->copyPtr();
  targetMean = targetList[numAtoms-1]->copyPtr();

  //printf("Atom i is used as the COG!\n");

  // 1) Get N that will be an axis of rotation for later further alignment 
  N = targetMean->getN();

  // 2) Get COG of 'ref' and 'target'
  targetCOG = targetMean->getX();
  refCOG = refMean->getX();


  // 3) Get axis of rotation and angle to align normal of mean of 'ref' to that of mean  of 'target'
  Vector3D temp = refMean->getN();
  Vector3D axis = temp.cross(targetMean->getN());
  double angle = acos(temp*(targetMean->getN()));
  //printf("angle = %f\n",angle);

  NQuat.setAxisAngle(axis, angle);
  NQuat.normalize();

  //clean
  delete [] targetList;
  delete [] refList;
  delete refMean;
  delete targetMean;
}


int NormalAlignedShapeSpace::parameterCount() { 
	return 1;  // Angle(1) 
}


SimilarityTransform3D NormalAlignedShapeSpace::createSimTrans(const Vector &v) {

	//2D Rotation with respect to N
	Vector3D axis(N.getX(), N.getY(), N.getZ());
	//double angle(exp(v(0)) - 1);
	double angle(v(0)/64.0);
	Quat temp;
	temp.setAxisAngle(axis, angle);
	temp.normalize();

	//Rotation combined to align N first, then rotate relative to origin
	Quat rot;
//	rot = temp*NQuat;
    rot = NQuat; //the rotation within the tagent plane is not stable to use

	//Translation needs to be recomputed after rotation relative to origin applied
	Vector3D rotRefCOG(refCOG);
	rot.rotateVector(rotRefCOG);

	Vector3D trans;
	trans = targetCOG - rotRefCOG;
	double scale(1);

	Vector3D refCOGClone(refCOG);
	rot.rotateVector(refCOGClone);
	refCOGClone = refCOGClone+trans;
	/* cout << "**** Debug rotated " << endl << "    "; refCOG.print();
	   cout << endl <<  "to: " << endl << "    "; refCOGClone.print(); 
	   cout << endl <<  "should be: " << endl << "    "; targetCOG.print();  
	   */  return SimilarityTransform3D(trans, rot, scale);
}
