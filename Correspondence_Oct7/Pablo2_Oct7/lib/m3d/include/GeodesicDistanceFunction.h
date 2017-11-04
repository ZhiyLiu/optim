#ifndef _GEODESIC_DISTANCE_FUNCTION_H
#define _GEODESIC_DISTANCE_FUNCTION_H

#include "optima.h"
#include "M3DObject.h"
#include "SimTransShapeSpace.h"
#include <vector>

using std::vector;

class GeodesicDistanceFunction : public Function {

typedef vector<M3DPrimitive*> primList;

public:
    GeodesicDistanceFunction(primList ref, primList target, SimTransShapeSpace* space);
	Vector doOptimization(); // this really invokes conjugate gradient

	virtual double evaluate(const Vector &v);
    primList createCandidateAtoms(const Vector& v);
	SimilarityTransform3D createSimTrans(const Vector &v);
	int transWeight;
    int scaleWeight;	
	int quatWeight;
	int thetaWeight;

protected:
     primList refPrims;
	 primList targetPrims;
	 SimTransShapeSpace* shapeSpace;
};


#endif /* _GEODESIC_DISTANCE_FUNCTION_H */
