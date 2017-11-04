#include "M3DPrimitiveCorrector.h"

// Predicts a primitive's position based on a neighbor
// Given the change in the neighbor primitive, we predict how the primitive
// will change
//
// newPredictee - the predicted primitive is returned in this variable
// referencePredictor - the initial state of the neighbor primitive
// referencePredictee - the initial state of the primitive we predict
// newPredictor - the changed neighbor primitive
void M3DPrimitiveCorrector::predictPrimitive(M3DPrimitive & newPredictee,
        const M3DPrimitive & referencePredictor,
        const M3DPrimitive & referencePredictee,
        const M3DPrimitive & newPredictor)
{
    Vector3D linkVector;
    Quat inversePredictorQ,
         newQ,
         linkRotation;
    double scale;

    linkVector = referencePredictee.getX() - referencePredictor.getX();

    inversePredictorQ = referencePredictor.getQ().conj();
    newQ = newPredictor.getQ() * inversePredictorQ * referencePredictee.getQ();
    linkRotation = newPredictor.getQ() * inversePredictorQ;

    linkRotation.rotateVector(linkVector);

    scale = newPredictor.getR() / referencePredictor.getR();

    newPredictee.setX(newPredictor.getX() + linkVector);
    newPredictee.setQ(newQ);
    newPredictee.setR(scale * referencePredictee.getR());
    newPredictee.setTheta(referencePredictee.getTheta());
}

// Weighted average of two primitives (linear interpolation)
// 0 weight for prim1, 1 weight for prim2
//
// avePrim - averaged primitive returned here
// prim1, prim2 - primitives to average
// weight - parameter in [0,1], goes from prim1 to prim2
void M3DPrimitiveCorrector::averagePrimitives(M3DPrimitive & avePrim,
        const M3DPrimitive & prim1,
        const M3DPrimitive & prim2, double weight)
{
    Quat q;

	if(weight == 0.0)
		avePrim = prim2;
	else if(weight == 1.0)
		avePrim = prim1;
	else
	{
	    avePrim.setX(weight * prim1.getX() + (1.0 - weight) * prim2.getX());
		avePrim.setQ(q.slerp(weight, prim1.getQ(), prim2.getQ()));
		avePrim.setR(weight * prim1.getR() + (1.0 - weight) * prim2.getR());
		avePrim.setTheta(weight * prim1.getTheta() +
						(1.0 - weight) * prim2.getTheta());
	}
}

