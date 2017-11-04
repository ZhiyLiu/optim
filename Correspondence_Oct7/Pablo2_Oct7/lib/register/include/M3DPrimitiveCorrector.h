#ifndef M3D_PRIMITIVE_CORRECTOR_H
#define M3D_PRIMITIVE_CORRECTOR_H

#include "M3DPrimitive.h"

class M3DPrimitiveCorrector
{
public:
    M3DPrimitiveCorrector() {}

    static void predictPrimitive(M3DPrimitive & newPredictee,
        const M3DPrimitive & referencePredictor, const M3DPrimitive & referencePredictee,
        const M3DPrimitive & newPredictor);

    static void averagePrimitives(M3DPrimitive & avePrim, const M3DPrimitive & prim1,
        const M3DPrimitive & prim2, double weight);
};

#endif

