#ifndef M3D_REGISTRATION_SIMILARITY_OPTIMIZER_H
#define M3D_REGISTRATION_SIMILARITY_OPTIMIZER_H

#include "M3DRegistrationOptimizer.h"
#include "M3DRegistrationSimilarityProblem.h"


class M3DRegistrationSimilarityOptimizer: public M3DRegistrationOptimizer
{
public:
    M3DRegistrationSimilarityOptimizer();
    ~M3DRegistrationSimilarityOptimizer() { }
};


#endif

