#ifndef M3D_ADAPTIVE_REGISTRATION_PGA_OPTIMIZER_H
#define M3D_ADAPTIVE_REGISTRATION_PGA_OPTIMIZER_H

#include "M3DAdaptiveRegistrationPGAProblem.h"
#include "M3DRegistrationOptimizer.h"
#include <vector>


class M3DAdaptiveRegistrationPGAOptimizer : public M3DRegistrationOptimizer
{
public:
    M3DAdaptiveRegistrationPGAOptimizer(M3DPGAStats * pgaPtr, int order);
    ~M3DAdaptiveRegistrationPGAOptimizer() { }

private:

};


#endif

