#ifndef M3D_ADAPTIVE_REGISTRATION_PGA_PROBLEM_H
#define M3D_ADAPTIVE_REGISTRATION_PGA_PROBLEM_H

#include "M3DRegistrationProblem.h"
#include "M3DPGAStats.h"
#include <vector>


class M3DAdaptiveRegistrationPGAProblemBoundsFunction : public BoundsFunction
{
public:
	M3DAdaptiveRegistrationPGAProblemBoundsFunction();
};

class M3DAdaptiveRegistrationPGAProblem : public M3DRegistrationProblem
{
public:
    M3DAdaptiveRegistrationPGAProblem(M3DPGAStats * pgaPtr, int _order);
    ~M3DAdaptiveRegistrationPGAProblem() {}

    int parameterCount() const;

    // The vector parameter is the weights for pga components
    virtual double evaluate(const Vector &x);
    virtual M3DObject * createTargetObject(const Vector & x);

    friend class M3DAdaptiveRegistrationPGAProblemBoundsFunction;

protected:
	M3DPGAStats * pga;
	int order;
};


#endif

