#ifndef M3D_REGISTRATION_PGA_PROBLEM_H
#define M3D_REGISTRATION_PGA_PROBLEM_H

#include "M3DRegistrationProblem.h"
#include "M3DPGAStats.h"

class M3DRegistrationPGAProblemBoundsFunction : public BoundsFunction
{
public:
	M3DRegistrationPGAProblemBoundsFunction();
};

class M3DRegistrationPGAProblem : public M3DRegistrationProblem
{
public:
    M3DRegistrationPGAProblem();
    ~M3DRegistrationPGAProblem() {}

    int parameterCount() const;

    // The vector parameter is evaluated as a (7+k)-tuple:
    // (translation in x,y,z; rotation around x,y,z; scale;
	//  then the k most significant eigenmodes of deformation)
    virtual double evaluate(const Vector &x);
    virtual M3DObject * createTargetObject(const Vector & x);

	void setPGA(M3DPGAStats * pgaPtr, int _order) { pga = pgaPtr; order = _order; } 

    friend class M3DRegistrationPGAProblemBoundsFunction;

protected:
	M3DPGAStats * pga;
	int order;
};

#endif

