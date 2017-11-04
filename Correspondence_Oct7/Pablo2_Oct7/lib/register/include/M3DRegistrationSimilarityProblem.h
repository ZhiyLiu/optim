#ifndef M3D_REGISTRATION_SIMILARITY_PROBLEM_H
#define M3D_REGISTRATION_SIMILARITY_PROBLEM_H

#include "M3DRegistrationProblem.h"

class M3DRegistrationSimilarityProblemBoundsFunction : public BoundsFunction
{
public:
	M3DRegistrationSimilarityProblemBoundsFunction();
};

class M3DRegistrationSimilarityProblem : public M3DRegistrationProblem
{
public:
    M3DRegistrationSimilarityProblem();
    ~M3DRegistrationSimilarityProblem() {}

    // The vector parameter is evaluated as a 7-tuple:
    // (translation in x,y,z; rotation around x,y,z; scale)
    virtual double evaluate(const Vector &x);

    virtual M3DObject * createTargetObject(const Vector & x);

    friend class M3DRegistrationSimilarityProblemBoundsFunction;

private:
	void computeTransformation(SimilarityTransform3D &transform, const Vector &x);
};

#endif

