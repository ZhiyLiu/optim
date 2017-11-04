#ifndef M3D_SIMILARITY_ELONGATION_PROBLEM_H
#define M3D_SIMILARITY_ELONGATION_PROBLEM_H

#include "M3DMainFigureProblem.h"

class M3DSimilarityElongationProblem;

class M3DSimilarityElongationProblemBoundsFunction : public BoundsFunction
{
	M3DSimilarityElongationProblem* problem;
public:
	M3DSimilarityElongationProblemBoundsFunction( M3DSimilarityElongationProblem* _problem )
		: problem(_problem) {}
	double bound(int i);
	double factor(int i);
};

class M3DSimilarityElongationProblem : public M3DMainFigureProblem
{
	bool tubePhiMode;
public:
    M3DSimilarityElongationProblem();
    ~M3DSimilarityElongationProblem() {}

	virtual void initialize(Match * m, M3DObject * referenceObject,
		int treeIndex);

    // The vector parameter is evaluated as a 7-tuple:
    // (translation in x,y,z; rotation around x,y,z; scale)
    virtual double evaluate(const Vector & x);

    virtual Components evaluateTerms(const Vector & x);

	// The predict argument is only used in Adaptive Pablo
    virtual M3DObject * createTargetObject(const Vector & x, int & figureId, bool predict = false);

    friend class M3DSimilarityElongationProblemBoundsFunction;

};


#endif

