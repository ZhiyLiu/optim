#ifndef M3D_SIMILARITY_PGA_PROBLEM_H
#define M3D_SIMILARITY_PGA_PROBLEM_H

#include "M3DMainFigureProblem.h"
#include "M3DPGAStats.h"
#include <vector>

class M3DSimilarityPGAProblemBoundsFunction : public BoundsFunction
{
public:
	M3DSimilarityPGAProblemBoundsFunction();
};

class M3DSimilarityPGAProblem : public M3DMainFigureProblem
{
public:
    M3DSimilarityPGAProblem();
    ~M3DSimilarityPGAProblem();

    int parameterCount() const;

    void initialize(Match * _match, M3DObject * referenceObject,
        int treeIndex);

    // The vector parameter is evaluated as a (7+k)-tuple:
    // (translation in x,y,z; rotation around x,y,z; scale;
	//  then the k most significant eigenmodes of deformation)
    virtual double evaluate(const Vector &x);
	virtual Components evaluateTerms( const Vector &x );

    M3DObject * createTargetObject(const Vector & x, int & figureId,
		bool predict = false);
	virtual std::vector<double> pgaCoefficients(const Vector &x );

	void setPGA(M3DPGAStats * pgaStat);

    friend class M3DSimilarityPGAProblemBoundsFunction;

protected:

    void initializePGA();

	M3DPGAStats * pga;
#ifdef BINARY
	M3DObject * similarityObject;
#endif
};

#endif

