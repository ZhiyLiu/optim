#ifndef M3D_DEFORM_PROBLEM
#define M3D_DEFORM_PROBLEM

#include "M3DObject.h"
#include "Image3D.h"
#include "optima.h"


class M3DSpokeProblemBoundsFunction : public BoundsFunction
{
public:
	M3DSpokeProblemBoundsFunction();
};


class Match;

class M3DSpokeProblem : public Function
{
public:
    M3DSpokeProblem(int appID = 0);
    M3DSpokeProblem(Match * _match, M3DObject * _referenceObject,
          M3DObject * _targetObject,int _figureId,int _atomId,int _spokeId,
		  int appID = 0);

    virtual ~M3DSpokeProblem();

    virtual double evaluate(const Vector & x);
	virtual double evaluateImageMatch(const Vector & x);	// Used by preview operation only

    virtual double getGeometricPenalty() { return geomPenalty; }

    // Applies a transformation, represented by x, to a spoke of a primitive
	void applyTransform(M3DPrimitive & prim, const double x);

protected:
    Match * match;

    const M3DObject * referenceObject;
    M3DObject * targetObject;

    int figureId;
    int atomId;
	int spokeId;

	double geomPenalty;

	int applicationID;
	// The default value of applicationID is 0, which means the optimizer is used by
	// Pablo.  If applicationID = 1, the optimizer is used by VSkelTool, which requires
	// a different method of computing atom penalties.

    friend class M3DSpokeProblemBoundsFunction;
	friend class M3DSpokeOptimizer;
};

#endif

