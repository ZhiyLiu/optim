// dibyendu

#ifndef M3D_SREP_PROBLEM
#define M3D_SREP_PROBLEM

#include "M3DObject.h"
#include "Image3D.h"
#include "optima.h"
#include "M3DSRepOptimizer.h"
#include "M3DPGAPrimitiveStats.h" 

class M3DSRepOptimizer;

class M3DSRepProblem;

class M3DSRepProblemBoundsFunction : public BoundsFunction
{
protected:
	const M3DSRepProblem* problem;
public:
	M3DSRepProblemBoundsFunction( const M3DSRepProblem* _problem )
		: problem(_problem) {}
	double bound(int i);
	double factor(int i);
};

class Match;

#define SREP_MAX_RETURN 1.7e308


class M3DSRepProblem : public Function
{
public:

    M3DSRepProblem(int appID = 0);

	// create an M3DSRepProblem for a particular spoke
    M3DSRepProblem(Match * _match, M3DObject * _referenceObject,
          M3DObject * _targetObject,int _figureId,int _atomId, int _spokeId, int appID);

	// create an M3DSRepProblem for a complete atom
    M3DSRepProblem(Match * _match, M3DObject * _referenceObject,
          M3DObject * _targetObject,int _figureId,int _atomId, int appID);

    virtual ~M3DSRepProblem();

    virtual double evaluate(const Vector & x);
	virtual Components evaluateTerms( const Vector &x );
	virtual double evaluateImageMatch(const Vector & x, M3DObject * _referenceObject = NULL);	// Used by preview operation only

    virtual double getGeometricPenalty() { return geomPenalty; }

	// creates a target primitive, an atom, by changing one spoke at a time
	// calls the applyVector() function to apply the spoke transformation to the atom
	void createTargetPrimitive(const Vector &x);

	// apply a transformation, represented by the vector (x), to a spoke of an atom (prim)
	void applyVector(M3DPrimitive & prim, const Vector & x);

	// update the target primitive (atom) after both the spokes have been transformed
	void updateTargetPrimitive(const Vector &x);

	// apply a transformation, represented by the vector (x), to an entire atom (prim)
	void applyVectorAtom(M3DPrimitive & prim, const Vector & x);

protected:

   	M3DFigure * figurePredictor;
    Match * match;

    const M3DObject * referenceObject;
    M3DObject * targetObject;

    int figureId ;
    int atomId ;
	int spokeId ;

	double geomPenalty;

	int applicationID;
	// The default value of applicationID is 0, which means the optimizer is used by
	// Pablo.  If applicationID = 1, the optimizer is used by VSkelTool, which requires
	// a different method of computing atom penalties.

    friend class M3DSRepProblemBoundsFunction;
};

#endif

