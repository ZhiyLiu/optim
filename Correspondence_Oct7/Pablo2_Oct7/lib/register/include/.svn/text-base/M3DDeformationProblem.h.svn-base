#ifndef M3D_DEFORM_PROBLEM
#define M3D_DEFORM_PROBLEM

#include "M3DObject.h"
#include "Image3D.h"
#include "optima.h"
#include "M3DDeformationOptimizer.h"
#include "M3DPGAPrimitiveStats.h" 
class M3DDeformationOptimizer;

class M3DDeformationProblem;

class M3DDeformationProblemBoundsFunction : public BoundsFunction
{
protected:
	const M3DDeformationProblem* problem;
public:
	M3DDeformationProblemBoundsFunction( const M3DDeformationProblem* _problem )
		: problem(_problem) {}
	double bound(int i);
	double factor(int i);
};

class Match;

#define DEFORMATION_MAX_RETURN 1.7e308


class M3DDeformationProblem : public Function
{
public:

    M3DDeformationProblem(int appID = 0);
    M3DDeformationProblem(Match * _match, M3DObject * _referenceObject,
          M3DObject * _targetObject,int _figureId,int _atomId,
		  M3DPGAPrimitiveStats * _primPga = NULL,int appID = 0);

    virtual ~M3DDeformationProblem();

    virtual double evaluate(const Vector & x);
	virtual Components evaluateTerms( const Vector &x );
	virtual double evaluateImageMatch(const Vector & x, M3DObject * _referenceObject = NULL);	// Used by preview operation only

    virtual double getGeometricPenalty() { return geomPenalty; }

  void createTargetPrimitive(const Vector &x);

  // apply a transformation, represented by the vector (x), to an atom (prim)
  void applyVector(M3DPrimitive & prim, const Vector & x);

  //deform the mean delta atom (prim) by the coefficient vector (x) in atom PGA space
  void applyVectorPGA(M3DPrimitive & prim, const Vector & x);

  void setPGA(M3DPGAPrimitiveStats * primPgaStats);

  //apply PGA deformation in Atom PGA stage. Including first applying the atom difference pga deformation
  //to interpolated atom and then adding up with the figure residual atom difference
  void doApplyVectorPGA(M3DFigure * currFigPtr, const Vector & x);


protected:

    M3DPGAPrimitiveStats * primPga;
   	M3DFigure * figurePredictor;
    Match * match;

    const M3DObject * referenceObject;
    M3DObject * targetObject;

    int figureId;
    int atomId;

	double geomPenalty;

	int applicationID;
	// The default value of applicationID is 0, which means the optimizer is used by
	// Pablo.  If applicationID = 1, the optimizer is used by VSkelTool, which requires
	// a different method of computing atom penalties.

    friend class M3DDeformationProblemBoundsFunction;
};

#endif

