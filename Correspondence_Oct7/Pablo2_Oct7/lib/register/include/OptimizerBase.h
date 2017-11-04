#ifndef _OPTIMIZER_BASE_H_
#define _OPTIMIZER_BASE_H_

#include "matrix.h"
#include "M3DObject.h"
#include "ConjugateGradientMethod.h"


/*
	This class is the base class for all optimizers.

	Having a common base class makes it easier for the objective function
	visualization tools to reset the optimization starting point.
*/

extern const int MAX_OPTIMIZER_ITERATIONS;

class OptimizerBase
{
public:

	// Is this used somewhere? - rrs
	enum optimizer_t { BoundaryDispl, Deformation, Spoke, MainFigure, Registration,
		Subfigure, SRep };
	//	Others that could be added:  AdaptRegistrPGA, MainFigResiduePGA,
	//	RegistrPGA, RegistrSimilarity, SimilarityElong, and SimilarityPGA

	optimizer_t optimizerType() { return type; }

	OptimizerBase();
	virtual ~OptimizerBase();

	virtual M3DObject * getReferenceObject() = 0;

	const Vector & getLastBestX(int i = 0) { 
		// M3DDeformationOptimizer has one Vector per atom.
		// All other optimizers return a single vector.
		return *(lastBestX[i]);
	}
	int getLastBestXSize() { return lastBestXSize; }

	double getLastBestVal() { return lastBestVal; }
	// Note: M3DDeformationOptimizer has another getLastBestVal() function

	double getLastPenalty() { return lastPenalty; }
	M3DObject * getLastBestObject() { return lastBestObject; }

	// Is the optimizer finished running?
	// M3DBoundaryDisplacementOptimizer has a different version
	virtual bool isFinished();

	void skipStage() { skipped = true; }
	bool stageSkipped() { return skipped; }

	// Perform up to iters iterations of optimization.  If true is
	// returned, then the optimization did something or is finished.
	// A return of false means that either an error occurred or the
	// iteration had no effect.  For example, compare the functions
	// in classes M3DDeformationOptimizer and M3DSubfigureOptimizer.
	// M3DDeformationOptimizer has another performIterations() function.
	virtual bool performIterations(int iters) = 0;

	// Get the target object from the problem for this position in space.
	virtual M3DObject * createTargetObject(const Vector & x) = 0;

	// Get the problem used by the optimizer.
	// BoundaryDisplacementOptimizer returns NULL.
	// The DeformationOptimizer presently returns NULL, but could
	// be changed to return an array of problems, one per atom.
	virtual Function * getProblem() = 0;

	// Get the Differentiable version of the problem 
	// used by the optimizer.
	// BoundaryDisplacementOptimizer returns NULL.
	// The DeformationOptimizer presently returns NULL, but could
	// be changed to return an array of problems, one per atom.
	virtual DifferentiableFunction * getDifferentiableProblem() = 0;

	// Tell the optimizer to start over from a new position in
	// the parameter space
	virtual void setOptimizerPosition(const Vector& x) = 0;

	virtual Vector* projectModel(M3DObject* model);

protected:

	// These variables must be set by the derived classes
	Vector ** lastBestX;	// The best the optimizer has done so far
	int lastBestXSize;		// Length of the lastBestX vector; may be 1
	double lastBestVal;		// Objective function value
	double lastPenalty;		// Geometric penalty value
	M3DObject * lastBestObject;
	optimizer_t type;
	Method * method;	// The method for optimization
	bool skipped;
};

class GeodesicDistanceProblem : public Function {
public:
	// need an optimizer base so that we can call createTargetObject
	GeodesicDistanceProblem(OptimizerBase* o, M3DObject* r) { opt = o; ref=r;}
	virtual double evaluate(const Vector& x);

private:
	M3DObject* ref;
	OptimizerBase* opt;
};

#endif	/* _OPTIMIZER_BASE_H_ */

