#include "OptimizerBase.h"
#include "LogManager.h"
#include "Tuning.h"

#ifndef BINARY
#endif /* NOT BINARY */


const int MAX_OPTIMIZER_ITERATIONS = 1 << 16;


// Default constructor
OptimizerBase::OptimizerBase()
{
	lastBestObject = NULL;
	lastPenalty = 0.0;
	lastBestVal = 0.0;
	lastBestX = NULL;
	lastBestXSize = 0;
	method = NULL;
	skipped = false;
}

OptimizerBase::~OptimizerBase()
{
	int i;

	if (lastBestObject != NULL)
		delete lastBestObject;

	if (method != NULL)
		delete method;

    if (lastBestX != NULL)
    {
        for (i = 0; i < lastBestXSize; i++) 
        {
            if (lastBestX[i] != NULL)
                delete lastBestX[i];
        }
        delete [] lastBestX;
    }
}

bool OptimizerBase::isFinished()
{
	if (! method)
		return false;
	if (skipped)
		return true;
	if (method->methodType() == Method::ConjugateGradient)
		return ((ConjugateGradientMethod *) method)->isFinished();
	else
		return true;
}

Vector* OptimizerBase::projectModel(M3DObject* model) {
#ifdef OPTIMIZATION_VISUALIZER
	if (model) {
		
		GeodesicDistanceProblem prob(this, model);
		NumericalFunction wrapper(prob, 1e-5);
		if (lastBestX && lastBestX[0]) {
			int cat = globalLogManager.getCurrentCategory();
			ConjugateGradientMethod opt(wrapper, 0 * (*lastBestX[0]));
			opt.setBrentLinearSearchBoundFactor( tuningWt(BrentLinearSearchBoundFactor) );
			while (!opt.isFinished() ) {
				opt.performIteration();
			}
			globalLogManager.setCategory(cat);
			return new Vector(opt.getBestEverX());
		} // lastBestX exists
	} // model exists
#endif
	return 0;
}

double GeodesicDistanceProblem::evaluate(const Vector& x) {
	if (opt && ref) {
		M3DObject* target = opt->createTargetObject(x);
		if (target) {
		    double distance = 0;
#ifdef BINARY
			int figureId = (int) tuningWt(BpFigureId);
			M3DFigure* targetFig = target->getFigurePtr(figureId);
			M3DFigure* refFig = ref->getFigurePtr(figureId);
			if (targetFig && refFig) {
			  distance = refFig->dist2FromFigure(targetFig);
            }
#else /*NOT BINARY*/
            // I don't know where to get FigureID from
			// but it shouldn't matter 
			int figureId = 0; // This is wrong for multifigure
			M3DFigure* targetFig = target->getFigurePtr(figureId);
			M3DFigure* refFig = ref->getFigurePtr(figureId);
			double atomDist = 0;
			if (targetFig && refFig) {
			  int a = targetFig->getPrimitiveCount();
              int b = refFig->getPrimitiveCount();
			  int numPrims = a < b ? a : b;
			  for (int i = 0; i<numPrims; i++) {
			    //atomDist = geo.atomDistance(targetFig->getPrimitivePtr(i), refFig->getPrimitivePtr(i));
			    targetFig->getPrimitivePtr(i)->atomDistance(refFig->getPrimitivePtr(i));
				distance += (atomDist * atomDist);
		      }
            }
#endif /*BINARY*/

            delete target;
			return distance;
		} // target exists
	} // optimizer and reference model exist
	return -1;
}
