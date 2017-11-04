#include <fstream>
#include <stdio.h>
#include <math.h>
#include "OptimizerBase.h"
#include "M3DMainFigResiduePGAOptimizer.h"

//#define DEBUG

using namespace std;

M3DMainFigResiduePGAOptimizer::M3DMainFigResiduePGAOptimizer(Match * m,
	M3DPGAStats * pgaPtr, int order) : M3DMainFigureOptimizer(m)
{
    boundsFunction = new M3DMainFigResiduePGAProblemBoundsFunction;
    problem = new M3DMainFigResiduePGAProblem(pgaPtr, order);
    numParameters = ((M3DMainFigResiduePGAProblem *) problem)->parameterCount();
}

bool M3DMainFigResiduePGAOptimizer::initialize(M3DObject * referenceObject,
    int treeIndex)
{
	M3DPGAStats * pgaPtr = ((M3DMainFigResiduePGAProblem *)problem)->getPGA();
	int order = ((M3DMainFigResiduePGAProblem *)problem)->getOrder();

	M3DObject * tempObj = referenceObject->assign();
	pgaPtr->applyMeanResidue(tempObj, order);
	delete tempObj;

#ifdef DEBUG
	cout << "treeIndex " << treeIndex << endl;
#endif 

	return M3DMainFigureOptimizer::initialize(referenceObject, treeIndex);
}


M3DObject * M3DMainFigResiduePGAOptimizer::createTargetObject(const Vector & v)
{
	int figureId;

	if (problem == NULL)
		return NULL;

	return problem->createTargetObject(v, figureId, true);
}

