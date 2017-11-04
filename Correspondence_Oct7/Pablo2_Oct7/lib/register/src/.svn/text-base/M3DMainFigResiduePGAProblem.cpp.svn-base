#include <math.h>
#include <iostream>
#include "M3DMainFigResiduePGAProblem.h"
#include "Tuning.h"

extern int globalVerbosity;


//#define CHECK_SIMILARITY_TRANSFORM /*JJ : to print the simialrity transforms */


using namespace std;


M3DMainFigResiduePGAProblem::M3DMainFigResiduePGAProblem(M3DPGAStats * pgaPtr, int _order)
{
	pga = pgaPtr;
	order = _order;
}

int M3DMainFigResiduePGAProblem::parameterCount() const
{
	return pga->getPGDataPtr(order)->numPGs;
}

double M3DMainFigResiduePGAProblem::evaluate(const Vector & x)
{
	double geomPenalty;
	double totalPenalty;
    M3DObject * targetObject;
    int figureId;

    evaluationCost += 1;

    if (match == NULL)
        return MAIN_FIGURE_MAX_RETURN;

	// Prediction is set.
	targetObject = createTargetObject(x, figureId, true);

	std::vector<double> vals;
	for(int i = 0; i != x.size(); ++i ) {
		vals.push_back(x(i)*tuningWt(FigurePGAFactor));
	}
	totalPenalty = match->computeFigureMatchAndPenalties(targetObject,
		figureId, NULL, vals, geomPenalty);
    lastPenalty = geomPenalty;		// Save for future queries

    delete targetObject;

#ifdef OPTIMIZATION_VISUALIZER
	copyResultsToVisualizer(&x, totalPenalty);
#endif

    return totalPenalty;
}

M3DObject * M3DMainFigResiduePGAProblem::createTargetObject(const Vector & x, int & figureId, bool predict)
{
    M3DFigureTreeNode * treeNode;

    if (refObject == NULL || x.size() != parameterCount())
        return NULL;

    treeNode = refObject->getFigureTreeRoot(tree_index);
    if(treeNode == NULL){
		std::cout<<"referenceObject in createTargetObject() is NULL"<<endl;
        return NULL;
	}

    figureId = treeNode->getFigureId();

    M3DObject * targetObject = refObject->assign();

    std::vector<double> vals;
    for(int i = 0; i < x.size(); i++)
        vals.push_back(x(i)*tuningWt(FigurePGAFactor));

	if(!pga->doPGADeform(targetObject, vals, order, predict)){
		std::cout<<"Error in pga deformation at M3DMainFigResiduePGAProblem"<<endl;
		return NULL;
	};

	// Apply initial transformation
	/*
		Ja-Yeon : 3/10/2005 
		targetObject should be selected to apply transformation.
		For display purpose, it should be unselected again.
	*/
	targetObject->select();
#ifdef UNFLIPPED
	Vector3D center(0.0, 1.0, 0.0);
#else
	Vector3D center(0.0, 0.0, 0.0);
#endif
	Quat qtmp = mainFigProblemSimTransform.getRotation(); 
	targetObject->rotateBy(qtmp, center);
    targetObject->scaleBy(mainFigProblemSimTransform.getScale(), center);
	targetObject->translateBy(mainFigProblemSimTransform.getTranslation());

#ifdef CHECK_SIMILARITY_TRANSFORM
	std::cout<<"order "<<order<<std::endl;

	std::cout<<"InitialTransform in M3DMainFigResiduePGAProblem::createTargetObject()" <<std::endl;
	initialTransform.getRotation().print();
	initialTransform.getTranslation().print();
	std::cout<<"scale "<<initialTransform.getScale()<<std::endl;
#endif

	// Apply hand placement.
	SimilarityTransform3D handTransform = pga->getHandPlacementTransform();
	qtmp = handTransform.getRotation();
	targetObject->rotateBy(qtmp, center);
    targetObject->scaleBy(handTransform.getScale(), center);
    targetObject->translateBy(handTransform.getTranslation());

#ifdef CHECK_SIMILARITY_TRANSFORM
	std::cout<<"handTransform in M3DMainFigResiduePGAProblem::createTargetObject() "<<std::endl;

	handTransform.getRotation().print();
	handTransform.getTranslation().print();
	std::cout<<"scale "<<handTransform.getScale()<<std::endl;
#endif

	// Apply best similarity from object_ensemble (or M3DRegistrationPGAOptimizer)...
	// JJ : 05/2005 This xform must be identity now because
	//      xform is not optimized at object ensemble stage.
	SimilarityTransform3D bestEnsembleTrans = pga->getBestTransform();
	qtmp = bestEnsembleTrans.getRotation();
	targetObject->rotateBy(qtmp, center);
    targetObject->scaleBy(bestEnsembleTrans.getScale(), center);
    targetObject->translateBy(bestEnsembleTrans.getTranslation());

#ifdef CHECK_SIMILARITY_TRANSFORM
	std::cout<<"bestEnsembleTrans in M3DMainFigResiduePGAProblem::createTargetObject()"<<std::endl;

	bestEnsembleTrans.getRotation().print();
	bestEnsembleTrans.getTranslation().print();
	std::cout<<"scale "<<bestEnsembleTrans.getScale()<<std::endl;
#endif

	/*
		Ja-Yeon : 3/10/2005 
		Unselect again except the figure of figureId
	*/
	targetObject->deselect();
	M3DFigure *figurePtr = targetObject->getFigurePtr(figureId);
	figurePtr->select();


	//REB REB REB Eli Robbie JVS: Need bone transformation here too!!!
	//This is for soft tissue.  If bone, then let's make sure the bone tranform is identity.


    return targetObject;
}

M3DMainFigResiduePGAProblemBoundsFunction::M3DMainFigResiduePGAProblemBoundsFunction()
{
    bounds	= Vector(1, tuningWt(FigurePGABound));
	factors	= Vector(1, tuningWt(FigurePGAFactor));
}

