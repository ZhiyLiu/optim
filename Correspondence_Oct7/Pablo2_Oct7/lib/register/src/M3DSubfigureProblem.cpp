#include <math.h>
#include "M3DSubfigureProblem.h"
#include "M3DSubfigureTransformation.h"
#include "Trackball.h"

const double SUBFIGURE_MAX_RETURN = 1.7e308;
const double SUBFIGURE_PARAMETER_SIZE = 5;

const double M3DSubfigureProblem::SUB_TRANSLATION_FACTOR = 1.0;
const double M3DSubfigureProblem::SUB_ROTATION_FACTOR = 2.0;
const double M3DSubfigureProblem::SUB_WIDEN_FACTOR = 4.0;
const double M3DSubfigureProblem::SUB_ELONGATE_FACTOR = 4.0;
const double M3DSubfigureProblem::SUB_HINGE_FACTOR = 8.0;
//#define DISPLAY_PENALTIES


using namespace std;

M3DSubfigureProblem::M3DSubfigureProblem(Match * _match,
                                         M3DObject * _referenceObject,
                                         int _figureId, double _penaltyWeight,
                                         double _constraintsPenaltyWeight)
{
    match = _match;
    figureId = _figureId;

    if(_referenceObject != NULL)
        referenceObject = _referenceObject->assign();
    else
        referenceObject = NULL;

    penaltyWeight = _penaltyWeight;
    constraintsPenaltyWeight = _constraintsPenaltyWeight;

    figureTreeNode = NULL;
    targetObject = NULL;

    numSurfacePoints = 0;
    surfacePoints = NULL;

    initialize();
#ifdef DISPLAY_PENALTIES
    std::cout << "Subfigure penalties:  \n";
#endif
}

M3DSubfigureProblem::~M3DSubfigureProblem()
{
    if(targetObject != NULL)
        delete targetObject;
}

double M3DSubfigureProblem::evaluate(const Vector &x)
{
    double matchVal;
    double penaltyVal;
    double totalPenalty;


    evaluationCost += 1;

    if(match == NULL)
        return SUBFIGURE_MAX_RETURN;

    if(x.size() != SUBFIGURE_PARAMETER_SIZE) 
        return SUBFIGURE_MAX_RETURN;

    applyVector(x);
#ifdef NO_MATCH_PENALTY
        matchVal = 0.0;
#else
        matchVal = -match->computeSubfigureMatch(targetObject, figureId);
#endif
    //matchVal = -match->computeSubfigureMatch(targetObject, figureId);
    penaltyVal = computePenalty(x);

    totalPenalty = penaltyVal - matchVal;
#ifdef DISPLAY_PENALTIES
    std::cout << "   " << penaltyVal << " model - " << matchVal << " match = "
		<< totalPenalty << std::endl;
#endif
    return totalPenalty;
}

double M3DSubfigureProblem::computePenalty(const Vector &x)
{
    if(x.size() != SUBFIGURE_PARAMETER_SIZE) 
        return SUBFIGURE_MAX_RETURN;

    return 0.0;
}

void M3DSubfigureProblem::initialize()
{
    M3DObject * object;
    M3DQuadFigure * figure;

    if(match == NULL)
        return;

    object = match->getReferenceObject();
    if(object == NULL)
        return;

    figure = dynamic_cast<M3DQuadFigure*>( (object->getFigurePtr(figureId)));

    if(figure == NULL)
        return;
}

void M3DSubfigureProblem::applyVector(const Vector &x)
{
    M3DSubfigureTransformation transform;

    if(referenceObject == NULL)
        return;

    if(targetObject != NULL)
        delete targetObject;

    targetObject = referenceObject->assign();
    figureTreeNode = targetObject->getFigureTreeNode(figureId);
    transform.init(targetObject, figureTreeNode);

    transform.translate(SUB_TRANSLATION_FACTOR*x(0), SUB_TRANSLATION_FACTOR*x(1));
    transform.rotate(SUB_ROTATION_FACTOR*x(2));
    transform.widen(pow(LOG_BASE, SUB_WIDEN_FACTOR*x(3)));			//(exp(x(3)));
    transform.elongate(pow(LOG_BASE, SUB_ELONGATE_FACTOR*x(4)));	//(exp(x(4)));
}

