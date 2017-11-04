#ifndef M3D_SUBFIGURE_PROBLEM_H
#define M3D_SUBFIGURE_PROBLEM_H

#include "Match.h"
#include "optima.h"

class M3DSubfigureProblem : public Function
{
public:
    M3DSubfigureProblem(Match * _match, M3DObject * _referenceObject,
                        int _figureId, double _penaltyWeight = 0.0,
                        double _constraintsPenaltyWeight = 0.0);

    ~M3DSubfigureProblem();

    Match * getMatch() { return match; }

    M3DObject * getReferenceObject() { return referenceObject; }
    void setPenaltyWeight(double w) { penaltyWeight = w; }
    void setConstraintsPenaltyWeight(double w) { constraintsPenaltyWeight = w; }

    // The vector parameter is evaluated as a 5-tuple:
    // (translation in u,v; rotation in u-v plane; scale; elongation)
    double evaluate(const Vector &x);

    double computePenalty(const Vector &x);

	static const double SUB_TRANSLATION_FACTOR;
	static const double SUB_ROTATION_FACTOR;
	static const double SUB_WIDEN_FACTOR;
	static const double SUB_ELONGATE_FACTOR;
	static const double SUB_HINGE_FACTOR;

private:
    void initialize();

    void applyVector(const Vector &x);


    Match * match;

    M3DObject * referenceObject;
    M3DObject * targetObject;

    int figureId;
    M3DFigureTreeNode * figureTreeNode;

    double penaltyWeight;
    double constraintsPenaltyWeight;

    int numSurfacePoints;
    double * surfacePoints;
};

#endif

