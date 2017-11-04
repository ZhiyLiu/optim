#ifndef CALCRESTREGULARITYFEATURES_H
#define CALCRESTREGULARITYFEATURES_H


#include "toolsfunc.h"
#include <vtksrepinterpolatecrestspokesquartic.h>
#include <vtksrepinterpolatemedialcrestcurve.h>
#include "vtkreadsrep.h"
#include <vtksrep.h>

#include "calregularityfeatures.h"

#include "time.h"

#include "uvmap.h"



class calcrestregularityfeatures
{
public:
    calcrestregularityfeatures();
    calcrestregularityfeatures(M3DQuadFigure* quadfig, DoubleVec subTVs, int interpolationLevel, int crestAtomNums);

    vtkSmartPointer<vtkSRep> getSrepFig(M3DQuadFigure* quadfig);

    void getSubQuadsPositionOfCrestRegion(vtkSmartPointer< vtkPoints > points_s, vtkSmartPointer< vtkPoints > points_b);

    void calculateCrestEdges_2(vtkSmartPointer< vtkPoints > points_s, vtkSmartPointer< vtkPoints > points_b,
                                                           MatrixType &verEdgeFeatures, MatrixType &horEdgeFeatures);
    double pointsDistance(double * p1, double * p2);
    void calculateCrestAngles_2(vtkSmartPointer< vtkPoints > points_b, MatrixType &angleFeatures);



private:

    DoubleVec subTVs;

    int subQuadVertexNum;
    int step;
    int exteriorAtomNum;

    M3DQuadFigure* quadfig;
    int interpolationLevel;
};

#endif // CALCRESTREGULARITYFEATURES_H
