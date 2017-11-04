#ifndef CALREGULARITYFEATURES_H
#define CALREGULARITYFEATURES_H


#include <iostream>

#include <vtkSmartPointer.h>

#include <vtkCellArray.h>
#include <vtkQuad.h>
#include <vtkLine.h>
#include <vtkLineSource.h>

#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkProperty.h>

#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderWindow.h>

#include "P3DControl.h"
#include "ControlParms.h"

#include "Vector3D.h"

#include "toolsfunc.h"

#include "uvmap.h"
#include "quadfigattribution.h"

typedef vnl_matrix<double> MatrixType;

class calregularityfeatures
{
public:
    calregularityfeatures();

    calregularityfeatures(M3DQuadFigure *quadFig,int interpolationLevel, int side, DoubleVec subTVs);

    ~calregularityfeatures();

    double lengthofedges(Vector3D point[2]);


    //get all the vetex hub position of subquad.
    vtkSmartPointer< vtkPoints > getSubQuadsPosition(int quadtype);
    void getSubQuadVertexesPosition(vtkSmartPointer< vtkPoints > bPoints, vtkSmartPointer< vtkPoints > sPoints);


    void calculateAnglesForStandardSide(vtkSmartPointer< vtkPoints > hubpos, MatrixType &angleFeatureMatrix, int quadType);
    void calculateEdges_method3(vtkSmartPointer< vtkPoints > hubpos, MatrixType &horEdgeFeatureMatrix,
                                                       MatrixType &verEdgeFeatureMatrix, int quadType);

    double pointsDistance(double * p1, double * p2);

    void calculateEdgesLength_3(vector<vtkSmartPointer< vtkPoints > > hubposition, MatrixType &horEdgeFeatureMatrix,
                                                     MatrixType &verEdgeFeatureMatrix, int quadType);


private:

    int side;

    DoubleVec subUVs;

    int interpolationLevel;
    M3DQuadFigure *curQuadFig;

    int quadNum;
    int subQuadVetexNum;
    int step;

};

#endif // CALREGULARITYFEATURES_H
