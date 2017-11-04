#ifndef VISUALIZATION_H
#define VISUALIZATION_H

#include <iostream>

#include <vtkSmartPointer.h>

#include <vtkCellArray.h>
#include <vtkQuad.h>
#include <vtkLine.h>

#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkProperty.h>

#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderWindow.h>

#include <vtkVersion.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkPointData.h>
#include <vtkUnsignedCharArray.h>
#include <vtkVertexGlyphFilter.h>
#include <vtkLookupTable.h>


#include <vtkCellData.h>
#include <vtkDoubleArray.h>

#include <vtkSphereSource.h>




#include "Vector3D.h"

#include "toolsfunc.h"
//#include "calregularityfeatures.h"


class visualization
{
public:
    visualization();// parameterless default constructor

    //Declaring a constructor. A constructor is similar to a function, but with the following differences: No return type; No return statement.
    //constructor with parameters
    visualization(M3DQuadFigure* quadfig, int interpolationLevel, int side, bool moved, double *quadColor, vtkSmartPointer<vtkRenderer> renderer,
                  double rX, double rY, double rZ, double opacity);

    visualization(M3DQuadFigure* quadfig, int interpolationLevel, int side, double *quadColor, vtkSmartPointer<vtkRenderer> renderer,
                  double rX, double rY, double rZ, double opacity);

    void drawEdgeOfQuads(bool moved);


    void drawQuads(vtkSmartPointer< vtkPoints > hubpos);
    void drawQuadByIndex(vtkSmartPointer< vtkPoints > hubpos, int quadIndex);
    void drawQuadByIndexAndHubpos(vtkSmartPointer< vtkPoints > hubpos, int quadIndex);
    void drawVolume(vtkSmartPointer< vtkPoints > boundaryHubpos, vtkSmartPointer< vtkPoints > skeletalHubpos, int quadIndex);

    void drawEdgeOfQuads_method2(vtkSmartPointer< vtkPoints > hubpos);
    void drawFrameOfQuads(vtkSmartPointer< vtkPoints > hubpos);
    void drawFrameOfQuadByIndex(vtkSmartPointer< vtkPoints > hubpos, int quadIndex);
    void drawSpoke(vtkSmartPointer< vtkPoints > boundaryHubpos, vtkSmartPointer< vtkPoints > skeletalHubpos, int quadIndex);
    void drawPointsOfTheUVCoordinate(VectorQuadPoint quadpoints_u, VectorQuadPoint quadpoints_v);

    void drawEdgeOfQuads_method3(vtkSmartPointer< vtkPoints > hubpos);

    void setColor(double *quadColor);
    void setRenderer(vtkSmartPointer<vtkRenderer> renderer);

    //get all the vetex hub position of subquad.
    vtkSmartPointer< vtkPoints > getSubQuadsPosition(int quadtype);
    void getSubQuadsUVCoordinateUsingDeltaUV();
    void getSubQuadsUVCoordinateNotUsingDeltaUV();

    double lengthofedges(Vector3D point[2]);



    void drawCurveAlongDiscreatPoints(vector<Vector3D> points, double *quadColor);

    void drawCrest(bool up, bool down, bool medial);
    void connectCrestCorrespondenceLines(vector<Vector3D> points_A, vector<Vector3D> points_a);
    void drawCrestCorrespondenceLines(bool crest_up, bool crest_down);
    void connectInteriorToCrestLines(vector<Vector3D> points);
    void drawUpOrDownBoundaryQuadLine_method1(bool up, bool down);
    void getAllSubQuadsVertexesOnBoundary(); // used for test the method getAllSubQuadsVertexesOnBoundary() in alignsrep.cpp
    void drawFrameOfQuadByIndex_m2(vector<Vector3D> points, int quadIndex);
    void drawUpOrDownBoundaryQuadLine_method2(bool up, int index1, bool down, int index2);


    void connectDiagnoal(vtkSmartPointer< vtkPoints > hubpos);

    void drawSubQuadVolume(vtkSmartPointer< vtkPoints > boundaryHubpos, vtkSmartPointer< vtkPoints > skeletalHubpos, int quadIndex,
                                   int subQuadIndex);

    void drawSpoke_2(vtkSmartPointer< vtkPoints > boundaryHubpos, vtkSmartPointer< vtkPoints > skeletalHubpos, int quadNums, int linewith,
                     bool showspoketail);
    void drawSpoke_3(vtkSmartPointer< vtkPoints > boundaryHubpos, vtkSmartPointer< vtkPoints > skeletalHubpos, int quadNums);
    void showPointBall(vtkSmartPointer< vtkPoints > points, double * color, double pointSize, vtkSmartPointer<vtkRenderer> renderer,
                                      double rX, double rY, double rZ, double opacity);

    void MakeLUT(size_t const & tableSize, vtkLookupTable *lut);

private:
    vtkSmartPointer<vtkRenderer> renderer;
    double *quadColor;

    int quadtype;
    int side;
    int interpolationLevel;
    int rowNum;
    int colNum;
    int quadNum;
    int subQuadPointsNum;
    int step;

    M3DQuadFigure* curQuadFig;
    toolsfunc tools;

    VectorQuadPoint quadpoints_u;   //store the u coordinate of all the sub-quad of quad[q].
    VectorQuadPoint quadpoints_v;   //store the v coordinate of all the sub-quad of quad[q].

    vtkSmartPointer< vtkPoints > hubPosition_b;   //all the 24*25 points xyz position on the boundary sheet.
    vtkSmartPointer< vtkPoints > hubPosition_s;   //all the 24*25 points xyz position on the skeletal sheet.

    //For the crest curves.
    vector<Vector3D> crest_b_points;
    vector<Vector3D> up_b_crest_points;
    vector<Vector3D> down_b_crest_points;
    vector<Vector3D> up_b_points;
    vector<Vector3D> down_b_points;

    // For compute the area of up and down boundary quads.
    vector<Vector3D> up_BP;//up boundary points list
    vector<Vector3D> down_BP;//down boundary points list

    double rX;
    double rY;
    double rZ;
    double opacity;

};

#endif // VISUALIZATION_H
