#ifndef VISUALIZESREPSURFACE_H
#define VISUALIZESREPSURFACE_H

#include "toolsfunc.h"

#include <vtksrepinterpolatecrestspokesquartic.h>
#include "vtkreadsrep.h"
#include <vtksrep.h>
#include "visualizecrest.h"

#include "slidestandardspokes.h"
#include "slidecrestspokes.h"




class visualizesrepsurface
{
public:
    visualizesrepsurface();
    visualizesrepsurface(double rX, double rY, double rZ, double opacity);

    typedef vnl_matrix<double> MatrixType;

    vtkSmartPointer< vtkPoints > getSurfacePointsSet(const char* srepfilename, int interpolationLevel);
    void getInterpolateUVCoordinate(VectorQuadPoint &interpolatedU, VectorQuadPoint &interpolatedV, int rowNum, int colNum, int step);
    void getUpAndDownSurfacePoints(M3DQuadFigure* quadFig, vtkSmartPointer< vtkPoints > pos, int interpolationLevel);
    void getCrestSurfacePoints(M3DQuadFigure* quadFig, vtkSmartPointer< vtkPoints > pos, int interpolationLevel);

    void drawPointsSetSurface_zcolor(vtkSmartPointer<vtkPoints> points, vtkSmartPointer<vtkRenderer> renderer);
    void drawPointsSetSurface(const char* srepfilename, int interpolationLevel, vtkSmartPointer<vtkRenderer> renderer);
    void drawPointsSetSurface2(vtkSmartPointer<vtkPoints> points, vtkSmartPointer<vtkRenderer> renderer);

    void drawCorrespondenceSpokes(const char* srepfilename, int interpolationLevel, vtkSmartPointer<vtkRenderer> renderer);
    void addSpokesToRender(vtkSmartPointer< vtkPoints > hubpos, vtkSmartPointer<vtkCellArray> cellarraypointsline,
                                                 vtkSmartPointer<vtkRenderer> renderer, int time);
    void drawQuads(vtkSmartPointer< vtkPoints > hubpos, int quadNum, int subQuadPointsNum, int step, vtkSmartPointer<vtkRenderer> renderer);

    void addSpokesToRender_2(vtkSmartPointer< vtkPoints > hubpos, vtkSmartPointer<vtkCellArray> cellarraypointsline,
                                                 vtkSmartPointer<vtkRenderer> renderer, double * quadColor);
    void drawSpokes(const char* srepfilename, int interpolationLevel, vtkSmartPointer<vtkRenderer> renderer);

    void getPoints(M3DQuadFigure* quadFig, int interpolationLevel, int side, vector<double> vars, int varNums, vtkSmartPointer< vtkPoints > surfacePoints);
    vtkSmartPointer< vtkPoints > getShiftedSurfacePoints(const char* srepfilename, int interpolationLevel, vector<double> vars_up, vector<double> vars_down,
                                                  vector<double> vars_crest, int varStand, int varCrest);

    void getUpOrDownSurfacePoints(M3DQuadFigure* quadFig, vtkSmartPointer< vtkPoints > pos, int interpolationLevel, int side);

    void getUVCoordinate(vector<double> &interpolatedU, vector<double> &interpolatedV, int rowNum, int colNum, int step);




private:
    toolsfunc tls;

    double rX;
    double rY;
    double rZ;
    double opacity;

};

#endif // VISUALIZESREPSURFACE_H
