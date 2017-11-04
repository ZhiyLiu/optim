#ifndef VISUALIZECREST_H
#define VISUALIZECREST_H


#include "toolsfunc.h"
#include <vtksrepinterpolatecrestspokesquartic.h>
#include "visualization.h"

#include "vtkreadsrep.h"
#include <vtksrep.h>

#include <vtkVersion.h>
#include <vtkSmartPointer.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkPointData.h>
#include <vtkCellArray.h>
#include <vtkUnsignedCharArray.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkVertexGlyphFilter.h>
#include <vtkProperty.h>
#include <vtkLookupTable.h>


class visualizecrest
{
public:
    visualizecrest();
    visualizecrest(M3DQuadFigure* quadfig, int interpolationLevel, double *quadColor, vtkSmartPointer<vtkRenderer> renderer,
                   double rX, double rY, double rZ, double opacity);

    void drawCrestFrames();
    void drawCrestQuads();
    void drawCrestSpokes();
    void drawCrestCurve(double linewidth, double * color);
    void drawMidOneCrestSpokes(double rX, double rY, double rZ, int linewidth, double * color, bool showpointball);
    vtkSmartPointer<vtkSRep> getSrepFig(M3DQuadFigure* quadfig);

    void connectDiagnoal();

    void drawCorrespondenceCrestSpokes();

    void drawCrestFramesByQuadIndex(int quadIndex, double *quadColor);
    void drawCrestQuadByQuadIndex(int quadIndex, double *quadColor);

    void drawCrestSpokesByQuadIndex(int quadIndex, double *quadColor);
    void drawCrestVolumesByQuadIndex(int quadIndex, double *quadColor);

    void drawCrestQuadWithDiffColor(int quadNum);
    void drawCrestBoundaryEquator(double rX, double rY, double rZ);

    void drawMiddleCrestSpokesByQuadIndex(int quadIndex, double *quadColor, bool showspoketail);


private:
    vtkSmartPointer<vtkRenderer> renderer;
    double *quadColor;
    int crestQuadNum;

    M3DQuadFigure* quadFig;

    toolsfunc tools;

    vector<double> subpoints_v; //store the v coordinate of all the sub-point along spoke tip connection curve.
    vector<double> subpoints_t; //store the u coordinate of all the sub-points along crest.
    vtkSmartPointer<vtkSRepInterpolateCrestSpokesQuartic> interpolatecrestspokes;

    double rX;
    double rY;
    double rZ;
    double opacity;

};

#endif // VISUALIZECREST_H
