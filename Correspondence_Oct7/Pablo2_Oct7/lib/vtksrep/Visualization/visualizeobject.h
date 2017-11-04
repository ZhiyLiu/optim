#ifndef VISUALIZEOBJECT_H
#define VISUALIZEOBJECT_H


#include "visualizecrest.h"
#include "visualization.h"

#include "slidestandardspokes.h"
#include "slidecrestspokes.h"

#include <vtkSphereSource.h>



class visualizeobject
{
public:
    visualizeobject();
    visualizeobject(M3DQuadFigure* quadfig, int interpolationLevel, double *quadColor, vtkSmartPointer<vtkRenderer> renderer,
                    bool moved, double rX, double rY, double rZ, double opacity);

    // Draw the frame(non-filled sub quads) of the boundary
    void showBoundaryFrames(M3DQuadFigure* quadFig, bool up, bool down, bool crest);

    // Draw the filled sub quads of the boundary
    void showBoundaryQuads(M3DQuadFigure* quadFig, bool up, bool down, bool crest);

    // Draw the skeletal sheet and spokes.
    void showStandardSrep(M3DQuadFigure* quadFig, int side, double *quadColor, int interpolationLevel, int linewith,bool showspoketail);
    //void showCrestSrep(M3DQuadFigure* quadFig, int interpolationLevel);
    void showSrep(M3DQuadFigure* quadFig, int interpolationLevel, int linewith);
    void showCrestSrep(M3DQuadFigure* quadFig, double * quadColor, int interpolationLevel, int linewidth, bool showspoketail);

    // Draw the old and new srep together.
    void showShiftedSrep(vector<double> coeff, bool onlyShifted, int varStand, int varCrest, int interpolationLevel, int linewith, bool showspoketail);
    void showShiftedCrestFrame(vector<double> crest_vars, int varCrest);
    void showShiftedCrestFrameAndSpokes(vector<double> crest_vars, int varCrest, int interpolationLevel);

    void showBoundaryTriangulars(M3DQuadFigure* quadFig, bool up, bool down, bool crest);
    void showSpokes(M3DQuadFigure* quadFig, bool up, bool down, bool crest);
    void showAtomPoints(M3DQuadFigure* quadFig);

    void showCrestFrameByQuadIndex(M3DQuadFigure* quadFig, int quadIndex, int interpolationLevel, double * quadColor);
    void showCrestQuadByQuadIndex(M3DQuadFigure* quadFig, int quadIndex, int interpolationLevel, double * quadColor);
    void showCrestSpokesByQuadIndex(M3DQuadFigure* quadFig, int quadIndex, int interpolationLevel, double * quadColor);
    void showCrestVolumesByQuadIndex(M3DQuadFigure* quadFig, int quadIndex, int interpolationLevel, double * quadColor);

    void drawSrep(M3DQuadFigure* quadFig, bool showup, bool showdown, bool showcrest,double *quadColor, int linewith);
    void showShiftedUpOrDownSideSrep(vector<double> coeff, double * color, int side, int linewith, bool showspoketail);

    void drawSrep_2(M3DQuadFigure* quadFig);

    void showShiftedCrestSideSrep(vector<double> coeff, double * color, int linewidth, bool showspoketail);
    void showShiftedUpOrDownSideSkeletalSheet(vector<double> coeff, double * color, int side, double opacity);

    void showShiftedCrestCurve(vector<double> coeff, double * color, double opacity);

    void showOriginalCrestSpokesByQuadIndex(M3DQuadFigure* quadFig, int quadIndex, int interpolationLevel, double * quadColor, bool showspoketail);

    void showSrep_2(M3DQuadFigure* quadFig, int interpolationLevel, int linewith);

    void showAtomPoints_2(M3DQuadFigure* quadFig, double radius, double*color);
    void showSpokes_2(M3DQuadFigure* quadFig, bool up, bool down, bool crest, double* color, double opacity);


private:
    vtkSmartPointer<vtkRenderer> renderer;
    double *quadColor;

    int interpolationLevel;
    int rowNums;
    int colNums;
    int quadNums;
    //int subQuadPointsNum;
    int step;
    bool moved;

    M3DQuadFigure* quadFig;
    toolsfunc tools;

    double rX;
    double rY;
    double rZ;
    double opacity;


};

#endif // VISUALIZEOBJECT_H
