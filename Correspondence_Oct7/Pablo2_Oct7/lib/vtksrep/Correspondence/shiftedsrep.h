#ifndef SHIFTEDSREP_H
#define SHIFTEDSREP_H

#include "toolsfunc.h"
#include "slidestandardspokes.h"
#include "slidecrestspokes.h"
#include "visualizecrest.h"
#include "visualizesrepsurface.h"

#include "vnl/algo/vnl_symmetric_eigensystem.h"
#include "itkMatrix.h"

#include <vtkPolyDataWriter.h>
#include <vtkPolyDataReader.h>
#include "vtkPolyData.h"



typedef vnl_matrix<double> MatrixType;


const int rowNum = 3;
const int colNum = 13;


class shiftedsrep
{
public:
    shiftedsrep();
    shiftedsrep(const char* rootDir, vector< std::string > srepNames, double rX, double rY, double rZ, double opacity);
    shiftedsrep(double rX, double rY, double rZ, double opacity);

    void readTrainingSreps();
    void getPoints(M3DQuadFigure* quadFig, int side, vector<double> vars, int varNums, vector<MatrixType> &boundaryPos,
                                vector<MatrixType> &skeletalPos);
    void meanSide(MatrixType mean_bpos, MatrixType mean_spos, int side);
    MatrixType meanSidePointPair(vector<MatrixType> points);
    void computeMeanSide(string vars_str, int side, int varNums);
    void computeMeanFig(string up_str, string down_str, string crest_str);

    void drawSpokes(bool showup, bool showdown, bool showcrest, int interpolationLevel, vtkSmartPointer<vtkRenderer> renderer);
    void drawMeanSrep();
    vtkSmartPointer< vtkPoints > getMeanSurfacePoints(int interpolationLevel);
    void drawPointsSetSurfaceOfTheMean(int interpolationLevel, vtkSmartPointer<vtkRenderer> renderer);

    void deformAlongEigMode(MatrixType X, int pointNum, int sampleNum);
    //int getPointsNum(int interpolationLevel);
    void saveMatrix(const char* filename, MatrixType matrix);
    void drawSurfacePoint(const char* filename, vtkSmartPointer<vtkRenderer> renderer);
    void saveToVTK(int pointNum, MatrixType mShape, string outputFile, int col);

    void deformOriginalSurfaceMean( int interpolationLevel);
    void surfacePointsEachSide(string vars_str, int side, int varNums, MatrixType &X, int interpolationLevel);

    int calPointsNum(int interpolationLevel);
    void getSurfacePoints(M3DQuadFigure* quadFig, int side, vector<double> vars, int varNums, vtkSmartPointer< vtkPoints > surfacePos,
                                int interpolationLevel);
    int getSSPN(int interpolationLevel);
    void deformShiftedSurfaceMean(string up_str, string down_str, string crest_str, int interpolationLevel);


private:
    vector<M3DQuadFigure *> quadFigList;//holding the input sreps.

    const char* rootDir;
    vector< std::string > srepNames;

    toolsfunc tls;

    int crestAtomNum; // 28. crest spoke number, each spoke has one variable.
    int interiorAtomNum; // 11
    int varStand; // 46. variables for standard spokes (up or down spokes).

    M3DQuadFigure* upMeanQuadFig; // storing the mean up side
    M3DQuadFigure* downMeanQuadFig; // storing the mean down side
    M3DQuadFigure* crestMeanQuadFig; // storing the mean crest side

    double rX;
    double rY;
    double rZ;
    double opacity;

    int srepNum;


};

#endif // SHIFTEDSREP_H
