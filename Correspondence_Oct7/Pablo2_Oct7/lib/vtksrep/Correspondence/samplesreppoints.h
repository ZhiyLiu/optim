#ifndef SAMPLESREPPOINTS_H
#define SAMPLESREPPOINTS_H


#include "toolsfunc.h"

#include <vtkPolyDataWriter.h>
#include <vtkPolyDataReader.h>
#include "vtkPolyData.h"

#include "slidestandardspokes.h"
#include "slidecrestspokes.h"
#include "visualizecrest.h"


typedef vnl_matrix<double> MatrixType;
//int row = 3;
//int col = 13;

class samplesreppoints
{
public:
    samplesreppoints();

    //vector<MatrixType> sampleSrepPoints(const char* srepFolder, const char* outputFolder, bool repeatSkeletalPoints);
    MatrixType samplePointsWithRepeatedSkeletal(M3DQuadFigure* quadFig);

    void saveToVTK(MatrixType X, string outputFile);

    void sampleSrepPoints(const char* srepFolder, vector< std::string > filename, const char* outputFolder, const char* coeffFolder, int rowNum,
                          int colNum, int type, int interpolationLevel);
    void sampleShiftedSrep(const char* srepFolder, vector< std::string > filename, int side, string vars_str_file, int varNums,
                                             vector<vtkSmartPointer<vtkPoints> > &sampledPoints,int type, int interpolationLevel);
    void getSurfacePoints(M3DQuadFigure* quadFig, int side, vector<double> vars, int varNums, vtkSmartPointer< vtkPoints > surfacePos, int type,
                          int interpolationLevel);
    void getSrepPoints(M3DQuadFigure* quadFig, int side, vtkSmartPointer< vtkPoints > surfacePos, int type, int interpolationLevel);

    void saveToText(vector<MatrixType> data, const char* filename);

    void getUpOrDownSidePoints(M3DQuadFigure* quadFig, vtkSmartPointer< vtkPoints > pos, int interpolationLevel, int side, int type);
    void getCrestSidePoints(M3DQuadFigure* quadFig, vtkSmartPointer< vtkPoints > pos, int interpolationLevel, int type);

    void saveVtkPointsToVTK(vtkSmartPointer<vtkPoints> points, string outputFile);
    void saveVtkPointsVectorToText(vector<vtkSmartPointer<vtkPoints> > data, const char* filename);

private:
    toolsfunc tls;


};

#endif // SAMPLESREPPOINTS_H
