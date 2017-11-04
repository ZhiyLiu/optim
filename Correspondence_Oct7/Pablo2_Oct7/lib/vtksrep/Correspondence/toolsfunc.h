#ifndef TOOLSFUNC_H
#define TOOLSFUNC_H

#include <vector>
#include "M3DQuadFigure.h"
#include "P3DControl.h"
#include "ControlParms.h"
#include "M3DQuadInterpolater.h"   //Jared's interpolation method to up and down spoke.
#include <iostream>
#include <fstream>
#include "time.h"

#include <vtkPolyDataWriter.h>
#include <vtkPolyDataReader.h>
#include "vtkPolyData.h"


typedef vector< double > VectorSRepFeaturesType;//store the features of each srep
typedef vector< VectorSRepFeaturesType > VectorTrainingSetFeaturesType;//sotre all the srep features over the training set. each row is a feature, each column is a case.

typedef vector<double> VectorDoublePoints;
typedef vector<VectorDoublePoints> VectorQuadPoint;//store the xyz position of the points.

typedef vector<Vector3D> Vector3DPosition;
typedef vector<Vector3DPosition> QuadPoint;


// itk typedefs
/*typedef itk::DefaultDynamicMeshTraits<float, 3, 3, float, float> MeshTraitsType;
typedef itk::Mesh<float, 3, MeshTraitsType>                      MeshType;

typedef vector<Vector3DPosition> Vector3DPositionList;
//Alignment to s-rep using Procrustes.
typedef itk::Mesh3DProcrustesAlignFilter<MeshType, MeshType> ProcrustesFilterType;
*/



class toolsfunc
{
public:
    toolsfunc();


    M3DQuadFigure* GetQuadFigure(const char* figfilename);

    string connectFilePath(string rootDir, const char* specificDir, int changingIndex, const char* extension);
    string conStr(string rootDir, const char* specificDir, string fileName, const char* extension) const;
    string connectStr(string rootDir, const char* specificDir) const;


    //clear all the files in a folder.
    void deleteFiles(string path);

    vector< std::string > getModelname(string filefolder);

    void savefeatures(string path, VectorTrainingSetFeaturesType featurematrix);

    void splitLine(double point1, double point2, int step, vector<double> &resultpoints);
    vector<double> splitQuad(double p0, double p1, double p2, double p3, int step);

    string getFileNameWithExtension(string path);
    string splitExtension(string fileName);
    string getPathName(string path);

    void saveMatrix3(string path, vector<double> A, vector<double> B,vector<double> C);
    void saveMatrix2(string path, vector<double> A, vector<double> B);

    double quadArea(Vector3D *point);

    string getBaseFileName(string path);
    bool quadAreaMin(Vector3D *point);

    double RunningTime(clock_t time1, clock_t time2);
    void save_vector_to_txt(const char* filename, vector<double> data);

    vector<double> splitQuad_2(double p0, double p1, double p2, double p3, int step);

    void getUVCoordinate(vector<double> &interpolatedU, vector<double> &interpolatedV, int rowNum, int colNum, int step);
    int getSSSN(int interpolationLevel, int rowNum, int colNum);
    int getCSSN(int interpolationLevel, int rowNum, int colNum);

    void saveVtkPointsToVTKFile(vtkSmartPointer<vtkPoints> points, string outputFile);

    double dotProductAngle(double * point_c, double * point_h, double * point_v);
    double dotProductAngle_2(Vector3D point_c, Vector3D point_h, Vector3D point_v);
    //Vector3D unitVector(Vector3D unitVec);
    Vector3D getVector(double * point1, double * point2);

    Vector3D trangularNormal_vec3D(Vector3D point_c, Vector3D point_h, Vector3D point_v);
    Vector3D trangularNormal_doublePot(double * point_c, double * point_h, double * point_v);
    Vector3D crossProduct(Vector3D a, Vector3D b);

    double lengthofedges(Vector3D point[2]);

};

#endif // TOOLSFUNC_H


