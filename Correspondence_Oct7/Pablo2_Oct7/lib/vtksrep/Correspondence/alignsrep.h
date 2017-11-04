/* Procrustes Alignment to sreps.
 * Liyun Tu
 * Feb 28, 2014
*/


/* First do alignment to the implied boundary of each s-rep, get the transformation and then
 * apply this transformation to its correspondence skeletal sheet.
*/

#ifndef ALIGNSREP_H
#define ALIGNSREP_H

#include <itkDefaultDynamicMeshTraits.h>
#include <itkMesh.h>

#include "M3DFigure.h"
#include "toolsfunc.h"
#include "itkMesh3DProcrustesAlignFilter.h"
#include <itkTransformMeshFilter.h>
#include <itkAffineTransform.h>

#include "M3DQuadFigure.h"
#include "P3DControl.h"
#include "ControlParms.h"
#include <vtksrepinterpolatecrestspokesquartic.h>
#include <vtksrepinterpolatemedialcrestcurve.h>

#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderWindow.h>

#include "vtkreadsrep.h"
#include <vtksrep.h>

//const int atomNums = 39;    //the total atoms on each srep.
//const int srepNums = 5;

class alignsrep
{
public:
    alignsrep();
    //alignsrep(const char* srepFolder);
    alignsrep(const char* filename, int interpolationLevel);




    //int prepareDataForMatlabProcr();
    //void saveGPA_Input_Matrix(const char* filename);

    void drawBoundary();
    void drawSkeletalCrestCurve(int interpolationlevel);    
    void drawCrestSpokes(int interpolationlevel);
    void getPointsOnBoundary();
    void getAllSubQuadsVertexesOnBoundary();
    vector<double> calculateAreas_method2(vector<Vector3D> surfacePoints);
    double quadArea(Vector3D *point);
    void getCrestQuadsArea();
    void getCrestQuadsArea_InteLevel0();
    void weightAreaToBoundaryPoints();
    void addBoundaryPointsAreas(double area1, double area2, double area3, double area4);
    void correPointAndArea(vector<double> quads_areas, vector<Vector3D> mid_points, vector<double> crest_areas,
                                               vector<Vector3D> crest_points);
    vector<Vector3D> weightPoints();
    void saveMatrix(vector<Vector3D> points,const char* filename);

    vector<Vector3D> getBoundaryPoints();
    vector<double> getBoundaryPointsAreas();

    vector<Vector3D> getPointsOnSkeletal();
    vector<double> getCrestAreas(int side);
    void drawCrestSpokes_method2();

    vector<double> getQuadsAreas(int side);    

    vector<double> getOnesWeight(int pointsNum);


    int getRows(); // return the srep's row number
    int getCols();
    int getAtomNums();
    int getInteriorAtomNums();
    int getCrestAtomNums();



private:

    toolsfunc tls;
    M3DQuadFigure* quadFig;
    int crestAtomNums;
    int crestQuadNum;
    int colNums;
    int rowNums;
    int quadNum;
    int subQuadPointsNum;
    int interiorAtomNums;
    int totalAtomNums;
    int interpolationLevel;//used only when calculating the quad area. Diff leve can got diff area precise.
    int step;
    //string figfilename;



    vtkSmartPointer<vtkRenderer> renderer;
    vtkSmartPointer<vtkSRep> srepfig;
    vtkSmartPointer<vtkReadSRep> readsrep;

    //vector<Vector3D> up_b; //store all the points on the top boundary.
   // vector<Vector3D> down_b; //store all the points on bottom boundary.
    vector<Vector3D> up_b_mid; //store the interior points on the top boundary.
    vector<Vector3D> down_b_mid; //store the interior points on bottom boundary.
    vector<Vector3D> up_b_crest; //store the crestAtomNums points on the top boundary.
    vector<Vector3D> down_b_crest; //store the crestAtomNums points on the bottom boundary.
    vector<Vector3D> medial_crest; //store the crestAtomNums points on the medial boundary.


    // For all the subquads' vertexes.
    vector<Vector3D> up_BP;//up boundary points list
    vector<Vector3D> down_BP;//down boundary points list


    vector<double> up_b_crest_areas, down_b_crest_areas;

    vector<Vector3D> boundaryPoints;
    VectorTrainingSetFeaturesType boundaryPointsAreas;
    int pointIndex;

    vector<double> up_quads_areas;
    vector<double> down_quads_areas;

};

#endif // ALIGNSREP_H
