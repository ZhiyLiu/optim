#ifndef THINPLATESPLINEPDMTOSREP_H
#define THINPLATESPLINEPDMTOSREP_H

#include "itkThinPlateSplineKernelTransform.h"
#include "itkPointSet.h"
#include <vtkPolyDataReader.h>
#include "vtkPolyData.h"
#include "toolsfunc.h"

#include <vtkPolyDataWriter.h>


class thinplatesplinepdmtosrep
{
public:
    thinplatesplinepdmtosrep();

    typedef double CoordinateRepType;
    typedef itk::ThinPlateSplineKernelTransform< CoordinateRepType,3> TransformType;
    typedef itk::Point< CoordinateRepType, 3 > PointType;
    typedef std::vector< PointType > PointArrayType;
    typedef TransformType::PointSetType PointSetType;
    typedef PointSetType::Pointer PointSetPointer;
    typedef PointSetType::PointIdentifier PointIdType;


    //typedef PointType::PointsContainer PointsContainer;

    typedef vnl_matrix<double> MatrixType;


    TransformType::Pointer get_tps(const char * sourcefilename, const char * targetfilename);
    PointType Vector3DtoPointType(Vector3D point);
    double calculateSpokeLength(PointType tail, PointType tip);
    Vector3D calculateSpokeDirection(PointType tail, PointType tip);
    int tps_to_srep(const char * sourcefilename, const char * targetfilename, const char* sourcesrepfilename, const char* targetsrepFolder);
    void saveModel(M3DQuadFigure* quadFig, const char * targetfilename, const char* sourcesrepfilename, const char* targetsrepFolder);

    int multiple_template_tps(const char * tamplatePDMFileList, const char * targetPDM, const char* tamplateSrepFileList, const char* targetSrepFolder);
    vector< std::string > readTempalte( const char* listName);

    PointType sumPoint(PointType point1, PointType point2);
    M3DQuadFigure* computeMeanSrep(vector<TransformType::Pointer > tpsList,  bool useTPS, vector<M3DQuadFigure* > quadFigList);
    PointType meanPoint(vector<TransformType::Pointer > tpsList, bool useTPS, vector<M3DQuadFigure* > quadFigList, int u, int v, int pointSide);

    void avarage_srep(vector< std::string > srepPathNames, const char* outFileName);

    void avarage_srep_new(vector< std::string > srepPathNames, const char* outFileName);
    PointType mean_point(vector<M3DQuadFigure* > quadFigList, int u, int v);
    Vector3D mean_direction(vector<M3DQuadFigure* > quadFigList, int u, int v, int side);
    double mean_radius(vector<M3DQuadFigure* > quadFigList, int u, int v, int side);
    M3DQuadFigure* compute_mean_srep_new(vector<M3DQuadFigure* > quadFigList);
    void cpns_avarage_srep(vector< std::string > srepPathNames, const char* outFileName);



private:
    toolsfunc tls;

};

#endif // THINPLATESPLINEPDMTOSREP_H
