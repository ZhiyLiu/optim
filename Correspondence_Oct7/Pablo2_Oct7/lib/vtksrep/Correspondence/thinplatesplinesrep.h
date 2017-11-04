#ifndef THINPLATESPLINESREP_H
#define THINPLATESPLINESREP_H


#include "itkThinPlateSplineKernelTransform.h"
#include "itkPointSet.h"
#include <fstream>
#include "alignsrep.h"
#include "itkMatrix.h"

#include <vtkPolyDataReader.h>
#include "vtkPolyData.h"


class thinplatesplinesrep
{
public:
    thinplatesplinesrep();   

    typedef double CoordinateRepType;
    typedef itk::ThinPlateSplineKernelTransform< CoordinateRepType,3> TransformType;
    typedef itk::Point< CoordinateRepType, 3 > PointType;
    typedef std::vector< PointType > PointArrayType;
    typedef TransformType::PointSetType PointSetType;
    typedef PointSetType::Pointer PointSetPointer;
    typedef PointSetType::PointIdentifier PointIdType;


    //typedef PointType::PointsContainer PointsContainer;

    typedef vnl_matrix<double> MatrixType;

    // itk typedefs
    typedef itk::DefaultDynamicMeshTraits<double, 3, 3, double, double> MeshTraitsType;
    typedef itk::Mesh<double, 3, MeshTraitsType>                      MeshType;

    // Access points
    typedef MeshType::PointsContainer::Iterator     PointsIterator;

    // Declare the type for PointsContainer
    typedef MeshType::PointsContainer     PointsContainerType;

    // Declare the type for PointsContainerPointer
    typedef MeshType::PointsContainerPointer PointsContainerPointer;

    // Declare the type for the filter
    typedef itk::TransformMeshFilter<MeshType, MeshType, TransformType>       FilterType;

    MeshType::Pointer  srepPointsSet(const char* sourcesrepfilename);
    int tps_srep(const char * sourcefilename, const char * targetfilename, const char* sourcesrepfilename);
    void generateSrep(PointSetType::Pointer outputMesh, const char* srepfilename);
    void saveGeneratedSreps(const char* filePathName, MatrixType newSpokeDirection, vector<double> newSpokeRadius, MatrixType transformedTails);
    void setCrestSpoke(int rowIndex, int colIndex, MatrixType newSpokeDirection,
                    vector<double> newSpokeRadius, int pIndex, int crestIndex, M3DQuadFigure* quadFig, MatrixType transformedTails);
    MatrixType correspondenceSpokeTailAndTip(MatrixType sPoints, int rowNums, int colNums,int totalAtomNums, int spokeNum);
    MatrixType getNewSpokesDirection(MatrixType bPoints, MatrixType sPoints);
    vector<double> getNewSpokesRadius(MatrixType bPoints, MatrixType sPoints);

};

#endif // THINPLATESPLINESREP_H
