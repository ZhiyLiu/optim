#ifndef SREPPROCRUSTES_H
#define SREPPROCRUSTES_H

#include "alignsrep.h"
//#include "procrustes.h"


class srepprocrustes
{
public:
    srepprocrustes();
    srepprocrustes(const char* srepFolder);

   // itk typedefs
    typedef itk::DefaultDynamicMeshTraits<double, 3, 3, double, double> MeshTraitsType;
    typedef itk::Mesh<double, 3, MeshTraitsType>                      MeshType;

    // Access points
    typedef MeshType::PointsContainer::Iterator     PointsIterator;

    // Define Procrustes alignment filter
    typedef itk::Mesh3DProcrustesAlignFilter<MeshType, MeshType> ProcrustesFilterType;

    // Declare the transform type
    //typedef double                                                        OutputMeshType;
    //typedef typename OutputMeshType::CoordRepType                         CoordRepType;
    //typedef AffineTransform<CoordRepType, 3>                              TransformType;
    typedef itk::AffineTransform<double,3>                              TransformType;
    // Declare the type for the filter
    typedef itk::TransformMeshFilter<MeshType, MeshType, TransformType>       FilterType;

    // Declare the type for PointsContainer
    typedef MeshType::PointsContainer     PointsContainerType;

    // Declare the type for PointsContainerPointer
    typedef MeshType::PointsContainerPointer PointsContainerPointer;

    // Declare the type for Points
    typedef MeshType::PointType           PointType;



    int procrustesAlignment( );
    void writeAlignToSkeletal(string filePathName, MeshType::Pointer outputMesh);


private:
    MeshType::Pointer mesh_b;//boundary mesh.
    MeshType::Pointer mesh_s;//skeletal mesh.
    vector<MeshType::Pointer> mesh_s_list;//store all the input srep's skeletal mesh.
    string srepFolder;//directory where a list of s-reps which we will performance alignment on.


    toolsfunc tls;

    TransformType::Pointer affineTransform;//a pointer, point to several matrixs, each matrix has 3 columns.

    //MatrixType_3 pMatrix;
    //int pRowIndex;
};

#endif // SREPPROCRUSTES_H
