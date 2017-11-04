#ifndef WEIGHTEDPROCRUSTES_H
#define WEIGHTEDPROCRUSTES_H

#include <vector>
#include "alignsrep.h"
#include "vnl/algo/vnl_symmetric_eigensystem.h"
#include "itkMatrix.h"



typedef vector<double> VectorDoublePoints;
typedef vector<VectorDoublePoints> VectorPoints;//store the xyz position of the points.


// itk typedefs
typedef itk::DefaultDynamicMeshTraits<double, 3, 3, double, double> MeshTraitsType;
typedef itk::Mesh<double, 3, MeshTraitsType>                      MeshType;

// Access points
typedef MeshType::PointsContainer::Iterator     PointsIterator;

// Declare the type for PointsContainer
typedef MeshType::PointsContainer     PointsContainerType;

// Declare the type for PointsContainerPointer
typedef MeshType::PointsContainerPointer PointsContainerPointer;

// Declare the type for Points
typedef MeshType::PointType           PointType;


typedef vnl_matrix<double> MatrixType;



class weightedprocrustes
{
public:
    weightedprocrustes();
    weightedprocrustes(bool realWeight);

    vector<MatrixType> weightedGPA(vector<MeshType::Pointer> X, VectorPoints W);

    // Apply the rotationMatrix, scaleVector and centroid to correspondence points set.
    vector<MatrixType> applyTransform(vector<MeshType::Pointer> X);

    double getFroDistance(MatrixType A, MatrixType B);
    double procrustesDistance(vector<MatrixType> srepList);

    // Return the new spokes direction between correspondence boundary and skeletal point.
    vector<MatrixType> getNewSpokesDirection(vector<MatrixType> bp, vector<MatrixType> sp);

    // Return the new spokes radius between correspondence boundary and skeletal point.
    VectorPoints getNewSpokesRadius(vector<MatrixType> bp, vector<MatrixType> sp);

    vector<MatrixType> getXp();
    vector<MatrixType> getQ();


private:
    vector<MatrixType> rotationMatrix; //rotationMatrix[i] holding the ith srep's rotation matrix.
    vector<double> scaleVector; // scaleVector[i] holding the ith srep's scale factor.
    MeshType::Pointer centroid; // Xbar[i] is centroid of ith srep.

    vector<MeshType::Pointer> mesh_Xp;
    vector<MeshType::Pointer> mesh_Q;

    bool realWeight; //if true, use correspondence area weight, if false, use all weight to 1;




};

#endif // WEIGHTEDPROCRUSTES_H
