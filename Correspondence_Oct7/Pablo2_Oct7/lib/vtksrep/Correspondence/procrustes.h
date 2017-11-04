#ifndef PROCRUSTES_H
#define PROCRUSTES_H

#include <vector>
#include "alignsrep.h"

class procrustes
{
public:
    procrustes();


    //procrustes(vector<MeshType::Pointer> mesh_b_list);


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

    MeshType::Pointer OPA(MeshType::Pointer MX, MeshType::Pointer MY);
    void GPA(vector<MeshType::Pointer> mesh_x_list);

    MatrixType getCenterOfMass(MatrixType srepPos);
    MatrixType centerThisSrep(MatrixType X, MatrixType COM);
    MeshType::Pointer calculateMean(vector<MeshType::Pointer> mesh_Xp);
    double getG(vector<MeshType::Pointer> mesh_x_list);
    double getFroDistance(MeshType::Pointer A, MeshType::Pointer B);
    vector<MeshType::Pointer> centerSrepsToOrigin(vector<MeshType::Pointer> mesh_x_list);

    vector<MeshType::Pointer> getTrans_X(); // return the transformed tanslation back (plus 0.5,0.5,0.5) sreps list.
    MatrixType getTrans_Rot(); // return the rotaion.
    double getTrans_Sca(); // return the scaling.
    MatrixType getTrans_Tra(); // return the translation.
    void saveMeshVector(const char* filename, vector<MeshType::Pointer> X);
    MatrixType getProcrustesMean();


private:
    MatrixType trans_T; // Translation. 3-by-1
    double trans_scale; // Scalar
    MatrixType trans_R; // Rotation matrix? or rotation transpose?

    //MeshType::Pointer mesh_Xn; // Current srep's boundary mesh.
    vector<MeshType::Pointer> mesh_Xp; // The centered sample sreps.
    //vector<MeshType::Pointer> mesh_Xpmean_input;  // The N-1 centered sreps used to compute the Xmean srep.
    //MeshType::Pointer mesh_Xpmean; // The mean srep.
    //MeshType::Pointer mesh_Xpn;  // The transformed Xn.
    //MeshType::Pointer mesh_b;//boundary mesh.
    //MeshType::Pointer mesh_s;//skeletal mesh.

    unsigned int srepNums; // sample of sreps.
    int pointsNum; // point number on each srep.


};

#endif // PROCRUSTES_H
