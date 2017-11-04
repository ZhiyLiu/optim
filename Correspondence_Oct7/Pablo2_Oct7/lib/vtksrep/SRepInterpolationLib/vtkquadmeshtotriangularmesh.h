#ifndef VTKQUADMESHTOTRIANGULARMESH_H
#define VTKQUADMESHTOTRIANGULARMESH_H

#include "vtkPolyDataAlgorithm.h"
#include <vector>
using namespace std;

class vtkQuadMeshToTriangularMesh : public vtkPolyDataAlgorithm
{    
public:
    static vtkQuadMeshToTriangularMesh* New();

    // Superclass method to update the pipeline
    virtual int RequestData(vtkInformation* request,
                            vtkInformationVector** inputVector,
                            vtkInformationVector* outputVector);
protected:
    vtkQuadMeshToTriangularMesh();
    ~vtkQuadMeshToTriangularMesh();
};

#endif // VTKQUADMESHTOTRIANGULARMESH_H
