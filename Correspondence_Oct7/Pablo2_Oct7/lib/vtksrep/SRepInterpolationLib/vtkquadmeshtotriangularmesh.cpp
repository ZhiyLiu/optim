#include "vtkquadmeshtotriangularmesh.h"

#include "vtkCellArray.h"
#include "vtkTriangle.h"
#include "vtkSmartPointer.h"

#include "vtkInformationVector.h"
#include "vtkInformation.h"

#include "vtkObjectFactory.h"
vtkStandardNewMacro(vtkQuadMeshToTriangularMesh);

vtkQuadMeshToTriangularMesh::vtkQuadMeshToTriangularMesh()
{
}

vtkQuadMeshToTriangularMesh::~vtkQuadMeshToTriangularMesh()
{
}

// Superclass method to update the pipeline
int vtkQuadMeshToTriangularMesh::RequestData(vtkInformation* request,
                        vtkInformationVector** inputVector,
                        vtkInformationVector* outputVector){



    // get the info objects
    vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
    vtkInformation *outInfo = outputVector->GetInformationObject(0);

    // get the input and output
    vtkPolyData *input = dynamic_cast<vtkPolyData*>(vtkPolyData::SafeDownCast(inInfo->Get(vtkPolyData::DATA_OBJECT())));
    vtkPolyData *output = vtkPolyData::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));

    vtkCellArray* outputcellarray = vtkCellArray::New();

    for(unsigned i = 0; i < input->GetNumberOfCells(); i++){
        vtkCell* quad = input->GetCell(i);

        vtkIdType id0 = quad->GetPointIds()->GetId(0);
        vtkIdType id1 = quad->GetPointIds()->GetId(1);
        vtkIdType id2 = quad->GetPointIds()->GetId(2);
        vtkIdType id3 = quad->GetPointIds()->GetId(3);

        vtkSmartPointer<vtkTriangle> cell0 = vtkSmartPointer<vtkTriangle>::New();
        cell0->GetPointIds()->SetId(0, id0);
        cell0->GetPointIds()->SetId(1, id1);
        cell0->GetPointIds()->SetId(2, id3);

        vtkSmartPointer<vtkTriangle> cell1 = vtkSmartPointer<vtkTriangle>::New();
        cell1->GetPointIds()->SetId(0, id1);
        cell1->GetPointIds()->SetId(1, id2);
        cell1->GetPointIds()->SetId(2, id3);

        outputcellarray->InsertNextCell(cell0);
        outputcellarray->InsertNextCell(cell1);
    }

    output->SetPoints(input->GetPoints());
    output->SetPolys(outputcellarray);

    return 1;
}
