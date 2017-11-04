#include "vtksreptogridlayout.h"


#include "vtkSmartPointer.h"
#include "vtksrep.h"

#include "vtkCellArray.h"

#include "vtkExecutive.h"
#include "vtkPointData.h"


#include "vtkObjectFactory.h"
vtkStandardNewMacro(vtkSRepToGridLayout);



vtkSRepToGridLayout::vtkSRepToGridLayout()
{

    m_Output = 0;

    // by default assume filters have one input and one output
    // subclasses that deviate should modify this setting
    this->SetNumberOfInputPorts(1);
    this->SetNumberOfOutputPorts(1);

}

vtkSRepToGridLayout::~vtkSRepToGridLayout()
{
}

void vtkSRepToGridLayout::SetInput(int index, vtkDataObject* input){
    if(input)
      {
      this->SetInputConnection(index, input->GetProducerPort());
      }
    else
      {
      // Setting a NULL input removes the connection.
      this->SetInputConnection(index, 0);
      }
}

void vtkSRepToGridLayout::SetInput(vtkDataObject* input){
    this->SetInput(0, input);
}



// Superclass method to update the pipeline
void vtkSRepToGridLayout::Update(){



    // get the input and output
    vtkSRep *input = dynamic_cast<vtkSRep*>(this->GetExecutive()->GetInputData(0, 0));

    m_Output = vtkSmartPointer<vtkSRep>::New();
    vtkSmartPointer< vtkPoints > hubpos = vtkSmartPointer< vtkPoints >::New();

    int numcols = input->GetNumColumns();

    for(int i = 0; i < input->GetNumberOfPoints(); i++){

        int nr = i/numcols - (numcols - 1)/2.0;
        int nc = (i%numcols)- (numcols - 1)/2.0;

        hubpos->InsertNextPoint(nr, nc, 0);
    }


    m_Output->SetPoints(hubpos);
    m_Output->GetPointData()->SetScalars(input->GetPointData()->GetScalars());
    m_Output->SetPolys(input->GetPolys());
    m_Output->SetAllSpokes(input->GetAllSpokes());
    m_Output->SetAllRadius(input->GetAllRadius());

}
