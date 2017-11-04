#include "vtksrepvisuprimitives.h"


#include "vtkSmartPointer.h"
#include "vtksrep.h"

#include "vtkCellArray.h"
#include "vtkQuad.h"
#include "vtkVertex.h"
#include "vtkLine.h"

#include "vtkExecutive.h"
#include "vtkPolyDataMapper.h"
#include "vtkActor.h"
#include "vtkProperty.h"
#include "vnl/vnl_cross.h"


#include "vtkObjectFactory.h"
vtkStandardNewMacro(vtkSRepVisuPrimitives);



vtkSRepVisuPrimitives::vtkSRepVisuPrimitives()
{

    // by default assume filters have one input and one output
    // subclasses that deviate should modify this setting
    this->SetNumberOfInputPorts(1);
    this->SetNumberOfOutputPorts(1);


    m_PrimitivesCollectionPolyData = vtkSmartPointer< vtkPolyDataCollection >::New();
    m_PrimitivesCollectionCells  = vtkSmartPointer< vtkCollection >::New();
    m_PrimitivesCollectionPoints = vtkSmartPointer< vtkCollection >::New();
    m_PrimitivesCollectionActors = vtkSmartPointer< vtkActorCollection >::New();


}

vtkSRepVisuPrimitives::~vtkSRepVisuPrimitives()
{
}

void vtkSRepVisuPrimitives::SetInput(int index, vtkDataObject* input){
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

void vtkSRepVisuPrimitives::SetInput(vtkDataObject* input){
    this->SetInput(0, input);
}



// Superclass method to update the pipeline
void vtkSRepVisuPrimitives::Update(){



    // get the input and output
    m_SRep = dynamic_cast<vtkSRep*>(this->GetExecutive()->GetInputData(0, 0));
    vtkSRep *input = m_SRep;

    int numprimitives = 8;


    for(int i = 0; i < numprimitives; i++){
        vtkSmartPointer<vtkPolyData> poly = vtkSmartPointer<vtkPolyData>::New();

        vtkSmartPointer< vtkPoints > polypoints = vtkSmartPointer< vtkPoints >::New();
        vtkSmartPointer<vtkCellArray> polycellarray = vtkSmartPointer<vtkCellArray>::New();

        m_PrimitivesCollectionPolyData->AddItem(poly);
        m_PrimitivesCollectionPoints->AddItem(polypoints);
        m_PrimitivesCollectionCells->AddItem(polycellarray);

    }

    for(unsigned i = 0; i < input->GetNumberOfPoints(); i++){

        double *temp = input->GetPoint(i);
        vtkSRep::VNLType p0(3);
        p0[0] = temp[0];
        p0[1] = temp[1];
        p0[2] = temp[2];

        for(int j = 0; j < numprimitives; j++){

            vtkPoints* polypoints = (vtkPoints*) m_PrimitivesCollectionPoints->GetItemAsObject(j);
            vtkIdType idp0 = polypoints->InsertNextPoint(p0[0], p0[1], p0[2]);

            if(j == 0){
                vtkSmartPointer<vtkVertex> vertex = vtkSmartPointer<vtkVertex>::New();
                vertex->GetPointIds()->SetId(0,idp0);

                vtkCellArray* cells = (vtkCellArray*)m_PrimitivesCollectionCells->GetItemAsObject(j);
                cells->InsertNextCell(vertex);
            }else{

                vtkSRep::VNLType pend0(3);
                pend0.fill(0);

                switch(j){
                    case 1://top spoke
                    {                        
                        pend0 = input->GetSpoke(i,vtkSRep::TOP_SPOKE) + p0;
                    }
                        break;
                    case 2://bottom spoke
                    {
                        pend0 = input->GetSpoke(i,vtkSRep::BOTTOM_SPOKE) + p0;
                    }
                        break;
                    case 3:
                    {
                        vtkSRep::VectorVNLType n0spokes = input->GetSpokes(i);                        

                        //crest spoke
                        if(n0spokes.size() > 2){

                            for(unsigned k = 2;  k < n0spokes.size(); k++){
                                pend0 = input->GetSpoke(i, k) + p0;

                                vtkIdType idp0end = polypoints->InsertNextPoint(pend0[0], pend0[1], pend0[2]);
                                vtkSmartPointer<vtkLine> line = vtkSmartPointer<vtkLine>::New();
                                line->GetPointIds()->SetId(0, idp0);
                                line->GetPointIds()->SetId(1, idp0end);
                                vtkCellArray* cells = (vtkCellArray*)m_PrimitivesCollectionCells->GetItemAsObject(j);
                                cells->InsertNextCell(line);
                            }
                        }
                    }
                        break;
                    case 4://normal to the sheet
                    {
                        //vtkSRep::VectorVNLType n0spokes = input->GetSpokes(i);
                        vtkSRep::VectorDoubleType n0radius = input->GetSpokesRadius(i);

                        //normal to the sheet

                        pend0 = input->GetMedialSheetNormal(i)*n0radius[0] + p0;


                        vtkSRep::VNLType norm(3);
                        norm.fill(0);
                        if(input->GetSpokes(i).size() > 0){
                            norm = vnl_cross_3d(input->GetSpoke(i, 1, true), input->GetSpoke(i, 0, true));
                            norm.normalize();//this normal is the normal of the wheel
                        }

                        pend0 = norm*n0radius[0] + p0;
                    }
                        break;
                    case 5://uderivatives
                    {

                        vtkSRep::VectorVNLType dp0 = input->GetDerivatives(i);
                        if(dp0.size() > 0)
                            pend0 = p0 + dp0[0];
                    }
                        break;
                    case 6://vderivatives
                    {
                        vtkSRep::VectorVNLType dp0 = input->GetDerivatives(i);
                        if(dp0.size() > 0)
                            pend0 = p0 + dp0[1];
                    }
                        break;

                }

                if(j!=3 && j!= 7){//different than the crest
                    vtkIdType idp0end = polypoints->InsertNextPoint(pend0[0], pend0[1], pend0[2]);
                    vtkSmartPointer<vtkLine> line = vtkSmartPointer<vtkLine>::New();
                    line->GetPointIds()->SetId(0, idp0);
                    line->GetPointIds()->SetId(1, idp0end);

                    vtkCellArray* cells = (vtkCellArray*)m_PrimitivesCollectionCells->GetItemAsObject(j);
                    cells->InsertNextCell(line);

                }

            }

        }

    }

    vtkSRep::VectorVNLType crestpositions = input->GetCrestMedialAtoms();
    vtkSRep::VectorVNLType crestderivatives = input->GetCrestMedialAtomsDerivatives();

    for(unsigned i = 0; i < crestderivatives.size(); i++){
        vtkPoints* polypoints = (vtkPoints*) m_PrimitivesCollectionPoints->GetItemAsObject(7);//7 is for the crestderivatives

        vtkSRep::VNLType p0 = crestpositions[i];
        vtkSRep::VNLType p0end = p0 + crestderivatives[i];

        vtkIdType idp0 = polypoints->InsertNextPoint(p0[0], p0[1], p0[2]);
        vtkIdType idp0end = polypoints->InsertNextPoint(p0end[0], p0end[1], p0end[2]);

        vtkSmartPointer<vtkLine> line = vtkSmartPointer<vtkLine>::New();
        line->GetPointIds()->SetId(0, idp0);
        line->GetPointIds()->SetId(1, idp0end);
        vtkCellArray* cells = (vtkCellArray*)m_PrimitivesCollectionCells->GetItemAsObject(7);
        cells->InsertNextCell(line);
    }


    for(int i = 0; i < numprimitives; i++){

        vtkPolyData* poly = (vtkPolyData*) m_PrimitivesCollectionPolyData->GetItemAsObject(i);

        vtkPoints* polypoints = (vtkPoints*) m_PrimitivesCollectionPoints->GetItemAsObject(i);
        vtkCellArray* polycellarray = (vtkCellArray*)m_PrimitivesCollectionCells->GetItemAsObject(i);

        poly->SetPoints(polypoints);

        switch(i){
            case 0://vertex
                poly->SetVerts(polycellarray);
            default://lines primitives
                poly->SetLines(polycellarray);
        }

        vtkSmartPointer<vtkPolyDataMapper> polymapper = vtkSmartPointer<vtkPolyDataMapper>::New();
        polymapper->SetInput(poly);

        vtkSmartPointer<vtkActor>  polyactor = vtkSmartPointer<vtkActor>::New();
        polyactor->SetMapper(polymapper);

        m_PrimitivesCollectionActors->AddItem(polyactor);


        switch(i){
            case 0://vertex
            {
                polyactor->GetProperty()->SetPointSize(15);
                polyactor->GetProperty()->SetColor(1,1,0);//yellow
            }
                break;
            case 1://top spoke
            {
                polyactor->GetProperty()->SetLineWidth(8);
                polyactor->GetProperty()->SetColor(1,0,1);//magenta

            }
                break;
            case 2://bottom spoke
            {
                polyactor->GetProperty()->SetLineWidth(8);
                polyactor->GetProperty()->SetColor(0,1,1);//cyan
            }
                break;
            case 3://crest spoke
            {
                polyactor->GetProperty()->SetLineWidth(8);
                polyactor->GetProperty()->SetColor(0,1,1);//red
            }
                break;
            case 4://sheet normal
            {
                polyactor->GetProperty()->SetLineWidth(8);
                polyactor->GetProperty()->SetColor(0,0,1);//blue
            }
                break;
            case 5://uderivatives
            {
                polyactor->GetProperty()->SetLineWidth(3);
                polyactor->GetProperty()->SetColor(0.8,0.5,0.2);//kindabrown
            }
                break;
            case 6://vderivatives
            {
                polyactor->GetProperty()->SetLineWidth(3);
                polyactor->GetProperty()->SetColor(0.3,0.8,0.4);//kinda green
            }
                break;
            case 7://crestderivatives
            {
                polyactor->GetProperty()->SetLineWidth(3);
                polyactor->GetProperty()->SetColor(0.8,0.2,0.5);
            }
                break;

        }

    }


}
