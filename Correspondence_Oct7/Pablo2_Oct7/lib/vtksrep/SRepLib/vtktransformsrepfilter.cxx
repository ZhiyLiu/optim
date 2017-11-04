#include "vtktransformsrepfilter.h"
#include "vtkObjectFactory.h"


#include "vtkInformationVector.h"
#include "vtkInformation.h"

#include "vtkAbstractTransform.h"
#include "vtkPointData.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkTransform.h"
#include "vnl/vnl_cross.h"

#define PI 3.14159265

vtkStandardNewMacro(vtkTransformSrepFilter);

vtkTransformSrepFilter::vtkTransformSrepFilter()
    : vtkTransformPolyDataFilter()
{
}

vtkTransformSrepFilter::~vtkTransformSrepFilter()
{

}

int vtkTransformSrepFilter::RequestData(
vtkInformation *request,
vtkInformationVector **inputVector,
vtkInformationVector *outputVector) {

    int ret = vtkTransformPolyDataFilter::RequestData(request, inputVector, outputVector);

    vtkInformation *outInfo = outputVector->GetInformationObject(0);

    // get the info objects
    vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);

    vtkPolyData *output = dynamic_cast<vtkPolyData*>(outInfo->Get(vtkDataObject::DATA_OBJECT()));
    // get the input and output
    vtkSRep *input = dynamic_cast<vtkSRep*>(inInfo->Get(vtkDataObject::DATA_OBJECT()));

    m_SRepOutput = vtkSmartPointer<vtkSRep>::New();
    m_SRepOutput->DeepCopy(input);//copy s-rep data
    m_SRepOutput->DeepCopy(output);//copy poly data

    for(unsigned i = 0; i < m_SRepOutput->GetNumberOfPoints(); i++){

        vtkSRep::VNLType d0 = m_SRepOutput->GetDerivatives(i)[0];
        d0.normalize();

        vtkSRep::VNLType norm(3);
        norm.fill(0);
        if(m_SRepOutput->GetSpokes(i).size() > 1){
            norm = vnl_cross_3d(m_SRepOutput->GetSpoke(i, 1, true), m_SRepOutput->GetSpoke(i, 0, true));
            norm.normalize();//this normal is the normal of the wheel
        }

        vtkSRep::VNLType rotvect = vnl_cross_3d(d0, norm);
        //rotvect.normalize();

        vtkSmartPointer<vtkTransform> transform = vtkSmartPointer<vtkTransform>::New();
        double angle = acos(dot_product(norm, d0))*180.0/PI;
        transform->RotateWXYZ(angle, rotvect[0], rotvect[1], rotvect[2]);

        vtkMatrix4x4* rotmatvtk = transform->GetMatrix();
        //rotmatvtk->Print(cout);
        vnl_matrix<double> rotmat((*rotmatvtk)[0], 4, 4);

        rotmat = rotmat.get_n_columns(0, 3).get_n_rows(0, 3);
        //cout<<"vnl: "<<rotmat<<endl;

        vtkSRep::VNLType sheetnorm = m_SRepOutput->GetMedialSheetNormal(i)*rotmat;
        m_SRepOutput->SetMedialSheetNormal(i, sheetnorm);

        for(unsigned j = 0; j < m_SRepOutput->GetSpokes(i).size(); j++){

            vtkSRep::VNLType rotspoke = m_SRepOutput->GetSpoke(i,j) * rotmat;
            m_SRepOutput->SetSpoke(i, j, rotspoke);

        }

    }

    /*for(unsigned i = 0; i < m_SRepOutput->GetNumberOfPoints(); i++){

        vtkSRep::VNLType p0 = m_SRepOutput->GetPointVNL(i);

        for(unsigned j = 0; j < m_SRepOutput->GetSpokes(i).size(); j++){

            vtkSRep::VNLType p1 = m_SRepOutput->GetSpoke(i, j) + p0;

            double *pout;
            pout = this->Transform->TransformDoublePoint(p1.data_block());

            p1 = vtkSRep::VNLType(pout, 3) - p0;

            m_SRepOutput->SetSpoke(i, j, p1);
        }

    }*/

    return ret;
}

void vtkTransformSrepFilter::PrintSelf(ostream& os, vtkIndent indent)
{
  vtkTransformPolyDataFilter::PrintSelf(os,indent);

  //os << indent << "Transform: " << this->Transform << "\n";
}
