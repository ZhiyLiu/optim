#include "vtksrepinterpolatemedialsheet.h"


#include "vtkSmartPointer.h"


#include "vtkCellArray.h"
#include "vtkQuad.h"
#include "vtkVertex.h"
#include "vtkGenericCell.h"

#include "vtkCleanPolyData.h"

#include "vtkInformationVector.h"
#include "vtkInformation.h"

#include "vnl/vnl_matrix.h"

#include "vtkObjectFactory.h"
vtkStandardNewMacro(vtkSRepInterpolateMedialSheet);



vtkSRepInterpolateMedialSheet::vtkSRepInterpolateMedialSheet()
{
    m_InterpolationLevel = 0;
    m_AtomId = -1;
    m_Gamma_u = 1;
    m_Gamma_v = 1;
    m_Input = 0;
}

vtkSRepInterpolateMedialSheet::~vtkSRepInterpolateMedialSheet()
{
    for(unsigned i = 0; i < m_HermitMatrices.size(); i++){
        VNLType** H = m_HermitMatrices[i];
        for(unsigned j = 0; j < 4; j++){
            delete[] H[j];
        }
        delete[] H;
    }
    m_HermitMatrices.clear();
}


// Superclass method to update the pipeline
int vtkSRepInterpolateMedialSheet::RequestData(vtkInformation* request,
                        vtkInformationVector** inputVector,
                        vtkInformationVector* outputVector){



    // get the info objects
    vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
    vtkInformation *outInfo = outputVector->GetInformationObject(0);

    // get the input and output
    vtkSRep *input = dynamic_cast<vtkSRep*>(vtkSRep::SafeDownCast(inInfo->Get(vtkSRep::DATA_OBJECT())));
    m_Input = input;
    vtkPolyData *output = vtkPolyData::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));


    vtkSmartPointer<vtkCellArray> interpolatedcellarray = vtkSmartPointer<vtkCellArray>::New();
    vtkSmartPointer<vtkPoints> interpolatedpoints = vtkSmartPointer<vtkPoints>::New();
    //vtkCellArray* interpolatedcellarray = vtkCellArray::New();
   // vtkPoints* interpolatedpoints = vtkPoints::New();

    vtkSmartPointer<vtkIdList> cellsid = vtkSmartPointer<vtkIdList>::New();
    if(m_AtomId != -1){
        input->GetPointCells(m_AtomId, cellsid);
    }else{
        for(unsigned i = 0; i < input->GetNumberOfCells(); i++){
            cellsid->InsertNextId(i);
        }
    }

    for(unsigned i = 0; i < cellsid->GetNumberOfIds(); i++){

        vtkSmartPointer<vtkGenericCell> quad = vtkSmartPointer<vtkGenericCell>::New();
        input->GetCell(cellsid->GetId(i), quad);

        if(quad->GetNumberOfPoints() == 4){

            vtkIdType id0 = quad->GetPointIds()->GetId(0);
            vtkIdType id1 = quad->GetPointIds()->GetId(1);
            vtkIdType id2 = quad->GetPointIds()->GetId(2);
            vtkIdType id3 = quad->GetPointIds()->GetId(3);

            double temp[3];
            input->GetPoint(id0, temp);
            VNLType p0(temp, 3);

            input->GetPoint(id1, temp);
            VNLType p1(temp, 3);

            input->GetPoint(id2, temp);
            VNLType p2(temp, 3);

            input->GetPoint(id3, temp);
            VNLType p3(temp, 3);

            //VectorVNLType test = input->GetDerivatives(15);
            //cout<<"test: "<<test[0]<<" "<<test[1]<<" "<<12<<endl;

            VectorVNLType dp0 = input->GetDerivatives(id0);
            VectorVNLType dp1 = input->GetDerivatives(id1);
            VectorVNLType dp2 = input->GetDerivatives(id2);
            VectorVNLType dp3 = input->GetDerivatives(id3);

            VNLType n0 = input->GetMedialSheetNormal(id0);

            VNLType n1 = input->GetMedialSheetNormal(id1);

            VNLType n2 = input->GetMedialSheetNormal(id2);

            VNLType n3 = input->GetMedialSheetNormal(id3);

            VNLType** H;
            H = new VNLType*[4];
            for(unsigned j = 0; j < 4; j++) H[j] = new VNLType[4];

            GetHermiteMatrix(p0, p1, p2, p3, dp0, dp1, dp2, dp3, n0, n1, n2, n3, H);
            m_HermitMatrices.push_back(H);

            /*cout<<"p: "<<endl<<p0[0]<<" "<<p0[1]<<" "<<p0[2]<<" "<<id0<<endl;
            cout<<p1[0]<<" "<<p1[1]<<" "<<p1[2]<<" "<<id1<<endl;
            cout<<p2[0]<<" "<<p2[1]<<" "<<p2[2]<<" "<<id2<<endl;
            cout<<p3[0]<<" "<<p3[1]<<" "<<p3[2]<<" "<<id3<<endl;
            cout<<"i: "<<endl;*/

            double step = pow((double)2, (double)m_InterpolationLevel);

            typedef vector< vtkIdType > VectorIdType;
            vector< VectorIdType > interpolatedids;

            double u0 = 0, uend = 1, stepu = 1.0/step;
            double v0 = 0, vend = 1, stepv = 1.0/step;

            if(m_AtomId != -1){

                for(unsigned j = 0; j < quad->GetPointIds()->GetNumberOfIds(); j++){
                    if(m_AtomId == quad->GetPointId(j)){
                        if(j == 0){
                            uend = m_Gamma_u;
                            vend = m_Gamma_v;
                        }else if(j == 1){
                            uend = m_Gamma_u;
                            v0 = 1 - m_Gamma_v;
                            vend = 1;
                        }else if(j == 2){
                            u0 = 1 - m_Gamma_u;
                            uend = 1;
                            v0 = 1 - m_Gamma_v;
                            vend = 1;
                        }else if(j == 3){
                            u0 = 1 - m_Gamma_u;
                            uend = 1;
                            vend = m_Gamma_v;
                        }
                        break;
                    }
                }
            }

            for(double u = u0, curru = 0; u <= uend; u += stepu, curru++){
                interpolatedids.push_back(VectorIdType());
                for(double v = v0 ; v <= vend; v += stepv){

                    VNLType interpolated = this->GetInterpolatedPoint(u, v, H);

                    vtkIdType tempid = interpolatedpoints->InsertNextPoint(interpolated[0], interpolated[1], interpolated[2]);

                    //cout<<interpolated[0]<<" "<<interpolated[1]<<" "<<interpolated[2]<<" id: "<<tempid<<endl;

                    interpolatedids[(int)curru].push_back(tempid);

                    /*vtkVertex* vertex = vtkVertex::New();
                    vertex->GetPointIds()->SetId(0, tempid);
                    interpolatedcellarray->InsertNextCell(vertex);*/

                }
            }

            for(unsigned i = 0; i < interpolatedids.size() - 1; i++){
                for(unsigned j = 0; j < interpolatedids[i].size() - 1; j++){

                    vtkSmartPointer<vtkQuad> quad = vtkSmartPointer<vtkQuad>::New();

                    vtkIdType idinterpolated0 = interpolatedids[i][j];
                    vtkIdType idinterpolated1 = interpolatedids[i+1][j];
                    vtkIdType idinterpolated2 = interpolatedids[i+1][j+1];
                    vtkIdType idinterpolated3 = interpolatedids[i][j+1];

                    //cout<<idinterpolated0<<" "<<idinterpolated1<<" "<<idinterpolated2<<" "<<idinterpolated3<<endl;

                    quad->GetPointIds()->SetId(0, idinterpolated0);
                    quad->GetPointIds()->SetId(1, idinterpolated1);
                    quad->GetPointIds()->SetId(2, idinterpolated2);
                    quad->GetPointIds()->SetId(3, idinterpolated3);

                    interpolatedcellarray->InsertNextCell(quad);

                }
            }

            for(unsigned i = 0; i < interpolatedids.size() - 1; i++) interpolatedids[i].clear();
            interpolatedids.clear();
        }
    }

    output->SetPoints(interpolatedpoints);
    output->SetPolys(interpolatedcellarray);
    //output->SetVerts(interpolatedcellarray);

    /*vtkSmartPointer<vtkCleanPolyData> clean = vtkSmartPointer<vtkCleanPolyData>::New();
    clean->SetInput(output);
    clean->PointMergingOn();
    clean->SetTolerance(0.1);
    clean->SetInput(output);
    output = clean->GetOutput();*/


    return 1;

}


/*
* \fn void GetInterpolatedPoint(double u, double v)
* \brief Computes the point using the cubic Hermite weight functions and the two scalar parameters
* \brief The cell id is calculated from u, v
* \param double u scalar for the weight functions
* \param double v second scalar for the weight functions
* \pre u and v should be between [0, 1]
*/
/*vtkSRepInterpolateMedialSheet::VNLType vtkSRepInterpolateMedialSheet::GetInterpolatedPoint(double u, double v){
    double ucoord = u*(m_Input->GetNumColumns()-1);
    double vcoord = v*(m_Input->GetNumRows()-1);
    vtkIdType cellid = m_Input->GetCellId(ucoord, vcoord);

    if(cellid != -1){
        return GetInterpolatedPoint(cellid, ucoord, vcoord);
    }
    VNLType s(3);
    s.fill(0);
    return s;

}*/

/*
* \fn void GetInterpolatedPoint(double u, double v, Vector3D H[4][4], Vector3D p)
* \brief Computes the point using the cubic Hermite weight functions and the two scalar parameters
* \param double u scalar for the weight functions
* \param double v second scalar for the weight functions
* \param cellid the id of the quad, see srep->GetNumberofCells()
* \pre u and v should be between [0, 1]
*/
vtkSRepInterpolateMedialSheet::VNLType vtkSRepInterpolateMedialSheet::GetInterpolatedPoint(vtkIdType cellid, double u, double v){
    return GetInterpolatedPoint(u, v, m_HermitMatrices[cellid]);
}

/* \brief Using Hermite interpolation calculates the derivative in the u
          direction
*/
vtkSRepInterpolateMedialSheet::VNLType vtkSRepInterpolateMedialSheet::GetInterpolatedDu(vtkIdType cellid, double u, double v){
    return GetInterpolatedDu(u, v, m_HermitMatrices[cellid]);
}

/* \brief Using Hermite interpolation calculates the derivative in the v
          direction
*/
vtkSRepInterpolateMedialSheet::VNLType vtkSRepInterpolateMedialSheet::GetInterpolatedDv(vtkIdType cellid, double u, double v){
    return GetInterpolatedDv(u, v, m_HermitMatrices[cellid]);
}


/*
*   \fn Vector3D* getHermiteMatrix(Vector3D p0, Vector3D p1, Vector3D p2, Vector3D p3, Vector3D n0, Vector3D n1, Vector3D n2, Vector3D n3)
*	\brief Calculates the corresponding hermiteMatrix from the 4 control points and the 4 given normals
*	\param Vector3D p0-p1 the control points
*	\param Vector3D n0-n1 the normals of the control points
*/
void vtkSRepInterpolateMedialSheet::GetHermiteMatrix(VNLType p0, VNLType p1, VNLType p2, VNLType p3, VectorVNLType dp0, VectorVNLType dp1, VectorVNLType dp2, VectorVNLType dp3, VNLType n0, VNLType n1, VNLType n2, VNLType n3, VNLType **H){

    //derivates from points using adjacent points
    VNLType dp0u = dp0[0];
    VNLType dp0v = dp0[1];

    VNLType dp1u = dp1[0];
    VNLType dp1v = dp1[1];

    VNLType dp2u = dp2[0];
    VNLType dp2v = dp2[1];

    VNLType dp3u = dp3[0];
    VNLType dp3v = dp3[1];

    //tangent Vector3Ds to normal plane described by the normal n
    VNLType tp0u = dp0u - dot_product(dp0u, n0) * n0;//v - ( v dot n ) n projection of point to plane
    VNLType tp0v = dp0v - dot_product(dp0v, n0) * n0;

    VNLType tp1u = dp1u - dot_product(dp1u, n1) * n1;
    VNLType tp1v = dp1v - dot_product(dp1v, n1) * n1;

    VNLType tp2u = dp2u - dot_product(dp2u, n2) * n2;
    VNLType tp2v = dp2v - dot_product(dp2v, n2) * n2;

    VNLType tp3u = dp3u - dot_product(dp3u, n3) * n3;
    VNLType tp3v = dp3v - dot_product(dp3v, n3) * n3;

    VNLType t2pxx(3);
    t2pxx.fill(0);

    //construction of the hermite matrix
    H[0][0] = p0;
    H[0][1] = p1;
    H[1][0] = p3;
    H[1][1] = p2;

    H[2][0] = tp0u;
    H[2][1] = tp1u;
    H[3][0] = tp3u;
    H[3][1] = tp2u;

    H[0][2] = tp0v;
    H[0][3] = tp1v;
    H[1][2] = tp3v;
    H[1][3] = tp2v;

    H[2][2] = t2pxx;
    H[2][3] = t2pxx;
    H[3][2] = t2pxx;
    H[3][3] = t2pxx;
}

/*
    * \fn void GetPoint(double u, double v, Vector3D H[4][4], Vector3D p)
    * \brief Computes the point using the cubic Hermite weight functions and the two scalar parameters
    * \param double u scalar for the weight functions
    * \param double v second scalar for the weight functions
    * \param Vector3D H[4][4] cubic hermite matrix
    * \param Vector3D p result of the interpolation
    */
vtkSRepInterpolateMedialSheet::VNLType vtkSRepInterpolateMedialSheet::GetInterpolatedPoint(double u, double v, VNLType **H){


    double Hu[4];
    double Hv[4];
    VNLType HuH[4];

    //weight functions
    Hu[0] = 2*pow(u, 3) - 3*pow(u, 2) + 1;
    Hu[1] = -2*pow(u, 3) + 3*pow(u, 2);
    Hu[2] = pow(u, 3) - 2*pow(u, 2) + u;
    Hu[3] = pow(u, 3) - pow(u, 2);

    Hv[0] = 2*pow(v, 3) - 3*pow(v, 2) + 1;
    Hv[1] = -2*pow(v, 3) + 3*pow(v, 2);
    Hv[2] = pow(v, 3) - 2*pow(v, 2) +v;
    Hv[3] = pow(v, 3) - pow(v, 2);

    //calculation of point using the weight functions
    HuH[0] = Hu[0] * H[0][0] + Hu[1] * H[1][0] + Hu[2] * H[2][0] + Hu[3] * H[3][0];
    HuH[1] = Hu[0] * H[0][1] + Hu[1] * H[1][1] + Hu[2] * H[2][1] + Hu[3] * H[3][1];
    HuH[2] = Hu[0] * H[0][2] + Hu[1] * H[1][2] + Hu[2] * H[2][2] + Hu[3] * H[3][2];
    HuH[3] = Hu[0] * H[0][3] + Hu[1] * H[1][3] + Hu[2] * H[2][3] + Hu[3] * H[3][3];

    return Hv[0] * HuH[0] + Hv[1] * HuH[1] + Hv[2] * HuH[2] + Hv[3] * HuH[3];


}
/*
    *
    * \brief Computes the derivative of the point in the u direction
    * \param double u scalar for the weight functions
    * \param double v second scalar for the weight functions
    * \param Vector3D H[4][4] cubic hermite matrix
    * \param Vector3D p result of the interpolation
    */
vtkSRepInterpolateMedialSheet::VNLType vtkSRepInterpolateMedialSheet::GetInterpolatedDu(double u, double v, VNLType **H){


    double Hu[4];
    double Hv[4];
    VNLType HuH[4];

    //weight functions
    Hu[0] = 6*u*u - 6*u;
    Hu[1] = -6*u*u + 6*u;
    Hu[2] = 3*u*u - 4*u + 1;
    Hu[3] = 3*u*u - 2*u;

    Hv[0] = 2*v*v*v-3*v*v+1;
    Hv[1] = -2*v*v*v+3*v*v;
    Hv[2] = v*v*v-2*v*v+v;
    Hv[3] = v*v*v-v*v;

    //calculation of point using the weight functions
    HuH[0] = Hu[0] * H[0][0] + Hu[1] * H[1][0] + Hu[2] * H[2][0] + Hu[3] * H[3][0];
    HuH[1] = Hu[0] * H[0][1] + Hu[1] * H[1][1] + Hu[2] * H[2][1] + Hu[3] * H[3][1];
    HuH[2] = Hu[0] * H[0][2] + Hu[1] * H[1][2] + Hu[2] * H[2][2] + Hu[3] * H[3][2];
    HuH[3] = Hu[0] * H[0][3] + Hu[1] * H[1][3] + Hu[2] * H[2][3] + Hu[3] * H[3][3];

    return Hv[0] * HuH[0] + Hv[1] * HuH[1] + Hv[2] * HuH[2] + Hv[3] * HuH[3];


}

/*
    * \brief Computes the derivative of the point in the v direction using the cubic Hermite weight
    * \param double u scalar for the weight functions
    * \param double v second scalar for the weight functions
    * \param Vector3D H[4][4] cubic hermite matrix
    * \param Vector3D p result of the interpolation
    */
vtkSRepInterpolateMedialSheet::VNLType vtkSRepInterpolateMedialSheet::GetInterpolatedDv(double u, double v, VNLType **H){


    double Hu[4];
    double Hv[4];
    VNLType HuH[4];

    //weight functions
    Hu[0] = 2*u*u*u-3*u*u+1;
    Hu[1] = -2*u*u*u+3*u*u;
    Hu[2] = u*u*u - 2*u*u + u;
    Hu[3] = u*u*u - u*u;

    //weight functions
    Hv[0] = 6*v*v - 6*v;
    Hv[1] = -6*v*v + 6*v;
    Hv[2] = 3*v*v - 4*v + 1;
    Hv[3] = 3*v*v - 2*v;

    //calculation of point using the weight functions
    HuH[0] = Hu[0] * H[0][0] + Hu[1] * H[1][0] + Hu[2] * H[2][0] + Hu[3] * H[3][0];
    HuH[1] = Hu[0] * H[0][1] + Hu[1] * H[1][1] + Hu[2] * H[2][1] + Hu[3] * H[3][1];
    HuH[2] = Hu[0] * H[0][2] + Hu[1] * H[1][2] + Hu[2] * H[2][2] + Hu[3] * H[3][2];
    HuH[3] = Hu[0] * H[0][3] + Hu[1] * H[1][3] + Hu[2] * H[2][3] + Hu[3] * H[3][3];

    return Hv[0] * HuH[0] + Hv[1] * HuH[1] + Hv[2] * HuH[2] + Hv[3] * HuH[3];


}
