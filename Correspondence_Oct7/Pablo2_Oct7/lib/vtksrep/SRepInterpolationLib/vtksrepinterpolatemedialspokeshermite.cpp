#include "vtksrepinterpolatemedialspokeshermite.h"


#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkQuad.h"
#include "vtkTriangle.h"

#include "vtkCellArray.h"
#include "vtkPoints.h"
#include "vtkPointData.h"
#include "vtkCleanPolyData.h"

#include "vtksrepinterpolatemedialsheet.h"

#include "vnl/algo/vnl_symmetric_eigensystem.h"
#include "vnl/algo/vnl_determinant.h"

#include "vtkTriangleFilter.h"
#include "vtkDecimatePro.h"
#include "vtkSmoothPolyDataFilter.h"
#include "vtkGenericCell.h"

#include "vtkObjectFactory.h"
vtkStandardNewMacro(vtkSRepInterpolateMedialSpokesHermite);

vtkSRepInterpolateMedialSpokesHermite::vtkSRepInterpolateMedialSpokesHermite()
{
    m_InterpolationLevel = 0;

    m_AtomId = -1;
    m_SpokeType = -1;
    m_Gamma_u = 1;
    m_Gamma_v = 1;
    m_Input = 0;
}


vtkSRepInterpolateMedialSpokesHermite::~vtkSRepInterpolateMedialSpokesHermite(){

    Clear();

}

void vtkSRepInterpolateMedialSpokesHermite::Clear(){
    for(unsigned i = 0; i < m_HermitMatrices.size(); i++){
        HermiteMatrixType H = m_HermitMatrices[i];
        for(unsigned j = 0; j < H.size(); j++){
            for(unsigned k = 0; k < 4; k++){
                delete[] H[j][k];
            }
            delete[] H[j];
        }
        H.clear();

        /*H = m_HermitMatricesRadius[i];
        for(unsigned j = 0; j < H.size(); j++){
            for(unsigned k = 0; k < 4; k++){
                delete[] H[j][k];
            }
            delete[] H[j];
        }
        H.clear();*/
    }

    m_HermitMatrices.clear();
    //m_HermitMatricesRadius.clear();
}

// Superclass method to update the pipeline
int vtkSRepInterpolateMedialSpokesHermite::RequestData(vtkInformation* request,
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


    vtkSmartPointer<vtkSRepInterpolateMedialSheet> medialinterpolator = vtkSmartPointer<vtkSRepInterpolateMedialSheet>::New();
    medialinterpolator->SetInput(input);
    medialinterpolator->SetInterpolationLevel(m_InterpolationLevel);
    medialinterpolator->Update();

    vtkDataArray* labelarray = input->GetPointData()->GetScalars();
    vtkSmartPointer< vtkDataArray > interpolatedlabelarraysrep = 0;
    vtkSmartPointer< vtkDataArray > interpolatedlabelarraysurface = 0;

    if(labelarray){
        interpolatedlabelarraysrep = (vtkDataArray*)vtkDataArray::CreateArray(VTK_DOUBLE);
        interpolatedlabelarraysurface = (vtkDataArray*)vtkDataArray::CreateArray(VTK_DOUBLE);
    }

    vector< unsigned > vectsides;
    if(m_SpokeType == vtkSRep::TOP_SPOKE){
        vectsides.push_back(vtkSRep::TOP_SPOKE);
    }else if(m_SpokeType == vtkSRep::BOTTOM_SPOKE){
        vectsides.push_back(vtkSRep::BOTTOM_SPOKE);
    }else{
        vectsides.push_back(vtkSRep::TOP_SPOKE);
        vectsides.push_back(vtkSRep::BOTTOM_SPOKE);
    }

    Clear();

    m_HermitMatrices.push_back(HermiteMatrixType());
    m_HermitMatrices.push_back(HermiteMatrixType());

    /*m_HermitMatricesRadius.push_back(HermiteMatrixType());
    m_HermitMatricesRadius.push_back(HermiteMatrixType());*/


    double step = pow((double)2, (double)m_InterpolationLevel);


    m_SRepOutput = vtkSmartPointer<vtkSRep>::New();
    vtkSmartPointer<vtkCellArray> srepoutinterpolatedcellarray = vtkSmartPointer<vtkCellArray>::New();
    vtkSmartPointer<vtkPoints> srepoutinterpolatedpoints = vtkSmartPointer<vtkPoints>::New();
    vtkSRep::RadiusVectorType srepoutallradius;
    vtkSRep::SpokesVectorType srepoutallspokes;
    typedef vector< vtkIdType > VectorIdType;

    vtkSRep::VectorIdsType vectcellsids;
    if(m_AtomId == -1){
        for(unsigned i = 0; i < input->GetNumberOfCells(); i++){
            vtkSmartPointer<vtkGenericCell> cell = vtkSmartPointer<vtkGenericCell>::New();
            input->GetCell(i, cell);
            if(cell->GetNumberOfPoints() == 4){
                vectcellsids.push_back(i);
            }
        }
    }else{
        vtkSmartPointer<vtkIdList> idlist = vtkSmartPointer<vtkIdList>::New();
        input->GetPointCells(m_AtomId, idlist);
        for(unsigned i = 0; i < idlist->GetNumberOfIds(); i++){
            vtkSmartPointer<vtkGenericCell> cell = vtkSmartPointer<vtkGenericCell>::New();
            input->GetCell(idlist->GetId(i), cell);
            if(cell->GetNumberOfPoints() == 4){
                vectcellsids.push_back(idlist->GetId(i));
            }
        }
    }

    for(unsigned i = 0; i < vectcellsids.size(); i++){

        vtkSmartPointer<vtkGenericCell> cell = vtkSmartPointer<vtkGenericCell>::New();
        input->GetCell(vectcellsids[i], cell);

        vtkCell* quad = cell;

        vtkIdType id0 = quad->GetPointIds()->GetId(0);
        vtkIdType id1 = quad->GetPointIds()->GetId(1);
        vtkIdType id2 = quad->GetPointIds()->GetId(2);
        vtkIdType id3 = quad->GetPointIds()->GetId(3);

        for(unsigned sides = 0; sides < vectsides.size(); sides++){

            VNLType** H;
            H = new VNLType*[4];
            for(unsigned j = 0; j < 4; j++) H[j] = new VNLType[4];

            GetHermiteMatrixSpokes(input, id0, id1, id2, id3, vectsides[sides], H);

            m_HermitMatrices[sides].push_back(H);

            VNLType** Hr;
            Hr = new VNLType*[4];
            for(unsigned j = 0; j < 4; j++) Hr[j] = new VNLType[4];

            //GetHermiteMatrixRadius(input, id0, id1, id2, id3, vectsides[sides], Hr);
            //m_HermitMatricesRadius[sides].push_back(Hr);
        }
    }

    vtkSRep::MapUVCoordToIdType mapuvtoid;

    for(unsigned i = 0; i < vectcellsids.size(); i++){

        vector< VectorIdType > srepinterpolatedids;

        double u0 = 0, uend = 1, stepu = 1.0/step;
        double v0 = 0, vend = 1, stepv = 1.0/step;        

        if(m_AtomId != -1){
            vtkSmartPointer<vtkGenericCell> cell = vtkSmartPointer<vtkGenericCell>::New();
            input->GetCell(vectcellsids[i], cell);

            for(unsigned j = 0; j < cell->GetPointIds()->GetNumberOfIds(); j++){
                if(m_AtomId == cell->GetPointId(j)){
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

        vtkSmartPointer<vtkGenericCell> cell = vtkSmartPointer<vtkGenericCell>::New();
        input->GetCell(vectcellsids[i], cell);

        vtkIdType id0 = cell->GetPointId(0);

        int x = id0%input->GetNumColumns();
        int y = id0/input->GetNumColumns();

        for(double u = u0, curru = 0; u <= uend; u += stepu, curru++){
            srepinterpolatedids.push_back(VectorIdType());
            for(double v = v0 ; v <= vend; v += stepv){

                pair<double, double> uvcoord((double)x + u, (double)y + v);

                if(mapuvtoid.find(uvcoord) == mapuvtoid.end()){


                    double interpolatedlabel = 0;

                    if(labelarray){                        
                        vtkSmartPointer<vtkGenericCell> cell = vtkSmartPointer<vtkGenericCell>::New();
                        input->GetCell(vectcellsids[i], cell);

                        vtkIdType id0 = cell->GetPointId(0);
                        vtkIdType id1 = cell->GetPointId(1);
                        vtkIdType id2 = cell->GetPointId(2);
                        vtkIdType id3 = cell->GetPointId(3);

                        double currentlabel = (1-u)*(1-v)*labelarray->GetTuple1(id0) + u*(1-v)*labelarray->GetTuple1(id3) + (1-u)*(v)*labelarray->GetTuple1(id1) + u*v*labelarray->GetTuple1(id2);

                        double mindist = VTK_DOUBLE_MAX;
                        for(unsigned j = 0; j < cell->GetNumberOfPoints(); j++){
                            vtkIdType id = cell->GetPointId(j);
                            double dist = fabs(labelarray->GetTuple1(id) - currentlabel);
                            if(mindist > dist){
                                mindist = dist;
                                interpolatedlabel = labelarray->GetTuple1(id);
                            }
                        }
                        interpolatedlabelarraysrep->InsertNextTuple1(round(interpolatedlabel));
                    }



                    VNLType interpolatedpoint = medialinterpolator->GetInterpolatedPoint(vectcellsids[i], u, v);

                    vtkIdType tempid = srepoutinterpolatedpoints->InsertNextPoint(interpolatedpoint[0], interpolatedpoint[1], interpolatedpoint[2]);
                    srepinterpolatedids[(int)(curru)].push_back(tempid);
                    mapuvtoid[uvcoord] = tempid;

                    vtkSRep::VectorDoubleType radius;
                    vtkSRep::VectorVNLType spokes;

                    for(unsigned sides = 0; sides < vectsides.size(); sides++){

                        VNLType interpolatedspoke = this->GetInterpolatedSpoke(i, sides, u, v);
                        radius.push_back(interpolatedspoke.magnitude());
                        spokes.push_back(interpolatedspoke.normalize());

                    }
                    srepoutallradius.push_back(radius);
                    srepoutallspokes.push_back(spokes);
                }else{
                    vtkIdType id = mapuvtoid[uvcoord];
                    srepinterpolatedids[(int)(curru)].push_back(id);
                }
            }
        }

        for(int j = 0; j < (int)srepinterpolatedids.size() - 1; j++){
            for(int k = 0; k < (int)srepinterpolatedids[j].size() - 1; k++){
                vtkSmartPointer<vtkQuad> quad = vtkSmartPointer<vtkQuad>::New();

                vtkIdType idinterpolated0 = srepinterpolatedids[j][k];
                vtkIdType idinterpolated1 = srepinterpolatedids[j+1][k];
                vtkIdType idinterpolated2 = srepinterpolatedids[j+1][k+1];
                vtkIdType idinterpolated3 = srepinterpolatedids[j][k+1];

                //cout<<idinterpolated0<<" "<<idinterpolated1<<" "<<idinterpolated2<<" "<<idinterpolated3<<endl;

                quad->GetPointIds()->SetId(0, idinterpolated0);
                quad->GetPointIds()->SetId(1, idinterpolated1);
                quad->GetPointIds()->SetId(2, idinterpolated2);
                quad->GetPointIds()->SetId(3, idinterpolated3);

                srepoutinterpolatedcellarray->InsertNextCell(quad);

            }
        }
        for(unsigned j = 0; j < srepinterpolatedids.size(); j++) srepinterpolatedids[j].clear();
        srepinterpolatedids.clear();
    }


    m_SRepOutput->SetPoints(srepoutinterpolatedpoints);
    m_SRepOutput->SetPolys(srepoutinterpolatedcellarray);
    m_SRepOutput->SetAllSpokes(srepoutallspokes);
    m_SRepOutput->SetAllRadius(srepoutallradius);
    if(labelarray){
        m_SRepOutput->GetPointData()->SetScalars(interpolatedlabelarraysrep);
    }
    //m_SRepOutput->SetUVToIdMap(mapuvtoid);
    //m_SRepOutput->GenerateIdToUVMap();



    for(unsigned i = 0; i < m_SRepOutput->GetNumberOfCells(); i++){
        vtkSmartPointer<vtkQuad> quad0 = vtkSmartPointer<vtkQuad>::New();
        m_SRepOutput->GetCell(i, (vtkGenericCell*)(quad0.GetPointer()));

        for(unsigned k = 0; k < vectsides.size(); k++){

            vtkSmartPointer<vtkQuad> quad = vtkSmartPointer<vtkQuad>::New();

            for(unsigned j = 0; j < quad0->GetPointIds()->GetNumberOfIds(); j++){
                vtkIdType id0 = quad0->GetPointIds()->GetId(j);

                double *temp = m_SRepOutput->GetPoint(id0);
                VNLType p0(temp, 3);


                VNLType sp0 = m_SRepOutput->GetSpoke(id0, k);

                p0 += sp0;

                vtkIdType idinterpolated = interpolatedpoints->InsertNextPoint(p0[0],p0[1],p0[2]);

                 quad->GetPointIds()->SetId(j, idinterpolated);

                 if(labelarray){
                    interpolatedlabelarraysurface->InsertNextTuple1(interpolatedlabelarraysrep->GetTuple1(id0));
                 }
            }

            if(false){
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

                interpolatedcellarray->InsertNextCell(cell0);
                interpolatedcellarray->InsertNextCell(cell1);
            }else{
                interpolatedcellarray->InsertNextCell(quad);
            }
        }
    }


    output->SetPoints(interpolatedpoints);
    output->SetPolys(interpolatedcellarray);
    output->GetPointData()->SetScalars(interpolatedlabelarraysurface);

    /*vtkSmartPointer<vtkDecimatePro> decimate = vtkSmartPointer<vtkDecimatePro>::New();
    decimate->SetInput(output);
    decimate->SetTargetReduction(.2);
    decimate->Update();
    output->DeepCopy(decimate->GetOutput());0

    vtkSmartPointer<vtkSmoothPolyDataFilter> smoothFilter = vtkSmartPointer<vtkSmoothPolyDataFilter>::New();
    //smoothFilter->SetInput(decimate->GetOutput());
    smoothFilter->SetInput(output);
    smoothFilter->SetNumberOfIterations(2);
    smoothFilter->GenerateErrorScalarsOff();
    smoothFilter->GenerateErrorVectorsOff();    
    smoothFilter->Update();

    output = smoothFilter->GetOutput();*/


    return 1;
}

double vtkSRepInterpolateMedialSpokesHermite::GetInterpolatedLabel(double u, double v, vtkIdType cellid){


    vtkSmartPointer<vtkGenericCell> cell = vtkSmartPointer<vtkGenericCell>::New();
    m_Input->GetCell(cellid, cell);
    vtkDataArray* labelarray = m_Input->GetPointData()->GetScalars();

    double interpolatedlabel = 0;

    if(labelarray && cell){
        vtkIdType id0 = cell->GetPointId(0);
        vtkIdType id1 = cell->GetPointId(1);
        vtkIdType id2 = cell->GetPointId(2);
        vtkIdType id3 = cell->GetPointId(3);

        double currentlabel = (1-u)*(1-v)*labelarray->GetTuple1(id0) + u*(1-v)*labelarray->GetTuple1(id3) + (1-u)*(v)*labelarray->GetTuple1(id1) + u*v*labelarray->GetTuple1(id2);

        double mindist = VTK_DOUBLE_MAX;
        for(unsigned j = 0; j < cell->GetNumberOfPoints(); j++){
            vtkIdType id = cell->GetPointId(j);
            double dist = fabs(labelarray->GetTuple1(id) - currentlabel);
            if(mindist > dist){
                mindist = dist;
                interpolatedlabel = labelarray->GetTuple1(id);
            }
        }
    }

    return interpolatedlabel;

}

/*
*   \fn Vector3D* getHermiteMatrix(Vector3D p0, Vector3D p1, Vector3D p2, Vector3D p3, Vector3D n0, Vector3D n1, Vector3D n2, Vector3D n3)
*	\brief Calculates the corresponding hermiteMatrix from the 4 control points and the 4 given normals
*	\param Vector3D p0-p1 the control points
*	\param Vector3D n0-n1 the normals of the control points
*/
void vtkSRepInterpolateMedialSpokesHermite::GetHermiteMatrix(VNLType **H, VNLType p0, VNLType p1, VNLType p2, VNLType p3, VectorVNLType dp0, VectorVNLType dp1, VectorVNLType dp2, VectorVNLType dp3, VNLType n0, VNLType n1, VNLType n2, VNLType n3, VNLType k0, VNLType k1, VNLType k2, VNLType k3){

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

    VNLType t2pxx(dp0u.size());
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

void vtkSRepInterpolateMedialSpokesHermite::GetHermiteMatrixSpokes(vtkSRep *input, vtkIdType id0, vtkIdType id1, vtkIdType id2, vtkIdType id3, vtkIdType side, VNLType **H){

    /*VNLType p0 = input->GetSpoke(id0, side, true);
    VNLType p1 = input->GetSpoke(id1, side, true);
    VNLType p2 = input->GetSpoke(id2, side, true);
    VNLType p3 = input->GetSpoke(id3, side, true);*/
    VNLType p0 = input->GetSpoke(id0, side);
    VNLType p1 = input->GetSpoke(id1, side);
    VNLType p2 = input->GetSpoke(id2, side);
    VNLType p3 = input->GetSpoke(id3, side);

    VectorVNLType dp0 = input->GetDerivativesSpoke(id0, side);
    VectorVNLType dp1 = input->GetDerivativesSpoke(id1, side);
    VectorVNLType dp2 = input->GetDerivativesSpoke(id2, side);
    VectorVNLType dp3 = input->GetDerivativesSpoke(id3, side);

    /*VNLType n0 = input->GetMedialSheetNormal(id0);
    VNLType n1 = input->GetMedialSheetNormal(id1);
    VNLType n2 = input->GetMedialSheetNormal(id2);
    VNLType n3 = input->GetMedialSheetNormal(id3);

    GetHermiteMatrix(H, p0, p1, p2, p3, dp0, dp1, dp2, dp3, n0, n1, n2, n3);*/



    H[0][0] = p0;
    H[0][1] = p1;
    H[1][0] = p3;
    H[1][1] = p2;

    H[2][0] = dp0[0];
    H[2][1] = dp1[0];
    H[3][0] = dp3[0];
    H[3][1] = dp2[0];

    H[0][2] = dp0[1];
    H[0][3] = dp1[1];
    H[1][2] = dp3[1];
    H[1][3] = dp2[1];

    VNLType t2pxx(3);
    t2pxx.fill(0);

    H[2][2] = t2pxx;
    H[2][3] = t2pxx;
    H[3][2] = t2pxx;
    H[3][3] = t2pxx;


    
}

/*void vtkSRepInterpolateMedialSpokesHermite::GetHermiteMatrixRadius(vtkSRep* input, vtkIdType id0, vtkIdType id1, vtkIdType id2, vtkIdType id3, vtkIdType side, VNLType **H){

    VNLType p0(1);
    p0[0] = input->GetSpokeRadius(id0, side);
    VNLType p1(1);
    p1[0] = input->GetSpokeRadius(id1, side);
    VNLType p2(1);
    p2[0] = input->GetSpokeRadius(id2, side);
    VNLType p3(1);
    p3[0] = input->GetSpokeRadius(id3, side);

    VectorVNLType dp0 = input->GetDerivativesSpokeRadius(id0, side);
    VectorVNLType dp1 = input->GetDerivativesSpokeRadius(id1, side);
    VectorVNLType dp2 = input->GetDerivativesSpokeRadius(id2, side);
    VectorVNLType dp3 = input->GetDerivativesSpokeRadius(id3, side);

    VNLType n0(1);
    n0.fill(1);
    VNLType n1(1);
    n1.fill(1);
    VNLType n2(1);
    n2.fill(1);
    VNLType n3(1);
    n3.fill(1);


    //GetHermiteMatrix(H, p0, p1, p2, p3, dp0, dp1, dp2, dp3, n0, n1, n2, n3, k0, k1, k2, k3);
    GetHermiteMatrix(H, p0, p1, p2, p3, dp0, dp1, dp2, dp3, n0, n1, n2, n3);

}*/

/*
* \fn void GetInterpolatedPoint(vtkIdType cellid, double u, double v)
* \brief Computes the point using the cubic Hermite weight functions and the two scalar parameters
* \param vtkIdType cellid, cellid in the srep
* \param unsigned side [0, 1] top or bottom
* \param double u scalar for the weight functions
* \param double v second scalar for the weight functions    n
* \pre u and v should be between [0, 1]
*/
vtkSRepInterpolateMedialSpokesHermite::VNLType vtkSRepInterpolateMedialSpokesHermite::GetInterpolatedSpoke(vtkIdType cellid, unsigned side, double u, double v){

    //return GetInterpolatedSpoke(u, v, m_HermitMatrices[side][cellid]).normalize()*GetInterpolatedSpoke(u, v, m_HermitMatricesRadius[side][cellid])[0];
    return GetInterpolatedSpoke(u, v, m_HermitMatrices[side][cellid]);
}

/*
* \fn void GetInterpolatedSpoke(vtkIdType cellid, double u, double v)
* \brief Computes the point using the cubic Hermite weight functions and the two scalar parameters
* \param vtkIdType cellid, cellid in the srep
* \param unsigned spoke position ex: 1-0 [top, bottom]
* \param double u scalar for the weight functions
* \param double v second scalar for the weight functions    n
* \pre u and v should be between [0, 1]
*/
/*vtkSRepInterpolateMedialSpokesHermite::VNLType vtkSRepInterpolateMedialSpokesHermite::GetInterpolatedSpoke(unsigned side, double u, double v){
    double ucoord = u*(m_Input->GetNumColumns()-1);
    double vcoord = v*(m_Input->GetNumRows()-1);
    vtkIdType cellid = m_Input->GetCellId(ucoord, vcoord);

    if(cellid != -1){
        return GetInterpolatedSpoke(cellid, side, ucoord, vcoord);
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
* \param Vector3D H[4][4] cubic hermite matrix
* \param Vector3D p result of the interpolation
* \pre u and v should be between [0, 1]
*/
vtkSRepInterpolateMedialSpokesHermite::VNLType vtkSRepInterpolateMedialSpokesHermite::GetInterpolatedSpoke(double u, double v, VNLType** H){
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
* \fn void GetInterpolatedPoint(vtkIdType cellid, double u, double v)
* \brief Computes the point using the cubic Hermite weight functions and the two scalar parameters
* \param vtkIdType cellid, cellid in the srep
* \param double u scalar for the weight functions
* \param double v second scalar for the weight functions    n
* \pre u and v should be between [0, 1]
*/
vtkSRepInterpolateMedialSpokesHermite::VNLType vtkSRepInterpolateMedialSpokesHermite::GetInterpolatedSpokeDerivative(vtkIdType cellid, unsigned side, double u, double v){

    return GetInterpolatedSpokeDerivative(u, v, m_HermitMatrices[side][cellid]);
    }

/*
* \fn void GetInterpolatedPoint(double u, double v, Vector3D H[4][4], Vector3D p)
* \brief Computes the point using the cubic Hermite weight functions and the two scalar parameters
* \param double u scalar for the weight functions
* \param double v second scalar for the weight functions
* \param Vector3D H[4][4] cubic hermite matrix
* \param Vector3D p result of the interpolation
* \pre u and v should be between [0, 1]
*/
vtkSRepInterpolateMedialSpokesHermite::VNLType vtkSRepInterpolateMedialSpokesHermite::GetInterpolatedSpokeDerivative(double u, double v, VNLType** H){
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
