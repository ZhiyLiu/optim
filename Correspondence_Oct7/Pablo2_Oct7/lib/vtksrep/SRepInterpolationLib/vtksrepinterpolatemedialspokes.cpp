#include "vtksrepinterpolatemedialspokes.h"


#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkQuad.h"

#include "vtkCellArray.h"
#include "vtkPoints.h"

#include "vnl/algo/vnl_symmetric_eigensystem.h"
#include "vnl/algo/vnl_determinant.h"

#include "vtkObjectFactory.h"
vtkStandardNewMacro(vtkSRepInterpolateMedialSpokes);

vtkSRepInterpolateMedialSpokes::vtkSRepInterpolateMedialSpokes()
{
    m_InterpolationLevel = 0;
}


vtkSRepInterpolateMedialSpokes::~vtkSRepInterpolateMedialSpokes(){

}

// Superclass method to update the pipeline
int vtkSRepInterpolateMedialSpokes::RequestData(vtkInformation* request,
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

    vector< unsigned > vectsides;
    //vectsides.push_back(vtkSRep::TOP_SPOKE);
    vectsides.push_back(vtkSRep::BOTTOM_SPOKE);

    vtkSmartPointer<vtkSRepInterpolateMedialSheet> medialinterpolator = vtkSmartPointer<vtkSRepInterpolateMedialSheet>::New();
    medialinterpolator->SetInput(input);
    medialinterpolator->SetInterpolationLevel(m_InterpolationLevel);
    medialinterpolator->Update();

    for(unsigned i = 0; i < input->GetNumberOfCells(); i++){
    //for(unsigned i = 0; i < 1; i++){

        vtkSmartPointer<vtkQuad> quad = vtkSmartPointer<vtkQuad>::New();
        quad->DeepCopy(input->GetCell(i));

        for(unsigned side = 0; side < vectsides.size(); side++){
            VectorVNLMatrixType srads;

            for(unsigned j = 0; j < quad->GetNumberOfPoints(); j++){
                vtkIdType pointid = quad->GetPointId(j);
                srads.push_back(input->GetSRadMatrix(pointid, vectsides[side]));
            }

            VectorVectorVNLType vectinterpolatedpoints;
            VectorVectorVNLType vectinterpolatedspokes;

            InterpolateSpokes(i, quad->GetPointId(0), srads, vectsides[side], medialinterpolator, vectinterpolatedpoints, vectinterpolatedspokes);

            VectorVectorIdsType vectpointsids;

            for(unsigned x = 0; x < vectinterpolatedpoints.size(); x++){
                vectpointsids.push_back(VectorIdsType());
                for(unsigned y = 0; y < vectinterpolatedpoints[x].size(); y++){

                    VNLType p = vectinterpolatedpoints[x][y] + vectinterpolatedspokes[x][y];
                    vtkIdType id = interpolatedpoints->InsertNextPoint(p[0], p[1], p[2]);
                    vectpointsids[x].push_back(id);

                }
            }

            for(unsigned x = 0; x < vectpointsids.size() - 1; x++){
                for(unsigned y = 0; y < vectpointsids[x].size() - 1; y++){

                    vtkQuad* quad = vtkQuad::New();
                    quad->GetPointIds()->SetId(0, vectpointsids[x][y]);
                    quad->GetPointIds()->SetId(1, vectpointsids[x+1][y]);
                    quad->GetPointIds()->SetId(2, vectpointsids[x+1][y+1]);
                    quad->GetPointIds()->SetId(3, vectpointsids[x][y+1]);

                    interpolatedcellarray->InsertNextCell(quad);

                }
            }

            for(unsigned x = 0; x < vectpointsids.size() - 1; x++){
                vectinterpolatedpoints[x].clear();
                vectinterpolatedspokes[x].clear();
                vectpointsids[x].clear();
            }

            vectinterpolatedpoints.clear();
            vectinterpolatedspokes.clear();
            vectpointsids.clear();
            srads.clear();
        }
    }

    output->SetPoints(interpolatedpoints);
    output->SetPolys(interpolatedcellarray);
    return 1;
}


/*
              3_______2
  quad =    u /      /
            0/______/1
                v

  */
void vtkSRepInterpolateMedialSpokes::InterpolateSpokes(vtkIdType cellid, vtkIdType pointid, VectorVNLMatrixType srads, vtkIdType side, vtkSmartPointer<vtkSRepInterpolateMedialSheet> medialinterpolator, VectorVectorVNLType& interpolatedpoints, VectorVectorVNLType &interpolatedspokes){

    // get the input and output
    vtkSRep *input = m_Input;

    VNLMatrixType rSrad0 = srads[0];
    VNLMatrixType rSrad1 = srads[1];
    VNLMatrixType rSrad2 = srads[2];
    VNLMatrixType rSrad3 = srads[3];

    cout<<"rSrad0 = "<<endl<<rSrad0<<endl;
    cout<<"rSrad1 = "<<endl<<rSrad1<<endl;
    cout<<"rSrad2 = "<<endl<<rSrad2<<endl;
    cout<<"rSrad3 = "<<endl<<rSrad3<<endl;

    vnl_symmetric_eigensystem<double> eigensystem0(rSrad0);
    vnl_symmetric_eigensystem<double> eigensystem1(rSrad1);
    vnl_symmetric_eigensystem<double> eigensystem2(rSrad2);
    vnl_symmetric_eigensystem<double> eigensystem3(rSrad3);

    //cout<<eigensystem0.get_eigenvalue(0)<<" "<<eigensystem0.get_eigenvalue(1)<<endl;
    //cout<<"o= "<<eigensystem0.get_eigenvector(0)<<" . "<<eigensystem0.get_eigenvector(1)<<endl;

    /*cout<<eigensystem1.get_eigenvalue(0)<<" "<<eigensystem1.get_eigenvalue(1)<<endl;
    cout<<eigensystem2.get_eigenvalue(0)<<" "<<eigensystem2.get_eigenvalue(1)<<endl;
    cout<<eigensystem3.get_eigenvalue(0)<<" "<<eigensystem3.get_eigenvalue(1)<<endl;*/

    double logLam0_0 = GetLogLambda(eigensystem0.get_eigenvalue(0));
    double logLam0_1 = GetLogLambda(eigensystem0.get_eigenvalue(1));
    double logLam1_0 = GetLogLambda(eigensystem1.get_eigenvalue(0));
    double logLam1_1 = GetLogLambda(eigensystem1.get_eigenvalue(1));
    double logLam2_0 = GetLogLambda(eigensystem2.get_eigenvalue(0));
    double logLam2_1 = GetLogLambda(eigensystem2.get_eigenvalue(1));
    double logLam3_0 = GetLogLambda(eigensystem3.get_eigenvalue(0));
    double logLam3_1 = GetLogLambda(eigensystem3.get_eigenvalue(1));

    VNLType e0_0 = eigensystem0.get_eigenvector(0);
    VNLType e0_1 = eigensystem0.get_eigenvector(1);
    VNLType e1_0 = eigensystem1.get_eigenvector(0);
    VNLType e1_1 = eigensystem1.get_eigenvector(1);
    VNLType e2_0 = eigensystem2.get_eigenvector(0);
    VNLType e2_1 = eigensystem2.get_eigenvector(1);
    VNLType e3_0 = eigensystem3.get_eigenvector(0);
    VNLType e3_1 = eigensystem3.get_eigenvector(1);


    double r = input->GetSpokeRadius(pointid, side);
    VNLType U = input->GetSpoke(pointid, side, false);
    VNLType x0 = medialinterpolator->GetInterpolatedPoint(cellid, 0,0);


    double step = pow((double)2, (double)m_InterpolationLevel);
    double stepsize = 1/step;


    //rowpos is the global position in the grid of the interpolated points
    for(double u = 0, rowpos = 0; u <= 1; u+=stepsize, rowpos++){

        interpolatedpoints.push_back(VectorVNLType());
        interpolatedspokes.push_back(VectorVNLType());

        if( u > 0){
            VectorVNLType interpolatedspoke = GetInterpolatedSpoke(logLam0_0, logLam0_1, logLam3_0, logLam3_1,
                                 e0_0, e0_1, e3_0, e3_1,
                                 U, u, 0, stepsize, medialinterpolator, cellid, true);

            x0 = interpolatedspoke[0];
            U = interpolatedspoke[1];
        }
        interpolatedpoints[(unsigned) rowpos].push_back(x0);
        interpolatedspokes[(unsigned) rowpos].push_back(U);

        cout<<"x0i = "<<x0<<endl;
        cout<<"Ui = "<<U<<endl;

        for(double v = stepsize; v <= 1; v+=stepsize){



            /*double stepu = stepsize/10;
            double stepv = stepsize/10;

            VNLType spokeintegral(3);//integration over a straight line to new interpolated spoke
            spokeintegral.fill(0);

            double curru = u;
            for(double currv = v - stepsize; currv <= v ; currv+=stepv){
                //weighted average
                // To find rSrad, need to calculate lambdas and eigenvectors
                double logAvg0 = (1-curru)*(1-currv)*logLam0_0 + (curru)*(1-currv)*logLam3_0 + (1-curru)*(currv)*logLam1_0 + (curru)*(currv)*logLam2_0;
                double logAvg1 = (1-curru)*(1-currv)*logLam0_1 + (curru)*(1-currv)*logLam3_1 + (1-curru)*(currv)*logLam1_1 + (curru)*(currv)*logLam2_1;

                double Lam0 = 1 - exp(logAvg0);
                double Lam1 = 1 - exp(logAvg1);

                vnl_diag_matrix<double> NewL(2);
                NewL[0] = Lam0;
                NewL[1] = Lam1;

                VNLType avge0, avge1;
                double avgtheta0, avgtheta1;

                avge0 = (1-curru)*(1-currv)*e0_0 + (curru)*(1-currv)*e3_0 + (1-curru)*(currv)*e1_0 + (curru)*(currv)*e2_0;
                avgtheta0 = atan2(avge0[1], avge0[0]);

                avge1 = (1-curru)*(1-currv)*e0_1 + (curru)*(1-currv)*e3_1 + (1-curru)*(currv)*e1_1 + (curru)*(currv)*e2_1;
                avgtheta1 = atan2(avge1[1], avge1[0]);


                VNLMatrixType NewV(2,2);
                NewV[0][0] = cos(avgtheta0);
                NewV[0][1] = cos(avgtheta1);
                NewV[1][0] = sin(avgtheta0);
                NewV[1][1] = sin(avgtheta1);

                //cout<<"i= "<<NewV.get_column(0)<<" . "<<NewV.get_column(1)<<endl;


                //VNLMatrixType NewrSrad;
                //if(vnl_determinant(NewV)!=0){
                vnl_symmetric_eigensystem<double> eigensystem(VNLMatrixType(2,2));
                eigensystem.D = NewL;
                eigensystem.V = NewV;
                VNLMatrixType NewrSrad = eigensystem.recompose();
                //}else{
                //    NewrSrad = rSrad0;
                //}
                //NewrSrad = rSrad0;

                //cout<<"Srad = "<<endl<<rSrad0<<endl;
                cout<<"newSrad = "<<endl<<NewrSrad<<endl;

                VNLType dxu = medialinterpolator->GetInterpolatedDu(cellid, curru, currv);
                VNLType dxv = medialinterpolator->GetInterpolatedDv(cellid, curru, currv);

                VNLMatrixType dSdun = GetdSdu(dxu,dxv, U, NewrSrad);

                VNLType dSdu = dSdun.get_row(0);
                VNLType dSdv = dSdun.get_row(1);

                U = U*r + (stepu * dSdu) + (stepv * dSdv);
                U = U.normalize();
            }*/

            VNLType x0 = medialinterpolator->GetInterpolatedPoint(cellid, u, v);
            VNLType newSpoke = U;

            //cout<<newSpoke + x0<<endl;

            interpolatedpoints[(unsigned) rowpos].push_back(x0);
            interpolatedspokes[(unsigned) rowpos].push_back(newSpoke);
        }
    }
}


vtkSRepInterpolateMedialSpokes::VectorVNLType vtkSRepInterpolateMedialSpokes::GetInterpolatedSpoke(double logLam0_0, double logLam0_1, double logLam1_0, double logLam1_1,
                                                                                                   VNLType e0_0, VNLType e0_1, VNLType e1_0, VNLType e1_1,
                                                                                                   VNLType U, double u, double v, double stepsize,
                                                                                                   vtkSmartPointer<vtkSRepInterpolateMedialSheet> medialinterpolator, vtkIdType cellid,
                                                                                                   bool moveu){

    double step = stepsize/10;

    double current = v - stepsize;
    double limit = v;

    if(moveu){
        current = u - stepsize;
        limit = u;
    }

    for(current; current <= limit; current+=step){
        //weighted average
        // To find rSrad, need to calculate lambdas and eigenvectors
        double logAvg0 = (1-current)*logLam0_0 + (current)*logLam1_0;
        double logAvg1 = (1-current)*logLam0_1 + (current)*logLam1_1;

        double Lam0 = 1 - exp(logAvg0);
        double Lam1 = 1 - exp(logAvg1);

        vnl_diag_matrix<double> NewL(2);
        NewL[0] = Lam0;
        NewL[1] = Lam1;

        VNLType avge0, avge1;
        double avgtheta0, avgtheta1;

        avge0 = (1-current)*e0_0 + (current)*e1_0;
        avgtheta0 = atan2(avge0[1], avge0[0]);

        avge1 = (1-current)*e0_1 + (current)*e1_1;
        avgtheta1 = atan2(avge1[1], avge1[0]);


        VNLMatrixType NewV(2,2);
        NewV[0][0] = cos(avgtheta0);
        NewV[0][1] = cos(avgtheta1);
        NewV[1][0] = sin(avgtheta0);
        NewV[1][1] = sin(avgtheta1);

        //cout<<"i= "<<NewV.get_column(0)<<" . "<<NewV.get_column(1)<<endl;


        //VNLMatrixType NewrSrad;
        //if(vnl_determinant(NewV)!=0){
        vnl_symmetric_eigensystem<double> eigensystem(VNLMatrixType(2,2));
        eigensystem.D = NewL;
        eigensystem.V = NewV;
        VNLMatrixType NewrSrad = eigensystem.recompose();
        //}else{
        //    NewrSrad = rSrad0;
        //}
        //NewrSrad = rSrad0;

        //cout<<"Srad = "<<endl<<rSrad0<<endl;

        VNLType dxu;
        VNLType dxv;
        VNLType x0;

        double stepu = step;
        double stepv = step;
        if(moveu){
            dxu = medialinterpolator->GetInterpolatedDu(cellid, current, v);
            dxv = medialinterpolator->GetInterpolatedDv(cellid, current, v);
            x0 = medialinterpolator->GetInterpolatedPoint(cellid, current, v);
            stepv = 0;
        }else{
            dxu = medialinterpolator->GetInterpolatedDu(cellid, u, current);
            dxv = medialinterpolator->GetInterpolatedDv(cellid, u, current);
            x0 = medialinterpolator->GetInterpolatedPoint(cellid, u, current);
            stepu = 0;
        }

        VNLMatrixType dSdun = GetdSdu(dxu,dxv, U, NewrSrad);

        VNLType dSdu = dSdun.get_row(0);
        VNLType dSdv = dSdun.get_row(1);

        U = U + (stepu * dSdu) + (stepv * dSdv);

        cout<<"newSrad = "<<endl<<NewrSrad<<endl;
        cout<<"dxu "<<dxu<<endl;
        cout<<"dxv "<<dxv<<endl;
        cout<<"x0 "<<x0<<endl;
        cout<<"step "<<step<<endl;
        cout<<"dSdu "<<dSdu<<endl;
        cout<<"dSdv "<<dSdv<<endl;
        cout<<"Ui= "<<U<<endl;


    }

    VNLType pos = medialinterpolator->GetInterpolatedPoint(cellid, u,v);
    VectorVNLType spoke;
    spoke.push_back(pos);
    spoke.push_back(U);
    return spoke;

}

/**
*	\fn double GetLogLambda(double lam)
*	\brief evaluates if the value lam < 1 and returns a value
*	\post if the value is greater than 1, returns a large negative value
*/
double vtkSRepInterpolateMedialSpokes::GetLogLambda(double lam){
    if(lam > 1){
        //return -1*exp(lam - 1);
        return -100;
    }else{
        return log(1 - lam);
    }
}

vtkSRepInterpolateMedialSpokes::VNLMatrixType vtkSRepInterpolateMedialSpokes::GetdSdu(VNLType dxu,VNLType dxv,VNLType U,VNLMatrixType rSrad){

    VNLMatrixType pUn(2, 3);
    pUn[0][0] = dxu[0];
    pUn[0][1] = dxu[1];
    pUn[0][2] = dxu[2];

    pUn[1][0] = dxv[0];
    pUn[1][1] = dxv[1];
    pUn[1][2] = dxv[2];

    VNLMatrixType Um(1,3);
    Um[0][0] = U[0];
    Um[0][1] = U[1];
    Um[0][2] = U[2];

    VNLMatrixType rUn(2,1);
    rUn[0][0] = 0;
    rUn[1][0] = 0;


    VNLMatrixType I(3, 3);//identity matrix
    I.fill(0);
    I.fill_diagonal(1);

    VNLMatrixType Q = pUn * ((Um.transpose() * Um) - I);

    return (rSrad.transpose() * Q) + (rUn * Um);
}
