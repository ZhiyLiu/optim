#include "vtksrepinterpolatecrestspokes.h"


#include "vtkSmartPointer.h"


#include "vtkCellArray.h"
#include "vtkLine.h"

#include "vtkInformationVector.h"
#include "vtkInformation.h"
#include "vnl/vnl_least_squares_function.h"
#include "vtkinterpolatecurve.h"
#include "vtkAppendPolyData.h"

#include "minimizecurvaturefunction.h"

#include "vnl/algo/vnl_levenberg_marquardt.h"

#include "vtkObjectFactory.h"
vtkStandardNewMacro(vtkSRepInterpolateCrestSpokes);



vtkSRepInterpolateCrestSpokes::vtkSRepInterpolateCrestSpokes()
{
    m_InterpolationLevel = 0;
}

vtkSRepInterpolateCrestSpokes::~vtkSRepInterpolateCrestSpokes()
{
}


// Superclass method to update the pipeline
int vtkSRepInterpolateCrestSpokes::RequestData(vtkInformation* request,
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


    vtkSRep::VectorIdsType crestids = input->GetCrestMedialAtomsIds();    

    double epsilon = 0.5;


    VectorSRepVectorVNLType epsilonpositionsvector;//this contains the inflated position at each spoke
    VectorSRepVectorVNLType epsilonpositionsderivativesvector;//contains the derivatives of the inflated positions
    VectorSRepVectorVNLType epsilonradiusvector;

    vtkSRep::VectorVNLType crestnormals = input->GetMedialSheetNormals(crestids);


    //Inflating the positions for every point by epsilon in the direction of the spokes vector to create the cylinder
    //around the crest curve
    //unsigned maxspokesnumber = 0;

    vector<vtkSRep::SPOKES_TYPE > spokestype;
    spokestype.push_back(vtkSRep::BOTTOM_SPOKE);
    spokestype.push_back(vtkSRep::CREST_SPOKE);
    spokestype.push_back(vtkSRep::TOP_SPOKE);


    for(unsigned i = 0; i < spokestype.size(); i++){

        vtkSRep::VectorVNLType epsilonpositions;
        vtkSRep::VectorVNLType epsilonradius;

        for(unsigned j = 0; j < crestids.size(); j++){

            vtkIdType crestatomid = crestids[j];

            vtkSRep::VectorVNLType spokes = input->GetSpokes(crestatomid);
            vtkSRep::VectorDoubleType radius = input->GetSpokesRadius(crestatomid);

            vtkSRep::VNLType eradius(1);
            eradius[0] = radius[spokestype[i]]*epsilon;

            double *temp = input->GetPoint(crestatomid);
            vtkSRep::VNLType atompos(3);
            atompos[0] = temp[0];
            atompos[1] = temp[1];
            atompos[2] = temp[2];
            vtkSRep::VNLType newpos = atompos + spokes[spokestype[i]]*eradius[0];

            epsilonpositions.push_back(newpos);
            epsilonradius.push_back(eradius);

        }

        epsilonpositionsvector.push_back(epsilonpositions);
        epsilonradiusvector.push_back(epsilonradius);

        epsilonpositions.clear();
    }


    for(unsigned i = 0; i < spokestype.size(); i++){

        vtkSRep::VectorVNLType crestderivatives = input->GetCrestMedialAtomsDerivatives(epsilonpositionsvector[i]);
        vtkSRep::VectorVNLType crestderivativesnormals;

        for(unsigned j = 0; j < crestderivatives.size(); j++){//recompute direction using the normal to the sheet

            VNLType vect = vnl_cross_3d(crestnormals[j], crestderivatives[j]);
            vect = vnl_cross_3d(vect, crestnormals[j]);

            crestderivativesnormals.push_back(vect);
        }

        epsilonpositionsderivativesvector.push_back(crestderivativesnormals);//this method works to calculate the derivatives

        crestderivativesnormals.clear();
    }

    //up to here epsilonpositionsvector and epsilonpositionsderivativesvector are organized by spokes i.e
    //epsilonpositionsvector[n] corresponds to all spokes_n in the crest, this would be simple to manipulate afterwards
    //to access a specific spoke epsilonpositionsvector[n][i] would give the vnl_vector of the spoke. the index [i] also
    //corresponds to crestids[i]


    //Start calculation of kappa and betta


    VectorSRepVectorVNLType kappabettavector = GetKappaBetta(crestids, epsilonpositionsderivativesvector, spokestype);


    vtkSRep::VectorVNLType crestpositions = input->GetCrestMedialAtoms(crestids);
    vtkSRep::VectorVNLType crestderivatives = input->GetCrestMedialAtomsDerivatives(crestpositions);

    for(unsigned j = 0; j < crestderivatives.size(); j++){//recompute direction using the normal to the sheet

        VNLType vect = vnl_cross_3d(crestnormals[j], crestderivatives[j]);
        vect = vnl_cross_3d(vect, crestnormals[j]);

        crestderivatives[j] = vect;
    }

    vtkSmartPointer< vtkInterpolateCurve > crestcurveinterpolator = vtkSmartPointer< vtkInterpolateCurve >::New();
    crestcurveinterpolator->SetInterpolationLevel(m_InterpolationLevel);
    crestcurveinterpolator->SetCurvePoints(crestpositions);
    crestcurveinterpolator->SetCurveDerivatives(crestderivatives);
    crestcurveinterpolator->Update();


    //start of interpolation
    vtkSmartPointer< vtkAppendPolyData > appendpoly = vtkSmartPointer< vtkAppendPolyData >::New();

    for(unsigned i = 1; i < kappabettavector.size() - 1; i++){

        vtkSmartPointer< vtkInterpolateCurve > curveinterpolator = vtkSmartPointer< vtkInterpolateCurve >::New();
        curveinterpolator->SetInterpolationLevel(m_InterpolationLevel);
        curveinterpolator->SetCurvePoints(epsilonpositionsvector[i]);
        curveinterpolator->SetCurveDerivatives(epsilonpositionsderivativesvector[i]);
        curveinterpolator->Update();

        appendpoly->AddInput(curveinterpolator->GetOutput());


        //for(unsigned j = 0; j < kappabettavector[i].size() - 1; j++){
        for(unsigned j = 0; j < 1; j++){

            VNLType kappabetta0 = kappabettavector[i][j];
            VNLType kappabetta1 = kappabettavector[i][j+1];
            double epsilonr = epsilonradiusvector[i][j][0];

            double kappa_rel0 = log(1 - kappabetta0[0]);
            double kappa_rel1 = log(1 - kappabetta1[0]);
            double betta0 = log(1 - kappabetta0[1]);
            double betta1 = log(1 - kappabetta1[1]);

            double step = pow((double)2, (double)m_InterpolationLevel);
            double stepsize = 1/step;
            double stept = stepsize/10;

            VNLType U = input->GetSpoke(crestids[j], spokestype[i], true);
            //VNLType x = crestcurveinterpolator->EvaluateFunction(j, 0);
            //U = x + U;

            for(double u = stepsize; u <= 1; u+=stepsize){

                for(double t = u-stepsize; t <= u; t+=stept){

                    double kappa_rel = (1-t)*kappa_rel0 + t*kappa_rel1;
                    kappa_rel = 1 - exp(kappa_rel);

                    double betta = (1-t)*betta0 + t*betta1;
                    betta = 1 - exp(betta);

                    vtkSRep::VNLMatrixType Srad(2,2);
                    double a1 = -1/epsilonr;
                    double a2 = kappa_rel/(1-epsilonr*kappa_rel);
                    Srad[0][0] = a1;
                    Srad[0][1] = betta * (a1 - a2);                    
                    Srad[1][0] = 0;
                    Srad[1][1] = a2;

                    VNLType dxu = curveinterpolator->EvaluateDerivativeFunction(j, t);
                    VNLType dxv(3);
                    dxv.fill(0);

                    //Srad = epsilonr*Srad;
                    VNLMatrixType dSdu = GetdSdu(dxu, dxv, epsilonr*U, Srad);

                    VNLType dUu = dSdu.get_row(0);
                    VNLType dUv = dSdu.get_row(1);

                    //cout<<dUu<<" "<<dUv<<endl;

                    U = U + dUu*stept + dUv*stept;
                    //U = U + dUu*stept + dUv*stept;
                    //U = U + dxu*stept;
                    //epsilonr = U.magnitude();
                    //U.normalize();

                }

                U = U*epsilonr/epsilon;


                VNLType x = crestcurveinterpolator->EvaluateFunction(j, u);
                vtkIdType id0 = interpolatedpoints->InsertNextPoint(x[0], x[1], x[2]);
                vtkIdType id1 = interpolatedpoints->InsertNextPoint(x[0] + U[0], x[1] + U[1], x[2] + U[2]);
                //vtkIdType id1 = interpolatedpoints->InsertNextPoint(U[0], U[1], U[2]);


                vtkLine* line = vtkLine::New();
                line->GetPointIds()->SetId(0, id0);
                line->GetPointIds()->SetId(1, id1);
                interpolatedcellarray->InsertNextCell(line);
            }
        }
    }





    //output->SetPoints(interpolatedpoints);
    //output->SetLines(interpolatedcellarray);


    vtkSmartPointer<vtkPolyData> outp = vtkSmartPointer<vtkPolyData>::New();
    outp->SetPoints(interpolatedpoints);
    outp->SetLines(interpolatedcellarray);

    appendpoly->AddInput(outp);
    appendpoly->Update();
    output->DeepCopy(appendpoly->GetOutput());



    return 1;

}


vtkSRepInterpolateCrestSpokes::VectorSRepVectorVNLType vtkSRepInterpolateCrestSpokes::GetKappaBetta(vtkSRep::VectorIdsType crestids, VectorSRepVectorVNLType epsilonpositionsderivativesvector, vector<vtkSRep::SPOKES_TYPE> spokestype){

    VectorSRepVectorVNLType kappabettavector;
    vtkSRep* input = m_Input;

    for(unsigned i = 0; i < spokestype.size(); i++){

        vtkSRep::VectorVNLType uvects = input->GetSpokes(crestids, spokestype[i], true);//get unit spokes
        vtkSRep::VectorVNLType uderivatives = input->GetCrestUDerivatives(crestids, false, spokestype[i]);//derivatives of U
        vtkSRep::VectorVNLType gammaderivatives = epsilonpositionsderivativesvector[i];//derivatives of the dilated points
        vtkSRep::VectorVNLType kappabetta;        

        for(unsigned j = 0; j < uvects.size(); j++){

            VNLType e_1 = input->GetMedialSheetNormal(crestids[j]);
            VNLType e_3 = input->GetMedialSheetVecte3(crestids[j]);//gammaderivatives[j].normalize();
            VNLType e_2 = vnl_cross_3d(e_1, e_3);
            e_2 = e_2.normalize();

            //cout<<dot_product(e_1, e_3)<<" "<<dot_product(e_1, e_2)<<endl;

            double kappa_rel = dot_product(uderivatives[j], e_3);
            kappa_rel = kappa_rel/dot_product(gammaderivatives[j], e_3);

            double dot = dot_product(e_1, uvects[j])/(e_1.magnitude()*uvects[j].magnitude());
            double theta = acos(dot);
            VNLType u_theta = -sin(theta)*e_1 + cos(theta)*e_2;//orthogonal complement
            u_theta = u_theta.normalize();

            VNLType v_theta = kappa_rel*(gammaderivatives[j] - (dot_product(gammaderivatives[j], uvects[j]) * uvects[j]) ) - uderivatives[j];
            double betta_2 = dot_product(v_theta, u_theta);
            //double betta_2 = dot_product(uderivatives[j] + kappa_rel*gammaderivatives[j], u_theta);

            //cout<<"kappa_rel = "<<kappa_rel<<" betta_2 = "<<betta_2<<endl;

            VNLType kb(2);
            kb[0] = kappa_rel;
            kb[1] = betta_2;

            kappabetta.push_back(kb);
        }
        kappabettavector.push_back(kappabetta);
    }
    return kappabettavector;
}


vtkSRepInterpolateCrestSpokes::VNLMatrixType vtkSRepInterpolateCrestSpokes::GetdSdu(VNLType dxu,VNLType dxv,VNLType U,VNLMatrixType Srad){

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

    /*VNLMatrixType Au = Srad.transpose()*pUn*Um.transpose();

    return  Au*Um - Srad.transpose()*pUn;*/

    //VNLMatrixType ruv = -pUn*Um.transpose();
    VNLMatrixType ruv(2,1);
    ruv.fill(0);

    VNLMatrixType I(3,3);
    I.fill(0);
    I.fill_diagonal(1);

    return Srad.transpose() * pUn *(Um.transpose() * Um - I) + ruv*Um;



}

