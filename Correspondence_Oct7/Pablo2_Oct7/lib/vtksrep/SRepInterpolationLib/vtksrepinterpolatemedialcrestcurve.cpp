#include "vtksrepinterpolatemedialcrestcurve.h"


#include "vtkSmartPointer.h"
#include "vtksrep.h"

#include "vtkCellArray.h"
#include "vtkLine.h"

#include "vtkInformationVector.h"
#include "vtkInformation.h"
#include "vnl/vnl_least_squares_function.h"

#include "minimizecurvaturefunction.h"

#include "vnl/algo/vnl_levenberg_marquardt.h"

#include "vtkObjectFactory.h"
vtkStandardNewMacro(vtkSRepInterpolateMedialCrestCurve);



vtkSRepInterpolateMedialCrestCurve::vtkSRepInterpolateMedialCrestCurve()
{
    m_InterpolationLevel = 0;
    m_Coeffs.clear();
}

vtkSRepInterpolateMedialCrestCurve::~vtkSRepInterpolateMedialCrestCurve()
{
    m_Coeffs.clear();
    m_CurveDerivatives.clear();
    m_CurvePoints.clear();
}


// Superclass method to update the pipeline
int vtkSRepInterpolateMedialCrestCurve::RequestData(vtkInformation* request,
                        vtkInformationVector** inputVector,
                        vtkInformationVector* outputVector){



    vtkInformation *outInfo = outputVector->GetInformationObject(0);

    vtkSRep::VectorVNLType crestpositions;
    vtkSRep::VectorIdsType crestids;


    // get the info objects
    vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
    // get the input and output
    vtkSRep *input = dynamic_cast<vtkSRep*>(vtkSRep::SafeDownCast(inInfo->Get(vtkSRep::DATA_OBJECT())));
    vtkPolyData *output = vtkPolyData::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));


    vtkSmartPointer<vtkCellArray> interpolatedcellarray = vtkSmartPointer<vtkCellArray>::New();
    vtkSmartPointer<vtkPoints> interpolatedpoints = vtkSmartPointer<vtkPoints>::New();

    vtkSRep::VectorVNLType crestderivatives;
    vtkSRep::VectorVNLType crestnormals;
    m_Coeffs.clear();

    if(m_CurvePoints.size() == 0){

        crestids = input->GetCrestMedialAtomsIds();
        crestpositions = input->GetCrestMedialAtoms(crestids);
        crestderivatives = input->GetCrestMedialAtomsDerivatives(crestpositions, false, m_CyclicCurve);
        crestnormals = input->GetMedialSheetNormals(crestids);

    }else{
        crestpositions = m_CurvePoints;
        if(m_CurveDerivatives.size() != 0){
            crestderivatives = m_CurveDerivatives;
        }else{
            if(input){
                crestderivatives = input->GetCrestMedialAtomsDerivatives(crestpositions);
            }else{
                for(unsigned i=0; i < crestpositions.size(); i++){
                    VNLType dp0(3);

                    if(i == 0){
                        dp0 = (crestpositions[i+1] - crestpositions[crestpositions.size()-1])/2.0;
                    }else if(i == crestpositions.size() - 1){
                        dp0 = (crestpositions[0] - crestpositions[i-1])/2.0;
                    }else{
                        dp0 = (crestpositions[i+1] - crestpositions[i-1])/2.0;
                    }
                    crestderivatives.push_back(dp0);
                }
            }
        }
    }


    int numpoints = pow(2.0,(double) m_InterpolationLevel);

    for(unsigned i=0; i < crestpositions.size(); i++){

        VNLType p0 = crestpositions[i];
        VNLType dp0;
        VNLType pend;
        VNLType dpend;

        if(crestnormals.size() > 0){
            //dp0 = vnl_cross_3d(crestvectors[i], crestnormals[i]);
            dp0 = vnl_cross_3d(crestnormals[i], crestderivatives[i]);
            dp0 = vnl_cross_3d(dp0, crestnormals[i]);

            if(i == crestpositions.size() - 1 && m_CyclicCurve){
                pend = crestpositions[0];

                //dpend = vnl_cross_3d(crestvectors[0], crestnormals[0]);
                dpend = vnl_cross_3d(crestnormals[0], crestderivatives[0]);
                dpend = vnl_cross_3d(dpend, crestnormals[0]);
            }else if(i == crestpositions.size() - 1 && !m_CyclicCurve){
                break;
            }else{
                pend = crestpositions[i+1];

                //dpend = vnl_cross_3d(crestvectors[i+1], crestnormals[i+1]);
                dpend = vnl_cross_3d(crestnormals[i+1], crestderivatives[i+1]);
                dpend = vnl_cross_3d(dpend, crestnormals[i+1]);
            }
        }else{
            dp0 = crestderivatives[i];

            if(i == crestpositions.size() - 1){
                if( m_CyclicCurve){
                    pend = crestpositions[0];
                    dpend = crestderivatives[0];
                }else{
                    break;
                }
            }else{
                pend = crestpositions[i+1];
                dpend = crestderivatives[i+1];
            }
        }

        VectorVNLType interpcrest = GetInterpolatedPoints(p0, pend, dp0, dpend, numpoints);

        for(unsigned j = 0; j < interpcrest.size(); j++){
            VNLType interp = interpcrest[j];
            interpolatedpoints->InsertNextPoint(interp[0],interp[1],interp[2]);
        }
    }

    for(unsigned i = 0; i < interpolatedpoints->GetNumberOfPoints() - 1; i++){
        vtkSmartPointer<vtkLine> line = vtkSmartPointer<vtkLine>::New();
        line->GetPointIds()->SetId(0, i);
        line->GetPointIds()->SetId(1, i + 1);

        interpolatedcellarray->InsertNextCell(line);
    }


    output->SetPoints(interpolatedpoints);
    output->SetLines(interpolatedcellarray);


    return 1;

}


vtkSRepInterpolateMedialCrestCurve::VectorVNLType vtkSRepInterpolateMedialCrestCurve::GetInterpolatedPoints(VNLType p0, VNLType pend, VNLType dp0, VNLType dpend, int numpoints){

    VectorVNLType interpolatedpoints;

    VectorVNLType coeffs;

    MinimizeCurvatureFunction quartic;

    for(unsigned i = 0; i < p0.size(); i++){

        quartic.SetR0(p0[i]);
        quartic.SetRend(pend[i]);
        quartic.SetdR0(dp0[i]);
        quartic.SetdRend(dpend[i]);

        vnl_levenberg_marquardt levenberg(quartic);

        vnl_vector<double> q0init(1);
        q0init[0] = -1;

        levenberg.minimize(q0init);

        coeffs.push_back(quartic.GetCoefficients(q0init[0]));

    }
    m_Coeffs.push_back(coeffs);

    for(int i = 0; i <= numpoints; i++){

        VNLType point(p0.size());
        double theta = ((double)i)/((double)numpoints);

        for(unsigned j = 0; j < p0.size(); j++){
            point[j] = quartic.EvaluateFunction(coeffs[j][0], coeffs[j][1], coeffs[j][2], p0[j], dp0[j], theta);
        }

        interpolatedpoints.push_back(point);

    }

    return interpolatedpoints;

}

vtkSRepInterpolateMedialCrestCurve::VNLType vtkSRepInterpolateMedialCrestCurve::GetInterpolatedPoint(vtkIdType idcoeffs, double t){

    if(idcoeffs < m_Coeffs.size()){
        MinimizeCurvatureFunction quartic;


        VNLType point(m_CurvePoints[idcoeffs].size());

        for(unsigned j = 0; j < point.size(); j++){
            point[j] = quartic.EvaluateFunction(m_Coeffs[idcoeffs][j][0], m_Coeffs[idcoeffs][j][1], m_Coeffs[idcoeffs][j][2], m_CurvePoints[idcoeffs][j], m_CurveDerivatives[idcoeffs][j], t);
        }

        return point;
    }



    return VNLType();
}

