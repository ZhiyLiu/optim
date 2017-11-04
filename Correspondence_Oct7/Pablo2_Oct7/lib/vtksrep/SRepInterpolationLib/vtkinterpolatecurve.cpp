#include "vtkinterpolatecurve.h"


#include "vtkSmartPointer.h"


#include "vtkCellArray.h"
#include "vtkLine.h"

#include "vtkInformationVector.h"
#include "vtkInformation.h"
#include "vnl/vnl_least_squares_function.h"


#include "vnl/algo/vnl_levenberg_marquardt.h"

#include "vtkObjectFactory.h"
vtkStandardNewMacro(vtkInterpolateCurve);



vtkInterpolateCurve::vtkInterpolateCurve()
{
    m_InterpolationLevel = 0;

    m_Cyclic = false;
    m_CurrCoef = 0;

    this->SetNumberOfInputPorts(0);
    this->SetNumberOfOutputPorts(1);
}

vtkInterpolateCurve::~vtkInterpolateCurve()
{
}


// Superclass method to update the pipeline
int vtkInterpolateCurve::RequestData(vtkInformation* request,
                        vtkInformationVector** inputVector,
                        vtkInformationVector* outputVector){



    vtkInformation *outInfo = outputVector->GetInformationObject(0);

    VectorVNLType curvepositions = m_CurvePoints;
    VectorVNLType curvederivatives = m_CurveDerivatives;

    if(curvederivatives.size() != curvepositions.size()){
        cout<<"curve derivatives not set"<<endl;
        curvederivatives.clear();
        for(unsigned i=0; i < curvepositions.size(); i++){

            VNLType dp0;

            if(i == 0){
                if(m_Cyclic){
                    dp0 = (curvepositions[i+1] - curvepositions[curvepositions.size()-1])/2.0;
                }else{
                    dp0 = curvepositions[i+1] - curvepositions[i];
                }
            }else if(i == curvepositions.size() - 1){
                if(m_Cyclic){
                    dp0 = (curvepositions[0] - curvepositions[i-1])/2.0;
                }else{
                    dp0 = (curvepositions[i] - curvepositions[i-1]);
                }
            }else{
                dp0 = (curvepositions[i+1] - curvepositions[i-1])/2.0;
            }

            curvederivatives.push_back(dp0);
        }
    }


    vtkPolyData *output = vtkPolyData::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));


    vtkSmartPointer<vtkCellArray> interpolatedcellarray = vtkSmartPointer<vtkCellArray>::New();
    vtkSmartPointer<vtkPoints> interpolatedpoints = vtkSmartPointer<vtkPoints>::New();

    int numpoints = pow(2.0,(double) m_InterpolationLevel);

    for(int i=0; i < curvepositions.size(); i++){

        VNLType p0 = curvepositions[i];
        VNLType dp0;
        VNLType pend;
        VNLType dpend;


        dp0 = curvederivatives[i];

        if(i == ((int)curvepositions.size()) - 1){
            if(m_Cyclic){
                pend = curvepositions[0];
                dpend = curvederivatives[0];
            }else{                
                break;
            }
        }else{
            pend = curvepositions[i+1];
            dpend = curvederivatives[i+1];
        }


        VectorVNLType interpcrest = GetInterpolatedPoints(p0, pend, dp0, dpend, numpoints);

        for(unsigned j = 0; j < interpcrest.size(); j++){
            VNLType interp = interpcrest[j];
            double temp[3] = {0,0,0};
            for(unsigned j = 0; j < interp.size() && j < 3; j++){
                temp[j] = interp.get(j);
            }
            interpolatedpoints->InsertNextPoint(temp[0],temp[1],temp[2]);
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


vtkInterpolateCurve::VectorVNLType vtkInterpolateCurve::GetInterpolatedPoints(VNLType p0, VNLType pend, VNLType dp0, VNLType dpend, int numpoints){


    VectorVNLType interpolatedpoints;

    VectorVNLType coeffs;    

    for(unsigned i = 0; i < p0.size(); i++){

        m_quartic.SetR0(p0[i]);
        m_quartic.SetRend(pend[i]);
        m_quartic.SetdR0(dp0[i]);
        m_quartic.SetdRend(dpend[i]);

        vnl_levenberg_marquardt levenberg(m_quartic);

        vnl_vector<double> q0init(1);
        q0init[0] = -1;

        levenberg.minimize(q0init);

        coeffs.push_back(m_quartic.GetCoefficients(q0init[0]));

    }

    for(int i = 0; i <= numpoints; i++){

        VNLType point(p0.size());
        double theta = ((double)i)/((double)numpoints);

        for(unsigned j = 0; j < p0.size(); j++){
            point[j] = m_quartic.EvaluateFunction(coeffs[j][0], coeffs[j][1], coeffs[j][2], p0[j], dp0[j], theta);
        }

        interpolatedpoints.push_back(point);

    }

    m_Coefficients.push_back(coeffs);
    m_P0.push_back(p0);
    m_dP0.push_back(dp0);



    return interpolatedpoints;
}

/*
 *  \brief Evaluates the function using the coefficients at position id
 *         The function is evaluated using the number of m_CurvePoints
 *         Therefore there are (m_CurvePoints.size() - 1) coefficients;
 *  \param vtkIdType id the position of the coefficients to use
 *  \param double t, theta parameter between [0, 1]
*/
vtkInterpolateCurve::VNLType vtkInterpolateCurve::EvaluateFunction(vtkIdType idcoeffs, double t){

    VNLType p0 = m_P0[idcoeffs];
    VNLType dp0 = m_dP0[idcoeffs];
    VectorVNLType coeffs = m_Coefficients[idcoeffs];

    VNLType point(p0.size());
    double theta = t;

    for(unsigned j = 0; j < p0.size(); j++){
        point[j] = m_quartic.EvaluateFunction(coeffs[j][0], coeffs[j][1], coeffs[j][2], p0[j], dp0[j], theta);
    }

    return point;
}

vtkInterpolateCurve::VNLType vtkInterpolateCurve::EvaluateDerivativeFunction(vtkIdType idcoeffs, double t){

    VNLType dp0 = m_dP0[idcoeffs];
    VectorVNLType coeffs = m_Coefficients[idcoeffs];

    VNLType dpoint(dp0.size());
    double theta = t;

    for(unsigned j = 0; j < dpoint.size(); j++){
        dpoint[j] = m_quartic.EvaluateDerivativeFunction(coeffs[j][0], coeffs[j][1], coeffs[j][2], dp0[j], theta);
    }

    return dpoint;
}


/*
 *  \brief Evaluates the function using the value of t, the whole curve is parametrized between
           0-1
           Warning!!!! use only when points are evenly spaced, the function calculates right
           the coefficients to use.
 *  \param double t, theta parameter between [0, 1]
*/
vtkInterpolateCurve::VNLType vtkInterpolateCurve::EvaluateFunction(double t){

    double tc = t*m_Coefficients.size();
    double tcfl = floor(tc);
    vtkIdType coef = tcfl;
    double tcoef = tc - tcfl;

    if(coef < m_Coefficients.size()){
        m_CurrCoef = coef;
    }else{
        m_CurrCoef = m_Coefficients.size()-1;
        tcoef = 1;
    }

    return EvaluateFunction(m_CurrCoef, tcoef);
}
