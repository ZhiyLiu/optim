#include "minimizecurvaturefunction.h"

#include "iostream"

using namespace std;

MinimizeCurvatureFunction::MinimizeCurvatureFunction()
    : vnl_least_squares_function(1,21,no_gradient)
{
    m_Rend = 0;
    m_R0 = 0;;
    m_dRend = 0;;
    m_dR0 = 0;
}

MinimizeCurvatureFunction::MinimizeCurvatureFunction(double rend, double r0, double drend, double dr0)
 : vnl_least_squares_function(1,21,no_gradient){
    SetRend(rend);
    SetR0(r0);
    SetdRend(drend);
    SetdR0(dr0);
}

MinimizeCurvatureFunction::~MinimizeCurvatureFunction()
{
}

void MinimizeCurvatureFunction::f(vnl_vector< double > const &x, vnl_vector< double > &fx){

    double thetamax = 1;//this theta can be another maximum
    double q0 = x[0];


    double q2 = 6.0/pow(thetamax,4)*(2*m_dRend*thetamax + q0*pow(thetamax,2) - 6*m_Rend + 6*m_R0 + 4*m_dR0*thetamax);
    double q1 = -6.0/pow(thetamax,3)*(-m_Rend+m_R0+m_dR0*thetamax+(q0*pow(thetamax,2))/2+(q2*pow(thetamax,4))/12);

    double theta = 0;

    for(unsigned i = 0; i <= 20; i++){
        theta = ((double)i)/20.0;
        fx[i] = q0+q1*theta + q2*pow(theta,2);
    }

}


vnl_vector< double > MinimizeCurvatureFunction::GetCoefficients(double q0, double thetamax){
    vnl_vector< double > coeffs(3);

    coeffs[0] = q0;
    double q2 = 6.0/pow(thetamax,4)*(2*m_dRend*thetamax + q0*pow(thetamax,2) - 6*m_Rend + 6*m_R0 + 4*m_dR0*thetamax);
    coeffs[2] = q2;
    coeffs[1] = -6.0/pow(thetamax,3)*(-m_Rend+m_R0+m_dR0*thetamax+(q0*pow(thetamax,2))/2+(q2*pow(thetamax,4))/12);

    return coeffs;

}

double MinimizeCurvatureFunction::EvaluateFunction(double q0, double q1, double q2, double r0, double dr0, double t){
    return r0 + dr0*t + q0*pow(t, 2.0)/2.0 + q1*pow(t,3.0)/6.0 + q2*pow(t,4.0)/12.0;
}

double MinimizeCurvatureFunction::EvaluateDerivativeFunction(double q0, double q1, double q2, double dr0, double t){
    return dr0 + q0*t + q1*pow(t,2.0)/2.0 + q2*pow(t,3)/3.0;
}
