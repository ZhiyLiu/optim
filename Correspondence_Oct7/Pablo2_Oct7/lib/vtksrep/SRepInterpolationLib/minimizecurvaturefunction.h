#ifndef MINIMIZECURVATUREFUNCTION_H
#define MINIMIZECURVATUREFUNCTION_H

#include "vnl/vnl_least_squares_function.h"

class MinimizeCurvatureFunction : public vnl_least_squares_function
{
public:
    MinimizeCurvatureFunction();
    MinimizeCurvatureFunction(double rend, double r0,double drend,double dr0);
    ~MinimizeCurvatureFunction();

    virtual void f(vnl_vector< double > const &x, vnl_vector< double > &fx);
    //virtual void gradf(vnl_vector< double > const &x, vnl_matrix< double > &jacobian);

    void SetRend(double rend){
        m_Rend = rend;
    }
    void SetR0(double r0){
        m_R0 = r0;
    }
    void SetdRend(double drend){
        m_dRend = drend;
    }
    void SetdR0(double dr0){
        m_dR0 = dr0;
    }

    vnl_vector< double > GetCoefficients(double q0, double thetamax = 1);

    double EvaluateFunction(double q0, double q1, double q2, double r0, double dr0, double t);
    double EvaluateDerivativeFunction(double q0, double q1, double q2, double dr0, double t);
private:
    double m_Rend;
    double m_R0;
    double m_dRend;
    double m_dR0;
};

#endif // MINIMIZECURVATUREFUNCTION_H
