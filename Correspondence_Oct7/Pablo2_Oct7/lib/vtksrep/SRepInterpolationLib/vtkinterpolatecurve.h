#ifndef vtkInterpolateCurve_H
#define vtkInterpolateCurve_H

#include "vtkPolyDataAlgorithm.h"
#include "vnl/vnl_vector.h"
#include "vnl/vnl_cross.h"
#include "vnl/vnl_numeric_traits.h"

#include "minimizecurvaturefunction.h"

#include <vector>

using namespace std;

class vtkInterpolateCurve : public vtkPolyDataAlgorithm
{
public:
    static vtkInterpolateCurve *New();


    typedef vnl_vector<double> VNLType;
    typedef vector< VNLType > VectorVNLType;
    typedef vector< VectorVNLType > VectorVectorVNLType;


    void SetInterpolationLevel(int level){
        m_InterpolationLevel = level;
    }
    int GetInterpolationLevel(){
        return m_InterpolationLevel;
    }

    // \brief optionally set the curvepoints to be interpolated
    void SetCurvePoints(VectorVNLType curvePoints){
        this->Modified();
        m_CurvePoints = curvePoints;
    }
    // \brief optionally set the curvepoints derivatives to be interpolated
    void SetCurveDerivatives(VectorVNLType curveDerivatives){
        this->Modified();
        m_CurveDerivatives = curveDerivatives;
    }    


    void SetCyclic(bool cyclic){
        m_Cyclic = cyclic;
    }

    /*
     *  \brief Evaluates the function using the coefficients at position id
     *         The function is evaluated using the number of m_CurvePoints
     *         Therefore there are (m_CurvePoints.size() - 1) coefficients;
     *  \param vtkIdType id the position of the coefficients to use
     *  \param double t, theta parameter between [0, 1]
    */
    VNLType EvaluateFunction(vtkIdType idcoeffs, double t);
    VNLType EvaluateDerivativeFunction(vtkIdType idcoeffs, double t);

    /*
     *  \brief Evaluates the function using the coefficients at position id
     *         The function is evaluated using the number of m_CurvePoints
     *         Therefore there are (m_CurvePoints.size() - 1) coefficients;
     *  \param vtkIdType id the position of the coefficients to use
     *  \param double t, theta parameter between [0, 1]
    */
    VNLType EvaluateFunction(double t);

    vtkIdType GetCurrCoef(){
        return m_CurrCoef;
    }

    int GetNumberOfPoints(){
        return m_Coefficients.size();
    }

protected:
    vtkInterpolateCurve();

    ~vtkInterpolateCurve();

    // Superclass method to update the pipeline
    virtual int RequestData(vtkInformation* request,
                            vtkInformationVector** inputVector,
                            vtkInformationVector* outputVector);

private:


    int m_InterpolationLevel;

    VectorVNLType m_CurvePoints;
    VectorVNLType m_CurveDerivatives;
    VectorVectorVNLType m_Coefficients;
    VectorVNLType m_P0;
    VectorVNLType m_dP0;
    bool m_Cyclic;

    MinimizeCurvatureFunction m_quartic;

    vtkIdType m_CurrCoef;

    VectorVNLType GetInterpolatedPoints(VNLType p0, VNLType pend, VNLType dp0, VNLType dpend, int numpoints = 20);

};

#endif
