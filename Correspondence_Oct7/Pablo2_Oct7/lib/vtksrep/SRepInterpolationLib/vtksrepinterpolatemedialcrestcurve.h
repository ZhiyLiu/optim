#ifndef VTKSREPINTERPOLATEMEDIALCRESTCURVE_H
#define VTKSREPINTERPOLATEMEDIALCRESTCURVE_H

#include "vtkPolyDataAlgorithm.h"
#include "vnl/vnl_vector.h"
#include "vnl/vnl_cross.h"
#include "vnl/vnl_numeric_traits.h"

#include <vector>
using namespace std;

class vtkSRepInterpolateMedialCrestCurve : public vtkPolyDataAlgorithm
{
public:
    static vtkSRepInterpolateMedialCrestCurve *New();


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
        m_CurvePoints = curvePoints;
    }
    // \brief optionally set the curvepoints derivatives to be interpolated
    void SetCurveDerivatives(VectorVNLType curveDerivatives){
        m_CurveDerivatives = curveDerivatives;
    }

    VNLType GetInterpolatedPoint(vtkIdType idcoeffs, double t);

    void SetCyclicCurve(bool cyclic){
        m_CyclicCurve = cyclic;
    }

protected:
    vtkSRepInterpolateMedialCrestCurve();

    ~vtkSRepInterpolateMedialCrestCurve();

    // Superclass method to update the pipeline
    virtual int RequestData(vtkInformation* request,
                            vtkInformationVector** inputVector,
                            vtkInformationVector* outputVector);

private:


    int m_InterpolationLevel;

    VectorVNLType m_CurvePoints;
    VectorVNLType m_CurveDerivatives;
    VectorVectorVNLType m_Coeffs;
    bool m_CyclicCurve;

    VectorVNLType GetInterpolatedPoints(VNLType p0, VNLType pend, VNLType dp0, VNLType dpend, int numpoints = 20);

};

#endif // VTKSREPINTERPOLATEMEDIALSHEET_H
