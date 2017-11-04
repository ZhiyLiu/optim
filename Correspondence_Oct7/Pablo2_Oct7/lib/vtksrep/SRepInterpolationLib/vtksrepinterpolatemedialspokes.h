#ifndef VTKSREPINTERPOLATEMEDIALSPOKES_H
#define VTKSREPINTERPOLATEMEDIALSPOKES_H

#include <vtkSmartPointer.h>

#include "vnl/vnl_vector.h"

#include "vtksrep.h"

#include "vtkPolyDataAlgorithm.h"
#include "vtksrepinterpolatemedialsheet.h"

#include "iostream"
using namespace std;


class vtkSRepInterpolateMedialSpokes : public vtkPolyDataAlgorithm
{
public:
    static vtkSRepInterpolateMedialSpokes *New();


    typedef vnl_vector<double> VNLType;
    typedef vector< VNLType > VectorVNLType;
    typedef vector< VectorVNLType > VectorVectorVNLType;

    typedef vtkSRep::VNLMatrixType VNLMatrixType;
    typedef vector< VNLMatrixType > VectorVNLMatrixType;

    typedef vtkSRep::VectorIdsType VectorIdsType;
    typedef vector< VectorIdsType > VectorVectorIdsType;


    void SetInterpolationLevel(int level){
        m_InterpolationLevel = level;
    }
    int GetInterpolationLevel(){
        return m_InterpolationLevel;
    }

protected:
    vtkSRepInterpolateMedialSpokes();
    ~vtkSRepInterpolateMedialSpokes();

    // Superclass method to update the pipeline
    virtual int RequestData(vtkInformation* request,
                            vtkInformationVector** inputVector,
                            vtkInformationVector* outputVector);

private:
    int m_InterpolationLevel;
    vtkSRep* m_Input;

    void InterpolateSpokes(vtkIdType cellid, vtkIdType pointid, VectorVNLMatrixType srads, vtkIdType side, vtkSmartPointer<vtkSRepInterpolateMedialSheet> medialinterpolator, VectorVectorVNLType& interpolatedpoints, VectorVectorVNLType& interpolatedspokes);

    VectorVNLType GetInterpolatedSpoke(double logLam0_0, double logLam0_1, double logLam1_0, double logLam1_1,
                           VNLType e0_0, VNLType e0_1, VNLType e1_0, VNLType e1_1,
                           VNLType U, double u, double v, double stepsize,
                           vtkSmartPointer<vtkSRepInterpolateMedialSheet> medialinterpolator, vtkIdType cellid,
                           bool moveu);

    /**
    *	\fn double GetLogLambda(double lam)
    *	\brief evaluates if the value lam < 1 and returns a value
    *	\post if the value is greater than 1, returns a large negative value
    */
    double GetLogLambda(double lam);

    double alpha(VNLType e);

    VNLMatrixType GetdSdu(VNLType dxu,VNLType dxv,VNLType U,VNLMatrixType rSrad);
};

#endif // VTKSREPINTERPOLATEMEDIALSPOKES_H
