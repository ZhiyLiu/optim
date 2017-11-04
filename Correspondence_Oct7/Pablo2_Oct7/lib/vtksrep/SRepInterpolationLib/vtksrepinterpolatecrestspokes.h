#ifndef VTKSREPINTERPOLATECRESTSPOKES_H
#define VTKSREPINTERPOLATECRESTSPOKES_H

#include "vtkPolyDataAlgorithm.h"
#include "vnl/vnl_vector.h"
#include "vnl/vnl_cross.h"
#include "vnl/vnl_numeric_traits.h"

#include "vtksrep.h"

using namespace std;

class VTK_FILTERING_EXPORT vtkSRepInterpolateCrestSpokes : public vtkPolyDataAlgorithm
{
public:
    static vtkSRepInterpolateCrestSpokes *New();


    typedef vnl_vector<double> VNLType;
    typedef vector< VNLType > VectorVNLType;
    typedef vector< vtkSRep::VectorVNLType > VectorSRepVectorVNLType;
    typedef vtkSRep::VNLMatrixType VNLMatrixType;



    void SetInterpolationLevel(int level){
        m_InterpolationLevel = level;
    }
    int GetInterpolationLevel(){
        return m_InterpolationLevel;
    }

protected:
    vtkSRepInterpolateCrestSpokes();

    ~vtkSRepInterpolateCrestSpokes();

    // Superclass method to update the pipeline
    virtual int RequestData(vtkInformation* request,
                            vtkInformationVector** inputVector,
                            vtkInformationVector* outputVector);

private:


    int m_InterpolationLevel;
    vtkSRep *m_Input;

    VectorSRepVectorVNLType GetKappaBetta(vtkSRep::VectorIdsType crestids, VectorSRepVectorVNLType epsilonpositionsderivativesvector, vector< vtkSRep::SPOKES_TYPE > spokestype);
    VNLMatrixType GetdSdu(VNLType dxu,VNLType dxv,VNLType U,VNLMatrixType Srad);

};

#endif // VTKSREPINTERPOLATEMEDIALSHEET_H
