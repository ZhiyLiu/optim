#ifndef VTKSREPINTERPOLATEMEDIALSHEET_H
#define VTKSREPINTERPOLATEMEDIALSHEET_H

#include "vtkPolyDataAlgorithm.h"
#include "vnl/vnl_vector.h"
#include "vnl/vnl_cross.h"
#include "vnl/vnl_numeric_traits.h"

#include "vtksrep.h"

#include <vector>
using namespace std;

class vtkSRepInterpolateMedialSheet : public vtkPolyDataAlgorithm
{
public:
    static vtkSRepInterpolateMedialSheet *New();


    typedef vnl_vector<double> VNLType;
    typedef vector< VNLType > VectorVNLType;

    void GetHermiteMatrix(VNLType p0, VNLType p1, VNLType p2, VNLType p3, VectorVNLType dp0, VectorVNLType dp1, VectorVNLType dp2, VectorVNLType dp3, VNLType n0, VNLType n1, VNLType n2, VNLType n3, VNLType **H);

    /*
    * \fn void GetInterpolatedPoint(double u, double v, Vector3D H[4][4], Vector3D p)
    * \brief Computes the point using the cubic Hermite weight functions and the two scalar parameters
    * \param double u scalar for the weight functions
    * \param double v second scalar for the weight functions
    * \param cellid the id of the quad, see srep->GetNumberofCells()
    * \pre u and v should be between [0, 1]
    */
    VNLType GetInterpolatedPoint(vtkIdType cellid, double u, double v);

    /*
    * \fn void GetInterpolatedPoint(double u, double v)
    * \brief Computes the point using the cubic Hermite weight functions and the two scalar parameters
    * \brief The cell id is calculated from u, v
    * \param double u scalar for the weight functions
    * \param double v second scalar for the weight functions
    * \pre u and v should be between [0, 1]
    */
    //VNLType GetInterpolatedPoint(double u, double v);

    /* \brief Using Hermite interpolation calculates the derivative in the u
              direction
    */
    VNLType GetInterpolatedDu(vtkIdType cellid, double u, double v);

    /* \brief Using Hermite interpolation calculates the derivative in the v
              direction
    */
    VNLType GetInterpolatedDv(vtkIdType cellid, double u, double v);

    void SetInterpolationLevel(int level){
        m_InterpolationLevel = level;
    }
    int GetInterpolationLevel(){
        return m_InterpolationLevel;
    }

    void SetAtomId(vtkIdType atomid){
        m_AtomId = atomid;
    }

    void SetGamma_u(double gammau){
        m_Gamma_u = gammau;
    }

    void SetGamma_v(double gammav){
        m_Gamma_v = gammav;
    }


protected:
    vtkSRepInterpolateMedialSheet();

    ~vtkSRepInterpolateMedialSheet();

    // Superclass method to update the pipeline
    virtual int RequestData(vtkInformation* request,
                            vtkInformationVector** inputVector,
                            vtkInformationVector* outputVector);

private:


    int m_InterpolationLevel;
    vector< VNLType** > m_HermitMatrices;
    vtkIdType m_AtomId;
    double m_Gamma_u;
    double m_Gamma_v;
    vtkSRep* m_Input;


    /*
    * \fn void GetInterpolatedPoint(double u, double v, Vector3D H[4][4], Vector3D p)
    * \brief Computes the point using the cubic Hermite weight functions and the two scalar parameters
    * \param double u scalar for the weight functions
    * \param double v second scalar for the weight functions
    * \param Vector3D H[4][4] cubic hermite matrix
    * \param Vector3D p result of the interpolation
    * \pre u and v should be between [0, 1]
    */
    VNLType GetInterpolatedPoint(double u, double v, VNLType** H);

    /* \brief Using Hermite interpolation calculates the derivative in the u
              direction
    */
    VNLType GetInterpolatedDu(double u, double v, VNLType** H);

    /* \brief Using Hermite interpolation calculates the derivative in the v
              direction
    */
    VNLType GetInterpolatedDv(double u, double v, VNLType** H);


};

#endif // VTKSREPINTERPOLATEMEDIALSHEET_H
