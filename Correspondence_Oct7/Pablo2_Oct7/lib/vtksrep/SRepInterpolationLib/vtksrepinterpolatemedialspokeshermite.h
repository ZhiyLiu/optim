#ifndef vtkSRepInterpolateMedialSpokesHermite_H
#define vtkSRepInterpolateMedialSpokesHermite_H

#include <vtkSmartPointer.h>

#include "vnl/vnl_vector.h"

#include "vtksrep.h"

#include "vtkPolyDataAlgorithm.h"

#include "iostream"
using namespace std;


class vtkSRepInterpolateMedialSpokesHermite : public vtkPolyDataAlgorithm
{
public:
    static vtkSRepInterpolateMedialSpokesHermite *New();


    typedef vnl_vector<double> VNLType;
    typedef vector< VNLType > VectorVNLType;
    typedef vector< VectorVNLType > VectorVectorVNLType;

    typedef vtkSRep::VNLMatrixType VNLMatrixType;
    typedef vector< VNLMatrixType > VectorVNLMatrixType;

    typedef vtkSRep::VectorIdsType VectorIdsType;
    typedef vector< VectorIdsType > VectorVectorIdsType;


    vtkSmartPointer<vtkSRep> GetSRepOutput(){
        return m_SRepOutput;
    }

    void SetInterpolationLevel(int level){
        m_InterpolationLevel = level;
    }
    int GetInterpolationLevel(){
        return m_InterpolationLevel;
    }

    /*
    * \fn void GetInterpolatedPoint(vtkIdType cellid, double u, double v)
    * \brief Computes the point using the cubic Hermite weight functions and the two scalar parameters
    * \param vtkIdType cellid, cellid in the srep
    * \param double u scalar for the weight functions
    * \param double v second scalar for the weight functions    n
    * \pre u and v should be between [0, 1]
    */
    VNLType GetInterpolatedSpoke(vtkIdType cellid, unsigned side, double u, double v);

    /*
    * \fn void GetInterpolatedPoint(vtkIdType cellid, double u, double v)
    * \brief Computes the point using the cubic Hermite weight functions and the two scalar parameters
    * \param vtkIdType cellid, cellid in the srep
    * \param double u scalar for the weight functions
    * \param double v second scalar for the weight functions    n
    * \pre u and v should be between [0, 1]
    */
    VNLType GetInterpolatedSpokeDerivative(vtkIdType cellid, unsigned side, double u, double v);

    /*
    * \fn void GetInterpolatedSpoke(vtkIdType cellid, double u, double v)
    * \brief Computes the point using the cubic Hermite weight functions and the two scalar parameters
    * \param vtkIdType cellid, cellid in the srep
    * \param unsigned spoke position ex: 1-0 [top, bottom]
    * \param double u scalar for the weight functions
    * \param double v second scalar for the weight functions    n
    * \pre u and v should be between [0, 1]
    */
    //VNLType GetInterpolatedSpoke(unsigned side, double u, double v);


    void SetAtomId(vtkIdType atomid){
        m_AtomId = atomid;
    }

    void SetGamma_u(double gammau){
        m_Gamma_u = gammau;
    }

    void SetGamma_v(double gammav){
        m_Gamma_v = gammav;
    }

    void SetSpokeType(vtkIdType spoketype){
        m_SpokeType = spoketype;
    }

    double GetInterpolatedLabel(double u, double v, vtkIdType cellid);

protected:
    vtkSRepInterpolateMedialSpokesHermite();
    ~vtkSRepInterpolateMedialSpokesHermite();

    // Superclass method to update the pipeline
    virtual int RequestData(vtkInformation* request,
                            vtkInformationVector** inputVector,
                            vtkInformationVector* outputVector);

private:
    int m_InterpolationLevel;
    vtkSRep* m_Input;


    vtkIdType m_SpokeType;
    vtkIdType m_AtomId;
    double m_Gamma_u;
    double m_Gamma_v;

    vtkSmartPointer<vtkSRep> m_SRepOutput;


    typedef vector< VNLType** >  HermiteMatrixType;
    vector< HermiteMatrixType > m_HermitMatrices;
    //vector< HermiteMatrixType > m_HermitMatricesRadius;


    /*
    * \fn void GetInterpolatedPoint(double u, double v, Vector3D H[4][4], Vector3D p)
    * \brief Computes the point using the cubic Hermite weight functions and the two scalar parameters
    * \param double u scalar for the weight functions
    * \param double v second scalar for the weight functions
    * \param Vector3D H[4][4] cubic hermite matrix
    * \param Vector3D p result of the interpolation
    * \pre u and v should be between [0, 1]
    */
    VNLType GetInterpolatedSpoke(double u, double v, VNLType** H);




    /*
    * \fn void GetInterpolatedPoint(double u, double v, Vector3D H[4][4], Vector3D p)
    * \brief Computes the point using the cubic Hermite weight functions and the two scalar parameters
    * \param double u scalar for the weight functions
    * \param double v second scalar for the weight functions
    * \param Vector3D H[4][4] cubic hermite matrix
    * \param Vector3D p result of the interpolation
    * \pre u and v should be between [0, 1]
    */
    VNLType GetInterpolatedSpokeDerivative(double u, double v, VNLType** H);


    void GetHermiteMatrix(VNLType **H, VNLType p0, VNLType p1, VNLType p2, VNLType p3, VectorVNLType dp0, VectorVNLType dp1, VectorVNLType dp2, VectorVNLType dp3, VNLType n0, VNLType n1, VNLType n2, VNLType n3, VNLType k0 = VNLType(), VNLType k1 = VNLType(), VNLType k2 = VNLType(), VNLType k3 = VNLType());
    //void GetHermiteMatrix(VNLType **H, VNLType p0, VNLType p1, VNLType p2, VNLType p3, VectorVNLType dp0, VectorVNLType dp1, VectorVNLType dp2, VectorVNLType dp3, VNLType n0, VNLType n1, VNLType n2, VNLType n3);
    void GetHermiteMatrixSpokes(vtkSRep* input, vtkIdType id0, vtkIdType id1, vtkIdType id2, vtkIdType id3, vtkIdType side, VNLType **H);
    //void GetHermiteMatrixRadius(vtkSRep* input, vtkIdType id0, vtkIdType id1, vtkIdType id2, vtkIdType id3, vtkIdType side, VNLType **H);

    void Clear();
};

#endif // vtkSRepInterpolateMedialSpokesHermite_H
