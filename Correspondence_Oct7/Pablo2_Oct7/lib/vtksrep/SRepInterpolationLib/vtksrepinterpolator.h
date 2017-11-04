#ifndef vtkSRepInterpolator_H
#define vtkSRepInterpolator_H

#include <vtkSmartPointer.h>

#include "vnl/vnl_vector.h"

#include "vtksrep.h"

#include "vtkAlgorithm.h"

#include "iostream"
using namespace std;


#include "vtksrepinterpolatemedialspokeshermite.h"
#include "vtksrepinterpolatemedialsheet.h"
#include "vtksrepinterpolatemedialcrestcurve.h"
#include "vtksrepinterpolatecrestspokesquartic.h"


class vtkSRepInterpolator : public vtkAlgorithm
{
public:
    static vtkSRepInterpolator *New();


    typedef vnl_vector<double> VNLType;
    typedef vector< VNLType > VectorVNLType;
    typedef vector< VectorVNLType > VectorVectorVNLType;

    typedef vtkSRep::VNLMatrixType VNLMatrixType;
    typedef vector< VNLMatrixType > VectorVNLMatrixType;

    typedef vtkSRep::VectorIdsType VectorIdsType;
    typedef vector< VectorIdsType > VectorVectorIdsType;

    /*
    * \fn void GetInterpolatedPoint(double u, double v, double t)
    * \brief Computes the point using the s-rep interpolators
    * \param double u [0, 1]
    * \param double v [0, 1]
    * \param double t [0, 1]
    * \pre u, v and t have to be between [0, 1]
    */
    VNLType GetInterpolatedPoint(double u, double v, double tau);

    /*
    * \fn void GetInterpolatedPoint(double u, double v, double t)
    * \brief Computes the point using the s-rep interpolators
    * \param double u [0, 1]
    * \param double v [0, 1]
    * \pre u, v and t have to be between [0, 1]
    */
    VNLType GetInterpolatedSpoke(double u, double v);

    //void SetInput(vtkDataObject* input);
    //void SetInput(int index, vtkDataObject* input);
    void SetInput(vtkSRep* srep){
        m_Input = srep;
    }

    vtkPolyData* TestPoly(){
        return m_TestPoly;
    }

    virtual void Update();

    double GetInterpolatedLabel(){
        return m_InterpolatedLabel;
    }


protected:
    vtkSRepInterpolator();
    ~vtkSRepInterpolator();

private:


    vtkSmartPointer<vtkPolyData> m_UVToCrestId;

    vector< vtkSmartPointer<vtkCellLocator> > m_CellLocators;

    vtkSRep* m_Input;
    vtkPolyData* m_TestPoly;

    void GenerateUVMap();
    vtkIdType GetCellId(double& u, double& v, int& loc);


    vtkSmartPointer< vtkSRepInterpolateMedialSpokesHermite > m_Medialspokesinterpolator;
    vtkSmartPointer< vtkSRepInterpolateMedialSheet > m_Medialsheetinterpolator;
    vtkSmartPointer<vtkSRepInterpolateCrestSpokesQuartic> m_Interpolatecrestspokes;
    vtkSmartPointer<vtkSRepInterpolateCrestSpokesQuartic> m_InterpolatecrestspokesBottom;
    double m_InterpolatedLabel;

};

#endif // vtkSRepInterpolator_H
