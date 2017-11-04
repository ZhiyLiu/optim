#ifndef vtkSRepInterpolateCrestSpokesQuartic_H
#define vtkSRepInterpolateCrestSpokesQuartic_H

#include "vtkPolyDataAlgorithm.h"
#include "vnl/vnl_vector.h"
#include "vnl/vnl_cross.h"
#include "vnl/vnl_numeric_traits.h"

#include "vtksrep.h"
#include "vtkinterpolatecurve.h"
#include "vtkSmartPointer.h"
#include "vtksrepinterpolatemedialspokeshermite.h"

using namespace std;

class vtkSRepInterpolateCrestSpokesQuartic : public vtkPolyDataAlgorithm
{
public:
    static vtkSRepInterpolateCrestSpokesQuartic *New();


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

    vtkSmartPointer<vtkSRep> GetSRepOutput(){
        return m_SRepOutput;
    }

    void SetAtomId(vtkIdType atomid){
        m_AtomId = atomid;
    }

    void SetGamma_t(double gammat){
        m_Gamma_t = gammat;
    }

    void SetGamma_theta(double gammatheta){
        m_Gamma_theta = gammatheta;
    }

    void SetSpokeType(vtkIdType spoketype){
        m_SpokeType = spoketype;
    }

    /*
      Set to true for quasi tube interpolation
    */
    void SetUseAllSpokes(bool setall){
        m_UseAllSpokes = setall;
    }
    void SetCyclicSpokes(bool cyclic){
        m_CyclicSpokes = cyclic;
    }
    void SetCyclicCurve(bool cyclic){
        m_CyclicCurve = cyclic;
    }
    void SetSphericalCaps(bool sphericalcaps){
        m_SphericalCaps = sphericalcaps;
    }

    void SetDoubleEndCrestAtoms(bool sides){
        m_DoubleEndCrestAtoms = sides;
    }

    vtkSRep::VNLType GetInterpolatedSpoke(vtkIdType cellid, double t, double v);
    vtkSRep::VNLType GetInterpolatedSpoke(double t, double v);
    vtkSRep::VNLType GetInterpolatedSpoke(double v);
    vtkSRep::VNLType GetInterpolatedPoint(vtkIdType cellid, double t);
    vtkSRep::VNLType GetInterpolatedPoint(double t);


protected:
    vtkSRepInterpolateCrestSpokesQuartic();

    ~vtkSRepInterpolateCrestSpokesQuartic();

    // Superclass method to update the pipeline
    virtual int RequestData(vtkInformation* request,
                            vtkInformationVector** inputVector,
                            vtkInformationVector* outputVector);

private:


    int m_InterpolationLevel;
    vtkSRep *m_Input;
    vtkSmartPointer<vtkSRep> m_SRepOutput;

    vtkIdType m_AtomId;
    vtkIdType m_SpokeType;
    double m_Gamma_t;
    double m_Gamma_theta;
    bool m_UseAllSpokes;
    bool m_CyclicSpokes;
    bool m_DoubleEndCrestAtoms;
    bool m_CyclicCurve;
    bool m_SphericalCaps;

    vtkSmartPointer< vtkInterpolateCurve > m_MedialCrestCurveInterpolator;
    vtkSRep::VNLType m_P0;
    vector< vtkSmartPointer< vtkInterpolateCurve > > m_Spokesinterpolator;    
    vtkSmartPointer< vtkInterpolateCurve > m_Inplanenterpolate;

    void InterpolateCrest(vector<vtkSRep::SPOKES_TYPE> spokestype,
                          vtkSmartPointer<vtkInterpolateCurve> medialcrestcurveinterpolator,
                          vector<vtkSmartPointer<vtkInterpolateCurve> > spokesinterpolator,                          
                          int crestpos, int spokepos);


};

#endif // VTKSREPINTERPOLATEMEDIALSHEET_H
