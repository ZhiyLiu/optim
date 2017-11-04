#ifndef VTKREADSREP_H
#define VTKREADSREP_H

#include "vtksrep.h"
#include "vtkDataReader.h"
#include "vtkSmartPointer.h"

#include "P3DControl.h"
#include "ControlParms.h"

using namespace std;

class vtkReadSRep : public vtkDataReader
{
public:
    static vtkReadSRep *New();

    virtual void Update();

    vtkSmartPointer<vtkSRep> GetOutput(){
        return m_Srep;
    }

    M3DQuadFigure* GetQuadFigure(){
        return dynamic_cast<M3DQuadFigure*>(GetFigure());
    }

protected:
    vtkReadSRep();

    ~vtkReadSRep();

private:

    vtkSmartPointer<vtkSRep> m_Srep;


    int globalVerbosity;			// Current verbosity level of Pablo
    P3DControl* p3d;
    M3DFigure* GetFigure();

};

#endif // VTKSREP_H
