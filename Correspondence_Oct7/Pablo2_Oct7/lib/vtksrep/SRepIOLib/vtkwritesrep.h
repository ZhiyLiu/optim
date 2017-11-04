#ifndef vtkWriteSRep_H
#define vtkWriteSRep_H

#include "vtksrep.h"
#include "vtkDataWriter.h"
#include "vtkSmartPointer.h"

#include "P3DControl.h"
#include "ControlParms.h"

using namespace std;

class vtkWriteSRep : public vtkDataWriter
{
public:
    static vtkWriteSRep *New();

    virtual int Write();

    M3DFigure* NewQuadFigure(const int numr, int numc);

protected:
    vtkWriteSRep();

    ~vtkWriteSRep();

private:

    vtkSmartPointer<vtkSRep> m_Srep;

    int globalVerbosity;			// Current verbosity level of Pablo
    P3DControl* p3d;

};

#endif // VTKSREP_H
