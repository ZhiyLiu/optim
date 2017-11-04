#ifndef VTKTRANSFORMSREPFILTER_H
#define VTKTRANSFORMSREPFILTER_H

#include "vtkTransformPolyDataFilter.h"
#include "vtksrep.h"

class vtkTransformSrepFilter : public vtkTransformPolyDataFilter
{
public:
    static vtkTransformSrepFilter *New();
    vtkTypeMacro(vtkTransformSrepFilter,vtkTransformPolyDataFilter);
    void PrintSelf(ostream& os, vtkIndent indent);

    vtkSmartPointer<vtkSRep> GetOutput(){
        return m_SRepOutput;
    }

protected:

    vtkTransformSrepFilter();
    ~vtkTransformSrepFilter();

    int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);

private:
    vtkSmartPointer<vtkSRep> m_SRepOutput;

};

#endif // VTKTRANSFORMSREPFILTER_H
