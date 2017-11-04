#ifndef vtkSRepToGridLayout_H
#define vtkSRepToGridLayout_H

#include "vtkAlgorithm.h"
#include "vtksrep.h"

#include <vector>

#include "vnl/vnl_vector.h"


#include "vtkPolyDataCollection.h"
#include "vtkMapperCollection.h"
#include "vtkActorCollection.h"
#include "vtkDataArrayCollection.h"

#include "vtkSmartPointer.h"

using namespace std;

class vtkSRepToGridLayout : public vtkAlgorithm
{
public:
    static vtkSRepToGridLayout *New();


    typedef vnl_vector<double> VNLType;
    typedef vector< VNLType > VectorVNLType;


    void SetInput(vtkDataObject* input);
    void SetInput(int index, vtkDataObject* input);

    vtkSRep* GetOutput(){
        return m_Output;
    }

    // Superclass method to update the pipeline
    virtual void Update();


protected:
    vtkSRepToGridLayout();

    ~vtkSRepToGridLayout();



private:


    vtkSmartPointer<vtkSRep> m_Output;


};

#endif // vtkSRepToGridLayout_H
