#ifndef VTKSREPVISUPRIMITIVES_H
#define VTKSREPVISUPRIMITIVES_H

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

class vtkSRepVisuPrimitives : public vtkAlgorithm
{
public:
    static vtkSRepVisuPrimitives *New();


    typedef vnl_vector<double> VNLType;
    typedef vector< VNLType > VectorVNLType;


    void SetInput(vtkDataObject* input);
    void SetInput(int index, vtkDataObject* input);

    // Superclass method to update the pipeline
    virtual void Update();

    vtkSmartPointer< vtkActorCollection > GetOuputActors(){
        return m_PrimitivesCollectionActors;
    }

protected:
    vtkSRepVisuPrimitives();

    ~vtkSRepVisuPrimitives();



private:


    vtkSRep* m_SRep;
    vtkSmartPointer< vtkPolyDataCollection > m_PrimitivesCollectionPolyData;
    vtkSmartPointer< vtkCollection > m_PrimitivesCollectionCells;
    vtkSmartPointer< vtkCollection > m_PrimitivesCollectionPoints;
    vtkSmartPointer< vtkActorCollection >  m_PrimitivesCollectionActors;

};

#endif // VTKSREPVISUPRIMITIVES_H
