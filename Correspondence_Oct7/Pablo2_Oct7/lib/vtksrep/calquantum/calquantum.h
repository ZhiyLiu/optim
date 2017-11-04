#include <vtkSmartPointer.h>

#include <vtksrep.h>
#include <vtksrepinterpolatemedialsheet.h>
#include <vtksrepvisuprimitives.h>

#include <vtkCellArray.h>
#include <vtkQuad.h>
#include <vtkLine.h>

#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkProperty.h>

#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderWindow.h>

#include <vtksrepinterpolatemedialcrestcurve.h>
#include <vtksrepinterpolatemedialspokes.h>
#include <vtksrepinterpolatemedialspokeshermite.h>
#include <vtksrepinterpolatecrestspokes.h>
#include <vtksrepinterpolatecrestspokesquartic.h>




double calquantum(vtkSmartPointer<vtkSRep> srepfig, vtkSmartPointer<vtkSRepInterpolateMedialSheet> medialsheetinterpolator, vtkSmartPointer<vtkSRep> srepcrest, vtkSmartPointer<vtkSRepInterpolateMedialSpokesHermite> medialspokesinterpolator);
