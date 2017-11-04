#include <vtkVersion.h>
#include <vtkSmartPointer.h>
#include <vtkActor.h>
#include <vtkSphereSource.h>
#include <vtkRendererCollection.h>
#include <vtkCellArray.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkObjectFactory.h>
#include <vtkPlaneSource.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkPropPicker.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>

// Handle mouse events
class MouseInteractorStyle2 : public vtkInteractorStyleTrackballCamera
{
  public:
    static MouseInteractorStyle2* New();
    vtkTypeMacro(MouseInteractorStyle2, vtkInteractorStyleTrackballCamera);

    virtual void OnLeftButtonDown()
    {
      int* clickPos = this->GetInteractor()->GetEventPosition();

      // Pick from this location.
      vtkSmartPointer<vtkPropPicker>  picker =
        vtkSmartPointer<vtkPropPicker>::New();
      picker->Pick(clickPos[0], clickPos[1], 0, this->GetDefaultRenderer());

      double* pos = picker->GetPickPosition();
      std::cout << "Pick position (world coordinates) is: "
                << pos[0] << " " << pos[1]
                << " " << pos[2] << std::endl;

      std::cout << "Picked actor: " << picker->GetActor() << std::endl;
      //Create a sphere
      vtkSmartPointer<vtkSphereSource> sphereSource =
        vtkSmartPointer<vtkSphereSource>::New();
      sphereSource->SetCenter(pos[0], pos[1], pos[2]);
      sphereSource->SetRadius(0.1);

      //Create a mapper and actor
      vtkSmartPointer<vtkPolyDataMapper> mapper =
        vtkSmartPointer<vtkPolyDataMapper>::New();
      mapper->SetInputConnection(sphereSource->GetOutputPort());

      vtkSmartPointer<vtkActor> actor =
        vtkSmartPointer<vtkActor>::New();
      actor->SetMapper(mapper);


      //this->GetInteractor()->GetRenderWindow()->GetRenderers()->GetDefaultRenderer()->AddActor(actor);
      this->GetDefaultRenderer()->AddActor(actor);
      // Forward events
      vtkInteractorStyleTrackballCamera::OnLeftButtonDown();
    }

  private:

};

vtkStandardNewMacro(MouseInteractorStyle2);

// Execute application.
int main(int, char *[])
{
  vtkSmartPointer<vtkPlaneSource> planeSource =
    vtkSmartPointer<vtkPlaneSource>::New();
  planeSource->Update();

  // Create a polydata object
  vtkPolyData* polydata = planeSource->GetOutput();

  // Create a mapper
  vtkSmartPointer<vtkPolyDataMapper> mapper =
    vtkSmartPointer<vtkPolyDataMapper>::New();
#if VTK_MAJOR_VERSION <= 5
  mapper->SetInput ( polydata );
#else
  mapper->SetInputData ( polydata );
#endif

  // Create an actor
  vtkSmartPointer<vtkActor> actor =
    vtkSmartPointer<vtkActor>::New();
  actor->SetMapper ( mapper );

  std::cout << "Actor address: " << actor << std::endl;

  // A renderer and render window
  vtkSmartPointer<vtkRenderer> renderer =
    vtkSmartPointer<vtkRenderer>::New();
  vtkSmartPointer<vtkRenderWindow> renderWindow =
    vtkSmartPointer<vtkRenderWindow>::New();
  renderWindow->AddRenderer ( renderer );

  // An interactor
  vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor =
    vtkSmartPointer<vtkRenderWindowInteractor>::New();
  renderWindowInteractor->SetRenderWindow ( renderWindow );

  // Set the custom stype to use for interaction.
  vtkSmartPointer<MouseInteractorStyle2> style =
    vtkSmartPointer<MouseInteractorStyle2>::New();
  style->SetDefaultRenderer(renderer);

  renderWindowInteractor->SetInteractorStyle( style );

  // Add the actors to the scene
  renderer->AddActor ( actor );
  renderer->SetBackground ( 0,0,1 );

  // Render and interact
  renderWindow->Render();
  renderWindowInteractor->Initialize();
  renderWindowInteractor->Start();

  return EXIT_SUCCESS;
}
/*
#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkPolyDataReader.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkArrowSource.h>
#include <vtkPoints.h>

#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkProperty.h>

#include <vtkLookupTable.h>
#include <vtkColorTransferFunction.h>
#include <vtkCurvatures.h>
#include <vtkPointData.h>

#include <vtkSmoothPolyDataFilter.h>

#include <vnl/vnl_vector.h>

#define PI 3.14159265

using namespace std;

int main(int argc, char *argv[])
{


    vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
    renderer->SetBackground(1,1,1);
    vtkSmartPointer<vtkRenderWindow> renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
    renderWindow->AddRenderer(renderer);
    vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
    renderWindowInteractor->SetRenderWindow(renderWindow);


    vtkSmartPointer<vtkPolyDataReader> reader =
    vtkSmartPointer<vtkPolyDataReader>::New();

    reader->SetFileName(argv[1]);
    reader->Update();
    vtkPolyData* orig = reader->GetOutput();


    vtkSmartPointer<vtkPolyDataReader> reader1 =
    vtkSmartPointer<vtkPolyDataReader>::New();

    reader1->SetFileName(argv[2]);
    reader1->Update();
    vtkPolyData* inflated = reader1->GetOutput();

    vtkPoints* pointsorig = orig->GetPoints();
    vtkPoints* pointsinflated = inflated->GetPoints();


    vtkSmartPointer<vtkPolyDataMapper> inflatedmapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    inflatedmapper->SetInput(inflated);
    vtkSmartPointer<vtkActor>  vtkactorinflated = vtkActor::New();
    vtkactorinflated->SetMapper(inflatedmapper);
    //vtkactorinflated->GetProperty()->SetOpacity(0.2);



    vtkSmartPointer<vtkPolyDataMapper> origmapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    origmapper->SetInput(orig);
    vtkSmartPointer<vtkActor>  vtkactororig = vtkActor::New();
    vtkactororig->SetMapper(origmapper);




    vtkSmartPointer<vtkSmoothPolyDataFilter> smoothpoly = vtkSmartPointer<vtkSmoothPolyDataFilter>::New();
    smoothpoly->SetInput(orig);
    smoothpoly->SetNumberOfIterations(100);
    smoothpoly->SetRelaxationFactor(0.025);
    smoothpoly->Update();


    vtkSmartPointer<vtkCurvatures> curv = vtkSmartPointer<vtkCurvatures>::New();
    curv->SetInput(smoothpoly->GetOutput());
    //curv->SetInput(orig);
    curv->SetCurvatureTypeToMean();
    curv->Update();

    vtkSmartPointer<vtkColorTransferFunction> color = vtkSmartPointer<vtkColorTransferFunction>::New();
    color->SetColorSpaceToRGB();

    //double* curvrange = curv->GetOutput()->GetPointData()->GetScalars()->GetRange();

    color->AddRGBPoint(-0.25,1,0,0);
    color->AddRGBPoint(0.25,0,1,0);

    color->Build();


    renderer->AddActor(vtkactorinflated);
    //renderer->AddActor(vtkactororig);

    orig->GetPointData()->SetScalars(curv->GetOutput()->GetPointData()->GetScalars());
    origmapper->SetLookupTable(color);

    inflated->GetPointData()->SetScalars(curv->GetOutput()->GetPointData()->GetScalars());
    inflatedmapper->SetLookupTable(color);

    for(int i = 0; i < orig->GetNumberOfPoints(); i+=20){
    //for(int i = 1000; i < 1020; i++){

        double* porig = pointsorig->GetPoint(i);
        double* pinflated = pointsinflated->GetPoint(i);


        vnl_vector<double> vorig(3);
        vorig.set(porig);

        vnl_vector<double> vinflated(3);
        vinflated.set(pinflated);


        vnl_vector<double> vdir = vinflated - vorig;
        double magnitud = vdir.magnitude();
        vdir.normalize();


        vtkSmartPointer<vtkPolyDataMapper> actormapper = vtkSmartPointer<vtkPolyDataMapper>::New();
        vtkSmartPointer<vtkActor>  vtkactor = vtkActor::New();
        vtkactor->SetMapper(actormapper);

        vtkSmartPointer<vtkArrowSource> vtkarrow		= vtkSmartPointer<vtkArrowSource>::New();
        vtkarrow->SetTipResolution(30);
        vtkarrow->SetShaftResolution( 30  );
        //vtkarrow->SetTipLength(magnitud);
        //vtkarrow->SetShaftRadius(2);
        actormapper->SetInput( vtkarrow->GetOutput() );

        double scale[3];
        scale[0] = 5;
        scale[1] = 5;
        scale[2] = 5;
        vtkactor->SetScale(scale);


        double orientation[3] = {0,0,0};
        vtkactor->SetOrientation(orientation);

        double valxy = 0, valxz = 0, valxyz = 0;
        double angle0 = atan2(vdir[2], vdir[0]) * 180.0 / PI;

        valxz = sqrt(pow(vdir[0],2) + pow(vdir[2],2));
        valxyz = sqrt(pow(vdir[0],2) + pow(vdir[1],2) + pow(vdir[2],2));
        valxy = sqrt(pow(valxyz,2) - pow(valxz,2));
        double angle1 = atan2(valxy, valxz) * 180.0 / PI;

        vtkactor->RotateY(-angle0);
        if(vdir[1] < 0){
                vtkactor->RotateZ(-angle1);
        }else{
                vtkactor->RotateZ(angle1);
        }


        vtkactor->SetPosition(vorig[0], vorig[1], vorig[2]);
        //renderer->AddActor(vtkactor);

    }


    renderWindow->Render();

    renderWindowInteractor->Start();

    return 0;
}*/
