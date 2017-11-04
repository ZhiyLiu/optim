/* This class draw suface points in correspondence color.
 * The input is vtk point file which storing all the surface point of a s-rep.
 * ./showPointsCluster /NIRAL/work/ltu/WorkSpace/May15/Pablo2/lib/vtksrep/Correspondence/deform_mean/meanShape.vtk
 * Liyun Tu
 * Jun 26, 2014
*/



#include "shiftedsrep.h"

using namespace std;

//const int rowNum = 3;
//const int colNum = 13;


int main(int argc, char *argv[])
{
    vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
    renderer->SetBackground(1,1,1);
    vtkSmartPointer<vtkRenderWindow> renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
    renderWindow->AddRenderer(renderer);
    vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
    renderWindowInteractor->SetRenderWindow(renderWindow);

    // rotate the picture, up or down. for hipcampus
    double rX = 0;//195;
    double rY = 0;//55;
    double rZ = 0;
    double opacity = 1;

    shiftedsrep shift_srep_obj(rX, rY, rZ, opacity);

    shift_srep_obj.drawSurfacePoint(argv[1], renderer);



    renderWindow->Render();
    renderWindowInteractor->Start();

    return 0;
}




