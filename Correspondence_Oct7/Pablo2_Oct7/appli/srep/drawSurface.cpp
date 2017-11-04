/* Draw the implied surface of given s-rep.
 * Command: ./drawSurface /NIRAL/work/ltu/WorkSpace/Test/30AllAtomsFromAtomStage/input/afterProcrustes/T0068-1-2-02_segmented_vent_pp_surfSPHARM_scale_atomStage_S.m3d
 * Liyun Tu
 * May 20, 2014
*/

#include "visualizesrepsurface.h"


int main(int argc, char *argv[])
{
    const char* srepfilename = argv[1];

    vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
    renderer->SetBackground(1,1,1);
    vtkSmartPointer<vtkRenderWindow> renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
    renderWindow->AddRenderer(renderer);
    vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
    renderWindowInteractor->SetRenderWindow(renderWindow);

    toolsfunc tf;
    int interpolationLevel = 2;

    double quadColor[3] = {1,1,1};

    // rotate the picture, up or down. for hipcampus
    double rX = 0;//195;
    double rY = 0;//55;
    double rZ = 0;

    visualizesrepsurface obj;
    //obj.drawPointsSetSurface(srepfilename, interpolationLevel, renderer);

    // Draw the spokes and medial sheet in color
    //obj.drawCorrespondenceSpokes(srepfilename, interpolationLevel, renderer);

    vtkSmartPointer< vtkPoints > ps = obj.getSurfacePointsSet(srepfilename,interpolationLevel);
    obj.drawPointsSetSurface_zcolor(ps, renderer);


    //obj.drawSpokes(srepfilename, interpolationLevel, renderer);





    renderWindow->Render();
    renderWindowInteractor->Start();

    return 0;
}











