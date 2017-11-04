// Command as: ./drawSrep -m /NIRAL/work/ltu/WorkSpace/Test/Feb20_forpaper/test1/dataSet/moved_sreps/1.m3d



#include "visualization.h"

using namespace std;


typedef vtkSmartPointer< vtkPoints > hp;
typedef vector<hp> hubposVector;

int main(int argc, char *argv[])
{
    vector< std::string > modelname;

    int modelnums = argc - 2;//arg[0] is exe name, arg[1] is "-m", the following is the model names.

    for(int i = 1; i <= modelnums; i++){
        modelname.push_back(argv[1+i]);
    }


    vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
    renderer->SetBackground(1,1,1);
    vtkSmartPointer<vtkRenderWindow> renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
    renderWindow->AddRenderer(renderer);
    vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
    renderWindowInteractor->SetRenderWindow(renderWindow);

    toolsfunc tf;

    int side =0;
    int interpolationLevel =0;
    int quadtype = 1;
    bool moved = true;
    double quadColor[3] = {1,0.5,0};
    int quadNums;



    //vtkSmartPointer< vtkPoints > hubPosition_b_up,hubPosition_s_up,hubPosition_b_moved_up,hubPosition_s_moved_up;

    //Loop each of the m3d files.
    cout<<"You input "<<modelname.size()<<" models"<<endl;
    for(int q=0; q< modelname.size();q++){

        M3DQuadFigure* quadFig =  tf.GetQuadFigure(modelname[q].c_str());

        /*vtkSmartPointer< vtkPoints > hubPosition_b_up = vtkSmartPointer< vtkPoints >::New();
        //vtkSmartPointer< vtkPoints > hubPosition_s_up = vtkSmartPointer< vtkPoints >::New();

        // Up side
        visualization visualObject_up(quadFig, interpolationLevel, side, moved, quadColor, renderer);

        hubPosition_b_up = visualObject_up.getSubQuadsPosition(0);
        //hubPosition_s_up = visualObject.getSubQuadsPosition(1);

        // Draw up boundary frame.
        visualObject_up.drawFrameOfQuads(hubPosition_b_up);


        // Down side
        vtkSmartPointer< vtkPoints > hubPosition_b_down = vtkSmartPointer< vtkPoints >::New();
        side = 1;
        quadColor[0] = 0.8; quadColor[1] = 0; quadColor[2] = 0;
        visualization visualObject_down(quadFig, interpolationLevel, side, moved, quadColor, renderer);

        hubPosition_b_down = visualObject_down.getSubQuadsPosition(0);

        // Draw down boundary frame.
        visualObject_down.drawFrameOfQuads(hubPosition_b_down);*/



        // Crest
        visualization visualObject_crest(quadFig, interpolationLevel, side, quadColor, renderer);

        //visualObject_crest.drawCrest(true, true, true);

        // Draw the line between up-boundary-crest points and medial-crest points; and down-boundary-crest and medial-crest.
       // visualObject_crest.drawCrestCorrespondenceLines(true,true);

        //visualObject_crest.drawUpOrDownBoundaryQuadLine_method1(true, false);

        // getAllSubQuadsVertexesOnBoundary() should be called before drawUpOrDownBoundaryQuadLine_method2.
        //visualObject_crest.getAllSubQuadsVertexesOnBoundary();
        //visualObject_crest.drawUpOrDownBoundaryQuadLine_method2(true, 50, false, 50);




    }






    renderWindow->Render();
    renderWindowInteractor->Start();


    return 0;
}






