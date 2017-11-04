/* Command this:
 * ltu@orion(53) loop-1/input_sreps$ /NIRAL/work/ltu/WorkSpace/data_pablo/binPablo2/drawFigure -m 10.m3d 16.m3d 21.m3d 27.m3d
 * 32.m3d 38.m3d 43.m3d 49.m3d 7.m3d
 **/

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

    hubposVector b;
    hubposVector s;
    vector<visualization> viObjList;

    //vtkSmartPointer< vtkPoints > hubPosition_b_up,hubPosition_s_up,hubPosition_b_moved_up,hubPosition_s_moved_up;

    //Loop each of the m3d files.
    cout<<"You input "<<modelname.size()<<" models"<<endl;
    for(int q=0; q< modelname.size();q++){

        M3DQuadFigure* quadFig =  tf.GetQuadFigure(modelname[q].c_str());
        quadNums = (quadFig->getRowCount()-1)*(quadFig->getColumnCount()-1);//16
        vtkSmartPointer< vtkPoints > hubPosition_b = vtkSmartPointer< vtkPoints >::New();
        vtkSmartPointer< vtkPoints > hubPosition_s = vtkSmartPointer< vtkPoints >::New();

        visualization visualObject(quadFig, interpolationLevel, 0, moved, quadColor, renderer,0,0,0,0);

        viObjList.push_back(visualObject);

        //Calculate quads areas.
        hubPosition_b = visualObject.getSubQuadsPosition(0);
        hubPosition_s = visualObject.getSubQuadsPosition(1);

        b.push_back(hubPosition_b);
        s.push_back(hubPosition_s);

/*   int q =0;
        cout<<"The following is quadfig "<< q << endl;
        M3DQuadFigure* quadFig =  tf.GetQuadFigure(modelname[q].c_str());
        visualization visualObject(quadFig, interpolationLevel, side, moved, quadColor, renderer);

        quadFig =  tf.GetQuadFigure(modelname[q+1].c_str());
        visualization visualObject_up(quadFig, interpolationLevel, 0, moved, quadColor, renderer);

        //input as a 3th m3d, this want to show the moved_up spokes.
        quadFig =  tf.GetQuadFigure(modelname[q+2].c_str());
        visualization visualObject_moved_up(quadFig, interpolationLevel, 0, moved, quadColor, renderer);

        //Calculate quads areas.
        hubPosition_b = visualObject.getSubQuadsPosition(0);
        hubPosition_s = visualObject.getSubQuadsPosition(1);

        hubPosition_b_up = visualObject_up.getSubQuadsPosition(0);
        hubPosition_s_up = visualObject_up.getSubQuadsPosition(1);

        hubPosition_b_moved_up = visualObject_moved_up.getSubQuadsPosition(0);
        hubPosition_s_moved_up = visualObject_moved_up.getSubQuadsPosition(1);



        //double quadColor[3] = {1,0.5,0};
        quadColor[0] =1+0.4*q;
        quadColor[1] =0.6*q;
        quadColor[2] =0;
        visualObject.setColor(quadColor);
        //visualObject.drawEdgeOfQuads_method3(hubPosition_b);

        quadColor[0] =0.35;
        quadColor[1] =0.16;
        quadColor[2] =0.02;
        visualObject.setColor(quadColor);
        visualObject.drawEdgeOfQuads_method3(hubPosition_s_up);
        //visualObject_up.drawEdgeOfQuads_method3(hubPosition_s_up);
        //visualObject_moved_up.drawEdgeOfQuads_method3(hubPosition_s_moved_up);

        //(1,0.6,0) dark yellow;
    for(int i =0;i<16;i++){
       //cout<<"Drawing quad: "<<i<<endl;
        //quadColor[0] =1;
        //quadColor[1] =0;
        //quadColor[2] =0.4+0.1*i; //fen se
        quadColor[0] =0;
        quadColor[1] =1;
        quadColor[2] =0.4;
        visualObject.drawSpoke(hubPosition_b,hubPosition_s,i);
        //visualObject.drawFrameOfQuadByIndex(hubPosition_s,i);
        quadColor[0] =0;
        quadColor[1] =1;
        quadColor[2] =0;
        //visualObject.drawFrameOfQuadByIndex(hubPosition_b,i);
        quadColor[0] =0;
        quadColor[1] =0;
        quadColor[2] =1;
        //visualObject_moved_up.drawSpoke(hubPosition_b_moved_up,hubPosition_s_moved_up,i);
        //qian hui (0.5,0.5,0.4)
        visualObject_moved_up.setColor(quadColor);
        //visualObject_moved_up.drawFrameOfQuadByIndex(hubPosition_s_moved_up,i);
        //visualObject.drawVolume(hubPosition_b,hubPosition_s,i);
        //visualObject.drawFrameOfQuadByIndex(hubPosition_b,i);
        //drawSpoke(quadfig, renderer, side, quadColor, interpolationLevel, i, true);
        quadColor[0] =1;
        quadColor[1] =0;
        quadColor[2] =0.5;
        visualObject_up.drawSpoke(hubPosition_b_up,hubPosition_s_up,i);
        //green color(0.3,0.8,0)
        visualObject_up.setColor(quadColor);
        //visualObject_up.drawFrameOfQuadByIndex(hubPosition_s_up,i);
        //visualObject.drawSpoke(hubPosition_b,hubPosition_s,i);

        // drawFrameOfQuadByIndex(quadfig, renderer, side, 0, quadColor, interpolationLevel, true, i);
        //drawVolume(quadfig, renderer, side, quadColor,interpolationLevel,i,true);
    }*/


   }

    for(int i =0;i<quadNums;i++){
       //cout<<"Drawing quad: "<<i<<endl;
        //quadColor[0] =1;
        //quadColor[1] =0;
        //quadColor[2] =0.4+0.1*i; //fen se

        for(int q=0; q < modelname.size();q++){
            quadColor[0] =1;//+0.02*q;//0.05+i*0.06;
            quadColor[1] =0.1;
            quadColor[2] =1;
            viObjList[q].setColor(quadColor);
            viObjList[q].drawSpoke(b[q],s[q],i);

            quadColor[0] =0;
            quadColor[1] =0.5;        //(0,0.2,0) dark green.
            quadColor[2] =0;
            viObjList[q].setColor(quadColor);            
            viObjList[q].drawFrameOfQuads(b[q]);
            quadColor[0] =0;
            quadColor[1] =0;        //(0,0.2,0) dark green.
            quadColor[2] =0;
            viObjList[q].setColor(quadColor);
            // Only draw the quads (without the subquads, no matter what's the interpolation level)
            viObjList[q].drawEdgeOfQuads_method3(b[q]);

            quadColor[0] =1;
            quadColor[1] =0;  //(0,0.2,0) dark green.(0,0.8,0) bright green.
            quadColor[2] =0;
            viObjList[q].setColor(quadColor);            
            viObjList[q].drawFrameOfQuads(s[q]);
            quadColor[0] =0;
            quadColor[1] =0;        //(0,0.2,0) dark green.
            quadColor[2] =1;
            viObjList[q].setColor(quadColor);
            // Only draw the quads (without the subquads, no matter what's the interpolation level)
            viObjList[q].drawEdgeOfQuads_method3(s[q]);
        }
   }

    //draw each atom's spoke in different color. one atom index one color.
    /*for(int q=0; q < modelname.size();q++){
        for(int i =0;i<quadNums;i++){
            int j = (i+1)%6 ;
            if(j==1){
                quadColor[0] =0.8;
                quadColor[1] =0;//(0,0.2,0) dark green.(0,0.4,0) bright green.
                quadColor[2] =0;
            }
            if(j==5){
                quadColor[0] =0.3;
                quadColor[1] =0;//(0,0.2,0) dark green.(0,0.4,0) bright green.
                quadColor[2] =0;
            }
            if(j==4){
                quadColor[0] =0;
                quadColor[1] =0.6;//(0,0.2,0) dark green.(0,0.4,0) bright green.
                quadColor[2] =0;
            }
            if(j==2){
                quadColor[0] =0;
                quadColor[1] =0.9;//(0,0.2,0) dark green.(0,0.4,0) bright green.
                quadColor[2] =0;
            }
            if(j==0){
                quadColor[0] =0;
                quadColor[1] =0;//(0,0.2,0) dark green.(0,0.4,0) bright green.
                quadColor[2] =0.7;
            }
            if(j==3){
                quadColor[0] =0;
                quadColor[1] =0;//(0,0.2,0) dark green.(0,0.4,0) bright green.
                quadColor[2] =0.4;
            }

            viObjList[q].setColor(quadColor);
            viObjList[q].drawSpoke(b[q],s[q],i);
        }
    }*/

    //only draw the middle 4 quads of the ellipsoid.
    /*for(int i =0;i<quadNums;i++){
        if(i==3||i==4||i==11||i==12){
            for(int q=0; q < modelname.size();q++){
                quadColor[0] =0;
                quadColor[1] =0.6;
                quadColor[2] =0;
                viObjList[q].setColor(quadColor);
                viObjList[q].drawSpoke(b[q],s[q],i);

                quadColor[0] =0;
                quadColor[1] =0.7;//(0,0.2,0) dark green.
                quadColor[2] =0;
                viObjList[q].setColor(quadColor);
                viObjList[q].drawFrameOfQuadByIndex(b[q],i);

                quadColor[0] =0;
                quadColor[1] =0;//(0,0.2,0) dark green.
                quadColor[2] =0.5;
                viObjList[q].setColor(quadColor);
                viObjList[q].drawFrameOfQuadByIndex(s[q],i);
            }
        }
    }*/


    //For test 2: lateral ventricles. Draw spoke's in diff colors.
    /*for(int q=0; q < modelname.size();q++){
        for(int i =0;i<quadNums;i++){
            int j = (i+1)%8;
            if(j==1){
                quadColor[0] =0.8;
                quadColor[1] =0;//(0,0.2,0) dark green.(0,0.4,0) bright green.
                quadColor[2] =0;
            }
            if(j==5){
                quadColor[0] =0.3;
                quadColor[1] =0;//(0,0.2,0) dark green.(0,0.4,0) bright green.
                quadColor[2] =0;
            }
            if(j==4){
                quadColor[0] =0;
                quadColor[1] =0.4;//(0,0.2,0) dark green.(0,0.4,0) bright green.
                quadColor[2] =0;
            }
            if(j==2){
                quadColor[0] =0;
                quadColor[1] =0.7;//(0,0.2,0) dark green.(0,0.4,0) bright green.
                quadColor[2] =0;
            }
            if(j==0){
                quadColor[0] =0;
                quadColor[1] =0;//(0,0.2,0) dark green.(0,0.4,0) bright green.
                quadColor[2] =0.7;
            }
            if(j==3){
                quadColor[0] =0.7;
                quadColor[1] =0;//(0,0.2,0) dark green.(0,0.4,0) bright green.
                quadColor[2] =0;
            }
            if(j==7){
                quadColor[0] =0.2;
                quadColor[1] =0.7;//(0,0.2,0) dark green.(0,0.4,0) bright green.
                quadColor[2] =0;
            }
            if(j==6){
                quadColor[0] =0;
                quadColor[1] =0;//(0,0.2,0) dark green.(0,0.4,0) bright green.
                quadColor[2] =0.2;
            }

            viObjList[q].setColor(quadColor);
            viObjList[q].drawSpoke(b[q],s[q],i);
            //viObjList[q].drawFrameOfQuadByIndex(b[q],i);
        }
    }*/






    renderWindow->Render();
    renderWindowInteractor->Start();


    return 0;
}





