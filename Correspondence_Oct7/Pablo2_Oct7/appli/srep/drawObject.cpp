/* Command this: ./drawObject -m /NIRAL/work/ltu/WorkSpace/data_pablo/figure-srep/3004_082201_3d_hippocampusLeft_c_pp-figure.m3d
 * This function show a shifted srep by given spokes' move variables (whole_coeffs).
 * whole_coeffs is gotten by combine three parts together: variables for up spokes, variables for down spokes and variables for crest spokes.
 *
 *
 *
 *
 * To call showShiftedSrep function, should set coeff first, which need set varCrest and varStand.
 *
 **/

#include "visualizeobject.h"

using namespace std;

const int rowNum = 3;
const int colNum = 13;


/* Connect three parts variables together.
 * up_vars:    size shoud be "crestAtomNums*1 + interiorAtomNums*2 - 4";    The four corners don't have variable.
 * down_vars:  size shoud be "crestAtomNums*1 + interiorAtomNums*2 - 4";    The four corners don't have variable.
 * crest_vars: size shoud be "crestAtomNums*1";                             The four corners each has one variable(for compute convenience.).
 */
vector<double> combineVariables(vector<double> up_vars, vector<double> down_vars, vector<double> crest_vars){

    vector<double> srep_coeff;

    // Connect variables
    for(unsigned int i = 0; i<up_vars.size(); i++){
        srep_coeff.push_back(up_vars[i]);
    }
    for(unsigned int i = 0; i<down_vars.size(); i++){
        srep_coeff.push_back(down_vars[i]);
    }
    for(unsigned int i = 0; i<crest_vars.size(); i++){
        srep_coeff.push_back(crest_vars[i]);
    }

    return srep_coeff;
}


/* Generate variables using given condition and inputs.
 * shift_type:  '1' means all moveable spokes can shift, vars size equals to arrayLength, its "crestAtomNums*1 - 4 + interiorAtomNums*2";
 *              '5' means only move given spokes, vars size shoud be ""; (This haven't be implement.....)
 * arrayLength: the length of the variables array to be generated. For up and down spokes (crestAtomNums*1 - 4 + interiorAtomNums*2),
 * for crest spoke (crestAtomNums).
 * moveDis: default move distance for every spokes. Only used when vars set to NULL!!!
 * vars: variables list. While shift_type set to '1': If vars set to NULL, use movdDis for all spokes;
 * If vars not NULL, use variables in vars.
 */
vector<double>  generateVariables(char shift_type, const char* vars, double moveDis, int arrayLength){
    vector<double> coeff;

    switch(shift_type){
      case '1': // Shift all moveable spokes.
        if(vars==NULL)
            cout<<"Set shift distance: "<<moveDis<<" to all moveable spokes."<<endl;
        for(unsigned int i = 0; i < arrayLength; i++) {
            coeff.push_back(moveDis);
            if(vars/*!=NULL*/){ //If given vars, use vars as shift variables
                cout<<vars[i]<<"  ";
                if(vars[i]>0 && vars[i]<1)
                    coeff[i] = vars[i];
                else {
                    cout<<"Error: Invalid vars set to this move, please check!!!"<<endl;
                    EXIT_FAILURE;
                }
            }
        }
        break;
    case '5': // Shift given spokes???? Left here, Need implement if needed.............
        break;
    //default:
        //cout<<"Invalid shift_type. Please check it in func generateVariables()!!"<<endl;
    }

    cout<<"===============the coeff is: "<<endl;
    for(unsigned int i = 0; i < coeff.size(); i++) {
        cout <<coeff[i]<<"  ";
    }
    cout<<endl;

    return coeff;
}




int main(int argc, char *argv[])
{
    vector< std::string > modelname;

    int modelnums = argc - 2;//arg[0] is exe name, arg[1] is "-m", the following is the model names.

    for(int i = 1; i <= modelnums; i++){
        modelname.push_back(argv[1+i]);
    }

    vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
    renderer->SetBackground(0,0,0);
    vtkSmartPointer<vtkRenderWindow> renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
    renderWindow->AddRenderer(renderer);
    vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
    renderWindowInteractor->SetRenderWindow(renderWindow);

    // rotate the picture, up or down. for hipcampus
    /*double rX = 0;//195;
    double rY = 0;//55;
    double rZ = 180;*/

    double rX = 0;//195;
    double rY = 0;//55;
    double rZ = 0;
    double opacity = 1;

    // rotate the picture, crest. for hipcampus
    /*double rX = 0;
    double rY = 0;
    double rZ = 180;*/


    // rotate the regularity example, Fig-5
    /*double rX = 110;
    double rY = 0;
    double rZ = 220;*/

    // Set para
    int crestAtomNum = rowNum*2 + colNum*2 - 4; // 28. crest spoke number, each spoke has one variable.
    int interiorAtomNum = rowNum*colNum - crestAtomNum; // 11
    int varStand = crestAtomNum*1 - 4 + interiorAtomNum*2; // 46. variables for standard spokes (up or down spokes).
    //int srepNum = 1;
    //int varNumEachSrep = varStand *2 + crestAtomNum; // variables for current srep.
    //int varNums = srepNum*varNumEachSrep;

    toolsfunc tf;
    int interpolationLevel = 2;
    bool moved = true;
    double quadColor[3] = {0,0.8,0};//Green.{0.3,0.3,0.3};//gray //

    vector<double> up_vars, down_vars , crest_vars; // set 0 for crest spokes.

    //Loop each of the m3d files.
    cout<<"You input "<<modelname.size()<<" models"<<endl;
    for(int q = 0; q < modelname.size(); q++){

        M3DQuadFigure* quadFig =  tf.GetQuadFigure(modelname[q].c_str());

        visualizeobject vObj(quadFig, interpolationLevel, quadColor, renderer, moved, rX, rY, rZ,opacity);
        //vObj.showBoundaryQuads(quadFig,true, true, true); //show the whole object boundary filled sub quads.

        //vObj.showBoundaryTriangulars(quadFig,true, false, true);
        //vObj.showBoundaryFrames(quadFig, true, false, true); //show the whole object boundary frame non-filled sub quads.
        //vObj.showSpokes(quadFig, true, true, true);
        vObj.showSrep(quadFig, 0, 1);
        vObj.showAtomPoints(quadFig);


        //For model q, its deltaU and deltaV are varArray[q*varNumEachSrep + j], j from 0 to varNumEachSrep.

        //const char* vars_crest = "";
        up_vars = generateVariables( '1', NULL, 0.2, varStand); // set 0.2 for up spokes.
        down_vars = generateVariables( '1', NULL, 0, varStand); // set 0 for down spokes.
        // set variables for crest spokes. 0.2 can be any value, not use if vars not NULL.
        crest_vars = generateVariables( '1', NULL, -0.3, crestAtomNum);

        vector<double> srep_coeff = combineVariables(up_vars, down_vars, crest_vars);
        /*cout<<"--------------------------------the variables for this sreps in this movement: "<<endl;
        for(unsigned int j = 0; j<varStand*2 + crestAtomNum; j++){
            cout<<srep_coeff[j]<<"  ";
        }
        cout<<endl;*/

        // Show a shifted srep
        //vObj.showShiftedSrep(srep_coeff, true, varStand, crestAtomNum, 0);

        // Show a shifted srep's crest frame.
        //vObj.showShiftedCrestFrameAndSpokes(crest_vars, crestAtomNum, 0);


        // Draw specific quads
/*        int quadNums = (quadFig->getRowCount()-1)*(quadFig->getColumnCount()-1);//16
        vtkSmartPointer< vtkPoints > hubPosition_b = vtkSmartPointer< vtkPoints >::New();
        vtkSmartPointer< vtkPoints > hubPosition_s = vtkSmartPointer< vtkPoints >::New();

        visualization visualObject(quadFig, interpolationLevel, 0, moved, quadColor, renderer, rX, rY,rZ);

        //Calculate quads areas.
        hubPosition_b = visualObject.getSubQuadsPosition(0);
        hubPosition_s = visualObject.getSubQuadsPosition(1);*/


           //cout<<"Drawing quad: "<<i<<endl;
            //quadColor[0] =1;
            //quadColor[1] =0;
            //quadColor[2] =0.4+0.1*i; //fen se





/*                quadColor[0] =0.9;
                quadColor[1] =0;        //red
                quadColor[2] =0;
                visualObject.setColor(quadColor);
                //visualObject.drawQuadByIndexAndHubpos(hubPosition_b,1);
                //visualObject.drawQuadByIndexAndHubpos(hubPosition_b,13);
                //visualObject.drawQuadByIndexAndHubpos(hubPosition_b,11);
                //visualObject.drawQuadByIndexAndHubpos(hubPosition_b,7);

                visualObject.drawVolume(hubPosition_b,hubPosition_s,1);
                visualObject.drawVolume(hubPosition_b,hubPosition_s,13);
                visualObject.drawVolume(hubPosition_b,hubPosition_s,11);
                visualObject.drawVolume(hubPosition_b,hubPosition_s,7);


                quadColor[0] =0.1;
                quadColor[1] =0;        //magenta
                quadColor[2] =0;
                visualObject.setColor(quadColor);
                visualObject.drawFrameOfQuadByIndex(hubPosition_s,1);
                visualObject.drawFrameOfQuadByIndex(hubPosition_s,13);
                visualObject.drawFrameOfQuadByIndex(hubPosition_s,11);
                visualObject.drawFrameOfQuadByIndex(hubPosition_s,7);

                visualObject.drawSpoke(hubPosition_b,hubPosition_s,1);
                visualObject.drawSpoke(hubPosition_b,hubPosition_s,13);
                visualObject.drawSpoke(hubPosition_b,hubPosition_s,11);
                visualObject.drawSpoke(hubPosition_b,hubPosition_s,7);

                visualObject.drawFrameOfQuadByIndex(hubPosition_b,1);
                visualObject.drawFrameOfQuadByIndex(hubPosition_b,13);
                visualObject.drawFrameOfQuadByIndex(hubPosition_b,11);
                visualObject.drawFrameOfQuadByIndex(hubPosition_b,7);*/

//                quadColor[0] =0;
//                quadColor[1] =0;        //(0,0.2,0) dark green.
//                quadColor[2] =0;
//                visualObject.setColor(quadColor);
                // Only draw the quads (without the subquads, no matter what's the interpolation level)
                //visualObject.drawEdgeOfQuads_method3(hubPosition_b);

                /*quadColor[0] =1;
                quadColor[1] =0;  //(1,0,0) red //(0,0.2,0) dark green.(0,0.8,0) bright green.
                quadColor[2] =0;
                visualObject.setColor(quadColor);
                //visualObject.drawFrameOfQuadByIndex(hubPosition_s,i);
                visualObject.drawQuadByIndexAndHubpos(hubPosition_s,i);
                quadColor[0] =0;
                quadColor[1] =0;        //(0,0.2,0) dark green.
                quadColor[2] =1;
                visualObject.setColor(quadColor);*/
                // Only draw the quads (without the subquads, no matter what's the interpolation level)
                //visualObject.drawEdgeOfQuads_method3(hubPosition_s);

      // }

/*      visualization visualObject_down(quadFig, interpolationLevel, 1, moved, quadColor, renderer, rX, rY,rZ);
      hubPosition_b = visualObject_down.getSubQuadsPosition(0);
      hubPosition_s = visualObject_down.getSubQuadsPosition(1);

/*      quadColor[0] =0;
      quadColor[1] =0.7; // green
      quadColor[2] =0;
      visualObject.setColor(quadColor);
      visualObject.drawVolume(hubPosition_b,hubPosition_s,4);
      visualObject.drawVolume(hubPosition_b,hubPosition_s,11);
      visualObject.drawVolume(hubPosition_b,hubPosition_s,1);
      visualObject.drawVolume(hubPosition_b,hubPosition_s,6);

      //visualObject.drawQuadByIndexAndHubpos(hubPosition_b,4);
      //visualObject.drawQuadByIndexAndHubpos(hubPosition_b,11);
      //visualObject.drawQuadByIndexAndHubpos(hubPosition_b,1);
      //visualObject.drawQuadByIndexAndHubpos(hubPosition_b,6);


    quadColor[0] =0;
     quadColor[1] =0.1;        //(0,0.2,0) dark green.
     quadColor[2] =0;
     visualObject.setColor(quadColor);
     visualObject.drawSpoke(hubPosition_b,hubPosition_s,4);
     visualObject.drawSpoke(hubPosition_b,hubPosition_s,11);
     visualObject.drawSpoke(hubPosition_b,hubPosition_s,1);
    visualObject.drawSpoke(hubPosition_b,hubPosition_s,6);

      visualObject.drawFrameOfQuadByIndex(hubPosition_b,4);
     visualObject.drawFrameOfQuadByIndex(hubPosition_b,11);
      visualObject.drawFrameOfQuadByIndex(hubPosition_b,1);
     visualObject.drawFrameOfQuadByIndex(hubPosition_b,6);

      visualObject.drawFrameOfQuadByIndex(hubPosition_s,4);
            visualObject.drawFrameOfQuadByIndex(hubPosition_s,11);
            visualObject.drawFrameOfQuadByIndex(hubPosition_s,1);
           visualObject.drawFrameOfQuadByIndex(hubPosition_s,6);*/


/*      quadColor[0] =0;
      quadColor[1] =0.7;        //(0,0.2,0) dark green.
      quadColor[2] =0;
      visualObject.setColor(quadColor);
      visualObject.drawSpoke(hubPosition_b,hubPosition_s,4);
      quadColor[0] =1;
      quadColor[1] =1;        // light yellow
      quadColor[2] =0;
      visualObject.setColor(quadColor);
      visualObject.drawFrameOfQuadByIndex(hubPosition_b,4);

      quadColor[0] =0;
      quadColor[1] =1;        //(0,0.2,0) dark green.
      quadColor[2] =1;
      visualObject.setColor(quadColor);
      visualObject.drawFrameOfQuadByIndex(hubPosition_s,4);
      quadColor[0] =0;
      quadColor[1] =0.7; // green
      quadColor[2] =0;
      visualObject.setColor(quadColor);
      visualObject.drawVolume(hubPosition_b,hubPosition_s,4);



      quadColor[0] =0.1;
      quadColor[1] =0.1;  //(1,0,0) red //(0,0.2,0) dark green.(0,0.8,0) bright green.
      quadColor[2] =1;
      visualObject.setColor(quadColor);
      visualObject.drawSubQuadVolume(hubPosition_b,hubPosition_s,4,28);
*/


        // Visualize the crest part.
/*            quadColor[0] =1;
            quadColor[1] =0.1;  //(1,0,0) red //(0,0.2,0) dark green.(0,0.8,0) bright green.
            quadColor[2] =1;
            visualObject.setColor(quadColor);
            visualObject.drawFrameOfQuadByIndex(hubPosition_s,7);
            visualObject.drawSpoke(hubPosition_b,hubPosition_s,7);
            visualObject.drawFrameOfQuadByIndex(hubPosition_b,7);

            quadColor[0] =0;
            quadColor[1] =0.9;  //(1,0,0) red //(0,0.2,0) dark green.(0,0.8,0) bright green.
            quadColor[2] =0;
            visualObject.setColor(quadColor);
            visualObject.drawVolume(hubPosition_b,hubPosition_s,7);
            //visualObject.drawQuadByIndexAndHubpos(hubPosition_s,i);
*/

        // Draw crest frame of quad 11
/*       quadColor[0] =0;
              quadColor[1] =1;        //green
              quadColor[2] =0;
              //vObj.showCrestFrameByQuadIndex(quadFig,11, interpolationLevel, quadColor);
              vObj.showCrestSpokesByQuadIndex(quadFig,11, interpolationLevel, quadColor);
              quadColor[0] =0.9;
              quadColor[1] =0;        //dark red
              quadColor[2] =0;
              //vObj.showCrestFrameByQuadIndex(quadFig,5, interpolationLevel, quadColor);
              vObj.showCrestSpokesByQuadIndex(quadFig,5, interpolationLevel, quadColor);
              quadColor[0] =0.2;
              quadColor[1] =0.7;        //light blue
              quadColor[2] =0.9;
              //vObj.showCrestFrameByQuadIndex(quadFig,3, interpolationLevel, quadColor);
              vObj.showCrestSpokesByQuadIndex(quadFig,3, interpolationLevel, quadColor);
              quadColor[0] =1;
              quadColor[1] =0.9;        //dark yellow
              quadColor[2] =0;
              //vObj.showCrestFrameByQuadIndex(quadFig,2, interpolationLevel, quadColor);
              vObj.showCrestSpokesByQuadIndex(quadFig,2, interpolationLevel, quadColor);
              quadColor[0] =0.9;
              quadColor[1] =0.5;        //dark red
              quadColor[2] =0.1;
              //vObj.showCrestFrameByQuadIndex(quadFig,14, interpolationLevel, quadColor);
              vObj.showCrestSpokesByQuadIndex(quadFig,14, interpolationLevel, quadColor);
              quadColor[0] =0.9;
              quadColor[1] =0;        //(0,0.2,0) dark green.
              quadColor[2] =0.9;
             // vObj.showCrestFrameByQuadIndex(quadFig,16, interpolationLevel, quadColor);
              vObj.showCrestSpokesByQuadIndex(quadFig,16, interpolationLevel, quadColor);
              quadColor[0] =0;
              quadColor[1] =0.3;        //light yellow.
              quadColor[2] =1;
              //vObj.showCrestFrameByQuadIndex(quadFig,9, interpolationLevel, quadColor);
              vObj.showCrestSpokesByQuadIndex(quadFig,9, interpolationLevel, quadColor);
              quadColor[0] =0;
              quadColor[1] =1;        //cyan.
              quadColor[2] =1;
              //vObj.showCrestFrameByQuadIndex(quadFig,7, interpolationLevel, quadColor);
              vObj.showCrestSpokesByQuadIndex(quadFig,7, interpolationLevel, quadColor);

        quadColor[0] =0;
        quadColor[1] =0.8;    //(1,0.4,0) dark yellow(orange)  //(1.0.6,0) light orange  //(0,0.2,0) dark green.
        quadColor[2] =0;
        vObj.showCrestQuadByQuadIndex(quadFig,11, interpolationLevel, quadColor);
        vObj.showCrestQuadByQuadIndex(quadFig,5, interpolationLevel, quadColor);
        vObj.showCrestQuadByQuadIndex(quadFig,3, interpolationLevel, quadColor);
        vObj.showCrestQuadByQuadIndex(quadFig,2, interpolationLevel, quadColor);
        vObj.showCrestQuadByQuadIndex(quadFig,14, interpolationLevel, quadColor);
        vObj.showCrestQuadByQuadIndex(quadFig,16, interpolationLevel, quadColor);
        vObj.showCrestQuadByQuadIndex(quadFig,9, interpolationLevel, quadColor);
        vObj.showCrestQuadByQuadIndex(quadFig,7, interpolationLevel, quadColor);*/

        // Draw the 5th, bigger figure.
/*        quadColor[0] =0.9;
        quadColor[1] =0;        //dark red
        quadColor[2] =0;
        vObj.showCrestSpokesByQuadIndex(quadFig,5, interpolationLevel, quadColor);
        quadColor[0] =0;
        quadColor[1] =0.8;    //(1,0.4,0) dark yellow(orange)  //(1.0.6,0) light orange  //(0,0.2,0) dark green.
        quadColor[2] =0;
        vObj.showCrestVolumesByQuadIndex(quadFig,5, interpolationLevel, quadColor);
        quadColor[0] =0;
        quadColor[1] =0.8;    //(1,0.4,0) dark yellow(orange)  //(1.0.6,0) light orange  //(0,0.2,0) dark green.
        quadColor[2] =0;
vObj.showCrestQuadByQuadIndex(quadFig,5, interpolationLevel, quadColor);
        quadColor[0] =0;
        quadColor[1] =0;
        quadColor[2] =0;
        vObj.showCrestFrameByQuadIndex(quadFig,5, interpolationLevel, quadColor);
        //vObj.showCrestFrameByQuadIndex(quadFig,16, interpolationLevel, quadColor);
        //vObj.showCrestFrameByQuadIndex(quadFig,7, interpolationLevel, quadColor);*/
//vObj.showSrep(quadFig, 0);

       /* quadColor[0] =1;
        quadColor[1] =0.4;    //(1,0.4,0) dark yellow(orange)  //(1.0.6,0) light orange  //(0,0.2,0) dark green.
        quadColor[2] =0;
        vObj.showCrestSpokesByQuadIndex(quadFig,11, interpolationLevel, quadColor);
        vObj.showCrestSpokesByQuadIndex(quadFig,5, interpolationLevel, quadColor);
        vObj.showCrestSpokesByQuadIndex(quadFig,3, interpolationLevel, quadColor);
        vObj.showCrestSpokesByQuadIndex(quadFig,2, interpolationLevel, quadColor);
        vObj.showCrestSpokesByQuadIndex(quadFig,14, interpolationLevel, quadColor);
        vObj.showCrestSpokesByQuadIndex(quadFig,16, interpolationLevel, quadColor);
        vObj.showCrestSpokesByQuadIndex(quadFig,9, interpolationLevel, quadColor);*/




/*        quadColor[0] =1;
        quadColor[1] =0.4;    //(1,0.4,0) dark yellow(orange)  //(1.0.6,0) light orange  //(0,0.2,0) dark green.
        quadColor[2] =0;
        vObj.showCrestSpokesByQuadIndex(quadFig,5, interpolationLevel, quadColor);
        quadColor[0] =0;
        quadColor[1] =0.9;        //(0,0.2,0) dark green.
        quadColor[2] =0;
        //vObj.showCrestVolumesByQuadIndex(quadFig,5, interpolationLevel, quadColor);
        quadColor[0] =0;
        quadColor[1] =0.9;        //(0,0.2,0) dark green.
        quadColor[2] =0;
        vObj.showCrestFrameByQuadIndex(quadFig,5, interpolationLevel, quadColor);*/


/*              for (unsigned int i =0; i<quadNums; i++){
                  quadColor[0] =0;
                  quadColor[1] =0.7;  //(1,0,0) red //(0,0.2,0) dark green.(0,0.8,0) bright green.
                  quadColor[2] =0;
                  visualObject.setColor(quadColor);
                visualObject.drawSpoke(hubPosition_b,hubPosition_s,i);
                quadColor[0] =0;
                quadColor[1] =0.9;  //(1,0,0) red //(0,0.2,0) dark green.(0,0.8,0) bright green.
                quadColor[2] =0;
                visualObject.setColor(quadColor);
                //visualObject.drawVolume(hubPosition_b,hubPosition_s,i);
                quadColor[0] =0.1;
                quadColor[1] =1;  //(1,0,0) red //(0,0.2,0) dark green.(0,0.8,0) bright green.
                quadColor[2] =1;
                visualObject.setColor(quadColor);
                visualObject.drawFrameOfQuadByIndex(hubPosition_s,i);
                quadColor[0] =1;
                quadColor[1] =0.1;  //(1,0,0) red //(0,0.2,0) dark green.(0,0.8,0) bright green.
                quadColor[2] =1;
                visualObject.setColor(quadColor);
                visualObject.drawFrameOfQuadByIndex(hubPosition_b,i);
              }
*/
/*              quadColor[0] =1;
              quadColor[1] =1;  //(1,0,0) red //(0,0.2,0) dark green.(0,0.8,0) bright green.
              quadColor[2] =0.01;
              visualObject.setColor(quadColor);
              visualObject.drawSubQuadVolume(hubPosition_b,hubPosition_s,0,10);

              quadColor[0] =0;
              quadColor[1] =0;  //(1,0,0) red //(0,0.2,0) dark green.(0,0.8,0) bright green.
              quadColor[2] =1;
              visualObject.setColor(quadColor);
              visualObject.drawSubQuadVolume(hubPosition_b,hubPosition_s,0,42);

              quadColor[0] =1;
              quadColor[1] =0.1;  //(1,0,0) red //(0,0.2,0) dark green.(0,0.8,0) bright green.
              quadColor[2] =1;
              visualObject.setColor(quadColor);
              visualObject.drawSubQuadVolume(hubPosition_b,hubPosition_s,0,30);
*/


    }

    renderWindow->Render();
    renderWindowInteractor->Start();

    return 0;
}



