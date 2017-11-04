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
    renderer->SetBackground(1,1,1);
    vtkSmartPointer<vtkRenderWindow> renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
    renderWindow->AddRenderer(renderer);
    vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
    renderWindowInteractor->SetRenderWindow(renderWindow);

    // rotate the picture, crest, for Fig 13, overlap sreps.
    double rX = 45;//195;
    double rY = -40;//55;
    double rZ = -90;
    double opacity = 1;

    // rotate the picture, crest, for Fig 14, only show spokes in different color,without srep. view 1.
    /*double rX = 90;//195;
    double rY = 0;//55;
    double rZ = 0;*/

    // rotate the picture, crest, for Fig 14, only show spokes in different color,without srep. view 2.
    /*double rX = 0;
    double rY = 85;
    double rZ = 0;*/


    // rotate the picture, crest, for Fig 14, only show spokes in different color,without srep. view 3.
    /*double rX = 90;
    double rY = 0;
    double rZ = 180;*/

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

    double minz = 0;
    double maxz = 80;
    // Create the color map
    vtkSmartPointer<vtkLookupTable> colorLookupTable = vtkSmartPointer<vtkLookupTable>::New();
    colorLookupTable->SetTableRange(minz, maxz);
    colorLookupTable->Build();



    M3DQuadFigure* quadFig =  tf.GetQuadFigure(modelname[0].c_str());
    visualizeobject vObj_initi(quadFig, interpolationLevel, quadColor, renderer, moved, rX, rY, rZ, opacity);
    vObj_initi.showSrep(quadFig, 0, 1);

    //Loop each of the m3d files.
    cout<<"You input "<<modelname.size()<<" models"<<endl;
    for(int q = 0; q < modelname.size(); q++){

        M3DQuadFigure* quadFig =  tf.GetQuadFigure(modelname[q].c_str());

        visualizeobject vObj(quadFig, interpolationLevel, quadColor, renderer, moved, rX, rY, rZ, opacity);
        //vObj.showBoundaryQuads(quadFig,true, true, true); //show the whole object boundary filled sub quads.

        //vObj.showBoundaryTriangulars(quadFig,true, false, true);
        //vObj.showBoundaryFrames(quadFig, true, false, true); //show the whole object boundary frame non-filled sub quads.
        //vObj.showSpokes(quadFig, true, true, true);
        //vObj.showSrep(quadFig, 0);
        //vObj.showAtomPoints(quadFig);



        /*cout<<"--------------------------------the variables for this sreps in this movement: "<<endl;
        for(unsigned int j = 0; j<varStand*2 + crestAtomNum; j++){
            cout<<srep_coeff[j]<<"  ";
        }
        cout<<endl;*/

        // Show a shifted srep
        //vObj.showShiftedSrep(srep_coeff, true, varStand, crestAtomNum, 0);

        // Show a shifted srep's crest frame.
        //vObj.showShiftedCrestFrameAndSpokes(crest_vars, crestAtomNum, 0);





        // Get a color correspondence to the z coordinate value from the color table
        double dcolor[3];
        colorLookupTable->GetColor(q, dcolor);

        // Draw crest frame of quad 11
        //quadColor[0] =1*q/20 + 0; //0+q*0.05;
        //quadColor[1] =0.5+q/20*0.5;    //(1,0.4,0) dark yellow(orange)  //(1.0.6,0) light orange  //(0,0.2,0) dark green.
        //quadColor[2] =0+q/20;
           //for(unsigned int i=0; i<28; i++){
        vObj.showCrestSpokesByQuadIndex(quadFig,17, interpolationLevel, dcolor);
                /*quadColor[0] =1;
                quadColor[1] =0.6;    //(1,0.4,0) dark yellow(orange)  //(1.0.6,0) light orange  //(0,0.2,0) dark green.
                quadColor[2] =0;
                   //for(unsigned int i=0; i<28; i++){
                        vObj.showCrestSpokesByQuadIndex(quadFig,16, interpolationLevel, quadColor);
           // }

        /*quadColor[0] =0;
        quadColor[1] =0.8;    //(1,0.4,0) dark yellow(orange)  //(1.0.6,0) light orange  //(0,0.2,0) dark green.
        quadColor[2] =0;
        vObj.showCrestQuadByQuadIndex(quadFig,17, interpolationLevel, quadColor);
        quadColor[0] =1;
        quadColor[1] =0.2;    //(1,0.4,0) dark yellow(orange)  //(1.0.6,0) light orange  //(0,0.2,0) dark green.
        quadColor[2] =0;
        vObj.showCrestQuadByQuadIndex(quadFig,16, interpolationLevel, quadColor);
        /*vObj.showCrestQuadByQuadIndex(quadFig,3, interpolationLevel, quadColor);
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





    }

    renderWindow->Render();
    renderWindowInteractor->Start();

    return 0;
}




