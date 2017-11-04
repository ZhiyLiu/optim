/* This class draw a whole object, which use old class visualization(which only can draw up or down side).
 * Based on visualization, here only add crest quads in to render.
 *
 * Liyun Tu
 * Apr 17, 2014
*/


#include "visualizeobject.h"

visualizeobject::visualizeobject()
{
}

visualizeobject::visualizeobject(M3DQuadFigure* quadfig, int interpolationLevel, double *quadColor, vtkSmartPointer<vtkRenderer> renderer,
                                 bool moved, double rX, double rY, double rZ, double opacity){
    this->quadFig = quadfig;
    this->interpolationLevel = interpolationLevel;
    this->quadColor = quadColor;
    this->renderer = renderer;
    this->step = pow((double)2, (double)interpolationLevel);
    this->moved = true;
    this->rowNums = this->quadFig->getRowCount();
    this->colNums = this->quadFig->getColumnCount();
    this->quadNums = (this->rowNums-1)*(this->colNums-1);
    this->rX = rX;
    this->rY = rY;
    this->rZ = rZ;
    this->opacity = opacity;
}


/* Draw the boundary surface of object in quads.
 * This function draw the object surface by non-filled subquads. Only shows the quad frame of the boundary surface.
 * We can choose to show which side..
 * up(/down/crest): if true, shows up(/down/crest) boundary quads.
*/
void visualizeobject::showBoundaryFrames(M3DQuadFigure* quadFig, bool up, bool down, bool crest){

    vtkSmartPointer< vtkPoints > hubPosition_b = vtkSmartPointer< vtkPoints >::New();

    // Add up boundary quads to render
    if(up){ // Show up quads
        //double quadColor[3]={0.2,0.2,0.2};
        visualization visualObject_up(quadFig, this->interpolationLevel, 0, this->moved, quadColor, this->renderer,this->rX, this->rY,this->rZ,
                                      this->opacity);

        hubPosition_b = visualObject_up.getSubQuadsPosition(0); //0(boundary), 1(skeletal)        

        /*double quadColor[3] = {1,0.5,0};
        for(int i =0; i<24; i++){
            quadColor[0] =1;//0.05+i*0.06;
            quadColor[1] =0.2+i*0.1;
            quadColor[2] =0;
            visualObject_up.setColor(quadColor);
            visualObject_up.drawQuadByIndexAndHubpos(hubPosition_b,i);
        }*/

        // Draw up boundary frame.
        visualObject_up.drawFrameOfQuads(hubPosition_b);
    }

    // Add down boundary quads to render
    if(down){ // Show down quads
        //quadColor[0] = 0.8; quadColor[1] = 0; quadColor[2] = 0;
        //double quadColor[3]={0.2,0.2,0.2};
        visualization visualObject_down(quadFig, this->interpolationLevel, 1, this->moved, quadColor, this->renderer, this->rX,
                                        this->rY, this->rZ, this->opacity);

        hubPosition_b = visualObject_down.getSubQuadsPosition(0); //0(boundary), 1(skeletal)

        /*double quadColor[3] = {1,0.5,0};
        for(int i =0; i<24; i++){
            quadColor[0] =0;//0.05+i*0.06;
            quadColor[1] =0.2+i*0.1;
            quadColor[2] =1;
            visualObject_down.setColor(quadColor);
            visualObject_down.drawQuadByIndexAndHubpos(hubPosition_b,i);
        }*/

        // Draw down boundary frame.
        visualObject_down.drawFrameOfQuads(hubPosition_b);
    }

    // Add crest boundary quads to render
    if(crest){ // Show crest quads.
        //double quadColor[3] = {0.1,0.1,0.1}; //0.5,0.2,0.9
        //double quadColor[3]={0.2,0.2,0.2};
        visualizecrest visualObject_crest(quadFig, this->interpolationLevel, quadColor, this->renderer, this->rX, this->rY, this->rZ, this->opacity);
        visualObject_crest.drawCrestFrames();
    }
}


void visualizeobject::showSpokes(M3DQuadFigure* quadFig, bool up, bool down, bool crest){

    vtkSmartPointer< vtkPoints > hubPosition_b = vtkSmartPointer< vtkPoints >::New();
    vtkSmartPointer< vtkPoints > hubPosition_s = vtkSmartPointer< vtkPoints >::New();

    // Add up boundary quads to render
    if(up){ // Show up quads
        //double quadColor[3]={0.2,0.2,0.2};
        visualization visualObject_up(quadFig, this->interpolationLevel, 0, this->moved, quadColor, this->renderer, this->rX,
                                      this->rY, this->rZ, this->opacity);
        cout<<"-----------interpolationLevel: "<<interpolationLevel<<endl;
        hubPosition_b = visualObject_up.getSubQuadsPosition(0); //0(boundary), 1(skeletal)
        hubPosition_s = visualObject_up.getSubQuadsPosition(1); //0(boundary), 1(skeletal)

        /*double quadColor[3] = {1,0.5,0};
        for(int i =0; i<24; i++){
            quadColor[0] =1;//0.05+i*0.06;
            quadColor[1] =0.2+i*0.1;
            quadColor[2] =0;
            visualObject_up.setColor(quadColor);
            visualObject_up.drawQuadByIndexAndHubpos(hubPosition_b,i);
        }*/

        double quadColor1[3]={0.1,1,1};//{0.1,0.1,0.1};;//
        visualObject_up.setColor(quadColor1);
        for(int i =0; i<24; i++){
        visualObject_up.drawSpoke(hubPosition_s, hubPosition_b, i);
        }

        // Draw up boundary frame.
        //visualObject_up.drawFrameOfQuads(hubPosition_b);
    }

    // Add down boundary quads to render
    if(down){ // Show down quads
        //quadColor[0] = 0.8; quadColor[1] = 0; quadColor[2] = 0;
        //double quadColor[3]={0.2,0.2,0.2};;
        visualization visualObject_down(quadFig, this->interpolationLevel, 1, this->moved, quadColor, this->renderer, this->rX,
                                        this->rY, this->rZ, this->opacity);

        hubPosition_b = visualObject_down.getSubQuadsPosition(0); //0(boundary), 1(skeletal)
        hubPosition_s = visualObject_down.getSubQuadsPosition(1); //0(boundary), 1(skeletal)

        /*double quadColor[3] = {1,0.5,0};
        for(int i =0; i<24; i++){
            quadColor[0] =0;//0.05+i*0.06;
            quadColor[1] =0.2+i*0.1;
            quadColor[2] =1;
            visualObject_down.setColor(quadColor);
            visualObject_down.drawQuadByIndexAndHubpos(hubPosition_b,i);
        }*/

        double quadColor1[3]={1,0.1,1};//{0.1,0.1,0.1};
        visualObject_down.setColor(quadColor1);
        for(int i =0; i<24; i++){
        visualObject_down.drawSpoke(hubPosition_s, hubPosition_b, i);
        }

        // Draw down boundary frame.
        //visualObject_down.drawFrameOfQuads(hubPosition_b);
    }

    // Add crest boundary quads to render
    if(crest){ // Show crest quads.
        //double quadColor[3] = {0.1,0.1,0.1}; //0.5,0.2,0.9
        double quadColor[3]={1,1,0.1};//{0.1,0.1,0.1};//
        visualizecrest visualObject_crest(quadFig, this->interpolationLevel, quadColor, this->renderer, this->rX, this->rY, this->rZ, this->opacity);
        //visualObject_crest.drawCrestFrames();
        visualObject_crest.drawCrestSpokes();
    }
}


void visualizeobject::showSpokes_2(M3DQuadFigure* quadFig, bool up, bool down, bool crest, double* color, double opacity){

    vtkSmartPointer< vtkPoints > hubPosition_b = vtkSmartPointer< vtkPoints >::New();
    vtkSmartPointer< vtkPoints > hubPosition_s = vtkSmartPointer< vtkPoints >::New();

    // Add up boundary quads to render
    if(up){ // Show up quads
        //double quadColor[3]={0.2,0.2,0.2};
        visualization visualObject_up(quadFig, this->interpolationLevel, 0, this->moved, quadColor, this->renderer, this->rX,
                                      this->rY, this->rZ, opacity);
        cout<<"-----------interpolationLevel: "<<interpolationLevel<<endl;
        hubPosition_b = visualObject_up.getSubQuadsPosition(0); //0(boundary), 1(skeletal)
        hubPosition_s = visualObject_up.getSubQuadsPosition(1); //0(boundary), 1(skeletal)

        visualObject_up.setColor(color);
        for(int i =0; i<24; i++){
        visualObject_up.drawSpoke(hubPosition_s, hubPosition_b, i);
        }
    }

    // Add down boundary quads to render
    if(down){ // Show down quads
        visualization visualObject_down(quadFig, this->interpolationLevel, 1, this->moved, quadColor, this->renderer, this->rX,
                                        this->rY, this->rZ, opacity);

        hubPosition_b = visualObject_down.getSubQuadsPosition(0); //0(boundary), 1(skeletal)
        hubPosition_s = visualObject_down.getSubQuadsPosition(1); //0(boundary), 1(skeletal)

        visualObject_down.setColor(color);
        for(int i =0; i<24; i++){
        visualObject_down.drawSpoke(hubPosition_s, hubPosition_b, i);
        }
    }

    // Add crest boundary quads to render
    if(crest){ // Show crest quads.
        visualizecrest visualObject_crest(quadFig, this->interpolationLevel, color, this->renderer, this->rX, this->rY, this->rZ, opacity);
        //visualObject_crest.drawCrestFrames();
        visualObject_crest.drawCrestSpokes();
    }
}


/* Draw the boundary surface of object in triangular.
 *
*/
void visualizeobject::showBoundaryTriangulars(M3DQuadFigure* quadFig, bool up, bool down, bool crest){

    vtkSmartPointer< vtkPoints > hubPosition_b = vtkSmartPointer< vtkPoints >::New();
    double quadColor[3]={0.3,0.3,0.3};

    // Add up boundary quads to render
    if(up){ // Show up quads
        /*quadColor[0] = 1;
        quadColor[1] = 0.1;
        quadColor[2] = 1; //set up quads green.*/
        visualization visualObject_up(quadFig, this->interpolationLevel, 0, this->moved, quadColor, this->renderer, this->rX,
                                      this->rY, this->rZ, this->opacity);

        hubPosition_b = visualObject_up.getSubQuadsPosition(0); //0(boundary), 1(skeletal)

        // Draw up boundary frame.
        //visualObject_up.drawQuads(hubPosition_b);

        visualObject_up.connectDiagnoal(hubPosition_b);
    }

    // Add down boundary quads to render
    if(down){ // Show down quads
       /* quadColor[0] = 0.1;
        quadColor[1] = 1;
        quadColor[2] = 1; // set color*/
        visualization visualObject_down(quadFig, this->interpolationLevel, 1, this->moved, quadColor, this->renderer, this->rX,
                                        this->rY, this->rZ, this->opacity);

        hubPosition_b = visualObject_down.getSubQuadsPosition(0); //0(boundary), 1(skeletal)

        // Draw down boundary frame.
        //visualObject_down.drawQuads(hubPosition_b);
        visualObject_down.connectDiagnoal(hubPosition_b);
    }

    // Add crest boundary quads to render
    if(crest){ // Show crest quads.
        /*quadColor[0] = 1;
        quadColor[1] = 1;
        quadColor[2] = 0.1; // set color, purpul*/
        visualizecrest visualObject_crest(quadFig, this->interpolationLevel, quadColor, this->renderer, this->rX, this->rY, this->rZ, this->opacity);
        //visualObject_crest.drawCrestQuads();
        visualObject_crest.connectDiagnoal();
    }
}



/* Draw the boundary surface of object in color (smoothly).
 * This function draw the object surface by filled subquads gotten from higher level interpolation, without show the quad frame.
 *
*/
void visualizeobject::showBoundaryQuads(M3DQuadFigure* quadFig, bool up, bool down, bool crest){

    vtkSmartPointer< vtkPoints > hubPosition_b = vtkSmartPointer< vtkPoints >::New();
    double quadColor[3];

    // Add up boundary quads to render
    if(up){ // Show up quads
        quadColor[0] = this->quadColor[0] + 0;
        quadColor[1] = this->quadColor[0] + 1;
        quadColor[2] = this->quadColor[0] + 0; //set up quads green.
        visualization visualObject_up(quadFig, this->interpolationLevel, 0, this->moved, quadColor, this->renderer, this->rX,
                                      this->rY, this->rZ, this->opacity);

        hubPosition_b = visualObject_up.getSubQuadsPosition(0); //0(boundary), 1(skeletal)

        // Draw up boundary frame.
        visualObject_up.drawQuads(hubPosition_b);
    }

    // Add down boundary quads to render
    if(down){ // Show down quads
        quadColor[0] = this->quadColor[0] + 0.8;
        quadColor[1] = this->quadColor[0] + 0;
        quadColor[2] = this->quadColor[0] + 0; // set color
        visualization visualObject_down(quadFig, this->interpolationLevel, 1, this->moved, quadColor, this->renderer, this->rX,
                                        this->rY, this->rZ, this->opacity);

        hubPosition_b = visualObject_down.getSubQuadsPosition(0); //0(boundary), 1(skeletal)

        // Draw down boundary frame.
        visualObject_down.drawQuads(hubPosition_b);
    }

    // Add crest boundary quads to render
    if(crest){ // Show crest quads.
        quadColor[0] = this->quadColor[0] + 0;
        quadColor[1] = this->quadColor[0] + 0;
        quadColor[2] = this->quadColor[0] + 0.7; // set color
        visualizecrest visualObject_crest(quadFig, this->interpolationLevel, quadColor, this->renderer, this->rX, this->rY, this->rZ, this->opacity);
        visualObject_crest.drawCrestQuads();
    }
}



/* Draw a srep's up or down skeletal sheet and spokes.
 * side: 0(up), 1(down).
 *
*/
void visualizeobject::showStandardSrep(M3DQuadFigure* quadFig, int side, double *quadColor, int interpolationLevel, int linewith, bool showspoketail){

    vtkSmartPointer< vtkPoints > hubPosition_s = vtkSmartPointer< vtkPoints >::New();
    vtkSmartPointer< vtkPoints > hubPosition_b = vtkSmartPointer< vtkPoints >::New();

    visualization visualObject(quadFig, interpolationLevel, side, this->moved, this->quadColor, this->renderer, rX, rY, rZ, this->opacity);

    hubPosition_s = visualObject.getSubQuadsPosition(1); //0(boundary), 1(skeletal)
    hubPosition_b = visualObject.getSubQuadsPosition(0); //0(boundary), 1(skeletal)

    visualObject.setColor(quadColor);

    double quadColor1[3] = {0,0.6,0};
    visualObject.setColor(quadColor1);
    //visualObject.drawQuads(hubPosition_s);

    // Add skeletal sheet frame to render
    visualObject.drawFrameOfQuads(hubPosition_s);


    // Set up spoke to yellow, down spoke to Cyan
    double color[3] = {0.1,1,1}; //Cyan //{0,0,1};//Megneta
    if(side==0)  { // {1,1,0}Yellow
        color[0] = 1;
        color[1] = 0.1;
        color[2] = 1;
    }
    visualObject.setColor(color);

    // Add spokes to render    
    visualObject.drawSpoke_2(hubPosition_b, hubPosition_s, this->quadNums, linewith, showspoketail);
    // if want to show spokes with different color. But there are spokes repeated, not reconmanded to use.
    //visualObject.drawSpoke_3(hubPosition_b, hubPosition_s, this->quadNums);

}


/* Draw a srep's crest spokes.*/
void visualizeobject::showCrestSrep(M3DQuadFigure* quadFig, double * quadColor, int interpolationLevel, int linewidth, bool showspoketail){

    // Draw the 0 level, in default.
    //visualizecrest visualObject_crest(quadFig, 0, quadColor, this->renderer, this->rX, this->rY, this->rZ);
    // Draw the true level.

    visualizecrest visualObject_crest(quadFig, 0, quadColor, this->renderer, this->rX, this->rY, this->rZ, this->opacity);

    // Add crest line on the skeletal to render
    quadColor[0]=0;
    quadColor[0]=0.6;
    quadColor[0]=0;
    //visualObject_crest.drawCrestCurve(1,quadColor);

    // Add crest spokes to render
    quadColor[0] = 1;
    quadColor[1] = 1;
    quadColor[2] = 0.1;
    visualObject_crest.drawMidOneCrestSpokes(rX, rY, rZ, linewidth, quadColor, showspoketail);

    // Add the crest line on the boundary to render
    //visualObject_crest.drawCrestBoundaryEquator(rX, rY, rZ);
}



void visualizeobject::drawSrep(M3DQuadFigure* quadFig, bool showup, bool showdown, bool showcrest, double *quadColor, int linewith){

    int rowNum = quadFig->getRowCount();
    int colNum = quadFig->getColumnCount();

    M3DQuadPrimitive* prim;
    M3DQuadEndPrimitive* endPrim;

    Vector3D spoketail, spoketip;

    vtkSmartPointer< vtkPoints > hubpos = vtkSmartPointer< vtkPoints >::New();
    vtkSmartPointer<vtkCellArray> lineCellArray = vtkSmartPointer<vtkCellArray>::New();
    vtkSmartPointer<vtkPolyData> linePolyData = vtkSmartPointer<vtkPolyData>::New();

    // Loop each atom of the quadfig
    for(unsigned u = 0; u < rowNum; u++){
        for(unsigned v = 0; v < colNum; v++){

            prim = dynamic_cast<M3DQuadPrimitive*>(quadFig->getPrimitivePtr(u,v));
            spoketail = prim->getX();

            // Add up spokes
            if(showup) {
                // Get up spoke tail and tip
                spoketip = prim->getX() + prim->getR0()*prim->getU0();

                vtkSmartPointer<vtkLine> medialsheetline = vtkSmartPointer<vtkLine>::New();

                vtkIdType id0 = hubpos->InsertNextPoint(spoketail.getX(), spoketail.getY(), spoketail.getZ());
                vtkIdType id1 = hubpos->InsertNextPoint(spoketip.getX(), spoketip.getY(), spoketip.getZ());

                // Add to render
                medialsheetline->GetPointIds()->SetId(0, id0);
                medialsheetline->GetPointIds()->SetId(1, id1);
                lineCellArray->InsertNextCell(medialsheetline);
            }

            // Add down spokes
            if(showdown){
                // Get down spoke tail and tip
                spoketip = prim->getX() + prim->getR1()*prim->getU1();

                vtkSmartPointer<vtkLine> medialsheetline = vtkSmartPointer<vtkLine>::New();

                vtkIdType id0 = hubpos->InsertNextPoint(spoketail.getX(), spoketail.getY(), spoketail.getZ());
                vtkIdType id1 = hubpos->InsertNextPoint(spoketip.getX(), spoketip.getY(), spoketip.getZ());

                // Add to render
                medialsheetline->GetPointIds()->SetId(0, id0);
                medialsheetline->GetPointIds()->SetId(1, id1);
                lineCellArray->InsertNextCell(medialsheetline);
            }

            // Add crest spokes
            if(u==0 || u==rowNum-1 || v==0 || v==colNum-1){
                if(showcrest){
                    endPrim = dynamic_cast<M3DQuadEndPrimitive*>(quadFig->getPrimitivePtr(u, v));

                    spoketip = endPrim->getX() + endPrim->getREnd()*endPrim->getUEnd();

                    vtkSmartPointer<vtkLine> medialsheetline = vtkSmartPointer<vtkLine>::New();

                    vtkIdType id0 = hubpos->InsertNextPoint(spoketail.getX(), spoketail.getY(), spoketail.getZ());
                    vtkIdType id1 = hubpos->InsertNextPoint(spoketip.getX(), spoketip.getY(), spoketip.getZ());

                    // Add to render
                    medialsheetline->GetPointIds()->SetId(0, id0);
                    medialsheetline->GetPointIds()->SetId(1, id1);
                    lineCellArray->InsertNextCell(medialsheetline);
                }
            }
        }
    }

    linePolyData->SetPoints(hubpos);
    linePolyData->SetLines(lineCellArray);

    // Create a lookup table to map cell data to colors
     vtkSmartPointer<vtkLookupTable> lut = vtkSmartPointer<vtkLookupTable>::New();
     int tableSize = 39;//std::max(3*3 + 1, 10);
     lut->SetNumberOfTableValues(tableSize);
     lut->Build();

     // Fill in a few known colors, the rest will be generated if needed
     lut->SetTableValue(0, 0.98,  0.98, 0.4300, 1);  //Blue
     lut->SetTableValue(1, 0.8900, 0.8100, 0.3400, 1); // Banana
     lut->SetTableValue(2, 1.0000, 0.3882, 0.2784, 1); // Tomato
     lut->SetTableValue(3, 0.9608, 0.8706, 0.7020, 1); // Wheat
     lut->SetTableValue(4, 0.9020, 0.9020, 0.9804, 1); // Lavender
     lut->SetTableValue(5, 1.0000, 0.4900, 0.2500, 1); // Flesh
     lut->SetTableValue(6, 0.5300, 0.1500, 0.3400, 1); // Raspberry
     lut->SetTableValue(7, 0.9804, 0.5020, 0.4471, 1); // Salmon
     lut->SetTableValue(8, 0.7400, 0.9900, 0.7900, 1); // Mint
     lut->SetTableValue(9, 0.2000, 0.6300, 0.7900, 1); // Peacock
     lut->SetTableValue(10,   0.92,  0.2100, 0.5300, 1);  //Black
     lut->SetTableValue(11, 1, 1, 1, 1); // Banana
     lut->SetTableValue(12, 1.0000, 0.7882, 0.2784, 1); // Tomato
     lut->SetTableValue(13, 0.0608, 0.2706, 0.9020, 1); // Wheat
     lut->SetTableValue(14, 0.9020, 0.9020, 0.1804, 1); // Lavender
     lut->SetTableValue(15, 1.0000, 0.4900, 0.7500, 1); // Flesh
     lut->SetTableValue(16, 0.5300, 0.1500, 1, 1); // Raspberry
     lut->SetTableValue(17, 1, 0, 0, 1); // Salmon
     lut->SetTableValue(18, 0.7400, 0.0900, 0.7900, 1); // Mint
     lut->SetTableValue(19, 0.2000, 0.6300, 0.2900, 1); // Peacock
     lut->SetTableValue(20  , 0.3,  0.8100, 0.3400, 1);  //Black
     lut->SetTableValue(21, 0.8900, 0.0100, 0.3400, 1); // Banana
     lut->SetTableValue(22, 1.0000, 1, 0.2784, 1); // Tomato
     lut->SetTableValue(23, 0.9608, 0.1706, 0.7020, 1); // Wheat
     lut->SetTableValue(24, 0, 0.9020, 0.9804, 1); // Lavender
     lut->SetTableValue(25, 0.7000, 0.6900, 0.2500, 1); // Flesh
     lut->SetTableValue(26, 0.5300, 0.1500, 0.1400, 1); // Raspberry
     lut->SetTableValue(27, 0.9804, 0.5020, 0.1471, 1); // Salmon
     lut->SetTableValue(28, 0.7400, 0.5900, 0.7900, 1); // Mint
     lut->SetTableValue(29, 0.2000, 0.1300, 0.7900, 1); // Peacock
     lut->SetTableValue(30, 0.7,  0.5300, 0.3400, 1);  //Black
     lut->SetTableValue(31, 0.5900, 0.3100, 0.3400, 1); // Banana
     lut->SetTableValue(32, 0, 0.3882, 0, 1); // Tomato
     lut->SetTableValue(33, 0, 0.8, 0, 1); // Wheat
     lut->SetTableValue(34, 0.5020, 0.9020, 0.1804, 1); // Lavender
     lut->SetTableValue(35, 1.0000, 0.4900, 0.2500, 1); // Flesh
     lut->SetTableValue(36, 0.5300, 0.1500, 0.3400, 1); // Raspberry
     lut->SetTableValue(37, 0.9804, 0.5020, 0.4471, 1); // Salmon
     lut->SetTableValue(38, 0.9, 0.8, 0.1900, 1); // Mint

       // Generate the colors for each point based on the color map
       vtkSmartPointer<vtkUnsignedCharArray> colors = vtkSmartPointer<vtkUnsignedCharArray>::New();
       colors->SetNumberOfComponents(3);
       colors->SetName("Colors");


       for(int i = 0; i < tableSize; i++) {
           double rgb[3];
           lut->GetColor(static_cast<double>(i)/(tableSize-1), rgb);

           unsigned char ucrgb[3];
           for(unsigned int j = 0; j < 3; j++) {
               ucrgb[j] = static_cast<unsigned char>(255.0 * rgb[j]);
           }

           colors->InsertNextTuple3(ucrgb[0], ucrgb[1], ucrgb[2]);
       }

    //linePolyData->GetCellData()->SetScalars(colors);

    vtkSmartPointer<vtkPolyDataMapper> pointlinemapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    pointlinemapper->SetInput(linePolyData);
    vtkSmartPointer<vtkActor>  lineactor = vtkActor::New();
    lineactor->SetMapper(pointlinemapper);
    lineactor->GetProperty()->SetLineWidth(linewith);
    lineactor->GetProperty()->SetColor(quadColor);

    lineactor->RotateX(rX);
    lineactor->RotateY(rY);
    lineactor->RotateZ(rZ);
    lineactor->GetProperty()->SetOpacity(this->opacity);
    this->renderer->AddActor(lineactor);
}


void visualizeobject::drawSrep_2(M3DQuadFigure* quadFig){

    int rowNum = quadFig->getRowCount();
    int colNum = quadFig->getColumnCount();

    M3DQuadPrimitive* prim;
    M3DQuadEndPrimitive* endPrim;

    Vector3D spoketail, spoketip;

    vtkSmartPointer< vtkPoints > hubpos = vtkSmartPointer< vtkPoints >::New();
    vtkSmartPointer<vtkCellArray> lineCellArray = vtkSmartPointer<vtkCellArray>::New();
    vtkSmartPointer<vtkPolyData> linePolyData = vtkSmartPointer<vtkPolyData>::New();

    // Loop each atom of the quadfig
    for(unsigned u = 0; u < rowNum; u++){
        for(unsigned v = 0; v < colNum; v++){

            prim = dynamic_cast<M3DQuadPrimitive*>(quadFig->getPrimitivePtr(u,v));
            spoketail = prim->getX();

            // Add up spokes

                // Get up spoke tail and tip
                spoketip = prim->getX() + prim->getR0()*prim->getU0();

                vtkSmartPointer<vtkLine> medialsheetline1 = vtkSmartPointer<vtkLine>::New();

                vtkIdType id0 = hubpos->InsertNextPoint(spoketail.getX(), spoketail.getY(), spoketail.getZ());
                vtkIdType id1 = hubpos->InsertNextPoint(spoketip.getX(), spoketip.getY(), spoketip.getZ());

                // Add to render
                medialsheetline1->GetPointIds()->SetId(0, id0);
                medialsheetline1->GetPointIds()->SetId(1, id1);
                lineCellArray->InsertNextCell(medialsheetline1);

            // Add down spokes
            // Get down spoke tail and tip
            spoketip = prim->getX() + prim->getR1()*prim->getU1();

            vtkSmartPointer<vtkLine> medialsheetline2 = vtkSmartPointer<vtkLine>::New();

            id0 = hubpos->InsertNextPoint(spoketail.getX(), spoketail.getY(), spoketail.getZ());
            id1 = hubpos->InsertNextPoint(spoketip.getX(), spoketip.getY(), spoketip.getZ());

            // Add to render
            medialsheetline2->GetPointIds()->SetId(0, id0);
            medialsheetline2->GetPointIds()->SetId(1, id1);
            lineCellArray->InsertNextCell(medialsheetline2);

            // Add crest spokes
            if(u==0 || u==rowNum-1 || v==0 || v==colNum-1){
                    endPrim = dynamic_cast<M3DQuadEndPrimitive*>(quadFig->getPrimitivePtr(u, v));

                    spoketip = endPrim->getX() + endPrim->getREnd()*endPrim->getUEnd();

                    vtkSmartPointer<vtkLine> medialsheetline3 = vtkSmartPointer<vtkLine>::New();

                    id0 = hubpos->InsertNextPoint(spoketail.getX(), spoketail.getY(), spoketail.getZ());
                    id1 = hubpos->InsertNextPoint(spoketip.getX(), spoketip.getY(), spoketip.getZ());

                    // Add to render
                    medialsheetline3->GetPointIds()->SetId(0, id0);
                    medialsheetline3->GetPointIds()->SetId(1, id1);
                    lineCellArray->InsertNextCell(medialsheetline3);
            }
        }
    }

    linePolyData->SetPoints(hubpos);
    linePolyData->SetLines(lineCellArray);

    // Create a lookup table to map cell data to colors
     vtkSmartPointer<vtkLookupTable> lut = vtkSmartPointer<vtkLookupTable>::New();
     int tableSize = 106;//std::max(3*3 + 1, 10);
     lut->SetNumberOfTableValues(tableSize);
     lut->Build();

     // Fill in a few known colors, the rest will be generated if needed
     lut->SetTableValue(0, 0,  0.5, 1, 1);
     lut->SetTableValue(1, 0.8900, 0.8100, 0.3400, 1); // Banana
     lut->SetTableValue(2, 1.0000, 0.3882, 0.2784, 1); // Tomato
     lut->SetTableValue(3, 0.9608, 0.8706, 0.7020, 1); // Wheat
     lut->SetTableValue(4, 0.9020, 0.9020, 0.9804, 1); // Lavender
     lut->SetTableValue(5, 1.0000, 0.4900, 0.2500, 1); // Flesh
     lut->SetTableValue(6, 0.5300, 0.1500, 0.3400, 1); // Raspberry
     lut->SetTableValue(7, 0.9804, 0.5020, 0.4471, 1); // Salmon
     lut->SetTableValue(8, 0.7400, 0.9900, 0.7900, 1); // Mint
     lut->SetTableValue(9, 0.2000, 0.6300, 0.7900, 1); // Peacock
     lut->SetTableValue(10,   0.92,  0.2100, 0.5300, 1);  //Black
     lut->SetTableValue(11, 1, 0.2, 0, 1); // Banana
     lut->SetTableValue(12, 1.0000, 0.7882, 0.2784, 1); // Tomato
     lut->SetTableValue(13, 0, 1, 0, 1); // Wheat
     lut->SetTableValue(14, 0.9020, 0.9020, 0.1804, 1); // Lavender
     lut->SetTableValue(15, 1.0000, 0.4900, 0.7500, 1); // Flesh
     lut->SetTableValue(16, 0.5300, 0.1500, 1, 1); // Raspberry
     lut->SetTableValue(17, 1, 1, 0, 1); // Salmon
     lut->SetTableValue(18, 0.7400, 0.0900, 0.7900, 1); // Mint
     lut->SetTableValue(19, 0.2000, 0.6300, 0.2900, 1); // Peacock
     lut->SetTableValue(20  , 0.3,  0.8100, 0.3400, 1);  //Black
     lut->SetTableValue(21, 0.8900, 0.0100, 0.3400, 1); // Banana
     lut->SetTableValue(22, 1.0000, 1, 0.2784, 1); // Tomato
     lut->SetTableValue(23, 0.9608, 0.1706, 0.7020, 1); // Wheat
     lut->SetTableValue(24, 0, 0.9020, 0.9804, 1); // Lavender
     lut->SetTableValue(25, 0.7000, 0.6900, 0.2500, 1); // Flesh
     lut->SetTableValue(26, 0.5300, 0.1500, 0.1400, 1); // Raspberry
     lut->SetTableValue(27, 0.9804, 0.5020, 0.1471, 1); // Salmon
     lut->SetTableValue(28, 0.7400, 0.5900, 0.7900, 1); // Mint
     lut->SetTableValue(29, 0.2000, 0.1300, 0.7900, 1); // Peacock
     lut->SetTableValue(30, 0.7,  0.5300, 0.3400, 1);  //Black
     lut->SetTableValue(31, 0.5900, 0.3100, 0.3400, 1); // Banana
     lut->SetTableValue(32, 0, 0.3882, 0, 1); // Tomato
     lut->SetTableValue(33, 0, 0.8, 0, 1); // Wheat
     lut->SetTableValue(34, 0.5020, 0.9020, 0.1804, 1); // Lavender
     lut->SetTableValue(35, 1.0000, 0.4900, 0.2500, 1); // Flesh
     lut->SetTableValue(36, 0.5300, 0.1500, 0.3400, 1); // Raspberry
     lut->SetTableValue(37, 0.9804, 0.5020, 0.4471, 1); // Salmon
     lut->SetTableValue(38, 1, 0.1, 1, 1); // Mint
     lut->SetTableValue(39, 0.9, 0.1, 0.1, 1); // Mint
     lut->SetTableValue(40, 0.9000, 0.6900, 0.2500, 1); // Flesh
     lut->SetTableValue(41, 0.7,  0.5300, 0.1400, 1);  //Black
     lut->SetTableValue(42, 0.98,  0.98, 0.0300, 1);  //Blue
     lut->SetTableValue(43, 1.0000, 1, 0.2784, 1); // Tomato
     lut->SetTableValue(44, 0.2900, 0.8100, 0.3400, 1); // Banana
     lut->SetTableValue(45, 1.0000, 0.9882, 0.2784, 1); // Tomato
     lut->SetTableValue(46, 0.9608, 0.4706, 0.7020, 1); // Wheat
     lut->SetTableValue(47, 0, 0.9, 1, 1); // Lavender
     lut->SetTableValue(48, 1.0000, 0.4900, 0.5500, 1); // Flesh
     lut->SetTableValue(49, 0.5300, 0.2500, 0.3400, 1); // Raspberry
     lut->SetTableValue(50, 0.9804, 0.5920, 0.4971, 1); // Salmon
     lut->SetTableValue(51, 0.7400, 0.9900, 0.0900, 1); // Mint
     lut->SetTableValue(52, 0.2000, 1.300, 0.7900, 1); // Peacock
     lut->SetTableValue(53,   0.2,  0.2100, 0.8900, 1);  //Black
     lut->SetTableValue(54, 0, 0.3, 0, 1); // Banana
     lut->SetTableValue(55, 1.0000, 0.8, 0, 1); // Tomato
     lut->SetTableValue(56, 1, 1, 0, 1); // Wheat
     lut->SetTableValue(57, 1, 0, 1, 1); // Lavender
     lut->SetTableValue(58, 1.0000, 0.1, 0.8, 1); // Flesh
     lut->SetTableValue(59, 0.5300, 0.1500, 1, 1); // Raspberry
     lut->SetTableValue(60, 0, 0.8, 1, 1); // Salmon
     lut->SetTableValue(61, 0.7400, 0.3900, 0.7900, 1); // Mint
     lut->SetTableValue(62, 0.2000, 0.9300, 0.2900, 1); // Peacock
     lut->SetTableValue(63  , 0.3,  0.6100, 0.3400, 1);  //Black
     lut->SetTableValue(64, 0.8900, 0.6100, 0.3400, 1); // Banana
     lut->SetTableValue(65, 0.9608, 0.4706, 0.7020, 1); // Wheat
     lut->SetTableValue(66, 0, 0.9020, 0.0804, 1); // Lavender
     lut->SetTableValue(67, 1, 1, 1, 1);
     lut->SetTableValue(68, 0.2300, 0.1500, 0.1400, 1); // Raspberry
     lut->SetTableValue(69, 0.9804, 0.5020, 0.6471, 1); // Salmon
     lut->SetTableValue(70, 0.7400, 0.5900, 0.9900, 1); // Mint
     lut->SetTableValue(71, 0.2000, 0.1900, 0.7900, 1); // Peacock
     lut->SetTableValue(72, 0.5900, 0.4500, 0.3400, 1); // Banana
     lut->SetTableValue(73, 0, 0.6682, 0, 1); // Tomato
     lut->SetTableValue(74, 0, 0.82, 0, 1); // Wheat
     lut->SetTableValue(75, 1, 0, 0, 1); // Lavender
     lut->SetTableValue(76, 0.5900, 0.9500, 0.3400, 1); // Raspberry
     lut->SetTableValue(77, 1, 0, 1, 1); // Salmon

     lut->SetTableValue(78, 0.9020, 0.7020, 0.3804, 1);
     lut->SetTableValue(79, 0.8900, 0.3100, 0.3400, 1);
     lut->SetTableValue(80, 1.0000, 0.1882, 0.2784, 1);
     lut->SetTableValue(81, 0.1608, 0.8706, 0.7020, 1);
     lut->SetTableValue(82, 1, 0.4, 0, 1);
     lut->SetTableValue(83, 1.0000, 0.1, 0.1, 1);
     lut->SetTableValue(84, 0.3, 1, 1, 1);
     lut->SetTableValue(85, 1, 1, 0.1, 1);
     lut->SetTableValue(86, 0.7400, 0.9900, 0.7900, 1);
     lut->SetTableValue(87, 0.2000, 0.6300, 0.7900, 1);
     lut->SetTableValue(88,   0.92,  0.2100, 0.5300, 1);
     lut->SetTableValue(89, 1, 0.1, 1, 1);
     lut->SetTableValue(90, 1.0000, 0.7882, 0.2784, 1);
     lut->SetTableValue(91, 0.0608, 0.2706, 0.9020, 1);
     lut->SetTableValue(92, 0.9020, 0.9020, 0.1804, 1);
     lut->SetTableValue(93, 1.0000, 0.4900, 0.7500, 1);
     lut->SetTableValue(94, 0.5300, 0.1500, 1, 1);
     lut->SetTableValue(95, 1, 0, 0, 1);
     lut->SetTableValue(96, 0.7400, 0.0900, 0.7900, 1);
     lut->SetTableValue(97, 0.2000, 0.6300, 0.2900, 1);
     lut->SetTableValue(98  , 0.3,  0.8100, 0.3400, 1);
     lut->SetTableValue(99, 0.8900, 0.0100, 0.3400, 1);
     lut->SetTableValue(100, 1, 1, 1, 1);
     lut->SetTableValue(101, 0.9608, 0.1706, 0.7020, 1);
     lut->SetTableValue(102, 0, 0.9020, 0.9804, 1);
     lut->SetTableValue(103, 0.7000, 0.6900, 0.2500, 1);
     lut->SetTableValue(104, 0.5300, 0.1500, 0.1400, 1);
     lut->SetTableValue(105, 0.9020, 0.7020, 0.1804, 1); //

       // Generate the colors for each point based on the color map
       vtkSmartPointer<vtkUnsignedCharArray> colors = vtkSmartPointer<vtkUnsignedCharArray>::New();
       colors->SetNumberOfComponents(3);
       colors->SetName("Colors");


       for(int i = 0; i < tableSize; i++) {
           double rgb[3];
           lut->GetColor(static_cast<double>(i)/(tableSize-1), rgb);

           unsigned char ucrgb[3];
           for(unsigned int j = 0; j < 3; j++) {
               ucrgb[j] = static_cast<unsigned char>(255.0 * rgb[j]);
           }

           colors->InsertNextTuple3(ucrgb[0], ucrgb[1], ucrgb[2]);
       }

    linePolyData->GetCellData()->SetScalars(colors);

    vtkSmartPointer<vtkPolyDataMapper> pointlinemapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    pointlinemapper->SetInput(linePolyData);
    vtkSmartPointer<vtkActor>  lineactor = vtkActor::New();
    lineactor->SetMapper(pointlinemapper);
    lineactor->GetProperty()->SetLineWidth(3);
    //lineactor->GetProperty()->SetColor(1,1,1);

    lineactor->RotateX(rX);
    lineactor->RotateY(rY);
    lineactor->RotateZ(rZ);
    lineactor->GetProperty()->SetOpacity(this->opacity);
    this->renderer->AddActor(lineactor);
}



/* Draw the srep. */
void visualizeobject::showSrep(M3DQuadFigure* quadFig, int interpolationLevel, int linewith){
    double quadColor[3] = {0,0.8,0}; // {1,0.6,1}; magenta  //{0,0.8,0};//Green     //{1,0.5,0} orange
    showStandardSrep(quadFig, 0, quadColor, interpolationLevel,linewith, true); //original srep has up and down skeletal sheet stick together, here can be 0 or 1.
    double quadColor1[3] = {0,0.8,0};//{0.6,1,1}; //{0,0,1} pulper
    showStandardSrep(quadFig, 1, quadColor, interpolationLevel, linewith, true); //original srep has up and down skeletal sheet stick together, here can be 0 or 1.
    showCrestSrep(quadFig, quadColor, interpolationLevel,linewith, true);
}

/* Draw the srep. */
void visualizeobject::showSrep_2(M3DQuadFigure* quadFig, int interpolationLevel, int linewith){
    double quadColor[3] = {1,0,1}; // {1,0.6,1}; magenta  //{0,0.8,0};//Green     //{1,0.5,0} orange
    //showStandardSrep(quadFig, 0, quadColor, interpolationLevel,linewith, false); //original srep has up and down skeletal sheet stick together, here can be 0 or 1.
    double quadColor1[3] = {0,1,1}; //cyan
    //showStandardSrep(quadFig, 1, quadColor1, interpolationLevel, linewith, false); //original srep has up and down skeletal sheet stick together, here can be 0 or 1.
    double quadColor2[3] = {1,0,0}; //pulper
    showCrestSrep(quadFig, quadColor2, interpolationLevel,linewith, true);
}


/* Draw a srep's skeletal sheet and spokes, and the shifted srep together.
 * We can set the deltau and deltav to set the shift distance.
 * The logic to set variables for each spoke is the same as class movespokes.cpp.
 * coeff: deltau and deltav variables for this shift. Its size is: varStand *2 + varCrest.
 * coeff: holding all the variables for each srep, the first varStand is for up spokes,
 * the next varStand for down spokes; the rest varCrest is for crest (atom's) crest spokes.
 * onlyShifted: if set to true, will only show the shifted new srep, without show the original one.
*/
void visualizeobject::showShiftedSrep(vector<double> coeff, bool onlyShifted, int varStand, int varCrest, int interpolationLevel, int linewidth
                                      , bool showspoketail){
    this->interpolationLevel = interpolationLevel;

    if(!onlyShifted){
        // Add original srep into render
        showSrep(this->quadFig, 0, linewidth);
        cout<<"onlyShifted is ture, should not output this!!!"<<endl;
    }

    slidestandardspokes mStandardSpoke;
    double *coeffStand = new double [varStand]; // storing up or down spokes variables.
    double *coeffCrest = new double [varCrest];

    // The first varStand is for up spokes.
    for(unsigned int i = 0; i < varStand; i++){
        coeffStand[i] = coeff[i];
    }
    // Shift up spokes
    M3DQuadFigure* upQuadFig = dynamic_cast<M3DQuadFigure*>(this->quadFig->clone()); //Move all base on original srep.
    mStandardSpoke.updateSpokeInfo(upQuadFig, coeffStand, 0);

    // Add shifted srep's up skeletal into render
    double quadColor_up[3] = {1,0.6,1};//{1,0.7,0}; orange     //{0.5,0.8,1};
    // Add shifted up srep into render
    showStandardSrep(upQuadFig, 0, quadColor_up, this->interpolationLevel, linewidth, showspoketail);

    // The following varStand variables are for down spokes;
    for(unsigned int i = 0; i < varStand; i++){
        coeffStand[i] = coeff[varStand + i];
    }
    // Shift down spokes
    M3DQuadFigure* downQuadFig = dynamic_cast<M3DQuadFigure*>(this->quadFig->clone()); //Move all base on original srep.
    mStandardSpoke.updateSpokeInfo(downQuadFig, coeffStand, 1); //Move base on original srep.
    // Add shifted srep's down skeletal into render
    double quadColor_down[3] = {0.6,1,1};//{0,0.7,1};dark cyan //{0.1,1,0};
    // Add shifted down srep into render
    showStandardSrep(downQuadFig, 1, quadColor_down, this->interpolationLevel, linewidth, showspoketail);

    // The last 28 variables in coeff is the variables for crest spokes.    
    for(unsigned int i = 0; i < varCrest; i++){
        coeffCrest[i] = coeff[varStand*2 + i];        
    }    

    // Shift crest spokes
    slidecrestspokes mCrestSpoke;
    M3DQuadFigure* crestQuadFig = dynamic_cast<M3DQuadFigure*>(this->quadFig->clone()); //Move all base on original srep.
    visualizecrest vcrest;
    vtkSmartPointer<vtkSRep> srepFig = vcrest.getSrepFig(this->quadFig);
    mCrestSpoke.updateCrestSpokeInfo(srepFig, crestQuadFig, coeffCrest); //Move base on original srep.

    // Add shifted crest srep
    showCrestSrep(crestQuadFig, this->quadColor, this->interpolationLevel, linewidth, showspoketail);

    delete coeffStand;
    delete coeffCrest;
}


void visualizeobject::showShiftedUpOrDownSideSkeletalSheet(vector<double> coeff, double * color, int side, double opacity){

    this->interpolationLevel = 0;

    int varNum = coeff.size();

    slidestandardspokes mStandardSpoke;
    double *coeffStand = new double [varNum]; // storing up or down spokes variables.

    // The first varStand is for up spokes.
    for(unsigned int i = 0; i < varNum; i++){
        coeffStand[i] = coeff[i];
    }

    // Shift up spokes
    M3DQuadFigure* quadFig = dynamic_cast<M3DQuadFigure*>(this->quadFig->clone()); //Move all base on original srep.
    mStandardSpoke.updateSpokeInfo(quadFig, coeffStand, side);

    vtkSmartPointer< vtkPoints > hubPosition_s = vtkSmartPointer< vtkPoints >::New();

    // draw the sub quad frame, interpolation level 2
    visualization vo(quadFig, 2, side, true, color, renderer, rX, rY, rZ, opacity);

    hubPosition_s = vo.getSubQuadsPosition(1); //0(boundary), 1(skeletal)

    vo.drawQuads(hubPosition_s);


    delete coeffStand;
}


void visualizeobject::showShiftedCrestCurve(vector<double> coeff, double * color, double opacity){

    this->interpolationLevel = 1;

    int varNum = coeff.size();

    double *coeffCrest = new double [varNum]; // storing up or down spokes variables.

    // The first varStand is for up spokes.
    for(unsigned int i = 0; i < varNum; i++){
        coeffCrest[i] = coeff[i];
    }

    // Shift crest spokes
    slidecrestspokes mCrestSpoke;
    M3DQuadFigure* crestQuadFig = dynamic_cast<M3DQuadFigure*>(this->quadFig->clone()); //Move all base on original srep.
    visualizecrest vcrest;
    vtkSmartPointer<vtkSRep> srepFig = vcrest.getSrepFig(this->quadFig);
    mCrestSpoke.updateCrestSpokeInfo(srepFig, crestQuadFig, coeffCrest); //Move base on original srep.

    delete coeffCrest;

    visualizecrest visualObject_crest(crestQuadFig, this->interpolationLevel, color, this->renderer, this->rX, this->rY, this->rZ, opacity);

    // Add crest spokes to render
    visualObject_crest.drawCrestCurve(2, color);

}




/* Draw a srep's skeletal sheet and spokes, and the shifted srep together.
 * We can set the deltau and deltav to set the shift distance.
 * The logic to set variables for each spoke is the same as class movespokes.cpp.
 * coeff: deltau and deltav variables for this shift. Its size is: varStand *2 + varCrest.
 * coeff: holding all the variables for each srep, the first varStand is for up spokes,
 * the next varStand for down spokes; the rest varCrest is for crest (atom's) crest spokes.
 * onlyShifted: if set to true, will only show the shifted new srep, without show the original one.
*/
void visualizeobject::showShiftedUpOrDownSideSrep(vector<double> coeff, double * color, int side, int linewith, bool showspoketail){

    this->interpolationLevel = 0;

    int varNum = coeff.size();

    slidestandardspokes mStandardSpoke;
    double *coeffStand = new double [varNum]; // storing up or down spokes variables.

    // The first varStand is for up spokes.
    for(unsigned int i = 0; i < varNum; i++){
        coeffStand[i] = coeff[i];
    }

    // Shift up spokes
    M3DQuadFigure* quadFig = dynamic_cast<M3DQuadFigure*>(this->quadFig->clone()); //Move all base on original srep.
    mStandardSpoke.updateSpokeInfo(quadFig, coeffStand, side);

    // Add shifted srep's up skeletal into render
    //double quadColor[3] = {1,0.8,0};//{1,0.7,0}; orange     //{0.5,0.8,1};

    // Add shifted up srep into render
    showStandardSrep(quadFig, side, color, this->interpolationLevel, linewith, showspoketail);

    delete coeffStand;
}



void visualizeobject::showShiftedCrestSideSrep(vector<double> coeff, double * color, int linewidth, bool showspoketail){

    this->interpolationLevel = 1;

    int varNum = coeff.size();

    double *coeffCrest = new double [varNum]; // storing up or down spokes variables.

    // The first varStand is for up spokes.
    for(unsigned int i = 0; i < varNum; i++){
        coeffCrest[i] = coeff[i];
    }

    // Shift crest spokes
    slidecrestspokes mCrestSpoke;
    M3DQuadFigure* crestQuadFig = dynamic_cast<M3DQuadFigure*>(this->quadFig->clone()); //Move all base on original srep.
    visualizecrest vcrest;
    vtkSmartPointer<vtkSRep> srepFig = vcrest.getSrepFig(this->quadFig);
    mCrestSpoke.updateCrestSpokeInfo(srepFig, crestQuadFig, coeffCrest); //Move base on original srep.

    // Add shifted crest srep
    showCrestSrep(crestQuadFig, color, this->interpolationLevel, linewidth, showspoketail);

    delete coeffCrest;

}





/* Show the crest frame.*/
void visualizeobject::showShiftedCrestFrame(vector<double> crest_vars, int varCrest){

    double *coeffCrest = new double [varCrest];

    if(crest_vars.size() != varCrest){
        cout<<"Error from visualizeobject::showShiftedCrestFrame: the crest variables is not set correctly!"<<endl;
        EXIT_FAILURE;
    }

    // The last 28 variables in coeff is the variables for crest spokes.
    for(unsigned int i = 0; i < varCrest; i++){
        coeffCrest[i] = crest_vars[i];
    }

    // Shift crest spokes
    slidecrestspokes mCrestSpoke;
    M3DQuadFigure* crestQuadFig = dynamic_cast<M3DQuadFigure*>(this->quadFig->clone()); //Move all base on original srep.
    visualizecrest vcrest;
    vtkSmartPointer<vtkSRep> srepFig = vcrest.getSrepFig(this->quadFig);
    mCrestSpoke.updateCrestSpokeInfo(srepFig, crestQuadFig, coeffCrest); //Move base on original srep.

    // Add shifted crest frame.
    //showBoundaryFrames(crestQuadFig, false, false, true);

    delete coeffCrest;
}


/*show the crest frame and spokes together.*/
void visualizeobject::showShiftedCrestFrameAndSpokes(vector<double> crest_vars, int varCrest, int interpolationLevel){
    this->interpolationLevel = interpolationLevel;

    double *coeffCrest = new double [varCrest];

    if(crest_vars.size() != varCrest){
        cout<<"Error from visualizeobject::showShiftedCrestFrame: the crest variables is not set correctly!"<<endl;
        EXIT_FAILURE;
    }

    // The last 28 variables in coeff is the variables for crest spokes.
    for(unsigned int i = 0; i < varCrest; i++){
        coeffCrest[i] = crest_vars[i];
    }

    // Shift crest spokes
    slidecrestspokes mCrestSpoke;
    M3DQuadFigure* crestQuadFig = dynamic_cast<M3DQuadFigure*>(this->quadFig->clone()); //Move all base on original srep.

    visualizecrest vcrest;
    vtkSmartPointer<vtkSRep> srepFig = vcrest.getSrepFig(this->quadFig);
    mCrestSpoke.updateCrestSpokeInfo(srepFig, crestQuadFig, coeffCrest); //Move base on original srep.

    // Add shifted crest frame.
    //showBoundaryFrames(crestQuadFig, false, false, true);

    visualizecrest visualObject_crest(crestQuadFig, this->interpolationLevel, quadColor, this->renderer, this->rX, this->rY, this->rZ, this->opacity);

    // Add crest spokes to render
    visualObject_crest.drawCrestSpokes();
    visualObject_crest.drawCrestCurve(0.5, quadColor);

    delete coeffCrest;
}


/* Show the atoms on the medial sheet in green sphere.
   (.3, .6, .3); // Background color green
*/
void visualizeobject::showAtomPoints(M3DQuadFigure* quadFig){
    // Get the hub position of each atom
    M3DQuadPrimitive* prim;
    for(unsigned int i = 0; i < quadFig->getRowCount(); i++){
        for(unsigned int j = 0; j < quadFig->getColumnCount(); j++){
            prim = dynamic_cast<M3DQuadPrimitive*>(quadFig->getPrimitivePtr(i,j));

            Vector3D hubPosition = prim->getX();

            // Create a sphere
            vtkSmartPointer<vtkSphereSource> sphereSource = vtkSmartPointer<vtkSphereSource>::New();
            sphereSource->SetCenter(hubPosition.getX(), hubPosition.getY(), hubPosition.getZ());
            sphereSource->SetRadius(0.0012); // for hipcampus //0.0010  // 0.0041

            vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
            mapper->SetInputConnection(sphereSource->GetOutputPort());

            vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
            actor->SetMapper(mapper);
            actor->GetProperty()->SetColor(1,1,1);//1,1,0.1 yellow    // 0,0.9,0 green
            actor->RotateX(rX);
            actor->RotateY(rY);
            actor->RotateZ(rZ);
            this->renderer->AddActor(actor);
        }
    }
}


void visualizeobject::showAtomPoints_2(M3DQuadFigure* quadFig, double radius, double*color){
    // Get the hub position of each atom
    M3DQuadPrimitive* prim;
    for(unsigned int i = 0; i < quadFig->getRowCount(); i++){
        for(unsigned int j = 0; j < quadFig->getColumnCount(); j++){
            prim = dynamic_cast<M3DQuadPrimitive*>(quadFig->getPrimitivePtr(i,j));

            Vector3D hubPosition = prim->getX();

            // Create a sphere
            vtkSmartPointer<vtkSphereSource> sphereSource = vtkSmartPointer<vtkSphereSource>::New();
            sphereSource->SetCenter(hubPosition.getX(), hubPosition.getY(), hubPosition.getZ());
            sphereSource->SetRadius(radius); // for hipcampus //0.0010  // 0.0041

            vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
            mapper->SetInputConnection(sphereSource->GetOutputPort());

            vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
            actor->SetMapper(mapper);
            actor->GetProperty()->SetColor(color);//1,1,0.1 yellow    // 0,0.9,0 green
            actor->RotateX(rX);
            actor->RotateY(rY);
            actor->RotateZ(rZ);
            this->renderer->AddActor(actor);
        }
    }
}



void visualizeobject::showCrestFrameByQuadIndex(M3DQuadFigure* quadFig, int quadIndex, int interpolationLevel, double * quadColor){
    visualizecrest visualObject_crest(quadFig, interpolationLevel, quadColor, this->renderer, this->rX, this->rY, this->rZ, this->opacity);
    visualObject_crest.drawCrestFramesByQuadIndex(quadIndex, quadColor);
}

void visualizeobject::showCrestQuadByQuadIndex(M3DQuadFigure* quadFig, int quadIndex, int interpolationLevel, double * quadColor){
    visualizecrest visualObject_crest(quadFig, interpolationLevel, quadColor, this->renderer, this->rX, this->rY, this->rZ, this->opacity);
    visualObject_crest.drawCrestQuadByQuadIndex(quadIndex, quadColor);
    //visualObject_crest.drawCrestQuadWithDiffColor(18);
}

void visualizeobject::showCrestSpokesByQuadIndex(M3DQuadFigure* quadFig, int quadIndex, int interpolationLevel, double * quadColor){
    visualizecrest visualObject_crest(quadFig, interpolationLevel, quadColor, this->renderer, this->rX, this->rY, this->rZ, this->opacity);
    visualObject_crest.drawCrestSpokesByQuadIndex(quadIndex, quadColor);
}


void visualizeobject::showCrestVolumesByQuadIndex(M3DQuadFigure* quadFig, int quadIndex, int interpolationLevel, double * quadColor){
    visualizecrest visualObject_crest(quadFig, interpolationLevel, quadColor, this->renderer, this->rX, this->rY, this->rZ, this->opacity);
    visualObject_crest.drawCrestVolumesByQuadIndex(quadIndex, quadColor);
}


void visualizeobject::showOriginalCrestSpokesByQuadIndex(M3DQuadFigure* quadFig, int quadIndex, int interpolationLevel, double * quadColor, bool showspoketail){
    visualizecrest visualObject_crest(quadFig, interpolationLevel, quadColor, this->renderer, this->rX, this->rY, this->rZ, this->opacity);
    visualObject_crest.drawMiddleCrestSpokesByQuadIndex(quadIndex, quadColor, showspoketail);
}

