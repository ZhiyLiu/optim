#include "visualization.h"



visualization::visualization()
{
}


/* Draw the edges(made up of many sub-lines) of quads. This function same as getEdgeOfQuadsUVCoordinateWithoutDeltaUV.
 * bool moved: true if draw the volume using the original r u p values in the input m3d file; false if draw the volume
 * using new values by adding deltaU and deltaV to primitive in the input m3d file.
*/
void visualization::drawEdgeOfQuads(bool moved){

    //get the correspondence delta u , v of each of the four points.
    M3DQuadPrimitive* prim0;
    M3DQuadInterpolater *tpm = new M3DQuadInterpolater(curQuadFig);

    vector<double> horizonaledgepointsu, horizonaledgepointsv, verticaledgepointsu, verticaledgepointsv;

    //Loop the left point of each quad in the srep.
    for(unsigned i = 0; i < rowNum; i++){
        for(unsigned j = 0; j < colNum; j++){

           //store the delta u and v of up and down spokes of Primitive[u][v].
           //deltau[0]: current point; deltau[1]: next horizonal point; deltau[2]: next vertical point
           double deltau[3]={0,0,0};
           double deltav[3]={0,0,0};

           if(!moved){
               switch(side){
                   case 0: //top side
                       //get the original medial sheet atom's delta u, v.
                       //considering horizonal edges & vertical edges.
                       //primitive[i][j]
                       prim0= dynamic_cast<M3DQuadPrimitive*>(curQuadFig->getPrimitivePtr(i, j));
                       deltau[0] = prim0->getDeltaU0();
                       deltav[0] = prim0->getDeltaV0();
                       //primitive[i][j+1]
                       if(j!=colNum-1){   //not the right side
                           prim0= dynamic_cast<M3DQuadPrimitive*>(curQuadFig->getPrimitivePtr(i, j+1));
                           deltau[1] = prim0->getDeltaU0();
                           deltav[1] = prim0->getDeltaV0();
                       }
                       //primitive[i+1][j]
                       if(i!=rowNum-1){    //not the bottom side
                           prim0= dynamic_cast<M3DQuadPrimitive*>(curQuadFig->getPrimitivePtr(i+1, j));
                           deltau[2] = prim0->getDeltaU0();
                           deltav[2] = prim0->getDeltaV0();
                       }
                       break;
                   case 1:  //down side
                       //get the original medial sheet atom's delta u, v.
                       //primitive[i][j]
                       prim0= dynamic_cast<M3DQuadPrimitive*>(curQuadFig->getPrimitivePtr(i, j));
                       deltau[0] = prim0->getDeltaU1();
                       deltav[0] = prim0->getDeltaV1();
                       //primitive[i][j+1]
                       if(j!=colNum-1){
                           prim0= dynamic_cast<M3DQuadPrimitive*>(curQuadFig->getPrimitivePtr(i, j+1));
                           deltau[1] = prim0->getDeltaU1();
                           deltav[1] = prim0->getDeltaV1();
                       }
                       //primitive[i+1][j]
                       if(i!=rowNum-1){
                           prim0= dynamic_cast<M3DQuadPrimitive*>(curQuadFig->getPrimitivePtr(i+1, j));
                           deltau[2] = prim0->getDeltaU1();
                           deltav[2] = prim0->getDeltaV1();
                       }
                       break;
                   default:
                       break;
               }
           }

           //new points position with delta u, v. In u, v coordinate.
           double quadpointsu[3], quadpointsv[3];
           //current point
           quadpointsu[0] = i + deltau[0];
           quadpointsv[0] = j + deltav[0];

           //next horizonal point
           quadpointsu[1] = i + deltau[1];
           quadpointsv[1] = j + 1 + deltav[1];

           //next vertical point
           quadpointsu[2] = i + 1 + deltau[2];
           quadpointsv[2] = j + deltav[2];

           cout<<"MSG from drawEdgeOfQuads: For pirmitve["<<i<<"]["<<j<<"], the sub-line-points are: "<<endl;
           if(j!=colNum-1){
           //calculate each of the horizonal edge, and store it in topboundaryquadedges
           tools.splitLine(quadpointsu[0], quadpointsu[1],step, horizonaledgepointsu);
           tools.splitLine(quadpointsv[0], quadpointsv[1],step, horizonaledgepointsv);

           for(unsigned m =0;m<horizonaledgepointsu.size();m++){
               cout<<"MSG from drawEdgeOfQuads: horizonaledgepointsu["<<m<<"] is: "<<horizonaledgepointsu[m]<<endl;
               cout<<"MSG from drawEdgeOfQuads: horizonaledgepointsv["<<m<<"] is: "<<horizonaledgepointsv[m]<<endl;
           }
           }

           if(i!=rowNum-1){
           //calculate each of the vertical edge, also store it in topboundaryquadedges
           tools.splitLine(quadpointsu[0], quadpointsu[2],step, verticaledgepointsu);
           tools.splitLine(quadpointsv[0], quadpointsv[2],step, verticaledgepointsv);
           for(unsigned m =0;m<verticaledgepointsu.size();m++){
               cout<<"MSG from drawEdgeOfQuads: verticaledgepointsu["<<m<<"] is: "<<verticaledgepointsu[m]<<endl;
               cout<<"MSG from drawEdgeOfQuads: verticaledgepointsv["<<m<<"] is: "<<verticaledgepointsv[m]<<endl;
           }
           }

           //store the horizonal two points, each point is a 3D vector, which store the x, y, z coordinate of this point.
           Vector3D horizonalpoint[2];
           //store the vertical two points, each point is a 3D vector, which store the x, y, z coordinate of this point.
           Vector3D verticalpoint[2];

           vtkSmartPointer<vtkPoints> hubpoints  = vtkSmartPointer<vtkPoints>::New();
           vtkSmartPointer<vtkCellArray> cellarraypointsline = vtkSmartPointer<vtkCellArray>::New();
           vtkSmartPointer<vtkPolyData> polypointsline = vtkSmartPointer<vtkPolyData>::New();

           //quadtype: Two kinds of quads, boundary quads(0), skeletal quads(1)
           switch(quadtype){
               case 0:
                   //boundary quads(0)
                   if(j!=colNum-1){
                       cout<<"MSG from drawEdgeOfQuads: For pirmitve["<<i<<"]["<<j<<"], has horizonal edges: "<<endl;
                       for(unsigned v = 0; v < step; v++){
                           horizonalpoint[0] = tpm->interpolateQuadSpoke(curQuadFig,horizonaledgepointsu[v],horizonaledgepointsv[v],side)->getB();
                           horizonalpoint[1] = tpm->interpolateQuadSpoke(curQuadFig,horizonaledgepointsu[v+1],horizonaledgepointsv[v+1],side)->getB();
                           cout<<"horizonaledgepointsu["<<v<<"] is: "<<horizonaledgepointsu[v]<<" horizonaledgepointsv["<<v<<"] is: "<<horizonaledgepointsv[v]<<endl;
                           cout<<"horizonalpoint[0].getX() is: "<<horizonalpoint[0].getX()<<" horizonalpoint[0].getY() is: "<<horizonalpoint[0].getY()<<" horizonalpoint[0].getZ() is: "<<horizonalpoint[0].getZ()<<endl;

                           // add the 4 sub-line to array.
                           vtkIdType id0 = hubpoints->InsertNextPoint(horizonalpoint[0].getX(), horizonalpoint[0].getY(), horizonalpoint[0].getZ());
                           vtkIdType id1 = hubpoints->InsertNextPoint(horizonalpoint[1].getX(), horizonalpoint[1].getY(), horizonalpoint[1].getZ());

                           vtkSmartPointer<vtkLine> medialsheetline = vtkSmartPointer<vtkLine>::New();
                           medialsheetline->GetPointIds()->SetId(0, id0);
                           medialsheetline->GetPointIds()->SetId(1, id1);

                           cellarraypointsline->InsertNextCell(medialsheetline);
                        }
                     }
                     if(i!=rowNum-1){
                         cout<<"MSG from drawEdgeOfQuads: For pirmitve["<<i<<"]["<<j<<"], has vertical edges: "<<endl;
                           for(unsigned f = 0; f < step; f++){
                               verticalpoint[0] = tpm->interpolateQuadSpoke(curQuadFig,verticaledgepointsu[f],verticaledgepointsv[f],side)->getB();
                               verticalpoint[1] = tpm->interpolateQuadSpoke(curQuadFig,verticaledgepointsu[f+1],verticaledgepointsv[f+1],side)->getB();
                              // cout<<"verticaledgepointsu["<<f<<"] is: "<<verticaledgepointsu[f]<<" verticaledgepointsv["<<f<<"] is: "<<verticaledgepointsv[f]<<endl;
                              // cout<<"verticalpoint[0].getX() is: "<<verticalpoint[0].getX()<<" verticalpoint[0].getY() is: "<<verticalpoint[0].getY()<<" verticalpoint[0].getZ() is: "<<verticalpoint[0].getZ()<<endl;
                               //drawLine(verticalpoint[0], verticalpoint[1], renderer);
                               // add the 4 sub-line to array.
                               vtkIdType id0 = hubpoints->InsertNextPoint(verticalpoint[0].getX(), verticalpoint[0].getY(), verticalpoint[0].getZ());
                               vtkIdType id1 = hubpoints->InsertNextPoint(verticalpoint[1].getX(), verticalpoint[1].getY(), verticalpoint[1].getZ());

                               vtkSmartPointer<vtkLine> medialsheetline = vtkSmartPointer<vtkLine>::New();
                               medialsheetline->GetPointIds()->SetId(0, id0);
                               medialsheetline->GetPointIds()->SetId(1, id1);

                               cellarraypointsline->InsertNextCell(medialsheetline);
                           }
                       }
                       break;
                   case 1:
                       //skeletal quads(1)
                       if(j!=colNum-1){
                           for(unsigned v = 0; v < step; v++){
                               horizonalpoint[0] = tpm->interpolateQuadSpoke(curQuadFig,horizonaledgepointsu[v],horizonaledgepointsv[v],side)->getX();
                               horizonalpoint[1] = tpm->interpolateQuadSpoke(curQuadFig,horizonaledgepointsu[v+1],horizonaledgepointsv[v+1],side)->getX();

                               // add the 4 sub-line to array.
                               vtkIdType id0 = hubpoints->InsertNextPoint(horizonalpoint[0].getX(), horizonalpoint[0].getY(), horizonalpoint[0].getZ());
                               vtkIdType id1 = hubpoints->InsertNextPoint(horizonalpoint[1].getX(), horizonalpoint[1].getY(), horizonalpoint[1].getZ());

                               vtkSmartPointer<vtkLine> medialsheetline = vtkSmartPointer<vtkLine>::New();
                               medialsheetline->GetPointIds()->SetId(0, id0);
                               medialsheetline->GetPointIds()->SetId(1, id1);

                               cellarraypointsline->InsertNextCell(medialsheetline);
                           }
                       }
                       if(i!=rowNum-1){
                           for(unsigned f = 0; f < step; f++){
                               verticalpoint[0] = tpm->interpolateQuadSpoke(curQuadFig,verticaledgepointsu[f],verticaledgepointsv[f],side)->getX();
                               verticalpoint[1] = tpm->interpolateQuadSpoke(curQuadFig,verticaledgepointsu[f+1],verticaledgepointsv[f+1],side)->getX();

                               // add the 4 sub-line to array.
                               vtkIdType id0 = hubpoints->InsertNextPoint(verticalpoint[0].getX(), verticalpoint[0].getY(), verticalpoint[0].getZ());
                               vtkIdType id1 = hubpoints->InsertNextPoint(verticalpoint[1].getX(), verticalpoint[1].getY(), verticalpoint[1].getZ());

                               vtkSmartPointer<vtkLine> medialsheetline = vtkSmartPointer<vtkLine>::New();
                               medialsheetline->GetPointIds()->SetId(0, id0);
                               medialsheetline->GetPointIds()->SetId(1, id1);

                               cellarraypointsline->InsertNextCell(medialsheetline);
                           }
                       }
                       break;
                   default:
                       break;
               }

               polypointsline->SetPoints(hubpoints);
               polypointsline->SetLines(cellarraypointsline);

               vtkSmartPointer<vtkPolyDataMapper> pointlinemapper = vtkSmartPointer<vtkPolyDataMapper>::New();
               pointlinemapper->SetInput(polypointsline);
               vtkSmartPointer<vtkActor>  lineactor = vtkActor::New();
               lineactor->SetMapper(pointlinemapper);
               lineactor->GetProperty()->SetLineWidth(1);
               lineactor->GetProperty()->SetColor(quadColor);// (0,1,0) is green, (1,1,1) is white.
               //lineactor->RotateY(95);
               //lineactor->RotateX(30);
               //lineactor->RotateY(-35);(35);
               renderer->AddActor(lineactor);

               horizonaledgepointsu.clear();
               horizonaledgepointsv.clear();
               verticaledgepointsu.clear();
               verticaledgepointsv.clear();
            }
        }

    //delete prim0; /*when use this, it will crash and throw Segmentation fault. Why? maybe M3DQuadPrimitive class has its own delete way?*/
    delete tpm;
}




/* Connect the diagnoal line of each (sub)quad, use the smaller area one.
*/
void visualization::connectDiagnoal(vtkSmartPointer< vtkPoints > hubpos){

    vtkSmartPointer<vtkCellArray> cellarraypointsline = vtkSmartPointer<vtkCellArray>::New();
    vtkSmartPointer<vtkPolyData> polypointsline = vtkSmartPointer<vtkPolyData>::New();

    double p0[3];
    double p1[3];
    double p2[3];
    double p3[3];

    for(int i=0; i<quadNum;i++){
        for(int j=0; j<subQuadPointsNum-(step+1);j++){// minus (step+1) because do not need the right column points.
           if((j+1)%(step+1)!=0){//do not need the bottom line's points.
               //Create four points (must be in counter clockwise order)
               hubpos->GetPoint(i*subQuadPointsNum+j, p0);
               hubpos->GetPoint(i*subQuadPointsNum+j+1, p1);
               hubpos->GetPoint(i*subQuadPointsNum+j+step+2, p2);
               hubpos->GetPoint(i*subQuadPointsNum+j+step+1, p3);

               // Find the diagnoal which make min area of the quad.
               Vector3D point[4];
               point[0].set(p0[0],p0[1],p0[2]);
               point[1].set(p1[0],p1[1],p1[2]);
               point[2].set(p2[0],p2[1],p2[2]);
               point[3].set(p3[0],p3[1],p3[2]);

               // If true, use diagonalp0p2; if false use diagonalp1p3
               bool p0p2 = tools.quadAreaMin(point);
               vtkSmartPointer<vtkLine> medialsheetline = vtkSmartPointer<vtkLine>::New();
               vtkIdType id0, id1;
               if(p0p2){ // Connect p0p2 as diagnoal
                   id0 = i*subQuadPointsNum+j;
                   id1 = i*subQuadPointsNum+j+step+2;
               }
               else{ // Connect p1p3 as diagnoal
                   id0 = i*subQuadPointsNum+j+1;
                   id1 = i*subQuadPointsNum+j+step+1;
               }

               medialsheetline->GetPointIds()->SetId(0, id0);
               medialsheetline->GetPointIds()->SetId(1, id1);

               cellarraypointsline->InsertNextCell(medialsheetline);
            }
         }
    }

    polypointsline->SetPoints(hubpos);
    polypointsline->SetLines(cellarraypointsline);

    vtkSmartPointer<vtkPolyDataMapper> pointlinemapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    pointlinemapper->SetInput(polypointsline);
    vtkSmartPointer<vtkActor>  lineactor = vtkActor::New();
    lineactor->SetMapper(pointlinemapper);
    lineactor->GetProperty()->SetLineWidth(1);
    lineactor->GetProperty()->SetColor(this->quadColor);// (0,1,0) is green, (1,1,1) is white.
    lineactor->RotateY(this->rY);
    lineactor->RotateX(this->rX);
    lineactor->RotateZ(this->rZ);
    lineactor->GetProperty()->SetOpacity(this->opacity);
    renderer->AddActor(lineactor);
}




/* Draw the quads with filled inside.
 * quadtype: 0(boundary quads), 1(skeletal quads)
 * side: 0(up), 1(down)
 * bool moved: true if draw the volume using the original r u p values in the input m3d file; false if draw the volume
 * using new values by adding deltaU and deltaV to primitive in the input m3d file.
*/
void visualization::drawQuads(vtkSmartPointer< vtkPoints > hubpos){

    vtkSmartPointer<vtkCellArray> quads = vtkSmartPointer<vtkCellArray>::New();

    for(int i=0; i<quadNum;i++){
        for(int j=0; j<subQuadPointsNum-(step+1);j++){// minus (step+1) because do not need the right column points.
           if((j+1)%(step+1)!=0){//do not need the bottom line's points.
                //define a up quad, and stroe it in a cellarray
                //Create four points (must be in counter clockwise order)
                vtkSmartPointer<vtkQuad> quad = vtkSmartPointer<vtkQuad>::New();
                quad->GetPointIds()->SetId(0, i*subQuadPointsNum+j);
                quad->GetPointIds()->SetId(1, i*subQuadPointsNum+j+1);
                quad->GetPointIds()->SetId(2, i*subQuadPointsNum+j+step+2);
                quad->GetPointIds()->SetId(3, i*subQuadPointsNum+j+step+1);

                quads->InsertNextCell(quad);//the size of the cellarray is the number of the quads, its 24.
            }
         }
    }

    // Create a polydata to store everything in
    vtkSmartPointer<vtkPolyData> polydataQuad = vtkSmartPointer<vtkPolyData>::New();

    // Add the points and quads to the dataset
    polydataQuad->SetPoints(hubpos);
    polydataQuad->SetPolys(quads);

    // Setup actor and mapper
    vtkSmartPointer<vtkPolyDataMapper> quadMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    #if VTK_MAJOR_VERSION <= 5
    quadMapper->SetInput(polydataQuad);
    #else
    quadMapper->SetInputData(polydataQuad);
    #endif

    vtkSmartPointer<vtkActor> quadActor = vtkSmartPointer<vtkActor>::New();
    quadActor->SetMapper(quadMapper);
    quadActor->RotateX(rX);
    quadActor->RotateY(rY);
    quadActor->RotateZ(rZ);
    quadActor->GetProperty()->SetOpacity(this->opacity);

    quadActor->GetProperty()->SetColor(quadColor);
    renderer->AddActor(quadActor);
}


/* Draw a specific quad with filled inside by its index number.
 * quadtype: 0(boundary quads), 1(skeletal quads)
 * side: 0(up), 1(down)
 * bool moved: true if draw the volume using the original r u p values in the input m3d file; false if draw the volume
 * using new values by adding deltaU and deltaV to primitive in the input m3d file.
*/
void visualization::drawQuadByIndex(vtkSmartPointer< vtkPoints > hubpos, int quadIndex){

    vtkSmartPointer<vtkCellArray> quads = vtkSmartPointer<vtkCellArray>::New();

    for(int j=0; j<subQuadPointsNum-(step+1);j++){// minus (step+1) because do not need the right column points.
        if((j+1)%(step+1)!=0){//do not need the bottom line's points.
            //define a up quad, and stroe it in a cellarray
             //Create four points (must be in counter clockwise order)
             vtkSmartPointer<vtkQuad> quad = vtkSmartPointer<vtkQuad>::New();
             quad->GetPointIds()->SetId(0, quadIndex*subQuadPointsNum+j);
             quad->GetPointIds()->SetId(1, quadIndex*subQuadPointsNum+j+1);
             quad->GetPointIds()->SetId(2, quadIndex*subQuadPointsNum+j+step+2);
             quad->GetPointIds()->SetId(3, quadIndex*subQuadPointsNum+j+step+1);

             quads->InsertNextCell(quad);//the size of the cellarray is the number of the quads, its 24.
         }
     }

    // Create a polydata to store everything in
    vtkSmartPointer<vtkPolyData> polydataQuad = vtkSmartPointer<vtkPolyData>::New();

    // Add the points and quads to the dataset
    polydataQuad->SetPoints(hubpos);
    polydataQuad->SetPolys(quads);

    // Setup actor and mapper
    vtkSmartPointer<vtkPolyDataMapper> quadMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    #if VTK_MAJOR_VERSION <= 5
    quadMapper->SetInput(polydataQuad);
    #else
    quadMapper->SetInputData(polydataQuad);
    #endif

    vtkSmartPointer<vtkActor> quadActor = vtkSmartPointer<vtkActor>::New();
    quadActor->SetMapper(quadMapper);

    quadActor->GetProperty()->SetColor(quadColor);
    renderer->AddActor(quadActor);
}


/* Draw a specific quad with filled inside by its index number.
 * quadtype: 0(boundary quads), 1(skeletal quads)
 * side: 0(up), 1(down)
*/
void visualization::drawQuadByIndexAndHubpos(vtkSmartPointer< vtkPoints > hubpos, int quadIndex){

    vtkSmartPointer<vtkCellArray> quads = vtkSmartPointer<vtkCellArray>::New();

    for(int j=0; j<subQuadPointsNum-(step+1);j++){// minus (step+1) because do not need the right column points.
        if((j+1)%(step+1)!=0){//do not need the bottom line's points.
            //define a up quad, and stroe it in a cellarray
             //Create four points (must be in counter clockwise order)
             vtkSmartPointer<vtkQuad> quad = vtkSmartPointer<vtkQuad>::New();
             quad->GetPointIds()->SetId(0, quadIndex*subQuadPointsNum+j);
             quad->GetPointIds()->SetId(1, quadIndex*subQuadPointsNum+j+1);
             quad->GetPointIds()->SetId(2, quadIndex*subQuadPointsNum+j+step+2);
             quad->GetPointIds()->SetId(3, quadIndex*subQuadPointsNum+j+step+1);

             quads->InsertNextCell(quad);//the size of the cellarray is the number of the quads, its 24.
         }
     }

    // Create a polydata to store everything in
    vtkSmartPointer<vtkPolyData> polydataQuad = vtkSmartPointer<vtkPolyData>::New();

    // Add the points and quads to the dataset
    polydataQuad->SetPoints(hubpos);
    polydataQuad->SetPolys(quads);

    // Setup actor and mapper
    vtkSmartPointer<vtkPolyDataMapper> quadMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    #if VTK_MAJOR_VERSION <= 5
    quadMapper->SetInput(polydataQuad);
    #else
    quadMapper->SetInputData(polydataQuad);
    #endif

    vtkSmartPointer<vtkActor> quadActor = vtkSmartPointer<vtkActor>::New();
    quadActor->SetMapper(quadMapper);
    quadActor->RotateX(rX);
    quadActor->RotateY(rY);
    quadActor->RotateZ(rZ);
    quadActor->GetProperty()->SetColor(quadColor);
    quadActor->GetProperty()->SetOpacity(0.5);
    renderer->AddActor(quadActor);
}


/* draw the volume of a specific quads between boundary and skeletal.
 * bool moved: true if draw the volume using the original r u p values in the input m3d file; false if draw the volume
 * using new values by adding deltaU and deltaV to primitive in the input m3d file.
*/
void visualization::drawVolume(vtkSmartPointer< vtkPoints > boundaryHubpos, vtkSmartPointer< vtkPoints > skeletalHubpos, int quadIndex){

    // After divided, each quad is divided into step*step sub-quad, which contains (step+1)*(step+1) points.
    vtkSmartPointer< vtkPoints > hubpos = vtkSmartPointer< vtkPoints >::New();
    vtkSmartPointer<vtkCellArray> quads = vtkSmartPointer<vtkCellArray>::New();

    // Add the left-side-quad of this quad to the list.
    for(int n=0; n<step;n++){
        cout<<"n is: "<<n<<endl;
        //Get the 25 sub-points of this quad.
        double p[3];
        vtkSmartPointer<vtkQuad> quad = vtkSmartPointer<vtkQuad>::New();

        boundaryHubpos->GetPoint(quadIndex*subQuadPointsNum+n,p);
        hubpos->InsertNextPoint(p[0],p[1],p[2]);
        quad->GetPointIds()->SetId(0, 0+4*n);

        skeletalHubpos->GetPoint(quadIndex*subQuadPointsNum+n,p);
        hubpos->InsertNextPoint(p[0],p[1],p[2]);
        quad->GetPointIds()->SetId(1, 1+4*n);

        skeletalHubpos->GetPoint(quadIndex*subQuadPointsNum+n+1,p);
        hubpos->InsertNextPoint(p[0],p[1],p[2]);
        quad->GetPointIds()->SetId(2, 2+4*n);

        boundaryHubpos->GetPoint(quadIndex*subQuadPointsNum+n+1,p);
        hubpos->InsertNextPoint(p[0],p[1],p[2]);
        quad->GetPointIds()->SetId(3, 3+4*n);

        quads->InsertNextCell(quad);
    }

    // Add the right-side-quad of this quad to the list.
    int pointIndex = 0;
    for(int n=subQuadPointsNum-(step+1); n<subQuadPointsNum-1;n++){
        cout<<"n is: "<<n<<endl;
        double p[3];
        vtkSmartPointer<vtkQuad> quad = vtkSmartPointer<vtkQuad>::New();

        boundaryHubpos->GetPoint(quadIndex*subQuadPointsNum+n,p);
        hubpos->InsertNextPoint(p[0],p[1],p[2]);
        quad->GetPointIds()->SetId(0, 0+4*(step+pointIndex));

        skeletalHubpos->GetPoint(quadIndex*subQuadPointsNum+n,p);
        hubpos->InsertNextPoint(p[0],p[1],p[2]);
        quad->GetPointIds()->SetId(1, 1+4*(step+pointIndex));

        skeletalHubpos->GetPoint(quadIndex*subQuadPointsNum+n+1,p);
        hubpos->InsertNextPoint(p[0],p[1],p[2]);
        quad->GetPointIds()->SetId(2, 2+4*(step+pointIndex));

        boundaryHubpos->GetPoint(quadIndex*subQuadPointsNum+n+1,p);
        hubpos->InsertNextPoint(p[0],p[1],p[2]);
        quad->GetPointIds()->SetId(3, 3+4*(step+pointIndex));

        quads->InsertNextCell(quad);
        pointIndex++;
    }

    // Add the bottom-side-quad of this quad to the list.
    int nn =0;
    for(int pointIndex = 0; pointIndex<step&&nn<subQuadPointsNum-1;pointIndex++){
        nn=step+(step+1)*pointIndex;
        cout<<"n is: "<<nn<<endl;
        double p[3];
        vtkSmartPointer<vtkQuad> quad = vtkSmartPointer<vtkQuad>::New();

        boundaryHubpos->GetPoint(quadIndex*subQuadPointsNum+nn,p);
        hubpos->InsertNextPoint(p[0],p[1],p[2]);
        quad->GetPointIds()->SetId(0, 0+4*(2*step+pointIndex));

        skeletalHubpos->GetPoint(quadIndex*subQuadPointsNum+nn,p);
        hubpos->InsertNextPoint(p[0],p[1],p[2]);
        quad->GetPointIds()->SetId(1, 1+4*(2*step+pointIndex));

        skeletalHubpos->GetPoint(quadIndex*subQuadPointsNum+nn+(step+1),p);
        hubpos->InsertNextPoint(p[0],p[1],p[2]);
        quad->GetPointIds()->SetId(2, 2+4*(2*step+pointIndex));

        boundaryHubpos->GetPoint(quadIndex*subQuadPointsNum+nn+(step+1),p);
        hubpos->InsertNextPoint(p[0],p[1],p[2]);
        quad->GetPointIds()->SetId(3, 3+4*(2*step+pointIndex));

        quads->InsertNextCell(quad);
    }

    // Add the top-side-quad of this quad to the list.
    nn =0;
    for(int pointIndex = 0; pointIndex<step&&nn<subQuadPointsNum-(step+1);pointIndex++){
        nn=0+(step+1)*pointIndex;
        cout<<"n is: "<<nn<<endl;
        double p[3];
        vtkSmartPointer<vtkQuad> quad = vtkSmartPointer<vtkQuad>::New();

        boundaryHubpos->GetPoint(quadIndex*subQuadPointsNum+nn,p);
        hubpos->InsertNextPoint(p[0],p[1],p[2]);
        quad->GetPointIds()->SetId(0, 0+4*(3*step+pointIndex));

        skeletalHubpos->GetPoint(quadIndex*subQuadPointsNum+nn,p);
        hubpos->InsertNextPoint(p[0],p[1],p[2]);
        quad->GetPointIds()->SetId(1, 1+4*(3*step+pointIndex));

        skeletalHubpos->GetPoint(quadIndex*subQuadPointsNum+nn+(step+1),p);
        hubpos->InsertNextPoint(p[0],p[1],p[2]);
        quad->GetPointIds()->SetId(2, 2+4*(3*step+pointIndex));

        boundaryHubpos->GetPoint(quadIndex*subQuadPointsNum+nn+(step+1),p);
        hubpos->InsertNextPoint(p[0],p[1],p[2]);
        quad->GetPointIds()->SetId(3, 3+4*(3*step+pointIndex));

        quads->InsertNextCell(quad);
    }

    // Create a polydata to store everything in
    vtkSmartPointer<vtkPolyData> polydataQuad = vtkSmartPointer<vtkPolyData>::New();

    // Add the points and quads to the dataset
    polydataQuad->SetPoints(hubpos);
    polydataQuad->SetPolys(quads);

    // Setup actor and mapper
    vtkSmartPointer<vtkPolyDataMapper> quadMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    #if VTK_MAJOR_VERSION <= 5
    quadMapper->SetInput(polydataQuad);
    #else
    quadMapper->SetInputData(polydataQuad);
    #endif

    vtkSmartPointer<vtkActor> quadActor = vtkSmartPointer<vtkActor>::New();
    quadActor->SetMapper(quadMapper);
    quadActor->RotateX(rX);
    quadActor->RotateY(rY);
    quadActor->RotateZ(rZ);
    quadActor->GetProperty()->SetColor(quadColor);
    quadActor->GetProperty()->SetOpacity(0.5);
    renderer->AddActor(quadActor);


    // Draw the quadIndex's boundary quads.
    drawQuadByIndexAndHubpos(boundaryHubpos, quadIndex);

    // Draw the quadIndex's skeletal quad.
    drawQuadByIndexAndHubpos(skeletalHubpos, quadIndex);
}



/* Draw the subquad volume by quad and subquad index.
 * subquad index is: 0, 1, 2, ..., step*step
*/
void visualization::drawSubQuadVolume(vtkSmartPointer< vtkPoints > boundaryHubpos, vtkSmartPointer< vtkPoints > skeletalHubpos, int quadIndex,
                               int subQuadIndex){
    // Given subQuadIndex, the quad's top-left corner point is in  (subQuadIndex%step) row.
    int m = subQuadIndex % step;

    // The point in column
    int n = (subQuadIndex-m)/step;

    int pindex = n*(step+1) + m;

    double p[3];

    // After divided, each quad is divided into step*step sub-quad, which contains (step+1)*(step+1) points.
    vtkSmartPointer< vtkPoints > hubpos = vtkSmartPointer< vtkPoints >::New();
    vtkSmartPointer<vtkCellArray> quads = vtkSmartPointer<vtkCellArray>::New();

    // Get the four corner of this sub-quad
    boundaryHubpos->GetPoint(quadIndex*subQuadPointsNum+pindex,p); // up-left
    vtkIdType id1_b = hubpos->InsertNextPoint(p[0],p[1],p[2]);

    skeletalHubpos->GetPoint(quadIndex*subQuadPointsNum+pindex,p);
    vtkIdType id1_s = hubpos->InsertNextPoint(p[0],p[1],p[2]);

    boundaryHubpos->GetPoint(quadIndex*subQuadPointsNum+pindex + step+1,p); // up-right
    vtkIdType id2_b = hubpos->InsertNextPoint(p[0],p[1],p[2]);

    skeletalHubpos->GetPoint(quadIndex*subQuadPointsNum+pindex + step+1,p);
    vtkIdType id2_s = hubpos->InsertNextPoint(p[0],p[1],p[2]);

    boundaryHubpos->GetPoint(quadIndex*subQuadPointsNum+pindex+1,p); // bottom-left
    vtkIdType id3_b = hubpos->InsertNextPoint(p[0],p[1],p[2]);

    skeletalHubpos->GetPoint(quadIndex*subQuadPointsNum+pindex+1,p);
    vtkIdType id3_s = hubpos->InsertNextPoint(p[0],p[1],p[2]);

    boundaryHubpos->GetPoint(quadIndex*subQuadPointsNum+pindex+1 + step+1,p); // bottom-right
    vtkIdType id4_b = hubpos->InsertNextPoint(p[0],p[1],p[2]);

    skeletalHubpos->GetPoint(quadIndex*subQuadPointsNum+pindex+1 + step+1,p);
    vtkIdType id4_s = hubpos->InsertNextPoint(p[0],p[1],p[2]);


    // Add the left-side-quad of this quad to the list.
    vtkSmartPointer<vtkQuad> quad_left = vtkSmartPointer<vtkQuad>::New();
    quad_left->GetPointIds()->SetId(0, id1_b);
    quad_left->GetPointIds()->SetId(1, id1_s);
    quad_left->GetPointIds()->SetId(2, id3_s);
    quad_left->GetPointIds()->SetId(3, id3_b);

    quads->InsertNextCell(quad_left);


    // Add the right-side-quad of this quad to the list.
    vtkSmartPointer<vtkQuad> quad_right = vtkSmartPointer<vtkQuad>::New();
    quad_right->GetPointIds()->SetId(0, id2_b);
    quad_right->GetPointIds()->SetId(1, id2_s);
    quad_right->GetPointIds()->SetId(2, id4_s);
    quad_right->GetPointIds()->SetId(3, id4_b);

    quads->InsertNextCell(quad_right);


    // Add the top-side-quad of this quad to the list.
    vtkSmartPointer<vtkQuad> quad_top = vtkSmartPointer<vtkQuad>::New();
    quad_top->GetPointIds()->SetId(0, id2_b);
    quad_top->GetPointIds()->SetId(1, id2_s);
    quad_top->GetPointIds()->SetId(2, id1_s);
    quad_top->GetPointIds()->SetId(3, id1_b);

    quads->InsertNextCell(quad_top);


    // Add the bottom-side-quad of this quad to the list.
    vtkSmartPointer<vtkQuad> quad_bottom = vtkSmartPointer<vtkQuad>::New();
    quad_bottom->GetPointIds()->SetId(0, id4_b);
    quad_bottom->GetPointIds()->SetId(1, id4_s);
    quad_bottom->GetPointIds()->SetId(2, id3_s);
    quad_bottom->GetPointIds()->SetId(3, id3_b);

    quads->InsertNextCell(quad_bottom);


    // Add the top quad
    vtkSmartPointer<vtkQuad> top = vtkSmartPointer<vtkQuad>::New();
    top->GetPointIds()->SetId(0, id4_b);
    top->GetPointIds()->SetId(1, id2_b);
    top->GetPointIds()->SetId(2, id1_b);
    top->GetPointIds()->SetId(3, id3_b);

    quads->InsertNextCell(top);


    // Add the bottom quad
    vtkSmartPointer<vtkQuad> bottom = vtkSmartPointer<vtkQuad>::New();
    bottom->GetPointIds()->SetId(0, id4_s);
    bottom->GetPointIds()->SetId(1, id2_s);
    bottom->GetPointIds()->SetId(2, id1_s);
    bottom->GetPointIds()->SetId(3, id3_s);

    quads->InsertNextCell(bottom);


    // Create a polydata to store everything in
    vtkSmartPointer<vtkPolyData> polydataQuad = vtkSmartPointer<vtkPolyData>::New();

    // Add the points and quads to the dataset
    polydataQuad->SetPoints(hubpos);
    polydataQuad->SetPolys(quads);

    // Setup actor and mapper
    vtkSmartPointer<vtkPolyDataMapper> quadMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    #if VTK_MAJOR_VERSION <= 5
    quadMapper->SetInput(polydataQuad);
    #else
    quadMapper->SetInputData(polydataQuad);
    #endif

    vtkSmartPointer<vtkActor> quadActor = vtkSmartPointer<vtkActor>::New();
    quadActor->SetMapper(quadMapper);
    quadActor->RotateX(rX);
    quadActor->RotateY(rY);
    quadActor->RotateZ(rZ);
    quadActor->GetProperty()->SetColor(quadColor);
    quadActor->GetProperty()->SetOpacity(0.5);
    renderer->AddActor(quadActor);
}



/* Draw the edge of each quad.
 * This function do the same thing as drawEdgeOfQuads.
*/
void visualization::drawEdgeOfQuads_method2(vtkSmartPointer< vtkPoints > hubpos){

    vtkSmartPointer<vtkCellArray> cellarraypointsline = vtkSmartPointer<vtkCellArray>::New();
    vtkSmartPointer<vtkPolyData> polypointsline = vtkSmartPointer<vtkPolyData>::New();

    for(int i=0; i<quadNum;i++){
        // Add the left side of each quad to the line-list.
        for(int k=0; k<step;k++){
            double p[3];
            vtkSmartPointer<vtkLine> medialsheetline = vtkSmartPointer<vtkLine>::New();

            hubpos->GetPoint(i*subQuadPointsNum+k,p);
            vtkIdType id0 = hubpos->InsertNextPoint(p[0],p[1],p[2]);

            hubpos->GetPoint(i*subQuadPointsNum+k+1,p);
            vtkIdType id1 = hubpos->InsertNextPoint(p[0],p[1],p[2]);

            medialsheetline->GetPointIds()->SetId(0, id0);
            medialsheetline->GetPointIds()->SetId(1, id1);

            cellarraypointsline->InsertNextCell(medialsheetline);
        }

        // Add the right side of each quad to the line-list.
        for(int k=subQuadPointsNum-(step+1); k<subQuadPointsNum-1;k++){
            double p[3];
            vtkSmartPointer<vtkLine> medialsheetline = vtkSmartPointer<vtkLine>::New();

            hubpos->GetPoint(i*subQuadPointsNum+k,p);
            vtkIdType id0 = hubpos->InsertNextPoint(p[0],p[1],p[2]);

            hubpos->GetPoint(i*subQuadPointsNum+k+1,p);
            vtkIdType id1 = hubpos->InsertNextPoint(p[0],p[1],p[2]);

            medialsheetline->GetPointIds()->SetId(0, id0);
            medialsheetline->GetPointIds()->SetId(1, id1);

            cellarraypointsline->InsertNextCell(medialsheetline);
        }

        // Add the bottom side of each quad to the line-list.
        int pointIndex =0;
        for(int k = 0; k<step&&pointIndex<subQuadPointsNum-1;k++){
            pointIndex = step + k*(step+1);
            double p[3];
            vtkSmartPointer<vtkLine> medialsheetline = vtkSmartPointer<vtkLine>::New();

            hubpos->GetPoint(i*subQuadPointsNum+pointIndex,p);
            vtkIdType id0 = hubpos->InsertNextPoint(p[0],p[1],p[2]);

            hubpos->GetPoint(i*subQuadPointsNum+pointIndex+step+1,p);
            vtkIdType id1 = hubpos->InsertNextPoint(p[0],p[1],p[2]);

            medialsheetline->GetPointIds()->SetId(0, id0);
            medialsheetline->GetPointIds()->SetId(1, id1);

            cellarraypointsline->InsertNextCell(medialsheetline);
        }

        // Add the top side of each quad to the line-list.
        pointIndex =0;
        for(int k = 0; k<step&&pointIndex<subQuadPointsNum-1;k++){
            pointIndex = 0 + k*(step+1);
            double p[3];
            vtkSmartPointer<vtkLine> medialsheetline = vtkSmartPointer<vtkLine>::New();

            hubpos->GetPoint(i*subQuadPointsNum+pointIndex,p);
            vtkIdType id0 = hubpos->InsertNextPoint(p[0],p[1],p[2]);

            hubpos->GetPoint(i*subQuadPointsNum+pointIndex+step+1,p);
            vtkIdType id1 = hubpos->InsertNextPoint(p[0],p[1],p[2]);

            medialsheetline->GetPointIds()->SetId(0, id0);
            medialsheetline->GetPointIds()->SetId(1, id1);

            cellarraypointsline->InsertNextCell(medialsheetline);
        }
    }

    polypointsline->SetPoints(hubpos);
    polypointsline->SetLines(cellarraypointsline);

    vtkSmartPointer<vtkPolyDataMapper> pointlinemapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    pointlinemapper->SetInput(polypointsline);
    vtkSmartPointer<vtkActor>  lineactor = vtkActor::New();
    lineactor->SetMapper(pointlinemapper);
    lineactor->GetProperty()->SetLineWidth(1);
    lineactor->GetProperty()->SetColor(quadColor);// (0,1,0) is green, (1,1,1) is white.
    renderer->AddActor(lineactor);
}





/* Draw the frame of each quad.
 * So we can see each quad was divided into how many sub-quad.
*/
void visualization::drawFrameOfQuads(vtkSmartPointer< vtkPoints > hubpos){

    vtkSmartPointer<vtkCellArray> cellarraypointsline = vtkSmartPointer<vtkCellArray>::New();
    vtkSmartPointer<vtkPolyData> polypointsline = vtkSmartPointer<vtkPolyData>::New();

    vtkSmartPointer< vtkPoints > hubPosNewSequence = vtkSmartPointer< vtkPoints >::New();

    //for each quad
    for(int i=0; i<quadNum;i++){
        //for each column's (step+1) points.
        //draw vertical lines.
        for(int n=0; n<step+1;n++){
            int pointIndex =0;
            for(int k = 0; k<step&&pointIndex<subQuadPointsNum-1;k++){
                pointIndex = n + k*(step+1);
                double p[3];
                vtkSmartPointer<vtkLine> medialsheetline = vtkSmartPointer<vtkLine>::New();

                hubpos->GetPoint(i*subQuadPointsNum+pointIndex,p);
                vtkIdType id0 = hubPosNewSequence->InsertNextPoint(p[0],p[1],p[2]);
                //cout<<"p1: "<<p[0]<<"  "<<p[1]<<"  "<<p[2]<<endl;

                hubpos->GetPoint(i*subQuadPointsNum+pointIndex+step+1,p);
                vtkIdType id1 = hubPosNewSequence->InsertNextPoint(p[0],p[1],p[2]);
                //cout<<"p2: "<<p[0]<<"  "<<p[1]<<"  "<<p[2]<<endl;

                medialsheetline->GetPointIds()->SetId(0, id0);
                medialsheetline->GetPointIds()->SetId(1, id1);

                cellarraypointsline->InsertNextCell(medialsheetline);
            }
        }

        //draw horizonal lines.
        for(int n=0; n<step+1;n++){
            for(int k = 0; k<step;k++){
                int currentpoint = k + n*(step+1);
                double p[3];
                vtkSmartPointer<vtkLine> medialsheetline = vtkSmartPointer<vtkLine>::New();

                hubpos->GetPoint(i*subQuadPointsNum+currentpoint,p);
                vtkIdType id0 = hubPosNewSequence->InsertNextPoint(p[0],p[1],p[2]);

                hubpos->GetPoint(i*subQuadPointsNum+currentpoint+1,p);
                vtkIdType id1 = hubPosNewSequence->InsertNextPoint(p[0],p[1],p[2]);

                medialsheetline->GetPointIds()->SetId(0, id0);
                medialsheetline->GetPointIds()->SetId(1, id1);

                cellarraypointsline->InsertNextCell(medialsheetline);
            }
        }
    }

    polypointsline->SetPoints(hubPosNewSequence);
    polypointsline->SetLines(cellarraypointsline);

    vtkSmartPointer<vtkPolyDataMapper> pointlinemapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    pointlinemapper->SetInput(polypointsline);
    vtkSmartPointer<vtkActor>  lineactor = vtkActor::New();
    lineactor->SetMapper(pointlinemapper);
    lineactor->GetProperty()->SetLineWidth(1);
    lineactor->GetProperty()->SetColor(this->quadColor);// (0,1,0) is green, (1,1,1) is white.
    lineactor->RotateX(rX);
    lineactor->RotateY(rY);
    lineactor->RotateZ(rZ);
    lineactor->GetProperty()->SetOpacity(this->opacity);
    renderer->AddActor(lineactor);
}






/* Draw the frame of quad by quad index.
*/
void visualization::drawFrameOfQuadByIndex(vtkSmartPointer< vtkPoints > hubpos, int quadIndex){

    vtkSmartPointer<vtkCellArray> cellarraypointsline = vtkSmartPointer<vtkCellArray>::New();
    vtkSmartPointer<vtkPolyData> polypointsline = vtkSmartPointer<vtkPolyData>::New();

    vtkSmartPointer< vtkPoints > hubPosNewSequence = vtkSmartPointer< vtkPoints >::New();

    //draw vertical lines.
    for(int n=0; n<step+1;n++){
        int pointIndex =0;
        for(int k = 0; k<step&&pointIndex<subQuadPointsNum-1;k++){
            pointIndex = n + k*(step+1);
            double p[3];
            vtkSmartPointer<vtkLine> medialsheetline = vtkSmartPointer<vtkLine>::New();

            hubpos->GetPoint(quadIndex*subQuadPointsNum+pointIndex,p);
            vtkIdType id0 = hubPosNewSequence->InsertNextPoint(p[0],p[1],p[2]);
            cout<<"----vertical lines-----p1: "<<p[0]<<"  "<<p[1]<<"  "<<p[2]<<endl;

            hubpos->GetPoint(quadIndex*subQuadPointsNum+pointIndex+step+1,p);
            vtkIdType id1 = hubPosNewSequence->InsertNextPoint(p[0],p[1],p[2]);
            cout<<"----vertical lines-----p2: "<<p[0]<<"  "<<p[1]<<"  "<<p[2]<<endl;

            medialsheetline->GetPointIds()->SetId(0, id0);
            medialsheetline->GetPointIds()->SetId(1, id1);

            cellarraypointsline->InsertNextCell(medialsheetline);
        }
    }

    //draw horizonal lines.
    for(int n=0; n<step+1;n++){
        for(int k = 0; k<step;k++){
            int currentpoint = k + n*(step+1);
            double p[3];
            vtkSmartPointer<vtkLine> medialsheetline = vtkSmartPointer<vtkLine>::New();

            hubpos->GetPoint(quadIndex*subQuadPointsNum+currentpoint,p);
            vtkIdType id0 = hubPosNewSequence->InsertNextPoint(p[0],p[1],p[2]);
            cout<<"----horizonal lines-----p1: "<<p[0]<<"  "<<p[1]<<"  "<<p[2]<<endl;

            hubpos->GetPoint(quadIndex*subQuadPointsNum+currentpoint+1,p);
            vtkIdType id1 = hubPosNewSequence->InsertNextPoint(p[0],p[1],p[2]);
            cout<<"----horizonal lines-----p1: "<<p[0]<<"  "<<p[1]<<"  "<<p[2]<<endl;

            medialsheetline->GetPointIds()->SetId(0, id0);
            medialsheetline->GetPointIds()->SetId(1, id1);

            cellarraypointsline->InsertNextCell(medialsheetline);
        }
    }

    polypointsline->SetPoints(hubPosNewSequence);
    polypointsline->SetLines(cellarraypointsline);

    vtkSmartPointer<vtkPolyDataMapper> pointlinemapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    pointlinemapper->SetInput(polypointsline);
    vtkSmartPointer<vtkActor>  lineactor = vtkActor::New();
    lineactor->SetMapper(pointlinemapper);
    lineactor->GetProperty()->SetLineWidth(1);
    lineactor->GetProperty()->SetColor(quadColor);// (0,1,0) is green, (1,1,1) is white.
    lineactor->RotateX(rX);
    lineactor->RotateY(rY);
    lineactor->RotateZ(rZ);
    //lineactor->GetProperty()->SetOpacity(0.5);
    renderer->AddActor(lineactor);
}



/* draw the spoke by its quadIndex.
 * bool moved: true if draw the spoke using the original r u p values in the input m3d file; false if draw the spoke
 * using new values by adding deltaU and deltaV to primitive in the input m3d file.
 * side: 0(up spokes); 1(down spokes).
*/
void visualization::drawSpoke(vtkSmartPointer< vtkPoints > boundaryHubpos, vtkSmartPointer< vtkPoints > skeletalHubpos, int quadIndex){

    vtkSmartPointer< vtkPoints > hubpos = vtkSmartPointer< vtkPoints >::New();
    vtkSmartPointer<vtkCellArray> cellarraypointsline = vtkSmartPointer<vtkCellArray>::New();
    vtkSmartPointer<vtkPolyData> polypointsline = vtkSmartPointer<vtkPolyData>::New();

    for(int n=0; n<subQuadPointsNum;n++){
        double p[3];
        vtkSmartPointer<vtkLine> medialsheetline = vtkSmartPointer<vtkLine>::New();

        skeletalHubpos->GetPoint(quadIndex*subQuadPointsNum+n,p);
        vtkIdType id0 = hubpos->InsertNextPoint(p[0],p[1],p[2]);

        boundaryHubpos->GetPoint(quadIndex*subQuadPointsNum+n,p);
        vtkIdType id1 = hubpos->InsertNextPoint(p[0],p[1],p[2]);

        medialsheetline->GetPointIds()->SetId(0, id0);
        medialsheetline->GetPointIds()->SetId(1, id1);

        cellarraypointsline->InsertNextCell(medialsheetline);
    }

    polypointsline->SetPoints(hubpos);
    polypointsline->SetLines(cellarraypointsline);

    vtkSmartPointer<vtkPolyDataMapper> pointlinemapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    pointlinemapper->SetInput(polypointsline);
    vtkSmartPointer<vtkActor>  lineactor = vtkActor::New();
    lineactor->SetMapper(pointlinemapper);
    lineactor->GetProperty()->SetLineWidth(2);
    lineactor->GetProperty()->SetColor(quadColor);// (0,1,0) is green, (1,1,1) is white.
    lineactor->RotateX(rX);
    lineactor->RotateY(rY);
    lineactor->RotateZ(rZ);
    //lineactor->GetProperty()->SetOpacity(0.8);
    renderer->AddActor(lineactor);
}


///* draw the spoke by its quadIndex.
// * bool moved: true if draw the spoke using the original r u p values in the input m3d file; false if draw the spoke
// * using new values by adding deltaU and deltaV to primitive in the input m3d file.
// * side: 0(up spokes); 1(down spokes).
//*/
void visualization::drawSpoke_2(vtkSmartPointer< vtkPoints > boundaryHubpos, vtkSmartPointer< vtkPoints > skeletalHubpos, int quadNums,
                                int linewith, bool showspoketail){

    vtkSmartPointer< vtkPoints > hubpos = vtkSmartPointer< vtkPoints >::New();
    vtkSmartPointer<vtkCellArray> cellarraypointsline = vtkSmartPointer<vtkCellArray>::New();
    vtkSmartPointer<vtkPolyData> polypointsline = vtkSmartPointer<vtkPolyData>::New();

    for(unsigned int i = 0; i < quadNums; i++){
        for(int n=0; n<subQuadPointsNum;n++){
            double p[3];
            vtkSmartPointer<vtkLine> medialsheetline = vtkSmartPointer<vtkLine>::New();

            skeletalHubpos->GetPoint(i*subQuadPointsNum+n,p);
            vtkIdType id0 = hubpos->InsertNextPoint(p[0],p[1],p[2]);

            boundaryHubpos->GetPoint(i*subQuadPointsNum+n,p);
            vtkIdType id1 = hubpos->InsertNextPoint(p[0],p[1],p[2]);

            medialsheetline->GetPointIds()->SetId(0, id0);
            medialsheetline->GetPointIds()->SetId(1, id1);

            cellarraypointsline->InsertNextCell(medialsheetline);
        }
    }

    polypointsline->SetPoints(hubpos);
    polypointsline->SetLines(cellarraypointsline);

/*    vtkSmartPointer<vtkVertexGlyphFilter> vertexFilter = vtkSmartPointer<vtkVertexGlyphFilter>::New();
    #if VTK_MAJOR_VERSION <= 5
    vertexFilter->SetInputConnection(polypointsline->GetProducerPort());
    #else
    vertexFilter->SetInputData(polypointsline);
    #endif
    vertexFilter->Update();

    vtkSmartPointer<vtkPolyData> polydata = vtkSmartPointer<vtkPolyData>::New();
    polydata->ShallowCopy(vertexFilter->GetOutput());

    double bounds[6];
    polydata->GetBounds(bounds);

    // Find min and max z
    double minz = bounds[4];
    double maxz = bounds[5];

    // Create the color map
    vtkSmartPointer<vtkLookupTable> colorLookupTable = vtkSmartPointer<vtkLookupTable>::New();
    colorLookupTable->SetTableRange(minz, maxz);
    colorLookupTable->Build();

    // Generate the colors for each point based on the color map
    vtkSmartPointer<vtkUnsignedCharArray> colors = vtkSmartPointer<vtkUnsignedCharArray>::New();
    colors->SetNumberOfComponents(3);
    colors->SetName("Colors");

    for(int i = 0; i < polypointsline->GetNumberOfPoints(); i++) {
        double p[3];
        polypointsline->GetPoint(i,p);

        // Get a color correspondence to the z coordinate value from the color table
        double dcolor[3];
        colorLookupTable->GetColor(p[1], dcolor); //p[2] is the z coordinate

        unsigned char color[3];
        for(unsigned int j = 0; j < 3; j++) {
            color[j] = static_cast<unsigned char>(255.0 * dcolor[j]);
        }

        colors->InsertNextTupleValue(color);
    }




    polypointsline->GetPointData()->SetScalars(colors);*/

    vtkSmartPointer<vtkPolyDataMapper> pointlinemapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    pointlinemapper->SetInput(polypointsline);
    vtkSmartPointer<vtkActor>  lineactor = vtkActor::New();
    lineactor->SetMapper(pointlinemapper);
    lineactor->GetProperty()->SetLineWidth(linewith);
    lineactor->GetProperty()->SetColor(quadColor);// (0,1,0) is green, (1,1,1) is white.
    lineactor->RotateX(rX);
    lineactor->RotateY(rY);
    lineactor->RotateZ(rZ);
    //lineactor->GetProperty()->SetOpacity(0.2);

    renderer->AddActor(lineactor);


    // If choosen to show the spoke tail
    if(showspoketail) {
        showPointBall(skeletalHubpos, quadColor, 0.0016, this->renderer, this->rX, this->rY, this->rZ, this->opacity);
    }
}


void visualization::showPointBall(vtkSmartPointer< vtkPoints > points, double * color, double pointSize, vtkSmartPointer<vtkRenderer> renderer,
                                  double rX, double rY, double rZ, double opacity){
    // loop each point
    double p[3] = {0.0, 0.0, 0.0};

    for(unsigned int i = 0; i < points->GetNumberOfPoints(); i++){
        points->GetPoint(i, p);

        // Create a sphere
        vtkSmartPointer<vtkSphereSource> sphereSource = vtkSmartPointer<vtkSphereSource>::New();
        sphereSource->SetCenter(p[0], p[1], p[2]);
        sphereSource->SetRadius(pointSize);

        vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
        mapper->SetInputConnection(sphereSource->GetOutputPort());

        vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
        actor->SetMapper(mapper);
        actor->GetProperty()->SetColor(color);
        actor->RotateX(rX);
        actor->RotateY(rY);
        actor->RotateZ(rZ);
        actor->GetProperty()->SetOpacity(opacity);
        renderer->AddActor(actor);
    }
}






/* draw the spoke of all the quad at one time.
 * bool moved: true if draw the spoke using the original r u p values in the input m3d file; false if draw the spoke
 * using new values by adding deltaU and deltaV to primitive in the input m3d file.
 * side: 0(up spokes); 1(down spokes).
*/
void visualization::drawSpoke_3(vtkSmartPointer< vtkPoints > boundaryHubpos, vtkSmartPointer< vtkPoints > skeletalHubpos, int quadNums){

    vtkSmartPointer< vtkPoints > hubpos = vtkSmartPointer< vtkPoints >::New();
    vtkSmartPointer<vtkCellArray> cellarraypointsline = vtkSmartPointer<vtkCellArray>::New();
    vtkSmartPointer<vtkPolyData> linepolydata = vtkSmartPointer<vtkPolyData>::New();

    int lineNum = 0;
    for(unsigned int i = 0; i < quadNums; i++){
        for(int n=0; n<subQuadPointsNum;n++){
            double p[3];
            vtkSmartPointer<vtkLine> medialsheetline = vtkSmartPointer<vtkLine>::New();

            skeletalHubpos->GetPoint(i*subQuadPointsNum+n,p);
            vtkIdType id0 = hubpos->InsertNextPoint(p[0],p[1],p[2]);

            boundaryHubpos->GetPoint(i*subQuadPointsNum+n,p);
            vtkIdType id1 = hubpos->InsertNextPoint(p[0],p[1],p[2]);

            medialsheetline->GetPointIds()->SetId(0, id0);
            medialsheetline->GetPointIds()->SetId(1, id1);

            cellarraypointsline->InsertNextCell(medialsheetline);
            lineNum++;
        }
    }
    cout<<"-----lineNum: "<<lineNum<<endl; // shoulde be 39 while level 0

    linepolydata->SetPoints(hubpos);
    linepolydata->SetLines(cellarraypointsline);

    cout<<"----line numbers---"<<linepolydata->GetNumberOfPoints()<<endl;

    // Create a lookup table to map cell data to colors
     vtkSmartPointer<vtkLookupTable> lut = vtkSmartPointer<vtkLookupTable>::New();
     int tableSize = 39;//std::max(3*3 + 1, 10);
     lut->SetNumberOfTableValues(tableSize);
     lut->Build();

     // Fill in a few known colors, the rest will be generated if needed
     lut->SetTableValue(0     , 0     , 0     , 0, 1);  //Black
     lut->SetTableValue(1, 0.8900, 0.8100, 0.3400, 1); // Banana
     lut->SetTableValue(2, 1.0000, 0.3882, 0.2784, 1); // Tomato
     lut->SetTableValue(3, 0.9608, 0.8706, 0.7020, 1); // Wheat
     lut->SetTableValue(4, 0.9020, 0.9020, 0.9804, 1); // Lavender
     lut->SetTableValue(5, 1.0000, 0.4900, 0.2500, 1); // Flesh
     lut->SetTableValue(6, 0.5300, 0.1500, 0.3400, 1); // Raspberry
     lut->SetTableValue(7, 0.9804, 0.5020, 0.4471, 1); // Salmon
     lut->SetTableValue(8, 0.7400, 0.9900, 0.7900, 1); // Mint
     lut->SetTableValue(9, 0.2000, 0.6300, 0.7900, 1); // Peacock
     lut->SetTableValue(10     , 0     , 0     , 0, 1);  //Black
     lut->SetTableValue(11, 0.8900, 0.8100, 0.3400, 1); // Banana
     lut->SetTableValue(12, 1.0000, 0.3882, 0.2784, 1); // Tomato
     lut->SetTableValue(13, 0.9608, 0.8706, 0.7020, 1); // Wheat
     lut->SetTableValue(14, 0.9020, 0.9020, 0.9804, 1); // Lavender
     lut->SetTableValue(15, 1.0000, 0.4900, 0.2500, 1); // Flesh
     lut->SetTableValue(16, 0.5300, 0.1500, 0.3400, 1); // Raspberry
     lut->SetTableValue(17, 0.9804, 0.5020, 0.4471, 1); // Salmon
     lut->SetTableValue(18, 0.7400, 0.9900, 0.7900, 1); // Mint
     lut->SetTableValue(19, 0.2000, 0.6300, 0.7900, 1); // Peacock
     lut->SetTableValue(20     , 0     , 0     , 0, 1);  //Black
     lut->SetTableValue(21, 0.8900, 0.8100, 0.3400, 1); // Banana
     lut->SetTableValue(22, 1.0000, 0.3882, 0.2784, 1); // Tomato
     lut->SetTableValue(23, 0.9608, 0.8706, 0.7020, 1); // Wheat
     lut->SetTableValue(24, 0.9020, 0.9020, 0.9804, 1); // Lavender
     lut->SetTableValue(25, 1.0000, 0.4900, 0.2500, 1); // Flesh
     lut->SetTableValue(26, 0.5300, 0.1500, 0.3400, 1); // Raspberry
     lut->SetTableValue(27, 0.9804, 0.5020, 0.4471, 1); // Salmon
     lut->SetTableValue(28, 0.7400, 0.9900, 0.7900, 1); // Mint
     lut->SetTableValue(29, 0.2000, 0.6300, 0.7900, 1); // Peacock
     lut->SetTableValue(30     , 0     , 0     , 0, 1);  //Black
     lut->SetTableValue(31, 0.8900, 0.8100, 0.3400, 1); // Banana
     lut->SetTableValue(32, 1.0000, 0.3882, 0.2784, 1); // Tomato
     lut->SetTableValue(33, 0.9608, 0.8706, 0.7020, 1); // Wheat
     lut->SetTableValue(34, 0.9020, 0.9020, 0.9804, 1); // Lavender
     lut->SetTableValue(35, 1.0000, 0.4900, 0.2500, 1); // Flesh
     lut->SetTableValue(36, 0.5300, 0.1500, 0.3400, 1); // Raspberry
     lut->SetTableValue(37, 0.9804, 0.5020, 0.4471, 1); // Salmon
     lut->SetTableValue(38, 0.7400, 0.9900, 0.7900, 1); // Mint

//       // Setup two colors - one for each line
//       unsigned char red[3] = {255, 0, 0};
//       unsigned char green[3] = {0, 255, 0};

       // Generate the colors for each point based on the color map
       vtkSmartPointer<vtkUnsignedCharArray> colors = vtkSmartPointer<vtkUnsignedCharArray>::New();
       colors->SetNumberOfComponents(3);
       colors->SetName("Colors");


       for(int i = 1; i < tableSize; i++) {
           double rgb[3];
           lut->GetColor(static_cast<double>(i)/(tableSize-1), rgb);

           unsigned char ucrgb[3];
           for(unsigned int j = 0; j < 3; j++) {
               ucrgb[j] = static_cast<unsigned char>(255.0 * rgb[j]);
           }

           colors->InsertNextTuple3(ucrgb[0], ucrgb[1], ucrgb[2]);

           // Print out what we have.
           /*std::cout << "(";
           PrintColour<double[3]>(rgb);
           std::cout << ") (";
           PrintColour<unsigned char[3]>(ucrgb);
           std::cout << ")" << std::endl;*/
       }

    //polydata->GetPointData()->SetScalars(colors);
    linepolydata->GetCellData()->SetScalars(colors);

    vtkSmartPointer<vtkPolyDataMapper> pointlinemapper = vtkSmartPointer<vtkPolyDataMapper>::New();
#if VTK_MAJOR_VERSION <= 5
  pointlinemapper->SetInput(linepolydata);
#else
  pointlinemapper->SetInputData(linepolydata);
#endif



    vtkSmartPointer<vtkActor>  lineactor = vtkActor::New();
    lineactor->SetMapper(pointlinemapper);
    lineactor->GetProperty()->SetLineWidth(4);
    //lineactor->GetProperty()->SetColor(quadColor);// (0,1,0) is green, (1,1,1) is white.
    lineactor->RotateX(rX);
    lineactor->RotateY(rY);
    lineactor->RotateZ(rZ);
    //lineactor->GetProperty()->SetOpacity(0.2);





    renderer->AddActor(lineactor);
}

//! Make a lookup table from a set of named colors.
/*
 * See: http://www.vtk.org/doc/nightly/html/classvtkColorTransferFunction.html
 */
void visualization::MakeLUT(size_t const & tableSize, vtkLookupTable *lut) {
//  vtkSmartPointer<vtkNamedColors> nc = vtkSmartPointer<vtkNamedColors>::New();

//  lut->SetNumberOfTableValues(tableSize);
//  lut->Build();

////  // Fill in a few known colors, the rest will be generated if needed
////  lut->SetTableValue(0, nc->GetColor4d("Black").GetData());
////  lut->SetTableValue(1, nc->GetColor4d("Banana").GetData());
////  lut->SetTableValue(2, nc->GetColor4d("Tomato").GetData());
////  lut->SetTableValue(3, nc->GetColor4d("Wheat").GetData());
////  lut->SetTableValue(4, nc->GetColor4d("Lavender").GetData());
////  lut->SetTableValue(5, nc->GetColor4d("Flesh").GetData());
////  lut->SetTableValue(6, nc->GetColor4d("Raspberry").GetData());
////  lut->SetTableValue(7, nc->GetColor4d("Salmon").GetData());
////  lut->SetTableValue(8, nc->GetColor4d("Mint").GetData());
////  lut->SetTableValue(9, nc->GetColor4d("Peacock").GetData());

//  // Fill in a few known colors, the rest will be generated if needed
//    lut->SetTableValue(0     , 0     , 0     , 0, 1);  //Black
//    lut->SetTableValue(1, 0.8900, 0.8100, 0.3400, 1); // Banana
//    lut->SetTableValue(2, 1.0000, 0.3882, 0.2784, 1); // Tomato
//    lut->SetTableValue(3, 0.9608, 0.8706, 0.7020, 1); // Wheat
//    lut->SetTableValue(4, 0.9020, 0.9020, 0.9804, 1); // Lavender
//    lut->SetTableValue(5, 1.0000, 0.4900, 0.2500, 1); // Flesh
//    lut->SetTableValue(6, 0.5300, 0.1500, 0.3400, 1); // Raspberry
//    lut->SetTableValue(7, 0.9804, 0.5020, 0.4471, 1); // Salmon
//    lut->SetTableValue(8, 0.7400, 0.9900, 0.7900, 1); // Mint
//    lut->SetTableValue(9, 0.2000, 0.6300, 0.7900, 1); // Peacock
//    lut->SetTableValue(0     , 0     , 0     , 0, 1);  //Black
//    lut->SetTableValue(1, 0.8900, 0.8100, 0.3400, 1); // Banana
//    lut->SetTableValue(2, 1.0000, 0.3882, 0.2784, 1); // Tomato
//    lut->SetTableValue(3, 0.9608, 0.8706, 0.7020, 1); // Wheat
//    lut->SetTableValue(4, 0.9020, 0.9020, 0.9804, 1); // Lavender
//    lut->SetTableValue(5, 1.0000, 0.4900, 0.2500, 1); // Flesh
//    lut->SetTableValue(6, 0.5300, 0.1500, 0.3400, 1); // Raspberry
//    lut->SetTableValue(7, 0.9804, 0.5020, 0.4471, 1); // Salmon
//    lut->SetTableValue(8, 0.7400, 0.9900, 0.7900, 1); // Mint
//    lut->SetTableValue(9, 0.2000, 0.6300, 0.7900, 1); // Peacock
//    lut->SetTableValue(0     , 0     , 0     , 0, 1);  //Black
//    lut->SetTableValue(1, 0.8900, 0.8100, 0.3400, 1); // Banana
//    lut->SetTableValue(2, 1.0000, 0.3882, 0.2784, 1); // Tomato
//    lut->SetTableValue(3, 0.9608, 0.8706, 0.7020, 1); // Wheat
//    lut->SetTableValue(4, 0.9020, 0.9020, 0.9804, 1); // Lavender
//    lut->SetTableValue(5, 1.0000, 0.4900, 0.2500, 1); // Flesh
//    lut->SetTableValue(6, 0.5300, 0.1500, 0.3400, 1); // Raspberry
//    lut->SetTableValue(7, 0.9804, 0.5020, 0.4471, 1); // Salmon
//    lut->SetTableValue(8, 0.7400, 0.9900, 0.7900, 1); // Mint
//    lut->SetTableValue(9, 0.2000, 0.6300, 0.7900, 1); // Peacock
//    lut->SetTableValue(0     , 0     , 0     , 0, 1);  //Black
//    lut->SetTableValue(1, 0.8900, 0.8100, 0.3400, 1); // Banana
//    lut->SetTableValue(2, 1.0000, 0.3882, 0.2784, 1); // Tomato
//    lut->SetTableValue(3, 0.9608, 0.8706, 0.7020, 1); // Wheat
//    lut->SetTableValue(4, 0.9020, 0.9020, 0.9804, 1); // Lavender
//    lut->SetTableValue(5, 1.0000, 0.4900, 0.2500, 1); // Flesh
//    lut->SetTableValue(6, 0.5300, 0.1500, 0.3400, 1); // Raspberry
//    lut->SetTableValue(7, 0.9804, 0.5020, 0.4471, 1); // Salmon
//    lut->SetTableValue(8, 0.7400, 0.9900, 0.7900, 1); // Mint
}



/* draw the points of u & v coordinate, changing with the interpolation level.
 * bool moved: true if draw the original u & v in the input m3d file; false if draw the new values by adding deltaU and deltaV to primitive.
 * quadpoints_u[q]: store the u coordinate of all the sub-quad of quad[q].
 * quadpoints_v[q]: store the v coordinate of all the sub-quad of quad[q].
*/
void visualization::drawPointsOfTheUVCoordinate(VectorQuadPoint quadpoints_u, VectorQuadPoint quadpoints_v){

    int quadNums = quadpoints_u.size();
    int pointNums = quadpoints_u[0].size();
    // Create the topology of the point (a vertex)
    vtkSmartPointer<vtkCellArray> vertices = vtkSmartPointer<vtkCellArray>::New();
    // Create the geometry of a point (the coordinate)
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    for(int quad=0; quad<quadNums; quad++){
        for(int pointNum=0; pointNum<pointNums; pointNum++){
            //cout<<"-----------------------------quad is "<<quad<<"-----pointNum is: "<<pointNum<<endl;
            vtkIdType pid[1];
            pid[0] = points->InsertNextPoint(quadpoints_u[quad][pointNum], quadpoints_v[quad][pointNum], 0);
            //cout<<""<<quadpoints_u[quad][pointNum]<<"  "<<quadpoints_v[quad][pointNum]<<endl;
            vertices->InsertNextCell(1,pid);
        }
    }

    // Create a polydata object
    vtkSmartPointer<vtkPolyData> point = vtkSmartPointer<vtkPolyData>::New();

    // Set the points and vertices we created as the geometry and topology of the polydata
    point->SetPoints(points);
    point->SetVerts(vertices);

    // Visualize
    vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    #if VTK_MAJOR_VERSION <= 5
      mapper->SetInput(point);
    #else
      mapper->SetInputData(point);
    #endif

      vtkSmartPointer<vtkActor> pointactor = vtkSmartPointer<vtkActor>::New();
      pointactor->SetMapper(mapper);
      pointactor->GetProperty()->SetPointSize(5);
      pointactor->GetProperty()->SetColor(quadColor);// (0,1,0) is green, (1,1,1) is white.
      pointactor->RotateZ(100);
      pointactor->RotateX(-115);
      renderer->AddActor(pointactor);
}


//calculate the distance between two points.
double visualization::lengthofedges(Vector3D point[2]){
    double x0 = point[0].getX();
    double y0 = point[0].getY();
    double z0 = point[0].getZ();
    //cout<<"point[0] is:"<<x0<<"; "<<y0<<"; "<<z0<<endl;

    double x1 = point[1].getX();
    double y1 = point[1].getY();
    double z1 = point[1].getZ();
    //cout<<"point[1] is:"<<x1<<"; "<<y1<<"; "<<z1<<endl;

    double distance = sqrt(pow(x1-x0, 2) + pow(y1-y0, 2) + pow(z1-z0, 2));

    //cout<<"distance between two edges is: "<<distance<<endl;

    return distance;
}


/* Draw the horizonal and vertical line of each quads. Each side of these quads only drawn once, no overlap.
*/
void visualization::drawEdgeOfQuads_method3(vtkSmartPointer< vtkPoints > hubpos){

    vtkSmartPointer<vtkCellArray> cellarraypointsline = vtkSmartPointer<vtkCellArray>::New();
    vtkSmartPointer<vtkPolyData> polypointsline = vtkSmartPointer<vtkPolyData>::New();

    double p[3];

    Vector3D point[2];
    vector<double> horizonaledges;

    //For each quad, its left side is horizonal lines.
    for(int m=0; m<quadNum;m++){
        cout<<"H"<<m+1<< " subline postion is:"<<endl;
        double sumOfSubLines = 0;
        //left side. For each sub-line.
        for(int n=0; n<step;n++){
            vtkSmartPointer<vtkLine> medialsheetline = vtkSmartPointer<vtkLine>::New();
            hubpos->GetPoint(m*subQuadPointsNum+n,p);
            vtkIdType id0 = hubpos->InsertNextPoint(p[0],p[1],p[2]);
            point[0].setX(p[0]);
            point[0].setY(p[1]);
            point[0].setZ(p[2]);

            hubpos->GetPoint(m*subQuadPointsNum+n+1,p);
            vtkIdType id1 = hubpos->InsertNextPoint(p[0],p[1],p[2]);
            point[1].setX(p[0]);
            point[1].setY(p[1]);
            point[1].setZ(p[2]);

            medialsheetline->GetPointIds()->SetId(0, id0);
            medialsheetline->GetPointIds()->SetId(1, id1);

            cellarraypointsline->InsertNextCell(medialsheetline);

            sumOfSubLines += lengthofedges(point);
        }
        horizonaledges.push_back(sumOfSubLines);
    }

    //For the quads (colNum-1)*colNumIndex -1, its right side is also horizonal lines.
    for(int m=0; m<quadNum;m++){
        double sumOfSubLines = 0;
        if((m+1)%(colNum-1)==0){
            cout<<"H"<<quadNum+(m+1)/(colNum-1)<< " subline position is:"<<endl;
            for(int n=subQuadPointsNum-(step+1); n<subQuadPointsNum-1;n++){
                vtkSmartPointer<vtkLine> medialsheetline = vtkSmartPointer<vtkLine>::New();
                hubpos->GetPoint(m*subQuadPointsNum+n,p);
                vtkIdType id0 = hubpos->InsertNextPoint(p[0],p[1],p[2]);
                point[0].setX(p[0]);
                point[0].setY(p[1]);
                point[0].setZ(p[2]);

                hubpos->GetPoint(m*subQuadPointsNum+n+1,p);
                vtkIdType id1 = hubpos->InsertNextPoint(p[0],p[1],p[2]);
                point[1].setX(p[0]);
                point[1].setY(p[1]);
                point[1].setZ(p[2]);

                medialsheetline->GetPointIds()->SetId(0, id0);
                medialsheetline->GetPointIds()->SetId(1, id1);

                cellarraypointsline->InsertNextCell(medialsheetline);
                sumOfSubLines += lengthofedges(point);
            }
            horizonaledges.push_back(sumOfSubLines);
        }
    }

    //print the horizonal lines length;
    for(unsigned o =0; o<horizonaledges.size();o++){
        cout<<"------------------line length is: "<<horizonaledges[o]<< "  ";
    }
    cout<<endl;

    //Vertical lines.
    //For each quad, its top side is vertical lines.
    for(int m=0; m<quadNum;m++){
        for(int k = 0; k<step;k++){
            vtkSmartPointer<vtkLine> medialsheetline = vtkSmartPointer<vtkLine>::New();
            hubpos->GetPoint(m*subQuadPointsNum+0+k*(step+1),p);
            vtkIdType id0 = hubpos->InsertNextPoint(p[0],p[1],p[2]);

            hubpos->GetPoint(m*subQuadPointsNum+(k+1)*(step+1),p);
            vtkIdType id1 = hubpos->InsertNextPoint(p[0],p[1],p[2]);

            medialsheetline->GetPointIds()->SetId(0, id0);
            medialsheetline->GetPointIds()->SetId(1, id1);

            cellarraypointsline->InsertNextCell(medialsheetline);
        }
    }

    //For the last colNum-1 quads, its bottom side is also vertical lines.
    for(int m=quadNum-(colNum-1); m<quadNum;m++){
        for(int k = 0; k<step;k++){
            vtkSmartPointer<vtkLine> medialsheetline = vtkSmartPointer<vtkLine>::New();
            hubpos->GetPoint(m*subQuadPointsNum+step + k*(step+1),p);
            vtkIdType id0 = hubpos->InsertNextPoint(p[0],p[1],p[2]);

            hubpos->GetPoint(m*subQuadPointsNum+step + (k+1)*(step+1),p);
            vtkIdType id1 = hubpos->InsertNextPoint(p[0],p[1],p[2]);

            medialsheetline->GetPointIds()->SetId(0, id0);
            medialsheetline->GetPointIds()->SetId(1, id1);

            cellarraypointsline->InsertNextCell(medialsheetline);
        }
    }

    polypointsline->SetPoints(hubpos);
    polypointsline->SetLines(cellarraypointsline);

    vtkSmartPointer<vtkPolyDataMapper> pointlinemapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    pointlinemapper->SetInput(polypointsline);
    vtkSmartPointer<vtkActor>  lineactor = vtkActor::New();
    lineactor->SetMapper(pointlinemapper);
    lineactor->GetProperty()->SetLineWidth(1);
    lineactor->GetProperty()->SetColor(quadColor);// (0,1,0) is green, (1,1,1) is white.
    //lineactor->RotateY(95);
    //lineactor->RotateX(30);
    //lineactor->RotateY(-35);(35);
    renderer->AddActor(lineactor);
}



void visualization::setRenderer(vtkSmartPointer<vtkRenderer> renderer){
    this->renderer = renderer;
}



/* Return all the sub quads's u & v coordinate into vector quadpoints_u & quadpoints_v.
 * The u, v is the moved value with deltaU and deltaV.
 * quadpoints_u[q]: store the u coordinate of all the sub-quad of quad[q].
 * quadpoints_v[q]: store the v coordinate of all the sub-quad of quad[q].
*/
void visualization::getSubQuadsUVCoordinateUsingDeltaUV(){

    int quadIndex = 0;

    M3DQuadPrimitive* prim0;
    vector<double> subquadpoint_v, subquadpoint_u;

    for(unsigned i = 0; i < rowNum -1; i++){ //the row number of the quads. its 3.
        for(unsigned j = 0; j < colNum -1; j++){//coloums, its 13.
                quadpoints_u.push_back(VectorDoublePoints());
                quadpoints_v.push_back(VectorDoublePoints());

                //deltau[0]: left-top point; deltau[1]: left-bottom point; deltau[2]: right-bottom; deltau[3]: right-top
                double deltau[4], deltav[4];

                switch(side){
                case 0: //top side
                    //get the primitive's delta u, v
                    //primitive[i][j]
                    prim0= dynamic_cast<M3DQuadPrimitive*>(curQuadFig->getPrimitivePtr(i, j));
                    deltau[0] = prim0->getDeltaU0();
                    deltav[0] = prim0->getDeltaV0();

                    //primitive[i+1][j]
                    prim0= dynamic_cast<M3DQuadPrimitive*>(curQuadFig->getPrimitivePtr(i+1, j));
                    deltau[1] = prim0->getDeltaU0();
                    deltav[1] = prim0->getDeltaV0();

                    //primitive[i+1][j+1]
                    prim0= dynamic_cast<M3DQuadPrimitive*>(curQuadFig->getPrimitivePtr(i+1, j+1));
                    deltau[2] = prim0->getDeltaU0();
                    deltav[2] = prim0->getDeltaV0();

                    //primitive[i][j+1]
                    prim0= dynamic_cast<M3DQuadPrimitive*>(curQuadFig->getPrimitivePtr(i, j+1));
                    deltau[3] = prim0->getDeltaU0();
                    deltav[3] = prim0->getDeltaV0();
                    //delete prim0;
                    break;
                case 1:  //down side
                    //get the primitive's delta u, v
                    //primitive[i][j]
                    prim0= dynamic_cast<M3DQuadPrimitive*>(curQuadFig->getPrimitivePtr(i, j));
                    deltau[0] = prim0->getDeltaU1();
                    deltav[0] = prim0->getDeltaV1();

                    //primitive[i+1][j]
                    prim0= dynamic_cast<M3DQuadPrimitive*>(curQuadFig->getPrimitivePtr(i+1, j));
                    deltau[1] = prim0->getDeltaU1();
                    deltav[1] = prim0->getDeltaV1();

                    //primitive[i+1][j+1]
                    prim0= dynamic_cast<M3DQuadPrimitive*>(curQuadFig->getPrimitivePtr(i+1, j+1));
                    deltau[2] = prim0->getDeltaU1();
                    deltav[2] = prim0->getDeltaV1();

                    //primitive[i][j+1]
                    prim0= dynamic_cast<M3DQuadPrimitive*>(curQuadFig->getPrimitivePtr(i, j+1));
                    deltau[3] = prim0->getDeltaU1();
                    deltav[3] = prim0->getDeltaV1();
                    //delete prim0;
                    break;
                default:
                    break;
                }

                //four points of a quad with delta u, v. counter clockwise.
                double quadpointsu[4], quadpointsv[4];
                //left-top point
                quadpointsu[0] = i + deltau[0];
                quadpointsv[0] = j + deltav[0];
                //left-bottom point
                quadpointsu[1] = i + 1 + deltau[1];
                quadpointsv[1] = j + deltav[1];
                //right-bottom point
                quadpointsu[2] = i + 1 + deltau[2];
                quadpointsv[2] = j + 1 + deltav[2];
                //right-top point
                quadpointsu[3] = i + deltau[3];
                quadpointsv[3] = j + 1 + deltav[3];

                //given four points of a quad (p0, p1, p2, p3), split each side of the quad into 2^interpolationLevel sub-line.
                //first, get the subquads's u coordinate, in column first order.
                subquadpoint_u = tools.splitQuad(quadpointsu[0],quadpointsu[1],quadpointsu[2],quadpointsu[3],step);

                //second, get the 25 subquads's v coordinate, in column first order.
                subquadpoint_v = tools.splitQuad(quadpointsv[0],quadpointsv[1],quadpointsv[2],quadpointsv[3],step);

                for(unsigned m =0; m<subquadpoint_u.size();m++){
                    quadpoints_u[quadIndex].push_back(subquadpoint_u[m]);
                    quadpoints_v[quadIndex].push_back(subquadpoint_v[m]);
                    //cout<<"MSG from getSubQuadsUVCoordinateUsingDeltaUV: quadpoints_u["<<quadIndex<<"]---------m is: "<<m<<"---subquadpoint_u["<<m<<"] is: "<<subquadpoint_u[m]<<endl;
                    //cout<<"MSG from getSubQuadsUVCoordinateUsingDeltaUV: quadpoints_v["<<quadIndex<<"]---------m is: "<<m<<"---subquadpoint_v["<<m<<"] is: "<<subquadpoint_v[m]<<endl;
                }

                quadIndex++;

                subquadpoint_v.clear();
                subquadpoint_u.clear();
        }
    }
}




/* Return all the sub quads's u & v coordinate into vector quadpoints_u & quadpoints_v.
 * Without considering the deltaU and deltaV in the m3d file(because the primitive has already been moved
 * and updated using them during write to m3d in privious steps).
 * quadpoints_u[q]: store the u coordinate of all the sub-quad of quad[q].
 * quadpoints_v[q]: store the v coordinate of all the sub-quad of quad[q].
*/
void visualization::getSubQuadsUVCoordinateNotUsingDeltaUV(){

    int quadIndex = 0;

    vector<double> subquadpoint_v, subquadpoint_u;

    for(unsigned i = 0; i < rowNum -1; i++){ //the row number of the quads. its 3.
        for(unsigned j = 0; j < colNum -1; j++){//coloums, its 13.
                quadpoints_u.push_back(VectorDoublePoints());
                quadpoints_v.push_back(VectorDoublePoints());

                //four points of a quad with delta u, v. counter clockwise.
                double quadpointsu[4], quadpointsv[4];
                //left-top point
                quadpointsu[0] = i;
                quadpointsv[0] = j;
                //left-bottom point
                quadpointsu[1] = i + 1;
                quadpointsv[1] = j;
                //right-bottom point
                quadpointsu[2] = i + 1;
                quadpointsv[2] = j + 1;
                //right-top point
                quadpointsu[3] = i;
                quadpointsv[3] = j + 1;

                //given four points of a quad (p0, p1, p2, p3), split each side of the quad into 2^interpolationLevel sub-line.
                //first, get the subquads's u coordinate, in column first order.
                subquadpoint_u = tools.splitQuad(quadpointsu[0],quadpointsu[1],quadpointsu[2],quadpointsu[3], step);

                //second, get the 25 subquads's v coordinate, in column first order.
                subquadpoint_v = tools.splitQuad(quadpointsv[0],quadpointsv[1],quadpointsv[2],quadpointsv[3], step);

                for(unsigned m =0; m<subquadpoint_u.size();m++){
                    quadpoints_u[quadIndex].push_back(subquadpoint_u[m]);
                    quadpoints_v[quadIndex].push_back(subquadpoint_v[m]);

                    //cout<<"MSG from getSubQuadsUVCoordinateNotUsingDeltaUV: quadpoints_u["<<quadIndex<<"]---------m is: "<<m<<"---subquadpoint_u["<<m<<"] is: "<<subquadpoint_u[m]<<endl;
                    //cout<<"MSG from getSubQuadsUVCoordinateNotUsingDeltaUV: quadpoints_v["<<quadIndex<<"]---------m is: "<<m<<"---subquadpoint_v["<<m<<"] is: "<<subquadpoint_v[m]<<endl;
                }

                quadIndex++;

                subquadpoint_v.clear();
                subquadpoint_u.clear();
        }
    }
}






/* Get all the sub-quad position by given u v coordinates.
 * quadtype: boundary quads(0), skeletal quads(1).
*/
vtkSmartPointer< vtkPoints > visualization::getSubQuadsPosition(int quadtype){

    vtkSmartPointer< vtkPoints > hubpos = vtkSmartPointer< vtkPoints >::New();

    M3DQuadInterpolater *tpm = new M3DQuadInterpolater(curQuadFig);

    //store the all the position of points for a quad, each point is a 3D vector, which store the x, y, z coordinate of this point.
    Vector3D point;

    for(int i =0; i<quadNum; i++){
        //cout<<"Currently drawing quad: "<<i<<endl;
        for(int j =0; j<subQuadPointsNum;j++){
            //quadtype: Two kinds of quads, boundary quads(0), skeletal quads(1)
            switch(quadtype){
                case 0:
                    //boundary quads(0)
                    //cout<<"quad: "<<i<<", point: "<<j<<" u is: "<<quadpoints_u[i][j]<< " v is: "<<quadpoints_v[i][j]<<endl;
                    point = tpm->interpolateQuadSpoke(curQuadFig,quadpoints_u[i][j],quadpoints_v[i][j],side)->getB();
                    //cout<<"point position: x is: "<<point.getX()<<" y is: "<<point.getY()<<" z is: "<<point.getZ()<<endl;

                    hubpos->InsertNextPoint(point.getX(),point.getY(),point.getZ());
                    break;
                case 1:
                    //skeletal quads(1)
                    //cout<<"quad: "<<i<<", point: "<<j<<" u is: "<<quadpoints_u[i][j]<< " v is: "<<quadpoints_v[i][j]<<endl;
                    point = tpm->interpolateQuadSpoke(curQuadFig,quadpoints_u[i][j],quadpoints_v[i][j],side)->getX();
                    //cout<<"point position: x is: "<<point.getX()<<" y is: "<<point.getY()<<" z is: "<<point.getZ()<<endl;
                    hubpos->InsertNextPoint(point.getX(),point.getY(),point.getZ());
                    break;
                default:
                    break;
            }
        }
    }

    delete tpm;

    return hubpos;
}




void visualization::setColor(double *quadColor){
    this->quadColor = quadColor;
}



/* Constructor for up and down
*/
visualization::visualization(M3DQuadFigure* quadfig, int interpolationLevel, int side, bool moved, double *quadColor,
                             vtkSmartPointer<vtkRenderer> renderer, double rX, double rY, double rZ, double opacity){
    this->curQuadFig = quadfig;
    this->rowNum = quadfig->getRowCount();
    this->colNum = quadfig->getColumnCount();
    this->quadNum = (rowNum-1)*(colNum-1);
    this->interpolationLevel = interpolationLevel;
    this->side = side;
    this->step = pow((double)2, (double)interpolationLevel);
    this->subQuadPointsNum = (step+1)*(step+1);

    //initialize the quadpoints_u, quadpoints_v vector of all the subquads. Each srep has diff quadpoints_u, quadpoints_v vector.
    if(moved){//use original primitive's u p r in the input m3d file.
        getSubQuadsUVCoordinateNotUsingDeltaUV();
    }
    else   {//use new primitive's u p r by add deltaU & deltaV.
        getSubQuadsUVCoordinateUsingDeltaUV();
    }

    this->quadColor = quadColor;
    this->renderer = renderer;

    this->rX = rX;
    this->rY = rY;
    this->rZ = rZ;
    this->opacity = opacity;
}



/* Constructor for crest
*/
visualization::visualization(M3DQuadFigure* quadfig, int interpolationLevel, int side, double *quadColor,
                             vtkSmartPointer<vtkRenderer> renderer, double rX, double rY, double rZ, double opacity){
    this->curQuadFig = quadfig;
    this->rowNum = quadfig->getRowCount();
    this->colNum = quadfig->getColumnCount();
    this->quadNum = (rowNum-1)*(colNum-1);//quadNum means the quad numbers on skeletal sheet.
    this->interpolationLevel = interpolationLevel;
    this->side = side;
    this->step = pow((double)2, (double)interpolationLevel);
    this->subQuadPointsNum = (step+1)*(step+1);

    this->rX = rX;
    this->rY = rY;
    this->rZ = rZ;
    this->opacity = opacity;

    M3DQuadPrimitive* prim;
    M3DQuadEndPrimitive* endPrim;
    Vector3D up_b_mid_point, down_b_mid_point, up_b_crest_point, down_b_crest_point, crest_b_point;

    // Loop each primitive's up and down spokes.
    for(unsigned u =0; u<this->rowNum; u++){
        for(unsigned v =0; v<this->colNum; v++){
            prim = dynamic_cast<M3DQuadPrimitive*>(this->curQuadFig->getPrimitivePtr(u,v));
            // Up boundary points
            up_b_mid_point = prim->getX() + prim->getR0()*prim->getU0();
            this->up_b_points.push_back(up_b_mid_point);

            // Down boundary points
            down_b_mid_point = prim->getX() + prim->getR1()*prim->getU1();
            this->down_b_points.push_back(down_b_mid_point);
        }
    }

    // Get the crest curve along the crest in counter clockwise
    for(unsigned v =0; v<this->colNum; v++){
        // Up boundary crest points. (End atom's up spoke's tip on the boundary). 2 quads use IM1, 2 quads use IM2.
        prim = dynamic_cast<M3DQuadPrimitive*>(this->curQuadFig->getPrimitivePtr(0,v));
        up_b_crest_point = prim->getX() + prim->getR0()*prim->getU0();
        this->up_b_crest_points.push_back(up_b_crest_point);

        // Down boundary crest points. (End atom's down spoke's tip on the boundary). 2 quads use IM1, 2 quads use IM2.
        down_b_crest_point = prim->getX() + prim->getR1()*prim->getU1();
        this->down_b_crest_points.push_back(down_b_crest_point);

        // Crest boundary points. (crest spokes' tip on the boundary.)
        endPrim = dynamic_cast<M3DQuadEndPrimitive*>(this->curQuadFig->getPrimitivePtr(0, v));
        crest_b_point = endPrim->getX() + endPrim->getREnd()*endPrim->getUEnd();
        this->crest_b_points.push_back(crest_b_point);
    }
    for(unsigned u =1; u<this->rowNum; u++){
        // Up boundary crest points. (End atom's up spoke's tip on the boundary). 2 quads use IM1, 2 quads use IM2.
        prim = dynamic_cast<M3DQuadPrimitive*>(this->curQuadFig->getPrimitivePtr(u, this->colNum - 1));
        up_b_crest_point = prim->getX() + prim->getR0()*prim->getU0();
        this->up_b_crest_points.push_back(up_b_crest_point);

        // Down boundary crest points. (End atom's down spoke's tip on the boundary). 2 quads use IM1, 2 quads use IM2.
        down_b_crest_point = prim->getX() + prim->getR1()*prim->getU1();
        this->down_b_crest_points.push_back(down_b_crest_point);

        // Crest boundary points. (crest spokes' tip on the boundary.)
        endPrim = dynamic_cast<M3DQuadEndPrimitive*>(this->curQuadFig->getPrimitivePtr(u, this->colNum - 1));
        crest_b_point = endPrim->getX() + endPrim->getREnd()*endPrim->getUEnd();
        this->crest_b_points.push_back(crest_b_point);
    }
    for(unsigned v =this->colNum-2; v>0; v--){
        // Up boundary crest points. (End atom's up spoke's tip on the boundary). 2 quads use IM1, 2 quads use IM2.
        prim = dynamic_cast<M3DQuadPrimitive*>(this->curQuadFig->getPrimitivePtr(this->rowNum-1, v));
        up_b_crest_point = prim->getX() + prim->getR0()*prim->getU0();
        this->up_b_crest_points.push_back(up_b_crest_point);

        // Down boundary crest points. (End atom's down spoke's tip on the boundary). 2 quads use IM1, 2 quads use IM2.
        down_b_crest_point = prim->getX() + prim->getR1()*prim->getU1();
        this->down_b_crest_points.push_back(down_b_crest_point);

        // Crest boundary points. (crest spokes' tip on the boundary.)
        endPrim = dynamic_cast<M3DQuadEndPrimitive*>(this->curQuadFig->getPrimitivePtr(this->rowNum-1, v));
        crest_b_point = endPrim->getX() + endPrim->getREnd()*endPrim->getUEnd();
        this->crest_b_points.push_back(crest_b_point);
    }
    for(unsigned u =this->rowNum-1; u>0; u--){
        // Up boundary crest points. (End atom's up spoke's tip on the boundary). 2 quads use IM1, 2 quads use IM2.
        prim = dynamic_cast<M3DQuadPrimitive*>(this->curQuadFig->getPrimitivePtr(u, 0));
        up_b_crest_point = prim->getX() + prim->getR0()*prim->getU0();
        this->up_b_crest_points.push_back(up_b_crest_point);

        // Down boundary crest points. (End atom's down spoke's tip on the boundary). 2 quads use IM1, 2 quads use IM2.
        down_b_crest_point = prim->getX() + prim->getR1()*prim->getU1();
        this->down_b_crest_points.push_back(down_b_crest_point);

        // Crest boundary points. (crest spokes' tip on the boundary.)
        endPrim = dynamic_cast<M3DQuadEndPrimitive*>(this->curQuadFig->getPrimitivePtr(u, 0));
        crest_b_point = endPrim->getX() + endPrim->getREnd()*endPrim->getUEnd();
        this->crest_b_points.push_back(crest_b_point);
    }

    this->quadColor = quadColor;
    this->renderer = renderer;
}





/* Input: a set of points counter clockwise along crest.
 * Draw curve between these discreat points.
 * These points in a sequence along the crest.
*/
void visualization::drawCurveAlongDiscreatPoints(vector<Vector3D> points, double *quadColor){
    vtkSmartPointer<vtkPoints> pointscrestspokes  = vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkCellArray> cellarraylines = vtkSmartPointer<vtkCellArray>::New();
    vtkSmartPointer<vtkPolyData> polycrestspokes = vtkSmartPointer<vtkPolyData>::New();

    unsigned n = points.size();
    //connect this point and its next point.
    for(unsigned i = 0; i < n-1; i++){
        vtkSmartPointer<vtkLine> line = vtkSmartPointer<vtkLine>::New();
        line->GetPointIds()->SetId(0, pointscrestspokes->InsertNextPoint(points[i].getX(),points[i].getY(),points[i].getZ()));
        line->GetPointIds()->SetId(1, pointscrestspokes->InsertNextPoint(points[i+1].getX(),points[i+1].getY(),points[i+1].getZ()));
        cellarraylines->InsertNextCell(line);
    }

    //connect the first point and the last point, make a closed curve.
    vtkSmartPointer<vtkLine> line = vtkSmartPointer<vtkLine>::New();
    line->GetPointIds()->SetId(0, pointscrestspokes->InsertNextPoint(points[0].getX(),points[0].getY(),points[0].getZ()));
    line->GetPointIds()->SetId(1, pointscrestspokes->InsertNextPoint(points[n-1].getX(),points[n-1].getY(),points[n-1].getZ()));
    cellarraylines->InsertNextCell(line);


    polycrestspokes->SetPoints(pointscrestspokes);
    polycrestspokes->SetLines(cellarraylines);

    vtkSmartPointer<vtkPolyDataMapper> crestspokesmapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    crestspokesmapper->SetInput(polycrestspokes);
    vtkSmartPointer<vtkActor>  spokesactor = vtkActor::New();
    spokesactor->SetMapper(crestspokesmapper);
    spokesactor->GetProperty()->SetLineWidth(1);
    spokesactor->GetProperty()->SetColor(quadColor);
    spokesactor->RotateZ(90);
    spokesactor->RotateX(-115);
    renderer->AddActor(spokesactor);
}


/* up: if ture, draw up boundary crest curve.
 * down: if ture, draw down boundary crest curve.
 * medial: if ture, draw medial sheet crest curve.(Composite from the skeletal's crest spoke tips.)
*/
void visualization::drawCrest(bool up, bool down, bool medial){

    double quadColor[3] = {1,0.5,0};
    if(medial)
        drawCurveAlongDiscreatPoints(this->crest_b_points, quadColor);//draw the boundary crest curve.
    if(up){
        quadColor[0] = 1;
        quadColor[1] = 0;
        quadColor[2] = 0;
        drawCurveAlongDiscreatPoints(this->up_b_crest_points, quadColor);//draw the up boundary crest curve.
    }

    if(down){
        quadColor[0] = 0;
        quadColor[1] = 1;
        quadColor[2] = 0;
        drawCurveAlongDiscreatPoints(this->down_b_crest_points, quadColor);//draw the down boundary crest curve.
    }
}





/* Connect two correspondence points. (e.g. one on the skeletal sheet, one on the boundary)
 * points_A, points_a has same size. Store a same sequence correspondence points, this function will connect A and a.
*/
void visualization::connectCrestCorrespondenceLines(vector<Vector3D> points_A, vector<Vector3D> points_a){
    vtkSmartPointer<vtkPoints> pointscrestspokes  = vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkCellArray> cellarraylines = vtkSmartPointer<vtkCellArray>::New();
    vtkSmartPointer<vtkPolyData> polycrestspokes = vtkSmartPointer<vtkPolyData>::New();

    //connect two correspondence points.
    for(unsigned i = 0; i < points_A.size(); i++){
        vtkSmartPointer<vtkLine> line = vtkSmartPointer<vtkLine>::New();
        line->GetPointIds()->SetId(0, pointscrestspokes->InsertNextPoint(points_A[i].getX(),points_A[i].getY(),points_A[i].getZ()));
        line->GetPointIds()->SetId(1, pointscrestspokes->InsertNextPoint(points_a[i].getX(),points_a[i].getY(),points_a[i].getZ()));
        cellarraylines->InsertNextCell(line);
    }

    polycrestspokes->SetPoints(pointscrestspokes);
    polycrestspokes->SetLines(cellarraylines);

    vtkSmartPointer<vtkPolyDataMapper> crestspokesmapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    crestspokesmapper->SetInput(polycrestspokes);
    vtkSmartPointer<vtkActor>  spokesactor = vtkActor::New();
    spokesactor->SetMapper(crestspokesmapper);
    spokesactor->GetProperty()->SetLineWidth(1);
    spokesactor->GetProperty()->SetColor(quadColor);
    spokesactor->RotateZ(90);
    spokesactor->RotateX(-115);
    renderer->AddActor(spokesactor);
}



/* draw the line between two correspondence points set.
 * crest_up: if true, will draw the lines of correspondence points between up boundary and medial sheet crest.
 * crest_down: if true, draw the lines of correspondence points between down boundary and medial sheet crest.
*/
void visualization::drawCrestCorrespondenceLines(bool crest_up, bool crest_down){
    bool c_up = false;
    bool c_down = false;

    c_up = crest_up;
    c_down = crest_down;

    // Draw lines between up boundary and crest
    if(c_up){
        this->quadColor[0]=0.8; this->quadColor[1]=0.5; this->quadColor[2]=0.5;
        connectCrestCorrespondenceLines(this->crest_b_points, this->up_b_crest_points);
    }
    if(c_down){
        this->quadColor[0]=0; this->quadColor[1]=0.5; this->quadColor[2]=0.5;
        connectCrestCorrespondenceLines(this->crest_b_points, this->down_b_crest_points);
    }
}



/* Connect the lines on the boundary, form quads on the top or down boundary.
 * Used to show the top boundary quads line or bottom boundary quad lines.
 * points: top or down boundary points set.
 * points contains the subquads points after interpolate.
*/
void visualization::connectInteriorToCrestLines(vector<Vector3D> points){
    vtkSmartPointer<vtkPoints> pointscrestspokes  = vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkCellArray> cellarraylines = vtkSmartPointer<vtkCellArray>::New();
    vtkSmartPointer<vtkPolyData> polycrestspokes = vtkSmartPointer<vtkPolyData>::New();

    // In the v direction(column direction). Each interi
    // If each side of the quad interpolate step points, the new row num will be:
    int rowNum = (this->rowNum-1)*(this->step-1) + this->rowNum;
    int colNum = (this->colNum-1)*(this->step-1) + this->colNum;

    // Split points into interiorRowNum pieces, each piece is an interior row points.
    // Connect line between row points, horizonal.
    for(unsigned m=0; m<rowNum; m++){
        for(unsigned i =0; i<colNum-1; i++){
            vtkSmartPointer<vtkLine> line = vtkSmartPointer<vtkLine>::New();
            line->GetPointIds()->SetId(0, pointscrestspokes->InsertNextPoint(points[m*colNum+i].getX(),points[m*colNum+i].getY(),points[m*colNum+i].getZ()));
            line->GetPointIds()->SetId(1, pointscrestspokes->InsertNextPoint(points[m*colNum+i+1].getX(),points[m*colNum+i+1].getY(),points[m*colNum+i+1].getZ()));
            cellarraylines->InsertNextCell(line);
        }
    }
    // Connect lines between vertical points.
    for(unsigned m =0; m<rowNum-1; m++){
        for(unsigned n=0; n<colNum;n++){
            vtkSmartPointer<vtkLine> line = vtkSmartPointer<vtkLine>::New();
            line->GetPointIds()->SetId(0, pointscrestspokes->InsertNextPoint(points[m*colNum+n].getX(),points[m*colNum+n].getY(),points[m*colNum+n].getZ()));
            line->GetPointIds()->SetId(1, pointscrestspokes->InsertNextPoint(points[m*colNum+n+colNum].getX(),points[m*colNum+n+colNum].getY(),points[m*colNum+n+colNum].getZ()));
            cellarraylines->InsertNextCell(line);
        }
    }

    polycrestspokes->SetPoints(pointscrestspokes);
    polycrestspokes->SetLines(cellarraylines);

    vtkSmartPointer<vtkPolyDataMapper> crestspokesmapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    crestspokesmapper->SetInput(polycrestspokes);
    vtkSmartPointer<vtkActor>  spokesactor = vtkActor::New();
    spokesactor->SetMapper(crestspokesmapper);
    spokesactor->GetProperty()->SetLineWidth(1);
    spokesactor->GetProperty()->SetColor(quadColor);
    spokesactor->RotateZ(90);
    spokesactor->RotateX(-115);
    renderer->AddActor(spokesactor);
}



/* Draw the quads of the top boundary or down boundary.
 * Compared to drawUpOrDownBoundaryQuadLine_method2, this method can only handle interpolate level 0, that is the original grid.
 * Because the up_b_points and down_b_points only got by primitive...not by interpolate.
*/
void visualization::drawUpOrDownBoundaryQuadLine_method1(bool up, bool down){
    if(up){
        connectInteriorToCrestLines(this->up_b_points);
    }
    if(down){
        connectInteriorToCrestLines(this->down_b_points);
    }

}




/* Get all the sub-quad position by interpolate to u v coordinate. *
 * quadtype: boundary quads(0), skeletal quads(1).
*/
void visualization::getAllSubQuadsVertexesOnBoundary(){

    int quadIndex = 0;

    vector<double> subquadpoint_v, subquadpoint_u;
    VectorQuadPoint quadpoints_u;   //store the u coordinate of all the sub-quad of quad[q].
    VectorQuadPoint quadpoints_v;   //store the v coordinate of all the sub-quad of quad[q].

    // Interpolate to the u v coordinate.
    for(unsigned i = 0; i < this->rowNum -1; i++){ //the row number of the quads. its 3.
        for(unsigned j = 0; j < this->colNum -1; j++){//coloums, its 13.
                quadpoints_u.push_back(VectorDoublePoints());
                quadpoints_v.push_back(VectorDoublePoints());

                //four points of a quad with delta u, v. counter clockwise.
                double quadpointsu[4], quadpointsv[4];
                //left-top point
                quadpointsu[0] = i;
                quadpointsv[0] = j;
                //left-bottom point
                quadpointsu[1] = i + 1;
                quadpointsv[1] = j;
                //right-bottom point
                quadpointsu[2] = i + 1;
                quadpointsv[2] = j + 1;
                //right-top point
                quadpointsu[3] = i;
                quadpointsv[3] = j + 1;

                //given four points of a quad (p0, p1, p2, p3), split each side of the quad into 2^interpolationLevel sub-line.
                //first, get the subquads's u coordinate, in column first order.
                subquadpoint_u = tools.splitQuad(quadpointsu[0],quadpointsu[1],quadpointsu[2],quadpointsu[3], this->step);

                //second, get the 25 subquads's v coordinate, in column first order.
                subquadpoint_v = tools.splitQuad(quadpointsv[0],quadpointsv[1],quadpointsv[2],quadpointsv[3], this->step);

                for(unsigned m =0; m<subquadpoint_u.size();m++){
                    quadpoints_u[quadIndex].push_back(subquadpoint_u[m]);
                    quadpoints_v[quadIndex].push_back(subquadpoint_v[m]);

                    //cout<<"MSG from getSubQuadsUVCoordinateNotUsingDeltaUV: quadpoints_u["<<quadIndex<<"]---------m is: "<<m<<"---subquadpoint_u["<<m<<"] is: "<<subquadpoint_u[m]<<endl;
                    //cout<<"MSG from getSubQuadsUVCoordinateNotUsingDeltaUV: quadpoints_v["<<quadIndex<<"]---------m is: "<<m<<"---subquadpoint_v["<<m<<"] is: "<<subquadpoint_v[m]<<endl;
                }

                quadIndex++;

                subquadpoint_v.clear();
                subquadpoint_u.clear();
        }
    }

    M3DQuadInterpolater *tpm = new M3DQuadInterpolater(this->curQuadFig);

    //store all the positions of points for a quad, each point is a 3D vector, which store the x, y, z coordinate of this point.
    Vector3D point;

    for(int i =0; i< this->quadNum; i++){
        for(int j =0; j< this->subQuadPointsNum; j++){
            //cout<<"quad: "<<i<<", point: "<<j<<" u is: "<<quadpoints_u[i][j]<< " v is: "<<quadpoints_v[i][j]<<endl;
            point = tpm->interpolateQuadSpoke(this->curQuadFig,quadpoints_u[i][j],quadpoints_v[i][j],0)->getB();
            //cout<<"point position: x is: "<<point.getX()<<" y is: "<<point.getY()<<" z is: "<<point.getZ()<<endl;
            this->up_BP.push_back(point);

            // Down boundary
            point = tpm->interpolateQuadSpoke(this->curQuadFig,quadpoints_u[i][j],quadpoints_v[i][j],1)->getB();
            this->down_BP.push_back(point);
        }
    }

    delete tpm;
}




/* Draw the quad frame, by index.
 * The only difference of this function with drawFrameOfQuadByIndex() is: input is a vectore of Vector3D.
 * The points in the vector are same order as drawFrameOfQuadByIndex.
*/
void visualization::drawFrameOfQuadByIndex_m2(vector<Vector3D> points, int quadIndex){

    vtkSmartPointer<vtkCellArray> cellarraypointsline = vtkSmartPointer<vtkCellArray>::New();
    vtkSmartPointer<vtkPolyData> polypointsline = vtkSmartPointer<vtkPolyData>::New();
    vtkSmartPointer< vtkPoints > hubPosNewSequence = vtkSmartPointer< vtkPoints >::New();
    Vector3D p;

    //draw horizonal lines.
    for(int n=0; n<step+1;n++){
        int pointIndex =0;
        for(int k = 0; k<step&&pointIndex<subQuadPointsNum-1;k++){
            pointIndex = n + k*(step+1);
            vtkSmartPointer<vtkLine> medialsheetline = vtkSmartPointer<vtkLine>::New();

            p = points[quadIndex*subQuadPointsNum+pointIndex];
            vtkIdType id0 = hubPosNewSequence->InsertNextPoint(p.getX(), p.getY(), p.getZ());

            p = points[quadIndex*subQuadPointsNum+pointIndex+step+1];
            vtkIdType id1 = hubPosNewSequence->InsertNextPoint(p.getX(), p.getY(), p.getZ());

            medialsheetline->GetPointIds()->SetId(0, id0);
            medialsheetline->GetPointIds()->SetId(1, id1);

            cellarraypointsline->InsertNextCell(medialsheetline);
        }
    }

    //draw vertical lines.
    for(int n=0; n<step+1;n++){
        for(int k = 0; k<step;k++){
            int currentpoint = k + n*(step+1);
            vtkSmartPointer<vtkLine> medialsheetline = vtkSmartPointer<vtkLine>::New();

            p = points[quadIndex*subQuadPointsNum+currentpoint];
            vtkIdType id0 = hubPosNewSequence->InsertNextPoint(p.getX(), p.getY(), p.getZ());

            p = points[quadIndex*subQuadPointsNum+currentpoint+1];
            vtkIdType id1 = hubPosNewSequence->InsertNextPoint(p.getX(), p.getY(), p.getZ());

            medialsheetline->GetPointIds()->SetId(0, id0);
            medialsheetline->GetPointIds()->SetId(1, id1);

            cellarraypointsline->InsertNextCell(medialsheetline);
        }
    }

    polypointsline->SetPoints(hubPosNewSequence);
    polypointsline->SetLines(cellarraypointsline);

    vtkSmartPointer<vtkPolyDataMapper> pointlinemapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    pointlinemapper->SetInput(polypointsline);
    vtkSmartPointer<vtkActor>  lineactor = vtkActor::New();
    lineactor->SetMapper(pointlinemapper);
    lineactor->GetProperty()->SetLineWidth(1);
    lineactor->GetProperty()->SetColor(quadColor);// (0,1,0) is green, (1,1,1) is white.
    lineactor->RotateY(95);
    lineactor->RotateX(30);
    //lineactor->RotateY(-35);(35);
    renderer->AddActor(lineactor);
}



/* Draw the quads of the top boundary or down boundary.
 * Compared to drawUpOrDownBoundaryQuadLine_method1, this method can handle arbitrary interpolate level.
 * Method getAllSubQuadsVertexesOnBoundary() should be called before drawUpOrDownBoundaryQuadLine_method2.
 * up: set to true if want the up boundary shown.
 * index1: set the index of quad want to shown, if this value beyond [0,quadNum], all the up boundary quads will be shown.
 * down: set to true if want the down boundary shown.
 * index2: set the index of quad want to shown, if this value beyond [0,quadNum], all the down boundary quads will be shown.
*/
void visualization::drawUpOrDownBoundaryQuadLine_method2(bool up, int index1, bool down, int index2){
    if(up){
        if(index1>=0&&index1<this->quadNum)
            drawFrameOfQuadByIndex_m2(this->up_BP, index1);
        else{
            for(unsigned i=0; i<this->quadNum; i++){
                drawFrameOfQuadByIndex_m2(this->up_BP, i);
            }
        }
    }
    if(down){
        if(index2>=0&&index2<this->quadNum)
            drawFrameOfQuadByIndex_m2(this->down_BP, index2);
        else{
            for(unsigned i=0; i<this->quadNum; i++){
                drawFrameOfQuadByIndex_m2(this->down_BP, i);
            }
        }
    }
}






