#include "visualizecrest.h"

visualizecrest::visualizecrest()
{
}


visualizecrest::visualizecrest(M3DQuadFigure* quadfig, int interpolationLevel, double *quadColor, vtkSmartPointer<vtkRenderer> renderer,
                               double rX, double rY, double rZ, double opacity){
    this->quadFig = quadfig;
    this->quadColor = quadColor;
    this->renderer = renderer;    
    int rowNums = quadFig->getRowCount();
    int colNums = quadFig->getColumnCount();
    this->crestQuadNum = rowNums*2 + (colNums-2)*2;

    toolsfunc tls;
    int step = pow((double)2, (double)interpolationLevel);

    // Along the curve between two correspondence up and down spokes.(v direction in crest interpolate method)
    tls.splitLine(0,1,step, this->subpoints_v); // Up spoke tip is consider as 0, down spoke 1. The medial crest is 0.5.

    // Along the crest curve in counter clockwise. the start from (0,0).
    tls.splitLine(0,1,step, this->subpoints_t);

    vtkSmartPointer<vtkSRep> srepfig = getSrepFig(quadfig);

    // New a crest interplator....
    this->interpolatecrestspokes = vtkSmartPointer<vtkSRepInterpolateCrestSpokesQuartic>::New();
    this->interpolatecrestspokes->SetInput(srepfig);
    this->interpolatecrestspokes->SetInterpolationLevel(interpolationLevel);
    this->interpolatecrestspokes->Update();

    this->rX = rX;
    this->rY = rY;
    this->rZ = rZ;
    this->opacity = opacity;
}



/* Draw the boundary crest in non-filled sub quads. */
void visualizecrest::drawCrestFrames(){

    vtkSmartPointer< vtkPoints > hubpos = vtkSmartPointer< vtkPoints >::New();
    vtkSmartPointer<vtkCellArray> cellarraypointsline = vtkSmartPointer<vtkCellArray>::New();
    vtkSmartPointer<vtkPolyData> polypointsline = vtkSmartPointer<vtkPolyData>::New();

    // Loop each crest quad.
    for(unsigned int i = 0; i< this->crestQuadNum; i++){
        for(unsigned int m = 0; m < subpoints_t.size()-1; m++){ //loop each sub points in t direction
            // Get the interpolated skeletal position at this point
            vtkSRep::VNLType s1 = this->interpolatecrestspokes->GetInterpolatedPoint(i,subpoints_t[m]);
            // Get the interpolated skeletal position at this point
            vtkSRep::VNLType s2 = this->interpolatecrestspokes->GetInterpolatedPoint(i,subpoints_t[m+1]);

            for(unsigned int n = 0; n < subpoints_v.size()-1; n++){ //loop each sub points in v direction
                vtkSmartPointer<vtkLine> medialsheetline = vtkSmartPointer<vtkLine>::New();

                // Get each interpolated boundary spoke direction at this point.
                vtkSRep::VNLType p0 = this->interpolatecrestspokes->GetInterpolatedSpoke(i,subpoints_t[m],subpoints_v[n]);
                vtkSRep::VNLType p1 = this->interpolatecrestspokes->GetInterpolatedSpoke(i,subpoints_t[m],subpoints_v[n+1]);
                vtkSRep::VNLType p2 = this->interpolatecrestspokes->GetInterpolatedSpoke(i,subpoints_t[m+1],subpoints_v[n+1]);
                vtkSRep::VNLType p3 = this->interpolatecrestspokes->GetInterpolatedSpoke(i,subpoints_t[m+1],subpoints_v[n]);

                vtkIdType id0 = hubpos->InsertNextPoint(s1[0] + p0[0], s1[1] + p0[1], s1[2] + p0[2]);
                vtkIdType id1 = hubpos->InsertNextPoint(s1[0] + p1[0], s1[1] + p1[1], s1[2] + p1[2]);
                vtkIdType id2 = hubpos->InsertNextPoint(s2[0] + p2[0], s2[1] + p2[1], s2[2] + p2[2]);
                vtkIdType id3 = hubpos->InsertNextPoint(s2[0] + p3[0], s2[1] + p3[1], s2[2] + p3[2]);

                // Add p0p1 into render.
                medialsheetline->GetPointIds()->SetId(0, id0);
                medialsheetline->GetPointIds()->SetId(1, id1);
                cellarraypointsline->InsertNextCell(medialsheetline);

                // Add p1p2 into render.
                medialsheetline->GetPointIds()->SetId(0, id1);
                medialsheetline->GetPointIds()->SetId(1, id2);
                cellarraypointsline->InsertNextCell(medialsheetline);

                // Add p2p3 into render.
                medialsheetline->GetPointIds()->SetId(0, id2);
                medialsheetline->GetPointIds()->SetId(1, id3);
                cellarraypointsline->InsertNextCell(medialsheetline);

                // Add p0p3 into render.
                medialsheetline->GetPointIds()->SetId(0, id0);
                medialsheetline->GetPointIds()->SetId(1, id3);
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
    //lineactor->RotateY(95);
    //lineactor->RotateX(30);
    //lineactor->RotateY(-65);
    lineactor->RotateX(rX);
    lineactor->RotateY(rY);
    lineactor->RotateZ(rZ);
    lineactor->GetProperty()->SetOpacity(this->opacity);
    this->renderer->AddActor(lineactor);
}



/* Draw the boundary crest in non-filled sub quads. */
void visualizecrest::drawCrestFramesByQuadIndex(int quadIndex, double *quadColor){

    vtkSmartPointer< vtkPoints > hubpos = vtkSmartPointer< vtkPoints >::New();
    vtkSmartPointer<vtkCellArray> cellarraypointsline = vtkSmartPointer<vtkCellArray>::New();
    vtkSmartPointer<vtkPolyData> polypointsline = vtkSmartPointer<vtkPolyData>::New();

    // Loop each crest quad.
    int i = quadIndex;
    for(unsigned int m = 0; m < subpoints_t.size()-1; m++){ //loop each sub points in t direction
        // Get the interpolated skeletal position at this point
        vtkSRep::VNLType s1 = this->interpolatecrestspokes->GetInterpolatedPoint(quadIndex,subpoints_t[m]);
        // Get the interpolated skeletal position at this point
        vtkSRep::VNLType s2 = this->interpolatecrestspokes->GetInterpolatedPoint(quadIndex,subpoints_t[m+1]);

        for(unsigned int n = 0; n < subpoints_v.size()-1; n++){ //loop each sub points in v direction
            vtkSmartPointer<vtkLine> medialsheetline = vtkSmartPointer<vtkLine>::New();

            // Get each interpolated boundary spoke direction at this point.
            vtkSRep::VNLType p0 = this->interpolatecrestspokes->GetInterpolatedSpoke(quadIndex,subpoints_t[m],subpoints_v[n]);
            vtkSRep::VNLType p1 = this->interpolatecrestspokes->GetInterpolatedSpoke(quadIndex,subpoints_t[m],subpoints_v[n+1]);
            vtkSRep::VNLType p2 = this->interpolatecrestspokes->GetInterpolatedSpoke(quadIndex,subpoints_t[m+1],subpoints_v[n+1]);
            vtkSRep::VNLType p3 = this->interpolatecrestspokes->GetInterpolatedSpoke(quadIndex,subpoints_t[m+1],subpoints_v[n]);

            vtkIdType id0 = hubpos->InsertNextPoint(s1[0] + p0[0], s1[1] + p0[1], s1[2] + p0[2]);
            vtkIdType id1 = hubpos->InsertNextPoint(s1[0] + p1[0], s1[1] + p1[1], s1[2] + p1[2]);
            vtkIdType id2 = hubpos->InsertNextPoint(s2[0] + p2[0], s2[1] + p2[1], s2[2] + p2[2]);
            vtkIdType id3 = hubpos->InsertNextPoint(s2[0] + p3[0], s2[1] + p3[1], s2[2] + p3[2]);

            // Add p0p1 into render.
            medialsheetline->GetPointIds()->SetId(0, id0);
            medialsheetline->GetPointIds()->SetId(1, id1);
            cellarraypointsline->InsertNextCell(medialsheetline);

            // Add p1p2 into render.
            medialsheetline->GetPointIds()->SetId(0, id1);
            medialsheetline->GetPointIds()->SetId(1, id2);
            cellarraypointsline->InsertNextCell(medialsheetline);

            // Add p2p3 into render.
            medialsheetline->GetPointIds()->SetId(0, id2);
            medialsheetline->GetPointIds()->SetId(1, id3);
            cellarraypointsline->InsertNextCell(medialsheetline);

            // Add p0p3 into render.
            medialsheetline->GetPointIds()->SetId(0, id0);
            medialsheetline->GetPointIds()->SetId(1, id3);
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
    lineactor->RotateX(rX);
    lineactor->RotateY(rY);
    lineactor->RotateZ(rZ);
    //lineactor->GetProperty()->SetOpacity(0.15);
    this->renderer->AddActor(lineactor);
}


/* Draw the boundary crest in non-filled sub quads. */
void visualizecrest::drawCrestQuadByQuadIndex(int quadIndex, double *quadColor){

    vtkSmartPointer< vtkPoints > hubpos = vtkSmartPointer< vtkPoints >::New();
    vtkSmartPointer<vtkCellArray> quads = vtkSmartPointer<vtkCellArray>::New();

    // For quad i.
    int i = quadIndex;
    for(unsigned int m = 0; m < subpoints_t.size()-1; m++){ //loop each sub points in t direction
        // Get the interpolated skeletal position at this point
        vtkSRep::VNLType s1 = this->interpolatecrestspokes->GetInterpolatedPoint(quadIndex,subpoints_t[m]);
        // Get the interpolated skeletal position at this point
        vtkSRep::VNLType s2 = this->interpolatecrestspokes->GetInterpolatedPoint(quadIndex,subpoints_t[m+1]);

        for(unsigned int n = 0; n < subpoints_v.size()-1; n++){ //loop each sub points in v direction
            //Create four points (must be in counter clockwise order)
            vtkSmartPointer<vtkQuad> quad = vtkSmartPointer<vtkQuad>::New();

            // Get each interpolated boundary spoke direction at this point.
            vtkSRep::VNLType p0 = this->interpolatecrestspokes->GetInterpolatedSpoke(quadIndex,subpoints_t[m],subpoints_v[n]);
            vtkSRep::VNLType p1 = this->interpolatecrestspokes->GetInterpolatedSpoke(quadIndex,subpoints_t[m],subpoints_v[n+1]);
            vtkSRep::VNLType p2 = this->interpolatecrestspokes->GetInterpolatedSpoke(quadIndex,subpoints_t[m+1],subpoints_v[n+1]);
            vtkSRep::VNLType p3 = this->interpolatecrestspokes->GetInterpolatedSpoke(quadIndex,subpoints_t[m+1],subpoints_v[n]);

            vtkIdType id0 = hubpos->InsertNextPoint(s1[0] + p0[0], s1[1] + p0[1], s1[2] + p0[2]);
            vtkIdType id1 = hubpos->InsertNextPoint(s1[0] + p1[0], s1[1] + p1[1], s1[2] + p1[2]);
            vtkIdType id2 = hubpos->InsertNextPoint(s2[0] + p2[0], s2[1] + p2[1], s2[2] + p2[2]);
            vtkIdType id3 = hubpos->InsertNextPoint(s2[0] + p3[0], s2[1] + p3[1], s2[2] + p3[2]);

            quad->GetPointIds()->SetId(0, id0);
            quad->GetPointIds()->SetId(1, id1);
            quad->GetPointIds()->SetId(2, id2);
            quad->GetPointIds()->SetId(3, id3);

            quads->InsertNextCell(quad);
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

    quadActor->RotateX(rX);
    quadActor->RotateY(rY);
    quadActor->RotateZ(rZ);
    //quadActor->GetProperty()->SetOpacity(0.15);
    quadActor->GetProperty()->SetOpacity(0.7);

    renderer->AddActor(quadActor);
}



/* Draw the boundary crest in non-filled sub quads. */
void visualizecrest::drawCrestQuadWithDiffColor(int quadNum){

    vtkSmartPointer< vtkPoints > hubpos = vtkSmartPointer< vtkPoints >::New();
    vtkSmartPointer<vtkCellArray> quads = vtkSmartPointer<vtkCellArray>::New();

    // For quad i.
    for (unsigned int i = 0; i < quadNum; i++){
    for(unsigned int m = 0; m < subpoints_t.size()-1; m++){ //loop each sub points in t direction
        // Get the interpolated skeletal position at this point
        vtkSRep::VNLType s1 = this->interpolatecrestspokes->GetInterpolatedPoint(i,subpoints_t[m]);
        // Get the interpolated skeletal position at this point
        vtkSRep::VNLType s2 = this->interpolatecrestspokes->GetInterpolatedPoint(i,subpoints_t[m+1]);

        for(unsigned int n = 0; n < subpoints_v.size()-1; n++){ //loop each sub points in v direction
            //Create four points (must be in counter clockwise order)
            vtkSmartPointer<vtkQuad> quad = vtkSmartPointer<vtkQuad>::New();

            // Get each interpolated boundary spoke direction at this point.
            vtkSRep::VNLType p0 = this->interpolatecrestspokes->GetInterpolatedSpoke(i,subpoints_t[m],subpoints_v[n]);
            vtkSRep::VNLType p1 = this->interpolatecrestspokes->GetInterpolatedSpoke(i,subpoints_t[m],subpoints_v[n+1]);
            vtkSRep::VNLType p2 = this->interpolatecrestspokes->GetInterpolatedSpoke(i,subpoints_t[m+1],subpoints_v[n+1]);
            vtkSRep::VNLType p3 = this->interpolatecrestspokes->GetInterpolatedSpoke(i,subpoints_t[m+1],subpoints_v[n]);

            vtkIdType id0 = hubpos->InsertNextPoint(s1[0] + p0[0], s1[1] + p0[1], s1[2] + p0[2]);
            vtkIdType id1 = hubpos->InsertNextPoint(s1[0] + p1[0], s1[1] + p1[1], s1[2] + p1[2]);
            vtkIdType id2 = hubpos->InsertNextPoint(s2[0] + p2[0], s2[1] + p2[1], s2[2] + p2[2]);
            vtkIdType id3 = hubpos->InsertNextPoint(s2[0] + p3[0], s2[1] + p3[1], s2[2] + p3[2]);

            quad->GetPointIds()->SetId(0, id0);
            quad->GetPointIds()->SetId(1, id1);
            quad->GetPointIds()->SetId(2, id2);
            quad->GetPointIds()->SetId(3, id3);

            quads->InsertNextCell(quad);
        }
    }
    }

    // Create a polydata to store everything in
    vtkSmartPointer<vtkPolyData> polydataQuad = vtkSmartPointer<vtkPolyData>::New();

    // Add the points and quads to the dataset
    polydataQuad->SetPoints(hubpos);
    polydataQuad->SetPolys(quads);

    // Find min and max z
    double minz = 0;
    double maxz = quadNum;

    // Create the color map
    vtkSmartPointer<vtkLookupTable> colorLookupTable = vtkSmartPointer<vtkLookupTable>::New();
    colorLookupTable->SetTableRange(minz, maxz);
    colorLookupTable->Build();

    // Generate the colors for each point based on the color map
    vtkSmartPointer<vtkUnsignedCharArray> colors = vtkSmartPointer<vtkUnsignedCharArray>::New();
    colors->SetNumberOfComponents(3);
    colors->SetName("Colors");

    for(int k = 0; k < quadNum; k++) {

        // Get a color correspondence to the z coordinate value from the color table
        double dcolor[3];
        colorLookupTable->GetColor(k, dcolor); //p[2] is the z coordinate

        unsigned char color[3];
        for(unsigned int j = 0; j < 3; j++) {
            color[j] = static_cast<unsigned char>(255.0 * dcolor[j]);
        }

        colors->InsertNextTupleValue(color);
    }

    polydataQuad->GetPointData()->SetScalars(colors);


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
    //quadActor->GetProperty()->SetOpacity(0.15);
    //quadActor->GetProperty()->SetOpacity(0.7);

    renderer->AddActor(quadActor);
}



/* Draw the boundary crest in filled sub quads. Looks like the surface. */
void visualizecrest::drawCrestQuads(){

    vtkSmartPointer< vtkPoints > hubpos = vtkSmartPointer< vtkPoints >::New();
    vtkSmartPointer<vtkCellArray> quads = vtkSmartPointer<vtkCellArray>::New();

    // Loop each crest quad.
    for(unsigned int i = 0; i< this->crestQuadNum; i++){
        for(unsigned int m = 0; m < subpoints_t.size()-1; m++){ //loop each sub points in t direction
            // Get the interpolated skeletal position at this point
            vtkSRep::VNLType s1 = this->interpolatecrestspokes->GetInterpolatedPoint(i,subpoints_t[m]);
            // Get the interpolated skeletal position at this point
            vtkSRep::VNLType s2 = this->interpolatecrestspokes->GetInterpolatedPoint(i,subpoints_t[m+1]);

            for(unsigned int n = 0; n < subpoints_v.size()-1; n++){ //loop each sub points in v direction

                vtkSmartPointer<vtkQuad> quad = vtkSmartPointer<vtkQuad>::New();

                // Get each interpolated boundary spoke direction at this point.
                vtkSRep::VNLType p0 = this->interpolatecrestspokes->GetInterpolatedSpoke(i,subpoints_t[m],subpoints_v[n]);
                vtkSRep::VNLType p1 = this->interpolatecrestspokes->GetInterpolatedSpoke(i,subpoints_t[m],subpoints_v[n+1]);
                vtkSRep::VNLType p2 = this->interpolatecrestspokes->GetInterpolatedSpoke(i,subpoints_t[m+1],subpoints_v[n+1]);
                vtkSRep::VNLType p3 = this->interpolatecrestspokes->GetInterpolatedSpoke(i,subpoints_t[m+1],subpoints_v[n]);

                vtkIdType id0 = hubpos->InsertNextPoint(s1[0] + p0[0], s1[1] + p0[1], s1[2] + p0[2]);
                vtkIdType id1 = hubpos->InsertNextPoint(s1[0] + p1[0], s1[1] + p1[1], s1[2] + p1[2]);
                vtkIdType id2 = hubpos->InsertNextPoint(s2[0] + p2[0], s2[1] + p2[1], s2[2] + p2[2]);
                vtkIdType id3 = hubpos->InsertNextPoint(s2[0] + p3[0], s2[1] + p3[1], s2[2] + p3[2]);

                quad->GetPointIds()->SetId(0, id0);
                quad->GetPointIds()->SetId(1, id1);
                quad->GetPointIds()->SetId(2, id2);
                quad->GetPointIds()->SetId(3, id3);

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

    quadActor->GetProperty()->SetColor(this->quadColor);
    quadActor->RotateZ(90);
    quadActor->RotateX(-115);
    quadActor->RotateY(35);
    this->renderer->AddActor(quadActor);
}





/* Connect the diagnoal line of each (sub)quad, use the smaller area one.
*/
void visualizecrest::connectDiagnoal(){

    vtkSmartPointer< vtkPoints > hubpos = vtkSmartPointer< vtkPoints >::New();
    vtkSmartPointer<vtkCellArray> cellarraypointsline = vtkSmartPointer<vtkCellArray>::New();
    vtkSmartPointer<vtkPolyData> polypointsline = vtkSmartPointer<vtkPolyData>::New();

    // Loop each crest quad.
    for(unsigned int i = 0; i< this->crestQuadNum; i++){
        for(unsigned int m = 0; m < subpoints_t.size()-1; m++){ //loop each sub points in t direction
            // Get the interpolated skeletal position at this point
            vtkSRep::VNLType s1 = this->interpolatecrestspokes->GetInterpolatedPoint(i,subpoints_t[m]);
            // Get the interpolated skeletal position at this point
            vtkSRep::VNLType s2 = this->interpolatecrestspokes->GetInterpolatedPoint(i,subpoints_t[m+1]);

            for(unsigned int n = 0; n < subpoints_v.size()-1; n++){ //loop each sub points in v direction

                // Get each interpolated boundary spoke direction at this point.
                vtkSRep::VNLType p0 = this->interpolatecrestspokes->GetInterpolatedSpoke(i,subpoints_t[m],subpoints_v[n]);
                vtkSRep::VNLType p1 = this->interpolatecrestspokes->GetInterpolatedSpoke(i,subpoints_t[m],subpoints_v[n+1]);
                vtkSRep::VNLType p2 = this->interpolatecrestspokes->GetInterpolatedSpoke(i,subpoints_t[m+1],subpoints_v[n+1]);
                vtkSRep::VNLType p3 = this->interpolatecrestspokes->GetInterpolatedSpoke(i,subpoints_t[m+1],subpoints_v[n]);

                // Find the diagnoal which make min area of the quad.
                Vector3D point[4];
                point[0].set(s1[0] + p0[0], s1[1] + p0[1], s1[2] + p0[2]);
                point[1].set(s1[0] + p1[0], s1[1] + p1[1], s1[2] + p1[2]);
                point[2].set(s2[0] + p2[0], s2[1] + p2[1], s2[2] + p2[2]);
                point[3].set(s2[0] + p3[0], s2[1] + p3[1], s2[2] + p3[2]);

                // If true, use diagonalp0p2; if false use diagonalp1p3
                bool p0p2 = tools.quadAreaMin(point);
                vtkSmartPointer<vtkLine> medialsheetline = vtkSmartPointer<vtkLine>::New();
                vtkIdType id0, id1;
                if(p0p2){ // Connect p0p2 as diagnoal
                    id0 = hubpos->InsertNextPoint(s1[0] + p0[0], s1[1] + p0[1], s1[2] + p0[2]);
                    id1 = hubpos->InsertNextPoint(s2[0] + p2[0], s2[1] + p2[1], s2[2] + p2[2]);
                }
                else{ // Connect p1p3 as diagnoal
                    id0 = hubpos->InsertNextPoint(s1[0] + p1[0], s1[1] + p1[1], s1[2] + p1[2]);
                    id1 = hubpos->InsertNextPoint(s2[0] + p3[0], s2[1] + p3[1], s2[2] + p3[2]);
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
    lineactor->RotateY(this->rZ);
    lineactor->GetProperty()->SetOpacity(this->opacity);
    renderer->AddActor(lineactor);
}



/* Only draw the middle one crest spokes.
*/
void visualizecrest::drawMidOneCrestSpokes(double rX, double rY, double rZ, int linewidth, double * color, bool showspoketail){

    vtkSmartPointer< vtkPoints > hubpos = vtkSmartPointer< vtkPoints >::New();
    vtkSmartPointer<vtkCellArray> cellarraypointsline = vtkSmartPointer<vtkCellArray>::New();
    vtkSmartPointer<vtkPolyData> polypointsline = vtkSmartPointer<vtkPolyData>::New();

    vtkSmartPointer< vtkPoints > skeletalpoint = vtkSmartPointer< vtkPoints >::New();

    // Loop each crest quad.
    for(unsigned int i = 0; i< this->crestQuadNum; i++){
        // Get the interpolated skeletal position at this point
        vtkSRep::VNLType s1 = this->interpolatecrestspokes->GetInterpolatedPoint(i,0);

        vtkSmartPointer<vtkLine> medialsheetline = vtkSmartPointer<vtkLine>::New();

        // Get each interpolated boundary spoke direction at this point.
        vtkSRep::VNLType bv = this->interpolatecrestspokes->GetInterpolatedSpoke(i, 0, 0.5);

        vtkIdType id0 = hubpos->InsertNextPoint(s1[0], s1[1], s1[2]);
        vtkIdType id1 = hubpos->InsertNextPoint(s1[0]+bv[0], s1[1]+bv[1], s1[2]+bv[2]);

        // Add to render
        medialsheetline->GetPointIds()->SetId(0, id0);
        medialsheetline->GetPointIds()->SetId(1, id1);
        cellarraypointsline->InsertNextCell(medialsheetline);

        skeletalpoint->InsertNextPoint(s1[0], s1[1], s1[2]);
    }

    polypointsline->SetPoints(hubpos);
    polypointsline->SetLines(cellarraypointsline);

    vtkSmartPointer<vtkPolyDataMapper> pointlinemapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    pointlinemapper->SetInput(polypointsline);
    vtkSmartPointer<vtkActor>  lineactor = vtkActor::New();
    lineactor->SetMapper(pointlinemapper);
    lineactor->GetProperty()->SetLineWidth(linewidth);
    lineactor->GetProperty()->SetColor(color);//1,1,0.1//red    // (1,0,0) is red, (0,1,0) is green, (1,1,1) is white.
    lineactor->RotateX(rX);
    lineactor->RotateY(rY);
    lineactor->RotateZ(rZ);
    //lineactor->GetProperty()->SetOpacity(0.15);
    this->renderer->AddActor(lineactor);

    //If choosen to show the spoke tails
    if(showspoketail){
        visualization vs;
        color[0]= 0;
        color[1]= 0.9;
        color[2] = 0;
        vs.showPointBall(skeletalpoint, color, 0.0012, this->renderer, this->rX, this->rY, this->rZ, this->opacity); //0.0016
    }
}




/* Draw the crest spokes.
 * showCrestLine: if true, will show the crest line.
*/
void visualizecrest::drawCrestSpokes(){

    vtkSmartPointer< vtkPoints > hubpos = vtkSmartPointer< vtkPoints >::New();
    vtkSmartPointer<vtkCellArray> cellarraypointsline = vtkSmartPointer<vtkCellArray>::New();
    vtkSmartPointer<vtkPolyData> polypointsline = vtkSmartPointer<vtkPolyData>::New();

    // Loop each crest quad.
    for(unsigned int i = 0; i< this->crestQuadNum; i++){
        for(unsigned int m = 0; m < subpoints_t.size(); m++){ //loop each sub points in t direction
            // Get the interpolated skeletal position at this point
            vtkSRep::VNLType s1 = this->interpolatecrestspokes->GetInterpolatedPoint(i,subpoints_t[m]);

            for(unsigned int n = 0; n < subpoints_v.size(); n++){ //loop each sub points in v direction
                vtkSmartPointer<vtkLine> medialsheetline = vtkSmartPointer<vtkLine>::New();

                // Get each interpolated boundary spoke direction at this point.
                vtkSRep::VNLType bv = this->interpolatecrestspokes->GetInterpolatedSpoke(i,subpoints_t[m],subpoints_v[n]);

                vtkIdType id0 = hubpos->InsertNextPoint(s1[0], s1[1], s1[2]);
                vtkIdType id1 = hubpos->InsertNextPoint(s1[0]+bv[0], s1[1]+bv[1], s1[2]+bv[2]);

                // Add to render
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
    lineactor->GetProperty()->SetColor(this->quadColor);// 1,1,0.1  (0,1,0) is green, (1,1,1) is white.
    lineactor->RotateX(rX);
    lineactor->RotateY(rY);
    lineactor->RotateZ(rZ);
    lineactor->GetProperty()->SetOpacity(this->opacity);
    this->renderer->AddActor(lineactor);
}



/* Draw the crest spokes.
 * showCrestLine: if true, will show the crest line.
*/
void visualizecrest::drawCrestSpokesByQuadIndex(int quadIndex, double *quadColor){

    vtkSmartPointer< vtkPoints > hubpos = vtkSmartPointer< vtkPoints >::New();
    vtkSmartPointer<vtkCellArray> cellarraypointsline = vtkSmartPointer<vtkCellArray>::New();
    vtkSmartPointer<vtkPolyData> polypointsline = vtkSmartPointer<vtkPolyData>::New();

    // Loop each crest quad.
    int i = quadIndex;
        for(unsigned int m = 0; m < subpoints_t.size(); m++){ //loop each sub points in t direction
            // Get the interpolated skeletal position at this point
            vtkSRep::VNLType s1 = this->interpolatecrestspokes->GetInterpolatedPoint(i,subpoints_t[m]);

            for(unsigned int n = 0; n < subpoints_v.size(); n++){ //loop each sub points in v direction
                vtkSmartPointer<vtkLine> medialsheetline = vtkSmartPointer<vtkLine>::New();

                // Get each interpolated boundary spoke direction at this point.
                vtkSRep::VNLType bv = this->interpolatecrestspokes->GetInterpolatedSpoke(i,subpoints_t[m],subpoints_v[n]);

                vtkIdType id0 = hubpos->InsertNextPoint(s1[0], s1[1], s1[2]);
                vtkIdType id1 = hubpos->InsertNextPoint(s1[0]+bv[0], s1[1]+bv[1], s1[2]+bv[2]);

                // Add to render
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
    lineactor->RotateX(rX);
    lineactor->RotateY(rY);
    lineactor->RotateZ(rZ);
    lineactor->GetProperty()->SetOpacity(this->opacity);
    this->renderer->AddActor(lineactor);
}

/* Only draw the middle crest spokes. Without showing the up and down crest spoke.
*/
void visualizecrest::drawMiddleCrestSpokesByQuadIndex(int quadIndex, double *quadColor, bool showspoketail){

    vtkSmartPointer< vtkPoints > hubpos = vtkSmartPointer< vtkPoints >::New();
    vtkSmartPointer<vtkCellArray> cellarraypointsline = vtkSmartPointer<vtkCellArray>::New();
    vtkSmartPointer<vtkPolyData> polypointsline = vtkSmartPointer<vtkPolyData>::New();

    // Loop each crest quad. Get the interpolated skeletal position at this point
    vtkSRep::VNLType s1 = this->interpolatecrestspokes->GetInterpolatedPoint(quadIndex,0);

    vtkSmartPointer<vtkLine> medialsheetline = vtkSmartPointer<vtkLine>::New();

    // Get each interpolated boundary spoke direction at this point.
    vtkSRep::VNLType bv = this->interpolatecrestspokes->GetInterpolatedSpoke(quadIndex,0, 0.5);

    vtkIdType id0 = hubpos->InsertNextPoint(s1[0], s1[1], s1[2]);
    vtkIdType id1 = hubpos->InsertNextPoint(s1[0]+bv[0], s1[1]+bv[1], s1[2]+bv[2]);

    // Add to render
    medialsheetline->GetPointIds()->SetId(0, id0);
    medialsheetline->GetPointIds()->SetId(1, id1);
    cellarraypointsline->InsertNextCell(medialsheetline);

    polypointsline->SetPoints(hubpos);
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
    lineactor->GetProperty()->SetOpacity(this->opacity);
    this->renderer->AddActor(lineactor);

}



/* Draw the volume of crest quad.
*/
void visualizecrest::drawCrestVolumesByQuadIndex(int quadIndex, double *quadColor){

    // Draw the crest quad
    drawCrestQuadByQuadIndex(quadIndex, quadColor);

    // Draw the quads of up and down boundary to skeletal point, (two side)
    vtkSmartPointer< vtkPoints > hubpos = vtkSmartPointer< vtkPoints >::New();
    vtkSmartPointer<vtkCellArray> quads = vtkSmartPointer<vtkCellArray>::New();

    // Loop each crest quad.
    for(unsigned int m = 0; m < subpoints_t.size()-1; m++){ //loop each sub points in t direction
        // Get the interpolated skeletal position at this point
        vtkSRep::VNLType s1 = this->interpolatecrestspokes->GetInterpolatedPoint(quadIndex,subpoints_t[m]);
        // Get the interpolated skeletal position at this point
        vtkSRep::VNLType s2 = this->interpolatecrestspokes->GetInterpolatedPoint(quadIndex,subpoints_t[m+1]);

        vtkSmartPointer<vtkQuad> quad = vtkSmartPointer<vtkQuad>::New();

        // Get each interpolated boundary spoke direction at this point.
        vtkSRep::VNLType p0 = this->interpolatecrestspokes->GetInterpolatedSpoke(quadIndex,subpoints_t[m],0);
        vtkSRep::VNLType p3 = this->interpolatecrestspokes->GetInterpolatedSpoke(quadIndex,subpoints_t[m+1],0);

        vtkIdType id0 = hubpos->InsertNextPoint(s1[0] + p0[0], s1[1] + p0[1], s1[2] + p0[2]);
        vtkIdType id1 = hubpos->InsertNextPoint(s1[0], s1[1], s1[2]);
        vtkIdType id2 = hubpos->InsertNextPoint(s2[0], s2[1], s2[2]);
        vtkIdType id3 = hubpos->InsertNextPoint(s2[0] + p3[0], s2[1] + p3[1], s2[2] + p3[2]);

        quad->GetPointIds()->SetId(0, id0);
        quad->GetPointIds()->SetId(1, id1);
        quad->GetPointIds()->SetId(2, id2);
        quad->GetPointIds()->SetId(3, id3);

        quads->InsertNextCell(quad);//the size of the cellarray is the number of the quads, its 24.


        p0 = this->interpolatecrestspokes->GetInterpolatedSpoke(quadIndex,subpoints_t[m],1);
        p3 = this->interpolatecrestspokes->GetInterpolatedSpoke(quadIndex,subpoints_t[m+1],1);

        id0 = hubpos->InsertNextPoint(s1[0] + p0[0], s1[1] + p0[1], s1[2] + p0[2]);
        id1 = hubpos->InsertNextPoint(s1[0], s1[1], s1[2]);
        id2 = hubpos->InsertNextPoint(s2[0], s2[1], s2[2]);
        id3 = hubpos->InsertNextPoint(s2[0] + p3[0], s2[1] + p3[1], s2[2] + p3[2]);

        quad->GetPointIds()->SetId(0, id0);
        quad->GetPointIds()->SetId(1, id1);
        quad->GetPointIds()->SetId(2, id2);
        quad->GetPointIds()->SetId(3, id3);

        quads->InsertNextCell(quad);//the size of the cellarray is the number of the quads, its 24.
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

    quadActor->GetProperty()->SetColor(this->quadColor);
    quadActor->RotateX(rX);
    quadActor->RotateY(rY);
    quadActor->RotateZ(rZ);
    quadActor->GetProperty()->SetOpacity(0.7);
    this->renderer->AddActor(quadActor);
}





/* Draw the crest spokes with beautiful different colors.
 * showCrestLine: if true, will show the crest line.
*/
void visualizecrest::drawCorrespondenceCrestSpokes(){

    vtkSmartPointer< vtkPoints > hubpos = vtkSmartPointer< vtkPoints >::New();
    vtkSmartPointer<vtkCellArray> cellarraypointsline = vtkSmartPointer<vtkCellArray>::New();
    vtkSmartPointer<vtkPolyData> linesPolydata = vtkSmartPointer<vtkPolyData>::New();

    // Loop each crest quad.
    for(unsigned int i = 0; i< this->crestQuadNum; i++){
        for(unsigned int m = 0; m < subpoints_t.size(); m++){ //loop each sub points in t direction
            // Get the interpolated skeletal position at this point
            vtkSRep::VNLType s1 = this->interpolatecrestspokes->GetInterpolatedPoint(i,subpoints_t[m]);

            for(unsigned int n = 0; n < subpoints_v.size(); n++){ //loop each sub points in v direction
                vtkSmartPointer<vtkLine> medialsheetline = vtkSmartPointer<vtkLine>::New();

                // Get each interpolated boundary spoke direction at this point.
                vtkSRep::VNLType bv = this->interpolatecrestspokes->GetInterpolatedSpoke(i,subpoints_t[m],subpoints_v[n]);

                vtkIdType id0 = hubpos->InsertNextPoint(s1[0], s1[1], s1[2]);
                vtkIdType id1 = hubpos->InsertNextPoint(s1[0]+bv[0], s1[1]+bv[1], s1[2]+bv[2]);

                // Add to render
                medialsheetline->GetPointIds()->SetId(0, id0);
                medialsheetline->GetPointIds()->SetId(1, id1);
                cellarraypointsline->InsertNextCell(medialsheetline);
            }
        }
    }

    linesPolydata->SetPoints(hubpos);
    linesPolydata->SetLines(cellarraypointsline);

    // Find min and max z
    double minz = 0;
    double maxz = 3800; //linesPolydata->GetNumberOfPoints();
    //cout<<"-------crest-------linesPolydata->GetNumberOfPoints() is:"<<linesPolydata->GetNumberOfPoints()<<endl; //leve=2, 1400

    // Create the color map
    vtkSmartPointer<vtkLookupTable> colorLookupTable = vtkSmartPointer<vtkLookupTable>::New();
    colorLookupTable->SetTableRange(minz, maxz);
    colorLookupTable->Build();

    // Generate the colors for each point based on the color map
    vtkSmartPointer<vtkUnsignedCharArray> colors = vtkSmartPointer<vtkUnsignedCharArray>::New();
    colors->SetNumberOfComponents(3);
    colors->SetName("Colors");

    for(int i = 2400; i<maxz; i++) {

        // Get a color correspondence to the z coordinate value from the color table
        double dcolor[3];
        colorLookupTable->GetColor(i, dcolor); //p[2] is the z coordinate

        unsigned char color[3];
        for(unsigned int j = 0; j < 3; j++) {
            color[j] = static_cast<unsigned char>(255 * dcolor[j]);
        }

        colors->InsertNextTupleValue(color);
    }

    linesPolydata->GetPointData()->SetScalars(colors);



    vtkSmartPointer<vtkPolyDataMapper> pointlinemapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    pointlinemapper->SetInput(linesPolydata);
    vtkSmartPointer<vtkActor>  lineactor = vtkActor::New();
    lineactor->SetMapper(pointlinemapper);
    lineactor->GetProperty()->SetLineWidth(2);
    //lineactor->RotateY(95);
    //lineactor->RotateX(30);
    //lineactor->RotateY(-35);(35);
    this->renderer->AddActor(lineactor);
}


/* Draw the crest curve along the skeletal sheet.
*/
void visualizecrest::drawCrestCurve(double linewidth, double * color){

    vtkSmartPointer< vtkPoints > hubpos = vtkSmartPointer< vtkPoints >::New();
    vtkSmartPointer<vtkCellArray> cellarraypointsline = vtkSmartPointer<vtkCellArray>::New();
    vtkSmartPointer<vtkPolyData> polypointsline = vtkSmartPointer<vtkPolyData>::New();

    // Loop each crest quad.
    for(unsigned int i = 0; i< this->crestQuadNum; i++){
        for(unsigned int m = 0; m < subpoints_t.size()-1; m++){ //loop each sub points in t direction
            // Get the interpolated skeletal position at this point
            vtkSRep::VNLType s1 = this->interpolatecrestspokes->GetInterpolatedPoint(i,subpoints_t[m]);
            // Get the interpolated skeletal position at this point
            vtkSRep::VNLType s2 = this->interpolatecrestspokes->GetInterpolatedPoint(i,subpoints_t[m+1]);

            vtkIdType id0 = hubpos->InsertNextPoint(s1[0], s1[1], s1[2]);
            vtkIdType id1 = hubpos->InsertNextPoint(s2[0], s2[1], s2[2]);

            vtkSmartPointer<vtkLine> medialsheetline = vtkSmartPointer<vtkLine>::New();

            // Add to render.
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
    lineactor->GetProperty()->SetLineWidth(linewidth);
    lineactor->GetProperty()->SetColor(0,0.6,0);// (0,1,0) is green, (1,1,1) is white.
    lineactor->RotateX(rX);
    lineactor->RotateY(rY);
    lineactor->RotateZ(rZ);
    lineactor->GetProperty()->SetOpacity(this->opacity);
    this->renderer->AddActor(lineactor);
}



/* Draw the middle crest spokes tip lines, medial crest boundary.
*/
void visualizecrest::drawCrestBoundaryEquator(double rX, double rY, double rZ){

    vtkSmartPointer< vtkPoints > hubpos = vtkSmartPointer< vtkPoints >::New();
    vtkSmartPointer<vtkCellArray> cellarraypointsline = vtkSmartPointer<vtkCellArray>::New();
    vtkSmartPointer<vtkPolyData> polypointsline = vtkSmartPointer<vtkPolyData>::New();

    // Loop each crest quad.
    for(unsigned int i = 0; i< this->crestQuadNum; i++){
        for(unsigned int m = 0; m < subpoints_t.size()-1; m++){ //loop each sub points in t direction
            // Get the interpolated skeletal position at this point
            vtkSRep::VNLType s1 = this->interpolatecrestspokes->GetInterpolatedPoint(i,subpoints_t[m]);
            // Get the interpolated skeletal position at this point
            vtkSRep::VNLType s2 = this->interpolatecrestspokes->GetInterpolatedPoint(i,subpoints_t[m+1]);

            // Get each interpolated boundary spoke direction at this point.
            vtkSRep::VNLType s1b = this->interpolatecrestspokes->GetInterpolatedSpoke(i, subpoints_t[m], 0.5);
            vtkSRep::VNLType s2b = this->interpolatecrestspokes->GetInterpolatedSpoke(i, subpoints_t[m+1], 0.5);

            vtkIdType id0 = hubpos->InsertNextPoint(s1[0]+s1b[0], s1[1]+s1b[1], s1[2]+s1b[2]);
            vtkIdType id1 = hubpos->InsertNextPoint(s2[0]+s2b[0], s2[1]+s2b[1], s2[2]+s2b[2]);


            vtkSmartPointer<vtkLine> medialsheetline = vtkSmartPointer<vtkLine>::New();

            // Add to render.
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
    lineactor->GetProperty()->SetLineWidth(2);
    lineactor->GetProperty()->SetColor(0,0,1);// (0,1,0) is green, (1,1,1) is white.
    lineactor->RotateX(rX);
    lineactor->RotateY(rY);
    lineactor->RotateZ(rZ);
    this->renderer->AddActor(lineactor);
}



/* Given a quadfig, generate its srepfig.*/
vtkSmartPointer<vtkSRep> visualizecrest::getSrepFig(M3DQuadFigure* quadfig){
    vtkSmartPointer<vtkSRep> m_Srep = vtkSmartPointer<vtkSRep>::New();

    vtkSmartPointer< vtkPoints > hubpos = vtkSmartPointer< vtkPoints >::New();
    vtkSmartPointer<vtkCellArray> cellarray = vtkSmartPointer<vtkCellArray>::New();

    vtkSRep::VectorSRepIdsType pointsIds;

    vtkSRep::RadiusVectorType allradius;
    vtkSRep::SpokesVectorType allspokes;

    if(quadfig){

        for(int u = 0; u < quadfig->getRowCount(); u++){
            pointsIds.push_back(vtkSRep::VectorIdsType());
            for(int v = 0; v < quadfig->getColumnCount(); v++){

                M3DQuadPrimitive* prim0 = dynamic_cast<M3DQuadPrimitive*>(quadfig->getPrimitivePtr(u, v));

                Vector3D x = prim0->getX();
                Vector3D u0 = prim0->getU0();
                Vector3D u1 = prim0->getU1();

                vtkSRep::VectorVNLType vnlspokes;
                vtkSRep::VNLType s(3);
                s[0] = u0.getX();
                s[1] = u0.getY();
                s[2] = u0.getZ();
                s = s.normalize();
                vnlspokes.push_back(s);

                s[0] = u1.getX();
                s[1] = u1.getY();
                s[2] = u1.getZ();
                s = s.normalize();
                vnlspokes.push_back(s);

                vtkSRep::VectorDoubleType radius;
                radius.push_back(prim0->getR0());
                radius.push_back(prim0->getR1());

                if(u == 0 || u == quadfig->getRowCount() - 1 || v == 0 || v == quadfig->getColumnCount() - 1){

                    M3DQuadEndPrimitive* prim0 = dynamic_cast<M3DQuadEndPrimitive*>(quadfig->getPrimitivePtr(u, v));
                    Vector3D uend = prim0->getUEnd();

                    s[0] = uend.getX();
                    s[1] = uend.getY();
                    s[2] = uend.getZ();
                    s = s.normalize();
                    vnlspokes.push_back(s);

                    radius.push_back(prim0->getREnd());

                }

                vtkIdType id = hubpos->InsertNextPoint(x.getX(), x.getY(), x.getZ());
                //cout<<quadfig->getPrimitiveID(u,v)<<"-"<<id<<endl;
                pointsIds[u].push_back(id);

                allspokes.push_back(vnlspokes);
                allradius.push_back(radius);
            }
        }

        for(unsigned i = 0; i < pointsIds.size() - 1; i++){
             for(unsigned j = 0; j < pointsIds[i].size() - 1; j++){

                 vtkSmartPointer<vtkQuad> quad = vtkSmartPointer<vtkQuad>::New();
                 quad->GetPointIds()->SetId(0, pointsIds[i][j]);
                 quad->GetPointIds()->SetId(1, pointsIds[i+1][j]);
                 quad->GetPointIds()->SetId(2, pointsIds[i+1][j+1]);
                 quad->GetPointIds()->SetId(3, pointsIds[i][j+1]);

                 //quad->Print(cout);

                 cellarray->InsertNextCell(quad);
             }
         }

        m_Srep->SetPoints(hubpos);
        m_Srep->SetPolys(cellarray);
        m_Srep->SetAllSpokes(allspokes);
        m_Srep->SetAllRadius(allradius);
        const float *color = quadfig->getColor();

        m_Srep->SetColor(color[0], color[1], color[2]);

        allspokes.clear();
        allradius.clear();
        pointsIds.clear();
    }

    return m_Srep;
}
