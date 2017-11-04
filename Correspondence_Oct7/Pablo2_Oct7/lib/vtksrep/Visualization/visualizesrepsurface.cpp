/* Draw the implied surface of given s-rep. Draw the boundary surface of object in vertexes (points cloud).
 * Liyun Tu
 * May 20, 2014
*/

#include "visualizesrepsurface.h"

visualizesrepsurface::visualizesrepsurface()
{
}

visualizesrepsurface::visualizesrepsurface(double rX, double rY, double rZ, double opacity)
{
    this->rX = rX;
    this->rY = rY;
    this->rZ = rZ;
    this->opacity = opacity;
}


/* Get all the points on this srep's implied surface using given interpolate degree. */
vtkSmartPointer< vtkPoints > visualizesrepsurface::getSurfacePointsSet(const char* srepfilename, int interpolationLevel){

    // Read in this s-rep
    M3DQuadFigure* quadFig = tls.GetQuadFigure(srepfilename);

    vtkSmartPointer< vtkPoints > surfacePoints = vtkSmartPointer< vtkPoints >::New();

    // Get up and down surface points
    getUpAndDownSurfacePoints(quadFig, surfacePoints, interpolationLevel);

    // Get crest surface points
    getCrestSurfacePoints(quadFig, surfacePoints, interpolationLevel);

    return surfacePoints;
}



void visualizesrepsurface::getUpAndDownSurfacePoints(M3DQuadFigure* quadFig, vtkSmartPointer< vtkPoints > pos, int interpolationLevel){
    int rowNums = quadFig->getRowCount();
    int colNums = quadFig->getColumnCount();
    int quadNum = (rowNums-1)*(colNums-1);//quadNum means the quad numbers on skeletal sheet.
    int step = pow((double)2, (double)interpolationLevel);
    int subQuadPointsNum = (step+1)*(step+1);

    // Get interpolated u & v coordinate, storing in interpolatedU & interpolatedV
    VectorQuadPoint interpolatedU;
    VectorQuadPoint interpolatedV;
    getInterpolateUVCoordinate(interpolatedU, interpolatedV, rowNums, colNums, step);

    // Get all the interpolated boundary points (on the surface)
    M3DQuadInterpolater *tpm = new M3DQuadInterpolater(quadFig);
    Vector3D point;
    for(unsigned int i = 0; i < quadNum; i++){
        //cout<<"Currently drawing quad: "<<i<<endl;
        for(unsigned int j = 0; j < subQuadPointsNum; j++){
            // Get up boundary points
            point = tpm->interpolateQuadSpoke(quadFig,interpolatedU[i][j],interpolatedV[i][j],0)->getB();
            pos->InsertNextPoint(point.getX(),point.getY(),point.getZ());

            // Get down boundary points
            point = tpm->interpolateQuadSpoke(quadFig,interpolatedU[i][j],interpolatedV[i][j],1)->getB();
            pos->InsertNextPoint(point.getX(),point.getY(),point.getZ());
        }
    }
    delete tpm;
}


void visualizesrepsurface::getUpOrDownSurfacePoints(M3DQuadFigure* quadFig, vtkSmartPointer< vtkPoints > pos, int interpolationLevel, int side){
    int rowNums = quadFig->getRowCount();
    int colNums = quadFig->getColumnCount();
    int step = pow((double)2, (double)interpolationLevel);

    // Get interpolated u & v coordinate, storing in interpolatedU & interpolatedV
    vector<double> interpolatedU;
    vector<double> interpolatedV;
    getUVCoordinate(interpolatedU, interpolatedV, rowNums, colNums, step);

    // Get all the interpolated boundary points (on the surface)
    M3DQuadInterpolater *tpm = new M3DQuadInterpolater(quadFig);
    Vector3D point;
    for(unsigned int i = 0; i < interpolatedU.size(); i++){
        // Get boundary points
        point = tpm->interpolateQuadSpoke(quadFig,interpolatedU[i],interpolatedV[i], side)->getB();
        pos->InsertNextPoint(point.getX(),point.getY(),point.getZ());
    }
    delete tpm;
}

void visualizesrepsurface::getCrestSurfacePoints(M3DQuadFigure* quadFig, vtkSmartPointer< vtkPoints > pos, int interpolationLevel){

    // New a crest interplator....
    vtkSmartPointer<vtkSRepInterpolateCrestSpokesQuartic> interpolatecrestspokes = vtkSmartPointer<vtkSRepInterpolateCrestSpokesQuartic>::New();
    visualizecrest visCrest;
    vtkSmartPointer<vtkSRep> srepfig = visCrest.getSrepFig(quadFig);
    interpolatecrestspokes->SetInput(srepfig);
    interpolatecrestspokes->SetInterpolationLevel(interpolationLevel);
    interpolatecrestspokes->Update();

    int rowNums = quadFig->getRowCount();
    int colNums = quadFig->getColumnCount();
    int crestQuadNum = rowNums*2 + (colNums-2)*2;
    int step = pow((double)2, (double)interpolationLevel);

    vector<double> subpoints; //store the v coordinate of all the sub-point along spoke tip connection curve.

    // Along the curve between two correspondence up and down spokes.(v direction in crest interpolate method)
    tls.splitLine(0,1,step, subpoints); // Up spoke tip is consider as 0, down spoke 1. The medial crest is 0.5.

    // Loop each crest quad.
    for(unsigned int i = 0; i< crestQuadNum; i++){
        for(unsigned int m = 0; m < subpoints.size()-1; m++){ //loop each sub points in v direction
            // Get the interpolated skeletal position at this point
            vtkSRep::VNLType s1 = interpolatecrestspokes->GetInterpolatedPoint(i,subpoints[m]);

            for(unsigned int n = 1; n < subpoints.size()-1; n++){ //loop each sub points in t direction, discard the up and down tip pos.
                // Get each interpolated boundary spoke direction at this point.
                vtkSRep::VNLType p0 = interpolatecrestspokes->GetInterpolatedSpoke(i,subpoints[m],subpoints[n]);

                pos->InsertNextPoint(s1[0] + p0[0], s1[1] + p0[1], s1[2] + p0[2]);
            }
        }
    }
}

/* Return all the sub quads's u & v coordinate into vector interpolateU & interpolateV. The shared edge repeat!!! Only can used for display.
 * interpolatedU: store the u coordinate of all the sub-quad of quad[q].
 * interpolatedV: store the v coordinate of all the sub-quad of quad[q].
*/
void visualizesrepsurface::getInterpolateUVCoordinate(VectorQuadPoint &interpolatedU, VectorQuadPoint &interpolatedV, int rowNum,
                                                      int colNum, int step){

    int quadIndex = 0;

    vector<double> subquadpoint_v, subquadpoint_u;

    for(unsigned int i = 0; i < rowNum -1; i++){ //the row number of the quads. its 3.
        for(unsigned int j = 0; j < colNum -1; j++){//coloums, its 13.
                interpolatedU.push_back(VectorDoublePoints());
                interpolatedV.push_back(VectorDoublePoints());

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
                subquadpoint_u = tls.splitQuad(quadpointsu[0],quadpointsu[1],quadpointsu[2],quadpointsu[3], step);

                //second, get the 25 subquads's v coordinate, in column first order.
                subquadpoint_v = tls.splitQuad(quadpointsv[0],quadpointsv[1],quadpointsv[2],quadpointsv[3], step);

                for(unsigned int m = 0; m < subquadpoint_u.size(); m++){
                    interpolatedU[quadIndex].push_back(subquadpoint_u[m]);
                    interpolatedV[quadIndex].push_back(subquadpoint_v[m]);
                }

                quadIndex++;

                subquadpoint_v.clear();
                subquadpoint_u.clear();
        }
    }
}


/* Return all the u v coordinate, no repeat!!
 * interpolatedU: store the u coordinate of all the sub-quad of quad[q].
 * interpolatedV: store the v coordinate of all the sub-quad of quad[q].
*/
void visualizesrepsurface::getUVCoordinate(vector<double> &interpolatedU, vector<double> &interpolatedV, int rowNum, int colNum, int step){
    vector<double> subquadpoint_v, subquadpoint_u;
    interpolatedU.clear();
    interpolatedV.clear();

    // Loop each quad, get the top-left subpositions
    for(unsigned int i = 0; i < rowNum -1; i++){ //the row number of the quads. its 3.
        for(unsigned int j = 0; j < colNum -1; j++){//coloums, its 13.
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
            subquadpoint_u = tls.splitQuad_2(quadpointsu[0],quadpointsu[1],quadpointsu[2],quadpointsu[3], step);

            // Get the 16 subquads's v coordinate, in column first order.
            subquadpoint_v = tls.splitQuad_2(quadpointsv[0],quadpointsv[1],quadpointsv[2],quadpointsv[3], step);

            for(unsigned int m = 0; m < subquadpoint_u.size(); m++){
                interpolatedU.push_back(subquadpoint_u[m]);
                interpolatedV.push_back(subquadpoint_v[m]);
            }

            subquadpoint_v.clear();
            subquadpoint_u.clear();
        }
    }

    // Get the bottom edge (rowNum -1) sub pieces.
    for(unsigned int j = 0; j < colNum -1; j++){//coloums, its 13.
        vector<double> subpieces;
        tls.splitLine(j, j+1, step, subpieces);

        for(unsigned int m = 0; m < subpieces.size()-1; m++){
            interpolatedU.push_back(rowNum -1);
            interpolatedV.push_back(subpieces[m]);
        }
    }

    // Get the the right edge (colNum -1) sub pieces
    for(unsigned int i = 0; i < rowNum -1; i++){//coloums, its 13.
        vector<double> subpieces;
        tls.splitLine(i, i+1, step, subpieces);

        for(unsigned int m = 0; m < subpieces.size()-1; m++){
            interpolatedU.push_back(subpieces[m]);
            interpolatedV.push_back(colNum -1);
        }
    }

    // Add the top-right corner.
    interpolatedU.push_back(rowNum -1);
    interpolatedV.push_back(colNum -1);
}


/* Draw the points set of the surface. */
void visualizesrepsurface::drawPointsSetSurface2(vtkSmartPointer<vtkPoints> points, vtkSmartPointer<vtkRenderer> renderer){

    vtkSmartPointer<vtkPolyData> pointsPolydata = vtkSmartPointer<vtkPolyData>::New();

    pointsPolydata->SetPoints(points);

    vtkSmartPointer<vtkVertexGlyphFilter> vertexFilter = vtkSmartPointer<vtkVertexGlyphFilter>::New();
    #if VTK_MAJOR_VERSION <= 5
    vertexFilter->SetInputConnection(pointsPolydata->GetProducerPort());
    #else
    vertexFilter->SetInputData(pointsPolydata);
    #endif
    vertexFilter->Update();

    vtkSmartPointer<vtkPolyData> polydata = vtkSmartPointer<vtkPolyData>::New();
    polydata->ShallowCopy(vertexFilter->GetOutput());

    // Setup colors
    unsigned char red[3] = {255, 0, 0};
    unsigned char green[3] = {0, 255, 0};
    unsigned char blue[3] = {0, 0, 255};
    int pointNums = points->GetNumberOfPoints();

    vtkSmartPointer<vtkUnsignedCharArray> colors = vtkSmartPointer<vtkUnsignedCharArray>::New();
    colors->SetNumberOfComponents(3);
    colors->SetName ("Colors");
    /*for (int i = 0; i < pointNums; ++i) {
        unsigned char tempColor[3] = {(int)c[i], (int)c[i+nV],
                                      (int)c[i+2*nV]};
        colors->InsertNextTupleValue (tempColor);
      }*/

    colors->InsertNextTupleValue(red);
    colors->InsertNextTupleValue(green);
    colors->InsertNextTupleValue(blue);

    polydata->GetPointData()->SetScalars(colors);

    // Visualization
    vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    #if VTK_MAJOR_VERSION <= 5
    mapper->SetInputConnection(polydata->GetProducerPort());
    #else
    mapper->SetInputData(polydata);
    #endif

    vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
    actor->SetMapper(mapper);
    actor->GetProperty()->SetPointSize(5);

    renderer->AddActor(actor);
}


/* Draw the points set of the surface. Use z coordinate as index search in the color map. */
void visualizesrepsurface::drawPointsSetSurface_zcolor(vtkSmartPointer<vtkPoints> points, vtkSmartPointer<vtkRenderer> renderer){

    vtkSmartPointer<vtkPolyData> pointsPolydata = vtkSmartPointer<vtkPolyData>::New();

    pointsPolydata->SetPoints(points);

    vtkSmartPointer<vtkVertexGlyphFilter> vertexFilter = vtkSmartPointer<vtkVertexGlyphFilter>::New();
    #if VTK_MAJOR_VERSION <= 5
    vertexFilter->SetInputConnection(pointsPolydata->GetProducerPort());
    #else
    vertexFilter->SetInputData(pointsPolydata);
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

    for(int i = 0; i < polydata->GetNumberOfPoints(); i++) {
        double p[3];
        polydata->GetPoint(i,p);

        // Get a color correspondence to the z coordinate value from the color table
        double dcolor[3];
        colorLookupTable->GetColor(p[1], dcolor); //p[2] is the z coordinate

        unsigned char color[3];
        for(unsigned int j = 0; j < 3; j++) {
            color[j] = static_cast<unsigned char>(255.0 * dcolor[j]);
        }

        colors->InsertNextTupleValue(color);
    }

    polydata->GetPointData()->SetScalars(colors);

    // Visualization
    vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    #if VTK_MAJOR_VERSION <= 5
    mapper->SetInputConnection(polydata->GetProducerPort());
    #else
    mapper->SetInputData(polydata);
    #endif

    vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
    actor->SetMapper(mapper);
    actor->GetProperty()->SetPointSize(5);

    actor->RotateX(rX);
    actor->RotateY(rY);
    actor->RotateZ(rZ);

    renderer->AddActor(actor);
}



/* Draw the points set of the surface. Use point index as index search in the color map. */
void visualizesrepsurface::drawPointsSetSurface(const char* srepfilename, int interpolationLevel, vtkSmartPointer<vtkRenderer> renderer){

    vtkSmartPointer< vtkPoints > points = getSurfacePointsSet(srepfilename, interpolationLevel);

    vtkSmartPointer<vtkPolyData> pointsPolydata = vtkSmartPointer<vtkPolyData>::New();

    pointsPolydata->SetPoints(points);

    vtkSmartPointer<vtkVertexGlyphFilter> vertexFilter = vtkSmartPointer<vtkVertexGlyphFilter>::New();
    #if VTK_MAJOR_VERSION <= 5
    vertexFilter->SetInputConnection(pointsPolydata->GetProducerPort());
    #else
    vertexFilter->SetInputData(pointsPolydata);
    #endif
    vertexFilter->Update();

    vtkSmartPointer<vtkPolyData> polydata = vtkSmartPointer<vtkPolyData>::New();
    polydata->ShallowCopy(vertexFilter->GetOutput());

    // Find min and max z
    double minz = 0;
    double maxz = polydata->GetNumberOfPoints();

    // Create the color map
    vtkSmartPointer<vtkLookupTable> colorLookupTable = vtkSmartPointer<vtkLookupTable>::New();
    colorLookupTable->SetTableRange(minz, maxz);
    colorLookupTable->Build();

    // Generate the colors for each point based on the color map
    vtkSmartPointer<vtkUnsignedCharArray> colors = vtkSmartPointer<vtkUnsignedCharArray>::New();
    colors->SetNumberOfComponents(3);
    colors->SetName("Colors");

    for(int i = 0; i < polydata->GetNumberOfPoints(); i++) {

        // Get a color correspondence to the z coordinate value from the color table
        double dcolor[3];
        colorLookupTable->GetColor(i, dcolor); //p[2] is the z coordinate

        unsigned char color[3];
        for(unsigned int j = 0; j < 3; j++) {
            color[j] = static_cast<unsigned char>(255.0 * dcolor[j]);
        }

        colors->InsertNextTupleValue(color);
    }

    polydata->GetPointData()->SetScalars(colors);

    // Visualization
    vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    #if VTK_MAJOR_VERSION <= 5
    mapper->SetInputConnection(polydata->GetProducerPort());
    #else
    mapper->SetInputData(polydata);
    #endif

    vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
    actor->SetMapper(mapper);
    actor->GetProperty()->SetPointSize(5);

    renderer->AddActor(actor);
}



/* Draw correspondence spokes. Use point index as index search in the color map. */
void visualizesrepsurface::drawCorrespondenceSpokes(const char* srepfilename, int interpolationLevel, vtkSmartPointer<vtkRenderer> renderer){

    vtkSmartPointer< vtkPoints > hubpos_medialsheet = vtkSmartPointer< vtkPoints >::New();

    // Read in this s-rep
    M3DQuadFigure* quadFig = tls.GetQuadFigure(srepfilename);

    int rowNums = quadFig->getRowCount();
    int colNums = quadFig->getColumnCount();
    int quadNum = (rowNums-1)*(colNums-1);//quadNum means the quad numbers on skeletal sheet.
    int step = pow((double)2, (double)interpolationLevel);
    int subQuadPointsNum = (step+1)*(step+1);

    // Get interpolated u & v coordinate, storing in interpolatedU & interpolatedV
    VectorQuadPoint interpolatedU;
    VectorQuadPoint interpolatedV;
    getInterpolateUVCoordinate(interpolatedU, interpolatedV, rowNums, colNums, step);

    // Get all the interpolated boundary points
    M3DQuadInterpolater *tpm = new M3DQuadInterpolater(quadFig);
    Vector3D point_s, point_b;
    // Get up spokes
    vtkSmartPointer<vtkCellArray> cellarraypointsline_up = vtkSmartPointer<vtkCellArray>::New();
    vtkSmartPointer< vtkPoints > hubpos_up = vtkSmartPointer< vtkPoints >::New();
    for(unsigned int i = 0; i < quadNum; i++){
        for(unsigned int j = 0; j < subQuadPointsNum; j++){
            point_b = tpm->interpolateQuadSpoke(quadFig,interpolatedU[i][j],interpolatedV[i][j],0)->getB();
            point_s = tpm->interpolateQuadSpoke(quadFig,interpolatedU[i][j],interpolatedV[i][j],0)->getX();

            hubpos_medialsheet->InsertNextPoint(point_s.getX(),point_s.getY(),point_s.getZ());

            vtkIdType id0 = hubpos_up->InsertNextPoint(point_s.getX(),point_s.getY(),point_s.getZ());
            vtkIdType id1 = hubpos_up->InsertNextPoint(point_b.getX(),point_b.getY(),point_b.getZ());

            vtkSmartPointer<vtkLine> medialsheetline = vtkSmartPointer<vtkLine>::New();

            medialsheetline->GetPointIds()->SetId(0, id0);
            medialsheetline->GetPointIds()->SetId(1, id1);

            cellarraypointsline_up->InsertNextCell(medialsheetline);
        }
    }

    addSpokesToRender(hubpos_up, cellarraypointsline_up, renderer, 0);

    // Get down spokes
    vtkSmartPointer<vtkCellArray> cellarraypointsline_down = vtkSmartPointer<vtkCellArray>::New();
    vtkSmartPointer< vtkPoints > hubpos_down = vtkSmartPointer< vtkPoints >::New();
    for(unsigned int i = 0; i < quadNum; i++){
        for(unsigned int j = 0; j < subQuadPointsNum; j++){
            point_b = tpm->interpolateQuadSpoke(quadFig,interpolatedU[i][j],interpolatedV[i][j],1)->getB();
            point_s = tpm->interpolateQuadSpoke(quadFig,interpolatedU[i][j],interpolatedV[i][j],1)->getX();

            vtkIdType id0 = hubpos_down->InsertNextPoint(point_s.getX(),point_s.getY(),point_s.getZ());
            vtkIdType id1 = hubpos_down->InsertNextPoint(point_b.getX(),point_b.getY(),point_b.getZ());

            vtkSmartPointer<vtkLine> medialsheetline = vtkSmartPointer<vtkLine>::New();

            medialsheetline->GetPointIds()->SetId(0, id0);
            medialsheetline->GetPointIds()->SetId(1, id1);

            cellarraypointsline_down->InsertNextCell(medialsheetline);
        }
    }

    addSpokesToRender(hubpos_down, cellarraypointsline_down, renderer, 1);

    // Draw the medial sheet.
    drawQuads(hubpos_medialsheet, quadNum, subQuadPointsNum, step, renderer);

    // Draw crest spokes
    double quadColor[3] = {0.3,0.3,0.3};//gray //{0,0.8,0};//Green.
    visualizecrest visualObject_crest(quadFig, interpolationLevel, quadColor, renderer,this->rX, this->rY, this->rZ,this->opacity);
    visualObject_crest.drawCorrespondenceCrestSpokes();

    delete tpm;
}



void visualizesrepsurface::addSpokesToRender(vtkSmartPointer< vtkPoints > hubpos, vtkSmartPointer<vtkCellArray> cellarraypointsline,
                                             vtkSmartPointer<vtkRenderer> renderer, int time){
    vtkSmartPointer<vtkPolyData> linesPolydata = vtkSmartPointer<vtkPolyData>::New();

    linesPolydata->SetPoints(hubpos);
    linesPolydata->SetLines(cellarraypointsline);

    // Find min and max z
    double minz = 0;
    double maxz = 3800; //linesPolydata->GetNumberOfPoints();
    //cout<<"--------------linesPolydata->GetNumberOfPoints() is:"<<linesPolydata->GetNumberOfPoints()<<endl; //level=2, 1200

    // Create the color map
    vtkSmartPointer<vtkLookupTable> colorLookupTable = vtkSmartPointer<vtkLookupTable>::New();
    colorLookupTable->SetTableRange(minz, maxz);
    colorLookupTable->Build();

    // Generate the colors for each point based on the color map
    vtkSmartPointer<vtkUnsignedCharArray> colors = vtkSmartPointer<vtkUnsignedCharArray>::New();
    colors->SetNumberOfComponents(3);
    colors->SetName("Colors");

    //double temp = (maxz*1.3 - maxz)/2;

    //for(int i = temp; i < maxz-temp; i++) {
    for(int i = 0*time; i < 1200*(time+1); i++) {

        // Get a color correspondence to the z coordinate value from the color table
        double dcolor[3];
        colorLookupTable->GetColor(i, dcolor); //p[2] is the z coordinate

        unsigned char color[3];
        for(unsigned int j = 0; j < 3; j++) {
            color[j] = static_cast<unsigned char>(255.0 * dcolor[j]);
        }

        colors->InsertNextTupleValue(color);
    }

    linesPolydata->GetPointData()->SetScalars(colors);

    // Visualization
    vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    #if VTK_MAJOR_VERSION <= 5
    mapper->SetInputConnection(linesPolydata->GetProducerPort());
    #else
    mapper->SetInputData(linesPolydata);
    #endif

    vtkSmartPointer<vtkActor> lineactor = vtkSmartPointer<vtkActor>::New();
    lineactor->SetMapper(mapper);
    lineactor->GetProperty()->SetLineWidth(2);
    //lineactor->RotateY(95);
    //lineactor->RotateX(30);
    //lineactor->RotateY(-35);(35);
    renderer->AddActor(lineactor);
}



void visualizesrepsurface::drawQuads(vtkSmartPointer< vtkPoints > hubpos, int quadNum, int subQuadPointsNum, int step,
                                     vtkSmartPointer<vtkRenderer> renderer){

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

    vtkSmartPointer<vtkActor> quadactor = vtkSmartPointer<vtkActor>::New();
    quadactor->SetMapper(quadMapper);

    quadactor->GetProperty()->SetColor(0.9,0.9,0.9);//1,1,0.1 yellow
    //quadactor->RotateY(95);
    //quadactor->RotateX(30);
    //quadActor->RotateY(35);
    renderer->AddActor(quadactor);
}





/* Draw correspondence spokes. Use point index as index search in the color map. */
void visualizesrepsurface::drawSpokes(const char* srepfilename, int interpolationLevel, vtkSmartPointer<vtkRenderer> renderer){

    vtkSmartPointer< vtkPoints > hubpos_medialsheet = vtkSmartPointer< vtkPoints >::New();

    // Read in this s-rep
    M3DQuadFigure* quadFig = tls.GetQuadFigure(srepfilename);

    int rowNums = quadFig->getRowCount();
    int colNums = quadFig->getColumnCount();
    int quadNum = (rowNums-1)*(colNums-1);//quadNum means the quad numbers on skeletal sheet.
    int step = pow((double)2, (double)interpolationLevel);
    int subQuadPointsNum = (step+1)*(step+1);

    // Get interpolated u & v coordinate, storing in interpolatedU & interpolatedV
    VectorQuadPoint interpolatedU;
    VectorQuadPoint interpolatedV;
    getInterpolateUVCoordinate(interpolatedU, interpolatedV, rowNums, colNums, step);

    // Get all the interpolated boundary points
    M3DQuadInterpolater *tpm = new M3DQuadInterpolater(quadFig);
    Vector3D point_s, point_b;
    // Get up spokes
    vtkSmartPointer<vtkCellArray> cellarraypointsline_up = vtkSmartPointer<vtkCellArray>::New();
    vtkSmartPointer< vtkPoints > hubpos_up = vtkSmartPointer< vtkPoints >::New();
    for(unsigned int i = 0; i < quadNum; i++){
        for(unsigned int j = 0; j < subQuadPointsNum; j++){
            point_b = tpm->interpolateQuadSpoke(quadFig,interpolatedU[i][j],interpolatedV[i][j],0)->getB();
            point_s = tpm->interpolateQuadSpoke(quadFig,interpolatedU[i][j],interpolatedV[i][j],0)->getX();

            hubpos_medialsheet->InsertNextPoint(point_s.getX(),point_s.getY(),point_s.getZ());

            vtkIdType id0 = hubpos_up->InsertNextPoint(point_s.getX(),point_s.getY(),point_s.getZ());
            vtkIdType id1 = hubpos_up->InsertNextPoint(point_b.getX(),point_b.getY(),point_b.getZ());

            vtkSmartPointer<vtkLine> medialsheetline = vtkSmartPointer<vtkLine>::New();

            medialsheetline->GetPointIds()->SetId(0, id0);
            medialsheetline->GetPointIds()->SetId(1, id1);

            cellarraypointsline_up->InsertNextCell(medialsheetline);
        }
    }

    double quadColor[3] = {0.3,0.3,0.3};//gray //{0,0.8,0};//Green.
    addSpokesToRender_2(hubpos_up, cellarraypointsline_up, renderer, quadColor);

    // Get down spokes
    vtkSmartPointer<vtkCellArray> cellarraypointsline_down = vtkSmartPointer<vtkCellArray>::New();
    vtkSmartPointer< vtkPoints > hubpos_down = vtkSmartPointer< vtkPoints >::New();
    for(unsigned int i = 0; i < quadNum; i++){
        for(unsigned int j = 0; j < subQuadPointsNum; j++){
            point_b = tpm->interpolateQuadSpoke(quadFig,interpolatedU[i][j],interpolatedV[i][j],1)->getB();
            point_s = tpm->interpolateQuadSpoke(quadFig,interpolatedU[i][j],interpolatedV[i][j],1)->getX();

            vtkIdType id0 = hubpos_down->InsertNextPoint(point_s.getX(),point_s.getY(),point_s.getZ());
            vtkIdType id1 = hubpos_down->InsertNextPoint(point_b.getX(),point_b.getY(),point_b.getZ());

            vtkSmartPointer<vtkLine> medialsheetline = vtkSmartPointer<vtkLine>::New();

            medialsheetline->GetPointIds()->SetId(0, id0);
            medialsheetline->GetPointIds()->SetId(1, id1);

            cellarraypointsline_down->InsertNextCell(medialsheetline);
        }
    }

    quadColor[0] = 0.3;
    quadColor[1] = 0.3;
    quadColor[2] = 0.3;
    addSpokesToRender_2(hubpos_down, cellarraypointsline_down, renderer, quadColor);

    // Draw the medial sheet.
    drawQuads(hubpos_medialsheet, quadNum, subQuadPointsNum, step, renderer);

    // Draw crest spokes
    quadColor[0] = 0.3;
    quadColor[1] = 0.3;
    quadColor[2] = 0.3;
    visualizecrest visualObject_crest(quadFig, interpolationLevel, quadColor, renderer,this->rX, this->rY,this->rZ,this->opacity);
    visualObject_crest.drawCrestSpokes();

    delete tpm;
}



void visualizesrepsurface::addSpokesToRender_2(vtkSmartPointer< vtkPoints > hubpos, vtkSmartPointer<vtkCellArray> cellarraypointsline,
                                             vtkSmartPointer<vtkRenderer> renderer, double * quadColor){
    vtkSmartPointer<vtkPolyData> linesPolydata = vtkSmartPointer<vtkPolyData>::New();

    linesPolydata->SetPoints(hubpos);
    linesPolydata->SetLines(cellarraypointsline);

    // Visualization
    vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    #if VTK_MAJOR_VERSION <= 5
    mapper->SetInputConnection(linesPolydata->GetProducerPort());
    #else
    mapper->SetInputData(linesPolydata);
    #endif

    vtkSmartPointer<vtkActor> lineactor = vtkSmartPointer<vtkActor>::New();
    lineactor->SetMapper(mapper);
    lineactor->GetProperty()->SetLineWidth(2);
    lineactor->GetProperty()->SetColor(quadColor);//1,1,0.1 yellow
    //lineactor->RotateY(95);
    //lineactor->RotateX(30);
    //lineactor->RotateY(-35);(35);
    renderer->AddActor(lineactor);
}






vtkSmartPointer< vtkPoints > visualizesrepsurface::getShiftedSurfacePoints(const char* srepfilename, int interpolationLevel, vector<double> vars_up,
                                                                           vector<double> vars_down, vector<double> vars_crest, int varStand, int varCrest){

    vtkSmartPointer< vtkPoints > pts = vtkSmartPointer< vtkPoints >::New();

    M3DQuadFigure* quadFig =  this->tls.GetQuadFigure(srepfilename);

    // Get up surface points, save in pts
    getPoints(quadFig, interpolationLevel, 0, vars_up, varStand, pts);

    // Get down surface points, save in pts
    getPoints(quadFig, interpolationLevel, 1, vars_down, varStand, pts);

    // Get crest surface points, save in pts
    getPoints(quadFig, interpolationLevel, 2, vars_crest, varCrest, pts);

    return pts;
}




/* Get surface point and save to the reference passing para surfacePoints
*/
void visualizesrepsurface::getPoints(M3DQuadFigure* quadFig, int interpolationLevel, int side, vector<double> vars, int varNums,
                                              vtkSmartPointer< vtkPoints > surfacePoints){


    // Convert coeffs from vector to pointer
    double *coeff = new double [varNums]; // storing up or down spokes variables.
    for(unsigned int i = 0; i < varNums; i++){
        coeff[i] = vars[i];
    }

    int rowNums = quadFig->getRowCount();
    int colNums = quadFig->getColumnCount();
    int quadNum = (rowNums-1)*(colNums-1);//quadNum means the quad numbers on skeletal sheet.
    int step = pow((double)2, (double)interpolationLevel);
    int subQuadPointsNum = (step+1)*(step+1);

    M3DQuadFigure* tempQuadFig = dynamic_cast<M3DQuadFigure*>(quadFig->clone()); //Move all base on original srep.

    switch (side){
    case 0 : // up or down
    case 1 :
    {   slidestandardspokes mStandardSpoke;
        // Shift spokes by given vars.
        mStandardSpoke.updateSpokeInfo(tempQuadFig, coeff, side);

        // Get interpolated u & v coordinate, storing in interpolatedU & interpolatedV
        VectorQuadPoint interpolatedU;
        VectorQuadPoint interpolatedV;
        getInterpolateUVCoordinate(interpolatedU, interpolatedV, rowNums, colNums, step);

        // Get all the interpolated boundary points (on the surface)
        M3DQuadInterpolater *tpm = new M3DQuadInterpolater(tempQuadFig);
        Vector3D point;
        for(unsigned int i = 0; i < quadNum; i++){
            for(unsigned int j = 0; j < subQuadPointsNum; j++){
                // Get boundary points
                point = tpm->interpolateQuadSpoke(tempQuadFig,interpolatedU[i][j],interpolatedV[i][j],side)->getB();
                surfacePoints->InsertNextPoint(point.getX(),point.getY(),point.getZ());
            }
        }
        delete tpm;
    }
        break;
    case 2 : // crest
        slidecrestspokes mCrestSpoke;
        // Shift crest spokes
        visualizecrest vcrest;
        vtkSmartPointer<vtkSRep> srepFig = vcrest.getSrepFig(quadFig);
        mCrestSpoke.updateCrestSpokeInfo(srepFig, tempQuadFig, coeff); //Move base on original srep.

        // Get the points on new (moved) surface.
        visualizesrepsurface vss;
        vss.getCrestSurfacePoints(tempQuadFig, surfacePoints, interpolationLevel);
        break;
    }

    delete coeff;
}





