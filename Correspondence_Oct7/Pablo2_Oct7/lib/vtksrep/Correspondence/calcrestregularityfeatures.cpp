#include "calcrestregularityfeatures.h"

calcrestregularityfeatures::calcrestregularityfeatures()
{
}



calcrestregularityfeatures::calcrestregularityfeatures(M3DQuadFigure* quadfig, DoubleVec subTVs, int interpolationLevel, int crestAtomNums) {
    this->subTVs = subTVs;

    this->quadfig = quadfig;
    this->interpolationLevel = interpolationLevel;

    step = pow(2.0, (double)interpolationLevel);

    this->subQuadVertexNum = pow((step+1.0), 2.0); //(step+1)*(step+1);

    this->exteriorAtomNum = crestAtomNums;
}



/* Get all the sub-quad position by given u v coordinates.
 * quadtype: boundary quads(0), skeletal quads(1).
 *
 * points_s: size is (step+1) *quadNum
 * points_b: size is (step+1)*(step+1) *quadNum
*/
void calcrestregularityfeatures::getSubQuadsPositionOfCrestRegion(vtkSmartPointer< vtkPoints > points_s, vtkSmartPointer< vtkPoints > points_b){

    vtkSmartPointer<vtkSRep> srepfig = getSrepFig(this->quadfig); // Use quadFig to compute its srepfig
    // New a crest interplator....
    vtkSmartPointer<vtkSRepInterpolateCrestSpokesQuartic> interpolatecrestspokes = vtkSmartPointer<vtkSRepInterpolateCrestSpokesQuartic>::New();
    interpolatecrestspokes->SetInput(srepfig);
    interpolatecrestspokes->SetInterpolationLevel(this->interpolationLevel);
    interpolatecrestspokes->Update();

    int linePosNum = this->subTVs.size(); // points number on each line, equals to step+1

    // for each crest quad
    for(int i = 0; i < this->exteriorAtomNum; i++){

        for(unsigned int m = 0; m < linePosNum; m++){ // In v direction (along the fold curve)

            // get the interpolated skeletal position at this point
            vtkSRep::VNLType s = interpolatecrestspokes->GetInterpolatedPoint(i, this->subTVs[m]);
            points_s->InsertNextPoint(s[0], s[1], s[2]);

            for(unsigned int n = 0; n < linePosNum; n++){ // In h direction (pendicular to the fold curve)
                // get the tips of the interpolated spokes pointing out from this point
                vtkSRep::VNLType b = interpolatecrestspokes->GetInterpolatedSpoke(i, this->subTVs[m], this->subTVs[n]);
                points_b->InsertNextPoint(s[0] + b[0], s[1] + b[1], s[2] + b[2]);
            }
        }
    }
}


void calcrestregularityfeatures::calculateCrestEdges_2(vtkSmartPointer< vtkPoints > points_s, vtkSmartPointer< vtkPoints > points_b,
                                                       MatrixType &verEdgeFeatures, MatrixType &horEdgeFeatures){
    double p_c[3]; // current point
    double p_n[3]; // next point

//    vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
//    renderer->SetBackground(0,0,0);
//    vtkSmartPointer<vtkRenderWindow> renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
//    renderWindow->AddRenderer(renderer);
//    vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
//    renderWindowInteractor->SetRenderWindow(renderWindow);

//    vtkSmartPointer<vtkPoints> hubpoints  = vtkSmartPointer<vtkPoints>::New();
//    vtkSmartPointer<vtkCellArray> cellarraypointsline = vtkSmartPointer<vtkCellArray>::New();
//    vtkSmartPointer<vtkPolyData> polypointsline = vtkSmartPointer<vtkPolyData>::New();


    // Loop each crest point.
    for(unsigned int i = 0; i< this->exteriorAtomNum; i++){

        double verticalEdge_s = 0.0; // The lines along medial sheet crest (in t direction). each crest point has one.
        double verticalEdge_b = 0.0; // The lines along crest spokes' tip. each crest point has one vertical edge.
        double horizonaledge_b =0.0; // The lines along up spoke and down spokes' tip. each crest point has one horizonal edge.

        for(unsigned int m = 0; m < step; m++){
            int currIndex = i*(step+1) + m; // i*(step+1) + m*(step+1);
            // The vertical edges on skeletal. (in v direction, along the fold curve)
            points_s->GetPoint(currIndex, p_c);
            points_s->GetPoint(currIndex+1, p_n);
            verticalEdge_s += pointsDistance(p_c, p_n);

            // The vertical edges on boundary. use the middel line (fold curve on boundary)
            int tempIndex = i*this->subQuadVertexNum + step/2 + m*(step+1);
            points_b->GetPoint(tempIndex, p_c);
            points_b->GetPoint(tempIndex+step+1, p_n);
            verticalEdge_b += pointsDistance(p_c, p_n);

             // The horizontal edges on boundary. each crest atom has one.
            currIndex = i*this->subQuadVertexNum + m;
            points_b->GetPoint(currIndex, p_c);
            points_b->GetPoint(currIndex+1, p_n);
            horizonaledge_b += pointsDistance(p_c, p_n);


//            // add sub-line to array.
//            vtkIdType id0 = hubpoints->InsertNextPoint(p_c[0], p_c[1], p_c[2]);
//            vtkIdType id1 = hubpoints->InsertNextPoint(p_n[0], p_n[1], p_n[2]);

//            vtkSmartPointer<vtkLine> medialsheetline = vtkSmartPointer<vtkLine>::New();
//            medialsheetline->GetPointIds()->SetId(0, id0);
//            medialsheetline->GetPointIds()->SetId(1, id1);

//            cellarraypointsline->InsertNextCell(medialsheetline);

//            polypointsline->SetPoints(hubpoints);
//            polypointsline->SetLines(cellarraypointsline);
        }

        // Add this vertical line on skeletal to list
        verEdgeFeatures[0][i] = verticalEdge_s;

        // Add this vertical line on boundary to list
        verEdgeFeatures[1][i] = verticalEdge_b;

        horEdgeFeatures[0][i] = horizonaledge_b;
    }

//    vtkSmartPointer<vtkPolyDataMapper> pointlinemapper = vtkSmartPointer<vtkPolyDataMapper>::New();
//    pointlinemapper->SetInput(polypointsline);
//    vtkSmartPointer<vtkActor>  lineactor = vtkActor::New();
//    lineactor->SetMapper(pointlinemapper);
//    lineactor->GetProperty()->SetLineWidth(1);
//    lineactor->GetProperty()->SetColor(0,1,0);// (0,1,0) is green, (1,1,1) is white.
//    renderer->AddActor(lineactor);

//    renderWindow->Render();
//    renderWindowInteractor->Start();

}


//calculate the distance between two points.
double calcrestregularityfeatures::pointsDistance(double * p1, double * p2){

    double distance = 0.0;

    for(unsigned int dim = 0; dim < 3; dim++){
        distance += pow((p2[dim] - p1[dim]), 2.0);
    }

    return sqrt(distance);
}




/* compute the avarage angle of the crest boundary.
 * This method do the same thing as calculateCrestAngles(), but much faster... (9920 ms VS 27710 ms)
*/
void calcrestregularityfeatures::calculateCrestAngles_2(vtkSmartPointer< vtkPoints > points_b, MatrixType &angleFeatures){

    double p0[3];
    double p1[3];
    double p2[3];
    double p3[3];

    toolsfunc tools;    

    // Loop each crest quad
    for(unsigned int n = 0; n < this->exteriorAtomNum; n++){
        double subAngle_topLeft = 0.0;
        double subAngle_botRight = 0.0;

        int index = n*this->subQuadVertexNum; // points number of each quad

        for(unsigned int m = 0; m < step; m++){ //loop each sub points along v direction
            for(unsigned int k = 0; k < step; k++){
                int currentpoint = k + m*(step+1);

                // top-left point.
                points_b->GetPoint(index+currentpoint,p0);  //m
                // bottom-left point
                points_b->GetPoint(index+currentpoint+1,p1);//m+1
                // bottom-right point
                points_b->GetPoint(index+currentpoint+1+step+1,p2);//m+6
                // top-right point
                points_b->GetPoint(index+currentpoint+step+1,p3);//m+5

                subAngle_topLeft += tools.dotProductAngle(p0, p3, p1);
                subAngle_botRight += tools.dotProductAngle(p2, p1, p3);
            }
        }

        // top-left corner trangular
        points_b->GetPoint(index,p0);
        points_b->GetPoint(index+(step+1),p3);
        points_b->GetPoint(index+1,p1);
        Vector3D normal_tl = tools.trangularNormal_doublePot(p0, p3, p1);

        // bottom-right corner trangular
        int p1_index = index + (step+1)*step - 1;
        points_b->GetPoint(p1_index+step+1,p2);
        points_b->GetPoint(p1_index,p1);
        points_b->GetPoint(p1_index+step,p3);
        Vector3D normal_br = tools.trangularNormal_doublePot(p2, p1, p3);

        normal_tl.normalize();
        normal_br.normalize();

        int subQuadNum = this->step * this->step;

        // up-left angle
        angleFeatures[0][n] = subAngle_topLeft/subQuadNum;

        // bot-right corner
        angleFeatures[1][n] = subAngle_botRight/subQuadNum;

        // normal swing
        angleFeatures[2][n] = normal_tl * normal_br;
    }
}



/* Given a quadfig, generate its srepfig.*/
vtkSmartPointer<vtkSRep> calcrestregularityfeatures::getSrepFig(M3DQuadFigure* quadfig){
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
        //const float *color = quadfig->getColor();

       //m_Srep->SetColor(color[0], color[1], color[2]);

        allspokes.clear();
        allradius.clear();
        pointsIds.clear();
    }

    return m_Srep;
}








