#include "thinplatesplinesrep.h"

thinplatesplinesrep::thinplatesplinesrep()
{
}


thinplatesplinesrep::MeshType::Pointer thinplatesplinesrep::srepPointsSet(const char* sourcesrepfilename){

    // Read in the target s-rep (template s-rep)
    int interpolationLevel = 2;
    alignsrep as(sourcesrepfilename, interpolationLevel);
    as.weightAreaToBoundaryPoints();

    int colNums = as.getCols();
    int rowNums = as.getRows();
    int crestAtomNums = rowNums*2 + colNums*2 - 4;
    int totalAtomNums = rowNums*colNums;
    int interiorAtomNums = totalAtomNums - crestAtomNums; // 11
    int spokeNum = 2*interiorAtomNums + 3*crestAtomNums; // should equals to 106

    vector<Vector3D> sPoints = as.getPointsOnSkeletal(); // Get skeletal points.
    vector<Vector3D> bPoints = as.getBoundaryPoints();

    // Create points
    PointType point;
    MeshType::Pointer mesh_b = MeshType::New();
    MatrixType s_matrix(3, sPoints.size()); // storing skeletal hub positions.
    s_matrix.fill(0);

    // Store each srep's boundary points into mesh_b.
    for(unsigned i = 0; i< bPoints.size(); i++){
        point[0]= bPoints[i].getX();
        point[1]= bPoints[i].getY();
        point[2]= bPoints[i].getZ();
        mesh_b->SetPoint(i, point);
    }
    // Draw out the boundary mesh points, in green.
    double color[3] = {0,1,0};
    //drawMeshPoints(renderer, mesh_b, color);

    // Store each srep's skeletal points into mesh_s.
    for(unsigned i = 0; i< sPoints.size(); i++){
        s_matrix[0][i] = sPoints[i].getX();
        s_matrix[1][i] = sPoints[i].getY();
        s_matrix[2][i] = sPoints[i].getZ();
    }
    // Draw out the skeletal mesh points, in white.
    color[0] =1;
    color[1] =1;
    color[2] =1;
    //drawMeshPoints(renderer, mesh_s, color);

    // Re-order skeletal points in the same sequence as boundary points.
    MatrixType spokeTails = correspondenceSpokeTailAndTip(s_matrix, rowNums, colNums, totalAtomNums, spokeNum);

    // Create a new points set holding skeletal (tail) points and boundary (tip) points(mesh_b).
    MeshType::Pointer mesh_bs = MeshType::New();

    // Add boundary points to mesh
    PointsContainerPointer pnt = mesh_b->GetPoints();
    PointsIterator pntIt;
    PointType pos;
    int pIndex =0; //point index
    for( pntIt = pnt->Begin(); pntIt != pnt->End(); ++pntIt ) {
        for( int dim = 0; dim < 3; dim++ ) {
            pos[dim] = pntIt.Value()[dim];
        }
        mesh_bs->SetPoint(pIndex, pos);
        pIndex++;
    }
    //cout<<"----------------The next point index is: "<<pIndex<<endl; //should be 106

    // Add tail points to mesh
    for(unsigned int m =0; m< spokeTails.columns(); m++){ // loop each point
        for( int dim = 0; dim < 3; dim++ ) {
            pos[dim] = spokeTails[dim][m];
        }
        mesh_bs->SetPoint(pIndex, pos);
        pIndex++;
    }

    // Draw out the boundary + skeletal mesh points, in yellow.
    //double color[3] = {0.9,0.8,0.7};
    //drawMeshPoints(renderer, mesh_bs, color);
    cout<<"-------------------------------mesh_bs->GetPoints() is: "<<endl;
    for( pntIt = mesh_bs->GetPoints()->Begin(); pntIt != mesh_bs->GetPoints()->End(); ++pntIt ) {
        cout<<pntIt.Value()<<endl;
    }

    return mesh_bs;
}


/* Get the tps transformation from the two input .vtk file (PDM).
 * Apply tps transformation to source s-rep, generate a s-rep for target image.
*/
int thinplatesplinesrep::tps_srep(const char * sourcefilename, const char * targetfilename, const char* sourcesrepfilename){

    // Landmarks correspondances may be associated with the SplineKernelTransforms
    // via Point Set containers. Let us define containers for the landmarks.
    // Define container for landmarks
    PointSetType::Pointer sourceLandMarks = PointSetType::New();
    PointSetType::Pointer targetLandMarks = PointSetType::New();
    PointType p1; PointType p2; // same as double p1[3];
    PointSetType::PointsContainer::Pointer sourceLandMarkContainer = sourceLandMarks->GetPoints();
    PointSetType::PointsContainer::Pointer targetLandMarkContainer = targetLandMarks->GetPoints();

    //PointSetType::PointsContainer::Iterator pntIt;
    PointIdType id_s = itk::NumericTraits< PointIdType >::Zero;
    PointIdType id_t = itk::NumericTraits< PointIdType >::Zero;

    // Read in the source points set
    vtkSmartPointer<vtkPolyData> polyData_source; // same as vtkPolyData* polyData_source;
    vtkSmartPointer<vtkPolyData> polyData_target;

    vtkSmartPointer<vtkPolyDataReader> reader_source = vtkSmartPointer<vtkPolyDataReader>::New();
    vtkSmartPointer<vtkPolyDataReader> reader_target = vtkSmartPointer<vtkPolyDataReader>::New();

    reader_source->SetFileName(sourcefilename);
    reader_source->Update();
    polyData_source = reader_source->GetOutput();
    cout<<"-----------------sourcefilename is: "<<sourcefilename<<endl;

    reader_target->SetFileName(targetfilename);
    reader_target->Update();
    polyData_target = reader_target->GetOutput();
    cout<<"-----------------targetfilename is: "<<targetfilename<<endl;

    // Read in the source points set
    for(unsigned int i = 0; i < polyData_source->GetNumberOfPoints(); i++){
        double p[3];
        polyData_source->GetPoint(i,p);
        // This is identical to:
        // polydata->GetPoints()->GetPoint(i,p);
        //std::cout << "Point " << i << " : (" << p[0] << " " << p[1] << " " << p[2] << ")" << std::endl;
        p1[0] = p[0];
        p1[1] = p[1];
        p1[2] = p[2];
        std::cout << "Point " << i << " : (" << p1[0] << " " << p1[1] << " " << p1[2] << ")" << std::endl;
        sourceLandMarkContainer->InsertElement(id_s++, p1);
    }

    // Read in the target points set
    for(unsigned int i = 0; i < polyData_target->GetNumberOfPoints(); i++){
        double p[3];
        polyData_target->GetPoint(i,p);
        //std::cout << "Point " << i << " : (" << p[0] << " " << p[1] << " " << p[2] << ")" << std::endl;
        p2[0] = p[0];
        p2[1] = p[1];
        p2[2] = p[2];
        std::cout << "Point " << i << " : (" << p2[0] << " " << p2[1] << " " << p2[2] << ")" << std::endl;
        targetLandMarkContainer->InsertElement(id_t++, p2);
    }
    //sourceLandMarks->SetPoints(sourceLandMarkContainer); // no need.
    //targetLandMarks->SetPoints(targetLandMarkContainer); // no need.

    TransformType::Pointer tps = TransformType::New();
    tps->SetSourceLandmarks(sourceLandMarks);
    tps->SetTargetLandmarks(targetLandMarks);
    // Output sourceLandMarks
    /*cout<<"-----------------------------------the sourceLandMarks input to tps is: "<<endl;
    for(unsigned int i = 0; i < sourceLandMarks->GetNumberOfPoints(); i++){
        cout<<sourceLandMarks->GetPoints()->GetPoint(i)<<endl;
    }*/

    cout<<"ComputeWMatrix... "<<endl;
    tps->ComputeWMatrix();

    cout<<"ComputeWMatrix finished!"<<endl;


    /*PointType inpos;
    inpos[0] = 0.603836;
    inpos[1] = 0.579571;
    inpos[2] = 0.342368;

    PointType outpos = tps->TransformPoint(inpos);
    cout<<"========================outpos: ("<<outpos[0]<<", "<<outpos[1]<<", "<<outpos[2]<<")"<<endl;*/


    // Get this srep's points set
    MeshType::Pointer pointsSet = srepPointsSet(sourcesrepfilename);

    // to save points after transform
    PointSetType::Pointer transformed_points = PointSetType::New();
    PointSetType::PointsContainer::Pointer transformed_points_container = transformed_points->GetPoints();
    TransformType::OutputPointType trackerPointNewPosition;

    PointsContainerPointer pnt = pointsSet->GetPoints();
    PointsIterator pntIt;
    PointType pos;
    int pointIndex = 0;
    for( pntIt = pnt->Begin(); pntIt != pnt->End(); ++pntIt ) {
        for( int dim = 0; dim < 3; dim++ ) {
            pos[dim] = pntIt.Value()[dim];
        }
        trackerPointNewPosition = tps->TransformPoint(pos);
        transformed_points_container->InsertElement(pointIndex++,trackerPointNewPosition);
    }

    // Output the transformed points
    cout<<"---------------------the transformed points is:"<<endl;
    typedef PointSetType::PointsContainer::Iterator PointsIterator;
    PointsIterator pointIterator = transformed_points_container->Begin();
    PointsIterator end = transformed_points_container->End();
    while( pointIterator != end )   {
        PointType p = pointIterator.Value(); // access the point
        std::cout << p << std::endl; // print the point
        ++pointIterator; // advance to next point
    }

    // Divide the transformedPoints to transformedTips and transformedTails, and re-generate new s-rep.
    generateSrep(transformed_points, sourcesrepfilename);

    return EXIT_SUCCESS;
}


/* Divide the transformed points to transformed tips and transformed tails.
 * srepfilename: is a arbitry m3d file template, used to write new srep info and generate to be a new srep.
*/
void thinplatesplinesrep::generateSrep(PointSetType::Pointer pointSet, const char* srepfilename){
    int spokeNum = pointSet->GetNumberOfPoints()/2;
    cout<<"----------------------spokeNum is: "<<spokeNum<<endl;

    MatrixType tipMatrix(3,spokeNum);
    tipMatrix.fill(0);
    MatrixType tailMatrix(3,spokeNum);
    tailMatrix.fill(0);

    //cout<<"------------------------ the passed in outputMesh is: "<<endl;
    typedef PointSetType::PointsContainer::Iterator PointsIterator;
    PointsIterator pointIterator = pointSet->GetPoints()->Begin();
    PointsIterator end = pointSet->GetPoints()->End();
    int tipIndex = 0;
    int tailIndex = 0;
    int pIndex = 0; // point index
    while( pointIterator != end )   {
        if(pIndex < spokeNum){ // 0-105 is the boundary points
            for( int dim = 0; dim < 3; dim++ ) {
                tipMatrix[dim][tipIndex] = pointIterator.Value()[dim];
            }
            tipIndex++;
        }
        else { // 106-211 is the tails points.
            for( int dim = 0; dim < 3; dim++ ) {
                tailMatrix[dim][tailIndex] = pointIterator.Value()[dim];
            }
            tailIndex++;
        }
        pIndex++;

        //PointType p = pointIterator.Value(); // access the point
        //std::cout << p << std::endl; // print the point
        ++pointIterator; // advance to next point
    }

    cout<<"--------------------The tipMatrix is: "<<endl;
    for(unsigned int o=0; o<tipMatrix.rows(); o++){
        for(unsigned int m=0; m<tipMatrix.columns();m++){
            cout<<tipMatrix[o][m]<<"  ";
        }
        cout<<endl;
    }

    cout<<"--------------------The tailMatrix is: "<<endl;
    for(unsigned int o=0; o<tailMatrix.rows(); o++){
        for(unsigned int m=0; m<tailMatrix.columns();m++){
            cout<<tailMatrix[o][m]<<"  ";
        }
        cout<<endl;
    }

    // Generate the new spokes using the transformed boundary and its correspondence skeletal points.
    MatrixType newSpokeDirections = getNewSpokesDirection(tipMatrix, tailMatrix);
    vector<double> newSpokeRadius = getNewSpokesRadius(tipMatrix, tailMatrix);

    // Center the transformed tail points to origin.
    int pointNum = tailMatrix.columns();
    double center[3] = {0,0,0};
    for(unsigned int m = 0; m< pointNum; m++){
        for(unsigned int dim=0; dim<3; dim++){
            center[dim] += tailMatrix[dim][m];
            //theta += tailMatrix[dim][m] * tailMatrix[dim][m];
        }
    }
    //
    for(unsigned int m = 0; m< pointNum; m++){
        for(unsigned int dim=0; dim<3; dim++){
            tailMatrix[dim][m] = tailMatrix[dim][m] - center[dim]/pointNum + 0.5; // plus 0.5, move to the center of pablo display window.
        }
    }
    cout<<"--------------------After centering, the tailMatrix is: "<<endl;
    for(unsigned int o=0; o<tailMatrix.rows(); o++){
        for(unsigned int m=0; m<tailMatrix.columns();m++){
            cout<<tailMatrix[o][m]<<"  ";
        }
        cout<<endl;
    }
    /*theta /= pointNum-1;
    for(unsigned int m = 0; m< pointNum; m++){
        for(unsigned int dim=0; dim<3; dim++){
            tailMatrix[dim][m] = tailMatrix[dim][m]/sqrt(theta);
        }
    }
    cout<<"--------------------After normalization, the tailMatrix is: "<<endl;
    for(unsigned int o=0; o<tailMatrix.rows(); o++){
        for(unsigned int m=0; m<tailMatrix.columns();m++){
            cout<<tailMatrix[o][m]<<"  ";
        }
        cout<<endl;
    }*/

    // Write each transformed srep to its own .m3d files.
    saveGeneratedSreps(srepfilename, newSpokeDirections, newSpokeRadius, tailMatrix);
}



/* Save the aligned points into its srep file(.m3d file).
 *
*/
void thinplatesplinesrep::saveGeneratedSreps(const char* filePathName, MatrixType newSpokeDirection, vector<double> newSpokeRadius,
                                             MatrixType transformedTails){
    Registry registry;
    registry.readFromFile(filePathName, false);
    toolsfunc tls;
    M3DQuadFigure* quadFig = tls.GetQuadFigure(filePathName);
    int rowNums = quadFig->getRowCount();
    int colNums = quadFig->getColumnCount();
    int totalAtomNums = rowNums*colNums;
    Vector3D spokeDir;
    double spokeRadiu;
    // For interior tails
    int pindex =0;
    int crestIndex=0;
    for(unsigned int u=1; u<rowNums-1;u++){
        for(unsigned int v=1; v<colNums-1; v++){
            // up spoke direction and radius
            spokeDir.setX(newSpokeDirection[0][pindex]);
            spokeDir.setY(newSpokeDirection[1][pindex]);
            spokeDir.setZ(newSpokeDirection[2][pindex]);
            dynamic_cast<M3DQuadPrimitive*>(quadFig->getPrimitivePtr(u,v))->setU0(spokeDir);

            spokeRadiu = newSpokeRadius[pindex];
            dynamic_cast<M3DQuadPrimitive*>(quadFig->getPrimitivePtr(u,v))->setR0(spokeRadiu);

            // down spoke direction and radius
            spokeDir.setX(newSpokeDirection[0][pindex+totalAtomNums]);
            spokeDir.setY(newSpokeDirection[1][pindex+totalAtomNums]);
            spokeDir.setZ(newSpokeDirection[2][pindex+totalAtomNums]);
            dynamic_cast<M3DQuadPrimitive*>(quadFig->getPrimitivePtr(u,v))->setU1(spokeDir);

            spokeRadiu = newSpokeRadius[pindex+totalAtomNums];
            dynamic_cast<M3DQuadPrimitive*>(quadFig->getPrimitivePtr(u,v))->setR1(spokeRadiu);

            // this atom's hub positions
            dynamic_cast<M3DQuadPrimitive*>(quadFig->getPrimitivePtr(u,v))->setX(transformedTails[0][pindex],transformedTails[1][pindex],transformedTails[2][pindex]);

            pindex++; // advance to next point
        }
    }

    // For up_crest tails
    // Get the crest curve along the crest in counter clockwise. the start point is (0,0).
    // u=0, v=0,1,...,colNums
    for(unsigned int v =0; v<colNums; v++){
        setCrestSpoke(0,v, newSpokeDirection, newSpokeRadius, pindex, crestIndex, quadFig, transformedTails);
        pindex++; // advance to next point
        crestIndex++;
    }
    // u=1,2,...,rowNums, v=colNums - 1
    for(unsigned int u =1; u<rowNums; u++){
        setCrestSpoke(u, colNums-1, newSpokeDirection, newSpokeRadius, pindex, crestIndex, quadFig, transformedTails);
        pindex++;
        crestIndex++;
    }
    // u=rowNums-1, v=1,2,...,colNums-2
    for(unsigned int v =colNums-2; v>0; v--){
        setCrestSpoke(rowNums-1, v, newSpokeDirection, newSpokeRadius, pindex, crestIndex, quadFig, transformedTails);
        pindex++;
        crestIndex++;
    }
    // u=rowNums-1,...,2,1, v=0
    for(unsigned u =rowNums-1; u>0; u--){
        setCrestSpoke(u,0, newSpokeDirection, newSpokeRadius, pindex, crestIndex, quadFig, transformedTails);
        pindex++;
        crestIndex++;
    }

    // Write the changed quadFig to m3d file.
    char figureStr[1024];
    //sprintf(figureStr, "model.figure[%d]", figNum);
    sprintf(figureStr, "model.figure[%d]", 0);//In m3d file, every figure is identified as figure[0], cannot change the index, why?
    //cout<<"Message from movespokes.cpp; model.figure[%d] is: "<<figureStr<<endl;
    quadFig->writeFigure(figureStr,registry); //write the figure into registry.

    string filename = tls.getPathName(filePathName) + string("/output/")+ tls.getFileNameWithExtension(filePathName);
    cout<<"-------------Write srep to: "<< filename << endl;
    //overwrite the input file to save the changes for next loop.
    registry.writeToFile(filename.c_str());

    cout<<"-------------Finished write!" << endl;
}



/* Set new spoke direction and radiu to srep*/
void thinplatesplinesrep::setCrestSpoke(int rowIndex, int colIndex, MatrixType newSpokeDirection,
                vector<double> newSpokeRadius, int pIndex, int crestIndex, M3DQuadFigure* quadFig, MatrixType transformedTails){
    int totalAtomNums = quadFig->getRowCount() * quadFig->getColumnCount();
    Vector3D spokeDir;
    double spokeRadiu;

    // up spoke direction and radius
    spokeDir.setX(newSpokeDirection[0][pIndex]);
    spokeDir.setY(newSpokeDirection[1][pIndex]);
    spokeDir.setZ(newSpokeDirection[2][pIndex]);
    dynamic_cast<M3DQuadPrimitive*>(quadFig->getPrimitivePtr(rowIndex,colIndex))->setU0(spokeDir);

    spokeRadiu = newSpokeRadius[pIndex];
    dynamic_cast<M3DQuadPrimitive*>(quadFig->getPrimitivePtr(rowIndex,colIndex))->setR0(spokeRadiu);

    // down spoke direction and radius
    spokeDir.setX(newSpokeDirection[0][pIndex+totalAtomNums]);
    spokeDir.setY(newSpokeDirection[1][pIndex+totalAtomNums]);
    spokeDir.setZ(newSpokeDirection[2][pIndex+totalAtomNums]);
    dynamic_cast<M3DQuadPrimitive*>(quadFig->getPrimitivePtr(rowIndex,colIndex))->setU1(spokeDir);

    spokeRadiu = newSpokeRadius[pIndex+totalAtomNums];
    dynamic_cast<M3DQuadPrimitive*>(quadFig->getPrimitivePtr(rowIndex,colIndex))->setR1(spokeRadiu);

    // crest spoke direction and radius
    spokeDir.setX(newSpokeDirection[0][crestIndex+2*totalAtomNums]);
    spokeDir.setY(newSpokeDirection[1][crestIndex+2*totalAtomNums]);
    spokeDir.setZ(newSpokeDirection[2][crestIndex+2*totalAtomNums]);
    dynamic_cast<M3DQuadEndPrimitive*>(quadFig->getPrimitivePtr(rowIndex,colIndex))->setUEnd(spokeDir);

    spokeRadiu = newSpokeRadius[crestIndex+2*totalAtomNums];
    dynamic_cast<M3DQuadEndPrimitive*>(quadFig->getPrimitivePtr(rowIndex,colIndex))->setREnd(spokeRadiu);

    // this atom's hub positions
    dynamic_cast<M3DQuadPrimitive*>(quadFig->getPrimitivePtr(rowIndex,colIndex))->setX(transformedTails[0][pIndex],transformedTails[1][pIndex],transformedTails[2][pIndex]);
}



/* Correspondence the boundary-skeletal point pair.
 * Boundary points are sampled in sequence: up_middle tips, up_crest_tips, down_middle tips, down_crest tips, medial crest tips.
 * Boundary points is more than skeletal points, to compute the new spokes, we need to correspondence the boundary points to its
 * skeletal points.*/
thinplatesplinesrep::MatrixType thinplatesplinesrep::correspondenceSpokeTailAndTip(MatrixType sPoints, int rowNums, int colNums,
                                                                                           int totalAtomNums, int spokeNum){

    MatrixType tail(3, spokeNum);
    tail.fill(0);

    // For up_middle(interior) tails
    int pindex =0;
    int crestIndex=0;
    for(unsigned m=1; m<rowNums-1;m++){
        for(unsigned n=1; n<colNums-1; n++){
            for(unsigned dim=0; dim<3; dim++){
                tail[dim][pindex] = sPoints[dim][m*colNums+n];

                // down middle tails, begin from totalAtomNums to 2*totalAtomNums.
                tail[dim][pindex+totalAtomNums] = sPoints[dim][m*colNums+n];
            }
            pindex++; // advance to next point
        }
    }

    // For up_crest tails
    // Get the crest curve along the crest in counter clockwise. the start point is (0,0).
    // u=0, v=0,1,...,colNums
    for(unsigned int v =0; v<colNums; v++){
        for(unsigned dim=0; dim<3; dim++){
            tail[dim][pindex] = sPoints[dim][v];

            // down_crest tails
            tail[dim][pindex+totalAtomNums] = sPoints[dim][v];

            // For medial_crest tails
            tail[dim][crestIndex+2*totalAtomNums] = sPoints[dim][v];
        }
        pindex++;
        crestIndex++;
    }

    // u=1,2,...,rowNums, v=colNums - 1
    for(unsigned int u =1; u<rowNums; u++){
        for(unsigned dim=0; dim<3; dim++){
            tail[dim][pindex] = sPoints[dim][u * colNums + colNums - 1];

            // down_crest tails
            tail[dim][pindex+totalAtomNums] = sPoints[dim][u * colNums + colNums - 1];

            // For medial_crest tips, begin from 2*totalAtomNums to spokeNum.
            tail[dim][crestIndex+2*totalAtomNums] = sPoints[dim][u * colNums + colNums - 1];
        }
        pindex++;
        crestIndex++;
    }

    // u=rowNums-1, v=1,2,...,colNums-2
    for(unsigned int v =colNums-2; v>0; v--){
        for(unsigned dim=0; dim<3; dim++){
            tail[dim][pindex] = sPoints[dim][(rowNums-1) * colNums + v];

            // down_crest tails
            tail[dim][pindex+totalAtomNums] = sPoints[dim][(rowNums-1) * colNums + v];

            // For medial_crest tails
            tail[dim][crestIndex+2*totalAtomNums] = sPoints[dim][(rowNums-1) * colNums + v];
        }
        pindex++;
        crestIndex++;
    }

    // u=rowNums-1,...,2,1, v=0
    for(unsigned u =rowNums-1; u>0; u--){
        for(unsigned dim=0; dim<3; dim++){
            tail[dim][pindex] = sPoints[dim][u * colNums];

            // down_crest tails
            tail[dim][pindex+totalAtomNums] = sPoints[dim][u * colNums];

            // For medial_crest tails
            tail[dim][crestIndex+2*totalAtomNums] = sPoints[dim][u * colNums];
        }
        pindex++;
        crestIndex++;
    }

    // Output tail
    /*cout<<"-----Finally---------------the tail(skeletal points correspondence to the boundary points) is: "<<endl;
    for(unsigned int o=0; o<tail.rows(); o++){
        for(unsigned int m=0; m<tail.columns();m++){
            cout<<tail[o][m]<<"  ";
        }
        cout<<endl;
    }*/

    return tail;
}



/* Compute the new spokes( r and u) between correspondence boundary point and skeletal point.
 * Return a matrix list, holding all the spokes' direction for each srep.
 * bPoints: this srep's boundary points. 3*k.
 * sPoints: this srep's skeletal points. 3*k.
 * newSpokes[i] is the ith srep's spoke information. each spoke has xyz direction.
*/
thinplatesplinesrep::MatrixType thinplatesplinesrep::getNewSpokesDirection(MatrixType bPoints, MatrixType sPoints){

    // This srep's spoke direction.
    MatrixType spokeVector = bPoints - sPoints;
    cout<<"--------------------the spokeVector is: "<<endl;
    for(unsigned int o=0; o<spokeVector.rows(); o++){
        for(unsigned int m=0; m<spokeVector.columns();m++){
            cout<<spokeVector[o][m]<<"  ";
        }
        cout<<endl;
    }

    // Normalize each spoke's direction to a unit vector
    for(unsigned int m=0; m< spokeVector.columns(); m++){
        double spokeRadiu=0;
        for(unsigned int dim=0; dim<3; dim++){
            spokeRadiu += spokeVector[dim][m] * spokeVector[dim][m];
        }
        for(unsigned int dim=0; dim<3; dim++){
            spokeVector[dim][m] /= sqrt(spokeRadiu);
        }
    }

    return spokeVector;
}


/* For each pair of correspondence points on boundary and skeletal
*/
vector<double> thinplatesplinesrep::getNewSpokesRadius(MatrixType bPoints, MatrixType sPoints){

    // This srep's spoke direction.
    MatrixType spokeVector = bPoints - sPoints;

    vector<double> spokeRadius;
    // Normalize each spoke's direction to a unit vector
    for(unsigned int m=0; m< spokeVector.columns(); m++){
        double spokeRadiu=0;
        for(unsigned int dim=0; dim<3; dim++){
            spokeRadiu += spokeVector[dim][m] * spokeVector[dim][m];
        }
        spokeRadius.push_back(sqrt(spokeRadiu));
    }

    return spokeRadius;
}
