/* Since the weightedprocrustes method do bad than weight=1, because the curverture of the boundary been weighted, which we don't want.
 * Here we take each srep's boundary points and skeletal points together, align the point set at a same time.
 *
 *
 * Liyun Tu
 * Apr 7, 2014
*/

/* srepFolder: The directory where s-reps stored.
 * For example: const char* srepFolder = "/NIRAL/work/ltu/WorkSpace/data_pablo/Pablo2/lib/vtksrep/Correspondence/models/alignment/";
 * command like: ./allPointsProcrustesTest /NIRAL/work/ltu/WorkSpace/May15/Pablo2/lib/vtksrep/Correspondence/models/alignment/
 * Liyun Tu
 * Apr 7, 2014
*/


#include "weightedprocrustes.h"
#include "toolsfunc.h"
#include "visualization.h"

#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderWindow.h>


/* Save matrix to .txt file.*/
void saveGPA_Input_Matrix(const char* filename, VectorTrainingSetFeaturesType matrix){

    std::ofstream fout;
    fout.open(filename);

    if(fout)  {
        for(int i =0; i< matrix.size();i++){
            for(int j=0; j<matrix[0].size(); j++){
                fout<< matrix[i][j]<<" ";
            }
            fout<<endl;
        }

        cout<<"Successfully saved matrix to: "<<filename<<endl;
    }
    else
        cerr<<"Write out failed, cannot open the file!"<<endl;

    fout.close();
}



/* Save matrix to .txt file.*/
void saveMeshMatrix(const char* filename, MatrixType matrix){

    std::ofstream fout;
    fout.open(filename);

    if(fout)  {
        for(int i =0; i< matrix.rows();i++){
            for(int j=0; j<matrix.columns(); j++){
                fout<< matrix[i][j]<<" ";
            }
            fout<<endl;
        }

        cout<<"Successfully saved matrix to: "<<filename<<endl;
    }
    else
        cerr<<"Write out failed, cannot open the file!"<<endl;

    fout.close();
}




/* Set new spoke direction and radiu to srep*/
void setCrestSpoke(int rowIndex, int colIndex, MatrixType newSpokeDirection,
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




/* Save the aligned points into its srep file(.m3d file).
 *
*/
void saveTransformedSreps(string filePathName, MatrixType newSpokeDirection, vector<double> newSpokeRadius,MatrixType transformedTails){

    Registry registry;
    registry.readFromFile(filePathName.c_str(), false);
    toolsfunc tls;
    M3DQuadFigure* quadFig = tls.GetQuadFigure(filePathName.c_str());
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

    //overwrite the input file to save the changes for next loop.
    registry.writeToFile(filename.c_str());

    cout<<"-------------Finished write srep to: "<<filename<<endl;
}


void drawMeshPoints(vtkSmartPointer<vtkRenderer> renderer, MeshType::Pointer mesh_p, double*color){

    // Create the topology of the point (a vertex)
    vtkSmartPointer<vtkCellArray> vertices = vtkSmartPointer<vtkCellArray>::New();

    // Create the geometry of a point (the coordinate)
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();

    // Loop points in mesh_p
    PointsContainerPointer pnt = mesh_p->GetPoints();
    PointsIterator pntIt;
    for( pntIt = pnt->Begin(); pntIt != pnt->End(); ++pntIt ) {
        vtkIdType pid[1];
        pid[0] = points->InsertNextPoint(pntIt.Value()[0], pntIt.Value()[1], pntIt.Value()[2]);
        vertices->InsertNextCell(1,pid);
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
    pointactor->GetProperty()->SetPointSize(3);
    pointactor->GetProperty()->SetColor(color);// (0,1,0) is green, (1,1,1) is white.

    renderer->AddActor(pointactor);
}


/*
 * matrix_p: 3*n matrix, column index is points number, row index is x y z.
*/
void drawMatrixPoints(vtkSmartPointer<vtkRenderer> renderer, MatrixType matrix_p, double*color){

    // Create the topology of the point (a vertex)
    vtkSmartPointer<vtkCellArray> vertices = vtkSmartPointer<vtkCellArray>::New();

    // Create the geometry of a point (the coordinate)
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();

    // Loop points in matrix_p
    for(unsigned int p = 0; p < matrix_p.columns();p++){ // loop each points
        for(unsigned int dim = 0; dim < 3; dim++ ) {
            vtkIdType pid[1];
            pid[0] = points->InsertNextPoint(matrix_p[0][p], matrix_p[1][p], matrix_p[2][p]);
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
    pointactor->GetProperty()->SetPointSize(3);
    pointactor->GetProperty()->SetColor(color);// (0,1,0) is green, (1,1,1) is white.

    renderer->AddActor(pointactor);
}




int main( int argc, char* argv[] )
{

    if( argc != 2 ) {
        std::cerr << "Usage: "<< std::endl;
        std::cerr << argv[0];
        std::cerr << " <srepFolder>";
        std::cerr << std::endl;
        return EXIT_FAILURE;
    }

    vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
    renderer->SetBackground(0,0,0);
    vtkSmartPointer<vtkRenderWindow> renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
    renderWindow->AddRenderer(renderer);
    vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
    renderWindowInteractor->SetRenderWindow(renderWindow);

    const char * srepFolder = argv[1];
    toolsfunc tls;

    int interpolationLevel =2;
    int crestAtomNums;
    int colNums;
    int rowNums;
    int interiorAtomNums;
    int totalAtomNums;
    bool needInitialSrepInfor = true;
    int spokeNum;

    vector<MeshType::Pointer> X_b; // a list holding the boundary points of all the sreps.
    VectorPoints W; // a vector holding all the area weight correspondence to each points.
    vector<MeshType::Pointer> X_s; // a list holding the skeletal points of all the sreps.

    //a list of s-reps path file names from srepFolder.
    vector< std::string > srepPathNames = tls.getModelname(srepFolder);
    cout<<"------------ Prepare data for Procrustes alignment: "<<srepFolder<<" --------------"<<endl;

    if(!srepPathNames.empty()){
        cout<<"-------------there are: "<<srepPathNames.size()<<"sreps."<<endl;
        for(int q =0; q< srepPathNames.size(); q++){

            const char* filename = srepPathNames[q].c_str();

            alignsrep as(filename, interpolationLevel);

            as.weightAreaToBoundaryPoints();

            // Initialize variables for following use.
            if(needInitialSrepInfor){// only initialize once.
                crestAtomNums = as.getCrestAtomNums();
                colNums = as.getCols();
                rowNums = as.getRows();
                interiorAtomNums = as.getInteriorAtomNums();
                totalAtomNums = as.getAtomNums();
                spokeNum = 2*interiorAtomNums + 3*crestAtomNums; // should equals to 106

                needInitialSrepInfor = false;
            }

            vector<Vector3D> sPoints = as.getPointsOnSkeletal(); // Get skeletal points.
            vector<Vector3D> bPoints = as.getBoundaryPoints();
            vector<double> wAreas = as.getOnesWeight(spokeNum*2); // Get weight. Currently, all points set to 1.

            // Create points
            PointType point;
            MeshType::Pointer mesh_b = MeshType::New();
            MeshType::Pointer mesh_s = MeshType::New();

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
                point[0]= sPoints[i].getX();
                point[1]= sPoints[i].getY();
                point[2]= sPoints[i].getZ();
                mesh_s->SetPoint(i, point);
            }
            // Draw out the skeletal mesh points, in white.
            color[0] =1;
            color[1] =1;
            color[2] =1;
            //drawMeshPoints(renderer, mesh_s, color);

            // Store this srep's skeletal points (mesh_s) into vector list.
            X_b.push_back(mesh_b);
            W.push_back(wAreas);
            X_s.push_back(mesh_s);
        }
    }
    else{
        cout<<"You input a invalid srep file, please check the path or name of this srep!"<<endl;
        return EXIT_FAILURE;
    }

    // Convert X_s from mesh to matrix
    vector<MatrixType> X_s_matrix;
    for(unsigned int i =0; i<X_s.size(); i++){//srep index
        // Add boundary points to mesh
        PointsContainerPointer pnt = X_s[i]->GetPoints();
        PointsIterator pntIt;
        MatrixType s_matrix(3, spokeNum);  // ------------------------------this shouldnt be spokenum, 39?????????????????????????
        int pIndex =0; //point index
        for( pntIt = pnt->Begin(); pntIt != pnt->End(); ++pntIt ) {
            for( int dim = 0; dim < 3; dim++ ) {
                s_matrix[dim][pIndex] = pntIt.Value()[dim];
            }
            pIndex++;
        }
        X_s_matrix.push_back(s_matrix);
    }

    // Correspondence the boundary-skeletal point pair.
    /* Boundary points are sampled in sequence: up_middle tips, up_crest_tips, down_middle tips, down_crest tips, medial crest tips.
     * Boundary points is more than skeletal points, to compute the new spokes, we need to correspondence the boundary points to its
     * skeletal points.*/
    vector<MatrixType> tails; // spoke tails.
    MatrixType tail(3, spokeNum);
    // For each srep, generate a new skeletal points(have repeat points) matrix according to the boundary points sequence.
    for(unsigned int i =0; i<X_s_matrix.size(); i++){//srep index
        MatrixType sPoints = X_s_matrix[i]; // this srep's tansformed skeletal points. 3*k.
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
        tails.push_back(tail);

        // Output tail
        /*cout<<"-----Finally---------------the tail(skeletal points correspondence to the boundary points) is: "<<endl;
        for(unsigned int o=0; o<tail.rows(); o++){
            for(unsigned int m=0; m<tail.columns();m++){
                cout<<tail[o][m]<<"  ";
            }
            cout<<endl;
        }*/
    }

    // Create a new points set for procrustes by adding tails points to boundary points(X_b).
    vector<MeshType::Pointer> X_bAnds; // a list holding the boundary points and tails points of all the sreps.
    // For each srep.
    for(unsigned int i =0; i< tails.size(); i++){//srep index
        MeshType::Pointer mesh_bs = MeshType::New();

        // Add boundary points to mesh
        PointsContainerPointer pnt = X_b[i]->GetPoints();
        PointsIterator pntIt;
        PointType point;
        int pIndex =0; //point index
        for( pntIt = pnt->Begin(); pntIt != pnt->End(); ++pntIt ) {
            for( int dim = 0; dim < 3; dim++ ) {
                point[dim] = pntIt.Value()[dim];
            }
            mesh_bs->SetPoint(pIndex, point);
            pIndex++;
        }
        //cout<<"----------------The next point index is: "<<pIndex<<endl; //should be 106

        // Add tail points to mesh
        MatrixType tailPoints = tails[i]; // this srep's tail points. 3*k.
        for(unsigned int m =0; m< tailPoints.columns(); m++){ // loop each point
            for( int dim = 0; dim < 3; dim++ ) {
                point[dim] = tailPoints[dim][m];
            }
            mesh_bs->SetPoint(pIndex, point);
            pIndex++;
        }

        // Save this mesh points
        X_bAnds.push_back(mesh_bs);

        // Draw out the boundary + skeletal mesh points, in yellow.
        double color[3] = {0.9,0.8,0.7};
        //drawMeshPoints(renderer, mesh_bs, color);
    }


    // Do procrustes alignment
    weightedprocrustes pro(false); // Set false to use all 1 weight.
    // Get procrustes transformation from srep's boundary + skeletal points.
    vector<MatrixType> transformedPoints = pro.weightedGPA(X_bAnds, W);

    // ----Test only
    //pro.weightedGPA(X_bAnds, W);
    //vector<MatrixType> transformedPoints = pro.getXp(); // Get the centered points set.
    //vector<MatrixType> transformedPoints = pro.getQ(); // Get the centered + scaled points set.
    // ----End of test only.

    // Divide the transformedPoints to transformedTips and transformedTails
    vector<MatrixType> transformedTips, transformedTails;
    for(unsigned int i =0; i< transformedPoints.size(); i++){//srep index
        MatrixType trans_Points = transformedPoints[i]; // each transformedPoints

        MatrixType tipMatrix(3,spokeNum);
        // 0-105 is the boundary points
        for(unsigned int p=0; p<spokeNum;p++){ // loop each points
            for( int dim = 0; dim < 3; dim++ ) {
                tipMatrix[dim][p] = trans_Points[dim][p];
            }
        }
        /*cout<<"--------------------The tipMatrix is: "<<endl;
        for(unsigned int o=0; o<tipMatrix.rows(); o++){
            for(unsigned int m=0; m<tipMatrix.columns();m++){
                cout<<tipMatrix[o][m]<<"  ";
            }
            cout<<endl;
        }*/

        transformedTips.push_back(tipMatrix);
        // Draw out the aligned boundary points, in green.
        double color[3] = {0,1,0};
        drawMatrixPoints(renderer, tipMatrix, color);

        MatrixType tailMatrix(3,spokeNum);
        // 106-211 is the tails points.
        int posIndex =0; //point index
        for(unsigned int p=spokeNum; p<trans_Points.columns();p++){ // loop each points
            for( int dim = 0; dim < 3; dim++ ) {
                tailMatrix[dim][posIndex] = trans_Points[dim][p];
            }
            posIndex++;
        }
        /*cout<<"--------------------The tailMatrix is: "<<endl;
        for(unsigned int o=0; o<tailMatrix.rows(); o++){
            for(unsigned int m=0; m<tailMatrix.columns();m++){
                cout<<tailMatrix[o][m]<<"  ";
            }
            cout<<endl;
        }*/
        transformedTails.push_back(tailMatrix);
        // Draw out the aligned skeletal points, in white.
        color[0] =1;
        color[1] =1;
        color[2] =1;
        drawMatrixPoints(renderer, tailMatrix, color);
    }

    // transformedTips Equals to pro.applyTransform(X_b);
    // transformedTails Equals to pro.applyTransform(tails);
    /*// transformedTips Equals to transformedTips_1;
    vector<MatrixType> transformedTips_1 = pro.applyTransform(X_b);
    // transformedTails Equals to transformedTails_1;
    //vector<MatrixType> transformedTails_1 = pro.applyTransform(tails);

    // Output transformedTips_1 and transformedTails_1
    for(unsigned int i =0; i< transformedTips_1.size(); i++){//srep index
        // output transformedTips_1
        cout<<"--------------the ttip is: "<<endl;
        MatrixType ttip = transformedTips_1[i];
        for(unsigned int j =0; j<ttip.columns();j++){
            for(int dim=0; dim<3; dim++){
                cout<<ttip[dim][j]<<"  ";
            }
            cout<<endl;
        }
        // output transformedTails_1
        cout<<"--------------the ttail is: "<<endl;
        MatrixType ttail = transformedTails_1[i];
        for(unsigned int j =0; j<ttail.columns();j++){
            for(int dim=0; dim<3; dim++){
                cout<<ttail[dim][j]<<"  ";
            }
            cout<<endl;
        }
    }*/



    // Generate the new spokes using the transformed boundary and its correspondence skeletal points.
    vector<MatrixType> newSpokeDirections = pro.getNewSpokesDirection(transformedTips, transformedTails);
    VectorPoints newSpokeRadius = pro.getNewSpokesRadius(transformedTips, transformedTails);

    // Apply the gotten transformation onto points set (skeletal points as example here).
    // value in alignedSkeletalPoints and transformedTails is same value, means the apply is correctly.
    //vector<MatrixType> alignedSkeletalPoints = pro.applyTransform(X_s);

    // Write each transformed srep to its own .m3d files.
    if(transformedTips.size() != srepPathNames.size()){
        cout<<"Something wrong while doing alignment. The transformedX.size() should equals to the input srep numbers."<<endl;
        return EXIT_FAILURE;
    }
    else{
        // For each srep
        for(int q =0; q< transformedTails.size(); q++){
            saveTransformedSreps(srepPathNames[q], newSpokeDirections[q], newSpokeRadius[q], transformedTails[q]);
        }
    }

    // Compute the procrustes distance over this transformed srep smaples.
    //double procrustesDis = pro.procrustesDistance(transformedBoundaryPoints);
    //cout<<"The procrustes distance over this transformed srep skeletal points is: "<<procrustesDis<<endl;



    renderWindow->Render();
    renderWindowInteractor->Start();


    return EXIT_SUCCESS;
}






