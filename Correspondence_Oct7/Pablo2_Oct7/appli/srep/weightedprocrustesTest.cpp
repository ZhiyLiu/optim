/* srepFolder: The directory where s-reps stored.
 * For example: const char* srepFolder = "/NIRAL/work/ltu/WorkSpace/data_pablo/Pablo2/lib/vtksrep/Correspondence/models/alignment/";
 * command like: ./weightedprocrustesTest /NIRAL/work/ltu/WorkSpace/data_pablo/Pablo2/lib/vtksrep/Correspondence/models/alignment/
 * Liyun Tu
 * Apr 1, 2014
*/



#include "weightedprocrustes.h"
#include "toolsfunc.h"



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
void setCrestSpokeUR(int rowIndex, int colIndex, MatrixType newSpokeDirection,
                vector<double> newSpokeRadius, int pIndex, int crestIndex, M3DQuadFigure* quadFig){
    int totalAtomNums = quadFig->getRowCount() * quadFig->getColumnCount();
    Vector3D spokeDir;
    double spokeRadiu;

    // up spoke
    spokeDir.setX(newSpokeDirection[0][pIndex]);
    spokeDir.setY(newSpokeDirection[1][pIndex]);
    spokeDir.setZ(newSpokeDirection[2][pIndex]);
    dynamic_cast<M3DQuadPrimitive*>(quadFig->getPrimitivePtr(rowIndex,colIndex))->setU0(spokeDir);

    spokeRadiu = newSpokeRadius[pIndex];
    dynamic_cast<M3DQuadPrimitive*>(quadFig->getPrimitivePtr(rowIndex,colIndex))->setR0(spokeRadiu);

    // down spoke
    spokeDir.setX(newSpokeDirection[0][pIndex+totalAtomNums]);
    spokeDir.setY(newSpokeDirection[1][pIndex+totalAtomNums]);
    spokeDir.setZ(newSpokeDirection[2][pIndex+totalAtomNums]);
    dynamic_cast<M3DQuadPrimitive*>(quadFig->getPrimitivePtr(rowIndex,colIndex))->setU1(spokeDir);

    spokeRadiu = newSpokeRadius[pIndex+totalAtomNums];
    dynamic_cast<M3DQuadPrimitive*>(quadFig->getPrimitivePtr(rowIndex,colIndex))->setR1(spokeRadiu);

    // crest spoke
    spokeDir.setX(newSpokeDirection[0][crestIndex+2*totalAtomNums]);
    spokeDir.setY(newSpokeDirection[1][crestIndex+2*totalAtomNums]);
    spokeDir.setZ(newSpokeDirection[2][crestIndex+2*totalAtomNums]);
    dynamic_cast<M3DQuadEndPrimitive*>(quadFig->getPrimitivePtr(rowIndex,colIndex))->setUEnd(spokeDir);

    spokeRadiu = newSpokeRadius[crestIndex+2*totalAtomNums];
    dynamic_cast<M3DQuadEndPrimitive*>(quadFig->getPrimitivePtr(rowIndex,colIndex))->setREnd(spokeRadiu);
}




/* Save the aligned points into its srep file(.m3d file).
 *
*/
void saveTransformedSreps(string filePathName, MatrixType transformedPoints,
                          MatrixType newSpokeDirection, vector<double> newSpokeRadius){

    Registry registry;
    registry.readFromFile(filePathName.c_str(), false);    
    toolsfunc tls;
    M3DQuadFigure* quadFig = tls.GetQuadFigure(filePathName.c_str());
    int rowNums = quadFig->getRowCount();
    int colNums = quadFig->getColumnCount();
    int totalAtomNums = rowNums*colNums;

    //Loop each primitive.
    int i =0; // point index
    for(unsigned u =0; u<rowNums; u++){
        for(unsigned v =0; v<colNums; v++){
            // Set spoke hub position to new value.
            quadFig->getPrimitivePtr(u,v)->setX(transformedPoints[0][i],transformedPoints[1][i],transformedPoints[2][i]);
            i++;
        }
    }

    // For spoke direction and radius.
/*    Vector3D spokeDir;
    double spokeRadiu;
    // For interior tails
    int pindex =0;
    int crestIndex=0;
    for(unsigned int u=1; u<rowNums-1;u++){
        for(unsigned int v=1; v<colNums-1; v++){
            // up spoke
            spokeDir.setX(newSpokeDirection[0][pindex]);
            spokeDir.setY(newSpokeDirection[1][pindex]);
            spokeDir.setZ(newSpokeDirection[2][pindex]);
            dynamic_cast<M3DQuadPrimitive*>(quadFig->getPrimitivePtr(u,v))->setU0(spokeDir);

            spokeRadiu = newSpokeRadius[pindex];
            dynamic_cast<M3DQuadPrimitive*>(quadFig->getPrimitivePtr(u,v))->setR0(spokeRadiu);

            // down spoke
            spokeDir.setX(newSpokeDirection[0][pindex+totalAtomNums]);
            spokeDir.setY(newSpokeDirection[1][pindex+totalAtomNums]);
            spokeDir.setZ(newSpokeDirection[2][pindex+totalAtomNums]);
            dynamic_cast<M3DQuadPrimitive*>(quadFig->getPrimitivePtr(u,v))->setU1(spokeDir);

            spokeRadiu = newSpokeRadius[pindex+totalAtomNums];
            dynamic_cast<M3DQuadPrimitive*>(quadFig->getPrimitivePtr(u,v))->setR1(spokeRadiu);

            pindex++; // advance to next point
        }
    }

    // For up_crest tails
    // Get the crest curve along the crest in counter clockwise. the start point is (0,0).
    // u=0, v=0,1,...,colNums
    for(unsigned int v =0; v<colNums; v++){
        setCrestSpokeUR(0,v, newSpokeDirection, newSpokeRadius, pindex, crestIndex, quadFig);
        pindex++; // advance to next point
        crestIndex++;
    }
    // u=1,2,...,rowNums, v=colNums - 1
    for(unsigned int u =1; u<rowNums; u++){
        setCrestSpokeUR(u, colNums-1, newSpokeDirection, newSpokeRadius, pindex, crestIndex, quadFig);
        pindex++;
        crestIndex++;
    }
    // u=rowNums-1, v=1,2,...,colNums-2
    for(unsigned int v =colNums-2; v>0; v--){
        setCrestSpokeUR(rowNums-1, v, newSpokeDirection, newSpokeRadius, pindex, crestIndex, quadFig);
        pindex++;
        crestIndex++;
    }
    // u=rowNums-1,...,2,1, v=0
    for(unsigned u =rowNums-1; u>0; u--){
        setCrestSpokeUR(u,0, newSpokeDirection, newSpokeRadius, pindex, crestIndex, quadFig);
        pindex++;
        crestIndex++;
    }
*/
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










int main( int argc, char* argv[] )
{

    if( argc != 2 ) {
        std::cerr << "Usage: "<< std::endl;
        std::cerr << argv[0];
        std::cerr << " <srepFolder>";
        std::cerr << std::endl;
        return EXIT_FAILURE;
    }

    const char * srepFolder = argv[1];
    toolsfunc tls;

    int interpolationLevel =2;
    int crestAtomNums;
    int colNums;
    int rowNums;
    int interiorAtomNums;
    int totalAtomNums;
    bool needInitialSrepInfor = true;

    vector<MeshType::Pointer> X_b; // a list holding the boundary points of all the sreps.
    VectorPoints W; // a vector holding all the area weight.
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

            vector<Vector3D> sPoints = as.getPointsOnSkeletal(); // Get skeletal points.
            vector<Vector3D> bPoints = as.getBoundaryPoints();
            // getOnesWeight() return the same points number of 1 weight. boundary and skeletal has diff points numbers.
            //vector<double> wAreas = as.getOnesWeight(106);//as.getBoundaryPointsAreas(); // Get weight.
            vector<double> wAreas = as.getOnesWeight(39);//as.getBoundaryPointsAreas(); // Get weight.

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

            // Store each srep's skeletal points into mesh_s.
            for(unsigned i = 0; i< sPoints.size(); i++){
                point[0]= sPoints[i].getX();
                point[1]= sPoints[i].getY();
                point[2]= sPoints[i].getZ();
                mesh_s->SetPoint(i, point);
            }

            //Output the points for check.
            /*std::cout << "Mesh points of boundary is: " << mesh_b->GetNumberOfPoints() << std::endl;
            PointsIterator pointIterator = mesh_b->GetPoints()->Begin();

            //Output these mesh points.
            PointsIterator end = mesh_b->GetPoints()->End();
            while( pointIterator != end )   {
             MeshType::PointType p = pointIterator.Value();  // access the point
              std::cout << p << std::endl;                    // print the point
              ++pointIterator;                                // advance to next point
            }*/

            // Store this srep's skeletal points (mesh_s) into vector list.
            X_b.push_back(mesh_b);
            W.push_back(wAreas);
            X_s.push_back(mesh_s);

            // Initialize variables for following use.
            if(needInitialSrepInfor){// only initialize once.
                crestAtomNums = as.getCrestAtomNums();
                colNums = as.getCols();
                rowNums = as.getRows();
                interiorAtomNums = as.getInteriorAtomNums();
                totalAtomNums = as.getAtomNums();

                needInitialSrepInfor = false;
            }
        }
    }
    else{
        cout<<"You input a invalid srep file, please check the path or name of this srep!"<<endl;
        return EXIT_FAILURE;
    }    

    int spokeNum = 2*interiorAtomNums + 3*crestAtomNums; // should equals to 106

    // Do procrustes alignment
    weightedprocrustes pro;

    // Get procrustes transformation from srep's boundary points.
    vector<MatrixType> transformedBoundaryPoints = pro.weightedGPA(X_b, W); //if use X_b,shoulde getOnesWeightForBoundary()

    // Apply the transformation onto its own skeletal points.
    vector<MatrixType> alignedSkeletalPoints = pro.applyTransform(X_s);

    // Test only!! input boundary point and correspondence skeletal point, the new generated m3d file is the same as original.
    /*vector<MatrixType> transformedBoundaryPoints, alignedSkeletalPoints;
    for(unsigned int j =0; j<X_s.size(); j++){//srep index
        // Convert MeshType X, Y to MatrixType
        MatrixType X(3, spokeNum);
        PointsIterator pntItX;
        int i =0;
        for( pntItX = X_s[j]->GetPoints()->Begin(); pntItX != X_s[j]->GetPoints()->End(); ++pntItX ) {
            for( int dim = 0; dim < 3; dim++ ) {
                X[dim][i] = pntItX.Value()[dim]; //X_s
            }
            i++;
        }
        alignedSkeletalPoints.push_back(X);
    }

    for(unsigned int j =0; j<X_s.size(); j++){//srep index
        // Convert MeshType X, Y to MatrixType
        MatrixType Y(3, spokeNum);
        PointsIterator pntItY;
        int i =0;
        for( pntItY = X_b[j]->GetPoints()->Begin(); pntItY != X_b[j]->GetPoints()->End(); ++pntItY ) {
            for( int dim = 0; dim < 3; dim++ ) {
                Y[dim][i] = pntItY.Value()[dim]; //X_b
            }
            i++;
        }

        transformedBoundaryPoints.push_back(Y);
        cout<<"--------------------The Boundary points is: "<<endl;
        for(unsigned int o=0; o<Y.rows(); o++){
            for(unsigned int m=0; m<Y.columns();m++){
                cout<<Y[o][m]<<"  ";
            }
            cout<<endl;
        }
    }*/

    // Correspondence the point pair.
    /* Boundary points are sampled in sequence: up_middle tips, up_crest_tips, down_middle tips, down_crest tips, medial crest tips.
     * Boundary points is more than skeletal points, to compute the new spokes, we need to correspondence the boundary points to its
     * skeletal points.*/
    vector<MatrixType> tails; // spoke tails.
    MatrixType tail(3, spokeNum);
    // For each srep, generate a new skeletal points(have repeat points) matrix according to the boundary points sequence.
    for(unsigned int i =0; i<X_s.size(); i++){//srep index
        MatrixType sPoints = alignedSkeletalPoints[i]; // this srep's tansformed skeletal points. 3*k.
       /* cout<<"--------------------the sPoints is: "<<endl;
        for(unsigned int o=0; o<sPoints.rows(); o++){
            for(unsigned int m=0; m<sPoints.columns();m++){
                cout<<sPoints[o][m]<<"  ";
            }
            cout<<endl;
        }*/

        /*MatrixType bPoints = transformedBoundaryPoints[i]; // this srep's tansformed skeletal points. 3*k.
        cout<<"--------------------the bPoints is: "<<endl;
        for(unsigned int o=0; o<bPoints.rows(); o++){
            for(unsigned int m=0; m<bPoints.columns();m++){
                cout<<bPoints[o][m]<<"  ";
            }
            cout<<endl;
        }*/

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

    // Generate the new spokes using the transformed boundary and its correspondence skeletal points.
    //vector<MatrixType> newSpokeDirections = pro.getNewSpokesDirection(transformedBoundaryPoints, tails);
    //VectorPoints newSpokeRadius = pro.getNewSpokesRadius(transformedBoundaryPoints, tails);

    // Test only. using the boundary point
    vector<MatrixType> newSpokeDirections = pro.getNewSpokesDirection(transformedBoundaryPoints, tails);
    VectorPoints newSpokeRadius = pro.getNewSpokesRadius(transformedBoundaryPoints, tails);

    // Write each transformed srep to its own .m3d files.
    if(transformedBoundaryPoints.size() != srepPathNames.size()){
        cout<<"Something wrong while doing alignment. The transformedX.size() should equals to the input srep numbers."<<endl;
        return EXIT_FAILURE;
    }
    else{
        // For each srep
        for(int q =0; q< alignedSkeletalPoints.size(); q++){
            saveTransformedSreps(srepPathNames[q], alignedSkeletalPoints[q], newSpokeDirections[q], newSpokeRadius[q]);
        }
    }

    // Compute the procrustes distance over this transformed srep smaples.
    //double procrustesDis = pro.procrustesDistance(transformedBoundaryPoints);
    //cout<<"The procrustes distance over this transformed srep skeletal points is: "<<procrustesDis<<endl;






    return EXIT_SUCCESS;
}





