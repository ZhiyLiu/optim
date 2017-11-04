/* srepFolder: The directory where s-reps stored.
 * For example: const char* srepFolder = "/NIRAL/work/ltu/WorkSpace/data_pablo/Pablo2/lib/vtksrep/Correspondence/models/alignment/";
 * command like: ./procrustes /NIRAL/work/ltu/WorkSpace/data_pablo/Pablo2/lib/vtksrep/Correspondence/models/alignment/
 * Liyun Tu
 * Mar 25, 2014
*/



#include "procrustes.h"
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
void saveMeshMatrix(const char* filename, procrustes::MatrixType matrix){

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




/* Save the aligned points into its srep file(.m3d file).
 *
*/
void saveTransformedSreps(string filePathName, procrustes::MeshType::Pointer outputMesh){
    //cout<<"----------write....."<<filePathName<<" to /models/alignment/"<<endl;//need change.

    Registry registry;
    registry.readFromFile(filePathName.c_str(), false);

    toolsfunc tls;
    M3DQuadFigure* quadFig = tls.GetQuadFigure(filePathName.c_str());

    //Get the point container
    procrustes::MeshType::PointsContainerPointer transformedPoints = outputMesh->GetPoints();
    procrustes::PointsContainerType::ConstIterator it = transformedPoints->Begin();

    //Loop each primitive, get the x,y,z coordinate of this srep's boundary.
    for(unsigned u =0; u<quadFig->getRowCount(); u++){
        for(unsigned v =0; v<quadFig->getColumnCount(); v++){
            if( it != transformedPoints->End() ) {
              procrustes::PointType p = it.Value();
              // Set spoke hub position to new value.
              quadFig->getPrimitivePtr(u,v)->setX(p[0],p[1],p[2]);
              ++it;
            }
        }
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

    VectorTrainingSetFeaturesType sMatrix; // storing the skeletal points.

    vector<procrustes::MeshType::Pointer> X; // a list holding all the sreps.

    //a list of s-reps path file names from srepFolder.
    vector< std::string > srepPathNames = tls.getModelname(srepFolder);
    cout<<"------------ Prepare data for Procrustes alignment: "<<srepFolder<<" --------------"<<endl;

    bool reSize = false;
    vector<procrustes::MeshType::Pointer> transformedX;

    if(!srepPathNames.empty()){
        cout<<"-------------there are: "<<srepPathNames.size()<<"sreps."<<endl;
        for(int q =0; q< srepPathNames.size(); q++){

            const char* filename = srepPathNames[q].c_str();

            alignsrep as(filename, interpolationLevel);

            as.weightAreaToBoundaryPoints();

            vector<Vector3D> sPoints = as.getPointsOnSkeletal(); // Get skeletal points.

            // Make room for a pointsNums*3 row vector. Each row will store x or y or z for the correspondence spoke on q's srep.
            if(!reSize){//Only need resize one time.
                for(unsigned i=0; i<sPoints.size()*3; i++){ // need multiply 3 because each point has xyz.
                    sMatrix.push_back(VectorSRepFeaturesType());
                }
                reSize=true;
            }

            // Create points
            procrustes::PointType point;
            procrustes::MeshType::Pointer mesh_s = procrustes::MeshType::New();

            // Store each srep's skeletal points into mesh_s.
            int sIndex =0;
            for(unsigned i = 0; i< sPoints.size(); i++){
                point[0]= sPoints[i].getX();
                point[1]= sPoints[i].getY();
                point[2]= sPoints[i].getZ();
                mesh_s->SetPoint(i, point);

                sMatrix[sIndex].push_back(sPoints[i].getX());
                sMatrix[sIndex+1].push_back(sPoints[i].getY());
                sMatrix[sIndex+2].push_back(sPoints[i].getZ());
                //cout<<"------------------------sPoints["<<i<<"]-----"<<sPoints[i]<<endl;
                sIndex = sIndex+3;
            }

            //Output the points for check.
           /* std::cout << "Mesh points of skeletal is: " << mesh_s->GetNumberOfPoints() << std::endl;
            procrustes::PointsIterator pointIterator = mesh_s->GetPoints()->Begin();

            //Output these mesh points.
            procrustes::PointsIterator end = mesh_s->GetPoints()->End();
            while( pointIterator != end )   {
             procrustes::MeshType::PointType p = pointIterator.Value();  // access the point
              std::cout << p << std::endl;                    // print the point
              ++pointIterator;                                // advance to next point
            }*/


            // Store this srep's skeletal points (mesh_s) into vector list.
            X.push_back(mesh_s);
        }

        // Save sMatrix to file
        string filename = string(srepFolder) + string("GPA_input_test_sMatrix.txt");
        saveGPA_Input_Matrix(filename.c_str(), sMatrix);

        // Do procrustes alignment
        procrustes pro;
        pro.GPA(X);
        transformedX = pro.getTrans_X();
        filename = string(srepFolder) + string("ske_mesh_Matrix_0.5.txt");
        pro.saveMeshVector(filename.c_str(), transformedX);//transformedX is the procrustes transformed X all coordinate plus 0.5(pablo center).

        procrustes::MatrixType procrustesMean = pro.getProcrustesMean();//Already plus 0.5(pablo center).
        // Save to file.
        filename = string(srepFolder) + string("procrustesMean_0.5.txt");
        saveMeshMatrix(filename.c_str(), procrustesMean);



    }
    else{
        cout<<"You input a invalid srep file, please check the path or name of this srep!"<<endl;
        return EXIT_FAILURE;
    }

    // Write each transformed srep to its own .m3d files.
    if(transformedX.size() != srepPathNames.size()){
        cout<<"Something wrong while doing alignment. The transformedX.size() should equals to the input srep numbers."<<endl;
        return EXIT_FAILURE;
    }
    else{
        for(int q =0; q< srepPathNames.size(); q++){
            saveTransformedSreps(srepPathNames[q], transformedX[q]);
        }
    }

    return EXIT_SUCCESS;
}





