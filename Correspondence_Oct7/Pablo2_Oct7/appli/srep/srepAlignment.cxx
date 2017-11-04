/* srepFolder: The directory where s-reps stored.
 * For example: const char* srepFolder = "/NIRAL/work/ltu/WorkSpace/data_pablo/Pablo2/lib/vtksrep/Correspondence/models/alignment/";
 * command like: ./srepAlignment /NIRAL/work/ltu/WorkSpace/data_pablo/Pablo2/lib/vtksrep/Correspondence/models/alignment/
*/




#include "alignsrep.h"



/* Save matrix to .txt file.*/
void saveGPA_Input_Matrix(const char* filename, VectorTrainingSetFeaturesType pMatrix){

    std::ofstream fout;
    fout.open(filename);
    //cout<<"---------------pMatrix.size() is: "<<pMatrix.size()<<endl<<"---------------pMatrix[0].size() is: "<<pMatrix[0].size()<<endl;

    if(fout)  {
        for(int i =0; i< pMatrix.size();i++){
            for(int j=0; j<pMatrix[0].size(); j++){
                fout<< pMatrix[i][j]<<" ";
            }
            fout<<endl;
        }

        cout<<"Successfully saved p values to: "<<filename<<endl;
    }
    else
        cerr<<"Write out failed, cannot open the file!"<<endl;

    fout.close();
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

    int side =0;
    int interpolationLevel =2;
    int quadtype = 1;
    bool moved = true;
    double quadColor[3] = {1,0.5,0};
    int quadNums;


    VectorTrainingSetFeaturesType pMatrix;
    VectorTrainingSetFeaturesType wMatrix; // storing the area weight of each point.
    VectorTrainingSetFeaturesType sMatrix; // storing the skeletal points.

    //a list of s-reps path file names from srepFolder.
    vector< std::string > srepPathNames = tls.getModelname(srepFolder);
    bool reSize = false;//used to make space for the first time.

    cout<<"------------ Prepare data for Procrustes alignment: "<<srepFolder<<" --------------"<<endl;
    if(!srepPathNames.empty()){
        cout<<"-------------there are: "<<srepPathNames.size()<<"sreps."<<endl;
        for(int q =0; q< srepPathNames.size(); q++){

            const char* filename = srepPathNames[q].c_str();

            alignsrep as(filename, interpolationLevel);

            as.weightAreaToBoundaryPoints();

            vector<Vector3D> bPoints = as.getBoundaryPoints();
            vector<Vector3D> sPoints = as.getPointsOnSkeletal(); // Get skeletal points.
            vector<double> wAreas = as.getBoundaryPointsAreas();

            int pointsNums =bPoints.size();
            // Make room for a pointsNums*3 row vector. Each row will store x or y or z for the correspondence spoke on q's srep.
            if(!reSize){//Only need resize one time.
                for(unsigned i=0; i< pointsNums*3; i++){
                    pMatrix.push_back(VectorSRepFeaturesType());
                }
                for(unsigned i=0; i<pointsNums; i++){
                    wMatrix.push_back(VectorSRepFeaturesType());
                }
                for(unsigned i=0; i<sPoints.size()*3; i++){ // need multiply 3 because each point has xyz.
                    sMatrix.push_back(VectorSRepFeaturesType());
                }
                reSize=true;
            }
            int pRowIndex = 0;            

            for(unsigned i = 0; i< bPoints.size(); i++){
                pMatrix[pRowIndex].push_back(bPoints[i].getX());
                pMatrix[pRowIndex+1].push_back(bPoints[i].getY());
                pMatrix[pRowIndex+2].push_back(bPoints[i].getZ());

                pRowIndex = pRowIndex+3;

                wMatrix[i].push_back(wAreas[i]);
            }

            // Store each srep's skeletal points into sMatrix.
            int sIndex =0;
            for(unsigned i = 0; i< sPoints.size(); i++){
                sMatrix[sIndex].push_back(sPoints[i].getX());
                sMatrix[sIndex+1].push_back(sPoints[i].getY());
                sMatrix[sIndex+2].push_back(sPoints[i].getZ());
                //cout<<"------------------------sPoints["<<i<<"]-----"<<sPoints[i]<<endl;
                sIndex = sIndex+3;
            }

        }

        // save pMatrix to file.
        string filename = string(srepFolder) + string("GPA_input_pMatrix.txt");
        saveGPA_Input_Matrix(filename.c_str(), pMatrix);

        // Save weight matrix to file.
        filename = string(srepFolder) + string("GPA_input_wMatrix.txt");
        saveGPA_Input_Matrix(filename.c_str(), wMatrix);

        // Save sMatrix to file
        filename = string(srepFolder) + string("GPA_input_sMatrix.txt");
        saveGPA_Input_Matrix(filename.c_str(), sMatrix);

        return EXIT_SUCCESS;
    }
    else{
        cout<<"You input a invalid srep file, please check the path or name of this srep!"<<endl;
        return EXIT_FAILURE;
    }

}


// command: ./srepAlignment /NIRAL/work/ltu/WorkSpace/data_pablo/Pablo2/lib/vtksrep/Correspondence/models/alignment/T0036-1-1-01_segmented_vent.m3d

/*int main( int argc, char* argv[] )
{

    if( argc != 2 ) {
        std::cerr << "Usage: "<< std::endl;
        std::cerr << argv[0];
        std::cerr << " <srepFolder>";
        std::cerr << std::endl;
        return EXIT_FAILURE;
    }

    const char * filename = argv[1];

    int interpolationLevel =5;

    alignsrep as(filename, interpolationLevel);

    as.weightAreaToBoundaryPoints();

    vector<double> upCrest = as.getCrestAreas(0);
    //vector<double> downCrest = as.getCrestAreas(1);


    // Save matrix to .txt file
    std::ofstream fout;
    const char * filePath = "/NIRAL/work/ltu/WorkSpace/data_pablo/Pablo2/lib/vtksrep/Correspondence/models/alignment/areas.txt";
    fout.open(filePath);

    if(fout)  {
        for(int i =0; i< upCrest.size();i++){
            fout<< upCrest[i]<<" ";
        }
        fout<<endl;

        cout<<"Successfully saved areas to: "<<filePath<<endl;
    }
    else
        cerr<<"Write out failed, cannot open the file!"<<endl;

    fout.close();

    // Draw the spokes
    as.drawCrestSpokes_method2();


    vector<double> upQuads = as.getQuadsAreas(0);
    //vector<double> downQuads = as.getQuadsAreas(1);

    // Save matrix to .txt file
    std::ofstream fout;
    const char * filePath = "/NIRAL/work/ltu/WorkSpace/data_pablo/Pablo2/lib/vtksrep/Correspondence/models/alignment/quad-areas.txt";
    fout.open(filePath);

    if(fout)  {
        for(int i =0; i< upQuads.size();i++){
            fout<< upQuads[i]<<" ";
        }
        fout<<endl;

        cout<<"Successfully saved areas to: "<<filePath<<endl;
    }
    else
        cerr<<"Write out failed, cannot open the file!"<<endl;

    fout.close();


    // Draw the up boundary and spokes.
    as.drawBoundary();

}*/
















  /**/

  //

/*    MeshType::Pointer  mesh = MeshType::New();

    // Create points
    MeshType::PointType p0,p1,p2,p3;

    p0[0]= -1.0; p0[1]= -1.0; p0[2]= 0.0; // first  point ( -1, -1, 0 )
    p1[0]=  1.0; p1[1]= -1.0; p1[2]= 0.0; // second point (  1, -1, 0 )
    p2[0]=  1.0; p2[1]=  1.0; p2[2]= 0.0; // third  point (  1,  1, 0 )
    p3[0]=  1.0; p3[1]=  1.0; p3[2]= 1.0; // fouth  point (  1,  1, 1 )

    mesh->SetPoint( 0, p0 );
    mesh->SetPoint( 1, p1 );
    mesh->SetPoint( 2, p2 );
    mesh->SetPoint( 3, p3 );

    std::cout << "Points = " << mesh->GetNumberOfPoints() << std::endl;

    // Access points
    typedef MeshType::PointsContainer::Iterator     PointsIterator;

    PointsIterator pointIterator = mesh->GetPoints()->Begin();

    PointsIterator end = mesh->GetPoints()->End();
    while( pointIterator != end )
      {
      MeshType::PointType p = pointIterator.Value();  // access the point
      std::cout << p << std::endl;                    // print the point
      ++pointIterator;                                // advance to next point
      }

*/






    /*typedef double PixelType;
      const unsigned int Dimension = 3;
      typedef itk::PointSet< PixelType, Dimension >   PointSetType;
      typedef PointSetType::PointType PointType;
      typedef PointSetType::PointsContainerPointer PointsContainerPointer;

      PointSetType::Pointer  PointSet = PointSetType::New();
      PointsContainerPointer  points = PointSet->GetPoints();

      // Create points
      PointType p0, p1, p2;

      p0[0]=  0.0; p0[1]= 0.0; p0[2]= 0.0;
      p1[0]=  0.1; p1[1]= 0.0; p1[2]= 0.0;
      p2[0]=  0.0; p2[1]= 0.1; p2[2]= 0.0;

      points->InsertElement(0, p0);
      points->InsertElement(1, p1);
      points->InsertElement(2, p2);

      return EXIT_SUCCESS;*/




  /*typedef MeshType::PointsContainer         PointsContainer;
  typedef MeshType::PointsContainerPointer  PointsContainerPointer;

  PointsContainerPointer points = PointsContainer::New();
  points->Reserve( 100 );

  typedef MeshType::PointType               PointType;
  PointType p;
  p[2] = 0.;

  typedef MeshType::PointIdentifier         PointIdentifier;
  PointIdentifier k = 0;

  for( int i = 0; i < 10; i++ )
    {
    p[0] = static_cast< CoordType >( i );

    for( int j = 0; j < 10; j++ )
      {
      p[1] = static_cast< CoordType >( j );
      points->SetElement( k, p );
      k++;
      }
    }

  mesh->SetPoints( points );

  k = 0;

  for( int i = 0; i < 9; i++ )
    {
    for( int j = 0; j < 9; j++ )
      {
//      mesh->AddFaceTriangle( k, k+1,  k+11 );
//      mesh->AddFaceTriangle( k, k+11, k+10 );
      k++;
      }
    k++;
    }

  typedef itk::MeshFileWriter< MeshType > WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( outputFileName );
  writer->SetInput( mesh );
  try
    {
    writer->Update();
    }
  catch( itk::ExceptionObject & error )
    {
    std::cerr << "Error: " << error << std::endl;
    return EXIT_FAILURE;
    }

  return EXIT_SUCCESS;*/





