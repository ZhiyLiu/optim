/* Prepare input shape location matrixs for correspondence evaluation method.
 * Compute the Generalization, Specificity and Compactness using Matin's tool named: CorrespondenceEvaluator, it require a input of:
 *
 * MeshListFile: List of shapes in the population
 * OutputFile: Where to store the evaluation results
 * -gen: Compute generalization
 * -spec: Compute specificity
 * -uniform: Use uniform distribution for random numbers (default is Gaussian)
 * Usage: CorrespondenceEvaluator MeshListFile OutputFile [-gen] [-spec] [-uniform]
 *
 * Description of these 3 statistic validation metrics refer to Ipek oguz's disertaion, P39 -P40,
 * and Bea's paper "Combined SPHARM-PDM and entropy-based particle systems shape analysis framework".
 * and Martin's paper "Evaluation of 3D Correspondence Methods for Model Building"
 *
 * This function input the folder of s-reps, and sample the skeletal + boundary points of the s-rep (leve 0) to represent the shape.
 * Since the original s-rep with their spokes share the tail, while the correspondence optimized s-rep have
 * their skeletal points seperate, do we need to keep the same point number for each object before and after if we want to compare them?
 * So we give two different choice: repeatSkeletalPoints, no repeatSkeletalPoints.
 * repeatSkeletalPoints means: count the shared tail in weighted 1, 2, 3, just as what we do in Procrustes alignment. (Briefly speaking,
 * we count each spoke's tail and tip, for the original s-rep, we count the skeletal point 3times)
 * no repeatSkeletalPoints means: the original s-rep have the shared skeletal points only counted once, this will cause a different nubmer of
 * points for each object before and after optimization.
 *
 * We can compare and see which one decreases much in the three metrics?
 *
 * Liyun Tu
 * Jul 10, 2014
 *
 * To run:
 * ./sampleSrepPointToVTK /NIRAL/work/ltu/WorkSpace/Test/30AllAtomsFromAtomStage/input/afterProcrustes/ /NIRAL/work/ltu/WorkSpace/Test/30AllAtomsFromAtomStage/input/afterProcrustes/M3DFileNameList.txt /NIRAL/work/ltu/WorkSpace/CorrespondenceEvaluation_SPIE/test/ /NIRAL/work/ltu/WorkSpace/Test/30AllAtomsFromAtomStage/after_optimized_coeffs 3 13 5 2
 *
 * To used this funtion, you must prepare:
 * 1) a M3DFileNameListFile which listing all the sreps' filename with its path.
 * 2) three string files, storing the variables for each side after shifted. (if want to sample the optimized srep).
 *
 *
 *
*/

#include "samplesreppoints.h"


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




bool readSreps(const char * listName, vector < std::string > &fileNames) {
  std::ifstream srepFileNameListFile;
  srepFileNameListFile.open(listName, std::ios_base::in);

  if (!srepFileNameListFile) {
    std::cerr << "Unable to open shape list \"" << srepFileNameListFile << "\"!" << std::endl;
    return false;
  }

  // parse the shape list file
  int numSamples = 0 ;
  std::string currentFileName;

  while (!srepFileNameListFile.eof()) {
    std::getline( srepFileNameListFile, currentFileName );
    if (currentFileName.empty() || currentFileName[0]=='#') { continue; }

    numSamples++;
    fileNames.resize( numSamples ) ;
    fileNames[numSamples - 1] = currentFileName ;
  }

  srepFileNameListFile.close();

  return true ;
}



/* pointType: 1 (only collect skeletal point), 3 (only boundary point), 5 (skeletal point + boundary point)
 *
 * If interpolationLevel = 2, pointType = 5: there will be 28*20(crest boundary points) + 28*4(fold curve points on skeletal) +
 * 441*2*2(up and down spoke tail and tip), totally 560+112+882*2 = 2436 points sampled for each srep;
 *
 * If interpolationLevel = 2, pointType = 3: there will be 28*20(crest boundary points) + 441*2(up and down spoke tip),
 * totally 560+882 = 1442 points sampled for each srep;
 *
 * If interpolationLevel = 2, pointType = 1: there will be 28*4(fold curve points on skeletal) + 441*2(up and down spoke tail),
 * totally 112+882 = 994 points sampled for each srep;
 *
 * If interpolationLevel = 0, pointType = 5: there will be 28(fold curve points on skeletal) + 39*2*2(up and down spoke tail),
 * 28(crest boundary points), totally 212 points sampled for each srep;

*/
int main( int argc, char* argv[] ){

    if (argc < 6) {
      std::cout << "Usage: " << argv[0] << " srepFolder M3DFileList outputFile coeffFilesPath rowNum colNum pointType interLevel" << std::endl << std::endl ;
      std::cout << "[1] SrepFolder: Where the srep stored" << std::endl;
      std::cout << "[2] M3DFileList: Where the file holding the srep file name list stored" << std::endl;
      std::cout << "[3] outputFileFolder: Where to store the output .vtk points files" << std::endl;
      std::cout << "[4] coeffFilesPath: where the three coefficient files stored. [For sampling before correspondence, set this to 0]" << std::endl ;
      std::cout << "[5] rowNum: Srep row count" << std::endl;
      std::cout << "[6] colNum: Srep column count" << std::endl;
      std::cout << "[7] pointType: 1 (only collect skeletal point), 3 (only boundary), 5 (skeletal + boundary)" << std::endl;
      std::cout << "[8] interLevel: interpolation level" << std::endl;
      std::cout << std::endl;
      return -1;
    }

    // parse command line arguments
    int rowNum = atoi(argv[5]);
    int colNum = atoi(argv[6]);
    int pointType = atoi(argv[7]);
    int interpolationLevel = atoi(argv[8]);

    // store the m3d file name to be sampled
    vector < string > srepNames;
    readSreps(argv[2], srepNames);

    samplesreppoints obj;
    obj.sampleSrepPoints(argv[1], srepNames, argv[3], argv[4], rowNum, colNum, pointType, interpolationLevel);

    return 0;
}

