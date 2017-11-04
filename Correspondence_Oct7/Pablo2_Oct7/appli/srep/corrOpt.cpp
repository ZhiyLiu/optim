/* Author: Liyun Tu
 * Date: Sep 27, 2013
 * e.g. /work/ltu/WorkSpace/May15/binPablo2/corrOpt /NIRAL/work/ltu/WorkSpace/input_and_temp_data/upSide 0 2 -useVarFile
*/

#include <iostream>

#include "optimizationusingnewuoa.h"
#include "oneplusonecostfunction.h"
#include "toolsfunc.h"

//using namespace std;


// Call the optimizer, this function only execute once.
int callOnePlusOneOptimizer(const char* rootDir, vector< std::string > inputSreps, int interpolationLevel, int regEntType, bool initialOpt,
                            const char* logFileName){

    typedef  itk::OnePlusOneEvolutionaryOptimizer  OptimizerType;
    // Declaration of a itkOptimizer
    OptimizerType::Pointer  itkOptimizer = OptimizerType::New();

    typedef itk::Statistics::NormalVariateGenerator GeneratorType;
    GeneratorType::Pointer generator = GeneratorType::New();
    generator->Initialize(12345); // This can be any value, get a same random value when we rerun the program.
    //generator->Initialize(time(NULL)); // The random value will not the same when we rerun the program each time.
    itkOptimizer->SetNormalVariateGenerator( generator );

    // Specify the initial radius, the shrink and grow mutation factors and termination criteria of the optimizer.
    double initRadius = 0.1; // search radius in parameter space
    double grow = 1.1; // search radius grow factor, should > 1.0
    double shrink = pow(grow, -0.25); // search radius shrink factor, should < 1.0
    itkOptimizer->Initialize(initRadius, grow, shrink); //Before running this optimizer, this initialize function should have been called.
    //Optimization will stop when meets either of two termination conditions: the maximum iteration limit or epsilon.
    itkOptimizer->SetEpsilon(1.0e-6); // minimal search radius.
    itkOptimizer->SetMaximumIteration(30000); //iteration times.
    itkOptimizer->SetMaximize(false); //should Minimize the metric. (We want to minimize the cost function.)

    std::ofstream out(logFileName, std::ios_base::out | std::ios_base::app);
    std::streambuf *coutbuf = std::cout.rdbuf(); //save old buf
    std::cout.rdbuf(out.rdbuf()); //redirect std::cout to out.txt!

    // Declaration of the Cost Function
    oneplusonecostfunction  costFunction(rootDir, inputSreps, interpolationLevel, regEntType);
    itkOptimizer->SetCostFunction( &costFunction); // Pass pointer.

    typedef oneplusonecostfunction::ParametersType    ParametersType;

    const unsigned int spaceDimension =  costFunction.GetNumberOfParameters(); //number of the deltau&deltav variables.
    //int variablesNum = spaceDimension * inputSreps.size();

    // Set initial position.
    ParametersType  initialPosition( spaceDimension );

    // for the first optimization, all the variables set to 0.
    if(initialOpt){
        for(unsigned int i = 0; i < spaceDimension; i++){
            initialPosition[i] = 0;
        }
    }
    else {
        string varFilePath = string(rootDir) + string("/vars.txt");
        ifstream ss(varFilePath.c_str());

        if (! ss.is_open()){
            std::cerr << "Msg from optimizationusingnewuoa::callNEWUOA: Unable to open the vars file \"" << varFilePath << "\"!" << std::endl;
            EXIT_FAILURE;
        }

        string token;

        int i = 0;
        while (ss >> token && i < spaceDimension) {
            initialPosition[i] = atof(token.c_str());
            i++;
        }
        std::cout<<"----Finished reading variables from "<<varFilePath<<std::endl;
    }

    itkOptimizer->SetInitialPosition( initialPosition );

    // Create the Command observer and register it with the optimizer.
    //CommandIterationUpdate::Pointer observer = CommandIterationUpdate::New();
    //itkOptimizer->AddObserver( itk::IterationEvent(), observer );

    typedef OptimizerType::ScalesType        ScalesType;
    ScalesType scales(spaceDimension); // variables numbers. scale diff types of paras to compensate for the differences in the dynamic ranges of the paras.
    scales.Fill( 1.0 );
    itkOptimizer->SetScales( scales );

    try  {
        itkOptimizer->StartOptimization();
    }
    catch( itk::ExceptionObject & e )  {
        std::cout << "Exception thrown ! " << std::endl;
        std::cout << "An error ocurred during Optimization" << std::endl;
        std::cout << "Location    = " << e.GetLocation()    << std::endl;
        std::cout << "Description = " << e.GetDescription() << std::endl;
        return EXIT_FAILURE;
    }

    // output the minimum position after optimization
    ParametersType finalPosition = itkOptimizer->GetCurrentPosition();
    for(unsigned int i = 0; i < spaceDimension; i++){
        std::cout << finalPosition[i] << " ";
    }
    std::cout<<std::endl;

    //reset to standard output
    std::cout.rdbuf(coutbuf);

    return EXIT_SUCCESS;
}



bool readSreps ( const char* listName,  vector< std::string > &inputSreps) {
  std::ifstream shapeListFile;
  shapeListFile.open( listName, std::ios_base::in );
  if (!shapeListFile) {
    std::cerr << "Unable to open shape list \"" << listName << "\"!" << std::endl;
    return false ;
  }

  // parse the shape list file
  std::string currentFileName;
  while (!shapeListFile.eof()) {
    std::getline( shapeListFile, currentFileName );
    if (currentFileName.empty() || currentFileName[0]=='#') { continue; }

    inputSreps.push_back(currentFileName);
  }

  shapeListFile.close();

  return true;
}


// call the optimizer
int main(int argc, char *argv[]){

    if(argc < 5) {
      std::cout << "Usage: " << argv[0] << " WorkingFolder SrepFileList regularityEntropyType logFileName [-interpLevel] [-useVarFile] [-useOnePlusOne]" << std::endl << std::endl ;
      std::cout << "WorkingFolder: where the 'moved_sreps', 'vars.txt', 'temp_data' stored" << std::endl;
      std::cout << "SrepFileList: the full path of the srep file name list file (input_sreps.txt)" << std::endl;
      std::cout << "regularityEntropyType: 0(up side), 1(down side), 2(crest side)" << std::endl;
      std::cout << "logFileName: the log file pathname that store the variables during iteration" << std::endl;
      std::cout << "interpLevel: interpolation level (default 2)" << std::endl;
      std::cout << "-useVarFile: if choose, read variables from given file (WorkingFolder/vars.txt)" << std::endl;
      std::cout << "-useOnePlusOne: if choose, will use oneplusone optimizer instead of NUEUOA (default use NEWUOA)" << std::endl;
      std::cout << std::endl;
      return -1;
    }

    int interpolationLevel = 2; // default use 2

    bool initialOpt = true;
    bool useNEWUOA = true;

    if ( argc > 5 )  {
      for ( short int i = 5 ; i < argc ; i++ )  {
        if ( !strcmp ( "-interpLevel", argv[i] ) ) {
            interpolationLevel = atoi(argv[i+1]);
        }
        if ( !strcmp ( "-useVarFile", argv[i] ) ) {
            initialOpt = false;
        }
        if ( !strcmp ( "-useOnePlusOne", argv[i] ) ) {
            useNEWUOA = false;
        }
      }
    }  

    const char* rootDir = argv[1];

    const char* SrepFileList = argv[2];
    vector< std::string > inputSreps;

    bool readfinish = readSreps(SrepFileList, inputSreps);
    if(!readfinish){
        std::cerr << "Something wrong while reading sreps from \"" << SrepFileList << "\"!" << std::endl;
        return EXIT_FAILURE;
    }

    int regularityEntropyType = atoi(argv[3]); // 0(up regularity entropy), 1(down regularity entropy), 2(crest regularity entropy).

    clock_t start_time = clock();

    if(useNEWUOA){
        optimizationusingnewuoa newuoa_instance(rootDir, inputSreps, interpolationLevel, regularityEntropyType);
        newuoa_instance.callNEWUOA(argv[4], initialOpt);
    }
    else{
        callOnePlusOneOptimizer(rootDir, inputSreps, interpolationLevel, regularityEntropyType, initialOpt, argv[4]);
    }

    clock_t end_time = clock();
    toolsfunc tls;
    std::cout << "Time elapsed for this iteration: " << double(tls.RunningTime(end_time, start_time)) << " ms";

    return 0;
}
