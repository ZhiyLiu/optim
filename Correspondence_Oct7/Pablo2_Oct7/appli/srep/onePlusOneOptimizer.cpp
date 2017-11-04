/* Author: Liyun Tu
 * Date: Apr, 2014
 * Call command: onePlusOneOptimizer -r /NIRAL/work/ltu/WorkSpace/data_pablo/Pablo2
 * can input many model at a same time, use whitespace between these models.
*/

#include <iostream>
#include "oneplusonecostfunction.h"

using namespace std;

// Call the optimizer, this function only execute once.
int callOnePlusOneOptimizer(const char* rootDir, vector< std::string > inputSreps, int interpolationLevel, int regEntType, bool initialOpt){

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
    itkOptimizer->SetMaximumIteration(25000); //iteration times.
    itkOptimizer->SetMaximize(false); //should Minimize the metric. (We want to minimize the cost function.)

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
        cout<<"----finished read variables from "<<varFilePath<<endl;
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
    cout<<endl;

    return EXIT_SUCCESS;
}



// call the optimizer
int main(int argc, char *argv[]){

    //arg[0] is exe name, arg[1] is "-r", the following is the model names.

    string dir = string(argv[2]) + string("/lib/vtksrep/Correspondence");
    const char* rootDir = dir.c_str();

    vector< std::string > inputSreps;

    // Set para
    int interpolationLevel = 2;
    int regularityEntropyType = 2; // 0(up regularity entropy), 1(down regularity entropy), 2(crest regularity entropy).
    bool initialOpt = true;

    //----------------------Comment this for test.-----------------------
    //Initialize the original modelname.
    //each time the deltaU, deltaV move based on "/NIRAL/work/ltu/WorkSpace/data_pablo/Pablo2/lib/vtksrep/Correspondence/models/testSet_Original".
    //string originalFilePath = "/NIRAL/work/ltu/WorkSpace/data_pablo/Pablo2/lib/vtksrep/Correspondence/models/testSet_Original/";
    //inputSreps = tools.getModelname(originalFilePath);

    /*------------Begin----------------For test only------------------------------------*/
    inputSreps.push_back(string("N0012-2-01_segmented_vent_pp_surfSPHARM_scale_atomStage_S.m3d"));
    inputSreps.push_back(string("N0052-2-02_segmented_vent_pp_surfSPHARM_scale_atomStage_S.m3d"));
    inputSreps.push_back(string("N0092-1-01_segmented_vent_pp_surfSPHARM_scale_atomStage_S.m3d"));
    inputSreps.push_back(string("N0092-1-02_segmented_vent_pp_surfSPHARM_scale_atomStage_S.m3d"));
    inputSreps.push_back(string("N0106-1-03_segmented_vent_pp_surfSPHARM_scale_atomStage_S.m3d"));
    inputSreps.push_back(string("N0147-2-02_segmented_vent_pp_surfSPHARM_scale_atomStage_S.m3d"));
    inputSreps.push_back(string("N0156-1-02_segmented_vent_pp_surfSPHARM_scale_atomStage_S.m3d"));
    inputSreps.push_back(string("N0156-1-03_segmented_vent_pp_surfSPHARM_scale_atomStage_S.m3d"));
    inputSreps.push_back(string("N0174-1-01_segmented_vent_pp_surfSPHARM_scale_atomStage_S.m3d"));
    inputSreps.push_back(string("N0174-1-02_segmented_vent_pp_surfSPHARM_scale_atomStage_S.m3d"));
    inputSreps.push_back(string("N0177-1-02_segmented_vent_pp_surfSPHARM_scale_atomStage_S.m3d"));
    inputSreps.push_back(string("N0188-1-03_segmented_vent_pp_surfSPHARM_scale_atomStage_S.m3d"));
    inputSreps.push_back(string("N0227-1-02_segmented_vent_pp_surfSPHARM_scale_atomStage_S.m3d"));
    inputSreps.push_back(string("N0227-1-03_segmented_vent_pp_surfSPHARM_scale_atomStage_S.m3d"));
    inputSreps.push_back(string("T0017-1-2-02_segmented_vent_pp_surfSPHARM_scale_atomStage_S.m3d"));
    inputSreps.push_back(string("T0039-1-2-02_segmented_vent_pp_surfSPHARM_scale_atomStage_S.m3d"));
    inputSreps.push_back(string("T0045-1-1-03_segmented_vent_pp_surfSPHARM_scale_atomStage_S.m3d"));
    inputSreps.push_back(string("T0058-1-1-01_segmented_vent_pp_surfSPHARM_scale_atomStage_S.m3d"));
    inputSreps.push_back(string("T0060-1-1-02_segmented_vent_pp_surfSPHARM_scale_atomStage_S.m3d"));
    inputSreps.push_back(string("T0060-1-1-03_segmented_vent_pp_surfSPHARM_scale_atomStage_S.m3d"));
    inputSreps.push_back(string("T0060-1-2-02_segmented_vent_pp_surfSPHARM_scale_atomStage_S.m3d"));
    inputSreps.push_back(string("T0062-1-1-02_segmented_vent_pp_surfSPHARM_scale_atomStage_S.m3d"));
    inputSreps.push_back(string("T0064-1-1-02_segmented_vent_pp_surfSPHARM_scale_atomStage_S.m3d"));
    inputSreps.push_back(string("T0068-1-2-02_segmented_vent_pp_surfSPHARM_scale_atomStage_S.m3d"));
    inputSreps.push_back(string("T0074-1-1-01_segmented_vent_pp_surfSPHARM_scale_atomStage_S.m3d"));
    inputSreps.push_back(string("T0074-1-1-02_segmented_vent_pp_surfSPHARM_scale_atomStage_S.m3d"));
    inputSreps.push_back(string("T0074-1-1-03_segmented_vent_pp_surfSPHARM_scale_atomStage_S.m3d"));
    inputSreps.push_back(string("T0081-1-1-01_segmented_vent_pp_surfSPHARM_scale_atomStage_S.m3d"));
    inputSreps.push_back(string("T0081-1-2-02_segmented_vent_pp_surfSPHARM_scale_atomStage_S.m3d"));

    cout << "-----------------Call One Plus One Evolutionary Optimizer... "<<endl;

    callOnePlusOneOptimizer(rootDir, inputSreps, interpolationLevel, regularityEntropyType, initialOpt);

    return 0;
}






