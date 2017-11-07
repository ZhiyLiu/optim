/* The purpose of this class is to test newuoa optimizer
 *
 * Zhiyuan Liu
 * 2017.11
 */
#include <string>
#include <vector>
#include <iostream>
#include "M3DNewuoaOptimizer.h"
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


int main(int argc, char *argv[]) {
    using namespace std;
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
        M3DNewuaOptimizer newuoa_instance(rootDir, inputSreps, interpolationLevel, regularityEntropyType);
        newuoa_instance.perform(argv[4], initialOpt);
    }
    else{
//        callOnePlusOneOptimizer(rootDir, inputSreps, interpolationLevel, regularityEntropyType, initialOpt, argv[4]);
    }

    clock_t end_time = clock();
//  toolsfunc tls;
//  std::cout << "Time elapsed for this iteration: " << double(tls.RunningTime(end_time, start_time)) << " ms";

    return 0;

}
