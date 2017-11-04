#include <itkMetaMeshConverter.h>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include "srepCorrespondenceEvaluator.h"

using namespace std;

//////////////////////////////////////////////////////////////////////////////////

typedef itk::DefaultStaticMeshTraits <double, 3, 3, double, double> traitsType ;
typedef itk::Mesh <double, 3, traitsType> meshType ;
typedef itk::MeshSpatialObject <meshType> meshSOType ;
typedef itk::MetaMeshConverter <3, double, traitsType> MeshConverterType ;
typedef itk::srepCorrespondenceEvaluator <meshType> evaluatorType ;

typedef vnl_matrix<double>                MatrixType;

//////////////////////////////////////////////////////////////////////////////////

bool readSrepProperties ( string filename, evaluatorType::Pointer evaluator, unsigned int &nSamples) {

    // open the property file
    std::ifstream propertyFile;
    propertyFile.open( filename.c_str(), std::ios_base::in );
    if (!propertyFile) {
        std::cerr << "Unable to open the s-rep property file \"" << filename << "\"!" << std::endl;
        return false ;
    }

    // parse the property file    
    vector<string> propertyVector;

    string currentLine;
    while (!propertyFile.eof()) {
        // read a line and store in "currentLine"
        std::getline( propertyFile, currentLine );

        if (!currentLine.empty()) {
            propertyVector.push_back(currentLine);
        }
    }
    int dimNum = propertyVector.size(); // property number
    cout<<"-------------------dimNum is: "<<dimNum<<endl;

    // parse string to get the smaple number
    string temp_str = propertyVector[0];
    stringstream ssin(temp_str);
    string token;
    nSamples = 0; // srep number
    vector<double> proRow;
    while (ssin >> token) {
        proRow.push_back(atof(token.c_str()));
        nSamples++;
    }
    cout<<"-------------------nSamples is: "<<nSamples<<endl;

    // property matrix of the training set
    MatrixType propertyMatrix(dimNum, nSamples);
    for(unsigned int i = 0; i < dimNum; i++){
        string temp_str = propertyVector[i];
        stringstream ss(temp_str);
        string val;
        int srepIndex = 0;
        while (ss >> val && srepIndex < nSamples) {
            propertyMatrix[i][srepIndex] = atof(val.c_str());

            srepIndex++;
        }
    }

    evaluator->SetNumberOfInputs(nSamples);

    evaluator->SetInput(propertyMatrix);

    return true ;
}

int main( int argc, char *argv[] ) {
    if (argc < 3)  {
        std::cout << "Usage: " << argv[0] << " MeshListFile OutputFile [-gen] [-spec] [-uniform]" << std::endl << std::endl ;
        std::cout << "MeshListFile: List of shapes in the population" << std::endl;
        std::cout << "OutputFile: Where to store the evaluation results" << std::endl;
        std::cout << "-gen: Compute generalization" << std::endl;
        std::cout << "-spec: Compute specificity" << std::endl;
        std::cout << "-uniform: Use uniform distribution for random numbers (default is Gaussian)" << std::endl ;
        std::cout << "-numSpec: number of samples created for specificity" << std::endl ;
        std::cout << std::endl ;
        return -1;
    }

    // parse command line arguments
    std::string meshListFileName = argv[1];
    std::string outputFileName = argv[2];

    bool computeGeneralization = false ;
    bool computeSpecificity = false ;
    bool gaussian = true ;

    int N = 1000; // Default create 1000 instance for specificity.

    if ( argc > 3 )  {
        for ( short int i = 3 ; i < argc ; i++ )  {
            if ( (!computeGeneralization) && ( !strcmp ( "-gen", argv[i] ) ) )  {
                computeGeneralization = true ;
            }
            if ( (!computeSpecificity) && ( !strcmp ( "-spec", argv[i] ) ) )  {
                computeSpecificity = true ;
            }
            if ( (gaussian) && ( !strcmp ( "-uniform", argv[i] ) ) ) {
                gaussian = false ;
            }
            if (( !strcmp ( "-numSpec", argv[i] ) ) )  {
                N = atoi(argv[i+1]);
            }
        }
    }

    // read input meshes
    unsigned int numSamples ;
    evaluatorType::Pointer evaluator = evaluatorType::New() ;

    if ( ! readSrepProperties ( meshListFileName, evaluator, numSamples ) ) {
        return -1 ;
    }
    std::cout << "reading done " << std::endl;

    // evaluate correspondence
    std::vector < double > generalization, specificity ;
    std::vector < double > generalizationError, specificityError ;
    generalization.resize ( numSamples - 1 ) ;
    specificity.resize ( numSamples ) ;
    generalizationError.resize ( numSamples - 1 ) ;
    specificityError.resize ( numSamples ) ;

    if ( computeGeneralization )  {
        for ( unsigned int i = 0  ; i < numSamples - 1 ; i++ )  {
            std::cout << "Starting generalization computation with " << i << " shape parameters... " ;
            generalization[i] = evaluator->GetGeneralization ( i, generalizationError[i] ) ;
            std::cout << "Finished." << std::endl ;
        }
    }
    if ( computeSpecificity )  {
        evaluator->SetDistribution ( gaussian ) ;
        for ( unsigned int i = 0 ; i < numSamples ; i++ ) {
            std::cout << "Starting specificity computation with " << i << " shape parameters... " ;
            specificity[i] = evaluator->GetSpecificity ( N, i, specificityError[i] ) ;
            std::cout << "Finished." << std::endl ;
        }
    }

    // write out evaluation results
    std::ofstream outputFile ;
    outputFile.open ( outputFileName.c_str() ) ;
    outputFile << "Correspondence evaluation for: " << argv[1] << std::endl << std::endl ;
    if ( computeGeneralization )  {
        outputFile << "Generalization: " << std::endl ;
        for ( unsigned int i = 0 ; i < numSamples - 1 ; i++ ) {
            outputFile << i << " " << generalization[i] << " " << generalizationError[i] << std::endl ;
        }
    }
    if ( computeSpecificity ) {
        outputFile << "Specificity: " << std::endl ;
        for ( unsigned int i = 0 ; i < numSamples ; i++ ) {
            outputFile << i << " " << specificity[i] << " " << specificityError[i] << std::endl ;
        }
    }
    outputFile.close ();

    return 0 ;
}



