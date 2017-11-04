/* This class implement the  Generalization, Specificity and Compactness metric described in Ipek's dissertaion 2.3.3.
 * The results should be the same as toolbox "CorrespondenceEvaluator"
*/


#include "correspondenceevaluation.h"

correspondenceevaluation::correspondenceevaluation()
{

}


/* Main logic of correspondence evaluation.
*/
//void correspondenceevaluation::corrEva(const char* rootDir){
//    // Read in data and put into matrix
//    MatrixType X = readVTKPointData(rootDir);

//    //





//}






/* Compute G(M)
*/
//void correspondenceevaluation::generalization(){

//}



///* Compute S(M)
//*/
//void correspondenceevaluation::specificity(){

//}



///* Compute C(M)
//*/
//void correspondenceevaluation::compactness(){

//}


/* Read point set into martix.
 * rootDir: contains the files of .vtk point.
*/
/*MatrixType correspondenceevaluation::readVTKPointData(const char* rootDir){
    vector<MatrixType> TS; // training set

    // Get the srep filenames from given folder
    vector< std::string > pointFiles;
    pointFiles = tls.getModelname(string(rootDir));

    // For each srep
    for(unsigned int i = 0; i < pointFiles.size(); i++){
        // Read in the source points set
        vtkSmartPointer<vtkPolyData> polyData = vtkSmartPointer<vtkPolyData>::New(); // same as vtkPolyData* polyData_source;

        vtkSmartPointer<vtkPolyDataReader> reader = vtkSmartPointer<vtkPolyDataReader>::New();

        reader->SetFileName(pointFiles[i].c_str());
        reader->Update();
        polyData = reader->GetOutput();

        int pointNum = polyData->GetNumberOfPoints();
        MatrixType X(3, pointNum);

        // Read in the points set
        for(unsigned int i = 0; i < pointNum; i++){
            double p[3];
            polyData->GetPoint(i,p);

            for(unsigned int dim = 0; dim < 3; dim++){
                X[dim][i] = p[dim];
            }

            //std::cout << "Point " << id_s << " : (" << p1[0] << " " << p1[1] << " " << p1[2] << ")" << std::endl;
        }

        TS.push_back(X);
    }

    // Convert TS from vector of size N to a matrix of 3*pointNum-by-N.
    int srepNum = TS.size();
    int pointNum = TS[0].colums();

    MatrixType X(3*pointNum, srepNum);

    // For each srep
    for(unsigned int i = 0; i<srepNum; i++){
        MatrixType temp = TS[i];
        if(temp.colums() == pointNum){
            // For each point
            for(unsigned int j = 0; j < pointNum; j++){
                for(unsigned int dim = 0; dim<3; dim++){
                    X[3*j+dim][i] = temp[dim][j];
                }
            }
        }
        else {
            cout<<"Msg from correspondenceevaluation::getDataMatrix: Wrong size!! Exit!"<<endl;
            exit(1);
        }
    }

    return X;
}*/






