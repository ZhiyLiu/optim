/* Compute the regularity entropy.
 * Input: feature matrix.
 * Output: the entropy of this feature matrix
 *
 * Liyun Tu
 * Apr 11, 2014
*/

#include "regularityentropy.h"


using namespace std;




regularityentropy::regularityentropy()
{
}





/* Return the entropy of the input feature matrix.
 * featureMatrix_A: m-by-n, holding the features for computing a specifice kind of entropy. m is feature dimesion, n is samples.
 * threshold: the contribution of an eigenvalue bigger than this threshold will be keep, otherwise will be omit.
 * Here we recommand set it to 0.01.
*/
double regularityentropy::calculateEntropy(MatrixType featureMatrix_A, double threshold){

    int m = featureMatrix_A.rows();
    int n = featureMatrix_A.columns();

    //Step 1: Compute the mean for each feature over the n samples.
    double rowMean[m];
    for(unsigned int i =0; i< m; i++){
        double sumRow = 0;
        for(unsigned int j =0; j< n; j++){
            sumRow += featureMatrix_A[i][j];
        }
        rowMean[i] = sumRow/n;
    }

    //Step 2: Center the data
    for(unsigned int j =0; j< n; j++){ // For each smaple
        for(unsigned int i =0; i< m; i++){ //For each feature
            featureMatrix_A[i][j] -= rowMean[i];
        }
    }

    //Step 3: Calculate the covariance matrix. m-by-m, symmetric.
    MatrixType covMatrix = (featureMatrix_A * featureMatrix_A.transpose()) / (n-1); //use n-1 instead of n.
    /*cout<<"------------covMatrix is: "<<endl;
    for(unsigned int i = 0; i < covMatrix.rows(); i++){ // equals to m
        for(unsigned int j =0; j< covMatrix.columns(); j++){ //equals to m
            cout<<covMatrix[i][j]<<"  ";
        }
        cout<<endl;
    }*/

    //Step 4: Normalized covariance matrix by standard deviation to get correlation matrix (dimensionless, symmetric)
    // Get the standard deviation of each dimesion based on the centered data.
    //double sd[m]; // standard deviation of each diemesion(each feature).
    vector<double> sd;
    sd.resize(m);
    for(unsigned int i =0; i< m; i++){ //For each feature
        double squaredSum =0;
        for(unsigned int j =0; j< n; j++){ // For each smaple
            squaredSum += featureMatrix_A[i][j] * featureMatrix_A[i][j]; // squared sum of each dimesion.
        }
        sd[i] = sqrt(squaredSum / (n-1)); //use n-1 instead of n. Exactly same with matlab function corrcov.
        //cout<<"--------------sd["<<i<<"] is: "<<sd[i]<<endl;
    }

    // entries in correlation matrix equals to correspondence entries in covMatrix divided by its two features's sd.
    MatrixType corMatrix(covMatrix.rows(),covMatrix.columns()); //m-by-m
    for(unsigned int i = 0; i < covMatrix.rows(); i++){ // equals to m
        for(unsigned int j =0; j< covMatrix.columns(); j++){ //equals to m
            corMatrix[i][j] = covMatrix[i][j]/(sd[i]*sd[j]); // divide by its correspondence two sd.            
        }       
    }

    //Step 5: Find the eigenvectors(PC, column vector, each column is a PC) and eigenvalues(diagonal of V) of corMatrix.
    vnl_symmetric_eigensystem<double> eig(corMatrix);
    vnl_diag_matrix<double> eigenValues(m, m);
    eigenValues = eig.D; //eig.D is default by increasing order.

    /*cout<<"eigenValues is: "<<endl;
    for(unsigned r=0; r<eigenValues.rows();r++){
        cout<<eigenValues(r,r)<<endl;
    }*/

    //Step 6: Select eigenvalues (only use eigenvalues bigger than theta) and Compute the entropy.
    // Sum all the m eigenvalues.
    double sumEig = 0;
    for(unsigned int i = 0; i< m; i++){
        sumEig += eigenValues(i,i);
    }

    double sumLogNubda = 0;
    int counter = 0;
    for(unsigned int i = 0; i< m; i++){
        double contribute = eigenValues(i,i) / sumEig;
        if(contribute > threshold) { /* usually threshold is set to 0.01*/
            //Compute the entropy
            sumLogNubda += log(eigenValues(i,i));

            counter++;
        }
    }

    double regEntropy = sumLogNubda*0.5 + counter * 1.4189; //(1 + log(2*PI))/2;

    return regEntropy;
}



