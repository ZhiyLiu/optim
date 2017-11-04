#ifndef REGULARITYENTROPY_H
#define REGULARITYENTROPY_H

#include <vector>
#include "vnl/algo/vnl_symmetric_eigensystem.h"
#include "itkMatrix.h"

//#include "macro_const_values.h"

//#define PI 3.1415926535

//double CONST_C = 1.4189; //(1 + log(2*PI))/2;


class regularityentropy
{
public:

    typedef vnl_matrix<double> MatrixType;


    regularityentropy();

    double calculateEntropy(MatrixType featureMatrix_A, double threshold);


private:
   // MatrixType featureMatrix_A; // holding the features for computing a specifice kind of entropy. Row is feature dimesion, column is samples.
};

#endif // REGULARITYENTROPY_H
