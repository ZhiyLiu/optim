
#ifndef _srepCorrespondenceEvaluator_txx
#define _srepCorrespondenceEvaluator_txx

#include "srepCorrespondenceEvaluator.h"
#include <fstream>
#include <stdlib.h>
#include <time.h>

namespace itk
{

  template<class TOutputMeshType>
  srepCorrespondenceEvaluator<TOutputMeshType>
  ::srepCorrespondenceEvaluator()
  {
    this->m_NumberOfInputs = 0 ;
    srand ( time(0) ) ;
    this->gaussianDistribution = true ;
  }

  template<class TOutputMeshType>
  srepCorrespondenceEvaluator<TOutputMeshType>
  ::~srepCorrespondenceEvaluator()
  {
  }

  template<class TOutputMeshType>
  double srepCorrespondenceEvaluator<TOutputMeshType>
  ::GetGeneralization(unsigned int numberOfShapeParameters, double &stdError) {
    double mean = 0, stddev = 0 ;
    std::vector < double > sample ;

    if ( numberOfShapeParameters > this->m_NumberOfInputs - 2 )
    {
      itkExceptionMacro ("Invalid number of shape parameters for generalization computation.") ;
      return -1 ;
    }

    sample.resize ( this->m_NumberOfInputs ) ;

    for ( int i = 0 ; i < this->m_NumberOfInputs ; i++ )
    {
      sample[i] = this->ComputeGeneralizationError ( i, numberOfShapeParameters ) ;
      mean += sample[i] ;
    }

    // compute sample mean
    mean /= this->m_NumberOfInputs ;

    // compute sample standard deviation
    for ( int i = 0 ; i < this->m_NumberOfInputs ; i++ )
    {
        sample[i] -= mean ;
        stddev += sample[i] * sample[i] ;
    }
    stddev = sqrt ( stddev / this->m_NumberOfInputs ) ;

    // compute standard error of generalization
    stdError = stddev / sqrt ( this->m_NumberOfInputs - 1 ) ;


    return mean ;
  }

  template<class TOutputMeshType>
  double srepCorrespondenceEvaluator<TOutputMeshType>
  ::GetSpecificity(unsigned int N, unsigned int numberOfShapeParameters, double &stdError)
  {
    double specificity = 0 ;

    if ( numberOfShapeParameters > this->m_NumberOfInputs - 1 ) {
      itkExceptionMacro ("Invalid number of shape parameters for generalization computation.") ;
      return -1 ;
    }

    // build model
    int proDim = this->featureMatrix.rows();
    MatrixType rowMean(proDim, 1);
    rowMean.fill(0);
    MatrixType pcaMatrix(proDim, this->m_NumberOfInputs);

    for (unsigned int row = 0; row < proDim; row++) {
        // compute mean of the training set
        for ( unsigned int sampleIdx = 0 ; sampleIdx < this->m_NumberOfInputs ; sampleIdx++ )  {
            rowMean[row][0] += this->featureMatrix[row][sampleIdx];
        }
        rowMean[row][0] /= ( this->m_NumberOfInputs );

        // build pca matrix, the matrix where we are going to perform pca to find principal modes of variation of the model
        for ( unsigned int sampleIdx = 0 ; sampleIdx < this->m_NumberOfInputs ; sampleIdx++ )  {
            (pcaMatrix)[row][sampleIdx] = this->featureMatrix[row][sampleIdx] - rowMean[row][0];
        }
    }

    std::vector<double> SqrtEigenValues;
    MatrixType PrincipalDirections;

    SqrtEigenValues.resize ( this->m_NumberOfInputs ) ;

    vnl_svd<double> svd( pcaMatrix );
    for (unsigned int m=0; m<this->m_NumberOfInputs-1; m++)  // why -1 ?
    {
      double sv = svd.W( m );
      SqrtEigenValues[m] = sqrt ( sv ) ;
    }
    PrincipalDirections = svd.U();

    std::vector < double > sampleErrors ;
      sampleErrors.resize ( N ) ;

    for ( int i = 0 ; i < N ; i++ ) {
        // generate random shape parameters
        MatrixType b(numberOfShapeParameters,1) ;
        b.fill(0);
        for ( unsigned int m = 0 ; m < numberOfShapeParameters ; m++ ) {
            //b[m] = 3 * ( ( rand () / (double) RAND_MAX ) * 2 * SqrtEigenValues[m] - SqrtEigenValues[m] ) ;
            if ( this->gaussianDistribution )
            {
                b[m][0] = this->GenerateGaussianRandomNumber() * SqrtEigenValues[m] ;
            }
            else
            {
                b[m][0] = this->GenerateUniformRandomNumber() * 3 * SqrtEigenValues[m] ;
            }
        }

        // construct random shape
        MatrixType reconstruction(proDim, 1);

        for ( unsigned int dim = 0 ; dim < proDim ; dim++ )  {
                reconstruction[dim][0] = rowMean[dim][0];

                for ( unsigned int m = 0 ; m < numberOfShapeParameters; m++ ){
                    // add the component associated with the mth mode of variation
                    reconstruction[dim][0] += PrincipalDirections[dim][m]*b[m][0];
                }
        }

        // find nearest shape in the training set and compute error
        int minJ = -1 ;
        double minDist = 999999999 ;
        double dist ;
        MatrixType currentMatrix;
        currentMatrix.set_size(proDim, 1);

        for ( int j = 0 ; j < this->m_NumberOfInputs ; j++ ){
            for(unsigned int dim =0; dim < proDim; dim++){
                currentMatrix[dim][0] = this->featureMatrix[dim][j];
            }

            dist = this->SurfaceDistance ( currentMatrix, reconstruction, proDim ) ;

            if ( dist < minDist ){
                minJ = j ;
                minDist = dist ;
            }
        }

        specificity += minDist ;
        sampleErrors[i] = minDist ;
    }

    // compute mean error
    specificity /= N ;

    // compute sample standard deviation
    double stdDev = 0 ;
    for ( int i = 0 ; i < N ; i++ )
    {
        sampleErrors[i] -= specificity ;
        stdDev += sampleErrors[i] * sampleErrors[i] ;
    }
    stdDev = sqrt ( stdDev / N ) ;

    // compute standard error of specificity
    stdError = stdDev / sqrt ( N ) ;

    return specificity ;
  }

  template<class TOutputMeshType>
  void srepCorrespondenceEvaluator<TOutputMeshType>
  ::SetNumberOfInputs( unsigned int number )
  {
    if ( number < 0 )
    {
      itkExceptionMacro( "Invalid number of inputs for correspondence evaluation." );
      return ;
    }
    this->m_NumberOfInputs = number ;
    this->m_Meshes.resize ( number ) ;
  }

  template<class TOutputMeshType>
  void srepCorrespondenceEvaluator<TOutputMeshType>
  ::SetInput( MatrixType featureMatrix ) {
    this->featureMatrix = featureMatrix;

      /*std::cout<<"-----------------------this->featureMatrix: "<<std::endl;
      for(unsigned int i = 0; i < this->featureMatrix.rows(); i++){
          for(unsigned int j = 0; j < this->featureMatrix.columns(); j++){
              std::cout<<this->featureMatrix[i][j]<<" ";
          }
          std::cout<<std::endl;
      }*/


    this->BringInputToOrigin() ;
  }

  // center the feature matrix to the origin
  template<class TOutputMeshType>
  void srepCorrespondenceEvaluator<TOutputMeshType>
  ::BringInputToOrigin() {
    double center = 0;

    // compute the feature mean (colum matrix, each row is the mean of a feature)
    for(unsigned int i = 0; i < this->featureMatrix.rows(); i++){
        for(unsigned int j = 0; j < this->featureMatrix.columns(); j++){
            center += this->featureMatrix[i][j];
        }

        center /= this->featureMatrix.columns();
        //std::cout<<"--------------center is: "<<center<<std::endl;

        // center each dimension
        for(unsigned int j = 0; j < this->featureMatrix.columns(); j++){
            this->featureMatrix[i][j] -= center;
        }
    }
  }

  template<class TOutputMeshType>
  double srepCorrespondenceEvaluator<TOutputMeshType>
  ::ComputeGeneralizationError ( unsigned int i, unsigned int M )  {
    if ( ( i < 0 ) || ( i >= this->m_NumberOfInputs ) )  {
      itkExceptionMacro ("Invalid input id for correspondence evalution.") ;
      return -1 ;
    }

    //////////////////////////////////////////////////////////////////////
    // 1. BUILD THE MODEL FROM THE TRAINING SET, WITH Xi REMOVED        //
    //////////////////////////////////////////////////////////////////////

    int proDim = this->featureMatrix.rows();

    MatrixType X(proDim, 1);
    for(unsigned int dim =0; dim < proDim; dim++){
        X[dim][0] = this->featureMatrix[dim][i];
    }

    MatrixType rowMean(proDim, 1);
    rowMean.fill(0);
    MatrixType pcaMatrix(proDim, this->m_NumberOfInputs-1);

    for (unsigned int row = 0; row < proDim; row++) {
      // compute mean of the training set, with Xi removed
      for ( unsigned int sampleIdx = 0 ; sampleIdx < this->m_NumberOfInputs ; sampleIdx++ ) {
        if ( sampleIdx == i )  {
          continue ;
        }

        rowMean[row][0] += this->featureMatrix[row][sampleIdx];
      }
      rowMean[row][0] /= ( this->m_NumberOfInputs - 1 );

      // build pca matrix, the matrix where we are going to perform pca to find principal modes of variation of the model      
      unsigned int skip = 0 ;
      for ( unsigned int sampleIdx = 0 ; sampleIdx < this->m_NumberOfInputs ; sampleIdx++ ) {
        if ( sampleIdx == i ) {
          skip = 1 ;
          continue ;
        }

        (pcaMatrix)[row][sampleIdx-skip] = this->featureMatrix[row][sampleIdx] - rowMean[row][0];
      }
    }

    std::vector<double> EigenValues;
    MatrixType PrincipalDirections;

    EigenValues.resize ( this->m_NumberOfInputs ) ;

    vnl_svd<double> svd( pcaMatrix );
    for (unsigned int m=0; m<this->m_NumberOfInputs-1; m++)
    {
      double sv = svd.W( m );
      EigenValues[m] = sv ;
    }
    PrincipalDirections = svd.U();
    //std::cout<<"--------------------------PrincipalDirections---:"<<PrincipalDirections.rows()<<"  "<<PrincipalDirections.columns()<<std::endl;
    /*for(unsigned int r =0; r<PrincipalDirections.rows(); r++){
        for(unsigned int c =0; c<PrincipalDirections.columns(); c++){
            std::cout<<PrincipalDirections[r][c]<<"  ";
        }
        std::cout<<std::endl;
    }*/

    //////////////////////////////////////////////////////////////////////
    // 2. ESTIMATE THE MODEL PARAMETERS FOR Xi                          //
    //////////////////////////////////////////////////////////////////////

    MatrixType b(M,1) ;
    b.fill(0);

    for ( unsigned int m = 0 ; m < M ; m++ ) {
        for ( unsigned int dim = 0 ; dim < proDim ; dim++ )  {
            b[m][0] += PrincipalDirections[dim][m] * ( X[dim][0] - rowMean[dim][0] ) ;
        }
    }


    //////////////////////////////////////////////////////////////////////
    // 3. RECONSTRUCT Xi USING M SHAPE PARAMETERS                       //
    //////////////////////////////////////////////////////////////////////

    MatrixType reconstruction(proDim, 1);

    for ( unsigned int dim = 0 ; dim < proDim ; dim++ )  {
        reconstruction[dim][0] = rowMean[dim][0];

        for ( unsigned int m = 0 ; m < M ; m++ ){
            reconstruction[dim][0] += PrincipalDirections[dim][m]*b[m][0];
        }
    }

    //////////////////////////////////////////////////////////////////////
    // 4. CALCULATE THE APPROXIMATION ERROR                             //
    //////////////////////////////////////////////////////////////////////

    return this->SurfaceDistance( X, reconstruction, proDim ) ;

  }

  template<class TOutputMeshType>
  double srepCorrespondenceEvaluator<TOutputMeshType>
  ::SurfaceDistance (MatrixType from, MatrixType to, unsigned int nDim){
      MatrixType distanceVector = from - to;
      double error = 0 ;

      for(unsigned int dim = 0; dim < nDim; dim++){
          error += sqrt(distanceVector[dim][0] * distanceVector[dim][0]);
      }

      return error / nDim;

  }

  template<class TOutputMeshType>
  void srepCorrespondenceEvaluator<TOutputMeshType>
  ::WriteSurfaceForDebug (MatrixType *surface, unsigned int nPoints, std::string filename)
  {
      std::ofstream out ;
    out.open (filename.c_str()) ;
    for ( unsigned int pt = 0 ; pt < nPoints ; pt++ )
      {
      out << (*surface)[0][pt] << " " ;
      out << (*surface)[1][pt] << " " ;
      out << (*surface)[2][pt] << std::endl  ;
      }
      out.close () ;
  }

  template<class TOutputMeshType>
  double srepCorrespondenceEvaluator<TOutputMeshType>
  ::GenerateUniformRandomNumber ()
  {
    // generate a random number uniformly distributed between [-1..1]
    double uniform = ( rand () / (double) RAND_MAX ) * 2 - 1.0 ;
    return uniform ;
  }

  template<class TOutputMeshType>
  double srepCorrespondenceEvaluator<TOutputMeshType>
  ::GenerateGaussianRandomNumber ()
  {
    // generate a random number normally distributed between [-1..1]

    // from http://www.taygeta.com/random/gaussian.html
    // Algorithm by Dr. Everett (Skip) Carter, Jr.

    float x1, x2, w, y1 ;

    do {
        x1 = this->GenerateUniformRandomNumber() ;
        x2 = this->GenerateUniformRandomNumber() ;
        w = x1 * x1 + x2 * x2;
    } while ( w >= 1.0 );

    w = sqrt( (-2.0 * log( w ) ) / w );
    y1 = x1 * w;

    return y1 ;
  }

}

#endif
