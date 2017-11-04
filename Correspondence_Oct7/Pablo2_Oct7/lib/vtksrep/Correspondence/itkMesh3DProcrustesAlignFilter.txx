#ifndef _itkMesh3DProcrustesAlignFilter_txx
#define _itkMesh3DProcrustesAlignFilter_txx

#include "itkMesh3DProcrustesAlignFilter.h"

namespace itk
{

template <class TInputMesh, class TOutputMesh>
Mesh3DProcrustesAlignFilter<TInputMesh, TOutputMesh>
::Mesh3DProcrustesAlignFilter()
{
  m_Convergence = 0.0001;
  m_Mean = OutputMeshType::New();
  m_OldMean = OutputMeshType::New();

  m_UseScaling = true;
  m_UseInitialAverage = false;
  m_UseSingleIteration = false;
  m_UseNormalization = true;
  m_AlignRotation = true;
}

template <class TInputMesh, class TOutputMesh>
Mesh3DProcrustesAlignFilter<TInputMesh, TOutputMesh>
::~Mesh3DProcrustesAlignFilter()
{
}

template <class TInputMesh, class TOutputMesh>
typename itk::ProcessObject::DataObjectPointer
Mesh3DProcrustesAlignFilter<TInputMesh, TOutputMesh>
::MakeOutput(itk::ProcessObject::DataObjectPointerArraySizeType /* idx */)
{
  return static_cast<itk::ProcessObject::DataObjectPointer>(TOutputMesh::New().GetPointer() );
}

template <class TInputMesh, class TOutputMesh>
typename Mesh3DProcrustesAlignFilter<TInputMesh, TOutputMesh>::OutputMeshType
* Mesh3DProcrustesAlignFilter<TInputMesh, TOutputMesh>
::GetOutput(unsigned int idx)
  {
  return m_MeshTransform[idx]->GetOutput();
  // return static_cast<TOutputMesh*>(this->ProcessObject::GetOutput(idx));
  }

template <class TInputMesh, class TOutputMesh>
void
Mesh3DProcrustesAlignFilter<TInputMesh, TOutputMesh>
::SetNumberOfInputs( unsigned int num )
{
  this->ProcessObject::SetNumberOfInputs( num );
  this->ProcessObject::SetNumberOfRequiredInputs( num );
  this->ProcessObject::SetNumberOfOutputs( num );
  this->ProcessObject::SetNumberOfRequiredOutputs( num );
  m_MeshTransform.resize( num );
  m_Center.resize( num );
  for( unsigned int i = 0; i < num; i++ )
    {
    OutputMeshPointer output = static_cast<TOutputMesh *>(this->MakeOutput(i).GetPointer() );
    this->ProcessObject::SetNthOutput( i, output.GetPointer() );
    m_MeshTransform[i] = TransformMeshType::New();
    m_MeshTransform[i]->GraftOutput( this->GetOutput( i ) );
    }
}

template <class TInputMesh, class TOutputMesh>
void
Mesh3DProcrustesAlignFilter<TInputMesh, TOutputMesh>
::SetInput( unsigned int idx, InputMeshPointer mesh )
{
  this->ProcessObject::SetNthInput( idx, mesh );
  m_MeshTransform[idx]->SetInput( this->GetInput( idx ) );
}

template <class TInputMesh, class TOutputMesh>
void
Mesh3DProcrustesAlignFilter<TInputMesh, TOutputMesh>
::GenerateData()
{
  // find mesh centers and store them
  for( unsigned int i = 0; i < this->GetNumberOfInputs(); i++ )
    {
    m_Center[i] = GetMeshCenter( i );
    }
  // check if mean and oldMean have to be initialized
  InputMeshPointer mesh = this->GetInput( 0 );
  unsigned int     numPoints = mesh->GetNumberOfPoints();
  if( m_Mean->GetNumberOfPoints() != numPoints )
    {
    OutputPointsContainerPointer MeanPoints = m_Mean->GetPoints();
    OutputPointsContainerPointer OldMeanPoints =  m_OldMean->GetPoints();
    MeanPoints->Reserve(numPoints);
    OldMeanPoints->Reserve(numPoints);

    // std::cout << m_Mean->GetNumberOfPoints() << "," << numPoints << std::endl;
    }

  OutputPointsContainer *meanPoints = m_Mean->GetPoints();
  OutputPointsContainer *oldMeanPoints = m_OldMean->GetPoints();
  typename OutputPointsContainer::Iterator meanIt, oldMeanIt;
  if( m_UseInitialAverage )
    {
     // m_MeanCenter.Fill( 0 );
    // initialize mean shape to the actual mean shape
    //  InitialMean();
    CalculateMean();
    }
  else
    {
    // initialize mean shape to first input mesh by copying the mesh
    typename InputMeshType::PointsContainer::ConstIterator meshIt;
    meshIt = mesh->GetPoints()->Begin();
    for( meanIt = meanPoints->Begin(); meanIt != meanPoints->End(); ++meanIt )
      {
      for( int dim = 0; dim < 3; dim++ )
        {
        meanIt.Value()[dim] = meshIt.Value()[dim];
        }
      ++meshIt;
      }
    m_MeanCenter = m_Center[0];
    }

  // iteratively update mean
  OutputPointType zeroPnt;
  for( int dim = 0; dim < 3; dim++ )
    {
    zeroPnt[dim] = 0;
    }
  CoordRepType diff, squaredDiff;

  do
    {
    // find current transformations to match the mean
    for( unsigned int i = 0; i < this->GetNumberOfInputs(); i++ )
      {
      TransformPointer transform = GetProcrustesMatch( i, m_Mean, m_MeanCenter);
      m_MeshTransform[i]->SetTransform( transform );
      PrintTransform( transform);
      }
    // copy mean to oldMean and set mean to 0
    oldMeanIt = oldMeanPoints->Begin();
    for( meanIt = meanPoints->Begin(); meanIt != meanPoints->End(); ++meanIt )
      {
      oldMeanIt.Value() = meanIt.Value();
      meanIt.Value() = zeroPnt;
      ++oldMeanIt;
      }
    // calculate new mean
    CalculateMean();
    // calculate average point distance between old mean and new mean
    squaredDiff = 0;
    oldMeanIt = oldMeanPoints->Begin();
    for( meanIt = meanPoints->Begin(); meanIt != meanPoints->End(); ++meanIt )
      {
      for( int dim = 0; dim < 3; dim++ )
        {
        diff = meanIt.Value()[dim] - oldMeanIt.Value()[dim];
        squaredDiff += diff * diff;
        }
      ++oldMeanIt;
      }
    diff = sqrt( squaredDiff );
    }
  while( diff > m_Convergence && !m_UseSingleIteration );
}

template <class TInputMesh, class TOutputMesh>
typename Mesh3DProcrustesAlignFilter<TInputMesh, TOutputMesh>::TranslationType
Mesh3DProcrustesAlignFilter<TInputMesh, TOutputMesh>
::GetMeshCenter( unsigned int idx )
{
  TranslationType center;

  center.Fill( 0 );
  InputPointType pnt;
  typename InputMeshType::PointsContainer * meshPoints = this->GetInput( idx )->GetPoints();
  typename InputMeshType::PointsContainer::ConstIterator pntIt;
  for( pntIt = meshPoints->Begin(); pntIt != meshPoints->End(); ++pntIt )
    {
    for( int dim = 0; dim < 3; dim++ )
      {
      center[dim] += pntIt.Value()[dim];
      }
    }
  center /= (double)meshPoints->Size();
  return center;
}

template <class TInputMesh, class TOutputMesh>
typename Mesh3DProcrustesAlignFilter<TInputMesh, TOutputMesh>::TransformPointer
Mesh3DProcrustesAlignFilter<TInputMesh, TOutputMesh>
::GetProcrustesMatch( unsigned int idx, OutputMeshPointer targetMesh, TranslationType targetCenter)
{
  TransformPointer result = TransformType::New();
  // copy source mesh coordinates to source matrix, translating to zero if necessary
  MatrixType       source;
  InputMeshPointer sourceMesh = this->GetInput( idx );

  source.set_size( sourceMesh->GetNumberOfPoints(), 3 );
  typename InputMeshType::PointsContainer::ConstIterator inputIt;
  unsigned int i = 0;
  for( inputIt = sourceMesh->GetPoints()->Begin(); inputIt != sourceMesh->GetPoints()->End(); ++inputIt )
    {
    for( int dim = 0; dim < 3; dim++ )
      {
      source[i][dim] = inputIt.Value()[dim] - m_Center[idx][dim];
      }
    i++;
    }
  // copy target mesh coordinates to target matrix
  MatrixType target;
  target.set_size( 3, targetMesh->GetNumberOfPoints() );
  typename OutputMeshType::PointsContainer::ConstIterator outputIt;
  i = 0;
  for( outputIt = targetMesh->GetPoints()->Begin(); outputIt != targetMesh->GetPoints()->End(); ++outputIt )
    {
    for( int dim = 0; dim < 3; dim++ )
      {
      target[dim][i] = outputIt.Value()[dim] - targetCenter[dim]; //targetCenter first iteration not 0, from 2th iteration on always 0...
      }
    i++;
    }
  cout<<"----------------------targetCenter is: "<<targetCenter<<endl;
  // do procrustes matching
  MatrixType            x1 = target * source / (target.fro_norm() * source.fro_norm() );
  vnl_svd<CoordRepType> svd( x1 );
  MatrixType            postTrans = svd.V() * svd.U().transpose();
  MatrixType            x2 = target * source * postTrans;
  CoordRepType          x2Trace = 0;
  for( unsigned int i_ = 0; i_ < x2.rows(); i_++ )
    {
    x2Trace += x2[i_][i_];
    }
  MatrixType   x3 = source.transpose() * source;
  CoordRepType x3Trace = 0;
  for( unsigned int i_ = 0; i_ < x3.rows(); i_++ )
    {
    x3Trace += x3[i_][i_];
    }
  CoordRepType scale = x2Trace / x3Trace;

  // set up transformation
  typename TransformType::InputPointType center;
  for( int dim = 0; dim < 3; dim++ ) {
      center[dim] = m_Center[idx][dim];
  }
  result->SetCenter( center ); //set the rotation center.

  //result->SetTranslation(targetCenter - m_Center[idx]); // This is not move to origin, accumulate deviation when large iteration.
  result->SetTranslation(- m_Center[idx]); // move the source center to the origin.

  if( m_AlignRotation )
    {

    typename TransformType::MatrixType rotMatrix;
    for( int r = 0; r < 3; r++ )
      {
      for( int c = 0; c < 3; c++ )
        {
        rotMatrix[r][c] = postTrans[c][r];
        }
      }
    result->SetMatrix( rotMatrix );
    }

  if( m_UseScaling )   {
    //result->Scale( scale); //The scaling is precomposed(A=RTS). Note that the scaling is applied centered at the origin.
    result->Scale( scale, 1); //The scaling is postcomposed(A=SRT). Note that the scaling is applied centered at the origin.
    }

  // Dont need. since targetCenter is already at the origin. If needed, not here, but put after update the transform(in the CalculateMean method)....
  /*typename TransformType::OutputVectorType targetOffset;
  for( int dim = 0; dim < 3; dim++ )
    {
    targetOffset[dim] = targetCenter[dim];
    cout<<"======================targetOffset is: "<<targetOffset[dim]<<endl;
    }

  result->Translate( targetOffset );*/

  return result;
}

template <class TInputMesh, class TOutputMesh>
void
Mesh3DProcrustesAlignFilter<TInputMesh, TOutputMesh>
::CalculateMean()
{
  typename OutputMeshType::PointsContainer * meanPoints = m_Mean->GetPoints();
  typename OutputMeshType::PointsContainer::Iterator meanIt;
  for( unsigned int i = 0; i < this->GetNumberOfInputs(); i++ )
    {
    m_MeshTransform[i]->Update();
    OutputMeshPointer transformedMesh = m_MeshTransform[i]->GetOutput();
    typename OutputMeshType::PointsContainer * transformedPoints = transformedMesh->GetPoints();
    typename OutputMeshType::PointsContainer::Iterator pntIt; //ConstIterator pntIt;//read only iterator

    // After transform, move the aligned points back to targetcenter
    /*for( pntIt = transformedPoints->Begin(); pntIt != transformedPoints->End(); ++pntIt ){
        for( int dim = 0; dim < 3; dim++ ) {
             pntIt.Value()[dim] +=  m_MeanCenter[dim];
        }
    }*/

    pntIt = transformedPoints->Begin();
    for( meanIt = meanPoints->Begin(); meanIt != meanPoints->End(); ++meanIt )
      {
      for( int dim = 0; dim < 3; dim++ )
        {
        meanIt.Value()[dim] += pntIt.Value()[dim];//+ m_MeanCenter[dim];
        }
      ++pntIt;
      }
     //cout<<"----------------------m_MeanCenter is: "<<m_MeanCenter<<endl;
    }

  for( meanIt = meanPoints->Begin(); meanIt != meanPoints->End(); ++meanIt )
    {
    for( int dim = 0; dim < 3; dim++ )
      {
      meanIt.Value()[dim] /= this->GetNumberOfInputs();
      }
    }
  if( m_UseNormalization ) // if not normalize, the aligned srep shrink to very small, why?
    {
    // scale mean to get a norm of 1
    double squaredNorm = 0;
    for( meanIt = meanPoints->Begin(); meanIt != meanPoints->End(); ++meanIt )
      {
      squaredNorm += meanIt.Value().GetVectorFromOrigin().GetSquaredNorm();
      }
    double scaleFactor = 1.0 / sqrt( squaredNorm );
    for( meanIt = meanPoints->Begin(); meanIt != meanPoints->End(); ++meanIt )
      {
      for( int dim = 0; dim < 3; dim++ )
        {
        meanIt.Value()[dim] *= scaleFactor;
        }
      }
    }

  m_MeanCenter.Fill( 0 );
  for( meanIt = meanPoints->Begin(); meanIt != meanPoints->End(); ++meanIt )
    {
    for( int dim = 0; dim < 3; dim++ )
      {
      m_MeanCenter[dim] += meanIt.Value()[dim];
      }
    }
  m_MeanCenter /= (double)meanPoints->Size();
  cout<<"========================================================m_MeanCenter is: "<<m_MeanCenter<<endl;
}



/*
template <class TInputMesh, class TOutputMesh>
void
Mesh3DProcrustesAlignFilter<TInputMesh, TOutputMesh>
::InitialMean()
{
  typename OutputMeshType::PointsContainer * meanPoints = m_Mean->GetPoints();
  typename OutputMeshType::PointsContainer::Iterator meanIt;

    // iteratively update mean
    OutputPointType zeroPnt;
    for( int dim = 0; dim < 3; dim++ )
      {
      zeroPnt[dim] = 0;
      }
    for( meanIt = meanPoints->Begin(); meanIt != meanPoints->End(); ++meanIt )  {
        // Initialize to 0
        meanIt.Value() = zeroPnt;
    }

  for( unsigned int i = 0; i < this->GetNumberOfInputs(); i++ )  {
      InputMeshPointer sourceMesh = this->GetInput( i );
      typename InputMeshType::PointsContainer::ConstIterator pntIt;
      cout<<"-----------------------the mesh used to compute initial mean is, before center: "<<endl;
      for( pntIt = sourceMesh->GetPoints()->Begin(); pntIt != sourceMesh->GetPoints()->End(); ++pntIt )
        {
          cout<<pntIt.Value()<<endl;
        }

    // Center the mesh
    double sum_x=0;
    double sum_y=0;
    double sum_z=0;
    int pointsNum =0;
    for( pntIt = sourceMesh->GetPoints()->Begin(); pntIt != sourceMesh->GetPoints()->End(); ++pntIt ) {
        sum_x += pntIt.Value()[0];
        sum_y += pntIt.Value()[1];
        sum_z += pntIt.Value()[2];
        pointsNum++;
    }
    double c_x=sum_x/pointsNum;
    double c_y=sum_y/pointsNum;
    double c_z=sum_z/pointsNum;
    //cout<<"=====================c_x,c_y,c_z  is: "<<c_x<<"  "<<c_y<<"  "<<c_z<<endl;

    pntIt = sourceMesh->GetPoints()->Begin();
    cout<<"-----------------------the mesh used to compute initial mean, after centroid: "<<endl;
    for( meanIt = meanPoints->Begin(); meanIt != meanPoints->End(); ++meanIt )  {
        meanIt.Value()[0] += pntIt.Value()[0] - c_x; //why no priority for the /
        meanIt.Value()[1] += pntIt.Value()[1] - c_y;
        meanIt.Value()[2] += pntIt.Value()[2] - c_z;

        cout << pntIt.Value()[0] - c_x<< "  "<< pntIt.Value()[1] - c_y<<"  "<<pntIt.Value()[2] - c_z<<endl;
        ++pntIt;
    }
  }

  cout<<"=================================the initial mean is: "<<endl;
  for( meanIt = meanPoints->Begin(); meanIt != meanPoints->End(); ++meanIt )
    {
    for( int dim = 0; dim < 3; dim++ )
      {
      meanIt.Value()[dim] /= this->GetNumberOfInputs();
      }
    cout<<meanIt.Value()<<endl;
    }

  m_MeanCenter.Fill( 0 );
  for( int dim = 0; dim < 3; dim++ )
    {
    m_MeanCenter[dim] =0;
    }

  for( meanIt = meanPoints->Begin(); meanIt != meanPoints->End(); ++meanIt )
    {
    for( int dim = 0; dim < 3; dim++ )
      {
      m_MeanCenter[dim] += meanIt.Value()[dim];
      }
    }
  m_MeanCenter /= (double)meanPoints->Size();
}

*/




}

#endif
