
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImage.h"
#include "itkResampleImageFilter.h"
#include "itkThinPlateSplineKernelTransform.h"
#include "itkPointSet.h"
#include <fstream>


int main( int argc, char * argv[] ) {
    if( argc < 4 )
    {
        std::cerr << "Missing Parameters " << std::endl;
        std::cerr << "Usage: " << argv[0];
        std::cerr << " landmarksFile inputImage ";
        std::cerr << "DeformedImage " << std::endl;
        std::cerr << "deformationField" << std::endl;
        return EXIT_FAILURE;
    }
    const unsigned int ImageDimension = 3;
    typedef unsigned char PixelType;
    typedef itk::Image< PixelType, ImageDimension > InputImageType;
    typedef itk::ImageFileReader< InputImageType > ReaderType;
    typedef itk::ImageFileWriter< InputImageType > DeformedImageWriterType;
    typedef itk::Point< float, ImageDimension > FieldPointType;
    typedef itk::Vector< float, ImageDimension > FieldVectorType;
    typedef itk::Image< FieldVectorType, ImageDimension > DisplacementFieldType;
    typedef itk::ImageFileWriter< DisplacementFieldType > FieldWriterType;
    typedef double CoordinateRepType;
    typedef itk::ThinPlateSplineKernelTransform< CoordinateRepType,ImageDimension> TransformType;
    typedef itk::Point< CoordinateRepType, ImageDimension > PointType;
    typedef std::vector< PointType > PointArrayType;
    typedef TransformType::PointSetType PointSetType;
    typedef PointSetType::Pointer PointSetPointer;
    typedef PointSetType::PointIdentifier PointIdType;
    typedef itk::ResampleImageFilter< InputImageType, InputImageType > ResamplerType;
    typedef itk::LinearInterpolateImageFunction< InputImageType, double > InterpolatorType;
    ReaderType::Pointer reader = ReaderType::New();
    reader->SetFileName( argv[2] );
    try   {
        reader->Update();
    }
    catch( itk::ExceptionObject & excp )  {
        std::cerr << "Exception thrown " << std::endl;
        std::cerr << excp << std::endl;
        return EXIT_FAILURE;
    }

    // Landmarks correspondances may be associated with the SplineKernelTransforms
    // via Point Set containers. Let us define containers for the landmarks.
    // Define container for landmarks
    PointSetType::Pointer sourceLandMarks = PointSetType::New();
    PointSetType::Pointer targetLandMarks = PointSetType::New();
    PointType p1; PointType p2;
    PointSetType::PointsContainer::Pointer sourceLandMarkContainer = sourceLandMarks->GetPoints();
    PointSetType::PointsContainer::Pointer targetLandMarkContainer = targetLandMarks->GetPoints();

    PointIdType id = itk::NumericTraits< PointIdType >::Zero;
    // Read in the list of landmarks
    std::ifstream infile;
    infile.open( argv[1] );
    while (!infile.eof())
    {
        infile >> p1[0] >> p1[1] >> p1[2] >> p2[0] >> p2[1] >> p2[2];
        sourceLandMarkContainer->InsertElement( id, p1 );
        targetLandMarkContainer->InsertElement( id++, p2 );
    }
    infile.close();

    TransformType::Pointer tps = TransformType::New();
    tps->SetSourceLandmarks(sourceLandMarks);
    tps->SetTargetLandmarks(targetLandMarks);
    tps->ComputeWMatrix();

    // The image is then resampled to produce an output image as defined by the transform. Here we use a LinearInterpolator.
    // Set the resampler params
    InputImageType::ConstPointer inputImage = reader->GetOutput();
    ResamplerType::Pointer resampler = ResamplerType::New();
    InterpolatorType::Pointer interpolator = InterpolatorType::New();
    resampler->SetInterpolator( interpolator );
    InputImageType::SpacingType spacing = inputImage->GetSpacing();
    InputImageType::PointType origin = inputImage->GetOrigin();
    InputImageType::DirectionType direction = inputImage->GetDirection();
    InputImageType::RegionType region = inputImage->GetBufferedRegion();
    InputImageType::SizeType size = region.GetSize();

    resampler->SetOutputSpacing( spacing );
    resampler->SetOutputDirection( direction );
    resampler->SetOutputOrigin( origin );
    resampler->SetSize( size );
    resampler->SetTransform( tps );

    resampler->SetOutputStartIndex( region.GetIndex() );
    resampler->SetInput( reader->GetOutput() );

    //Set and write deformed image
    DeformedImageWriterType::Pointer deformedImageWriter = DeformedImageWriterType::New();
    deformedImageWriter->SetInput( resampler->GetOutput() );
    deformedImageWriter->SetFileName( argv[3] );
    try  {
        deformedImageWriter->Update();
    }
    catch( itk::ExceptionObject & excp )  {
        std::cerr << "Exception thrown " << std::endl;
        std::cerr << excp << std::endl;
        return EXIT_FAILURE;
    }

    // The deformation field is computed as the difference between the input and the deformed image by using an iterator.
    // Compute the deformation field
    DisplacementFieldType::Pointer field = DisplacementFieldType::New();
    field->SetRegions( region );
    field->SetOrigin( origin );
    field->SetSpacing( spacing );
    field->Allocate();
    typedef itk::ImageRegionIterator< DisplacementFieldType > FieldIterator;
    FieldIterator fi( field, region );
    fi.GoToBegin();
    TransformType::InputPointType point1;
    TransformType::OutputPointType point2;
    DisplacementFieldType::IndexType index;
    FieldVectorType displacement;
    while( ! fi.IsAtEnd() ) {
        index = fi.GetIndex();
        field->TransformIndexToPhysicalPoint( index, point1 );
        point2 = tps->TransformPoint( point1 );
        for ( unsigned int i = 0;i < ImageDimension;i++)  {
            displacement[i] = point2[i] - point1[i];
        }
        fi.Set( displacement );
        ++fi;
    }
    //Write computed deformation field
    FieldWriterType::Pointer fieldWriter = FieldWriterType::New();
    fieldWriter->SetFileName( argv[4] );
    fieldWriter->SetInput( field );
    try  {
        fieldWriter->Update();
    }
    catch( itk::ExceptionObject & excp )
    {
        std::cerr << "Exception thrown " << std::endl;
        std::cerr << excp << std::endl;
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
