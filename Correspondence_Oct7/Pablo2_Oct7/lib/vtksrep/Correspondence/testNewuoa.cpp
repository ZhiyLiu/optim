#include <iostream>
#include "M3DNewuoaOptimizer.h"
//#include "itkImage.h"
//#include "itkImageFileReader.h"
#include "ImageDistanceMap.h"
#include "RAWImageFile.h"
#include "M3DObjectFile.h"
#include "WorldSystem.h"

int main(int argc, char** argv)
{
    if(argc < 3)
    {
        std::cout << "[Usage:]" << std::endl;
        return -1;
    }

    char* modelName = argv[1];
    char* imagePath = argv[2];

/*
    const unsigned int Dimension = 2;
    typedef unsigned char                      PixelType;
    typedef itk::Image< PixelType, Dimension > ImageType;
    typedef itk::ImageFileReader< ImageType >  ReaderType;
    typedef itk::Image<float, 3>               RealImage;
    ReaderType::Pointer reader = ReaderType::New();
    reader->SetFileName( argv[1] );
    reader->Update();

    ImageType::Pointer image = reader->GetOutput();
    itk::SignedDanielssonDistanceMapImageFilter<ImageType,RealImage>::Pointer ddmFilter =
    itk::SignedDanielssonDistanceMapImageFilter<ImageType,RealImage>::New();
    ddmFilter->SetInput(image);
    ddmFilter->Update();
*/
    // read image 
    RAWImageFile imageFile;
    Image3D* image3D = imageFile.read(imagePath,true);
    ImageDistanceMap* binaryDistanceMap = new ImageDistanceMap(image3D);
    binaryDistanceMap->fromImage3D(image3D);
    M3DNewuoaOptimizer optimizer;

    // read model
    M3DObjectFile objectFile;
    int primitiveId = 0;
    M3DObject* object = objectFile.read(modelName, primitiveId,
                                          new SimilarityTransform3D, new WorldSystem);

    optimizer.setImage(binaryDistanceMap);
    optimizer.setObject(object);

    M3DObject* modelAfterOptimizer = new M3DObject;
    optimizer.perform(modelAfterOptimizer);

    objectFile.write("out.m3d", *modelAfterOptimizer);
    delete modelAfterOptimizer;
    return 0;
}