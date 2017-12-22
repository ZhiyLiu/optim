#include <iostream>
#include "M3DNewuoaOptimizer.h"
#include "itkImage.h"
#include "itkImageFileReader.h"

int main(int argc, char** argv)
{
    if(argc < 2)
    {
        std::cout << "[Usage:]" << std::endl;
        return -1;
    }

    std::string modelName = argv[0];
    std::string imagePath = argv[1];

    const unsigned int Dimension = 2;
    typedef unsigned char                      PixelType;
    typedef itk::Image< PixelType, Dimension > ImageType;
    typedef itk::ImageFileReader< ImageType >  ReaderType;
    ReaderType::Pointer reader = ReaderType::New();
    reader->SetFileName( argv[1] );
    reader->Update();

    ImageType::Pointer image = reader->GetOutput();
    /*M3DNewuoaOptimizer optimizer;
    optimizer.setTargetImage();
    optimizer.setSrepModel();
    optimizer.perform();*/
    return 0;
}