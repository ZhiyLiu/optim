#ifndef ITKVectorDiscreteGaussianImageFilter_H
#define ITKVectorDiscreteGaussianImageFilter_H



#include "itkImageRegionIterator.h"
#include "itkImageToImageFilter.h"

namespace itk{

template<class TInputImage,class TOutputImage>
class VectorDiscreteGaussianImageFilter :
public ImageToImageFilter<TInputImage,TOutputImage>
{
public:
    /** Standard class typedefs. */
    typedef VectorDiscreteGaussianImageFilter             Self;
    typedef ImageToImageFilter<TInputImage,TOutputImage>      Superclass;
    typedef SmartPointer<Self>                                Pointer;
    typedef SmartPointer<const Self>                          ConstPointer;


    /** Pixel Type of the input image */
    typedef TInputImage                                       InputImageType;
    typedef TOutputImage                                      OutputImageType;
    typedef typename TInputImage::PixelType                   PixelType;

    /** Method for creation through the object factory. */
    itkNewMacro(Self);

    /** Run-time type information (and related methods). */
    itkTypeMacro(VectorDiscreteGaussianImageFilter, ImageToImageFilter);


    /** Some additional typedefs.  */
    typedef typename InputImageType::Pointer      InputImagePointer;
    typedef typename InputImageType::RegionType   InputImageRegionType;
    typedef typename InputImageType::PixelType    InputImagePixelType;

    /** Some additional typedefs.  */
    typedef typename OutputImageType::Pointer     OutputImagePointer;
    typedef typename OutputImageType::RegionType  OutputImageRegionType;
    typedef typename OutputImageType::PixelType   OutputImagePixelType;



    typedef itk::Image< double, InputImageType::ImageDimension > DiscreteImageType;
    typedef typename DiscreteImageType::Pointer  DiscreteImagePointerType;

    /** Iterators def**/
    typedef itk::ImageRegionIterator< InputImageType > RegionIteratorType;

    typedef itk::ImageRegionIterator< DiscreteImageType > RegionIteratorDiscreteType;


    /** Set the texel dimension the default is 3 for RGB textures
     * NumericTraits<PixelType>::Zero. */
    itkSetMacro(Variance,double);

    /** Get the texel dimension */
    itkGetConstMacro(Variance,double);

protected:
    VectorDiscreteGaussianImageFilter();
    ~VectorDiscreteGaussianImageFilter();
    virtual void GenerateData();

private:
    double m_Variance;
};

}


#ifndef ITK_MANUAL_INSTANTIATION
#include "itkvectordiscretegaussianimagefilter.txx"
#endif

#endif // ITKVectorDiscreteGaussianImageFilter_H
