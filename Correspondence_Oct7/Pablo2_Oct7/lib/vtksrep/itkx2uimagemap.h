#ifndef ITKX2UIMAGEMAP_H
#define ITKX2UIMAGEMAP_H

#include "itkImageSource.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "vector"
//#include "itksampleannconnector.h"
#include "itkListSample.h"
#include "vtkPolyData.h"
#include "vtkRenderer.h"

#include "vtksrep.h"
#include "vtksrepinterpolatemedialspokeshermite.h"
#include "vtksrepinterpolatemedialsheet.h"
#include "vtksrepinterpolatemedialcrestcurve.h"
#include "vtksrepinterpolatecrestspokesquartic.h"

#include <itkVectorInterpolateImageFunction.h>

#ifndef PI
#define PI 3.14159265
#endif

using namespace std;

namespace itk{

template<   class TX2UImage,
            class TU2XImage
            //class SampleANNType = Statistics::SampleANNConnector< Statistics::ListSample< Vector< double, 3 > > >
        >
class ITK_EXPORT X2UImageMap
        : public ImageSource< TX2UImage >
{
public:

    typedef X2UImageMap         Self;
    typedef ImageSource<TX2UImage>  Superclass;
    typedef SmartPointer<Self>         Pointer;
    typedef SmartPointer<const Self>   ConstPointer;


    /** Method for creation through the object factory. */
    itkNewMacro(Self);

    /** Run-time type information (and related methods). */
    itkTypeMacro(X2UImageMap,ImageSource);


    /** Some convenient typedefs. */
    typedef TX2UImage                             X2UImageType;
    typedef typename X2UImageType::Pointer        X2UImagePointerType;
    typedef typename X2UImageType::ConstPointer   X2UImageConstPointerType;
    typedef typename X2UImageType::RegionType     X2UImageRegionType;
    typedef typename X2UImageType::PixelType      X2UImagePixelType;
    typedef typename X2UImageType::IndexType      X2UImageIndexType;


    typedef itk::ImageRegionIteratorWithIndex< X2UImageType > X2URegionIteratorType;

    /** If an imaging filter can be implemented as a multithreaded
     * algorithm, the filter will provide an implementation of
     * ThreadedGenerateData().  This superclass will automatically split
     * the output image into a number of pieces, spawn multiple threads,
     * and call ThreadedGenerateData() in each thread. Prior to spawning
     * threads, the BeforeThreadedGenerateData() method is called. After
     * all the threads have completed, the AfterThreadedGenerateData()
     * method is called. If an image processing filter cannot support
     * threading, that filter should provide an implementation of the
     * GenerateData() method instead of providing an implementation of
     * ThreadedGenerateData().  If a filter provides a GenerateData()
     * method as its implementation, then the filter is responsible for
     * allocating the output data.  If a filter provides a
     * ThreadedGenerateData() method as its implementation, then the
     * output memory will allocated automatically by this superclass.
     * The ThreadedGenerateData() method should only produce the output
     * specified by "outputThreadRegion"
     * parameter. ThreadedGenerateData() cannot write to any other
     * portion of the output image (as this is responsibility of a
     * different thread).
     *
     * \sa GenerateData(), SplitRequestedRegion() */
    virtual void ThreadedGenerateData(const X2UImageRegionType& outputRegionForThread, int threadId );


    /** The GenerateData method normally allocates the buffers for all of the
     * outputs of a filter. Some filters may want to override this default
     * behavior. For example, a filter may have multiple outputs with
     * varying resolution. Or a filter may want to process data in place by
     * grafting its input to its output. */
    virtual void AllocateOutputs();

    /** If an imaging filter needs to perform processing after the buffer
     * has been allocated but before threads are spawned, the filter can
     * can provide an implementation for BeforeThreadedGenerateData(). The
     * execution flow in the default GenerateData() method will be:
     *      1) Allocate the output buffer
     *      2) Call BeforeThreadedGenerateData()
     *      3) Spawn threads, calling ThreadedGenerateData() in each thread.
     *      4) Call AfterThreadedGenerateData()
     * Note that this flow of control is only available if a filter provides
     * a ThreadedGenerateData() method and NOT a GenerateData() method. */
    virtual void BeforeThreadedGenerateData();

    /** If an imaging filter needs to perform processing after all
     * processing threads have completed, the filter can can provide an
     * implementation for AfterThreadedGenerateData(). The execution
     * flow in the default GenerateData() method will be:
     *      1) Allocate the output buffer
     *      2) Call BeforeThreadedGenerateData()
     *      3) Spawn threads, calling ThreadedGenerateData() in each thread.
     *      4) Call AfterThreadedGenerateData()
     * Note that this flow of control is only available if a filter provides
     * a ThreadedGenerateData() method and NOT a GenerateData() method. */
    virtual void AfterThreadedGenerateData();

    /** Split the output's RequestedRegion into "num" pieces, returning
     * region "i" as "splitRegion". This method is called "num" times. The
     * regions must not overlap. The method returns the number of pieces that
     * the routine is capable of splitting the output RequestedRegion,
     * i.e. return value is less than or equal to "num". */
    virtual int SplitRequestedRegion(int i, int num, X2UImageRegionType& splitRegion);

    /**
      Set the dimensions of the mrep
    */
    itkSetObjectMacro(Input, vtkSRep)
    itkGetMacro(Input, vtkSRep*)

    itkSetMacro(SizeX, int)
    itkGetMacro(SizeX,int)

    itkSetMacro(SizeY, int)
    itkGetMacro(SizeY,int)

    itkSetMacro(SizeZ, int)
    itkGetMacro(SizeZ,int)

    itkSetMacro(BinaryImage, bool)
    itkGetMacro(BinaryImage,bool)

    itkSetMacro(DistanceTransformImage, bool)
    itkGetMacro(DistanceTransformImage,bool)

    typedef TU2XImage                                  U2XImageType;
    typedef typename U2XImageType::Pointer             U2XImagePointerType;
    typedef typename U2XImageType::RegionType          U2XImageRegionType;
    typedef typename U2XImageType::PixelType           U2XImagePixelType;
    typedef typename U2XImageType::IndexType           U2XImageIndexType;

    typedef itk::ImageRegionIteratorWithIndex< U2XImageType > U2XRegionIteratorType;

    /**
      Get the U2XImage
    */
    itkSetMacro(U2XImage, U2XImagePointerType);
    itkGetMacro(U2XImage, U2XImagePointerType);

protected:

    X2UImageMap();
    ~X2UImageMap();


    private:

    double m_Bounds[6];

    int m_SizeX;
    int m_SizeY;
    int m_SizeZ;
    bool m_BinaryImage;
    bool m_DistanceTransformImage;

    U2XImagePointerType m_U2XImage;
    vtkSRep* m_Input;

    vtkSmartPointer< vtkSRepInterpolateMedialSpokesHermite > m_Medialspokesinterpolator;
    vtkSmartPointer< vtkSRepInterpolateMedialSheet > m_Medialsheetinterpolator;    
    vtkSmartPointer<vtkSRepInterpolateCrestSpokesQuartic> m_Interpolatecrestspokes;

};

}

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkx2uimagemap.txx"
#endif

#endif // ITKX2UIMAGEMAP_H







