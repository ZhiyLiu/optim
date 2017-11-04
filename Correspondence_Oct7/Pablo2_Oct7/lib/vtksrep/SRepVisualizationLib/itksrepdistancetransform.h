#ifndef ITKSrepDistanceTransform_H
#define ITKSrepDistanceTransform_H

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

template<  class TDistanceImage    >
class ITK_EXPORT SrepDistanceTransform
        : public ImageSource< TDistanceImage >
{
public:

    typedef SrepDistanceTransform         Self;
    typedef ImageSource<TDistanceImage>  Superclass;
    typedef SmartPointer<Self>         Pointer;
    typedef SmartPointer<const Self>   ConstPointer;


    /** Method for creation through the object factory. */
    itkNewMacro(Self);

    /** Run-time type information (and related methods). */
    itkTypeMacro(SrepDistanceTransform,ImageSource);


    /** Some convenient typedefs. */
    typedef TDistanceImage                        DistanceImageType;
    typedef typename DistanceImageType::Pointer        DistanceImagePointerType;
    typedef typename DistanceImageType::ConstPointer   DistanceImageConstPointerType;
    typedef typename DistanceImageType::RegionType     DistanceImageRegionType;
    typedef typename DistanceImageType::PixelType      DistanceImagePixelType;
    typedef typename DistanceImageType::IndexType      DistanceImageIndexType;


    typedef itk::ImageRegionIteratorWithIndex< DistanceImageType > DistanceRegionIteratorType;

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
    virtual void ThreadedGenerateData(const DistanceImageRegionType& outputRegionForThread, int threadId );


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
    virtual int SplitRequestedRegion(int i, int num, DistanceImageRegionType& splitRegion);

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


protected:

    SrepDistanceTransform();
    ~SrepDistanceTransform();


    private:

    double m_Bounds[6];

    int m_SizeX;
    int m_SizeY;
    int m_SizeZ;

    vtkSRep* m_Input;

    vtkSmartPointer< vtkSRepInterpolateMedialSpokesHermite > m_Medialspokesinterpolator;
    vtkSmartPointer< vtkSRepInterpolateMedialSheet > m_Medialsheetinterpolator;    
    vtkSmartPointer<vtkSRepInterpolateCrestSpokesQuartic> m_Interpolatecrestspokes;

};

}

#ifndef ITK_MANUAL_INSTANTIATION
#include "itksrepdistancetransform.txx"
#endif

#endif // ITKSrepDistanceTransform_H







