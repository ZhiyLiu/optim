#ifndef ITKSrepDistanceTransform_TXX
#define ITKSrepDistanceTransform_TXX

#include "itkImage.h"
#include "itkRGBToLuminanceImageFilter.h"
#include "vtkGenericCell.h"



namespace itk{

template< class TDistanceImage >
SrepDistanceTransform< TDistanceImage >::SrepDistanceTransform()
{

    m_SizeX = 512;
    m_SizeY = 512;
    m_SizeZ = 512;

}

template< class TDistanceImage >
SrepDistanceTransform< TDistanceImage >::~SrepDistanceTransform()
{
}

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
template< class TDistanceImage >
void SrepDistanceTransform< TDistanceImage >::ThreadedGenerateData(const DistanceImageRegionType& outputRegionForThread, int threadId ){

    if(this->GetDebug())
        cout<<"this is thread "<<threadId<<endl;


    //typename DistanceImageType::SpacingType spacing = this->GetOutput()->GetSpacing();
    //double xspc = spacing[0];
    //double yspc = spacing[1];
    //double zspc = spacing[2];
    //double minspc = min(min(xspc, yspc), zspc)/8.0;
    double minspc = 1/pow(2,5);//min(min(xspc, yspc), zspc)/2.0;
    //const double* bounds = this->m_Bounds;

    DistanceRegionIteratorType outit(this->GetOutput(), this->GetOutput()->GetLargestPossibleRegion());

    vtkSRep* srepmedial = m_Medialspokesinterpolator->GetSRepOutput();

    int nthreads = this->GetNumberOfThreads();
    int itersize = ceil((double)srepmedial->GetNumberOfPoints() / nthreads);

    int startpos = itersize * threadId;
    int endpos = startpos + itersize;

    //PlanesOperations pop;
    if(this->GetDebug()){
        cout<<"this is thread = "<<threadId;
    }    

    //if(threadId == 0)
    {
        for(int i = startpos; i < endpos && i < (int)m_Input->GetNumberOfCells(); i++){
        //for(int i = 0; i < (int)m_Input->GetNumberOfCells(); i++){            

            for(double u = 0; u <= 1; u += minspc){
                for(double v = 0; v <= 1; v += minspc){

                    vtkSRep::VNLType position = m_Medialsheetinterpolator->GetInterpolatedPoint(i, u, v);

                    for(unsigned side = 0; side < 2; side++){
                        vtkSRep::VNLType spoke = m_Medialspokesinterpolator->GetInterpolatedSpoke(i, side, u, v);
                        for(double t = 0; t < 2; t+= minspc){

                            vtkSRep::VNLType currentpoint = spoke * t + position;

                            DistanceImagePixelType outpixel;
                            outpixel = 1 - t;


                            typedef typename DistanceImageType::PointType PointType;
                            PointType point;
                            point[0] = currentpoint[0];
                            point[1] = currentpoint[1];
                            point[2] = currentpoint[2];

                            DistanceImageIndexType index;

                            if(this->GetOutput()->TransformPhysicalPointToIndex(point, index)){
                                outit.SetIndex(index);
                                if(outit.Get() == -1 || t == 0){
                                    outit.Set(outpixel);
                                }
                            }else{
                                t = 2;
                            }

                            //cout<<"scalarpoint "<<(int)scalarpoint[0] <<" "<<(int)scalarpoint[1] <<" "<<(int)scalarpoint[2] <<" "<<endl;
                        }
                    }
                }
            }
        }
    }


    //if(threadId == 0)
    {
        /*vtkSRep::VectorIdsType crestids = m_Input->GetCrestMedialAtomsIds();

        int itersizecrest = ceil((double)crestids.size() / nthreads);

        int startposcrest = itersizecrest * threadId;
        int endposcrest = startposcrest + itersizecrest;

        vtkSmartPointer<vtkSRepInterpolateCrestSpokesQuartic> interpolatecrestspokes = vtkSmartPointer<vtkSRepInterpolateCrestSpokesQuartic>::New();
        interpolatecrestspokes->SetSpokeType(vtkSRep::CREST_SPOKE);
        interpolatecrestspokes->SetInput(m_Input);
        interpolatecrestspokes->Update();


        for(int i = startposcrest; i < endposcrest && i < crestids.size(); i++){
        //for(int i = 0; i < crestids.size(); i++){

            for(double u = 0; u <= 1; u += minspc){
                vtkSRep::VNLType position = interpolatecrestspokes->GetInterpolatedPoint(i, u);
                for(double v = 0; v <= 1; v += minspc){

                    vtkSRep::VNLType spoke;
                    if(v == 0){
                        spoke = interpolatecrestspokes->GetInterpolatedSpoke(i, u, v);
                    }else{
                        spoke = interpolatecrestspokes->GetInterpolatedSpoke(v);
                    }

                    for(double t = 0; t <= 2; t+= minspc){

                        vtkSRep::VNLType currentpoint = spoke * t + position;

                        DistanceImagePixelType outpixel;
                        outpixel = 1 - t;

                        typedef typename DistanceImageType::PointType PointType;
                        PointType point;
                        point[0] = currentpoint[0];
                        point[1] = currentpoint[1];
                        point[2] = currentpoint[2];

                        DistanceImageIndexType index;

                        if(this->GetOutput()->TransformPhysicalPointToIndex(point, index)){
                            outit.SetIndex(index);
                            if(outit.Get() == -1 || t == 0){
                                outit.Set(outpixel);
                            }
                        }else{
                            t = 2;
                        }
                    }
                }
            }
        }*/
    }
}


/** The GenerateData method normally allocates the buffers for all of the
 * outputs of a filter. Some filters may want to override this default
 * behavior. For example, a filter may have multiple outputs with
 * varying resolution. Or a filter may want to process data in place by
 * grafting its input to its output. */
template< class TDistanceImage >
void SrepDistanceTransform< TDistanceImage >::AllocateOutputs(){

    if(m_Input){        

        m_Medialsheetinterpolator = vtkSmartPointer< vtkSRepInterpolateMedialSheet >::New();
        m_Medialsheetinterpolator->SetInput(m_Input);
        m_Medialsheetinterpolator->Update();

        m_Medialspokesinterpolator = vtkSmartPointer< vtkSRepInterpolateMedialSpokesHermite >::New();
        m_Medialspokesinterpolator->SetInput(m_Input);        
        m_Medialspokesinterpolator->Update();

        m_Medialspokesinterpolator->GetOutput()->GetBounds(m_Bounds);

        m_Bounds[0] -= 0.1;
        m_Bounds[1] += 0.1;
        m_Bounds[2] -= 0.1;
        m_Bounds[3] += 0.1;
        m_Bounds[4] -= 0.1;
        m_Bounds[5] += 0.1;


        double spacing[3];
        spacing[0] = (m_Bounds[1] - m_Bounds[0])/((double)m_SizeX);
        spacing[1] = (m_Bounds[3] - m_Bounds[2])/((double)m_SizeY);
        spacing[2] = (m_Bounds[5] - m_Bounds[4])/((double)m_SizeZ);

        double origin[3];
        origin[0] = m_Bounds[0];
        origin[1] = m_Bounds[2];
        origin[2] = m_Bounds[4];

        if(this->GetDebug()){
            cout<<" x "<<m_Bounds[0]<<" "<<m_Bounds[1]<<endl;
            cout<<" Y "<<m_Bounds[2]<<" "<<m_Bounds[3]<<endl;
            cout<<" Z "<<m_Bounds[4]<<" "<<m_Bounds[5]<<endl;
            cout<<" xspc "<<spacing[0]<<endl;
            cout<<" Yspc "<<spacing[1]<<endl;
            cout<<" Zspc "<<spacing[2]<<endl;
        }

        DistanceImagePointerType imageobject = this->GetOutput();
        DistanceImageRegionType outregion;

        outregion.SetSize(0, m_SizeX);
        outregion.SetSize(1, m_SizeY);
        outregion.SetSize(2, m_SizeZ);

        imageobject->SetRegions(outregion);
        imageobject->SetSpacing(spacing);
        imageobject->SetOrigin(origin);
        imageobject->Allocate();


        DistanceImagePixelType pix = -1;
        imageobject->FillBuffer(pix);


    }else{
        itkExceptionMacro("Missing the surface information of the mrep ->SetPolyData");
    }

}

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
template< class TDistanceImage >
void SrepDistanceTransform< TDistanceImage >::BeforeThreadedGenerateData(){


}

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
template< class TDistanceImage >
void SrepDistanceTransform< TDistanceImage >::AfterThreadedGenerateData(){


}

/** Split the output's RequestedRegion into "num" pieces, returning
 * region "i" as "splitRegion". This method is called "num" times. The
 * regions must not overlap. The method returns the number of pieces that
 * the routine is capable of splitting the output RequestedRegion,
 * i.e. return value is less than or equal to "num". */
template< class TDistanceImage >
int SrepDistanceTransform< TDistanceImage >::SplitRequestedRegion(int i, int num, DistanceImageRegionType& splitRegion){
    return num;
}

}//namespace

#endif
