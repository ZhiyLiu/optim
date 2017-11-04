#ifndef ITKX2UIMAGEMAP_TXX
#define ITKX2UIMAGEMAP_TXX

#include "itkImage.h"
#include "itkRGBToLuminanceImageFilter.h"
#include "vtkGenericCell.h"

namespace itk{

template< class TX2UImage, class TU2XImage >
X2UImageMap< TX2UImage, TU2XImage >::X2UImageMap()
{
    //m_SampleANN = 0;
    //m_ObjectSurface = 0;
    m_Input = 0;
    m_SizeX = 256.0;
    m_SizeY = 256.0;
    m_SizeZ = 256.0;
    //m_Step = 1;
    //m_Dimension = 0;
    m_U2XImage = 0;
    m_BinaryImage = false;
    m_DistanceTransformImage = false;
}

template< class TX2UImage, class TU2XImage >
X2UImageMap< TX2UImage, TU2XImage >::~X2UImageMap()
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
template< class TX2UImage, class TU2XImage >
void X2UImageMap< TX2UImage, TU2XImage >::ThreadedGenerateData(const X2UImageRegionType& outputRegionForThread, int threadId ){

    if(this->GetDebug())
        cout<<"this is thread "<<threadId<<endl;


    typename X2UImageType::SpacingType spacing = this->GetOutput()->GetSpacing();
    double xspc = spacing[0];
    double yspc = spacing[1];
    double zspc = spacing[2];
    //double minspc = min(min(xspc, yspc), zspc)/8.0;
    double minspc = min(min(xspc, yspc), zspc)/2.0;
    //const double* bounds = this->m_Bounds;

    X2URegionIteratorType outit(this->GetOutput(), this->GetOutput()->GetLargestPossibleRegion());

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
    /*{
        for(unsigned i = startpos; i < endpos && i < srepmedial->GetNumberOfPoints(); i++){
        //for(unsigned i = 0; i < srepmedial->GetNumberOfPoints(); i++){

            vtkSRep::VNLType uv = srepmedial->GetUVCoord(i);

            //cout<<uv<<endl;

            double r = 255.0*(uv[0])/(m_Input->GetNumColumns());
            double g = 255.0*(uv[1])/(m_Input->GetNumRows());

            vtkSRep::VectorVNLType spokes = srepmedial->GetSpokes(i);
            double point[3];
            srepmedial->GetPoint(i, point);
            vtkSRep::VNLType position(point, 3);

            for(unsigned j = 0; j < spokes.size(); j++){

                vtkSRep::VNLType spoke = spokes[j];

                for(double t = 0; t <= 1; t+= 0.01){

                    vtkSRep::VNLType currentpoint = spoke * t + position;
                    X2UImagePixelType outpixel;
                    outpixel[0] = (unsigned char)(r);
                    if(j == 0){

                        outpixel[1] = (unsigned char)(g/2.0);
                    }else if(j == 1){

                        outpixel[1] = (unsigned char)(255.0 - g/2.0);
                    }

                    double b = t*250.0;
                    if(b > 255){
                        b = 255;
                    }                    
                    outpixel[2] = (unsigned char) b;


                    typedef typename X2UImageType::PointType PointType;
                    PointType point;
                    point[0] = currentpoint[0];
                    point[1] = currentpoint[1];
                    point[2] = currentpoint[2];

                    X2UImageIndexType index;

                    if(this->GetOutput()->TransformPhysicalPointToIndex(point, index)){
                        outit.SetIndex(index);
                        outit.Set(outpixel);
                    }
                }
            }
        }
    }*/

    //if(threadId == 0)
    {
        for(int i = startpos; i < endpos && i < (int)m_Input->GetNumberOfCells(); i++){
        //for(int i = 0; i < (int)m_Input->GetNumberOfCells(); i++){
            vtkSmartPointer<vtkGenericCell> cell = vtkSmartPointer<vtkGenericCell>::New();
            m_Input->GetCell(i, cell);
            vtkIdType id0 = cell->GetPointIds()->GetId(0);

            int x = id0/m_Input->GetNumColumns();
            int y = id0%m_Input->GetNumColumns();


            for(double u = 0; u <= 1; u += minspc){

                double r = 255.0*(y + u)/(m_Input->GetNumColumns()-1);

                for(double v = 0; v <= 1; v += minspc){

                    double g = 255.0*(x + v)/(m_Input->GetNumRows()-1);


                    vtkSRep::VNLType position = m_Medialsheetinterpolator->GetInterpolatedPoint(i, u, v);

                    for(unsigned side = 0; side < 2; side++){
                        vtkSRep::VNLType spoke = m_Medialspokesinterpolator->GetInterpolatedSpoke(i, side, u, v);
                        for(double t = 0; t <= 1.05; t+= minspc){

                            vtkSRep::VNLType currentpoint = spoke * t + position;
                            X2UImagePixelType outpixel;
                            outpixel[0] = (unsigned char)(r);
                            if(side == 0){
                                outpixel[1] = (unsigned char)(g/2.0);
                            }else if(side == 1){

                                outpixel[1] = (unsigned char)(255.0 - g/2.0);
                            }

                            double b = t*250.0;
                            if(b > 255){
                                b = 255;
                            }
                            outpixel[2] = (unsigned char) b;


                            typedef typename X2UImageType::PointType PointType;
                            PointType point;
                            point[0] = currentpoint[0];
                            point[1] = currentpoint[1];
                            point[2] = currentpoint[2];

                            X2UImageIndexType index;

                            if(this->GetOutput()->TransformPhysicalPointToIndex(point, index)){
                                outit.SetIndex(index);
                                outit.Set(outpixel);
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
        vtkSRep::VectorIdsType crestids = m_Input->GetCrestMedialAtomsIds();

        int itersizecrest = ceil((double)crestids.size() / nthreads);

        int startposcrest = itersizecrest * threadId;
        int endposcrest = startposcrest + itersizecrest;

        vtkSmartPointer<vtkSRepInterpolateCrestSpokesQuartic> interpolatecrestspokes = vtkSmartPointer<vtkSRepInterpolateCrestSpokesQuartic>::New();
        interpolatecrestspokes->SetSpokeType(vtkSRep::CREST_SPOKE);
        interpolatecrestspokes->SetInput(m_Input);
        interpolatecrestspokes->Update();


        for(int i = startposcrest; i < endposcrest && i < crestids.size(); i++){
        //for(int i = 0; i < crestids.size(); i++){

            vtkIdType id0 = crestids[i];
            vtkIdType id1;
            if(i < crestids.size()-1){
                id1 = crestids[i+1];
            }else{
                id1 = crestids[0];
            }

            int x = id0/m_Input->GetNumColumns();
            int y = id0%m_Input->GetNumColumns();
            int x1 = id1/m_Input->GetNumColumns();
            int y1 = id1%m_Input->GetNumColumns();

            for(double u = 0; u <= 1; u += minspc){

                double r = 255.0/(m_Input->GetNumColumns()-1);
                if(y1 - y == 1){
                    r*=(y + (y1-y)*u);
                }else{
                    r*=(y + (y1-y)*u);
                }

                double g = 255.0/(m_Input->GetNumRows()-1);
                if(x1-x == 1){
                    g*=(x + (x1-x)*u);
                }else{
                    g*=(x + (x1-x)*u);
                }


                vtkSRep::VNLType position = interpolatecrestspokes->GetInterpolatedPoint(i, u);

                for(double v = 0; v <= 1; v += minspc){

                    vtkSRep::VNLType spoke;
                    if(v == 0){
                        spoke = interpolatecrestspokes->GetInterpolatedSpoke(i, u, v);
                    }else{
                        spoke = interpolatecrestspokes->GetInterpolatedSpoke(v);
                    }

                    X2UImagePixelType outpixel;
                    outpixel[0] = (unsigned char)(r);
                    if(v > 0.5){
                        outpixel[1] = (unsigned char)(g/2.0);
                    }else{
                        outpixel[1] = (unsigned char)(255.0 - g/2.0);
                    }

                    for(double t = 0; t <= 1; t+= minspc){

                        vtkSRep::VNLType currentpoint = spoke * t + position;

                        double b = t*255.0;
                        if(b > 255){
                            b = 255;
                        }
                        outpixel[2] = (unsigned char) b;


                        typedef typename X2UImageType::PointType PointType;
                        PointType point;
                        point[0] = currentpoint[0];
                        point[1] = currentpoint[1];
                        point[2] = currentpoint[2];

                        X2UImageIndexType index;

                        if(this->GetOutput()->TransformPhysicalPointToIndex(point, index)){
                            outit.SetIndex(index);
                            outit.Set(outpixel);
                        }

                        //cout<<"scalarpoint "<<(int)scalarpoint[0] <<" "<<(int)scalarpoint[1] <<" "<<(int)scalarpoint[2] <<" "<<endl;
                    }

                }
            }
        }
    }

}


/** The GenerateData method normally allocates the buffers for all of the
 * outputs of a filter. Some filters may want to override this default
 * behavior. For example, a filter may have multiple outputs with
 * varying resolution. Or a filter may want to process data in place by
 * grafting its input to its output. */
template< class TX2UImage, class TU2XImage >
void X2UImageMap< TX2UImage, TU2XImage >::AllocateOutputs(){

    if(m_Input){        

        m_Medialsheetinterpolator = vtkSmartPointer< vtkSRepInterpolateMedialSheet >::New();
        m_Medialsheetinterpolator->SetInput(m_Input);
        m_Medialsheetinterpolator->Update();

        m_Medialspokesinterpolator = vtkSmartPointer< vtkSRepInterpolateMedialSpokesHermite >::New();
        m_Medialspokesinterpolator->SetInput(m_Input);        
        m_Medialspokesinterpolator->Update();

        m_Medialspokesinterpolator->GetOutput()->GetBounds(m_Bounds);

        m_Bounds[0] -= 1;
        m_Bounds[1] += 1;
        m_Bounds[2] -= 1;
        m_Bounds[3] += 1;
        m_Bounds[4] -= 1;
        m_Bounds[5] += 1;


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

        X2UImagePointerType imageobject = this->GetOutput();
        X2UImageRegionType outregion;

        outregion.SetSize(0, m_SizeX);
        outregion.SetSize(1, m_SizeY);
        outregion.SetSize(2, m_SizeZ);

        imageobject->SetRegions(outregion);
        imageobject->SetSpacing(spacing);
        imageobject->SetOrigin(origin);
        imageobject->Allocate();


        X2UImagePixelType texel;
        texel.Fill(0);
        imageobject->FillBuffer(texel);



        /*m_U2XImage = U2XImageType::New();
        U2XImageRegionType regionU2X;
        regionU2X.SetSize(0,128);
        regionU2X.SetSize(1,128);
        regionU2X.SetSize(2,128);
        double u2xspacing[3] = {2,2,2};
        double u2xorigin[3] = {0,0,0};
        m_U2XImage->SetRegions(regionU2X);
        m_U2XImage->SetSpacing (u2xspacing);
        m_U2XImage->SetOrigin(u2xorigin);
        m_U2XImage->Allocate();

        U2XImagePixelType texelu2x;
        texelu2x.Fill(0);
        m_U2XImage->FillBuffer(texelu2x);*/



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
template< class TX2UImage, class TU2XImage >
void X2UImageMap< TX2UImage, TU2XImage >::BeforeThreadedGenerateData(){


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
template< class TX2UImage, class TU2XImage >
void X2UImageMap< TX2UImage, TU2XImage >::AfterThreadedGenerateData(){

    /*typedef itk::VectorDiscreteGaussianImageFilter< X2UImageType, X2UImageType > VectorDiscreteGaussianType;
    typename VectorDiscreteGaussianType::Pointer vectorgaussian = VectorDiscreteGaussianType::New();

    vectorgaussian->SetInput(this->GetOutput());
    vectorgaussian->SetVariance(1);
    vectorgaussian->Update();
    this->SetNthOutput(0, vectorgaussian->GetOutput());*/


}

/** Split the output's RequestedRegion into "num" pieces, returning
 * region "i" as "splitRegion". This method is called "num" times. The
 * regions must not overlap. The method returns the number of pieces that
 * the routine is capable of splitting the output RequestedRegion,
 * i.e. return value is less than or equal to "num". */
template< class TX2UImage, class TU2XImage >
int X2UImageMap< TX2UImage, TU2XImage >::SplitRequestedRegion(int i, int num, X2UImageRegionType& splitRegion){
    return num;
}

}//namespace

#endif
