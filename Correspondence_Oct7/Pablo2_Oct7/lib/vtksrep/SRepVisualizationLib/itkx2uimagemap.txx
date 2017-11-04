#ifndef ITKX2UIMAGEMAP_TXX
#define ITKX2UIMAGEMAP_TXX

#include "itkImage.h"
#include "itkRGBToLuminanceImageFilter.h"
#include "vtkGenericCell.h"
#include "vtkCellArray.h"
#include "vtkVertex.h"

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
}

template< class TX2UImage, class TU2XImage >
X2UImageMap< TX2UImage, TU2XImage >::~X2UImageMap()
{
    m_SrepInterpolate.clear();
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
void X2UImageMap< TX2UImage, TU2XImage >::ThreadedGenerateData(const X2UImageRegionType& outputRegionForThread, ThreadIdType threadId ){

    if(this->GetDebug())
        cout<<"this is thread "<<threadId<<endl;

    X2URegionIteratorType outit(this->GetOutput(), this->GetOutput()->GetLargestPossibleRegion());
    LabelRegionIteratorType labelit(this->GetLabelImage(), this->GetLabelImage()->GetLargestPossibleRegion());

    if(this->GetDebug()){
        cout<<"this is thread = "<<threadId;
    }

    //if(threadId == 0)
    {
        if(true){

            double max = 255;
            double stepr = 0.025;
            double stepg = 0.025;
            double stepb = 5;
            double coeffb = 3;
            double maxb = max*coeffb;

            int nthreads = this->GetNumberOfThreads();
            //int nthreads = 1;
            double itersize = max / nthreads;

            double startpos = itersize * threadId;
            double endpos = startpos + itersize;

            for(double r = startpos; r <= endpos; r+=stepr){
                for(double g = 0; g <= max; g+=stepg){

                    vtkSRep::VNLType position = m_SrepInterpolate[threadId]->GetInterpolatedPoint(r/max, g/max, 0);
                    double interpolatedlabel = m_SrepInterpolate[threadId]->GetInterpolatedLabel();
                    vtkSRep::VNLType spoke = m_SrepInterpolate[threadId]->GetInterpolatedSpoke(r/max, g/max);

                    for(double b = 0; b <= maxb; b+=stepb){

                        X2UImagePixelType outpixel;
                        outpixel[0] = r;
                        outpixel[1] = g;
                        outpixel[2] = b;

                        vtkSRep::VNLType p = spoke*b/max + position;

                        typedef typename X2UImageType::PointType PointType;
                        PointType point;
                        point[0] = p[0];
                        point[1] = p[1];
                        point[2] = p[2];

                        X2UImageIndexType index;

                        if(this->GetOutput()->TransformPhysicalPointToIndex(point, index)){
                            outit.SetIndex(index);
                            labelit.SetIndex(index);

                            if(outpixel[2] <= max  || outit.Get()[2] == 0 || outpixel[2] < outit.Get()[2]){
                                outit.Set(outpixel);
                                labelit.Set(interpolatedlabel);
                            }
                        }
                    }
                }
            }
        }else{
            double max = 256;
            double stepr = 0.1;
            double stepg = 0.1;
            double stepb = 1;
            double coeffb = 1.0;
            double maxb = max*coeffb;

            int nthreads = this->GetNumberOfThreads();
            //int nthreads = 1;
            double itersize = max / nthreads;

            double startpos = itersize * threadId;
            double endpos = startpos + itersize;

            for(double r = startpos; r < endpos; r+=stepr){
                for(int side = 0; side < 2; side++){
                    for(double g = 0; g < max/2.0; g+=stepg){
                        double tempg = g;
                        if(side == 1){
                            tempg = 256 - g;
                        }

                        vtkSRep::VNLType position = m_SrepInterpolate[threadId]->GetInterpolatedPoint(r/max, tempg/max, 0);
                        vtkSRep::VNLType spoke = m_SrepInterpolate[threadId]->GetInterpolatedSpoke(r/max, tempg/max);

                        for(double b = 0; b < maxb; b+=stepb){

                            X2UImagePixelType outpixel;
                            outpixel[0] = r;
                            outpixel[1] = g*2.0;

                            double tempb = b/2.0 + 128;
                            if(side == 1){
                                tempb = 128 - b/2.0;
                            }

                            outpixel[2] = tempb;

                            vtkSRep::VNLType p = spoke*b/max + position;

                            typedef typename X2UImageType::PointType PointType;
                            PointType point;
                            point[0] = p[0];
                            point[1] = p[1];
                            point[2] = p[2];

                            X2UImageIndexType index;

                            if(this->GetOutput()->TransformPhysicalPointToIndex(point, index)){
                                outit.SetIndex(index);
                                outit.Set(outpixel);
                            }
                        }
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


        /*vtkSmartPointer<vtkSRepInterpolator> srepinterpolate = vtkSmartPointer<vtkSRepInterpolator>::New();
        srepinterpolate->SetInput(m_Input);
        srepinterpolate->Update();

        vtkSmartPointer<vtkUnsignedCharArray> cols = vtkSmartPointer<vtkUnsignedCharArray>::New();
        cols->SetName("cols");
        cols->SetNumberOfComponents(3);

        vtkSmartPointer<vtkPoints> sgp = vtkSmartPointer<vtkPoints>::New();
        vtkSmartPointer<vtkCellArray> scells = vtkSmartPointer<vtkCellArray>::New();

        double max = 256;
        //double maxb = 2;

        for(int r = 0; r <= max; r++){
            for(int g = 0; g <= max; g++){

                vtkSRep::VNLType position = srepinterpolate->GetInterpolatedPoint(r/max, g/max, 0);
                vtkSRep::VNLType spoke = srepinterpolate->GetInterpolatedSpoke(r/max, g/max);

                for(int b = 0; b <= 1; b++){

                    vtkSRep::VNLType p = spoke*b + position;


                    vtkIdType id = sgp->InsertNextPoint(p[0], p[1], p[2]);

                    cols->InsertNextTuple3(r/max*256, g/max*256, b*256);
                    vtkSmartPointer<vtkVertex> vert = vtkSmartPointer<vtkVertex>::New();
                    vert->GetPointIds()->SetId(0, id);
                    scells->InsertNextCell(vert);

                }
            }
        }

        m_Sgrid = vtkSmartPointer<vtkPolyData>::New();
        m_Sgrid->SetPoints(sgp);
        m_Sgrid->SetVerts(scells);
        m_Sgrid->BuildLinks();
        m_Sgrid->BuildCells();
        m_Sgrid->GetPointData()->AddArray(cols);*/

        vtkSRep::VNLType max(3);
        max.fill(VTK_DOUBLE_MIN);

        vtkSRep::VNLType min(3);
        min.fill(VTK_DOUBLE_MAX);

        for(unsigned i = 0; i < m_Input->GetNumberOfPoints(); i++){
            double p[3];
            m_Input->GetPoint(i, p);
            vtkSRep::VNLType p0(p, 3);
            vtkSRep::VectorVNLType spokes =  m_Input->GetSpokes(i);
            for(unsigned j = 0; j < spokes.size(); j++){
                vtkSRep::VNLType s = p0 + spokes[j];
                for(unsigned k = 0; k < s.size(); k++){
                    if(s[k] < min[k]){
                        min[k] = s[k];
                    }
                    if(s[k] > max[k]){
                        max[k] = s[k];
                    }
                }
            }
        }
        
        double dx = (max[0]-min[0]);
        double dy = (max[1]-min[1]);
        double dz = (max[2]-min[2]);

        m_Bounds[0] = min[0] - dx*0.1;
        m_Bounds[1] = max[0] + dx*0.1;
        m_Bounds[2] = min[1] - dy*0.1;
        m_Bounds[3] = max[1] + dy*0.1;
        m_Bounds[4] = min[2] - dz*0.1;
        m_Bounds[5] = max[2] + dz*0.1;


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


        m_LabelImage = LabelImageType::New();
        LabelImageRegionType regionlabel;

        regionlabel.SetSize(0,m_SizeX);
        regionlabel.SetSize(1,m_SizeY);
        regionlabel.SetSize(2,m_SizeZ);

        m_LabelImage->SetRegions(regionlabel);
        m_LabelImage->SetSpacing(spacing);
        m_LabelImage->SetOrigin(origin);
        m_LabelImage->Allocate();

        LabelImagePixelType label = 0;
        m_LabelImage->FillBuffer(label);


        /*m_U2XImage = U2XImageType::New();
        U2XImageRegionType regionU2X;

        regionU2X.SetSize(0,m_SizeX);
        regionU2X.SetSize(1,m_SizeY);
        regionU2X.SetSize(2,m_SizeZ);
        double u2xspacing[3] = {1,1,1};
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

    for(unsigned i = 0; i < this->GetNumberOfThreads();i++){
        m_SrepInterpolate.push_back(vtkSmartPointer<vtkSRepInterpolator>::New());
        m_SrepInterpolate[i]->SetInput(m_Input);
        m_SrepInterpolate[i]->Update();
    }
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
    vectorgaussian->SetVariance(0.5);
    vectorgaussian->Update();
    this->SetNthOutput(0, vectorgaussian->GetOutput());

    X2URegionIteratorType outit(this->GetOutput(), this->GetOutput()->GetLargestPossibleRegion());

    X2UImagePixelType texel;
    texel.Fill(0);
    outit.GoToBegin();
    while(!outit.IsAtEnd()){

        if(outit.Get()[2]>255 || outit.Get() == texel){
            outit.Set(texel);
        }

        ++outit;
    }*/

    /*m_Fillholes = X2UFillHolesType::New();
    m_Fillholes->SetInput(this->GetOutput());

    X2UImagePixelType pix;
    pix.Fill(0);
    m_Fillholes->SetBackgroundValue(pix);
    m_Fillholes->Update();
    this->SetNthOutput(0, m_Fillholes->GetOutput());*/




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
