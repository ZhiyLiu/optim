


#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkPolyDataReader.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkArrowSource.h>
#include <vtkPoints.h>

#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkProperty.h>

#include <vtkLookupTable.h>
#include <vtkColorTransferFunction.h>
#include <vtkCurvatures.h>
#include <vtkPointData.h>

#include <vtkSmoothPolyDataFilter.h>
#include <vtkPolyDataNormals.h>

#include <vtkDoubleArray.h>

#include <vtkImageData.h>
#include <vtkPolyDataToImageStencil.h>
#include <vtkImageStencil.h>

#include <vtkMetaImageWriter.h>
#include <vtkPolyDataWriter.h>

#include <vtkImageDilateErode3D.h>

#include <vnl/vnl_vector.h>

#include <map>
#include <vector>

#define PI 3.14159265

using namespace std;

void help(char* execname){
    cout<<"Calculates the mean curvature over a set of patients point wise and visualize the result on the first patient inflated surface."<<endl;
    cout<<"Usage: "<<execname<<" -d <patients dir> -p <patient name> [options]"<<endl;
    cout<<"options:"<<endl;
    cout<<"--h --help show help menu"<<endl;
    cout<<"-p <patient name> add another patient";

}


void PolyDataImageStencilExport(vtkPolyData* polydata, vtkImageData *imagein, double* bounds){



    //vtkPolyData* polydata = (vtkPolyData*) poly->GetInput();

    vtkSmartPointer<vtkPolyDataToImageStencil> polytostencil = vtkSmartPointer<vtkPolyDataToImageStencil>::New();

    polytostencil->SetInput(polydata);

    polytostencil->Update();



    //double *bounds = polydata->GetBounds();

   // vtkSmartPointer<vtkImageData> imagein = vtkSmartPointer<vtkImageData>::New();



    imagein->SetExtent(bounds[0] - 1, bounds[1] + 1, bounds[2] - 1, bounds[3] + 1, bounds[4] - 1, bounds[5] + 1);

    imagein->SetScalarTypeToUnsignedShort();

    imagein->AllocateScalars();





    int* extent = imagein->GetExtent();



    for (int x = extent[0]; x <= extent[1]; x++){

        for (int y = extent[2]; y <= extent[3]; y++){

            for (int z  =extent[4]; z <= extent[5]; z++){

                unsigned short* pixel = static_cast<unsigned short*>(imagein->GetScalarPointer(x,y,z));

                *pixel = 0;

            }

        }

    }



    vtkSmartPointer<vtkImageStencil> stencil = vtkSmartPointer<vtkImageStencil>::New();

    stencil->SetInput(imagein);

    stencil->SetStencil(polytostencil->GetOutput());

    stencil->ReverseStencilOn();

    stencil->SetBackgroundValue(128);

    stencil->Update();

    imagein->DeepCopy(stencil->GetOutput());



}


int main(int argc, char *argv[])
{


    if (argc == 1){
        help(argv[0]);
        return 0;
    }

    string dirname = "";
    vector< std::string > patientname;

    for(int i = 1; i < argc; i++){

        if(string(argv[i]) == "--h" || string(argv[i]) == "--help"){
            help(argv[0]);
            return 0;
        }else if (string(argv[i]) == "-d"){
            if(i + 1 >= argc){
                std::cout << "ERROR: Sample filename missing."
                          <<std::endl
                          <<"-sp <filename>"
                          << std::endl ;
                return EXIT_FAILURE;
            }

            dirname = argv[i + 1];


        }else if(string(argv[i]) == "-p"){

            //imagesize[0] = atoi(argv[i+1]);
            patientname.push_back(argv[i + 1]);

        }

    }


    vector< vnl_vector<float> > allthickness;

    for(unsigned i = 0; i < patientname.size(); i++){

        string thickfilename = dirname;
        thickfilename.append("/");
        thickfilename.append(patientname[i]);
        thickfilename.append("/surf/vtk/lh.thickness.asc");

        ifstream thickfile (thickfilename.c_str());

        if (thickfile.is_open())
        {

            vnl_vector<float> vectthickness;
            vector<float> tempvecthick;

            string line = "";

            while ( thickfile.good() ){
                getline (thickfile,line);
                int pos0 = 0;
                pos0 = line.find_last_of(" ");

                if(pos0 != -1){
                    string num = line.substr(pos0, line.size() - pos0);


                    tempvecthick.push_back(atof(num.c_str()));
                }



            }
            thickfile.close();

            if(vectthickness.size() == 0){
                vectthickness.set_size(tempvecthick.size());
                vectthickness.fill(0);
            }

            for(unsigned i = 0; i < tempvecthick.size(); i++){
                vectthickness[i] = tempvecthick[i];
            }

            allthickness.push_back(vectthickness);

            vectthickness.clear();
            tempvecthick.clear();
        }else{
             cout << "Unable to label intensities file from freesurfer. <subject>/mri/aseg.auto_noCCseg.label_intensities.txt";
             return 0;
        }

    }

    for(unsigned n = 0; n < patientname.size(); n++){

        vtkSmartPointer<vtkPolyDataReader> reader1 =  vtkSmartPointer<vtkPolyDataReader>::New();

        string inflatedfilename = dirname;
        inflatedfilename.append("/");
        inflatedfilename.append(patientname[n]);
        inflatedfilename.append("/surf/vtk/lh.inflated.vtk");

        reader1->SetFileName(inflatedfilename.c_str());
        reader1->Update();
        vtkPolyData* inflated = reader1->GetOutput();



        vtkDataArray* normalDataArray =  inflated->GetPointData()->GetArray("Normals");

        vtkSmartPointer<vtkPolyDataNormals> normalGenerator = vtkSmartPointer<vtkPolyDataNormals>::New();

        cout<<inflated->GetNumberOfPoints()<< endl;

        vtkPolyData* inflatednormals = 0;

        if(!normalDataArray ){
            normalGenerator->SetInput(inflated);
            normalGenerator->ComputePointNormalsOn();
            normalGenerator->ComputeCellNormalsOff();
            normalGenerator->SetConsistency(1);
            normalGenerator->Update();
            inflatednormals = normalGenerator->GetOutput();
            normalDataArray =  inflatednormals->GetPointData()->GetArray("Normals");
        }

        vtkSmartPointer<vtkPolyData> inflatedpial = vtkSmartPointer<vtkPolyData>::New();
        vtkSmartPointer<vtkPoints> inflatedpialpoints = vtkSmartPointer<vtkPoints>::New();


        vnl_vector<float> vectthick = allthickness[n];

        for(unsigned i = 0; i < inflated->GetNumberOfPoints(); i++){

            double normal[3];
            normalDataArray->GetTuple(i, normal);

            double thickness = vectthick[i];

            double point[3];
            inflated->GetPoint(i, point);

            inflatedpialpoints->InsertNextPoint(point[0] + thickness*normal[0], point[1]  + thickness*normal[1], point[2]  + thickness*normal[2]);

            int celltype = inflated->GetCellType(i);
            //cout<<celltype<<endl;

        }

        inflatedpial->SetPoints(inflatedpialpoints);
        inflatedpial->SetPolys(inflated->GetPolys());

        vtkSmartPointer<vtkImageData> inflatedimage = vtkSmartPointer<vtkImageData>::New();
        vtkSmartPointer<vtkImageData> inflatedpialimage = vtkSmartPointer<vtkImageData>::New();

        PolyDataImageStencilExport(inflated, inflatedimage, inflatedpial->GetBounds());
        PolyDataImageStencilExport(inflatedpial, inflatedpialimage, inflatedpial->GetBounds());


        int* extent = inflatedpialimage->GetExtent();

        for (int x = extent[0]; x <= extent[1]; x++){

            for (int y = extent[2]; y <= extent[3]; y++){

                for (int z  =extent[4]; z <= extent[5]; z++){

                    unsigned short* pixelpialinflated = static_cast<unsigned short*>(inflatedpialimage->GetScalarPointer(x,y,z));

                    unsigned short* pixelinflated = static_cast<unsigned short*>(inflatedimage->GetScalarPointer(x,y,z));

                    if(*pixelinflated !=0){
                        *pixelpialinflated = 0;
                    }

                }

            }

        }

        vtkSmartPointer<vtkImageDilateErode3D> dilate = vtkSmartPointer<vtkImageDilateErode3D>::New();
        dilate->SetInput(inflatedpialimage);
        dilate->SetDilateValue(128);
        dilate->SetErodeValue(0);
        dilate->SetKernelSize(2,2,2);
        dilate->Update();


        string cortexfilename = dirname;
        cortexfilename.append("/");
        cortexfilename.append(patientname[n]);
        cortexfilename.append("/surf/vtk/lh.inflatedSheet.mhd");

        cout<<"writing to: "<<cortexfilename<<endl;

        vtkSmartPointer<vtkMetaImageWriter> writer = vtkSmartPointer<vtkMetaImageWriter>::New();

        writer->SetCompression(false);
        writer->SetInput(dilate->GetOutput());
        writer->SetFileName(cortexfilename.c_str());
        writer->Write();
    }

    return 0;
}


