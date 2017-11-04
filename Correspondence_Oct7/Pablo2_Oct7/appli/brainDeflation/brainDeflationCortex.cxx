


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

#include <vnl/vnl_vector.h>

#include <map>
#include <vector>

#define PI 3.14159265

using namespace std;

void help(char* execname){
    cout<<"Calculates the thickness surface for the cortex and displays it."<<endl;
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

    vector< vector< float > > allcurvature;
    vnl_vector<float> vectcurvature;
    vector< vnl_vector< int > > alllabels;
    vector< map< string, int > > alllabelsmap;

    vector< vnl_vector<float> > allthickness;

    for(unsigned i = 0; i < patientname.size(); i++){
        string curvfilename = dirname;
        curvfilename.append("/");
        curvfilename.append(patientname[i]);
        curvfilename.append("/surf/vtk/lh.curv.asc");

        string line;
        ifstream myfile (curvfilename.c_str());

        if (myfile.is_open())
        {

            vector<float> tempveccurvature;

            while ( myfile.good() ){
                getline (myfile,line);
                int pos0 = 0;
                pos0 = line.find_last_of(" ");

                if(pos0 != -1){
                    string num = line.substr(pos0, line.size() - pos0);


                    tempveccurvature.push_back(atof(num.c_str()));
                }



            }
            myfile.close();

            if(vectcurvature.size() == 0){
                vectcurvature.set_size(tempveccurvature.size());
                vectcurvature.fill(0);
            }

            for(unsigned i = 0; i < vectcurvature.size(); i++){
                vectcurvature[i] += tempveccurvature[i];
            }

            allcurvature.push_back(tempveccurvature);

            tempveccurvature.clear();
        }else{
             cout << "Unable to label intensities file from freesurfer. <subject>/mri/aseg.auto_noCCseg.label_intensities.txt";
             return 0;
        }



        string annotfilename = dirname;
        annotfilename.append("/");
        annotfilename.append(patientname[i]);
        annotfilename.append("/surf/vtk/lh.aparc.stats.asc");


        ifstream annotfile (annotfilename.c_str());

        if (annotfile.is_open())
        {

            vector<int> tempvecannot;

            while ( annotfile.good() ){
                getline (annotfile,line);
                int pos0 = 0;
                pos0 = line.find_last_of(" ");

                if(pos0 != -1){
                    string num = line.substr(pos0, line.size() - pos0);


                    tempvecannot.push_back(atoi(num.c_str()));
                }



            }
            annotfile.close();

            vnl_vector<int> vnltemp;
            vnltemp.set_size(tempvecannot.size());
            for(unsigned i=0; i < tempvecannot.size(); i++){
                vnltemp[i] = tempvecannot[i];
            }

            alllabels.push_back(vnltemp);

            tempvecannot.clear();
        }else{
             cout << "Unable to get annotation file";
             return 0;
        }


        string labelannotfilename = dirname;
        labelannotfilename.append("/");
        labelannotfilename.append(patientname[i]);
        labelannotfilename.append("/surf/vtk/lh.aparc.stats.label.asc");


        ifstream labelfile (labelannotfilename.c_str());

        map< string, int > maplabels;

        if (labelfile.is_open()){

            while ( labelfile.good() ){

                getline (labelfile,line);
                int pos0 = 0;
                int pos1 = 0;
                pos0 = line.find(" ");
                pos1 = line.find(",", pos0+1);

                if(pos0 != -1 && pos1 != -1){

                    string label = line.substr(pos0, pos1-pos0);

                    pos0 = line.find_last_of(" ");
                    int numlabel = -1;

                    if(pos0 != -1){
                        numlabel = atoi(line.substr(pos0, line.size() - pos0).c_str());
                    }


                    maplabels[label] = numlabel;
                }

            }
            labelfile.close();

        }else{
             cout << "Unable to get annotation file";
             return 0;
        }


        alllabelsmap.push_back(maplabels);



        string thickfilename = dirname;
        thickfilename.append("/");
        thickfilename.append(patientname[i]);
        thickfilename.append("/surf/vtk/lh.thickness.asc");

        ifstream thickfile (thickfilename.c_str());

        if (thickfile.is_open())
        {

            vnl_vector<float> vectthickness;
            vector<float> tempvecthick;

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


    vectcurvature = vectcurvature / (patientname.size());


    vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
    renderer->SetBackground(1,1,1);
    vtkSmartPointer<vtkRenderWindow> renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
    renderWindow->AddRenderer(renderer);
    vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
    renderWindowInteractor->SetRenderWindow(renderWindow);


    vtkSmartPointer<vtkPolyDataReader> reader = vtkSmartPointer<vtkPolyDataReader>::New();

    string origfilename = dirname;
    origfilename.append("/");
    origfilename.append(patientname[0]);
    origfilename.append("/surf/vtk/lh.orig.vtk");

    reader->SetFileName(origfilename.c_str());
    reader->Update();
    vtkPolyData* orig = reader->GetOutput();


    vtkSmartPointer<vtkPolyDataReader> reader1 =  vtkSmartPointer<vtkPolyDataReader>::New();

    string inflatedfilename = dirname;
    inflatedfilename.append("/");
    inflatedfilename.append(patientname[0]);
    inflatedfilename.append("/surf/vtk/lh.inflated.vtk");

    reader1->SetFileName(inflatedfilename.c_str());
    reader1->Update();
    vtkPolyData* inflated = reader1->GetOutput();


    string pialfilename = dirname;
    pialfilename.append("/");
    pialfilename.append(patientname[0]);
    pialfilename.append("/surf/vtk/lh.pial.vtk");

    vtkSmartPointer<vtkPolyDataReader> reader2 =  vtkSmartPointer<vtkPolyDataReader>::New();
    reader2->SetFileName(pialfilename.c_str());
    reader2->Update();
    vtkPolyData* pial = reader2->GetOutput();


    vtkSmartPointer<vtkPolyDataMapper> inflatedmapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    inflatedmapper->SetInput(inflated);
    vtkSmartPointer<vtkActor>  vtkactorinflated = vtkActor::New();
    vtkactorinflated->SetMapper(inflatedmapper);


    vtkSmartPointer<vtkPolyDataMapper> origmapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    origmapper->SetInput(orig);
    vtkSmartPointer<vtkActor>  vtkactororig = vtkActor::New();
    vtkactororig->SetMapper(origmapper);


    vtkSmartPointer<vtkPolyDataMapper> pialmapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    pialmapper->SetInput(pial);
    vtkSmartPointer<vtkActor>  vtkactorpial = vtkActor::New();
    vtkactorpial->SetMapper(pialmapper);

    vtkSmartPointer<vtkColorTransferFunction> colorannot = vtkSmartPointer<vtkColorTransferFunction>::New();
    colorannot->SetColorSpaceToHSV();

    vtkDataArray* labelarray = (vtkDataArray*)vtkDataArray::CreateArray(VTK_DOUBLE);
    //double* curvrange = curv->GetOutput()->GetPointData()->GetScalars()->GetRange();


    for(unsigned i = 0; i < 1; i++){
        vnl_vector<int> vectlabel = alllabels[i];
        map< string, int > maplabels = alllabelsmap[i];


        double step = 1.0/((double)maplabels.size());
        double hue = 0;

        map< string, int >::iterator it;


        for ( it=maplabels.begin() ; it != maplabels.end(); it++ ){

            int label = (*it).second;
            colorannot->AddHSVPoint(label, hue, 1, 1);
            hue += step;
        }


        for(unsigned i = 0; i < vectlabel.size();i++){
            double val[1];
            val[0] = vectlabel[i];
            labelarray->InsertTuple(i, val);
        }
    }

    colorannot->Build();


    vtkSmartPointer<vtkColorTransferFunction> color = vtkSmartPointer<vtkColorTransferFunction>::New();
    color->SetColorSpaceToRGB();

    //double* curvrange = curv->GetOutput()->GetPointData()->GetScalars()->GetRange();

    color->AddRGBPoint(-0.1,1,0,0);
    color->AddRGBPoint(0.1,0,1,0);

    color->Build();

    vtkDataArray* curvaturearray = (vtkDataArray*)vtkDataArray::CreateArray(VTK_DOUBLE);

    for(unsigned i = 0; i < vectcurvature.size();i++){
        double val[1];
        val[0] = vectcurvature[i];
        curvaturearray->InsertTuple(i, val);
    }





    //thickness stuff
    vtkSmartPointer<vtkColorTransferFunction> colorthick = vtkSmartPointer<vtkColorTransferFunction>::New();
    colorthick->SetColorSpaceToHSV();

    vtkDataArray* thickarray = (vtkDataArray*)vtkDataArray::CreateArray(VTK_DOUBLE);
    //double* curvrange = curv->GetOutput()->GetPointData()->GetScalars()->GetRange();

    for(unsigned i = 0; i < 1; i++){
        vnl_vector<float> vectthick = allthickness[i];

        for(unsigned i = 0; i < vectthick.size();i++){
            double val[1];
            val[0] = vectthick[i];
            thickarray->InsertTuple(i, val);
        }

        double step = 1.0/5.0;
        double hue = 0;


        for ( hue; hue < 1; hue+=step ){
            colorthick->AddHSVPoint(hue*5, hue, 1, 1);
        }
    }

    colorthick->Build();

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


    vnl_vector<float> vectthick = allthickness[0];

    for(unsigned i = 0; i < inflated->GetNumberOfPoints(); i++){

        double normal[3];
        normalDataArray->GetTuple(i, normal);

        double thickness = vectthick[i];

        double point[3];
        inflated->GetPoint(i, point);

        inflatedpialpoints->InsertNextPoint(point[0] + thickness*normal[0], point[1]  + thickness*normal[1], point[2]  + thickness*normal[2]);

    }

    inflatedpial->SetPoints(inflatedpialpoints);
    inflatedpial->SetPolys(inflated->GetPolys());

    vtkSmartPointer<vtkImageData> inflatedimage = vtkSmartPointer<vtkImageData>::New();
    vtkSmartPointer<vtkImageData> inflatedpialimage = vtkSmartPointer<vtkImageData>::New();

    PolyDataImageStencilExport(inflated, inflatedimage, inflatedpial->GetBounds());
    PolyDataImageStencilExport(inflatedpial, inflatedpialimage, inflatedpial->GetBounds());


    vtkSmartPointer<vtkPolyDataMapper> inflatedpialmapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    inflatedpialmapper->SetInput(inflatedpial);
    vtkSmartPointer<vtkActor>  vtkactorinflatedpial = vtkActor::New();
    vtkactorinflatedpial->SetMapper(inflatedpialmapper);



    cout<<normalDataArray->GetNumberOfTuples()<<" "<<thickarray->GetNumberOfTuples()<<" "<<inflated->GetNumberOfPoints()<< endl;
    inflated->GetPointData()->SetNormals(normalDataArray);
    normalDataArray =  inflated->GetPointData()->GetArray("Normals");
    cout<<normalDataArray->GetNumberOfTuples()<<" "<<thickarray->GetNumberOfTuples()<<" "<<inflated->GetNumberOfPoints()<< endl;



    //thickness stuff end
    //vtkactorinflated->GetProperty()->SetOpacity(0.2);
    //vtkactorpial->GetProperty()->SetOpacity(0.2);
    //renderer->AddActor(vtkactororig);
    vtkactorinflated->GetProperty()->SetRepresentationToWireframe();
    renderer->AddActor(vtkactorinflated);   

    vtkactorinflatedpial->GetProperty()->SetOpacity(0.2);
    //renderer->AddActor(vtkactorinflatedpial);

    //renderer->AddActor(vtkactorpial);


    //orig->GetPointData()->SetScalars(curvaturearray);
    //origmapper->SetLookupTable(color);

    orig->GetPointData()->SetScalars(labelarray);
    origmapper->SetLookupTable(colorannot);

    //inflated->GetPointData()->SetScalars(curvaturearray);
    //inflatedmapper->SetLookupTable(color);

    //inflated->GetPointData()->SetScalars(labelarray);
    //inflatedmapper->SetLookupTable(colorannot);

    inflated->GetPointData()->SetScalars(thickarray);
    inflatedmapper->SetLookupTable(colorthick);

    inflatedpial->GetPointData()->SetScalars(labelarray);
    inflatedpialmapper->SetLookupTable(colorannot);

    pial->GetPointData()->SetScalars(labelarray);
    pialmapper->SetLookupTable(colorannot);

     //vtkactorinflated->GetProperty()->SetOpacity(0.2);

    renderWindow->Render();

    renderWindowInteractor->Start();


    vtkSmartPointer<vtkMetaImageWriter> writer = vtkSmartPointer<vtkMetaImageWriter>::New();
    writer->SetInput(inflatedimage);
    writer->SetFileName("inflatedImage.mhd");
    writer->SetCompression(false);
    writer->Write();

    writer->SetInput(inflatedpialimage);
    writer->SetFileName("inflatedPialImage.mhd");
    writer->Write();


    vtkSmartPointer<vtkPolyDataWriter> polywriter = vtkSmartPointer<vtkPolyDataWriter>::New();
    polywriter->SetInput(inflated);
    polywriter->SetFileName("inflated.vtk");
    polywriter->Write();

    polywriter->SetInput(inflatedpial);
    polywriter->SetFileName("inflatedpial.vtk");
    polywriter->Write();




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

    writer->SetInput(inflatedpialimage);
    writer->SetFileName("CortexImage.mhd");
    writer->Write();

    return 0;
}

