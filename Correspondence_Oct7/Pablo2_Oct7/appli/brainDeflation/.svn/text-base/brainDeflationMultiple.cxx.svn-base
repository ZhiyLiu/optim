


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

    }

    ofstream outfile;
    outfile.open("lhall.curv.asc", ios::out);

    for(unsigned i = 0; i < vectcurvature.size(); i++){
        string line = "";

        for(unsigned j = 0; j < allcurvature.size(); j++){

            char buf[50];
            sprintf(buf, " %f ", allcurvature[j][i]);
            line.append(buf);

        }

        outfile << line <<endl;
    }

    outfile.close();


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

    vtkPoints* pointsorig = orig->GetPoints();
    vtkPoints* pointsinflated = inflated->GetPoints();


    vtkSmartPointer<vtkPolyDataMapper> inflatedmapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    inflatedmapper->SetInput(inflated);
    vtkSmartPointer<vtkActor>  vtkactorinflated = vtkActor::New();
    vtkactorinflated->SetMapper(inflatedmapper);


    vtkSmartPointer<vtkPolyDataMapper> origmapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    origmapper->SetInput(orig);
    vtkSmartPointer<vtkActor>  vtkactororig = vtkActor::New();
    vtkactororig->SetMapper(origmapper);    

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


     //vtkactorinflated->GetProperty()->SetOpacity(0.2);
    renderer->AddActor(vtkactororig);
    renderer->AddActor(vtkactorinflated);

    //orig->GetPointData()->SetScalars(curvaturearray);
    //origmapper->SetLookupTable(color);

    orig->GetPointData()->SetScalars(labelarray);
    origmapper->SetLookupTable(colorannot);

    //inflated->GetPointData()->SetScalars(curvaturearray);
    //inflatedmapper->SetLookupTable(color);

    inflated->GetPointData()->SetScalars(labelarray);
    inflatedmapper->SetLookupTable(colorannot);

     //vtkactorinflated->GetProperty()->SetOpacity(0.2);

    renderWindow->Render();

    renderWindowInteractor->Start();



    if(alllabelsmap.size() > 0)    {

        map<string, int> labelmap = alllabelsmap[0];

        map<string, int>::iterator it;

        for(it = labelmap.begin(); it != labelmap.end(); ++it){

            string key = (*it).first;



            for(unsigned i = 0; i < alllabelsmap.size(); i++){

                map<string, int> currentmap = alllabelsmap[0];

                int val = currentmap[key];


                vector< float > vectcurv = allcurvature[i];
                vnl_vector< int > vectlabel = alllabels[i];

                for(unsigned j = 0; j < vectlabel.size(); j++){
                    if(vectlabel[j] == val){

                    }
                }
            }
        }
    }



    return 0;
}
