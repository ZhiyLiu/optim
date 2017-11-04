


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
#include <vtkCell.h>
#include <vtkCellArray.h>
#include <vtkQuad.h>


#include <vtkPointLocator.h>
#include <vtkSphereSource.h>


#include <vtksrep.h>
#include <vtksrepinterpolatemedialsheet.h>
#include <vtksrepvisuprimitives.h>
#include <vtksrepinterpolatemedialcrestcurve.h>
#include <vtksrepinterpolatemedialspokes.h>
#include <vtksrepinterpolatecrestspokes.h>
#include <vtksrepinterpolatecrestspokesquartic.h>
#include <vtkTransform.h>
#include <vtkTransformPolyDataFilter.h>

#include <vnl/vnl_vector.h>

#include <map>

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

        map<string,int>::iterator it;
        for ( it=maplabels.begin() ; it != maplabels.end(); it++ )
            cout << (*it).first << " => " << (*it).second << endl;


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

    vtkSmartPointer<vtkPolyDataMapper> origmapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    origmapper->SetInput(orig);
    vtkSmartPointer<vtkActor>  vtkactororig = vtkActor::New();
    vtkactororig->SetMapper(origmapper);
    renderer->AddActor(vtkactororig);

    string pialfilename = dirname;
    pialfilename.append("/");
    pialfilename.append(patientname[0]);
    pialfilename.append("/surf/vtk/lh.pial.vtk");

    vtkSmartPointer<vtkPolyDataReader> reader2 =  vtkSmartPointer<vtkPolyDataReader>::New();
    reader2->SetFileName(pialfilename.c_str());
    reader2->Update();
    vtkPolyData* pial = reader2->GetOutput();


    vtkSmartPointer<vtkPolyDataReader> reader3 =  vtkSmartPointer<vtkPolyDataReader>::New();

    string spherefilename = dirname;
    spherefilename.append("/");
    spherefilename.append(patientname[0]);
    spherefilename.append("/surf/vtk/lh.sphere.vtk");

    reader3->SetFileName(spherefilename.c_str());
    reader3->Update();
    vtkPolyData* sphere = reader3->GetOutput();

    /*vtkSmartPointer<vtkPolyDataMapper> spheremapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    spheremapper->SetInput(sphere);
    vtkSmartPointer<vtkActor>  vtkactorsphere = vtkActor::New();
    vtkactorsphere->SetMapper(spheremapper);
    //vtkactorsphere->GetProperty()->SetOpacity(0.5);
    renderer->AddActor(vtkactorsphere);*/


    vtkSmartPointer<vtkColorTransferFunction> colorannot = vtkSmartPointer<vtkColorTransferFunction>::New();
    colorannot->SetColorSpaceToHSV();

    vtkDataArray* labelarray = (vtkDataArray*)vtkDataArray::CreateArray(VTK_DOUBLE);
    //double* curvrange = curv->GetOutput()->GetPointData()->GetScalars()->GetRange();


    vnl_vector<double> northpole(3);
    northpole.fill(3);

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

        int unknownlabel = maplabels["unknown"];
        int numunkown = 0;

        for(unsigned i = 0; i < vectlabel.size();i++){
            double val[1];
            val[0] = vectlabel[i];
            labelarray->InsertTuple(i, val);

            if(val[0]==unknownlabel){
                numunkown++;
                double temp[3];
                sphere->GetPoint(i, temp);
                northpole[0] += temp[0];
                northpole[1] += temp[1];
                northpole[2] += temp[2];
            }
        }

        northpole /= numunkown;
        northpole.normalize();
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


    vnl_vector<double> north(3);
    north.fill(0);
    north[2] = 1;

    vnl_vector<double> cross = vnl_cross_3d(northpole, north);
    double angle = acos(dot_product(northpole, north))*180.0/PI;

    vtkSmartPointer<vtkTransform> transform = vtkSmartPointer<vtkTransform>::New();
    transform->RotateWXYZ(angle, cross[0], cross[1], cross[2]);

    vtkSmartPointer<vtkTransformPolyDataFilter> transformFilter =  vtkSmartPointer<vtkTransformPolyDataFilter>::New();
    transformFilter->SetTransform(transform);
    transformFilter->SetInputConnection(sphere->GetProducerPort());
    transformFilter->Update();
    sphere = transformFilter->GetOutput();


    vtkSmartPointer<vtkPointLocator> locator = vtkSmartPointer<vtkPointLocator>::New();
    locator->SetDataSet(sphere);
    locator->AutomaticOn();
    locator->SetNumberOfPointsPerBucket(3);
    locator->BuildLocator();


    vtkSmartPointer<vtkSRep> srepfig = vtkSmartPointer<vtkSRep>::New();
    vtkSmartPointer< vtkPoints > hubpos = vtkSmartPointer< vtkPoints >::New();
    vtkSmartPointer<vtkCellArray> cellarray = vtkSmartPointer<vtkCellArray>::New();

    vtkSRep::VectorSRepIdsType pointsIds;

    vtkSRep::RadiusVectorType allradius;
    vtkSRep::SpokesVectorType allspokes;

    double bounds[6];
    sphere->GetBounds(bounds);

    vector<double> steps;
    double logsum = 0;
    for(double step = 0.1; step < 0.8; step+=0.02){
        double temp = fabs(log(step));
        steps.push_back(temp);
        logsum += temp;
    }
    steps.insert(steps.end(), steps.rbegin(), steps.rend());


    double x = -logsum;
    for(unsigned xn = 0; xn < steps.size(); xn++){
        x += steps[xn];
        double y = -logsum;
        pointsIds.push_back(vtkSRep::VectorIdsType());
        for(unsigned yn = 0; yn < steps.size(); yn++){
            y += steps[yn];

            double point[3];
            double X = 4*x/logsum;
            double Y = 4*y/logsum;
            point[0] = (2*X/(1+pow(X,2)+pow(Y,2)));
            point[1] = (2*Y/(1+pow(X,2)+pow(Y,2)));
            point[2] = ((-1+pow(X,2)+pow(Y,2))/(1+pow(X,2)+pow(Y,2)));

            vtkSmartPointer<vtkIdList> pointsid = vtkSmartPointer<vtkIdList>::New();
            locator->FindClosestNPoints(6, point, pointsid);

            vtkSRep::VNLType avg(3);
            avg.fill(0);

            vtkSRep::VNLType pialavg(3);
            pialavg.fill(0);

            for(unsigned i = 0; i < pointsid->GetNumberOfIds(); i++){
                double pointinf[3];

                orig->GetPoint(pointsid->GetId(i), pointinf);
                sphere->GetPoint(pointsid->GetId(i), pointinf);
                avg[0]+= pointinf[0];
                avg[1]+= pointinf[1];
                avg[2]+= pointinf[2];


                pial->GetPoint(pointsid->GetId(i), pointinf);
                sphere->GetPoint(pointsid->GetId(i), pointinf);
                pialavg[0]+= pointinf[0];
                pialavg[1]+= pointinf[1];
                pialavg[2]+= pointinf[2];

            }

            avg /= pointsid->GetNumberOfIds();
            pialavg /= pointsid->GetNumberOfIds();

            vnl_vector<double> midpoint = (avg + pialavg)/2.0;

            pointsIds[(int)xn].push_back(hubpos->InsertNextPoint(midpoint[0], midpoint[1], midpoint[2]));

            vtkSRep::VNLType s(3);
            vtkSRep::VNLType s1(3);

            s = pialavg - midpoint;
            s1 = avg - midpoint;

            vtkSRep::VectorDoubleType radius;
            radius.push_back(s.magnitude());
            radius.push_back(s1.magnitude());

            vtkSRep::VectorVNLType vnlspokes;

            vnlspokes.push_back(s.normalize());
            vnlspokes.push_back(s1.normalize());


            allspokes.push_back(vnlspokes);
            allradius.push_back(radius);
        }

    }


    for(unsigned i = 0; i < pointsIds.size() - 1; i++){
         for(unsigned j = 0; j < pointsIds[i].size() - 1; j++){

             vtkSmartPointer<vtkQuad> quad = vtkSmartPointer<vtkQuad>::New();
             quad->GetPointIds()->SetId(0, pointsIds[i][j]);
             quad->GetPointIds()->SetId(1, pointsIds[i+1][j]);
             quad->GetPointIds()->SetId(2, pointsIds[i+1][j+1]);
             quad->GetPointIds()->SetId(3, pointsIds[i][j+1]);

             //quad->Print(cout);

             cellarray->InsertNextCell(quad);

         }
     }

    srepfig->SetPoints(hubpos);
    srepfig->SetPolys(cellarray);
    srepfig->SetAllSpokes(allspokes);
    //srepfig->SetAllSpokesRadius(allradius);
    //srepfig->SetGridTopolgyIds(pointsIds);


    int interpolationlevel = 3;


    vtkSmartPointer< vtkSRepInterpolateMedialSheet > medialsheetinterpolator = vtkSmartPointer< vtkSRepInterpolateMedialSheet >::New();
    medialsheetinterpolator->SetInput(srepfig);
    medialsheetinterpolator->SetInterpolationLevel(interpolationlevel);
    medialsheetinterpolator->Update();


    vtkSmartPointer< vtkPolyData > medialsheet = medialsheetinterpolator->GetOutput();

    vtkSmartPointer<vtkPolyDataMapper> medialsheetmapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    medialsheetmapper->SetInputConnection(medialsheet->GetProducerPort());
    vtkSmartPointer<vtkActor>  medialsheetactor = vtkActor::New();
    medialsheetactor->SetMapper(medialsheetmapper);
    //medialsheetactor->GetProperty()->SetOpacity(0.5);
    medialsheetactor->GetProperty()->SetPointSize(10);
    medialsheetactor->GetProperty()->SetColor(0.2,0.6,1);
    medialsheetactor->GetProperty()->SetRepresentationToWireframe();
    renderer->AddActor(medialsheetactor);




    vtkSmartPointer< vtkSRepVisuPrimitives > visuprimitives = vtkSmartPointer< vtkSRepVisuPrimitives >::New();
    visuprimitives->SetInput(srepfig);
    visuprimitives->Update();




    vtkActorCollection* actorcollection = visuprimitives->GetOuputActors();

    for(unsigned i = 0; i < 4; i++){

        vtkActor* actor = (vtkActor*) actorcollection->GetItemAsObject(i);

        renderer->AddActor(actor);
    }



    renderWindow->Render();

    renderWindowInteractor->Start();


    return 0;
}

