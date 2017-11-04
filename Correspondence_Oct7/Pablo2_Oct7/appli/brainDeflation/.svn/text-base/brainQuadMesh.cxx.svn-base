


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


    vtkSmartPointer<vtkPolyDataReader> reader3 =  vtkSmartPointer<vtkPolyDataReader>::New();

    string spherefilename = dirname;
    spherefilename.append("/");
    spherefilename.append(patientname[0]);
    spherefilename.append("/surf/vtk/lh.sphere.vtk");

    reader3->SetFileName(spherefilename.c_str());
    reader3->Update();
    vtkPolyData* sphere = reader3->GetOutput();


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


        for ( hue = 0; hue < 1; hue+=step ){
            colorthick->AddHSVPoint(hue*5, hue, 1, 1);
        }
    }

    colorthick->Build();

    vtkDataArray* normalDataArray =  inflated->GetPointData()->GetArray("Normals");

    vtkSmartPointer<vtkPolyDataNormals> normalGenerator = vtkSmartPointer<vtkPolyDataNormals>::New();

    cout<<inflated->GetNumberOfPoints()<< endl;

    vtkPolyData* inflatednormals = 0;



    if(!normalDataArray ){
        //normalGenerator->SetInput(inflated);
        normalGenerator->SetInput(orig);
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
    inflatedpial->BuildLinks();


    double bounds[6];
    inflated->GetBounds(bounds);
    double origin[3] ;
    origin[0] = (bounds[1] + bounds[0])/2.0;
    origin[1] = (bounds[3] + bounds[2])/2.0;
    origin[2] = (bounds[5] + bounds[4])/2.0;


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

    double step = PI/80;
    for(double phi = 0, n = 0; phi < 2*PI; phi+=step, n++){
        pointsIds.push_back(vtkSRep::VectorIdsType());
        for(double theta = 0; theta < 2*PI; theta+=step){

            double point[3];
            point[0] = 100*sin(theta)*cos(phi);
            point[1] = 100*sin(theta)*sin(phi);
            point[2] = 100*cos(theta);

            vtkSmartPointer<vtkIdList> pointsid = vtkSmartPointer<vtkIdList>::New();
            locator->FindClosestNPoints(9, point, pointsid);

            vtkSRep::VNLType avg(3);
            avg.fill(0);

            vtkSRep::VNLType normal(3);
            normal.fill(0);

            double thickness = 0;

            for(unsigned i = 0; i < pointsid->GetNumberOfIds(); i++){
                double pointinf[3];
                //inflated->GetPoint(pointsid->GetId(i), pointinf);
                orig->GetPoint(pointsid->GetId(i), pointinf);
                avg[0]+= pointinf[0];
                avg[1]+= pointinf[1];
                avg[2]+= pointinf[2];


                double normalinf[3];
                normalDataArray->GetTuple(pointsid->GetId(i), normalinf);
                normal[0] += normalinf[0];
                normal[1] += normalinf[1];
                normal[2] += normalinf[2];

                thickness += vectthick[pointsid->GetId(i)];
            }

            avg /= pointsid->GetNumberOfIds();
            normal /= pointsid->GetNumberOfIds();
            normal = normal.normalize();
            thickness /= pointsid->GetNumberOfIds();
            thickness/=2.0;

            avg = avg + normal*thickness;


            vtkSRep::VectorVNLType vnlspokes;
            vtkSRep::VNLType s(3);
            s = thickness*normal;
            s = s.normalize();
            vnlspokes.push_back(s);

            vtkSRep::VNLType s1(3);
            s1 = -thickness*normal;
            s1 = s1.normalize();
            vnlspokes.push_back(s1);

            vtkSRep::VectorDoubleType radius;
            radius.push_back(thickness);
            radius.push_back(thickness);

            pointsIds[(int)n].push_back(hubpos->InsertNextPoint(avg[0], avg[1], avg[2]));
            allspokes.push_back(vnlspokes);
            allradius.push_back(radius);


            /*vtkSmartPointer<vtkSphereSource> sphereSource = vtkSmartPointer<vtkSphereSource>::New();
            sphereSource->SetCenter(avg[0], avg[1], avg[2]);
             sphereSource->SetRadius(1.0);

             vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
             mapper->SetInputConnection(sphereSource->GetOutputPort());

             vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
             actor->SetMapper(mapper);
             renderer->AddActor(actor);*/
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
    srepfig->SetAllSpokesRadius(allradius);
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
    renderer->AddActor(medialsheetactor);




    vtkSmartPointer< vtkSRepVisuPrimitives > visuprimitives = vtkSmartPointer< vtkSRepVisuPrimitives >::New();
    visuprimitives->SetInput(srepfig);
    visuprimitives->Update();




    vtkActorCollection* actorcollection = visuprimitives->GetOuputActors();

    for(unsigned i = 1; i < 4; i++){

        vtkActor* actor = (vtkActor*) actorcollection->GetItemAsObject(i);

        renderer->AddActor(actor);
    }






    vtkSmartPointer<vtkPolyDataMapper> inflatedpialmapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    inflatedpialmapper->SetInput(inflatedpial);
    //inflatedpialmapper->SetInput(inflatedpialquad);
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
    //renderer->AddActor(vtkactorinflated);

    //vtkactorinflatedpial->GetProperty()->SetOpacity(0.2);
    vtkactorinflatedpial->GetProperty()->SetRepresentationToWireframe();
    //renderer->AddActor(vtkactorinflatedpial);

    //renderer->AddActor(vtkactorpial);


    //orig->GetPointData()->SetScalars(curvaturearray);
    //origmapper->SetLookupTable(color);

    orig->GetPointData()->SetScalars(labelarray);
    origmapper->SetLookupTable(colorannot);

    //inflated->GetPointData()->SetScalars(curvaturearray);
    //inflatedmapper->SetLookupTable(color);

    inflated->GetPointData()->SetScalars(labelarray);
    inflatedmapper->SetLookupTable(colorannot);
    //inflated->GetPointData()->SetScalars(thickarray);
    //inflatedmapper->SetLookupTable(colorthick);

    //inflatedpial->GetPointData()->SetScalars(labelarray);
    //inflatedpialmapper->SetLookupTable(colorannot);

    pial->GetPointData()->SetScalars(labelarray);
    pialmapper->SetLookupTable(colorannot);

     //vtkactorinflated->GetProperty()->SetOpacity(0.2);

    renderWindow->Render();

    renderWindowInteractor->Start();



    vtkSmartPointer<vtkPolyDataWriter> polywriter = vtkSmartPointer<vtkPolyDataWriter>::New();
    polywriter->SetInput(inflated);
    polywriter->SetFileName("inflated.vtk");
    polywriter->Write();

    polywriter->SetInput(inflatedpial);
    polywriter->SetFileName("inflatedpial.vtk");
    polywriter->Write();



    return 0;
}


/*vtkSmartPointer<vtkImageData> inflatedimage = vtkSmartPointer<vtkImageData>::New();


PolyDataImageStencilExport(inflated, inflatedimage, inflated->GetBounds());

vtkSmartPointer<vtkImageData> inflatedimagequad = vtkSmartPointer<vtkImageData>::New();
inflatedimagequad->DeepCopy(inflatedimage);


int* extent = inflatedimage->GetExtent();
vector< int* > offset;
int *a = new int[3];
a[0] = -1;
a[1] = 0;
a[2] = 0;
offset.push_back(a);
a = new int[3];
a[0] = 1;
a[1] = 0;
a[2] = 0;
offset.push_back(a);
a = new int[3];
a[0] = 0;
a[1] = -1;
a[2] = 0;
offset.push_back(a);
a = new int[3];
a[0] = 0;
a[1] = 1;
a[2] = 0;
offset.push_back(a);
a = new int[3];
a[0] = 0;
a[1] = 0;
a[2] = -1;
offset.push_back(a);
a = new int[3];
a[0] = 0;
a[1] = 0;
a[2] = 1;
offset.push_back(a);


for (int z  =extent[4] + 1; z < extent[5]; z++){
    for (int x = extent[0] + 1; x < extent[1]; x++){
        for (int y = extent[2] + 1; y < extent[3]; y++){

            unsigned short* inflatedpixel = static_cast<unsigned short*>(inflatedimage->GetScalarPointer(x,y,z));
            unsigned short* inflatedquadpixel = static_cast<unsigned short*>(inflatedimagequad->GetScalarPointer(x,y,z));

            *inflatedquadpixel = 0;

            if(*inflatedpixel != 0 ){//&& z%2 == 0 && x%2==0 && y%2==0){

                bool border = false;

                for(unsigned i = 0; i < offset.size() && !border; i++){

                    unsigned short* borderpixel = static_cast<unsigned short*>(inflatedimage->GetScalarPointer(offset[i][0] + x, offset[i][1] + y, offset[i][2] + z));
                    if(*borderpixel == 0){
                        border = true;
                    }

                }


                if(border){
                    *inflatedquadpixel = 128;
                }
            }
        }
    }
}

for (int z  =extent[4] + 1; z < extent[5]; z++){
    for (int x = extent[0] + 1; x < extent[1]; x++){
        for (int y = extent[2] + 1; y < extent[3]; y++){

            unsigned short* inflatedquadpixel = static_cast<unsigned short*>(inflatedimagequad->GetScalarPointer(x,y,z));

            if(x%4 != 0 || y%4 != 0){
                *inflatedquadpixel = 0;
            }
        }
    }
}



vtkSmartPointer<vtkMetaImageWriter> imagewriter = vtkSmartPointer<vtkMetaImageWriter>::New();
imagewriter->SetInput(inflatedimagequad);
imagewriter->SetFileName("output1.mhd");
imagewriter->Write();*/


/*vtkSmartPointer<vtkQuadricClustering> quadriccluster = vtkSmartPointer<vtkQuadricClustering>::New();
quadriccluster->SetInput(inflatedpial);
//quadriccluster->SetAutoAdjustNumberOfDivisions(1);
quadriccluster->Update();

inflatedpial = quadriccluster->GetOutput();





vtkSmartPointer<vtkPolyData> inflatedpialquad = vtkSmartPointer<vtkPolyData>::New();
vtkSmartPointer<vtkPoints> inflatedpialquadpoints = vtkSmartPointer<vtkPoints>::New();
vtkSmartPointer<vtkCellArray> inflatedpialquadcells = vtkSmartPointer<vtkCellArray>::New();

map<vtkIdType, vtkIdType> pointindex;

for(unsigned i = 0; i < inflatedpial->GetNumberOfCells(); i++){

    vtkSmartPointer<vtkIdList> pointsid = vtkSmartPointer<vtkIdList>::New();
    inflatedpial->GetCellPoints(i, pointsid);


    double point0[3];
    inflatedpial->GetPoint(pointsid->GetId(0), point0);

    double point1[3];
    inflatedpial->GetPoint(pointsid->GetId(1), point1);

    double point2[3];
    inflatedpial->GetPoint(pointsid->GetId(2), point2);

    double avg[3];
    avg[0]= (point0[0] + point1[0] + point2[0])/3.0;
    avg[1]= (point0[1] + point1[1] + point2[1])/3.0;
    avg[2]= (point0[2] + point1[2] + point2[2])/3.0;

    double avg01[3];
    avg01[0]= (point0[0] + point1[0])/2.0;
    avg01[1]= (point0[1] + point1[1])/2.0;
    avg01[2]= (point0[2] + point1[2])/2.0;

    double avg02[3];
    avg02[0]= (point0[0] + point2[0])/2.0;
    avg02[1]= (point0[1] + point2[1])/2.0;
    avg02[2]= (point0[2] + point2[2])/2.0;

    double avg12[3];
    avg12[0]= (point1[0] + point2[0])/2.0;
    avg12[1]= (point1[1] + point2[1])/2.0;
    avg12[2]= (point1[2] + point2[2])/2.0;


    vtkIdType id0 = pointsid->GetId(0);
    if(pointindex.find(id0) == pointindex.end()){
        pointindex[id0] = inflatedpialquadpoints->InsertNextPoint(point0[0], point0[1], point0[2]);
    }

    vtkIdType id1 = pointsid->GetId(1);
    if(pointindex.find(id1) == pointindex.end()){
        pointindex[id1] = inflatedpialquadpoints->InsertNextPoint(point1[0], point1[1], point2[2]);
    }

    vtkIdType id2 = pointsid->GetId(2);
    if(pointindex.find(id2) == pointindex.end()){
        pointindex[id2] = inflatedpialquadpoints->InsertNextPoint(point2[0], point2[1], point2[2]);
    }

    vtkIdType id01 = inflatedpialquadpoints->InsertNextPoint(avg01[0], avg01[1], avg01[2]);
    vtkIdType id02 = inflatedpialquadpoints->InsertNextPoint(avg02[0], avg02[1], avg02[2]);
    vtkIdType id12 = inflatedpialquadpoints->InsertNextPoint(avg12[0], avg12[1], avg12[2]);
    vtkIdType idavg = inflatedpialquadpoints->InsertNextPoint(avg[0], avg[1], avg[2]);

    vtkQuad* quad0 = vtkQuad::New();
    quad0->GetPointIds()->SetId(0, pointindex[id0]);
    quad0->GetPointIds()->SetId(1, id01);
    quad0->GetPointIds()->SetId(2, idavg);
    quad0->GetPointIds()->SetId(3, id02);


    vtkQuad* quad1 = vtkQuad::New();
    quad1->GetPointIds()->SetId(0, id01);
    quad1->GetPointIds()->SetId(1, pointindex[id1]);
    quad1->GetPointIds()->SetId(2, id12);
    quad1->GetPointIds()->SetId(3, idavg);


    vtkQuad* quad2 = vtkQuad::New();
    quad2->GetPointIds()->SetId(0, id02);
    quad2->GetPointIds()->SetId(1, idavg);
    quad2->GetPointIds()->SetId(2, id12);
    quad2->GetPointIds()->SetId(3, pointindex[id2]);

    inflatedpialquadcells->InsertNextCell(quad0);
    inflatedpialquadcells->InsertNextCell(quad1);
    inflatedpialquadcells->InsertNextCell(quad2);


}


inflatedpialquad->SetPoints(inflatedpialquadpoints);
inflatedpialquad->SetPolys(inflatedpialquadcells);*/
