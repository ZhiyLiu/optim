#include <vtkSmartPointer.h>

#include <vtksrep.h>
#include <vtksrepinterpolatemedialsheet.h>
#include <vtksrepvisuprimitives.h>

#include <vtkCellArray.h>
#include <vtkQuad.h>
#include <vtkLine.h>

#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkProperty.h>

#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderWindow.h>

#include <vtksrepinterpolatemedialcrestcurve.h>
#include <vtksrepinterpolatemedialspokes.h>
#include <vtksrepinterpolatecrestspokes.h>
#include <vtksrepinterpolatecrestspokesquartic.h>

#include "P3DControl.h"
#include "ControlParms.h"

ControlParms * globalControl;	// Read the user's preferences file
int globalVerbosity;			// Current verbosity level of Pablo

using namespace std;

void help(char* execname){
    cout<<"Srep test using vtk"<<endl;
    //cout<<"Usage: "<<execname<<" -d <patients dir> -p <patient name> [options]"<<endl;
    //cout<<"options:"<<endl;
    cout<<"--h --help show help menu"<<endl;
    cout<<"-m <model name> ";

}

M3DTubeFigure* GetTubeFigure(const char* figfilename){
    globalControl = new ControlParms(NULL, -1, false);	// Ignore user's preferences
    globalVerbosity = 0;
    globalControl->setDefault(OutputVerbosity, globalVerbosity);
    globalControl->setDefault(ReorderModels, 0);
    globalControl->setDefault(SmoothImages, false);
    globalControl->setDefault(ConvertImages, false);
    globalControl->setDefault(ByteOrder, 1);
    globalControl->setDefault(CompressImages, true);
    globalControl->setDefault(ShowLandmarks, true);

    P3DControl* p3d = new P3DControl(10);
    p3d->read(figfilename, false);
    M3DObject* m3dobject = p3d->getObjectPtr();
    M3DFigure * figure = m3dobject->getFigurePtr( 0 ) ;
   return dynamic_cast<M3DTubeFigure*>( figure );
}

int main(int argc, char *argv[])
{


    if (argc < 1){
        help(argv[0]);
        return 0;
    }


    vector< std::string > modelname;

    for(int i = 1; i < argc; i++){

        if(string(argv[i]) == "--h" || string(argv[i]) == "--help"){
            help(argv[0]);
            return 0;
        }else if(string(argv[i]) == "-m"){
            modelname.push_back(argv[i + 1]);
        }

    }

         vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
         renderer->SetBackground(1,1,1);
         vtkSmartPointer<vtkRenderWindow> renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
         renderWindow->AddRenderer(renderer);
         vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
         renderWindowInteractor->SetRenderWindow(renderWindow);



         M3DTubeFigure* tubefig = GetTubeFigure(modelname[0].c_str());



         vtkSmartPointer<vtkSRep> srepfig = vtkSmartPointer<vtkSRep>::New();
             vtkSmartPointer< vtkPoints > hubpos = vtkSmartPointer< vtkPoints >::New();
             vtkSmartPointer<vtkCellArray> cellarray = vtkSmartPointer<vtkCellArray>::New();

             vtkSRep::VectorSRepIdsType pointsIds;

             vtkSRep::RadiusVectorType allradius;
             vtkSRep::SpokesVectorType allspokes;



             pointsIds.push_back(vtkSRep::VectorIdsType());

             for(int u = 0; u < (int)tubefig->getPrimitiveCount(); u++){

                 M3DTubePrimitive* prim0 = dynamic_cast<M3DTubePrimitive*>(tubefig->getPrimitivePtr(u));
                 Vector3D x = prim0->getX();

                 pointsIds[0].push_back(hubpos->InsertNextPoint(x.getX(), x.getY(), x.getZ()));

                 vtkSRep::VectorVNLType vnlspokes;
                 vtkSRep::VectorDoubleType radius;

                 for(int v = 0; v < (int)tubefig->getNumberOfSpokes(); v++){


                     Vector3D u0 = prim0->getYN(v);


                     vtkSRep::VNLType s(3);
                     s[0] = u0.getX();
                     s[1] = u0.getY();
                     s[2] = u0.getZ();

                     radius.push_back(s.magnitude());
                     vnlspokes.push_back(s.normalize());


                 }

                 allspokes.push_back(vnlspokes);
                 allradius.push_back(radius);
             }


             for(unsigned i = 0; i < pointsIds.size(); i++){
                  for(unsigned j = 0; j < pointsIds[i].size() - 1; j++){

                      vtkSmartPointer<vtkLine> line = vtkSmartPointer<vtkLine>::New();
                      line->GetPointIds()->SetId(0, pointsIds[i][j]);
                      line->GetPointIds()->SetId(1, pointsIds[i][j+1]);

                      //quad->Print(cout);

                      cellarray->InsertNextCell(line);

                  }
              }

             srepfig->SetPoints(hubpos);
             srepfig->SetLines(cellarray);
             srepfig->SetAllSpokes(allspokes);
             srepfig->SetAllSpokesRadius(allradius);
             srepfig->SetNumberOfSpokes(tubefig->getNumberOfSpokes());
             //srepfig->SetGridTopolgyIds(pointsIds);



         int interpolationlevel = 3;




         vtkSmartPointer<vtkSRepInterpolateMedialCrestCurve> curveinterpolation = vtkSmartPointer<vtkSRepInterpolateMedialCrestCurve>::New();
         curveinterpolation->SetInput(srepfig);
         curveinterpolation->SetInterpolationLevel(interpolationlevel);
         curveinterpolation->Update();
         vtkPolyData* medialcrestcurve = curveinterpolation->GetOutput();

         vtkSmartPointer<vtkPolyDataMapper> medialcrestcurvemapper = vtkSmartPointer<vtkPolyDataMapper>::New();
         medialcrestcurvemapper->SetInputConnection(medialcrestcurve->GetProducerPort());
         vtkSmartPointer<vtkActor>  medialsheetcrestcurveactor = vtkActor::New();
         medialsheetcrestcurveactor->SetMapper(medialcrestcurvemapper);
         medialsheetcrestcurveactor->GetProperty()->SetLineWidth(5);
         medialsheetcrestcurveactor->GetProperty()->SetColor(1,1,0);
         renderer->AddActor(medialsheetcrestcurveactor);


         //vtkSRep::VectorIdsType crestids = hippoc->GetCrestMedialAtomsIds();

         //int n = 0;

         //while(n < 100){
         //    n++;
          //   for(unsigned i = 0; i < crestids.size(); i++){
                 //vtkSmartPointer<vtkSRepInterpolateCrestSpokes> interpolatecrestspokes = vtkSmartPointer<vtkSRepInterpolateCrestSpokes>::New();
                 vtkSmartPointer<vtkSRepInterpolateCrestSpokesQuartic> interpolatecrestspokes = vtkSmartPointer<vtkSRepInterpolateCrestSpokesQuartic>::New();
                 interpolatecrestspokes->SetInput(srepfig);
                 interpolatecrestspokes->SetInterpolationLevel(interpolationlevel);
                 interpolatecrestspokes->SetUseAllSpokes(true);
                 interpolatecrestspokes->SetCyclicSpokes(true);
                 //interpolatecrestspokes->SetAtomId(crestids[i]);
                 //interpolatecrestspokes->SetAtomId(8);
                 //interpolatecrestspokes->SetGamma_t(0.5);
                 //interpolatecrestspokes->SetGamma_theta(0.5);
                 //interpolatecrestspokes->SetSpokeType(vtkSRep::TOP_SPOKE);
                 interpolatecrestspokes->Update();




                 vtkPolyData* interpolatedcrest = interpolatecrestspokes->GetOutput();

                 vtkSmartPointer<vtkPolyDataMapper> crestspokescurvemapper = vtkSmartPointer<vtkPolyDataMapper>::New();
                 crestspokescurvemapper->SetInputConnection(interpolatedcrest->GetProducerPort());
                 vtkSmartPointer<vtkActor>  crestspokesactor = vtkActor::New();
                 crestspokesactor->SetMapper(crestspokescurvemapper);
                 crestspokesactor->GetProperty()->SetLineWidth(5);
                 crestspokesactor->GetProperty()->SetColor(0.2,0.8,0.3);
                 renderer->AddActor(crestspokesactor);




                 vtkSmartPointer<vtkSRep> srepcrest = interpolatecrestspokes->GetSRepOutput();

                 vtkSmartPointer<vtkPolyData> polycrestspokes = vtkSmartPointer<vtkPolyData>::New();
                 vtkSmartPointer<vtkCellArray> cellarraycrestspokes = vtkSmartPointer<vtkCellArray>::New();
                 vtkSmartPointer<vtkPoints> pointscrestspokes  = vtkSmartPointer<vtkPoints>::New();

                 for(unsigned i = 0; i < srepcrest->GetNumberOfPoints(); i++){

                     double point[3];

                     srepcrest->GetPoint(i, point);

                     vtkIdType id0 = pointscrestspokes->InsertNextPoint(point[0], point[1], point[2]);

                     vtkSRep::VectorVNLType currentspokes = srepcrest->GetSpokes(i);
                     vtkSRep::VectorDoubleType radius = srepcrest->GetSpokesRadius(i);

                     for(unsigned j = 0; j < currentspokes.size(); j++){

                         vtkSRep::VNLType p1 = currentspokes[j]*radius[j];

                         vtkSmartPointer<vtkLine> crestspokeline = vtkSmartPointer<vtkLine>::New();
                         crestspokeline->GetPointIds()->SetId(0, id0);
                         crestspokeline->GetPointIds()->SetId(1, pointscrestspokes->InsertNextPoint(point[0] + p1[0], point[1] + p1[1], point[2] + p1[2]));

                         cellarraycrestspokes->InsertNextCell(crestspokeline);

                     }

                 }

                 polycrestspokes->SetPoints(pointscrestspokes);
                 polycrestspokes->SetLines(cellarraycrestspokes);

                 vtkSmartPointer<vtkPolyDataMapper> crestspokesmapper = vtkSmartPointer<vtkPolyDataMapper>::New();
                 crestspokesmapper->SetInput(polycrestspokes);
                 vtkSmartPointer<vtkActor>  spokesactor = vtkActor::New();
                 spokesactor->SetMapper(crestspokesmapper);
                 spokesactor->GetProperty()->SetLineWidth(5);
                 spokesactor->GetProperty()->SetColor(0,1,0);
                 //renderer->AddActor(spokesactor);
                 //renderWindow->Render();


                 //renderer->RemoveActor(spokesactor);

           // }
         //}



                 vtkSmartPointer< vtkSRepVisuPrimitives > visuprimitives = vtkSmartPointer< vtkSRepVisuPrimitives >::New();
                 visuprimitives->SetInput(srepfig);
                 visuprimitives->Update();




                 vtkActorCollection* actorcollection = visuprimitives->GetOuputActors();

                 for(unsigned i = 0; i < 5; i++){

                     vtkActor* actor = (vtkActor*) actorcollection->GetItemAsObject(i);

                     renderer->AddActor(actor);
                 }




         renderWindow->Render();
         renderWindowInteractor->Start();



    return 0;
}




