#include "vtkwritesrep.h"

#include "vtkQuad.h"

#include "vtkObjectFactory.h"
#include "vnl/vnl_cross.h"

vtkStandardNewMacro(vtkWriteSRep);

vtkWriteSRep::vtkWriteSRep()
{

    globalControl = 0;	// Read the user's preferences file
    globalVerbosity = 0;			// Current verbosity level of Pablo
    p3d = 0;

    m_Srep = 0;

}

vtkWriteSRep::~vtkWriteSRep()
{

    if(globalControl){
        delete globalControl;
        globalControl = 0;
    }

    if(p3d){
        delete p3d;
        p3d = 0;
    }

}


int vtkWriteSRep::Write(){


    vtkSmartPointer<vtkSRep> srepfig = dynamic_cast<vtkSRep*>(vtkSRep::SafeDownCast(this->GetInput()));

    M3DFigure* m3dfig = NewQuadFigure(srepfig->GetNumRows(), srepfig->GetNumColumns());

    double* c = srepfig->GetColor();
    float cf[3];
    cf[0] = c[0];
    cf[1] = c[1];
    cf[2] = c[2];
    m3dfig->setColor(cf);

    for(unsigned l = 0; l < srepfig->GetNumRows(); l++){
         for(unsigned m = 0; m < srepfig->GetNumColumns(); m++){

             double point[3];
             vtkIdType srepspokeid = l*srepfig->GetNumColumns() + m;
             srepfig->GetPoint(srepspokeid, point);

             if(srepfig->GetNumRows() > 1){

                 M3DQuadFigure* quadfig = dynamic_cast<M3DQuadFigure*>(m3dfig);
                 M3DQuadPrimitive* prim = dynamic_cast<M3DQuadPrimitive*>(quadfig->getPrimitivePtr(l, m));
                 prim->setX(point[0], point[1], point[2]);

                 vtkSRep::VNLType vnlvect0 = srepfig->GetSpoke(srepspokeid, vtkSRep::TOP_SPOKE);
                 double r0 = vnlvect0.magnitude();
                 vnlvect0.normalize();
                 Vector3D vect0 = prim->getU0();
                 vect0.setX(vnlvect0[0]);
                 vect0.setY(vnlvect0[1]);
                 vect0.setZ(vnlvect0[2]);
                 prim->setU0(vect0);
                 prim->setR(0, r0);


                 vtkSRep::VNLType vnlvect1 = srepfig->GetSpoke(srepspokeid, vtkSRep::BOTTOM_SPOKE);
                 double r1 = vnlvect1.magnitude();
                 vnlvect1.normalize();

                 Vector3D vect1 = prim->getU1();
                 vect1.setX(vnlvect1[0]);
                 vect1.setY(vnlvect1[1]);
                 vect1.setZ(vnlvect1[2]);
                 prim->setU1(vect1);
                 prim->setR(1, r1);

                 if(srepfig->GetSpokes(srepspokeid).size() == 3){
                     vtkSRep::VNLType vnlvectcrest = srepfig->GetSpoke(srepspokeid, vtkSRep::CREST_SPOKE);
                     double rend = vnlvectcrest.magnitude();
                     vnlvectcrest.normalize();
                     Vector3D vectend = prim->getUEnd();
                     vectend.setX(vnlvectcrest[0]);
                     vectend.setY(vnlvectcrest[1]);
                     vectend.setZ(vnlvectcrest[2]);
                     prim->setUEnd(vectend);
                     prim->setR(2, rend);
                 }
             }else{
                 M3DTubeFigure* tubefig = dynamic_cast<M3DTubeFigure*>(m3dfig);
                 tubefig->setPositivePolarity(true);
                 tubefig->setPositiveSpace(true);
                 tubefig->setTolerance(50);

                 double *r = new double[srepfig->GetSpokes(srepspokeid).size()];
                 Vector3D vect0;
                 Vector3D vect2;

                 for(unsigned j = 0; j < srepfig->GetSpokes(srepspokeid).size(); j++){

                     vtkSRep::VNLType u0 = (srepfig->GetSpoke(srepspokeid, j))*-1;
                     r[j] = u0.magnitude();

                     if(j == 0){
                         u0.normalize();

                         vect0.setX(u0[0]);
                         vect0.setY(u0[1]);
                         vect0.setZ(u0[2]);

                         vtkSRep::VNLType u1 = srepfig->GetMedialSheetNormal(srepspokeid);

                         u1 = vnl_cross_3d(u0, u1);

                         vect2.setX(u1[0]);
                         vect2.setY(u1[1]);
                         vect2.setZ(u1[2]);

                     }

                 }

                 M3DTubePrimitive* prim = 0;
                 prim = new M3DTubePrimitive(point[0], point[1], point[2], 0, r, vect0, vect2, srepfig->GetSpokes(srepspokeid).size());
                 prim->setBaseAtom(false);
                 prim->select();


                 tubefig->setPrimitivePtr(m, prim);
             }
         }
     }

    p3d->write(this->GetFileName());
}
M3DFigure* vtkWriteSRep::NewQuadFigure(const int numr, int numc){
    globalControl = new ControlParms(NULL, -1, false);	// Ignore user's preferences
    globalVerbosity = 0;
    globalControl->setDefault(OutputVerbosity, globalVerbosity);
    globalControl->setDefault(ReorderModels, 0);
    globalControl->setDefault(SmoothImages, false);
    globalControl->setDefault(ConvertImages, false);
    globalControl->setDefault(ByteOrder, 1);
    globalControl->setDefault(CompressImages, true);
    globalControl->setDefault(ShowLandmarks, true);

    p3d = new P3DControl(10);
    float color[3] = {0.2,0.5,0.8};    
    p3d->addQuadFigure(numr, numc, "Cortex", color);    

    M3DObject* m3dobject = p3d->getObjectPtr();
    M3DFigure * figure = m3dobject->getFigurePtr( 0 ) ;
    //figure->setLandmarkName(0,0)
    return figure;
}
