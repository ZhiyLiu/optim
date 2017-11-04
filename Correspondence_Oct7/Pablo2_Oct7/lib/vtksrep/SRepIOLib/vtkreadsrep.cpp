#include "vtkreadsrep.h"

#include "vtkQuad.h"
#include "vtkLine.h"

#include "vtkTransform.h"
#include "vtkObjectFactory.h"
vtkStandardNewMacro(vtkReadSRep);

#include "vtkMatrix4x4.h"
#include "vnl/vnl_cross.h"
#include "vtkPointData.h"

//ControlParms * globalControl;	// Read the user's preferences file

vtkReadSRep::vtkReadSRep()
{

    globalControl = 0;	// Read the user's preferences file
    globalVerbosity = 0;			// Current verbosity level of Pablo
    p3d = 0;

    m_Srep = 0;

}

vtkReadSRep::~vtkReadSRep()
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


void vtkReadSRep::Update(){


    M3DFigure* figure = GetFigure();
    M3DQuadFigure* quadfig = dynamic_cast<M3DQuadFigure*>(figure);

    m_Srep = vtkSmartPointer<vtkSRep>::New();

    vtkSmartPointer< vtkPoints > hubpos = vtkSmartPointer< vtkPoints >::New();
    vtkSmartPointer<vtkCellArray> cellarray = vtkSmartPointer<vtkCellArray>::New();

    vtkSRep::VectorSRepIdsType pointsIds;

    vtkSRep::RadiusVectorType allradius;
    vtkSRep::SpokesVectorType allspokes;

    if(quadfig){

        for(int u = 0; u < quadfig->getRowCount(); u++){
            pointsIds.push_back(vtkSRep::VectorIdsType());
            for(int v = 0; v < quadfig->getColumnCount(); v++){

                M3DQuadPrimitive* prim0 = dynamic_cast<M3DQuadPrimitive*>(quadfig->getPrimitivePtr(u, v));

                Vector3D x = prim0->getX();
                Vector3D u0 = prim0->getU0();
                Vector3D u1 = prim0->getU1();

                vtkSRep::VectorVNLType vnlspokes;
                vtkSRep::VNLType s(3);
                s[0] = u0.getX();
                s[1] = u0.getY();
                s[2] = u0.getZ();
                s = s.normalize();
                vnlspokes.push_back(s);

                s[0] = u1.getX();
                s[1] = u1.getY();
                s[2] = u1.getZ();
                s = s.normalize();
                vnlspokes.push_back(s);

                vtkSRep::VectorDoubleType radius;
                radius.push_back(prim0->getR0());
                radius.push_back(prim0->getR1());

                if(u == 0 || u == quadfig->getRowCount() - 1 || v == 0 || v == quadfig->getColumnCount() - 1){

                    M3DQuadEndPrimitive* prim0 = dynamic_cast<M3DQuadEndPrimitive*>(quadfig->getPrimitivePtr(u, v));
                    Vector3D uend = prim0->getUEnd();

                    s[0] = uend.getX();
                    s[1] = uend.getY();
                    s[2] = uend.getZ();
                    s = s.normalize();
                    vnlspokes.push_back(s);

                    radius.push_back(prim0->getREnd());

                }


                vtkIdType id = hubpos->InsertNextPoint(x.getX(), x.getY(), x.getZ());
                //cout<<quadfig->getPrimitiveID(u,v)<<"-"<<id<<endl;
                pointsIds[u].push_back(id);



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

        m_Srep->SetPoints(hubpos);
        m_Srep->SetPolys(cellarray);
        m_Srep->SetAllSpokes(allspokes);
        m_Srep->SetAllRadius(allradius);
        const float *color = quadfig->getColor();

        m_Srep->SetColor(color[0], color[1], color[2]);

        allspokes.clear();
        allradius.clear();

        delete quadfig;
    }else{

        M3DTubeFigure* tubefig = dynamic_cast<M3DTubeFigure*>( figure );


        vtkSmartPointer<vtkDataArray> normals = (vtkDataArray*)vtkDataArray::CreateArray(VTK_DOUBLE);
        normals->SetNumberOfComponents(3);


        pointsIds.push_back(vtkSRep::VectorIdsType());

        for(int u = 0; u < (int)tubefig->getPrimitiveCount(); u++){

            M3DTubePrimitive* prim0 = dynamic_cast<M3DTubePrimitive*>(tubefig->getPrimitivePtr(u));
            Vector3D x = prim0->getX();

            pointsIds[0].push_back(hubpos->InsertNextPoint(x.getX(), x.getY(), x.getZ()));

            vtkSRep::VectorVNLType vnlspokes;
            vtkSRep::VectorDoubleType radius;

            vtkSRep::VNLType u0vnl(3);
            vtkSRep::VNLType u1vnl(3);
            /*vtkSRep::VNLType u2vnl(3);
            vtkSmartPointer<vtkTransform> transform = vtkSmartPointer<vtkTransform>::New();*/

            int numspokes = (int)tubefig->getNumberOfSpokes();

            for(int v = 0; v < numspokes; v++){

                Vector3D u0 = prim0->getYN(v);

                vtkSRep::VNLType s(3);
                s[0] = u0.getX();
                s[1] = u0.getY();
                s[2] = u0.getZ();

                radius.push_back(s.magnitude());
                vnlspokes.push_back(s.normalize());

                if( v == 0){
                    Vector3D u0 = prim0->getU0();
                    u0vnl[0] = u0.getX();
                    u0vnl[1] = u0.getY();
                    u0vnl[2] = u0.getZ();

                    Vector3D u2 = prim0->getB();
                    u1vnl[0] = u2.getX();
                    u1vnl[1] = u2.getY();
                    u1vnl[2] = u2.getZ();

                    u1vnl = vnl_cross_3d(u1vnl, u0vnl);
                    //cout<<u1vnl<<endl;
                    normals->InsertNextTuple3(u1vnl[0], u1vnl[1], u1vnl[2]);
                }

                /*vtkSRep::VNLType s;

                if( v == 0){
                    Vector3D u0 = prim0->getU0();
                    u0vnl[0] = u0.getX();
                    u0vnl[1] = u0.getY();
                    u0vnl[2] = u0.getZ();

                    //cout<<u0vnl<<endl;
                    //cout<<u0vnl.normalize()<<endl;

                    s = u0vnl*prim0->getRN(v);

                    Vector3D u2 = prim0->getB();
                    u1vnl[0] = u2.getX();
                    u1vnl[1] = u2.getY();
                    u1vnl[2] = u2.getZ();

                    u1vnl = vnl_cross_3d(u1vnl, u0vnl);
                    //cout<<u1vnl<<endl;
                    normals->InsertNextTuple3(u1vnl[0], u1vnl[1], u1vnl[2]);

                    u2vnl = vnl_cross_3d(u1vnl, u0vnl);
                    u2vnl.normalize();
                }else{
                    transform->Identity();
                    double angle = ((double)v)*360.0/((double)numspokes);
                    transform->RotateWXYZ(angle, u2vnl[0], u2vnl[1], u2vnl[2]);

                    vtkSmartPointer<vtkMatrix4x4> matrix = transform->GetMatrix();

                    vnl_matrix<double> mat(3, 3);
                    for(unsigned i = 0; i < 3; i++){
                        for(unsigned j = 0; j < 3; j++){
                            mat[i][j] = matrix->GetElement(i,j);
                        }
                    }

                    //cout<<mat<<endl;
                    //cout<<u0vnl<<endl;

                    s = (u0vnl*mat);
                    //cout<<s<<endl;
                    s *= prim0->getRN(v);

                }


                radius.push_back(s.magnitude());
                vnlspokes.push_back(s.normalize());*/


            }

            m_Srep->GetPointData()->SetNormals(normals);

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

        m_Srep->SetPoints(hubpos);
        m_Srep->SetPolys(cellarray);
        m_Srep->SetAllSpokes(allspokes);
        m_Srep->SetAllRadius(allradius);

        const float *color = tubefig->getColor();

        m_Srep->SetColor(color[0], color[1], color[2]);

        allspokes.clear();
        allradius.clear();

        delete tubefig;


    }
}

M3DFigure* vtkReadSRep::GetFigure(){

    const char* figfilename = GetFileName();

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
    if(m3dobject){
        M3DFigure * figure = m3dobject->getFigurePtr( 0 ) ;

       return figure;
    }

    return 0;
}
