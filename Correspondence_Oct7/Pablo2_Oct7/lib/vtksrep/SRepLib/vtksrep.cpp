#include "vtksrep.h"
#include "vtkPolyData.h"

#include "vtkSmartPointer.h"
#include "vtkIdList.h"
#include "vtkCell.h"
#include "vtkQuad.h"
#include "vtkLine.h"
#include "vtkGenericCell.h"
#include "vtkTriangle.h"

#include "vnl/vnl_cross.h"
#include "vtkPointData.h"

#include "vtkObjectFactory.h"
vtkStandardNewMacro(vtkSRep)



vtkSRep::vtkSRep()
{

    m_NumRows = -1;
    m_NumColumns = -1;
    m_UVToId = 0;
    m_UVToIdLocator = 0;

    m_UVToIdCrest = 0;
    m_UVToIdCrestLocator = 0;
    m_NumberOfSpokes = -1;
    m_Color[0] = 1;
    m_Color[1] = 0;
    m_Color[2] = 0;
}

vtkSRep::~vtkSRep()
{
    for(unsigned i = 0; i < m_Spokes.size(); i++){
        m_Spokes[i].clear();
        m_SpokesRadius[i].clear();
    }
    m_Spokes.clear();
    m_SpokesRadius.clear();    
}



/* \brief Given a vector of hub ids and the spoke position
          returns a list containing the hub spoke at the spokeposition
          (each hub contains 1..N spokes)
*/
vtkSRep::VectorVNLType vtkSRep::GetSpokes(VectorIdsType vecthubids, vtkIdType spokepos, bool unit){
    VectorVNLType spokesatpos;
    for(unsigned i = 0; i < vecthubids.size(); i++){
        VectorVNLType spokes_hub = GetSpokes(vecthubids[i], unit);
        if(spokepos < (int)spokes_hub.size()  ){            
            spokesatpos.push_back(spokes_hub[spokepos]);
        }else{
            VNLType zero(3);
            zero.fill(0);
            spokesatpos.push_back(zero);
        }
    }
    return spokesatpos;
}

vtkSRep::VectorVNLType vtkSRep::GetDerivatives(vtkIdType pointid){

    VectorVNLType derivatives;


    vtkSmartPointer< vtkIdList > cellsid = vtkSmartPointer< vtkIdList >::New();

    this->GetPointCells(pointid, cellsid);

    VNLType uderivative(3);
    uderivative.fill(0);
    VNLType vderivative(3);
    vderivative.fill(0);


    vtkIdType numids = cellsid->GetNumberOfIds();

    if(numids == 1){
        vtkSmartPointer<vtkGenericCell> quad0 = vtkSmartPointer<vtkGenericCell>::New();
        this->GetCell(cellsid->GetId(0), quad0);

        if(quad0->GetNumberOfPoints() == 4){//slab

            double pointnext[3] = {0, 0, 0};
            double pointprev[3] = {0, 0, 0};

            for(int i = 0; i < quad0->GetPointIds()->GetNumberOfIds(); i++){
                if(quad0->GetPointId(i) == pointid){

                    vtkIdType nextid = i + 1;
                    vtkIdType previd = i - 1;

                    if(nextid == 4){
                        nextid = 0;
                    }else if(previd < 0){
                        previd = 3;
                    }

                    this->GetPoint(quad0->GetPointIds()->GetId(nextid), pointnext);
                    this->GetPoint(quad0->GetPointIds()->GetId(previd), pointprev);

                    if(i == 0 || i == 1){

                        vderivative[0] = (pointnext[0] - pointprev[0])/2.0;
                        vderivative[1] = (pointnext[1] - pointprev[1])/2.0;
                        vderivative[2] = (pointnext[2] - pointprev[2])/2.0;

                        if(i == 0){
                            uderivative = -1.0*vderivative;
                        }else{
                            uderivative = vderivative;
                        }
                    }else{

                        uderivative[0] = (pointnext[0] - pointprev[0])/2.0;
                        uderivative[1] = (pointnext[1] - pointprev[1])/2.0;
                        uderivative[2] = (pointnext[2] - pointprev[2])/2.0;

                        if(i == 2){
                            vderivative = -1.0*uderivative;
                        }else{
                            vderivative = -1.0*uderivative;
                            uderivative = vderivative;
                        }

                    }

                }
            }

            derivatives.push_back(uderivative);
            derivatives.push_back(vderivative);
        }else{//tube
            double pointnext[3] = {0, 0, 0};
            double pointprev[3] = {0, 0, 0};

            this->GetPoint(quad0->GetPointIds()->GetId(0), pointprev);
            this->GetPoint(quad0->GetPointIds()->GetId(1), pointnext);

            uderivative[0] = pointnext[0] - pointprev[0];
            uderivative[1] = pointnext[1] - pointprev[1];
            uderivative[2] = pointnext[2] - pointprev[2];

            derivatives.push_back(uderivative);
            derivatives.push_back(vderivative);
        }

    }else{        


        vtkSmartPointer<vtkGenericCell> quad = vtkSmartPointer<vtkGenericCell>::New();
        vtkSmartPointer<vtkGenericCell> quad0 = vtkSmartPointer<vtkGenericCell>::New();
        vtkSmartPointer<vtkGenericCell> quad1 = vtkSmartPointer<vtkGenericCell>::New();
        bool vboolprev = false, vboolnext = false, uboolnext = false, uboolprev = false;

        double upointnext[3] = {0, 0, 0};
        double upointprev[3] = {0, 0, 0};

        double vpointprev[3] = {0, 0, 0};
        double vpointnext[3] = {0, 0, 0};

        double center[3] = {0, 0, 0};
        this->GetPoint(pointid, center);

        bool tube = false;

        if(numids == 2){
            this->GetCell(cellsid->GetId(0), quad);
            if(quad->GetNumberOfPoints() == 2){
                tube = true;
            }
        }

        if(tube){

            this->GetCell(cellsid->GetId(0), quad0);
            this->GetPoint(quad0->GetPointIds()->GetId(0), upointprev);

            this->GetCell(cellsid->GetId(1), quad1);
            this->GetPoint(quad1->GetPointIds()->GetId(1), upointnext);

            uboolnext = true;
            uboolprev = true;

        }else{


            vtkSmartPointer<vtkGenericCell> quad2 = vtkSmartPointer<vtkGenericCell>::New();
            vtkSmartPointer<vtkGenericCell> quad3 = vtkSmartPointer<vtkGenericCell>::New();

            //cout<<"center: "<<center[0]<<" "<<center[1]<<" "<<center[2]<<endl;

            for(unsigned i = 0; i < numids; i++){
                this->GetCell(cellsid->GetId(i), quad);

                if(quad->GetNumberOfPoints() == 4){

                    if(quad->GetPointIds()->GetId(0) == pointid){

                        this->GetCell(cellsid->GetId(i), quad0);

                        this->GetPoint(quad0->GetPointIds()->GetId(1), vpointnext);
                        this->GetPoint(quad0->GetPointIds()->GetId(3), upointnext);

                        vboolnext = true;
                        uboolnext = true;

                        //cout<<"right: "<<vpointnext[0]<<" "<<vpointnext[1]<<" "<<vpointnext[2]<<endl;
                        //cout<<"up: "<<upointnext[0]<<" "<<upointnext[1]<<" "<<upointnext[2]<<endl;

                    }else if(quad->GetPointIds()->GetId(1) == pointid){

                        this->GetCell(cellsid->GetId(i), quad1);

                        this->GetPoint(quad1->GetPointIds()->GetId(0), vpointprev);
                        vboolprev = true;

                        if(quad1->GetPointIds()->GetNumberOfIds() > 2){
                            this->GetPoint(quad1->GetPointIds()->GetId(2), upointnext);
                            uboolnext = true;
                        }

                        //cout<<"left: "<<vpointprev[0]<<" "<<vpointprev[1]<<" "<<vpointprev[2]<<endl;
                        //cout<<"up: "<<upointnext[0]<<" "<<upointnext[1]<<" "<<upointnext[2]<<endl;

                    }else if(quad->GetPointIds()->GetId(2) == pointid){

                        this->GetCell(cellsid->GetId(i), quad2);

                        this->GetPoint(quad2->GetPointIds()->GetId(3), vpointprev);
                        this->GetPoint(quad2->GetPointIds()->GetId(1), upointprev);

                        vboolprev = true;
                        uboolprev = true;

                        //cout<<"left: "<<vpointprev[0]<<" "<<vpointprev[1]<<" "<<vpointprev[2]<<endl;
                        //cout<<"down: "<<upointprev[0]<<" "<<upointprev[1]<<" "<<upointprev[2]<<endl;

                    }else if(quad->GetPointIds()->GetId(3) == pointid){
                        this->GetCell(cellsid->GetId(i), quad3);

                        this->GetPoint(quad3->GetPointIds()->GetId(0), upointprev);
                        this->GetPoint(quad3->GetPointIds()->GetId(2), vpointnext);

                        uboolprev = true;
                        vboolnext = true;

                        //cout<<"right: "<<vpointnext[0]<<" "<<vpointnext[1]<<" "<<vpointnext[2]<<endl;
                        //cout<<"down: "<<upointprev[0]<<" "<<upointprev[1]<<" "<<upointprev[2]<<endl;
                    }
                }
            }
        }

        if(uboolnext && uboolprev ){
            uderivative[0] = (upointnext[0] - upointprev[0])/2.0;
            uderivative[1] = (upointnext[1] - upointprev[1])/2.0;
            uderivative[2] = (upointnext[2] - upointprev[2])/2.0;
        }else if(uboolnext && !uboolprev ){
            uderivative[0] = (upointnext[0] - center[0]);
            uderivative[1] = (upointnext[1] - center[1]);
            uderivative[2] = (upointnext[2] - center[2]);
        }else if(!uboolnext && uboolprev ){
            uderivative[0] = (center[0] - upointprev[0]);
            uderivative[1] = (center[1] - upointprev[1]);
            uderivative[2] = (center[2] - upointprev[2]);
        }else{
            uderivative.fill(0);
        }

        if(vboolnext && vboolprev){
            vderivative[0] = (vpointnext[0] - vpointprev[0])/2.0;
            vderivative[1] = (vpointnext[1] - vpointprev[1])/2.0;
            vderivative[2] = (vpointnext[2] - vpointprev[2])/2.0;
        }else if(vboolnext && !vboolprev){
            vderivative[0] = (vpointnext[0] - center[0]);
            vderivative[1] = (vpointnext[1] - center[1]);
            vderivative[2] = (vpointnext[2] - center[2]);
        }else if(!vboolnext && vboolprev){
            vderivative[0] = (center[0] - vpointprev[0]);
            vderivative[1] = (center[1] - vpointprev[1]);
            vderivative[2] = (center[2] - vpointprev[2]);
        }else{
            vderivative.fill(0);
        }

        //cout<<"uder :" <<uderivative[0]<<" "<<uderivative[1]<<" "<<uderivative[2]<<endl;
        //cout<<"vder :" <<vderivative[0]<<" "<<vderivative[1]<<" "<<vderivative[2]<<endl;

        derivatives.push_back(uderivative);
        derivatives.push_back(vderivative);
    }


    return derivatives;
}


/*vtkSRep::VectorVNLType vtkSRep::GetDerivatives(vtkIdType rowid, vtkIdType colid){
    return GetDerivatives(m_GridTopolgyIds[rowid][colid]);
}*/

vtkSRep::VectorVNLType vtkSRep::GetDerivativesSpoke(vtkIdType pointid, vtkIdType spokepos, bool unit){

        VectorVNLType derivatives;


        vtkSmartPointer< vtkIdList > cellsid = vtkSmartPointer< vtkIdList >::New();

        this->GetPointCells(pointid, cellsid);

        VNLType uderivative;
        VNLType vderivative;


        vtkIdType numids = cellsid->GetNumberOfIds();

        if(numids == 1){
            vtkSmartPointer<vtkGenericCell> quad0 = vtkSmartPointer<vtkGenericCell>::New();
            this->GetCell(cellsid->GetId(0), quad0);

            VNLType pointnext;
            VNLType pointprev;

            for(int i = 0; i < quad0->GetPointIds()->GetNumberOfIds(); i++){
                if(quad0->GetPointId(i) == pointid){

                    vtkIdType nextid = i + 1;
                    vtkIdType previd = i - 1;

                    if(nextid == 4){
                        nextid = 0;
                    }else if(previd < 0){
                        previd = 3;
                    }

                    pointnext = this->GetSpoke(quad0->GetPointIds()->GetId(nextid), spokepos);
                    pointprev = this->GetSpoke(quad0->GetPointIds()->GetId(previd),spokepos);

                    if(i == 0 || i == 1){

                        vderivative = (pointnext - pointprev)/2.0;

                        if(i == 0){
                            uderivative = -1.0*vderivative;
                        }else{
                            uderivative = vderivative;
                        }
                    }else{

                        uderivative = (pointnext - pointprev)/2.0;

                        if(i == 2){
                            vderivative = -1.0*uderivative;
                        }else{
                            vderivative = -1.0*uderivative;
                            uderivative = vderivative;
                        }
                    }
                }
            }

            derivatives.push_back(uderivative);
            derivatives.push_back(vderivative);

        }else{

            vtkSmartPointer<vtkGenericCell> quad = vtkSmartPointer<vtkGenericCell>::New();
            vtkSmartPointer<vtkGenericCell> quad0 = vtkSmartPointer<vtkGenericCell>::New();
            vtkSmartPointer<vtkGenericCell> quad1 = vtkSmartPointer<vtkGenericCell>::New();
            vtkSmartPointer<vtkGenericCell> quad2 = vtkSmartPointer<vtkGenericCell>::New();
            vtkSmartPointer<vtkGenericCell> quad3 = vtkSmartPointer<vtkGenericCell>::New();

            VNLType center = this->GetSpoke(pointid, spokepos, unit);

            VNLType upointnext;
            VNLType upointprev;

            VNLType vpointprev;
            VNLType vpointnext;

            bool vboolprev = false, vboolnext = false, uboolnext = false, uboolprev = false;

            for(unsigned i = 0; i < numids; i++){

                this->GetCell(cellsid->GetId(i), quad);

                if(quad->GetPointIds()->GetId(0) == pointid){
                    this->GetCell(cellsid->GetId(i), quad0);

                    upointnext = this->GetSpoke(quad0->GetPointIds()->GetId(3), spokepos);
                    vpointnext = this->GetSpoke(quad0->GetPointIds()->GetId(1),spokepos);

                    if(unit){
                        upointnext.normalize();
                        vpointnext.normalize();
                    }

                    uboolnext = true;
                    vboolnext = true;

                    //cout<<"right: "<<vpointnext[0]<<" "<<vpointnext[1]<<" "<<vpointnext[2]<<end
                    //cout<<"up: "<<upointnext[0]<<" "<<upointnext[1]<<" "<<upointnext[2]<<endl;

                }else if(quad->GetPointIds()->GetId(1) == pointid){
                    this->GetCell(cellsid->GetId(i), quad1);

                    upointnext = this->GetSpoke(quad1->GetPointIds()->GetId(2),spokepos);
                    vpointprev = this->GetSpoke(quad1->GetPointIds()->GetId(0),spokepos);

                    if(unit){
                        upointnext.normalize();
                        vpointprev.normalize();
                    }

                    uboolnext = true;
                    vboolprev = true;

                    //cout<<"left: "<<vpointprev[0]<<" "<<vpointprev[1]<<" "<<vpointprev[2]<<endl;
                    //cout<<"up: "<<upointnext[0]<<" "<<upointnext[1]<<" "<<upointnext[2]<<endl;

                }else if(quad->GetPointIds()->GetId(2) == pointid){
                    this->GetCell(cellsid->GetId(i), quad2);

                    upointprev = this->GetSpoke(quad2->GetPointIds()->GetId(1), spokepos);
                    vpointprev = this->GetSpoke(quad2->GetPointIds()->GetId(3), spokepos);

                    if(unit){
                        upointprev.normalize();
                        vpointprev.normalize();
                    }

                    uboolprev = true;
                    vboolprev = true;

                    //cout<<"left: "<<vpointprev[0]<<" "<<vpointprev[1]<<" "<<vpointprev[2]<<endl;
                    //cout<<"down: "<<upointprev[0]<<" "<<upointprev[1]<<" "<<upointprev[2]<<endl;

                }else if(quad->GetPointIds()->GetId(3) == pointid){
                    this->GetCell(cellsid->GetId(i), quad3);

                    upointprev = this->GetSpoke(quad3->GetPointIds()->GetId(0),spokepos);
                    vpointnext = this->GetSpoke(quad3->GetPointIds()->GetId(2),spokepos);

                    if(unit){
                        upointprev.normalize();
                        vpointnext.normalize();
                    }

                    uboolprev = true;
                    vboolnext = true;

                    //cout<<"right: "<<vpointnext[0]<<" "<<vpointnext[1]<<" "<<vpointnext[2]<<endl;
                    //cout<<"down: "<<upointprev[0]<<" "<<upointprev[1]<<" "<<upointprev[2]<<endl;
                }
            }

            if(uboolnext && uboolprev ){
                uderivative = (upointnext - upointprev)/2.0;
            }else if(uboolnext && !uboolprev ){
                uderivative = (upointnext - center);                
                VNLType N = this->GetMedialSheetNormal(pointid);
                uderivative = uderivative - dot_product(uderivative, N) * N;
            }else if(!uboolnext && uboolprev ){
                uderivative = (center - upointprev);
                VNLType N = this->GetMedialSheetNormal(pointid);
                uderivative = uderivative - dot_product(uderivative, N) * N;
            }else{
                uderivative.fill(0);
            }

            if(vboolnext && vboolprev){
                vderivative = (vpointnext - vpointprev)/2.0;
            }else if(vboolnext && !vboolprev){
                vderivative = (vpointnext - center);
                VNLType N = this->GetMedialSheetNormal(pointid);
                vderivative = vderivative - dot_product(vderivative, N) * N;
            }else if(!vboolnext && vboolprev){
                vderivative = (center - vpointprev);
                VNLType N = this->GetMedialSheetNormal(pointid);
                vderivative = vderivative - dot_product(vderivative, N) * N;
            }else{
                vderivative.fill(0);
            }

            //cout<<"center :" <<center<<endl;
            //cout<<"un - up = "<<upointnext<<" - "<<upointprev<<endl;
            //cout<<"vn - vp = "<<vpointnext<<" - "<<vpointprev<<endl;
            //cout<<"uder :" <<uderivative<<endl;
            //cout<<"vder :" <<vderivative<<endl;

            derivatives.push_back(uderivative);
            derivatives.push_back(vderivative);

        }

        return derivatives;
}

/*vtkSRep::VectorVNLType vtkSRep::GetDerivativesSpokeRadius(vtkIdType pointid, vtkIdType spokepos){
    VectorVNLType derivatives;


    vtkSmartPointer< vtkIdList > cellsid = vtkSmartPointer< vtkIdList >::New();

    this->GetPointCells(pointid, cellsid);

    double uderivative = 0;
    double vderivative = 0;


    vtkIdType numids = cellsid->GetNumberOfIds();

    if(numids == 1){
        vtkSmartPointer<vtkGenericCell> quad0 = vtkSmartPointer<vtkGenericCell>::New();
        this->GetCell(cellsid->GetId(0), quad0);


        double pointnext = 0;
        double pointprev = 0;

        for(int i = 0; i < quad0->GetPointIds()->GetNumberOfIds(); i++){
            if(quad0->GetPointId(i) == pointid){

                vtkIdType nextid = i + 1;
                vtkIdType previd = i - 1;

                if(nextid == 4){
                    nextid = 0;
                }else if(previd < 0){
                    previd = 3;
                }

                pointnext = this->GetSpokeRadius(quad0->GetPointIds()->GetId(nextid), spokepos);
                pointprev = this->GetSpokeRadius(quad0->GetPointIds()->GetId(previd),spokepos);

                if(i == 0 || i == 1){

                    vderivative = (pointnext - pointprev)/2.0;

                    if(i == 0){
                        uderivative = -1.0*vderivative;
                    }else{
                        uderivative = vderivative;
                    }
                }else{

                    uderivative = (pointnext - pointprev)/2.0;

                    if(i == 2){
                        vderivative = -1.0*uderivative;
                    }else{
                        vderivative = -1.0*uderivative;
                        uderivative = vderivative;
                    }
                }
            }
        }

        VNLType uder(1);
        uder[0] = uderivative;

        VNLType vder(1);
        vder[0] = vderivative;

        derivatives.push_back(uder);
        derivatives.push_back(vder);

    }else{

        vtkSmartPointer<vtkGenericCell> quad = vtkSmartPointer<vtkGenericCell>::New();
        vtkSmartPointer<vtkGenericCell> quad0 = vtkSmartPointer<vtkGenericCell>::New();
        vtkSmartPointer<vtkGenericCell> quad1 = vtkSmartPointer<vtkGenericCell>::New();
        vtkSmartPointer<vtkGenericCell> quad2 = vtkSmartPointer<vtkGenericCell>::New();
        vtkSmartPointer<vtkGenericCell> quad3 = vtkSmartPointer<vtkGenericCell>::New();

        double center = this->GetSpokeRadius(pointid, spokepos);

        double upointnext = 0;
        double upointprev = 0;

        double vpointprev = 0;
        double vpointnext = 0;

        bool vboolprev = false, vboolnext = false, uboolnext = false, uboolprev = false;

        for(unsigned i = 0; i < numids; i++){
            this->GetCell(cellsid->GetId(i), quad);
            if(quad->GetPointIds()->GetId(0) == pointid){
                this->GetCell(cellsid->GetId(i), quad0);

                upointnext = this->GetSpokeRadius(quad0->GetPointIds()->GetId(3), spokepos);
                vpointnext = this->GetSpokeRadius(quad0->GetPointIds()->GetId(1),spokepos);

                uboolnext = true;
                vboolnext = true;

                //cout<<"right: "<<vpointnext[0]<<" "<<vpointnext[1]<<" "<<vpointnext[2]<<end
                //cout<<"up: "<<upointnext[0]<<" "<<upointnext[1]<<" "<<upointnext[2]<<endl;

            }else if(quad->GetPointIds()->GetId(1) == pointid){
                this->GetCell(cellsid->GetId(i), quad1);

                upointnext = this->GetSpokeRadius(quad1->GetPointIds()->GetId(2),spokepos);
                vpointprev = this->GetSpokeRadius(quad1->GetPointIds()->GetId(0),spokepos);

                uboolnext = true;
                vboolprev = true;

                //cout<<"left: "<<vpointprev[0]<<" "<<vpointprev[1]<<" "<<vpointprev[2]<<endl;
                //cout<<"up: "<<upointnext[0]<<" "<<upointnext[1]<<" "<<upointnext[2]<<endl;

            }else if(quad->GetPointIds()->GetId(2) == pointid){
                this->GetCell(cellsid->GetId(i), quad2);

                upointprev = this->GetSpokeRadius(quad2->GetPointIds()->GetId(1), spokepos);
                vpointprev = this->GetSpokeRadius(quad2->GetPointIds()->GetId(3), spokepos);

                uboolprev = true;
                vboolprev = true;

                //cout<<"left: "<<vpointprev[0]<<" "<<vpointprev[1]<<" "<<vpointprev[2]<<endl;
                //cout<<"down: "<<upointprev[0]<<" "<<upointprev[1]<<" "<<upointprev[2]<<endl;

            }else if(quad->GetPointIds()->GetId(3) == pointid){
                this->GetCell(cellsid->GetId(i), quad3);

                upointprev = this->GetSpokeRadius(quad3->GetPointIds()->GetId(0),spokepos);
                vpointnext = this->GetSpokeRadius(quad3->GetPointIds()->GetId(2),spokepos);

                uboolprev = true;
                vboolnext = true;

                //cout<<"right: "<<vpointnext[0]<<" "<<vpointnext[1]<<" "<<vpointnext[2]<<endl;
                //cout<<"down: "<<upointprev[0]<<" "<<upointprev[1]<<" "<<upointprev[2]<<endl;
            }
        }

        if(uboolnext && uboolprev ){
            uderivative = (upointnext - upointprev)/2.0;
        }else if(uboolnext && !uboolprev ){
            uderivative = (upointnext - center);
        }else if(!uboolnext && uboolprev ){
            uderivative = (center - upointprev);
        }

        if(vboolnext && vboolprev){
            vderivative = (vpointnext - vpointprev)/2.0;
        }else if(vboolnext && !vboolprev){
            vderivative = (vpointnext - center);
        }else if(!vboolnext && vboolprev){
            vderivative = (center - vpointprev);
        }

        VNLType uder(1);
        uder[0] = uderivative;

        VNLType vder(1);
        vder[0] = vderivative;

        derivatives.push_back(uder);
        derivatives.push_back(vder);

    }

    return derivatives;
}*/

vtkSRep::VectorVNLType vtkSRep::GetCrestMedialAtoms(){

    VectorIdsType medialcrestidsvector = this->GetCrestMedialAtomsIds();
    return GetCrestMedialAtoms(medialcrestidsvector);

}

// \brief given an id of a hub returns the e3 vector i.e. cross(topSpoke, bottomSpoke)
vtkSRep::VNLType vtkSRep::GetMedialSheetVecte3(vtkIdType hubid){
    VNLType e3 = vnl_cross_3d(GetTopSpoke(hubid), GetBottomSpoke(hubid));
    return e3.normalize();
}

// \brief given an id of a hub returns the normal to the sheet
vtkSRep::VNLType vtkSRep::GetMedialSheetNormal(vtkIdType hubid){
    VNLType normal(3);
    if(this->GetPointData()->GetNormals()){
        double tempnormal[3];
        this->GetPointData()->GetNormals()->GetTuple(hubid, tempnormal);
        normal[0] = tempnormal[0];
        normal[1] = tempnormal[1];
        normal[2] = tempnormal[2];
    }else{
        normal = GetTopSpoke(hubid) - GetBottomSpoke(hubid);
    }
    return normal.normalize();
}

// \brief given an id of a hub sets the normal to the sheet
void vtkSRep::SetMedialSheetNormal(vtkIdType hubid, VNLType normal){
    if(this->GetPointData()->GetNormals()){
        this->GetPointData()->GetNormals()->SetTuple3(hubid, normal[0], normal[1], normal[2]);
    }
}

// \brief given an id of a hub returns the normal to the sheet
vtkSRep::VectorVNLType vtkSRep::GetMedialSheetNormals(VectorIdsType hubids){
    VectorVNLType normals;
    for(unsigned i = 0; i < hubids.size(); i++){
        normals.push_back(GetMedialSheetNormal(hubids[i]));
    }
    return normals;
}

// \brief returns a list of vnl_vectors that contain the position of every crest medial atom
// \param crestids the list of the ids of crestmedialAtoms
vtkSRep::VectorVNLType vtkSRep::GetCrestMedialAtoms(VectorIdsType crestids){
    vtkSRep::VectorVNLType medialcrestvector;

    for(unsigned i=0; i < crestids.size();i++){

        double temp[3];
        this->GetPoint(crestids[i], temp);

        vtkSRep::VNLType vect(3);
        vect[0] = temp[0];
        vect[1] = temp[1];
        vect[2] = temp[2];

        medialcrestvector.push_back(vect);
    }
    return medialcrestvector;
}

/* \brief Calculates the next atom id that corresponds to a crest medial atom
   \params vtkIdType id  the current atom id that is to evaluate
           this id is use to recover all the cells that this point belong to
   \params vtkIdType inquadid corresponds to the current inquadid for the given id
           this is use to locate the correct cell
           if no cell is found to have the inquadid for the given id then
           inquadid is increased. This happens at the corners of the slab
    \params vectids returns all the ids of the points found
*/

vtkSRep::VectorIdsType vtkSRep::GetCrestMedialAtomsIds(){

    VectorIdsType medialcrestidsvector;

    if(this->GetNumberOfCells() > 0){

        vtkSmartPointer<vtkGenericCell> cell = vtkSmartPointer<vtkGenericCell>::New();
        this->GetCell(0, cell);

        vtkIdType id = cell->GetPointIds()->GetId(0);

        vtkIdType id1 = cell->GetPointIds()->GetId(1);

        medialcrestidsvector.push_back(id);

        int numids = cell->GetPointIds()->GetNumberOfIds();

        GetCrestMedialAtomsIds(id1, 0, medialcrestidsvector, numids);


    }

    return medialcrestidsvector;

}
void vtkSRep::GetCrestMedialAtomsIds(vtkIdType pointid, vtkIdType inquadid, VectorIdsType& vectids, int numids ){

    vtkSmartPointer< vtkIdList > cellsid = vtkSmartPointer< vtkIdList >::New();
    this->GetPointCells(pointid, cellsid);

    bool found = false;
    vtkIdType nextcellsid = 0;

    int numpoints = 0;

    if(cellsid->GetNumberOfIds() > 0){
        vtkSmartPointer<vtkGenericCell> cell = vtkSmartPointer<vtkGenericCell>::New();
        this->GetCell(cellsid->GetId(0), cell);
        numpoints = cell->GetNumberOfPoints();
    }

    for(unsigned i = 0; i < cellsid->GetNumberOfIds() && !found; i++){
        vtkSmartPointer<vtkGenericCell> cell = vtkSmartPointer<vtkGenericCell>::New();
        this->GetCell(cellsid->GetId(i), cell);

        if(cell->GetPointIds()->GetId(inquadid) == pointid){
            nextcellsid = cellsid->GetId(i);
            found = true;
        }
    }

    if(found){

        vtkSmartPointer<vtkGenericCell> cell = vtkSmartPointer<vtkGenericCell>::New();
        this->GetCell(nextcellsid, cell);

        vtkIdType nextid = inquadid + 1;

        if(inquadid + 1 == numpoints && numids > 2 ){
            nextid = 0;
        }

        vtkIdType id = cell->GetPointIds()->GetId(nextid);

        if(numids == 2){
            vectids.push_back(pointid);
            GetCrestMedialAtomsIds(id, inquadid, vectids, numids);
        }else{

            GetCrestMedialAtomsIds(id, inquadid, vectids, numids);
            vectids.push_back(pointid);
        }

    }else if(!found){//this means that is a corner
        if(inquadid + 1 < numpoints && numids > 2){
            GetCrestMedialAtomsIds(pointid, inquadid + 1, vectids, numids);
        }else if(numids == 2){
            vectids.push_back(pointid);
        }
    }
}

// \brief returns a list of vnl_vectors that contain the ids of every crest medial atom on the given slab side
//        side [0-3] as in the positions of a square.
vtkSRep::VectorIdsType vtkSRep::GetCrestMedialAtomsIds(int side){
    VectorIdsType allids = GetCrestMedialAtomsIds();
    if(allids.size() > 0)
        allids.push_back(allids[0]);

    VectorIdsType sideids;


    int currentside = 0;
    bool insert = false;
    bool last = false;
    for(unsigned i = 0; i < allids.size(); i++){
        vtkSmartPointer< vtkIdList > cellsid = vtkSmartPointer< vtkIdList >::New();
        this->GetPointCells(allids[i], cellsid);

        if(cellsid->GetNumberOfIds() == 1){
            if(currentside == side){
                insert = true;
            }else{
                if(insert){
                    last = true;
                }
                insert = false;
            }
            currentside++;
        }
        if(insert || last){
            if(last){
                last = false;
            }
            sideids.push_back(allids[i]);
        }
    }

    return sideids;
}

vtkSRep::VectorVNLType vtkSRep::GetCrestMedialAtomsDerivatives(bool unit){
    return GetCrestMedialAtomsDerivatives(GetCrestMedialAtoms(), unit);
}

vtkSRep::VectorVNLType vtkSRep::GetCrestMedialAtomsDerivatives(VectorIdsType crestids, bool unit){
    return GetCrestMedialAtomsDerivatives(GetCrestMedialAtoms(crestids), unit);
}

vtkSRep::VectorVNLType vtkSRep::GetCrestMedialAtomsDerivatives(VectorVNLType crestpositions, bool unit, bool cyclic){

    VectorVNLType crestderivatives;
    for(unsigned i=0; i < crestpositions.size(); i++){
        VNLType dp0(3);

        if(i == 0){
            if(cyclic){
                dp0 = (crestpositions[i+1] - crestpositions[crestpositions.size()-1])/2.0;
            }else{
                dp0 = (crestpositions[i+1] - crestpositions[0]);
            }
        }else if(i == crestpositions.size() - 1){
            if(cyclic){
                dp0 = (crestpositions[0] - crestpositions[i-1])/2.0;
            }else{
                dp0 = (crestpositions[i] - crestpositions[i-1]);
            }
        }else{
            dp0 = (crestpositions[i+1] - crestpositions[i-1])/2.0;
        }

        if(unit){
            dp0 = dp0.normalize();
        }

        crestderivatives.push_back(dp0);
    }

    return crestderivatives;
}

vtkSRep::VectorVNLType vtkSRep::GetCrestUDerivatives(bool unit, vtkIdType spokepos){
    return GetCrestUDerivatives(GetCrestMedialAtomsIds(), unit, spokepos);
}

vtkSRep::VectorVNLType vtkSRep::GetCrestUDerivatives(VectorIdsType crestids, bool unit, vtkIdType spokepos){
    VectorVNLType crestvectors;
    for(unsigned i=0; i < crestids.size(); i++){
        vtkIdType id = crestids[i];
        VectorVNLType spokes = this->GetSpokes(id);
        VNLType s = spokes[spokepos];

        VectorDoubleType radius = this->GetSpokesRadius(crestids[i]);
        double rad = radius[spokepos];

        VNLType crestspoke = s * rad;
        crestvectors.push_back(crestspoke);
    }
    return GetCrestUDerivatives(crestvectors, unit);
}

vtkSRep::VectorVNLType vtkSRep::GetCrestUDerivatives(VectorVNLType crestvectors, bool unit){
    return GetCrestMedialAtomsDerivatives(crestvectors, unit);
}

vtkSRep::VNLMatrixType vtkSRep::GetSRadMatrix(vtkIdType hubid, vtkIdType spokepos){

    VNLType Uvect = this->GetSpoke(hubid, spokepos, true);//unit vector
    VectorVNLType derivatives = this->GetDerivatives(hubid);//derivatives
    VectorVNLType Uderivatives = this->GetDerivativesSpoke(hubid, spokepos, false);//unit spoke derivatives
    double r = this->GetSpokeRadius(hubid, spokepos);


    VNLMatrixType U(1,3);
    U[0][0] = Uvect[0];
    U[0][1] = Uvect[1];
    U[0][2] = Uvect[2];

    VNLMatrixType pUn(2, 3);
    pUn[0][0] = derivatives[0][0];
    pUn[0][1] = derivatives[0][1];
    pUn[0][2] = derivatives[0][2];
    pUn[1][0] = derivatives[1][0];
    pUn[1][1] = derivatives[1][1];
    pUn[1][2] = derivatives[1][2];

    VNLMatrixType I(3, 3);
    I.fill(0);
    I.fill_diagonal(1);

    VNLMatrixType UtU = U.transpose() * U;
    VNLMatrixType Q = pUn * (UtU - I);

    VNLMatrixType dSdu(2, 3);
    dSdu[0][0] = Uderivatives[0][0];
    dSdu[0][1] = Uderivatives[0][1];
    dSdu[0][2] = Uderivatives[0][2];
    dSdu[1][0] = Uderivatives[1][0];
    dSdu[1][1] = Uderivatives[1][1];
    dSdu[1][2] = Uderivatives[1][2];

    dSdu *= r;


    VNLMatrixType QQtinv = vnl_inverse(Q*Q.transpose());

    VNLMatrixType rSrad = (dSdu * (Q.transpose()*QQtinv));//equ 2.13 on QuiogHan See def of dSdu
    return rSrad.transpose();


    //Jared's
    /*VNLMatrixType C11(3, 3);

    C11[0][0] = derivatives[0][0];
    C11[1][0] = derivatives[0][1];
    C11[2][0] = derivatives[0][2];

    C11[0][1] = derivatives[1][0];
    C11[1][1] = derivatives[1][1];
    C11[2][1] = derivatives[1][2];

    C11[0][2] = Uvect[0];
    C11[1][2] = Uvect[1];
    C11[2][2] = Uvect[2];


    VNLMatrixType Ci11 = vnl_inverse(C11);

    VNLMatrixType Au11(3,1);
    Au11[0][0] = Uderivatives[0][0];
    Au11[1][0] = Uderivatives[0][1];
    Au11[2][0] = Uderivatives[0][2];

    VNLMatrixType Av11(3,1);
    Av11[0][0] = Uderivatives[1][0];
    Av11[1][0] = Uderivatives[1][1];
    Av11[2][0] = Uderivatives[1][2];

    VNLMatrixType Bu11 = Ci11 * Au11;
    VNLMatrixType Bv11 = Ci11 * Av11;


    VNLMatrixType Srad11(2,2);
    Srad11[0][0] = -1*Bu11[0][0];
    Srad11[1][0] = -1*Bu11[1][0];
    Srad11[0][1] = -1*Bv11[0][0];
    Srad11[1][1] = -1*Bv11[1][0];

    VNLMatrixType rSrad11 = Srad11 * r;

    return rSrad11;*/


}

vtkSRep::VectorIdsType vtkSRep::GetInternalMedialAtomIds(vtkIdType startid, vtkIdType nextincellid){

    vtkSmartPointer<vtkIdList> cellsidlist = vtkSmartPointer<vtkIdList>::New();
    this->GetPointCells(startid, cellsidlist);
    VectorIdsType vectids;

    vectids.push_back(startid);

    for(unsigned i = 0; i < cellsidlist->GetNumberOfIds(); i++){
        vtkSmartPointer< vtkGenericCell > cell = vtkSmartPointer< vtkGenericCell >::New();

        this->GetCell(cellsidlist->GetId(i), cell);
        if(cell->GetPointId(0) == startid && nextincellid < cell->GetNumberOfPoints()){
            GetInternalMedialAtomIds(cell->GetPointId(nextincellid), nextincellid, vectids);
        }
    }


    return vectids;
}

void vtkSRep::GetInternalMedialAtomIds(vtkIdType startid, vtkIdType nextincellid, VectorIdsType& vectids ){
    vtkSmartPointer<vtkIdList> cellsidlist = vtkSmartPointer<vtkIdList>::New();
    this->GetPointCells(startid, cellsidlist);

    vectids.push_back(startid);

    for(unsigned i = 0; i < cellsidlist->GetNumberOfIds(); i++){
        vtkSmartPointer< vtkGenericCell > cell = vtkSmartPointer< vtkGenericCell >::New();
        this->GetCell(cellsidlist->GetId(i), cell);
        if(cell->GetPointId(0) == startid && nextincellid < cell->GetNumberOfPoints()){
            GetInternalMedialAtomIds(cell->GetPointId(nextincellid), nextincellid, vectids);
        }
    }

}

double vtkSRep::GetCellArea(vtkIdType cellid){
    vtkSmartPointer< vtkGenericCell > cell = vtkSmartPointer< vtkGenericCell >::New();
    this->GetCell(cellid, cell);

    double p0[3], p1[3], p2[3], p3[3];

    if(dynamic_cast<vtkQuad*>(cell.GetPointer()) != 0){
        cell->GetPoints()->GetPoint(0, p0);
        cell->GetPoints()->GetPoint(1, p1);
        cell->GetPoints()->GetPoint(2, p2);
        cell->GetPoints()->GetPoint(3, p3);

        vtkSmartPointer<vtkTriangle> cell0 = vtkSmartPointer<vtkTriangle>::New();
        cell0->GetPoints()->SetPoint(0, p0[0], p0[1], p0[2]);
        cell0->GetPoints()->SetPoint(1, p1[0], p1[1], p1[2]);
        cell0->GetPoints()->SetPoint(2, p3[0], p3[1], p3[2]);

        vtkSmartPointer<vtkTriangle> cell1 = vtkSmartPointer<vtkTriangle>::New();
        cell1->GetPoints()->SetPoint(0, p1[0], p1[1], p1[2]);
        cell1->GetPoints()->SetPoint(1, p2[0], p2[1], p2[2]);
        cell1->GetPoints()->SetPoint(2, p3[0], p3[1], p3[2]);


        return cell0->ComputeArea() + cell1->ComputeArea();
    }else if(dynamic_cast<vtkLine*>(cell.GetPointer()) != 0){
        cout<<"Cell type is a line. Returning lenght of cell!"<<endl;
        cell->GetPoints()->GetPoint(0, p0);
        cell->GetPoints()->GetPoint(1, p1);
        VNLType vnlp0(p0,3);
        VNLType vnlp1(p1,3);
        return (vnlp0 - vnlp1).magnitude();
    }else{
        cout<<"Cell type not recognized. Area not calculated!"<<endl;
        return -1;
    }
}

double vtkSRep::GetCellAreaAverage(){
    double avg = 0;
    for(unsigned i = 0; i < this->GetNumberOfCells(); i++){
        avg+= GetCellArea(i);
    }
    return avg/this->GetNumberOfCells();
}

int vtkSRep::GetNumRows(){
    if(m_NumRows == -1){
        int nc = this->GetNumColumns() - 1;
        if(nc != 0){
            int numr = (int)(this->GetNumberOfCells()/nc);
            if(numr == 1 && this->GetCell(0)->GetNumberOfPoints() == 2){
                //This is a tube figure
                m_NumRows = 0; //only one row
            }else{
                m_NumRows = numr; //this is forslab
            }

            m_NumRows += 1;

        }else{
            m_NumRows = 1.0;
        }
    }
    return m_NumRows;
}

int vtkSRep::GetNumColumns(){
    if(m_NumColumns == -1){        
        m_NumColumns = this->GetCrestMedialAtomsIds(0).size();
    }

    return m_NumColumns;
}


void vtkSRep::DeepCopy(vtkDataObject *src){

    vtkPolyData::DeepCopy(src);

    if(dynamic_cast<vtkSRep*>(src)){

        vtkSRep* srepfig = dynamic_cast<vtkSRep*>(src);
        double *color = srepfig->GetColor();
        this->SetColor(color[0], color[1], color[2]);

        this->SetAllSpokes(srepfig->GetAllSpokes());
        this->SetAllRadius(srepfig->GetAllRadius());
    }
}
