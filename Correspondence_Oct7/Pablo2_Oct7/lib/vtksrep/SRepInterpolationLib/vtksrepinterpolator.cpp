#include "vtksrepinterpolator.h"


#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkQuad.h"
#include "vtkTriangle.h"

#include "vtkCellArray.h"
#include "vtkPoints.h"
#include "vtkPointData.h"
#include "vtkCleanPolyData.h"

#include "vtksrepinterpolatemedialsheet.h"

#include "vnl/algo/vnl_symmetric_eigensystem.h"
#include "vnl/algo/vnl_determinant.h"

#include "vtkExecutive.h"

#include "vtkObjectFactory.h"
vtkStandardNewMacro(vtkSRepInterpolator);

vtkSRepInterpolator::vtkSRepInterpolator()
{
    m_Input = 0;
    m_Medialsheetinterpolator = 0;
    m_Interpolatecrestspokes = 0;
    m_InterpolatecrestspokesBottom = 0;
    m_InterpolatedLabel = 0;

}


vtkSRepInterpolator::~vtkSRepInterpolator(){

    m_CellLocators.clear();


}

/*void vtkSRepInterpolator::SetInput(int index, vtkDataObject* input){
    if(input)
      {
      this->SetInputConnection(index, input->GetProducerPort());
      }
    else
      {
      // Setting a NULL input removes the connection.
      this->SetInputConnection(index, 0);
      }
}
void vtkSRepInterpolator::SetInput(vtkDataObject* input){
    this->SetInput(0, input);
}*/

void vtkSRepInterpolator::Update(){

    if(!m_Input){
        return;
    }

    GenerateUVMap();

    m_Medialsheetinterpolator = vtkSmartPointer< vtkSRepInterpolateMedialSheet >::New();
    m_Medialsheetinterpolator->SetInput(m_Input);
    m_Medialsheetinterpolator->Update();

    m_Medialspokesinterpolator = vtkSmartPointer< vtkSRepInterpolateMedialSpokesHermite >::New();
    m_Medialspokesinterpolator->SetInput(m_Input);
    m_Medialspokesinterpolator->Update();

    m_Interpolatecrestspokes = vtkSmartPointer<vtkSRepInterpolateCrestSpokesQuartic>::New();    
    m_Interpolatecrestspokes->SetCyclicCurve(true);
    m_Interpolatecrestspokes->SetInput(m_Input);
    m_Interpolatecrestspokes->Update();

}

void vtkSRepInterpolator::GenerateUVMap(){

    int numrows = m_Input->GetNumRows();
    int numcols = m_Input->GetNumColumns();

    m_CellLocators.clear();

    for(unsigned side = 0; side < 2; side++){
        vtkSmartPointer<vtkPoints> uvpoints = vtkSmartPointer<vtkPoints>::New();

        for(unsigned i = 0; i < m_Input->GetNumberOfPoints(); i++){
            int r = i%numcols;
            int g = i/numcols;

            int cr = r+1;
            int cg = 0;

            if(side == 0){
                cg = g+1;
            }else{
                cg = g+1 + (numrows-g)*2;
            }

            uvpoints->InsertNextPoint(cr, cg, 0);
            //cout<<cr<<" "<<cg<<endl;
        }

        vtkSmartPointer<vtkPolyData> poly = vtkSmartPointer<vtkPolyData>::New();

        poly->SetPoints(uvpoints);
        vtkSmartPointer<vtkCellArray> polys = vtkSmartPointer<vtkCellArray>::New();
        polys->DeepCopy(m_Input->GetPolys());
        poly->SetPolys(polys);
        poly->BuildCells();
        poly->BuildLinks();
        vtkSmartPointer<vtkCellLocator> uvtoidlocator = vtkSmartPointer<vtkCellLocator>::New();
        uvtoidlocator->SetDataSet(poly);
        uvtoidlocator->SetTolerance(0);
        uvtoidlocator->AutomaticOn();
        //m_UVToIdLocator->SetNumberOfPointsPerBucket(3);
        uvtoidlocator->BuildLocator();
        m_CellLocators.push_back(uvtoidlocator);
    }

    vector< vector< vtkIdType > > pointsids;

    vtkSmartPointer<vtkPoints> uvpoints = vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkCellArray> cellarray = vtkSmartPointer<vtkCellArray>::New();

    for(int y = 0; y < numrows + 2; y++){
        pointsids.push_back(vector< vtkIdType >());
        for(int x = 0; x < numcols + 2; x++){
            if(x <= 1 || y <= 1 || x >= numcols || y >= numrows){
                pointsids[y].push_back(uvpoints->InsertNextPoint(x, y, 0));
                //cout<<x<<" "<<y<<endl;
            }else{
                pointsids[y].push_back(-1);
            }

        }
        //cout<<endl;
    }


    //crest interpolation
    for(unsigned j = 0; j < 4; j++){

        int startx = 0;
        int endx = 0;

        int starty = 0;
        int endy = 0;

        int dx = 1;
        int dy = 1;

        //bool start = true;

        if(j == 0){
            startx = 1;
            endx = numcols;
            starty = 1;
            endy = 0;
            dy = -1;
        }else if(j == 1){
            startx = numcols;
            endx = numcols+1;
            dx=1;
            starty = 1;
            endy = numrows;
        }else if(j == 2){
            startx = numcols;
            endx = 1;
            dx = -1;
            starty=numrows;
            endy=numrows+1;
        }else if(j == 3){
            startx = 1;
            endx = 0;
            starty=numrows;
            endy = 1;
            dy = -1;
            dx = -1;
        }

        for(int y = starty; y != endy; y += dy){
            for(int x = startx; x != endx; x += dx){

                vtkSmartPointer< vtkQuad > quad = 0;
                quad = vtkSmartPointer< vtkQuad >::New();
                quad->GetPointIds()->SetId(0, pointsids[y][x]);
                quad->GetPointIds()->SetId(1, pointsids[y][x+dx]);
                quad->GetPointIds()->SetId(2, pointsids[y+dy][x+dx]);
                quad->GetPointIds()->SetId(3, pointsids[y+dy][x]);

                cellarray->InsertNextCell(quad);

            }
            //cout<<endl;
        }

    }

    vtkSmartPointer<vtkPolyData> poly = vtkSmartPointer<vtkPolyData>::New();

    poly->SetPoints(uvpoints);
    poly->SetPolys(cellarray);

    //cout<<poly->GetNumberOfCells()<<endl;
    vtkSmartPointer<vtkCellLocator> uvtoidlocator = vtkSmartPointer<vtkCellLocator>::New();
    uvtoidlocator->SetDataSet(poly);
    uvtoidlocator->AutomaticOn();
    uvtoidlocator->SetTolerance(0);
    //m_UVToIdLocator->SetNumberOfPointsPerBucket(3);
    uvtoidlocator->BuildLocator();
    m_CellLocators.push_back(uvtoidlocator);

    vtkSmartPointer<vtkPoints> uvpointsside = vtkSmartPointer<vtkPoints>::New();
    for(unsigned i = 0; i < poly->GetNumberOfPoints(); i++){

        double p[3];
        poly->GetPoint(i, p);

        double cg = p[1]+2 + (numrows-p[1])*2;

        uvpointsside->InsertNextPoint(p[0], cg, p[2]);

        //cout<<p[0]<<" "<<cg<<endl;

    }

    poly = vtkSmartPointer<vtkPolyData>::New();

    poly->SetPoints(uvpointsside);
    poly->SetPolys(cellarray);

    uvtoidlocator = vtkSmartPointer<vtkCellLocator>::New();
    uvtoidlocator->SetDataSet(poly);
    uvtoidlocator->AutomaticOn();
    uvtoidlocator->SetTolerance(0);
    //m_UVToIdLocator->SetNumberOfPointsPerBucket(3);
    uvtoidlocator->BuildLocator();
    m_CellLocators.push_back(uvtoidlocator);


    //cout<<"cids "<<m_Input->GetCrestMedialAtomsIds().size()<<endl;


    //handle the corners
    cellarray = vtkSmartPointer<vtkCellArray>::New();
    for(unsigned j = 0; j < 4; j++){

        int startx = 0;
        int endx = 0;

        int starty = 0;
        int endy = 0;

        int dx = 1;
        int dy = 1;

        //bool start = true;

        if(j == 0){
            startx = 1;
            endx = 0;
            starty = 1;
            endy = 0;
            dx = -1;
            dy = -1;
        }else if(j == 1){
            startx = numcols;
            endx = numcols+1;
            starty = 1;
            endy = 0;
            dy = -1;
        }else if(j == 2){
            startx = numcols;
            endx = numcols+1;
            dx = 1;
            starty=numrows;
            endy=numrows+1;
        }else if(j == 3){
            startx = 1;
            endx = 0;
            dx = -1;

            starty=numrows;
            endy = numrows+1;
            dy = 1;

        }

        for(int y = starty; y != endy; y += dy){
            for(int x = startx; x != endx; x += dx){

                vtkSmartPointer< vtkQuad > quad = 0;
                quad = vtkSmartPointer< vtkQuad >::New();
                quad->GetPointIds()->SetId(0, pointsids[y][x]);
                quad->GetPointIds()->SetId(1, pointsids[y][x+dx]);
                quad->GetPointIds()->SetId(2, pointsids[y+dy][x+dx]);
                quad->GetPointIds()->SetId(3, pointsids[y+dy][x]);

                cellarray->InsertNextCell(quad);

            }
            //cout<<endl;
        }

    }

    poly = vtkSmartPointer<vtkPolyData>::New();

    poly->SetPoints(uvpoints);
    poly->SetPolys(cellarray);

    //cout<<poly->GetNumberOfCells()<<endl;
    uvtoidlocator = vtkSmartPointer<vtkCellLocator>::New();
    uvtoidlocator->SetDataSet(poly);
    uvtoidlocator->AutomaticOn();
    uvtoidlocator->SetTolerance(0);
    //m_UVToIdLocator->SetNumberOfPointsPerBucket(3);
    uvtoidlocator->BuildLocator();
    m_CellLocators.push_back(uvtoidlocator);

    uvpointsside = vtkSmartPointer<vtkPoints>::New();
    for(unsigned i = 0; i < poly->GetNumberOfPoints(); i++){

        double p[3];
        poly->GetPoint(i, p);

        double cg = p[1]+2 + (numrows-p[1])*2;

        uvpointsside->InsertNextPoint(p[0], cg, p[2]);

        //cout<<p[0]<<" "<<cg<<endl;

    }

    poly = vtkSmartPointer<vtkPolyData>::New();

    poly->SetPoints(uvpointsside);
    poly->SetPolys(cellarray);

    uvtoidlocator = vtkSmartPointer<vtkCellLocator>::New();
    uvtoidlocator->SetDataSet(poly);
    uvtoidlocator->AutomaticOn();
    uvtoidlocator->SetTolerance(0);
    //m_UVToIdLocator->SetNumberOfPointsPerBucket(3);
    uvtoidlocator->BuildLocator();
    m_CellLocators.push_back(uvtoidlocator);
}


vtkIdType vtkSRepInterpolator::GetCellId(double& u, double& v, int &loc){

    double point[3] = {u, v, -1};
    double point2[3] = {u, v, 1};    

    for(unsigned i = 0; i < m_CellLocators.size(); i++){
        vtkSmartPointer<vtkIdList> cellsids = vtkSmartPointer<vtkIdList>::New();
        m_CellLocators[i]->FindCellsAlongLine(point, point2, 0, cellsids);

        if(cellsids->GetNumberOfIds() != 0){

            vtkIdType cellid = cellsids->GetId(0);
            vtkSmartPointer<vtkIdList> ptids = vtkSmartPointer<vtkIdList>::New();

            vtkPolyData* poly = dynamic_cast<vtkPolyData*>(m_CellLocators[i]->GetDataSet());

            poly->GetCellPoints(cellid, ptids);
            double p0[3];
            poly->GetPoint(ptids->GetId(0), p0);

            u = fabs(u - p0[0]);
            v = fabs(v - p0[1]);
            loc = i;
            return cellid;

        }
    }

    return -1;
}

/*
* \fn void GetInterpolatedPoint(double u, double v, double t)
* \brief Computes the point using the s-rep interpolators
* \param double u [0, 1]
* \param double v [0, 1]
* \param double t [0, 1]
* \pre u, v and t have to be between [0, 1]
*/
vtkSRepInterpolator::VNLType vtkSRepInterpolator::GetInterpolatedPoint(double u, double v, double tau){

    double ucoord = u*(m_Input->GetNumColumns()+1.0);
    double vcoord = v*((m_Input->GetNumRows()+1.0)*2.0);
    int loc = -1;

    vtkIdType cellid = GetCellId(ucoord, vcoord, loc);

    if(cellid != -1){
        if(loc < 2){
            vtkSRep::VNLType position = m_Medialsheetinterpolator->GetInterpolatedPoint(cellid, ucoord, vcoord);
            vtkSRep::VNLType spoke = m_Medialspokesinterpolator->GetInterpolatedSpoke(cellid, loc, ucoord, vcoord);
            m_InterpolatedLabel = m_Medialspokesinterpolator->GetInterpolatedLabel(ucoord, vcoord, cellid);

            return spoke*tau + position;
        }else if(loc == 2){

            double tempu = ucoord;
            double tempv = vcoord;

            double tempucoord = u*(m_Input->GetNumColumns()+1);
            if(tempucoord < 1 || tempucoord >= m_Input->GetNumColumns()){
                tempu = vcoord;
                tempv = ucoord;
            }

            vtkSRep::VNLType position = m_Interpolatecrestspokes->GetInterpolatedPoint(cellid, tempu);
            vtkSRep::VNLType spoke = m_Interpolatecrestspokes->GetInterpolatedSpoke(cellid, tempu, tempv/2.0);
            return spoke*tau + position;
        }else if(loc == 3){

            double tempu = ucoord;
            double tempv = vcoord;

            double tempucoord = u*(m_Input->GetNumColumns()+1);
            if(tempucoord < 1 || tempucoord >= m_Input->GetNumColumns()){
                tempu = vcoord;
                tempv = ucoord;
            }


            vtkSRep::VNLType position = m_Interpolatecrestspokes->GetInterpolatedPoint(cellid, tempu);
            vtkSRep::VNLType spoke = m_Interpolatecrestspokes->GetInterpolatedSpoke(cellid, tempu, 1 - tempv/2.0);
            return spoke*tau + position;
        }else if(loc == 4){


            double tempucoord = u*(m_Input->GetNumColumns()+1);
            double tempu = ucoord;
            double tempv = vcoord;
            if(tempucoord < 1 || tempucoord >= m_Input->GetNumColumns()){
                tempu = vcoord;
                tempv = ucoord;
            }

            vtkSRep::VNLType position;
            vtkSRep::VNLType spoke;

            if(cellid == 0){

                position = m_Interpolatecrestspokes->GetInterpolatedPoint(0, 0);

                double temp = sqrt(pow(tempu, 2) + pow(tempv, 2));
                spoke = m_Interpolatecrestspokes->GetInterpolatedSpoke(0, 0, temp/2.0);

            }else if(cellid == 1){

                cellid = m_Input->GetNumColumns() - 1;

                position = m_Interpolatecrestspokes->GetInterpolatedPoint(cellid, 0);
                double temp = sqrt(pow(tempu, 2) + pow(tempv, 2));
                spoke = m_Interpolatecrestspokes->GetInterpolatedSpoke(cellid, 0, temp/2.0);
            }else if(cellid == 2){

                cellid = m_Input->GetNumColumns() - 1 + m_Input->GetNumRows() - 1;

                position = m_Interpolatecrestspokes->GetInterpolatedPoint(cellid, 0);
                double temp = sqrt(pow(tempu, 2) + pow(tempv, 2));
                spoke = m_Interpolatecrestspokes->GetInterpolatedSpoke(cellid, 0, temp/2.0);
            }else{
                cellid = 2*(m_Input->GetNumColumns() - 1) + m_Input->GetNumRows() - 1;
                position = m_Interpolatecrestspokes->GetInterpolatedPoint(cellid, 0);

                double temp = sqrt(pow(tempu, 2) + pow(tempv, 2));
                spoke = m_Interpolatecrestspokes->GetInterpolatedSpoke(cellid, 0, temp/2.0);
            }


            return spoke*tau + position;

        }else if(loc == 5){

            double tempucoord = u*(m_Input->GetNumColumns()+1);
            double tempu = ucoord;
            double tempv = vcoord;
            if(tempucoord < 1 || tempucoord >= m_Input->GetNumColumns()){
                tempu = vcoord;
                tempv = ucoord;
            }

            vtkSRep::VNLType position;
            vtkSRep::VNLType spoke;

            if(cellid == 0){

                position = m_Interpolatecrestspokes->GetInterpolatedPoint(0, 0);

                double temp = sqrt(pow(tempu, 2) + pow(tempv, 2));
                spoke = m_Interpolatecrestspokes->GetInterpolatedSpoke(0, 0, 1- temp/2.0);

            }else if(cellid == 1){


                cellid = m_Input->GetNumColumns() - 1;

                position = m_Interpolatecrestspokes->GetInterpolatedPoint(cellid, 0);
                double temp = sqrt(pow(tempu, 2) + pow(tempv, 2));
                spoke = m_Interpolatecrestspokes->GetInterpolatedSpoke(cellid, 0, 1 - temp/2.0);
            }else if(cellid == 2){

                cellid = m_Input->GetNumColumns() - 1 + m_Input->GetNumRows() -1;

                position = m_Interpolatecrestspokes->GetInterpolatedPoint(cellid, 0);
                double temp = sqrt(pow(tempu, 2) + pow(tempv, 2));
                spoke = m_Interpolatecrestspokes->GetInterpolatedSpoke(cellid, 0, 1 - temp/2.0);
            }else{
                cellid = m_Input->GetNumColumns() - 1 + m_Input->GetNumRows() - 1 + m_Input->GetNumColumns() - 1;
                position = m_Interpolatecrestspokes->GetInterpolatedPoint(cellid, 0);
                double temp = sqrt(pow(tempu, 2) + pow(tempv, 2));
                spoke = m_Interpolatecrestspokes->GetInterpolatedSpoke(cellid, 0, 1 - temp/2.0);
            }

            return spoke*tau + position;
        }
    }

    VNLType p(3);
    p.fill(0);

    cout<<"not"<<endl;

    return p;

}

vtkSRepInterpolator::VNLType vtkSRepInterpolator::GetInterpolatedSpoke(double u, double v){

    double ucoord = u*(m_Input->GetNumColumns()+1.0);
    double vcoord = v*((m_Input->GetNumRows()+1.0)*2.0);
    int loc = -1;

    vtkIdType cellid = GetCellId(ucoord, vcoord, loc);

    if(cellid != -1){
        if(loc < 2){

            vtkSRep::VNLType spoke = m_Medialspokesinterpolator->GetInterpolatedSpoke(cellid, loc, ucoord, vcoord);
            return spoke;

        }else if(loc == 2){

            double tempu = ucoord;
            double tempv = vcoord;

            double tempucoord = u*(m_Input->GetNumColumns()+1);
            if(tempucoord < 1 || tempucoord >= m_Input->GetNumColumns()){
                tempu = vcoord;
                tempv = ucoord;
            }

            vtkSRep::VNLType spoke = m_Interpolatecrestspokes->GetInterpolatedSpoke(cellid, tempu, tempv/2.0);
            return spoke;
        }else if(loc == 3){

            double tempu = ucoord;
            double tempv = vcoord;

            double tempucoord = u*(m_Input->GetNumColumns()+1);
            if(tempucoord < 1 || tempucoord >= m_Input->GetNumColumns()){
                tempu = vcoord;
                tempv = ucoord;
            }

            vtkSRep::VNLType spoke = m_Interpolatecrestspokes->GetInterpolatedSpoke(cellid, tempu, 1 - tempv/2.0);
            return spoke;
        }else if(loc == 4){


            double tempucoord = u*(m_Input->GetNumColumns()+1);
            double tempu = ucoord;
            double tempv = vcoord;
            if(tempucoord < 1 || tempucoord >= m_Input->GetNumColumns()){
                tempu = vcoord;
                tempv = ucoord;
            }

            vtkSRep::VNLType spoke;

            if(cellid == 0){

                double temp = sqrt(pow(tempu, 2) + pow(tempv, 2));
                spoke = m_Interpolatecrestspokes->GetInterpolatedSpoke(0, 0, temp/2.0);

            }else if(cellid == 1){

                cellid = m_Input->GetNumColumns() - 1;
                double temp = sqrt(pow(tempu, 2) + pow(tempv, 2));
                spoke = m_Interpolatecrestspokes->GetInterpolatedSpoke(cellid, 0, temp/2.0);
            }else if(cellid == 2){

                cellid = m_Input->GetNumColumns() - 1 + m_Input->GetNumRows() - 1;
                double temp = sqrt(pow(tempu, 2) + pow(tempv, 2));
                spoke = m_Interpolatecrestspokes->GetInterpolatedSpoke(cellid, 0, temp/2.0);
            }else{
                cellid = 2*(m_Input->GetNumColumns() - 1) + m_Input->GetNumRows() - 1;
                double temp = sqrt(pow(tempu, 2) + pow(tempv, 2));
                spoke = m_Interpolatecrestspokes->GetInterpolatedSpoke(cellid, 0, temp/2.0);
            }


            return spoke;

        }else if(loc == 5){

            double tempucoord = u*(m_Input->GetNumColumns()+1);
            double tempu = ucoord;
            double tempv = vcoord;
            if(tempucoord < 1 || tempucoord >= m_Input->GetNumColumns()){
                tempu = vcoord;
                tempv = ucoord;
            }

            vtkSRep::VNLType spoke;

            if(cellid == 0){

                double temp = sqrt(pow(tempu, 2) + pow(tempv, 2));
                spoke = m_Interpolatecrestspokes->GetInterpolatedSpoke(0, 0, 1- temp/2.0);

            }else if(cellid == 1){


                cellid = m_Input->GetNumColumns() - 1;
                double temp = sqrt(pow(tempu, 2) + pow(tempv, 2));
                spoke = m_Interpolatecrestspokes->GetInterpolatedSpoke(cellid, 0, 1 - temp/2.0);
            }else if(cellid == 2){

                cellid = m_Input->GetNumColumns() - 1 + m_Input->GetNumRows() -1;
                double temp = sqrt(pow(tempu, 2) + pow(tempv, 2));
                spoke = m_Interpolatecrestspokes->GetInterpolatedSpoke(cellid, 0, 1 - temp/2.0);
            }else{
                cellid = m_Input->GetNumColumns() - 1 + m_Input->GetNumRows() - 1 + m_Input->GetNumColumns() - 1;
                double temp = sqrt(pow(tempu, 2) + pow(tempv, 2));
                spoke = m_Interpolatecrestspokes->GetInterpolatedSpoke(cellid, 0, 1 - temp/2.0);
            }

            return spoke;
        }
    }

    VNLType p(3);
    p.fill(0);

    cout<<"not"<<endl;

    return p;

}
