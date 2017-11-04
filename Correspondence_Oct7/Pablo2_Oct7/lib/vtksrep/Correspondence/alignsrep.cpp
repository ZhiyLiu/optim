/* This is a preprocess of correspondence, because s-reps may not have a same size or position, we
 * align(translation, rotation and scaling) them to be similar.
 * We first do Procrustes alignment to the boundaries of the correspondence sreps, and then perform the same translation to
 * its skeletal sheet.
 * We sample each srep's up boundary, down boundary and crest boundary together, and do alignment at the same time.
 * Author: Liyun Tu
 * Date: Mar 1, 2014.
*/



#include "alignsrep.h"
#include <iostream>
#include "visualization.h"

using namespace std;


alignsrep::alignsrep()
{
}


/*alignsrep::alignsrep(const char* srepFolder)
{
    this->srepFolder = string(srepFolder);
}*/


alignsrep::alignsrep(const char* filename/*, M3DQuadFigure* quadfig, */,int interpolationLevel/*, string srepFolder*/){
    this->quadFig = tls.GetQuadFigure(filename);
    this->rowNums = this->quadFig->getRowCount();
    this->colNums = this->quadFig->getColumnCount();
    this->quadNum = (this->rowNums-1)*(this->colNums-1);
    this->crestAtomNums =  (this->rowNums-2)*2 + this->colNums*2;//number of crest atoms.
    this->crestQuadNum = (this->rowNums-1)*2 + (this->colNums-1)*2;//Equal to crestAtomNums ??
    this->totalAtomNums = this->rowNums * this->colNums;
    this->interiorAtomNums = this->totalAtomNums - this->crestAtomNums;
    this->interpolationLevel = interpolationLevel;
    this->step = pow((double)2, (double)interpolationLevel);
    this->subQuadPointsNum = (step+1)*(step+1);

    //this->srepFolder = string(srepFolder);


    // Initialize srepfig.
    this->readsrep = vtkSmartPointer<vtkReadSRep>::New();
    readsrep->SetFileName(filename);
    readsrep->Update();

    this->srepfig = readsrep->GetOutput();
}



/* srepFolder: directory where a list of s-reps which we will performance alignment on.
*/
/*int alignsrep::prepareDataForMatlabProcr(){
    cout<<"------------ Prepare data for Procrustes alignment: "<<srepFolder<<" --------------"<<endl;

    if(!srepPathNames.empty()){
        cout<<"-------------there are: "<<srepPathNames.size()<<"sreps."<<endl;
        for(int q =0; q<srepPathNames.size();q++){
            this->figfilename = srepPathNames[q];
            cout<<"Begain to read in srep: "<<this->figfilename<<endl;

            this->quadFig = tls.GetQuadFigure(this->figfilename.c_str());

            M3DQuadPrimitive* prim;

            this->pRowIndex = 0;

            //Loop each primitive, get the x,y,z coordinate of this srep's boundary.
            for(unsigned u =0; u<this->quadFig->getRowCount(); u++){
                for(unsigned v =0; v<this->quadFig->getColumnCount(); v++){
                    prim = dynamic_cast<M3DQuadPrimitive*>(this->quadFig->getPrimitivePtr(u,v));

                    //pMatrix holding 3*atomNum row X q columns data.
                    this->pMatrix(this->pRowIndex,q) = prim->getX().getX();
                    this->pMatrix(this->pRowIndex+1,q) = prim->getX().getY();
                    this->pMatrix(this->pRowIndex+2,q) = prim->getX().getZ();


                    this->pRowIndex = this->pRowIndex+3;
                }
            }
        }

        // save pMatrix to file.
        string filename = string(this->srepFolder) + string("GPA_input_pMatrix.txt");
        saveGPA_Input_Matrix(filename.c_str());

        return EXIT_SUCCESS;
    }
    else{
        cout<<"You input a invalid srep file, please check the path or name of this srep!"<<endl;
        return EXIT_FAILURE;
    }
}*/





/* Save matrix to .txt file.*/
void alignsrep::saveMatrix(vector<Vector3D> points,const char* filename){

    std::ofstream fout;
    fout.open(filename);

    if(fout)  {
        for(int i =0; i< points.size();i++){
            fout<< points[i].getX()<<" "<<points[i].getY()<<" "<<points[i].getZ()<<endl;
        }
        //fout<<endl;

        cout<<"Successfully saved matrix to: "<<filename<<endl;
    }
    else
        cerr<<"Write out failed, cannot open the file!"<<endl;

    fout.close();
}


/* Get all the points on boundary. Including: up boundary points, down boundary points, crest boundray points.
 * We store these points in 5 vectors:
 * up_b: store all the points on the top boundary.
 * down_b: store all the points on bottom boundary.
 * up_b_crest: store the crestAtomNums points on the top boundary.
 * down_b_crest: store the crestAtomNums points on the bottom boundary.
 * medial_crest: store the crestAtomNums points on the medial boundary.
 * For the crest points, the sequence is: primitive[0,0], [0,1], ..., [0,colNums-1], [1,colNums-1], [2,colNums-1], ..., [rowNums-1,colNums-1],
 * [rowNums-1,colNums-2], [rowNums-1,colNums-3], ..., [rowNums-1,2], [rowNums-1,1], [rowNums-1,0], [rowNums-2,0], ..., [1,0].
 * up_b_crest, down_b_crest and medial_crest have a same size of crestAtomNums. There points id are correspondence.
 * up_b_mid and down_b_mid have a same size of interiorAtomNums, there points id are correpondence.
 *
*/
void alignsrep::getPointsOnBoundary(){

    M3DQuadPrimitive* prim;
    M3DQuadEndPrimitive* endPrim;
    Vector3D up_b_point, down_b_point, up_b_crest_point, down_b_crest_point, crest_b_point;

    // Loop the interior atoms, which only have up and down spokes.
    for(unsigned u =1; u<this->rowNums-1; u++){
        for(unsigned v =1; v<this->colNums-1; v++){
            prim = dynamic_cast<M3DQuadPrimitive*>(this->quadFig->getPrimitivePtr(u,v));
            // Up boundary points.
            up_b_point = prim->getX() + prim->getR0()*prim->getU0();
            this->up_b_mid.push_back(up_b_point);

            // Down boundary points
            down_b_point = prim->getX() + prim->getR1()*prim->getU1();
            this->down_b_mid.push_back(down_b_point);
        }
    }

    // Get the crest curve along the crest in counter clockwise. the start point is (0,0).
    for(unsigned v =0; v<colNums; v++){
        // Up boundary crest points. (End atom's up spoke's tip on the boundary). 2 quads use IM1, 2 quads use IM2.
        prim = dynamic_cast<M3DQuadPrimitive*>(this->quadFig->getPrimitivePtr(0,v));
        up_b_crest_point = prim->getX() + prim->getR0()*prim->getU0();        
        this->up_b_crest.push_back(up_b_crest_point);

        // Down boundary crest points. (End atom's down spoke's tip on the boundary). 2 quads use IM1, 2 quads use IM2.
        down_b_crest_point = prim->getX() + prim->getR1()*prim->getU1();
        this->down_b_crest.push_back(down_b_crest_point);

        // Crest boundary points. (crest spokes' tip on the boundary.)
        endPrim = dynamic_cast<M3DQuadEndPrimitive*>(this->quadFig->getPrimitivePtr(0, v));
        crest_b_point = endPrim->getX() + endPrim->getREnd()*endPrim->getUEnd();        
        this->medial_crest.push_back(crest_b_point);
    }
    for(unsigned u =1; u<rowNums; u++){
        // Up boundary crest points. (End atom's up spoke's tip on the boundary). 2 quads use IM1, 2 quads use IM2.
        prim = dynamic_cast<M3DQuadPrimitive*>(this->quadFig->getPrimitivePtr(u, colNums - 1));
        up_b_crest_point = prim->getX() + prim->getR0()*prim->getU0();
        this->up_b_crest.push_back(up_b_crest_point);

        // Down boundary crest points. (End atom's down spoke's tip on the boundary). 2 quads use IM1, 2 quads use IM2.
        down_b_crest_point = prim->getX() + prim->getR1()*prim->getU1();
        this->down_b_crest.push_back(down_b_crest_point);

        // Crest boundary points. (crest spokes' tip on the boundary.)
        endPrim = dynamic_cast<M3DQuadEndPrimitive*>(this->quadFig->getPrimitivePtr(u, colNums - 1));
        crest_b_point = endPrim->getX() + endPrim->getREnd()*endPrim->getUEnd();
        this->medial_crest.push_back(crest_b_point);
    }
    for(unsigned v =colNums-2; v>0; v--){
        // Up boundary crest points. (End atom's up spoke's tip on the boundary). 2 quads use IM1, 2 quads use IM2.
        prim = dynamic_cast<M3DQuadPrimitive*>(this->quadFig->getPrimitivePtr(rowNums-1, v));
        up_b_crest_point = prim->getX() + prim->getR0()*prim->getU0();
        this->up_b_crest.push_back(up_b_crest_point);

        // Down boundary crest points. (End atom's down spoke's tip on the boundary). 2 quads use IM1, 2 quads use IM2.
        down_b_crest_point = prim->getX() + prim->getR1()*prim->getU1();
        this->down_b_crest.push_back(down_b_crest_point);

        // Crest boundary points. (crest spokes' tip on the boundary.)
        endPrim = dynamic_cast<M3DQuadEndPrimitive*>(this->quadFig->getPrimitivePtr(rowNums-1, v));
        crest_b_point = endPrim->getX() + endPrim->getREnd()*endPrim->getUEnd();
        this->medial_crest.push_back(crest_b_point);
    }
    for(unsigned u =rowNums-1; u>0; u--){
        // Up boundary crest points. (End atom's up spoke's tip on the boundary). 2 quads use IM1, 2 quads use IM2.
        prim = dynamic_cast<M3DQuadPrimitive*>(this->quadFig->getPrimitivePtr(u, 0));
        up_b_crest_point = prim->getX() + prim->getR0()*prim->getU0();
        this->up_b_crest.push_back(up_b_crest_point);

        // Down boundary crest points. (End atom's down spoke's tip on the boundary). 2 quads use IM1, 2 quads use IM2.
        down_b_crest_point = prim->getX() + prim->getR1()*prim->getU1();
        this->down_b_crest.push_back(down_b_crest_point);

        // Crest boundary points. (crest spokes' tip on the boundary.)
        endPrim = dynamic_cast<M3DQuadEndPrimitive*>(this->quadFig->getPrimitivePtr(u, 0));
        crest_b_point = endPrim->getX() + endPrim->getREnd()*endPrim->getUEnd();
        this->medial_crest.push_back(crest_b_point);
    }
}








/* Get each quad area and store them in tables.
 * Each quad's area is compute by given a interpolationLevel, and sum the subquads' areas together.
 * The interpolationLevel used in dividing quads into subquads and sum its area can be any level.
 * Get the interpolationLevel=0 boundary points and store them in tables.
 * Apply each boundary point to its neighbor 4 or 3 quads.
*/
void alignsrep::weightAreaToBoundaryPoints(){

    // Step 1: Get boundary points(The interpolationLevel=0 quads' vertexes).
    // Store the points into talbes: up_b, down_b, up_b_crest, down_b_crest and medial_crest.
    getPointsOnBoundary();


    // Step 2: Get areas for each quad, store in tables.
    // Get the up and down boundary quad vertexes positions, store them in up_BP and down_BP seperately.
    getAllSubQuadsVertexesOnBoundary();

    // Areas of up and down quads.    
    up_quads_areas = calculateAreas_method2(this->up_BP);
    down_quads_areas = calculateAreas_method2(this->down_BP);

    // Areas of crest quads. up_b_crest_areas store the area of quads between up boundary crest and medial crest;
    // down_b_crest_areas store area of quads between down boundary crest and medial crest.
    //getCrestQuadsArea_InteLevel0();
    getCrestQuadsArea();


    // Step 3: weight each boundary point with its neighbor 4 or 3 quad area.
    this->pointIndex = 0; // index of the boundary point.
    // Add the up boundary points and its correspondence areas.(Including up middle points and up crest points.)
    correPointAndArea(up_quads_areas, this->up_b_mid, this->up_b_crest_areas, this->up_b_crest);

    // Add the down boundary points and its correspondence areas.(Including down middle points and down crest points.)
    correPointAndArea(down_quads_areas, this->down_b_mid, this->down_b_crest_areas, this->down_b_crest);

    // Add the crest boundary points and its correspondence areas.
    // For the first point in crest boundary
    this->boundaryPoints.push_back(this->medial_crest[0]);
    double area1, area2, area3, area4;
    area1 = this->down_b_crest_areas[this->medial_crest.size()-1];
    area2 = this->down_b_crest_areas[0];
    area3 = this->up_b_crest_areas[0];
    area4 = this->up_b_crest_areas[this->medial_crest.size()-1];
    addBoundaryPointsAreas(area1, area2, area3, area4);
    // For the rest points in crest boundary.
    for(unsigned i =1; i<this->medial_crest.size();i++){
        this->boundaryPoints.push_back(this->medial_crest[i]);
        double area1, area2, area3, area4;
        area1 = this->down_b_crest_areas[i-1];
        area2 = this->down_b_crest_areas[i];
        area3 = this->up_b_crest_areas[i];
        area4 = this->up_b_crest_areas[i-1];
        addBoundaryPointsAreas(area1, area2, area3, area4);
    }

    // Output vector for check...
    /*cout<<"this->boundaryPoints.size() is: "<<this->boundaryPoints.size()<<endl;
    for(unsigned i=0; i<this->boundaryPoints.size();i++){
        cout<<"boundary point["<<i<<"] is: "<<this->boundaryPoints[i]<<endl;
    }
    cout<<"this->boundaryPointsAreas.size() is: "<<this->boundaryPointsAreas.size()<<endl;
    for(unsigned i=0; i<this->boundaryPointsAreas.size();i++){
        cout<<"boundaryPointsAreas["<<i<<"] is: "<<this->boundaryPointsAreas[i][0]<<" "<<this->boundaryPointsAreas[i][1]<<" "
           <<this->boundaryPointsAreas[i][2]<<" "<<this->boundaryPointsAreas[i][3]<<endl;
    }*/


    // Step 4: weight points using correspondence 4 areas.
    //vector<Vector3D> wPoints = weightPoints();//incorrect logic, don't use.
    // Step 4: call alignsrep class from arepAlignment.cxx, generate boundary points and area for each srep...


    // Step 5: Save weighted boundary points to file for each srep. This used just for test.
    /*toolsfunc tl;
    string filename = tl.getFileNameWithExtension(this->figfilename);
    cout<<"--------------filename is: "<<filename<<endl;
    string tmp = tl.splitExtension(filename);
    cout<<"--------------filename without exetion is: "<<tmp<<endl;
    string str = string(this->srepFolder) + tmp + string("_wPoints.txt");
    cout<<"--------------whole path of filename is: "<<str<<endl;
    saveMatrix(wPoints, str.c_str());*/

}



/* Save area into vector. */
void alignsrep::addBoundaryPointsAreas(double area1, double area2, double area3, double area4){

    this->boundaryPointsAreas.push_back(VectorSRepFeaturesType());

    this->boundaryPointsAreas[this->pointIndex].push_back(area1);
    this->boundaryPointsAreas[this->pointIndex].push_back(area2);
    this->boundaryPointsAreas[this->pointIndex].push_back(area3);
    this->boundaryPointsAreas[this->pointIndex].push_back(area4);

    this->pointIndex++;
}


/* Loop each boundary points and save in vector this->boundaryPoints.
 * Save each boundary points' 4 areas in vector this->boundaryPointsAreas, it has same index with the boundaryPoints vector.
 * quads_areas: a vector storing all the up or down boundary quads areas.
 * mid_points: a vector storing the middle points on up or down boundary.
 * crest_areas: storing up or down crest boudary quad areas.
 * crest_points: storing up or down crest boundary points.
*/
void alignsrep::correPointAndArea(vector<double> quads_areas, vector<Vector3D> mid_points, vector<double> crest_areas,
                                           vector<Vector3D> crest_points){
    int quadNumEachRow = this->colNums-1; // quad number each row
    double area1, area2, area3, area4;
    // For interior tips
    int pindex =0;
    for(unsigned n=0; n<this->rowNums-2;n++){
        for(unsigned m=0; m<this->colNums-2; m++){
            this->boundaryPoints.push_back(mid_points[pindex]);

            area1 = quads_areas[quadNumEachRow*n + m];
            area2 = quads_areas[quadNumEachRow*n + m+1];
            area3 = quads_areas[quadNumEachRow*n + m+1+quadNumEachRow];
            area4 = quads_areas[quadNumEachRow*n + m+quadNumEachRow];
            /*cout<<"-----"<<quadNumEachRow*n + m<<"-----  "<<area1<<endl;
            cout<<"-----"<<quadNumEachRow*n + m+1<<"-----  "<<area2<<endl;
            cout<<"-----"<<quadNumEachRow*n + m+1+quadNumEachRow<<"-----  "<<area3<<endl;
            cout<<"-----"<<quadNumEachRow*n + m+quadNumEachRow<<"-----  "<<area4<<endl;*/

            addBoundaryPointsAreas(area1, area2, area3, area4);

            pindex++;
        }
    }

    // For up boundary crest tips
    if(this->crestQuadNum!=crest_points.size() || this->crestQuadNum!= crest_areas.size()
            || crest_points.size()!= crest_areas.size()){
        cout<<"Message from alignsrep.cpp: the crest quad number should equals to the up_b_crest_areas vector's size"<<endl;
        return;
    }
    // For the first tip
    this->boundaryPoints.push_back(crest_points[0]);
    area1 = crest_areas[this->crestQuadNum -1];
    area2 = crest_areas[0];
    area3 = quads_areas[0];
    area4 = 0;
    addBoundaryPointsAreas(area1, area2, area3, area4);

    //For the following colNums-2 tips.
    int k; // Up crest index;
    for(k=1;k<quadNumEachRow;k++){
        this->boundaryPoints.push_back(crest_points[k]);

        area1 = crest_areas[k-1];
        area2 = crest_areas[k];
        area3 = quads_areas[k];
        area4 = quads_areas[k-1];

        addBoundaryPointsAreas(area1, area2, area3, area4);
    }

    // For the (colNums-1)th tip, the bottom-right corner.
    k = this->colNums-1;
    this->boundaryPoints.push_back(crest_points[k]); // k = this->colNums-1;
    area1 = crest_areas[k-1];
    area2 = crest_areas[k];
    area3 = quads_areas[k-1];
    area4 = 0;
    addBoundaryPointsAreas(area1, area2, area3, area4);
    k++;

    // For the following rowNum-2 tips
    for(unsigned j=1; j<this->rowNums-1;j++){
        this->boundaryPoints.push_back(crest_points[k]);

        area1 = crest_areas[k-1];
        area2 = crest_areas[k];
        area3 = quads_areas[quadNumEachRow*j+(this->colNums-2)];
        area4 = quads_areas[quadNumEachRow*(j-1)+(this->colNums-2)];

        addBoundaryPointsAreas(area1, area2, area3, area4);
        k++;
    }

    // For the (rowNum+colNums-1)th tip, the top-right corner.
    this->boundaryPoints.push_back(crest_points[k]);
    area1 = crest_areas[k-1];
    area2 = crest_areas[k];
    area3 = quads_areas[this->quadNum-1];
    area4 = 0;
    addBoundaryPointsAreas(area1, area2, area3, area4);
    k++;

    // For the following colNum-1 tips.
    for(int i =1; i<quadNumEachRow;i++){
        this->boundaryPoints.push_back(crest_points[k]);

        area1 = crest_areas[k-1];
        area2 = crest_areas[k];
        area3 = quads_areas[this->quadNum-1-(i-1)];
        area4 = quads_areas[this->quadNum-1-i];

        addBoundaryPointsAreas(area1, area2, area3, area4);
        k++;
    }

    // For the top-left corner tip
    this->boundaryPoints.push_back(crest_points[k]);
    area1 = crest_areas[k-1];
    area2 = crest_areas[k];
    area3 = quads_areas[this->quadNum-quadNumEachRow];
    area4 = 0;
    addBoundaryPointsAreas(area1, area2, area3, area4);
    k++;

    // For the following tips
    for(int i =this->rowNums-2; i>0; i--){
        this->boundaryPoints.push_back(crest_points[k]);

        area1 = crest_areas[k-1];
        area2 = crest_areas[k];
        area3 = quads_areas[(i-1)*quadNumEachRow];
        area4 = quads_areas[i*quadNumEachRow];

        addBoundaryPointsAreas(area1, area2, area3, area4);
        k++;
    }
    //cout<<"Here, k should equals to 20------: "<< k<<endl;
}



/* Weight each point with its correspondence neighbor 4 quads areas.
*/
vector<Vector3D> alignsrep::weightPoints(){

    // Sum all the boundary quads area.
    double totalArea;
    vector<double> weightArea;
    vector<Vector3D> wightedPoints;

    /*for(unsigned i=0; i<this->boundaryPointsAreas.size();i++){
        double sum4Areas =0;//sum the 4 quads area for each point.
        for(unsigned j=0; j<this->boundaryPointsAreas[0].size();j++){
            totalArea += this->boundaryPointsAreas[i][j];
            sum4Areas += this->boundaryPointsAreas[i][j];
        }
        weightArea.push_back(sum4Areas);
    }

    // Divid each points sum area by the total areas.
    if(weightArea.size()!=this->boundaryPoints.size() || weightArea.size()!=this->boundaryPointsAreas.size()){
        cout<<"Message from alignsrep::weightPoints: the points number should equals to the area weight!!"<<endl;
    }
    for(unsigned i=0; i<weightArea.size();i++){
        wightedPoints.push_back(Vector3D());
        weightArea[i] = weightArea[i]/totalArea;
        //cout<<"=============the area weight for each point is: "<<weightArea[i]<<endl;

        //cout<<"-------------the boundary point is: "<<this->boundaryPoints[i]<<endl;

        // Mutiply each point with its correspondence area weight
        wightedPoints[i] = this->boundaryPoints[i]*weightArea[i];
        //cout<<"************** the weighted point is: "<< wightedPoints[i]<<endl;
    }*/

    return wightedPoints;
}





/* Get all the sub-quad position by interpolate to u v coordinate. *
 * quadtype: boundary quads(0), skeletal quads(1).
*/
void alignsrep::getAllSubQuadsVertexesOnBoundary(){

    int quadIndex = 0;

    vector<double> subquadpoint_v, subquadpoint_u;
    VectorQuadPoint quadpoints_u;   //store the u coordinate of all the sub-quad of quad[q].
    VectorQuadPoint quadpoints_v;   //store the v coordinate of all the sub-quad of quad[q].

    // Interpolate to the u v coordinate.
    for(unsigned i = 0; i < this->rowNums -1; i++){ //the row number of the quads. its 3.
        for(unsigned j = 0; j < this->colNums -1; j++){//coloums, its 13.
                quadpoints_u.push_back(VectorDoublePoints());
                quadpoints_v.push_back(VectorDoublePoints());

                //four points of a quad with delta u, v. counter clockwise.
                double quadpointsu[4], quadpointsv[4];
                //left-top point
                quadpointsu[0] = i;
                quadpointsv[0] = j;
                //left-bottom point
                quadpointsu[1] = i + 1;
                quadpointsv[1] = j;
                //right-bottom point
                quadpointsu[2] = i + 1;
                quadpointsv[2] = j + 1;
                //right-top point
                quadpointsu[3] = i;
                quadpointsv[3] = j + 1;

                //given four points of a quad (p0, p1, p2, p3), split each side of the quad into 2^interpolationLevel sub-line.
                //first, get the subquads's u coordinate, in column first order.
                subquadpoint_u = this->tls.splitQuad(quadpointsu[0],quadpointsu[1],quadpointsu[2],quadpointsu[3], this->step);

                //second, get the 25 subquads's v coordinate, in column first order.
                subquadpoint_v = this->tls.splitQuad(quadpointsv[0],quadpointsv[1],quadpointsv[2],quadpointsv[3], this->step);

                for(unsigned m =0; m<subquadpoint_u.size();m++){
                    quadpoints_u[quadIndex].push_back(subquadpoint_u[m]);
                    quadpoints_v[quadIndex].push_back(subquadpoint_v[m]);

                    //cout<<"MSG from getSubQuadsUVCoordinateNotUsingDeltaUV: quadpoints_u["<<quadIndex<<"]---------m is: "<<m<<"---subquadpoint_u["<<m<<"] is: "<<subquadpoint_u[m]<<endl;
                    //cout<<"MSG from getSubQuadsUVCoordinateNotUsingDeltaUV: quadpoints_v["<<quadIndex<<"]---------m is: "<<m<<"---subquadpoint_v["<<m<<"] is: "<<subquadpoint_v[m]<<endl;
                }

                quadIndex++;

                subquadpoint_v.clear();
                subquadpoint_u.clear();
        }
    }

    M3DQuadInterpolater *tpm = new M3DQuadInterpolater(this->quadFig);

    //store all the positions of points for a quad, each point is a 3D vector, which store the x, y, z coordinate of this point.
    Vector3D point;

    for(int i =0; i< this->quadNum; i++){
        for(int j =0; j< this->subQuadPointsNum; j++){
            //cout<<"quad: "<<i<<", point: "<<j<<" u is: "<<quadpoints_u[i][j]<< " v is: "<<quadpoints_v[i][j]<<endl;
            point = tpm->interpolateQuadSpoke(this->quadFig,quadpoints_u[i][j],quadpoints_v[i][j],0)->getB();
            //cout<<"point position: x is: "<<point.getX()<<" y is: "<<point.getY()<<" z is: "<<point.getZ()<<endl;
            up_BP.push_back(point);

            // Down boundary
            point = tpm->interpolateQuadSpoke(this->quadFig,quadpoints_u[i][j],quadpoints_v[i][j],1)->getB();
            down_BP.push_back(point);
        }
    }

    delete tpm;
}



/* Given all the quads vertexes points on the surface, compute each quad area and store into a vector.
*/
vector<double> alignsrep::calculateAreas_method2(vector<Vector3D> surfacePoints){
    int subQuadNum = step*step;

    vector<double> subquadareas;//store the 16 subquads area.
    //store the sum of 16 subquad to get the quad area.
    vector<double> quadarea;

    Vector3D points[4];
    double subquadarea;

    //for each quad
    for(int n=0; n<quadNum;n++){
        //locate space for a new value, initial it to 0.
        quadarea.push_back(0);

        //For each subQuadNum's sub-quad.
        for(int m=0; m<step; m++){
            for(int k = 0; k<step; k++){
                int currentpoint = k + m*(step+1);
                //Left-top point.
                points[0] = surfacePoints[n*subQuadPointsNum+currentpoint]; //m
                //Bottom side
                points[1] = surfacePoints[n*subQuadPointsNum+currentpoint+1]; //m+1
                //Right side
                points[2] = surfacePoints[n*subQuadPointsNum+currentpoint+1+step+1]; //m+6
                //Top side
                points[3] = surfacePoints[n*subQuadPointsNum+currentpoint+step+1]; //m+5

                //each of the 16 subquads area.
                subquadarea = quadArea(points);

                //put all the 16 sub areas into a vector subquadareas.
                subquadareas.push_back(subquadarea);
            }
        }

        //sum the 16 subquads areas.
        for(unsigned i =0; i<subQuadNum;i++){
            quadarea[n] += subquadareas[i];
        }

        //clear vector for next quad.
        subquadareas.clear();
    }

    return quadarea;
}



/* Calculate the quad area by divide each of its subquads into two triangles.
 * Input: Four points of a quad, each point contains the x, y, z coordinate.
 * The four vetex point of the quad is in counter clockwise. whiche mean: point[0] is the top-left corner;
 * point[1] is the bottom-left corner; point[2] is the bottom-right corner; point[3] is the top-right corner.
 * There're two diagonal for each subquad, use each to calculate the sum of its two triangles. we can get two area,
 * we want the smaller one as the area of this subquad.
*/
double alignsrep::quadArea(Vector3D *point){
    //the first method: diagonal of p0p2
    double diagonalp0p2 = sqrt(pow(point[0].getX()-point[2].getX(), 2) + pow(point[0].getY()-point[2].getY(), 2) + pow(point[0].getZ()-point[2].getZ(), 2));
    //the second method: diagonal of p1p3
    double diagonalp1p3 = sqrt(pow(point[1].getX()-point[3].getX(), 2) + pow(point[1].getY()-point[3].getY(), 2) + pow(point[1].getZ()-point[3].getZ(), 2));

    //the four sides of each sub-quad.
    //left edge: p0p1
    double p0p1 = sqrt(pow(point[0].getX()-point[1].getX(), 2) + pow(point[0].getY()-point[1].getY(), 2) + pow(point[0].getZ()-point[1].getZ(), 2));
    //bottom edge: p1p2
    double p1p2 = sqrt(pow(point[1].getX()-point[2].getX(), 2) + pow(point[1].getY()-point[2].getY(), 2) + pow(point[1].getZ()-point[2].getZ(), 2));
    //right edge: p2p3
    double p2p3 = sqrt(pow(point[2].getX()-point[3].getX(), 2) + pow(point[2].getY()-point[3].getY(), 2) + pow(point[2].getZ()-point[3].getZ(), 2));
    //top edge: p0p3
    double p0p3 = sqrt(pow(point[0].getX()-point[3].getX(), 2) + pow(point[0].getY()-point[3].getY(), 2) + pow(point[0].getZ()-point[3].getZ(), 2));

    //calculate triangle area using Heron's formula
    double s1 = (diagonalp0p2 + p0p1 + p1p2)/2; //half of the sum of three edges.
    double s2 = (diagonalp0p2 + p0p3 + p2p3)/2;
    double s3 = (diagonalp1p3 + p0p1 + p0p3)/2;
    double s4 = (diagonalp1p3 + p1p2 + p2p3)/2;

    //two triangle's area using the first method, diagonalp0p2
    double size1 = sqrt(fabs(s1*(s1-p0p1)*(s1-p1p2)*(s1-diagonalp0p2)));
    double size2 = sqrt(fabs(s2*(s2-p0p3)*(s2-p2p3)*(s2-diagonalp0p2)));

    //two triangle's area using the second method, diagonalp1p3
    double size3 = sqrt(fabs(s3*(s3-p0p1)*(s3-p0p3)*(s3-diagonalp1p3)));
    double size4 = sqrt(fabs(s4*(s4-p1p2)*(s4-p2p3)*(s4-diagonalp1p3)));

    //the smaller one of the two method as the area for the sub-quad
    return min(size1+size2, size3+size4);
}



void alignsrep::drawBoundary(){
    vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
    renderer->SetBackground(1,1,1);
    vtkSmartPointer<vtkRenderWindow> renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
    renderWindow->AddRenderer(renderer);
    vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
    renderWindowInteractor->SetRenderWindow(renderWindow);

    int quadNums = (this->quadFig->getRowCount()-1)*(this->quadFig->getColumnCount()-1);//16    
    bool moved = true;
    double quadColor[3] = {1,0.5,0};


    vtkSmartPointer< vtkPoints > hubPosition_b_u = vtkSmartPointer< vtkPoints >::New();//up boundary points.
    vtkSmartPointer< vtkPoints > hubPosition_s = vtkSmartPointer< vtkPoints >::New();//skeletal points.
    //vtkSmartPointer< vtkPoints > hubPosition_b_d = vtkSmartPointer< vtkPoints >::New();//down boundary points.

    visualization visualObject_u(this->quadFig, this->interpolationLevel, 0, moved, quadColor, renderer,0,0,0,0.0);
    //visualization visualObject_d(this->quadFig, this->interpolationLevel, 1, moved, quadColor, renderer);//down boundary

    // Up boundary points.
    hubPosition_b_u = visualObject_u.getSubQuadsPosition(0);//0 is boundary
    hubPosition_s = visualObject_u.getSubQuadsPosition(1);//1 is skeletal

    // Down boundary points.
    //hubPosition_b_d = visualObject_d.getSubQuadsPosition(0);

    for(int i =0;i<quadNums;i++){
        visualObject_u.setColor(quadColor);
        visualObject_u.drawSpoke(hubPosition_b_u,hubPosition_s,i);
        visualObject_u.drawFrameOfQuadByIndex(hubPosition_b_u,i);

        //visualObject_d.drawSpoke(hubPosition_b_d,hubPosition_s,i);
        //visualObject_d.drawFrameOfQuadByIndex(hubPosition_b_d,i);
    }

    renderWindow->Render();
    renderWindowInteractor->Start();
}






/* Draw the crest spokes.*/
void alignsrep::drawCrestSpokes(int interpolationlevel){

    vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
    renderer->SetBackground(1,1,1);
    vtkSmartPointer<vtkRenderWindow> renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
    renderWindow->AddRenderer(renderer);
    vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
    renderWindowInteractor->SetRenderWindow(renderWindow);


    vtkSmartPointer<vtkSRepInterpolateCrestSpokesQuartic> interpolatecrestspokes = vtkSmartPointer<vtkSRepInterpolateCrestSpokesQuartic>::New();
    vtkSmartPointer<vtkPoints> pointscrestspokes  = vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkCellArray> cellarraycrestspokes = vtkSmartPointer<vtkCellArray>::New();
    vtkSmartPointer<vtkPolyData> polycrestspokes = vtkSmartPointer<vtkPolyData>::New();

    interpolatecrestspokes->SetInput(this->srepfig);
    interpolatecrestspokes->SetInterpolationLevel(interpolationlevel);
    interpolatecrestspokes->Update();


    // Method 1
    vtkSmartPointer<vtkSRep> srepcrest = interpolatecrestspokes->GetSRepOutput();
    //cout<<"---srepcrest->GetNumberOfPoints() is: "<<srepcrest->GetNumberOfPoints()<<endl; //Why this equal (crestQuadNum+1) not crestQuadNum?????
    cout<<"crestQuadNum is: "<<crestQuadNum<<endl;

    // From cellid 0 to 19
    for(unsigned i = 0; i < this->crestQuadNum-1; i++){

        vtkSRep::VNLType p0 = interpolatecrestspokes->GetInterpolatedPoint(i,0);
        vtkSRep::VNLType p0s = interpolatecrestspokes->GetInterpolatedSpoke(i,0,0.5);


        vtkSRep::VNLType p1 = interpolatecrestspokes->GetInterpolatedPoint(i+1,0);
        vtkSRep::VNLType p1s = interpolatecrestspokes->GetInterpolatedSpoke(i+1,0,0.5);

        vtkSmartPointer<vtkLine> crestspokeline = vtkSmartPointer<vtkLine>::New();
        crestspokeline->GetPointIds()->SetId(0, pointscrestspokes->InsertNextPoint(p0[0]+p0s[0], p0[1]+p0s[1], p0[2]+p0s[2]));
        crestspokeline->GetPointIds()->SetId(1, pointscrestspokes->InsertNextPoint(p1[0]+p1s[0], p1[1]+p1s[1], p1[2]+p1s[2]));

        cellarraycrestspokes->InsertNextCell(crestspokeline);

    }

    // For the last cellid (this->crestQuadNum-1), its connect to 0.
    vtkSRep::VNLType p0 = interpolatecrestspokes->GetInterpolatedPoint(this->crestQuadNum-1,0);
    vtkSRep::VNLType p0s = interpolatecrestspokes->GetInterpolatedSpoke(this->crestQuadNum-1,0,0.5);


    vtkSRep::VNLType p1 = interpolatecrestspokes->GetInterpolatedPoint(0,0);
    vtkSRep::VNLType p1s = interpolatecrestspokes->GetInterpolatedSpoke(0,0,0.5);

    vtkSmartPointer<vtkLine> crestspokeline = vtkSmartPointer<vtkLine>::New();
    crestspokeline->GetPointIds()->SetId(0, pointscrestspokes->InsertNextPoint(p0[0]+p0s[0], p0[1]+p0s[1], p0[2]+p0s[2]));
    crestspokeline->GetPointIds()->SetId(1, pointscrestspokes->InsertNextPoint(p1[0]+p1s[0], p1[1]+p1s[1], p1[2]+p1s[2]));

    cellarraycrestspokes->InsertNextCell(crestspokeline);

    polycrestspokes->SetPoints(pointscrestspokes);
    polycrestspokes->SetLines(cellarraycrestspokes);

    vtkSmartPointer<vtkPolyDataMapper> crestspokesmapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    crestspokesmapper->SetInput(polycrestspokes);
    vtkSmartPointer<vtkActor>  spokesactor = vtkActor::New();
    spokesactor->SetMapper(crestspokesmapper);
    spokesactor->GetProperty()->SetLineWidth(1);
    spokesactor->GetProperty()->SetColor(0,1,0);
    renderer->AddActor(spokesactor);


    // Method 2
   /* vtkPolyData* interpolatedcrest  = interpolatecrestspokes->GetOutput();

    vtkSmartPointer<vtkPolyDataMapper> crestspokescurvemapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    crestspokescurvemapper->SetInputConnection(interpolatedcrest->GetProducerPort());
    vtkSmartPointer<vtkActor>  crestspokesactor = vtkActor::New();
    crestspokesactor->SetMapper(crestspokescurvemapper);
    double *color = srepfig->GetColor();
    crestspokesactor->GetProperty()->SetColor(color[0],color[1],color[2]);
    renderer->AddActor(crestspokesactor);*/



    renderWindow->Render();
    renderWindowInteractor->Start();

}




/* Draw the crest spokes.*/
void alignsrep::drawCrestSpokes_method2(){

    vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
    renderer->SetBackground(1,1,1);
    vtkSmartPointer<vtkRenderWindow> renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
    renderWindow->AddRenderer(renderer);
    vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
    renderWindowInteractor->SetRenderWindow(renderWindow);


    vtkSmartPointer<vtkSRepInterpolateCrestSpokesQuartic> interpolatecrestspokes = vtkSmartPointer<vtkSRepInterpolateCrestSpokesQuartic>::New();
    vtkSmartPointer<vtkPoints> pointscrestspokes  = vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkCellArray> cellarraycrestspokes = vtkSmartPointer<vtkCellArray>::New();
    vtkSmartPointer<vtkPolyData> polycrestspokes = vtkSmartPointer<vtkPolyData>::New();

    interpolatecrestspokes->SetInput(this->srepfig);
    interpolatecrestspokes->SetInterpolationLevel(this->interpolationLevel);
    interpolatecrestspokes->Update();

    vector<double> subpoints_i_up; //store the v coordinate of all the sub-point along spoke tip connection curve.
    vector<double> subpoints_i_down; //store the v coordinate of all the sub-point along spoke tip connection curve.
    vector<double> subpoints_a; //store the u coordinate of all the sub-points along crest.
    Vector3D points[4];

    // Along the curve between two correspondence up and down spokes.(v direction in crest interpolate method)
    // Up spoke tip is consider as 0, down spoke 1. The medial crest is 0.5.
    // For up crest
    this->tls.splitLine(0,0.5,this->step,subpoints_i_up);
    // For down crest
    this->tls.splitLine(0.5,1,this->step,subpoints_i_down);

    // Along the crest curve in counter clockwise. the start from (0,0).
    this->tls.splitLine(0,1,this->step,subpoints_a);

    // Loop each crest quad.
    for(unsigned i =0; i< this->crestQuadNum; i++){
        for(unsigned m =0; m< subpoints_a.size(); m++){
            // Get the interpolated boundary position at this point
            vtkSRep::VNLType p0 = interpolatecrestspokes->GetInterpolatedPoint(i,subpoints_a[m]);

            // For up crest.
            for(unsigned n=0; n<subpoints_i_up.size(); n++){
                // Get each interpolated boundary spoke direction at this point.
                vtkSRep::VNLType p1 = interpolatecrestspokes->GetInterpolatedSpoke(i,subpoints_a[m],subpoints_i_up[n]);

                vtkSmartPointer<vtkLine> crestspokeline = vtkSmartPointer<vtkLine>::New();
                crestspokeline->GetPointIds()->SetId(0, pointscrestspokes->InsertNextPoint(p0[0], p0[1], p0[2]));
                crestspokeline->GetPointIds()->SetId(1, pointscrestspokes->InsertNextPoint(p0[0]+p1[0], p0[1]+p1[1], p0[2]+p1[2]));

                cellarraycrestspokes->InsertNextCell(crestspokeline);
            }

            // For down crest.
            for(unsigned n=0; n<subpoints_i_down.size(); n++){
                // Get each interpolated boundary spoke direction at this point.
                vtkSRep::VNLType p1 = interpolatecrestspokes->GetInterpolatedSpoke(i,subpoints_a[m],subpoints_i_down[n]);

                vtkSmartPointer<vtkLine> crestspokeline = vtkSmartPointer<vtkLine>::New();
                crestspokeline->GetPointIds()->SetId(0, pointscrestspokes->InsertNextPoint(p0[0], p0[1], p0[2]));
                crestspokeline->GetPointIds()->SetId(1, pointscrestspokes->InsertNextPoint(p0[0]+p1[0], p0[1]+p1[1], p0[2]+p1[2]));

                cellarraycrestspokes->InsertNextCell(crestspokeline);
            }
        }
    }

    polycrestspokes->SetPoints(pointscrestspokes);
    polycrestspokes->SetLines(cellarraycrestspokes);

    vtkSmartPointer<vtkPolyDataMapper> crestspokesmapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    crestspokesmapper->SetInput(polycrestspokes);
    vtkSmartPointer<vtkActor>  spokesactor = vtkActor::New();
    spokesactor->SetMapper(crestspokesmapper);
    spokesactor->GetProperty()->SetLineWidth(1);
    spokesactor->GetProperty()->SetColor(0,1,0);
    renderer->AddActor(spokesactor);


    renderWindow->Render();
    renderWindowInteractor->Start();

}


/* Draw the crest spokes. This is the old method, unit spoke???*/
/*void alignsrep::drawCrestSpokes_old(int interpolationlevel){

    vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
    renderer->SetBackground(1,1,1);
    vtkSmartPointer<vtkRenderWindow> renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
    renderWindow->AddRenderer(renderer);
    vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
    renderWindowInteractor->SetRenderWindow(renderWindow);


    vtkSmartPointer<vtkSRepInterpolateCrestSpokesQuartic> interpolatecrestspokes = vtkSmartPointer<vtkSRepInterpolateCrestSpokesQuartic>::New();
    vtkSmartPointer<vtkPoints> pointscrestspokes  = vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkCellArray> cellarraycrestspokes = vtkSmartPointer<vtkCellArray>::New();
    vtkSmartPointer<vtkPolyData> polycrestspokes = vtkSmartPointer<vtkPolyData>::New();

    interpolatecrestspokes->SetInput(this->srepfig);
    interpolatecrestspokes->SetInterpolationLevel(interpolationlevel);
    interpolatecrestspokes->Update();


    vtkSmartPointer<vtkSRep> srepcrest = interpolatecrestspokes->GetSRepOutput();
    int count =0;
    for(unsigned i = 0; i < srepcrest->GetNumberOfPoints(); i++){
        double point[3];

        srepcrest->GetPoint(i, point);

        vtkIdType id0 = pointscrestspokes->InsertNextPoint(point[0], point[1], point[2]);

        vtkSRep::VectorVNLType currentspokes = srepcrest->GetSpokes(i, false);
        vtkSRep::VectorDoubleType radius = srepcrest->GetSpokesRadius(i);
        for(unsigned j = 0; j < currentspokes.size(); j++){
            cout<<"j is: "<<j<<" , currentspokes[j] is: -------------------- "<<currentspokes[j]<<endl;
            vtkSRep::VNLType p1 = currentspokes[j]*radius[j];

            vtkSmartPointer<vtkLine> crestspokeline = vtkSmartPointer<vtkLine>::New();
            crestspokeline->GetPointIds()->SetId(0, id0);
            crestspokeline->GetPointIds()->SetId(1, pointscrestspokes->InsertNextPoint(point[0] + p1[0], point[1] + p1[1], point[2] + p1[2]));
            cout<<"-----------: "<<point[0] + p1[0]<< point[1] + p1[1]<< point[2] + p1[2]<<endl;
            cout<<"----------------line number is: "<<count<<endl;
            count++;
            cellarraycrestspokes->InsertNextCell(crestspokeline);
        }

        polycrestspokes->SetPoints(pointscrestspokes);
        polycrestspokes->SetLines(cellarraycrestspokes);

        vtkSmartPointer<vtkPolyDataMapper> crestspokesmapper = vtkSmartPointer<vtkPolyDataMapper>::New();
        crestspokesmapper->SetInput(polycrestspokes);
        vtkSmartPointer<vtkActor>  spokesactor = vtkActor::New();
        spokesactor->SetMapper(crestspokesmapper);
        spokesactor->GetProperty()->SetLineWidth(1);
        spokesactor->GetProperty()->SetColor(0,1,0);
        renderer->AddActor(spokesactor);
    }


    vtkPolyData* interpolatedcrest = 0;
    //vtkSmartPointer<vtkSRepInterpolateCrestSpokesQuartic> interpolatecrestspokes = vtkSmartPointer<vtkSRepInterpolateCrestSpokesQuartic>::New();
    //interpolatecrestspokes->SetInput(this->srepfig);
    //interpolatecrestspokes->SetInterpolationLevel(interpolationlevel);
    //interpolatecrestspokes->Update();
    interpolatedcrest = interpolatecrestspokes->GetOutput();

    vtkSmartPointer<vtkPolyDataMapper> crestspokescurvemapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    crestspokescurvemapper->SetInputConnection(interpolatedcrest->GetProducerPort());
    vtkSmartPointer<vtkActor>  crestspokesactor = vtkActor::New();
    crestspokesactor->SetMapper(crestspokescurvemapper);
    double *color = srepfig->GetColor();
    crestspokesactor->GetProperty()->SetColor(color[0],color[1],color[2]);
    renderer->AddActor(crestspokesactor);



    renderWindow->Render();
    renderWindowInteractor->Start();

}
*/

/* Draw the skeletal sheet crest curve.*/
void alignsrep::drawSkeletalCrestCurve(int interpolationlevel){
    vtkSmartPointer<vtkSRepInterpolateMedialCrestCurve> curveinterpolation = vtkSmartPointer<vtkSRepInterpolateMedialCrestCurve>::New();
    curveinterpolation->SetInput(this->srepfig);
    curveinterpolation->SetInterpolationLevel(interpolationlevel);
    curveinterpolation->Update();
    vtkPolyData* medialcrestcurve = curveinterpolation->GetOutput();

    vtkSmartPointer<vtkPolyDataMapper> medialcrestcurvemapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    medialcrestcurvemapper->SetInputConnection(medialcrestcurve->GetProducerPort());
    vtkSmartPointer<vtkActor>  medialsheetcrestcurveactor = vtkActor::New();
    medialsheetcrestcurveactor->SetMapper(medialcrestcurvemapper);
    medialsheetcrestcurveactor->GetProperty()->SetLineWidth(1);
    medialsheetcrestcurveactor->GetProperty()->SetColor(1,1,0);
    renderer->AddActor(medialsheetcrestcurveactor);
}



/* Get all the sub-quad position on the crest by interpolate to t v.
 * Here for the crest interpolate, the t is along the crest, v is along the curve between tip of up spoke and tip of down spoke.
 * GetInterpolatedSpoke(cellid,t,0) is the crest atoms' up spoke's tip on the boundary.
 * GetInterpolatedSpoke(cellid,t,0.5) is the crest spoke's tip on the boundary.
 * GetInterpolatedSpoke(cellid,t,1) is the crest atoms' down spoke's tip on the boundary.
 * For each position GetInterpolatedPoint(cellid,t) there are many spokes(depends on the interpolate level), each spoke's position on
 * the boundary can be get using GetInterpolatedPoint(cellid,t) + GetInterpolatedSpoke(cellid,t,v).
 * Cellid from (0,0), counter clockwise of the crest atoms. That is: (0,0),(0,1),(0,2),...
 * Four each cellid, the t alway changing from 0 to 1.
 * The v changing [0,0.5] is the curve between up spoke's tip to crest spoke's tip, named as up crest.
 * The v changing [0.5,1] is the curve between crest spoke's tip to down spoke's tip, named as up crest.
*/
void alignsrep::getCrestQuadsArea(){

    int crestQuadNum = (this->rowNums-1)*2 + (this->colNums-1)*2;
    if(this->up_b_crest.size() != crestQuadNum){
        cout<<"Message from:alignsrep::getAllSubQuadsVertexesOnCrest: something wrong with the up_b_crest!"<<endl;
        cout<<"this->up_b_crest.size() is: "<<this->up_b_crest.size()<<" while crestQuadNum is: "<<crestQuadNum<<" they should be same!"<<endl;
        return;
    }

    vector<double> subpoints_i_up; //store the v coordinate of all the sub-point along spoke tip connection curve.
    vector<double> subpoints_i_down; //store the v coordinate of all the sub-point along spoke tip connection curve.
    vector<double> subpoints_a; //store the u coordinate of all the sub-points along crest.
    Vector3D points[4];

    // Along the curve between two correspondence up and down spokes.(v direction in crest interpolate method)
    // Up spoke tip is consider as 0, down spoke 1. The medial crest is 0.5.
    // For up crest
    this->tls.splitLine(0,0.5,this->step,subpoints_i_up);
    // For down crest
    this->tls.splitLine(0.5,1,this->step,subpoints_i_down);

    // Along the crest curve in counter clockwise. the start from (0,0).
    this->tls.splitLine(0,1,this->step,subpoints_a);

    // New a crest interplator....
    vtkSmartPointer<vtkSRepInterpolateCrestSpokesQuartic> interpolatecrestspokes = vtkSmartPointer<vtkSRepInterpolateCrestSpokesQuartic>::New();
    interpolatecrestspokes->SetInput(this->srepfig);
    interpolatecrestspokes->SetInterpolationLevel(this->interpolationLevel);
    interpolatecrestspokes->Update();

    // Loop each crest quad.
    for(unsigned i =0; i< crestQuadNum; i++){
        this->up_b_crest_areas.push_back(0);
        this->down_b_crest_areas.push_back(0);
        double subquadarea_up =0;
        double subquadarea_down =0;

        for(unsigned m =0; m< subpoints_a.size()-1; m++){
            // Get the interpolated boundary position at this point
            vtkSRep::VNLType bp_c = interpolatecrestspokes->GetInterpolatedPoint(i,subpoints_a[m]);

            // Get the interpolated boundary position at this point
            vtkSRep::VNLType bp_c_next = interpolatecrestspokes->GetInterpolatedPoint(i,subpoints_a[m+1]);

            // For up crest.            
            for(unsigned n=0; n<subpoints_i_up.size()-1; n++){
                // Get each interpolated boundary spoke direction at this point.
                vtkSRep::VNLType p0 = interpolatecrestspokes->GetInterpolatedSpoke(i,subpoints_a[m],subpoints_i_up[n]);
                points[0].setX(bp_c[0] + p0[0]);
                points[0].setY(bp_c[1] + p0[1]);
                points[0].setZ(bp_c[2] + p0[2]);
                vtkSRep::VNLType p1 = interpolatecrestspokes->GetInterpolatedSpoke(i,subpoints_a[m],subpoints_i_up[n+1]);
                points[1].setX(bp_c[0] + p1[0]);
                points[1].setY(bp_c[1] + p1[1]);
                points[1].setZ(bp_c[2] + p1[2]);
                vtkSRep::VNLType p2 = interpolatecrestspokes->GetInterpolatedSpoke(i,subpoints_a[m+1],subpoints_i_up[n+1]);
                points[2].setX(bp_c_next[0] + p2[0]);
                points[2].setY(bp_c_next[1] + p2[1]);
                points[2].setZ(bp_c_next[2] + p2[2]);
                vtkSRep::VNLType p3 = interpolatecrestspokes->GetInterpolatedSpoke(i,subpoints_a[m+1],subpoints_i_up[n]);
                points[3].setX(bp_c_next[0] + p3[0]);
                points[3].setY(bp_c_next[1] + p3[1]);
                points[3].setZ(bp_c_next[2] + p3[2]);
                //cout<<points[0]<<points[1]<<points[2]<<points[3]<<endl;

                // Compute the sub quad area.
                subquadarea_up += quadArea(points);
            }

            // For down crest.            
            for(unsigned n=0; n<subpoints_i_down.size()-1; n++){
                // Get each interpolated boundary spoke direction at this point.
                vtkSRep::VNLType p0 = interpolatecrestspokes->GetInterpolatedSpoke(i,subpoints_a[m],subpoints_i_down[n]);
                points[0].setX(bp_c[0] + p0[0]);
                points[0].setY(bp_c[1] + p0[1]);
                points[0].setZ(bp_c[2] + p0[2]);
                vtkSRep::VNLType p1 = interpolatecrestspokes->GetInterpolatedSpoke(i,subpoints_a[m],subpoints_i_down[n+1]);
                points[1].setX(bp_c[0] + p1[0]);
                points[1].setY(bp_c[1] + p1[1]);
                points[1].setZ(bp_c[2] + p1[2]);
                vtkSRep::VNLType p2 = interpolatecrestspokes->GetInterpolatedSpoke(i,subpoints_a[m+1],subpoints_i_down[n+1]);
                points[2].setX(bp_c_next[0] + p2[0]);
                points[2].setY(bp_c_next[1] + p2[1]);
                points[2].setZ(bp_c_next[2] + p2[2]);
                vtkSRep::VNLType p3 = interpolatecrestspokes->GetInterpolatedSpoke(i,subpoints_a[m+1],subpoints_i_down[n]);
                points[3].setX(bp_c_next[0] + p3[0]);
                points[3].setY(bp_c_next[1] + p3[1]);
                points[3].setZ(bp_c_next[2] + p3[2]);
                //cout<<points[0]<<points[1]<<points[2]<<points[3]<<endl;

                // Compute the sub quad area.
                subquadarea_down += quadArea(points);
            }
        }

        // Compute each quad's area.
        this->up_b_crest_areas[i] = subquadarea_up;
        //cout<<this->up_b_crest_areas[i]<<endl;
        this->down_b_crest_areas[i] = subquadarea_down;
        //cout<<this->down_b_crest_areas[i]<<endl;

    }

    /*vtkSRep::VNLType p0 = interpolatecrestspokes->GetInterpolatedSpoke(0,0,0.5);
    cout<<"up crest point:["<<0<<","<<0<<"] is: "<<p0[0]<<p0[1]<<p0[2]<<endl;
    vtkSRep::VNLType pp = interpolatecrestspokes->GetInterpolatedPoint(0,0);
    cout<<"media crest boundary points:["<<0<<","<<0<<"] is: "<<pp[0]+p0[0]<<pp[1]+p0[1]<<pp[2]+p0[2]<<endl;

    // GetInterpolatedPoint(0,0)+GetInterpolatedSpoke(0,0,0) equals to getX()+getU0*getR0;
    // GetInterpolatedPoint(0,0)+GetInterpolatedSpoke(0,0,0.5) equals to getX()+getUEnd*getREnd;
    // GetInterpolatedPoint(0,0)+GetInterpolatedSpoke(0,0,1) equals to getX()+getU1*getR1;

    */
}




/* This method only compute the crest quads' area under interpolation level 0.
 * In fact, this method do not use interpolation method, it use the primitives. Just should be the same as do it use getCrestQuadsArea()
 * under interpolation level 0.
*/
void alignsrep::getCrestQuadsArea_InteLevel0(){

    int crestQuadNum = this->crestAtomNums;

    if(this->up_b_crest.size() != crestQuadNum){
        cout<<"Message from:alignsrep::getAllSubQuadsVertexesOnCrest: something wrong with the up_b_crest!"<<endl;
        cout<<"this->up_b_crest.size() is: "<<this->up_b_crest.size()<<" while crestQuadNum is: "<<crestQuadNum<<" they should be same!"<<endl;
        return;
    }

    Vector3D points[4];

    // Loop each crest quad.
    for(unsigned i =0; i< crestQuadNum-1; i++){
        this->up_b_crest_areas.push_back(0);
        this->down_b_crest_areas.push_back(0);

        // For up crest.
        points[0] = this->up_b_crest[i]; //p0
        points[1] = this->medial_crest[i]; //p1
        points[2] = this->medial_crest[i+1]; //p2
        points[3] = this->up_b_crest[i+1]; //p3
        //cout<<"up crest: "<<points[0]<<points[1]<<points[2]<<points[3]<<endl;

        // Compute each quad's area.
        this->up_b_crest_areas[i] = quadArea(points);
        //cout<<"--------this->up_b_crest_areas["<<i<<"] is: "<<this->up_b_crest_areas[i]<<endl;

        // For down crest.
        points[0] = this->medial_crest[i]; //p0
        points[1] = this->down_b_crest[i]; //p1
        points[2] = this->down_b_crest[i+1]; //p2
        points[3] = this->medial_crest[i+1]; //p3
        cout<<"down crest: "<<points[0]<<points[1]<<points[2]<<points[3]<<endl;

        // Compute each quad's area.
        this->down_b_crest_areas[i] = quadArea(points);
        //cout<<"--------this->down_b_crest_areas["<<i<<"] is: "<<this->down_b_crest_areas[i]<<endl;
    }

    int i = crestQuadNum-1;
    // The last quad is special, it compozite of the last point and the first point.
    // For the last quad of up crest.
    points[0] = this->up_b_crest[i]; //p0
    points[1] = this->medial_crest[i]; //p1
    points[2] = this->medial_crest[0]; //p2
    points[3] = this->up_b_crest[0]; //p3

    // Compute each quad's area.
    this->up_b_crest_areas[i] = quadArea(points);
    cout<<"--------this->up_b_crest_areas["<<i<<"] is: "<<this->up_b_crest_areas[i]<<endl;

    // For the last quad of down crest.
    points[0] = this->medial_crest[i]; //p0
    points[1] = this->down_b_crest[i]; //p1
    points[2] = this->down_b_crest[0]; //p2
    points[3] = this->medial_crest[0]; //p3

    // Compute each quad's area.
    this->down_b_crest_areas[i] = quadArea(points);
    cout<<"--------this->down_b_crest_areas["<<i<<"] is: "<<this->down_b_crest_areas[i]<<endl;
}


/* Return the boundary ponints vector.*/
vector<Vector3D> alignsrep::getBoundaryPoints(){
    return this->boundaryPoints;
}

/* Return sum of the boundary points' corresponding 4 or 3 areas.
 * wAreas is a points number row by 1 coulum matrix.
*/
vector<double> alignsrep::getBoundaryPointsAreas(){
    vector<double> wAreas;
    for(unsigned i =0; i<this->boundaryPointsAreas.size(); i++){
        double sumarea = 0;
        for(unsigned j=0; j<this->boundaryPointsAreas[0].size(); j++){
            sumarea += this->boundaryPointsAreas[i][j];
        }
        wAreas.push_back(sumarea);
    }

    return wAreas;
}


/* Get the skeletal points on this quadFig.
*/
vector<Vector3D> alignsrep::getPointsOnSkeletal(){

    M3DQuadPrimitive* prim;

    vector<Vector3D> sPoints;

    //Loop each primitive, get the x,y,z coordinate of this srep's boundary.
    for(unsigned u =0; u<this->quadFig->getRowCount(); u++){
        for(unsigned v =0; v<this->quadFig->getColumnCount(); v++){
            prim = dynamic_cast<M3DQuadPrimitive*>(this->quadFig->getPrimitivePtr(u,v));

            sPoints.push_back(prim->getX());
        }
    }

    /*for (unsigned i=0; i<sPoints.size();i++){
        cout<<"-------------------"<<sPoints[i]<<endl;
    }*/

    return sPoints;
}




/* Get the area of up or down boundary crest.*/
vector<double> alignsrep::getCrestAreas(int side){

    if(side==0){
        return this->up_b_crest_areas;
    }
    else
        return this->down_b_crest_areas;
}


vector<double> alignsrep::getQuadsAreas(int side){
    if(side==0){
        return this->up_quads_areas;
    }
    else
        return this->down_quads_areas;
}



/* For compare only..
 * Set all the wAreas to 1(as no weight).
 * Return a points number row by 1 coulum matrix.
*/
vector<double> alignsrep::getOnesWeight(int pointsNum){
    cout<<"----------------- alignesrep pointsNum is: "<<pointsNum<<endl;
    vector<double> wAreas;
    for(unsigned i =0; i< pointsNum; i++){   // point size
        wAreas.push_back(1);
    }

    return wAreas;
}



int alignsrep::getRows(){
    return this->rowNums;
}
int alignsrep::getCols(){
    return this->colNums;
}
int alignsrep::getAtomNums(){
    return this->totalAtomNums;
}
int alignsrep::getInteriorAtomNums(){
    return this->interiorAtomNums;
}
int alignsrep::getCrestAtomNums(){
    return this->crestAtomNums;
}





