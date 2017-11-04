#include <iostream>
#include "M3DQuadInterpolater.h"   //Jared's interpolation method to up and down spoke.
#include "calregularityfeatures.h"
#include "toolsfunc.h"
#include <math.h>
#include "visualization.h"



using namespace std;

calregularityfeatures::calregularityfeatures()
{
}


calregularityfeatures::calregularityfeatures(M3DQuadFigure *quadFig, int interpolationLevel, int side, DoubleVec subUVs) {

    this->side = side;

    this->subUVs = subUVs;

    this->curQuadFig = quadFig;
    this->interpolationLevel = interpolationLevel;

    this->quadNum = (quadFig->getRowCount()-1)*(quadFig->getColumnCount()-1);

    this->step = pow(2.0, (double)interpolationLevel);

    this->subQuadVetexNum = pow((step+1.0), 2.0); //(step+1)*(step+1);


}



calregularityfeatures::~calregularityfeatures()
{
}






/* Get all the sub-quad position by given u v coordinates.
 * quadtype: boundary quads(0), skeletal quads(1).
*/
vtkSmartPointer< vtkPoints > calregularityfeatures::getSubQuadsPosition(int quadtype){

    DoubleVecVec subUs;
    DoubleVecVec subVs;

    uvmap mp;
    mp.getSubUVCoordinateForStandSide(subUs, subVs, this->curQuadFig->getRowCount(), this->curQuadFig->getColumnCount(), this->step);

    vtkSmartPointer< vtkPoints > hubpos = vtkSmartPointer< vtkPoints >::New();    

    M3DQuadInterpolater standsideintp;

    //store the all the position of points for a quad, each point is a 3D vector, which store the x, y, z coordinate of this point.
    Vector3D point;

    for(int i = 0; i < this->quadNum; i++){
        //cout<<"Currently drawing quad: "<<i<<endl;
        for(int j =0; j < this->subQuadVetexNum;j++){
            M3DSpoke newspoke = standsideintp.interpolateQuadSpoke2(this->curQuadFig, subUs[i][j], subVs[i][j], side);
            std::cout<<"----subUs["<<i<<"]["<<j<<"] is: "<<subUs[i][j]<<"    subVs["<<i<<"]["<<j<<"] is: "<<subVs[i][j]<<std::endl;
            //quadtype: Two kinds of quads, boundary quads(0), skeletal quads(1)
            if(quadtype==0){ //boundary quads(0)
                point = newspoke.getB();
                //cout<<"point position on boundary: x is: "<<point.getX()<<" y is: "<<point.getY()<<" z is: "<<point.getZ()<<endl;
            }
            else { //skeletal quads(1)
                point = newspoke.getX();
            }

            hubpos->InsertNextPoint(point.getX(),point.getY(),point.getZ());
        }
    }

    return hubpos;
}



//calculate the distance between two points.
double calregularityfeatures::lengthofedges(Vector3D point[2]){
    double x0 = point[0].getX();
    double y0 = point[0].getY();
    double z0 = point[0].getZ();
    //cout<<"point[0] is:"<<x0<<"; "<<y0<<"; "<<z0<<endl;

    double x1 = point[1].getX();
    double y1 = point[1].getY();
    double z1 = point[1].getZ();
    //cout<<"point[1] is:"<<x1<<"; "<<y1<<"; "<<z1<<endl;

    double distance = sqrt(pow(x1-x0, 2) + pow(y1-y0, 2) + pow(z1-z0, 2));

    //cout<<"distance between two edges is: "<<distance<<endl;

    return distance;
}

////calculate the distance between two points.
//double calregularityfeatures::pointsDistance(double * p1, double * p2){

//    double distance = sqrt(pow(p2[0]-p1[0], 2) + pow(p2[1]-p1[1], 2) + pow(p2[2]-p1[2], 2));

//    return distance;
//}

//calculate the distance between two points.
double calregularityfeatures::pointsDistance(double * p1, double * p2){

    double distance = 0.0;

    for(unsigned int dim = 0; dim < 3; dim++){
        distance += pow((p2[dim] - p1[dim]), 2.0);
    }

    return sqrt(distance);
}




void calregularityfeatures::calculateEdges_method3(vtkSmartPointer< vtkPoints > hubpos, MatrixType &horEdgeFeatureMatrix,
                                                   MatrixType &verEdgeFeatureMatrix, int quadType){
    int colNum = this->curQuadFig->getColumnCount();

    double p1[3];
    double p2[3];

    int horEdgeNum = 0;
    //For each quad, its left side is horizonal lines.
    for(int m=0; m<this->quadNum; m++){
        double sumOfSubLines = 0;
        //left side. For each sub-line.
        for(int n=0; n<this->step;n++){
            int currIndex = m*this->subQuadVetexNum+n;

            hubpos->GetPoint(currIndex, p1);
            hubpos->GetPoint(currIndex+1, p2);

            sumOfSubLines += pointsDistance(p1, p2);
        }
        horEdgeFeatureMatrix[quadType][horEdgeNum] =sumOfSubLines;
        horEdgeNum++;
    }

    //For the quads (colNum-1)*colNumIndex -1, its right side is also horizonal lines.
    for(int m=0; m<this->quadNum; m++){
        if((m+1)%(colNum-1)==0){
            double sumOfSubLines = 0;
            for(int n=this->subQuadVetexNum-(this->step+1); n<this->subQuadVetexNum-1;n++){
                int currIndex = m*this->subQuadVetexNum+n;

                hubpos->GetPoint(currIndex, p1);
                hubpos->GetPoint(currIndex+1, p2);

                sumOfSubLines += pointsDistance(p1, p2);
            }
            horEdgeFeatureMatrix[quadType][horEdgeNum] =sumOfSubLines;
            horEdgeNum++;
        }
    }

    //Vertical lines.
    int verEdgeNum = 0;
    //For each quad, its top side is vertical lines.
    for(int m=0; m<this->quadNum;m++){
        double sumOfSubLines = 0;
        for(int k = 0; k < this->step; k++){
            int currIndex = m*this->subQuadVetexNum + k*(this->step+1);
            int nextIndex = currIndex + (this->step+1);

            hubpos->GetPoint(currIndex, p1);
            hubpos->GetPoint(nextIndex, p2);

           sumOfSubLines += pointsDistance(p1, p2);
        }
        verEdgeFeatureMatrix[quadType][verEdgeNum] = sumOfSubLines;
        verEdgeNum++;
    }

    //For the last colNum-1 quads, its bottom side is also vertical lines.
    for(int m=this->quadNum-(colNum-1); m<this->quadNum;m++){
        double sumOfSubLines = 0;
        for(int k = 0; k<this->step;k++){
            int currIndex = m*this->subQuadVetexNum+this->step + k*(this->step+1);
            int nextIndex = currIndex + (this->step+1);

            hubpos->GetPoint(currIndex, p1);
            hubpos->GetPoint(nextIndex, p2);

           sumOfSubLines += pointsDistance(p1, p2);
        }
        verEdgeFeatureMatrix[quadType][verEdgeNum] = sumOfSubLines;
        verEdgeNum++;
    }
}





/* Compute the angles on each quad.
 * use the avarage angle of the sub-quads.
 */
void calregularityfeatures::calculateAnglesForStandardSide(vtkSmartPointer< vtkPoints > hubpos, MatrixType &angleFeatureMatrix, int quadType) {

    double p0[3];
    double p1[3];
    double p2[3];
    double p3[3];

    toolsfunc tools;

    // For each quad
    for(int n = 0; n < this->quadNum; n++){
        double angle_subQuads_tl = 0.0; // top-left angle
        double angle_subQuads_br = 0.0; // bottom-right angle

        //For each sub-quad.
        for(int m=0; m<this->step;m++){
            for(int k = 0; k<this->step;k++){
                int currentpoint = k + m*(this->step+1);
                // top-left point.
                hubpos->GetPoint(n*this->subQuadVetexNum+currentpoint,p0);  //m
                // bottom-left point
                hubpos->GetPoint(n*this->subQuadVetexNum+currentpoint+1,p1);//m+1
                // bottom-right point
                hubpos->GetPoint(n*this->subQuadVetexNum+currentpoint+1+this->step+1,p2);//m+6
                // top-right point
                hubpos->GetPoint(n*this->subQuadVetexNum+currentpoint+this->step+1,p3);//m+5

                angle_subQuads_tl += tools.dotProductAngle(p0, p3, p1);
                angle_subQuads_br += tools.dotProductAngle(p2, p1, p3);
            }
        }

        // top-left corner trangular
        hubpos->GetPoint(n*this->subQuadVetexNum,p0);
        hubpos->GetPoint(n*this->subQuadVetexNum+(this->step+1),p3);
        hubpos->GetPoint(n*this->subQuadVetexNum+1,p1);
        Vector3D normal_tl = tools.trangularNormal_doublePot(p0, p3, p1);

        // bottom-right corner trangular
        int p1_index = (this->step+1)*this->step - 1;
        hubpos->GetPoint(n*this->subQuadVetexNum+p1_index+this->step+1,p2);
        hubpos->GetPoint(n*this->subQuadVetexNum+p1_index,p1);
        hubpos->GetPoint(n*this->subQuadVetexNum+p1_index+this->step,p3);
        Vector3D normal_br = tools.trangularNormal_doublePot(p2, p1, p3);

        normal_tl.normalize();
        normal_br.normalize();

        int subQuadNum = this->step * this->step;

        if(quadType == 0){ // on boundary
            // up-left angle
            angleFeatureMatrix[0][n] = angle_subQuads_tl/subQuadNum;

            // bot-right corner
            angleFeatureMatrix[1][n] = angle_subQuads_br/subQuadNum;

            // normal swing
            angleFeatureMatrix[2][n] = normal_tl * normal_br;
        }
        else { // on skeletal
            // up-left angle
            angleFeatureMatrix[3][n] = angle_subQuads_tl/subQuadNum;

            // bot-right corner
            angleFeatureMatrix[4][n] = angle_subQuads_br/subQuadNum;

            // normal swing
            angleFeatureMatrix[5][n] = normal_tl * normal_br;
        }
    }
}



/* Get all the sub-quad position by given u v coordinates.
 * quadtype: boundary quads(0), skeletal quads(1).
*/
void calregularityfeatures::getSubQuadVertexesPosition(vtkSmartPointer< vtkPoints > bPoints, vtkSmartPointer< vtkPoints > sPoints){

    M3DQuadInterpolater standsideintp;

    //store the all the position of points for a quad, each point is a 3D vector, which store the x, y, z coordinate of this point.
    Vector3D point;

    int linePosNum = this->subUVs.size(); // points number on each line
    // for each quad
    for(unsigned u = 0; u < curQuadFig->getRowCount() -1; u++){ //the row number of the quads. its 3.
        for(unsigned v = 0; v < curQuadFig->getColumnCount() -1; v++){ //colums, its 13.
            for(unsigned int m = 0; m < linePosNum; m++){ // subdivision in v direction (along the fold curve) ,v
                double subV = v + this->subUVs[m];
                for(unsigned int n = 0; n < linePosNum; n++){ // In h direction (u)
                    double subU = u + this->subUVs[n];

                    M3DSpoke newspoke = standsideintp.interpolateQuadSpoke2(this->curQuadFig, subU, subV, side);
                    //std::cout<<"----subUs["<<u<<"]["<<v<<"] is: "<<subU<<"    subVs["<<u<<"]["<<v<<"] is: "<<subV<<std::endl;

                    //boundary quads(0)
                    point = newspoke.getB();
                    bPoints->InsertNextPoint(point.getX(),point.getY(),point.getZ());

                    //skeletal quads(1)
                    point = newspoke.getX();
                    sPoints->InsertNextPoint(point.getX(),point.getY(),point.getZ());
                }
            }
        }
    }
}
