#include "uvmap.h"

uvmap::uvmap()
{
}



/* Return subdivided u & v coordinate based on the interpolation level.
 * Do the same thing as function: calregularityfeatures::getSubQuadsUVCoordinateNotUsingDeltaUV()
*/
void uvmap::getSubUVCoordinateForStandSide(DoubleVecVec &subUs, DoubleVecVec &subVs, int rowNum, int colNum, int step){

    int quadIndex = 0;

    for(unsigned i = 0; i < rowNum -1; i++){ //the row number of the quads. its 3.
        for(unsigned j = 0; j < colNum -1; j++){//coloums, its 13.
                subUs.push_back(DoubleVec());
                subVs.push_back(DoubleVec());

                //four points of a quad with delta u, v. counter clockwise.
                //double quadpointsu[4], quadpointsv[4];
                DoubleVec quadpointsu, quadpointsv;
                quadpointsv.resize(4);
                quadpointsu.resize(4);
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
                DoubleVec subquadpoint_u = splitQuad(quadpointsu[0],quadpointsu[1],quadpointsu[2],quadpointsu[3], step);

                //second, get the 25 subquads's v coordinate, in column first order.
                DoubleVec subquadpoint_v = splitQuad(quadpointsv[0],quadpointsv[1],quadpointsv[2],quadpointsv[3], step);

                for(unsigned m =0; m<subquadpoint_u.size();m++){
                    subUs[quadIndex].push_back(subquadpoint_u[m]);
                    subVs[quadIndex].push_back(subquadpoint_v[m]);
               }

                quadIndex++;
        }
    }
}




/* Return subdivided u & v coordinate based on the interpolation level.
 * Do the same thing as function: calregularityfeatures::getSubQuadsUVCoordinateNotUsingDeltaUV()
*/
DoubleVecVec uvmap::getSubUVCoordinateForCrestSide(int step){
    DoubleVecVec subUVs;
    subUVs.push_back(DoubleVec());

    // Along the curve between two correspondence up and down spokes.(v direction in crest interpolate method)
    splitLine(0, 1, step, subUVs[0]); // Up spoke tip is consider as 0, down spoke 1. The medial crest is 0.5.

    return subUVs;
}


DoubleVec uvmap::getUVSubdivision(int interpolationLevel){
    int step = pow(2.0, (double)interpolationLevel);

    DoubleVec subUVs;

    // Along the curve between two correspondence up and down spokes.(v direction in crest interpolate method)
    splitLine(0, 1, step, subUVs); // Up spoke tip is consider as 0, down spoke 1. The medial crest is 0.5.

    return subUVs;
}


/* Split a quad into many sub quads. By bilinear interpolating to u and v seperately for each quad point.
 * Return a double vector that store all the u or v coordinate of the sub-quad.
 * p0, p1, p2, p3 is u or v coordinate of the four vertex of the quad.
 * first consider only the u coordinate, then do the same to v coordinate.
 * The bilinear interpolating works this way:
 * Firstly, divide the row line p0p3 and p1p2 each into n(n equals to 2^interpolationLevel) pieces;
 * Secondly, divide the n+1 column line into n pieces.
 * It return all the points contained in the quad in a column first order.
 */
DoubleVec uvmap::splitQuad(double p0, double p1, double p2, double p3, int step){

    //do bilinear interpolate to the two rows.
    //for the first row, p0p3
    DoubleVec firstrow, secondrow;
    splitLine(p0, p3, step, firstrow);

    //for the second row, p1p2
    splitLine(p1, p2, step, secondrow);

    //do bilinear interpolate to the step+1 coloum line, each line is divided into 2^interpolationLevel pieces.
    DoubleVec columntemp;
    DoubleVec quadpoints;  //store the final subquads vertexes.

    //loop each pair of points in firstrow[i] and secondrow[i].
    for(unsigned i =0; i<firstrow.size(); i++){
        //interpolate to each pair of column
        splitLine(firstrow[i], secondrow[i], step, columntemp);

        //store the column's step+1 points into vector quadpoints.
        for(unsigned j=0; j<columntemp.size(); j++){
            //cout<<"columntemp ["<<j<<"] is: "<<columntemp[j]<<endl;
            quadpoints.push_back(columntemp[j]);
        }

        //clear the vector for next loop.
        columntemp.clear();
    }

    firstrow.clear();
    secondrow.clear();

    return quadpoints;
}



/* Bilinear interpolate between two points:
 * the line between point1 and point2 is divided into 2^interpolationLevel pieces.
 * step: the pieces a line will be divided into;
 * the resultpoints vector size will be step+1, because step sub-line has step+1 points.
*/
void uvmap::splitLine(double point1, double point2, int step, DoubleVec &resultpoints){
    // Clear vector
    resultpoints.clear();

    double stepdistance = (point2 - point1)/step;

    //resultpoints is an array, contain step+1 double value (step+1 points).
    for(unsigned i =0; i<step+1; i++){
        resultpoints.push_back(point1 + i*stepdistance);
    }
}




