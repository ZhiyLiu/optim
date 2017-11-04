#ifndef UVMAP_H
#define UVMAP_H

#include <vector>
#include "math.h"

typedef std::vector<double> DoubleVec;
typedef std::vector<DoubleVec> DoubleVecVec;//store the xyz position of the points.


class uvmap
{
public:
    uvmap();

    DoubleVec splitQuad(double p0, double p1, double p2, double p3, int step);
    void splitLine(double point1, double point2, int step, DoubleVec &resultpoints);
    void getSubUVCoordinateForStandSide(DoubleVecVec &subUs, DoubleVecVec &subVs, int rowNum, int colNum, int step);
    DoubleVecVec getSubUVCoordinateForCrestSide(int step);

    DoubleVec getUVSubdivision(int step);




};

#endif // UVMAP_H