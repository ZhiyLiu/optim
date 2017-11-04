#ifndef QUADFIGATTRIBUTION_H
#define QUADFIGATTRIBUTION_H


#include "M3DQuadFigure.h"


class quadfigattribution
{
public:
    quadfigattribution();
    quadfigattribution(M3DQuadFigure *quadFig);
    quadfigattribution(M3DQuadFigure *quadFig, int interpolationLevel);

    int getRowNum();
    int getColNum();
    int getExteriorAtomNum();
    int getInteriorAtomNum();

    int getStep();
    int getSubQuadNum();
    int getQuadNum();
    int getSubQuadVertexNum();
    M3DQuadFigure * getQuadFig();
    int getInterpolationLevel();


    void setFigAttribution(M3DQuadFigure *quadFig, int interpolationLevel);

    int getHorEdgeNumOfStandSide();
    int getVerEdgeNumOfStandSide();

    int getBotRigHorPosIndex();
    int getBotRigVerPosIndex();



private:
    M3DQuadFigure *quadFig;
    int interpolationLevel;
//    int rowNum;
//    int colNum;
//    int exteriorAtomNum;
//    int interiorAtomNum;

//    int step;
//    int subQuadNum;
//    int quadNum;
//    int subQuadVertexNum; //the subdivided sub-quads' vertexes


//    int horEdgeNum = (this->rowNum-1) * this->colNum; //26
//    int verEdgeNum = (this->colNum-1) * this->rowNum; //36

};

#endif // QUADFIGATTRIBUTION_H
