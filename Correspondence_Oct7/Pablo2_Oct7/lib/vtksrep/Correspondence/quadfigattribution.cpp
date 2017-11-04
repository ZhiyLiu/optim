#include "quadfigattribution.h"

quadfigattribution::quadfigattribution()
{
}

//quadfigattribution::quadfigattribution(M3DQuadFigure *quadFig, int interpolationLevel){
//    this->rowNum =

//    this->exteriorAtomNum = rowNum*2 + colNum*2 - 4; // 28
//    this->interiorAtomNum = rowNum*colNum - this->exteriorAtomNum; // 11

//    this->step = pow((double)2, (double)interpolationLevel);
//    this->subQuadNum = pow((double)2, (double)step);
//    this->quadNum = (rowNum-1)*(colNum-1);
//    this->subQuadVertexNum = (step+1)*(step+1);
//}

quadfigattribution::quadfigattribution(M3DQuadFigure *quadFig, int interpolationLevel){
    this->quadFig = quadFig;
    this->interpolationLevel = interpolationLevel;
}




int quadfigattribution::getRowNum(){
    return this->quadFig->getRowCount();
}

int quadfigattribution::getColNum(){
    return this->quadFig->getColumnCount();
}

int quadfigattribution::getExteriorAtomNum(){
    return (getRowNum()+getColNum())*2 - 4;
}

int quadfigattribution::getInteriorAtomNum(){
    return getRowNum()*getColNum() - getExteriorAtomNum();
}

int quadfigattribution::getStep(){
    return pow(2.0, (double)this->interpolationLevel);
}

int quadfigattribution::getSubQuadNum(){
    return pow((double)getStep(), 2.0);
}

int quadfigattribution::getQuadNum(){
    return (getRowNum()-1) * (getColNum()-1);
}

int quadfigattribution::getSubQuadVertexNum(){
    int temp = getStep() + 1;

    return temp*temp;
}

M3DQuadFigure * quadfigattribution::getQuadFig(){
    return this->quadFig;
}

int quadfigattribution::getInterpolationLevel(){
    return this->interpolationLevel;
}


void quadfigattribution::setFigAttribution(M3DQuadFigure *quadFig, int interpolationLevel){
    this->quadFig = quadFig;
    this->interpolationLevel = interpolationLevel;
}


/* return the horizontal edge number of the standard side*/
int quadfigattribution::getHorEdgeNumOfStandSide(){
    return (getRowNum()-1) * getColNum(); //26
}

int quadfigattribution::getVerEdgeNumOfStandSide(){
    return (getColNum()-1) * getRowNum(); //36
}


/* return the index of the bottom-right sub-point (p1, in horizontal direction) */
int quadfigattribution::getBotRigHorPosIndex(){
    return (getStep()+1)*getStep() - 1;
}

int quadfigattribution::getBotRigVerPosIndex(){
    return getSubQuadVertexNum()-2;
}



















