#include <math.h>
#include "M3DFigureElongater.h"

M3DFigureElongater::M3DFigureElongater(M3DQuadFigure * _figure) :
    elongationMethod(M3D_ELONGATE_LONG_AXIS), rowOrientedLinks(false)
{
    referenceFigure = NULL;
    initialize(_figure);
}

M3DFigureElongater::~M3DFigureElongater()
{
    if(referenceFigure != NULL)
        delete referenceFigure;
}

void M3DFigureElongater::initialize(M3DQuadFigure * _figure)
{
    figure = _figure;

	if (referenceFigure != NULL)
		delete referenceFigure;

    if(figure == NULL)
    {
        referenceFigure = NULL;
        return;
    }

    // FIXME: not sure what to do when a M3DTubeFigure is passed - elongation is not defined.
    // perhaps this class should be renamed M3DQuadFigureElongator
    referenceFigure = dynamic_cast<M3DQuadFigure *>( figure->assign() );	// Doesn't copy the name

    setElongationMethod(elongationMethod);
}

void M3DFigureElongater::update()
{
    *referenceFigure = *figure; 	// Doesn't copy the name
}

void M3DFigureElongater::setElongationMethod(M3DElongationMethod m)
{
    int rowCount,
        columnCount;


    elongationMethod = m;

    if(referenceFigure == NULL)
        return;

    rowCount = referenceFigure->getRowCount();
    columnCount = referenceFigure->getColumnCount();

    if((elongationMethod == M3D_ELONGATE_LONG_AXIS && rowCount >= columnCount) ||
       (elongationMethod == M3D_ELONGATE_SHORT_AXIS && rowCount <= columnCount) ||
       (elongationMethod == M3D_ELONGATE_ROW))

       rowOrientedLinks = true;

    else
        rowOrientedLinks = false;
}

void M3DFigureElongater::elongate(double value)
{
    M3DPrimitive * primitive;

    int linkCount,
        rowCount, columnCount,
        i, j,
        primitiveId,
        primitiveIndexDelta,
        neighborIndexDelta,
        numPrimitivesPerLink;

    int prevPrimitiveId,
        neighborPrimitiveId;


    if(figure == NULL || referenceFigure == NULL)
        return;

    (*figure) = (*referenceFigure);

    rowCount = figure->getRowCount();
    columnCount = figure->getColumnCount();

    if(rowOrientedLinks)
    {
        primitiveIndexDelta = columnCount;
        numPrimitivesPerLink = rowCount;
        neighborIndexDelta = 1;
        linkCount = columnCount;
    }
    else
    {
        primitiveIndexDelta = 1;
        numPrimitivesPerLink = columnCount;
        neighborIndexDelta = columnCount;
        linkCount = rowCount;
    }

    int middleIndex = (linkCount / 2) * neighborIndexDelta;
    int deltaIndex;

//    lastElongateValue = value * lastElongateValue;

    int startIndex = middleIndex;
    for(deltaIndex = -neighborIndexDelta; deltaIndex <= neighborIndexDelta;
        deltaIndex += (2*neighborIndexDelta))
    {
        for(i = startIndex; i < linkCount*neighborIndexDelta && i >= 0; i += deltaIndex)
        {
            primitiveId = i;
            primitive = figure->getPrimitivePtr(primitiveId);
            if(primitive == NULL)
                continue;

            if(i == middleIndex)
            {
                for(j = 1; j < numPrimitivesPerLink; j++)
                {
                    prevPrimitiveId = primitiveId;
                    primitiveId += primitiveIndexDelta;

                    setPrimitiveToPreviousPrediction(figure, primitiveId,
                        prevPrimitiveId, value);
                }
            }

            else
            {
                for(j = 1; j < numPrimitivesPerLink; j++)
                {
                    prevPrimitiveId = primitiveId;
                    primitiveId += primitiveIndexDelta;
                    neighborPrimitiveId = primitiveId - deltaIndex;

                    setPrimitiveToAveragePrediction(figure, primitiveId,
                        prevPrimitiveId, neighborPrimitiveId,
                        1.0, value);
                }
            }
        }

        startIndex = middleIndex + neighborIndexDelta;
    }

    middleIndex = (numPrimitivesPerLink / 2) * primitiveIndexDelta;
    M3DPrimitive * newPrimitive = figure->getPrimitivePtr(middleIndex);
    M3DPrimitive * oldPrimitive = referenceFigure->getPrimitivePtr(middleIndex);
    Vector3D trans = oldPrimitive->getX() - newPrimitive->getX();

    figure->translateBy(trans, true);
}

void M3DFigureElongater::setPrimitiveToPreviousPrediction(
    M3DQuadFigure * targetFigure, int primitiveId, int prevPrimitiveId, double scale)
{
    M3DPrimitive * primitive,
                 * prevPrimitive,
                 * refPrimitive,
                 * refPrevPrimitive;

    Vector3D prevDiff;

    Quat q, dq;


    primitive = targetFigure->getPrimitivePtr(primitiveId);
    prevPrimitive = targetFigure->getPrimitivePtr(prevPrimitiveId);

    refPrimitive = referenceFigure->getPrimitivePtr(primitiveId);
    refPrevPrimitive = referenceFigure->getPrimitivePtr(prevPrimitiveId);

    if(primitive == NULL || prevPrimitive == NULL ||
       refPrimitive == NULL || refPrevPrimitive == NULL)
        return;

    prevDiff = scale * (refPrimitive->getX() - refPrevPrimitive->getX());

    dq = prevPrimitive->getQ() * refPrevPrimitive->getQ().conj();
    dq.rotateVector(prevDiff);

    q = refPrevPrimitive->getQ().conj() * refPrimitive->getQ();
/*
    // Take quaternion to the "elongate" power
    double angle = q.getW();
    if(angle < -1.0) angle = -1.0;
    if(angle > 1.0) angle = 1.0;
    angle = 2.0 * acos(angle) * scale;
    q.setAxisAngle(q.getVector(), angle);
*/
    primitive->setX(prevPrimitive->getX() + prevDiff);
    primitive->setQ(prevPrimitive->getQ() * q);
}

void M3DFigureElongater::setPrimitiveToAveragePrediction(
    M3DQuadFigure * targetFigure, int primitiveId, int prevPrimitiveId,
    int neighborPrimitiveId, double scale, double elongation)
{
    M3DPrimitive * primitive,
                 * prevPrimitive,
                 * neighborPrimitive,
                 * refPrimitive,
                 * refPrevPrimitive,
                 * refNeighborPrimitive;

    M3DQuadPrimitive prevPrediction,
                 neighborPrediction;

    Vector3D prevDiff,
             neighborDiff;
    double prevDist,
           neighborDist,
           weight;

    Quat q, dq;

    double angle;

    primitive = targetFigure->getPrimitivePtr(primitiveId);
    prevPrimitive = targetFigure->getPrimitivePtr(prevPrimitiveId);
    neighborPrimitive = targetFigure->getPrimitivePtr(neighborPrimitiveId);

    refPrimitive = referenceFigure->getPrimitivePtr(primitiveId);
    refPrevPrimitive = referenceFigure->getPrimitivePtr(prevPrimitiveId);
    refNeighborPrimitive = referenceFigure->getPrimitivePtr(neighborPrimitiveId);

    if(primitive == NULL || prevPrimitive == NULL ||
       neighborPrimitive == NULL || refPrimitive == NULL ||
       refPrevPrimitive == NULL || refNeighborPrimitive == NULL)
        return;

    // Calculate the averaging weight based on prev/neighbor distances
    prevDiff = elongation * (refPrimitive->getX() - refPrevPrimitive->getX());
    neighborDiff = scale * (refPrimitive->getX() - refNeighborPrimitive->getX());
//    prevDist = prevDiff.norm();
//    neighborDist = neighborDiff.norm();
//    weight = neighborDist / (prevDist + neighborDist);

    dq = prevPrimitive->getQ() * refPrevPrimitive->getQ().conj();

    // Take quaternion to the "elongate" power
    angle = dq.getW();
    if(angle < -1.0) angle = -1.0;
    if(angle > 1.0) angle = 1.0;
    angle = 2.0 * acos(angle) * elongation;
    dq.setAxisAngle(dq.getVector(), angle);

    dq.rotateVector(prevDiff);

    q = refPrevPrimitive->getQ().conj() * refPrimitive->getQ();

    prevPrediction.setX(prevPrimitive->getX() + prevDiff);
    prevPrediction.setQ(prevPrimitive->getQ() * q);
    prevPrediction.setR(primitive->getR());
    prevPrediction.setTheta(primitive->getTheta());

    dq = neighborPrimitive->getQ() * refNeighborPrimitive->getQ().conj();

    // Take quaternion to the "scale" power
//    angle = 2.0 * acos(dq.getW()) * scale;
//    dq.setAxisAngle(dq.getVector(), angle);

    dq.rotateVector(neighborDiff);

    q = refNeighborPrimitive->getQ().conj() * refPrimitive->getQ();

    neighborPrediction.setX(neighborPrimitive->getX() + neighborDiff);
    neighborPrediction.setQ(neighborPrimitive->getQ() * q);
    neighborPrediction.setR(primitive->getR());
    neighborPrediction.setTheta(primitive->getTheta());

    prevDiff = primitive->getX() - prevPrimitive->getX();
    neighborDiff = primitive->getX() - neighborPrimitive->getX();
    prevDist = prevDiff.norm();
    neighborDist = neighborDiff.norm();
    weight = elongation * neighborDist / (scale * prevDist + elongation * neighborDist);

    averagePrimitives(*primitive, prevPrediction, neighborPrediction, weight);
    //(*primitive) = neighborPrediction;
}

void M3DFigureElongater::averagePrimitives(M3DPrimitive & avePrim,
    const M3DPrimitive & prim1, const M3DPrimitive & prim2, double weight)
{
    Quat q;

    double weight2 = 1.0 - weight;

    avePrim.setX(weight2 * prim1.getX() + weight * prim2.getX());
    avePrim.setQ(q.slerp(weight, prim1.getQ(), prim2.getQ()));
    avePrim.setR(weight2 * prim1.getR() + weight * prim2.getR());
    avePrim.setTheta(weight2 * prim1.getTheta() + weight * prim2.getTheta());
}

