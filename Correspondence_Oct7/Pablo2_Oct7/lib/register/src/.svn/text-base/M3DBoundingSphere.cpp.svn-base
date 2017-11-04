#include <math.h>
#include "M3DBoundingSphere.h"
#include "M3DFigure.h"

#define D_XFERLIST
#define D_POINTLIST_SERVER2
#include "Shapedepend.h"

M3DBoundingSphere::M3DBoundingSphere(const M3DObject & object, int subdivisions)
{
    int figureCount;
    int numPoints;
    int i, j;

    M3DFigure * figure;
    Xferlist * xferList;
    ThallCode::Pointlist_server2 * pList;

    Bpoint * bPoints;
    Bpoint * bPointPtr;
    Vector3D pos;
    Vector3D diff;
    double distSqr, maxDistSqr;


    center = object.getCOG();

    maxDistSqr = 0.0;

    figureCount = object.getFigureCount();
    for(i = 0; i < figureCount; i++)
    {
        figure = object.getFigurePtr(i);

        if(figure == NULL)
            continue;

        xferList = convertM3DtoXfer(figure);

        pList = new ThallCode::Pointlist_server2;
        pList->init(xferList);
        pList->ComputeSubdivPointCloud(subdivisions);
        pList->subdivboundaryinfo(&numPoints, &bPoints);

        bPointPtr = bPoints;
        for(j = 0; j < numPoints; j++, bPointPtr++)
        {
            pos.set(bPointPtr->pnt[0], bPointPtr->pnt[1], bPointPtr->pnt[2]);
            diff = (pos - center);
            distSqr = diff * diff;

            if(distSqr > maxDistSqr)
                maxDistSqr = distSqr;
        }
    }

    radius = sqrt(maxDistSqr);
}

M3DBoundingSphere::M3DBoundingSphere(const M3DFigure & figure, int subdivisions)
{
    center = figure.getCOG();
}

