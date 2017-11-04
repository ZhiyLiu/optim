#include "M3DQuadInterpolator.h"
#include <math.h>
#include <stdio.h>

M3DQuadInterpolator::M3DQuadInterpolator(
	M3DPrimitive * p1, M3DPrimitive * p2,
	M3DPrimitive * p3, M3DPrimitive * p4,
	int xmesh, int ymesh)
{
	bez2D = NULL;
    radiusBez2D = NULL;

	initMedialPatch(p1, p2, p3, p4, xmesh, ymesh);
}

M3DQuadInterpolator::M3DQuadInterpolator(
	const M3DQuadFigure &figure, int startPrimId,
	int xmesh, int ymesh)
{
	bez2D = NULL;
    radiusBez2D = NULL;

	init(figure, startPrimId, xmesh, ymesh);
}

M3DQuadInterpolator::~M3DQuadInterpolator()
{
	if(bez2D != NULL)
		delete bez2D;

    if(radiusBez2D == NULL)
        delete radiusBez2D;
}

bool M3DQuadInterpolator::init(const M3DQuadFigure &figure, int startPrimId,
								 int xmesh, int ymesh)
{
	M3DPrimitive * primPtr1,
                 * primPtr2,
                 * primPtr3,
                 * primPtr4;

    int primsPerRow;

	if(bez2D != NULL)
		delete bez2D;
    if(radiusBez2D != NULL)
        delete radiusBez2D;

	bez2D = NULL;
    radiusBez2D = NULL;

    primsPerRow = figure.getColumnCount();

	primPtr1 = figure.getPrimitivePtr(startPrimId);
	primPtr2 = figure.getPrimitivePtr(startPrimId + 1);
	primPtr3 = figure.getPrimitivePtr(startPrimId + primsPerRow);
	primPtr4 = figure.getPrimitivePtr(startPrimId + primsPerRow + 1);

    if(primPtr1 == NULL || primPtr2 == NULL ||
       primPtr3 == NULL || primPtr4 == NULL)
       return false;

	initMedialPatch(primPtr1, primPtr2, primPtr3, primPtr4, xmesh, ymesh);

	return true;
}

// Given a quad of primitives and a u and v, interpolates primitive
// along the interpolated, parameterized medial sheet (Bezier patch)
M3DPrimitive * M3DQuadInterpolator::interpolatePrim(double u, double v)
{
    M3DPrimitive * primPtr;

    Vector3D du, dv,
             b, bPerp, n;

	if(bez2D == NULL || radiusBez2D == NULL)
		return NULL;

    primPtr = new M3DQuadPrimitive;

	primPtr->setX(bez2D->calcPoint(u, v));
    primPtr->setR(radiusBez2D->calcPoint(u, v).getX());

    du = bez2D->calcDU(u,v);
    dv = bez2D->calcDV(u,v);

    n = du.cross(dv);
    n.normalize();

    b = -(du * radiusBez2D->calcDU(u, v).getX() +
          dv * radiusBez2D->calcDV(u, v).getX());
    double length = b.normalize();

    bPerp = b.cross(n);

    Quat q(b, n, bPerp);
    primPtr->setQ(q);

    primPtr->setTheta(acos(length));

    return primPtr;
}

M3DPrimitive * M3DQuadInterpolator::getPrimitive(int i, int j)
{
    M3DPrimitive * primPtr;

	double u, v;

    Vector3D du, dv,
             b, bPerp, n;

	if(bez2D == NULL || radiusBez2D == NULL)
		return NULL;

	u = (double) i / (double) bez2D->getXMesh();
	v = (double) j / (double) bez2D->getYMesh();

    primPtr = new M3DQuadPrimitive;

	primPtr->setX(bez2D->getPoint(i, j));
    primPtr->setR(radiusBez2D->getPoint(i, j).getX());

    du = bez2D->calcDU(u,v);
    dv = bez2D->calcDV(u,v);

    n = du.cross(dv);
    n.normalize();

    b = -(du * radiusBez2D->calcDU(u, v).getX() +
          dv * radiusBez2D->calcDV(u, v).getX());
    double length = b.normalize();

    bPerp = b.cross(n);

    Quat q(b, n, bPerp);
    primPtr->setQ(q);

    primPtr->setTheta(acos(length));

    return primPtr;
}

void M3DQuadInterpolator::initMedialPatch(M3DPrimitive * prim1,
                                          M3DPrimitive * prim2,
                                          M3DPrimitive * prim3,
                                          M3DPrimitive * prim4,
                                          int xmesh, int ymesh)
{
	Vector3D ctrlPts[4][4];

	ctrlPts[0][0] = prim1->getX();
	ctrlPts[0][3] = prim2->getX();
	ctrlPts[3][0] = prim3->getX();
	ctrlPts[3][3] = prim4->getX();

	const double third = 1.0 / 3.0;

	// Difference between two primitive positions
	Vector3D diff;
	// Projection of the difference vector onto the b-n-plane
	Vector3D proj;

	diff = prim2->getX() - prim1->getX();
	proj = (prim1->getB() * diff) * prim1->getB() +
		   (prim1->getBPerp() * diff) * prim1->getBPerp();
	ctrlPts[0][1] = prim1->getX() + (proj * third);
	proj = (prim2->getB() * diff) * prim2->getB() +
		   (prim2->getBPerp() * diff) * prim2->getBPerp();
	ctrlPts[0][2] = prim2->getX() - (proj * third);

	diff = prim3->getX() - prim1->getX();
	proj = (prim1->getB() * diff) * prim1->getB() +
		   (prim1->getBPerp() * diff) * prim1->getBPerp();
	ctrlPts[1][0] = prim1->getX() + (proj * third);
	proj = (prim3->getB() * diff) * prim3->getB() +
		   (prim3->getBPerp() * diff) * prim3->getBPerp();
	ctrlPts[2][0] = prim3->getX() - (proj * third);

	diff = prim4->getX() - prim2->getX();
	proj = (prim2->getB() * diff) * prim2->getB() +
		   (prim2->getBPerp() * diff) * prim2->getBPerp();
	ctrlPts[1][3] = prim2->getX() + (proj * third);
	proj = (prim4->getB() * diff) * prim4->getB() +
		   (prim4->getBPerp() * diff) * prim4->getBPerp();
	ctrlPts[2][3] = prim4->getX() - (proj * third);

	diff = prim4->getX() - prim3->getX();
	proj = (prim3->getB() * diff) * prim3->getB() +
		   (prim3->getBPerp() * diff) * prim3->getBPerp();
	ctrlPts[3][1] = prim3->getX() + (proj * third);
	proj = (prim4->getB() * diff) * prim4->getB() +
		   (prim4->getBPerp() * diff) * prim4->getBPerp();
	ctrlPts[3][2] = prim4->getX() - (proj * third);

	ctrlPts[1][1] = ctrlPts[0][1] + ctrlPts[1][0] - prim1->getX();
	ctrlPts[1][2] = ctrlPts[0][2] + ctrlPts[1][3] - prim2->getX();
	ctrlPts[2][1] = ctrlPts[2][0] + ctrlPts[3][1] - prim3->getX();
	ctrlPts[2][2] = ctrlPts[3][2] + ctrlPts[2][3] - prim4->getX();

	if(bez2D != NULL)
		delete bez2D;

	bez2D = new Bezier2D(ctrlPts, xmesh, ymesh);

    ctrlPts[0][0] = Vector3D(prim1->getR(), 0.0, 0.0);
    ctrlPts[0][3] = Vector3D(prim2->getR(), 0.0, 0.0);
    ctrlPts[3][0] = Vector3D(prim3->getR(), 0.0, 0.0);
    ctrlPts[3][3] = Vector3D(prim4->getR(), 0.0, 0.0);

    diff = prim2->getX() - prim1->getX();
    diff.normalize();
    ctrlPts[0][1] = Vector3D(prim1->getR() + (prim1->getB() * diff) * -cos(prim1->getTheta()), 0.0, 0.0);
    ctrlPts[0][2] = Vector3D(prim2->getR() + (prim2->getB() * -diff) * -cos(prim2->getTheta()), 0.0, 0.0);

    diff = prim3->getX() - prim1->getX();
    diff.normalize();
    ctrlPts[1][0] = Vector3D(prim1->getR() + (prim1->getB() * diff) * -cos(prim1->getTheta()), 0.0, 0.0);
    ctrlPts[2][0] = Vector3D(prim3->getR() + (prim3->getB() * -diff) * -cos(prim3->getTheta()), 0.0, 0.0);

    diff = prim4->getX() - prim2->getX();
    diff.normalize();
    ctrlPts[1][3] = Vector3D(prim2->getR() + (prim2->getB() * diff) * -cos(prim2->getTheta()), 0.0, 0.0);
    ctrlPts[2][3] = Vector3D(prim4->getR() + (prim4->getB() * -diff) * -cos(prim4->getTheta()), 0.0, 0.0);

    diff = prim4->getX() - prim3->getX();
    diff.normalize();
    ctrlPts[3][1] = Vector3D(prim3->getR() + (prim3->getB() * diff) * -cos(prim3->getTheta()), 0.0, 0.0);
    ctrlPts[3][2] = Vector3D(prim4->getR() + (prim4->getB() * -diff) * -cos(prim4->getTheta()), 0.0, 0.0);

	ctrlPts[1][1] = (ctrlPts[0][1] + ctrlPts[1][0]) * 0.5;// - ctrlPts[0][0];
	ctrlPts[1][2] = (ctrlPts[0][2] + ctrlPts[1][3]) * 0.5;// - ctrlPts[0][3];
	ctrlPts[2][1] = (ctrlPts[2][0] + ctrlPts[3][1]) * 0.5;// - ctrlPts[3][0];
	ctrlPts[2][2] = (ctrlPts[3][2] + ctrlPts[2][3]) * 0.5;// - ctrlPts[3][3];

    printf("%f %f %f %f\n", ctrlPts[0][0].getX(), ctrlPts[0][1].getX(), ctrlPts[0][2].getX(), ctrlPts[0][3].getX());
    printf("%f %f %f %f\n", ctrlPts[1][0].getX(), ctrlPts[1][1].getX(), ctrlPts[1][2].getX(), ctrlPts[1][3].getX());
    printf("%f %f %f %f\n", ctrlPts[2][0].getX(), ctrlPts[2][1].getX(), ctrlPts[2][2].getX(), ctrlPts[2][3].getX());
    printf("%f %f %f %f\n", ctrlPts[3][0].getX(), ctrlPts[3][1].getX(), ctrlPts[3][2].getX(), ctrlPts[3][3].getX());
	if(radiusBez2D != NULL)
		delete radiusBez2D;

	radiusBez2D = new Bezier2D(ctrlPts, xmesh, ymesh);
}


