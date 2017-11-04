#ifndef BEZIER2D_H
#define BEZIER2D_H

#include "Vector3D.h"

class Bezier2D
{
public:
	Bezier2D();
	Bezier2D(double * _ctrlPts, int xmesh = 0, int ymesh = 0);
	Bezier2D(Vector3D _ctrlPts[4][4], int xmesh = 0, int ymesh = 0);

	Vector3D &getPoint(int i, int j)
	{
		//assert(i <= xMesh + 1 && i >= 0 && j <= yMesh + 1 && j >= 0);
		return points[(yMesh + 1) * i + j];
	}

	Vector3D &getNormal(int i, int j)
	{
		//assert(i <= xMesh + 1 && i >= 0 && j <= yMesh + 1 && j >= 0);
		return normals[(yMesh + 1) * i + j];
	}

	Vector3D calcPoint(double u, double v);
	Vector3D calcDU(double u, double v);
	Vector3D calcDV(double u, double v);
	Vector3D calcNormal(double u, double v);

	int getXMesh() { return xMesh; }
	int getYMesh() { return yMesh; }

	~Bezier2D();

private:
	void calcCoeffMat();
	void calcPoints();
	void initMesh();

	int xMesh;
	int yMesh;

	Vector3D ctrlPts[4][4];
	Vector3D coeffMat[4][4];
	Vector3D * points;
	Vector3D * normals;
};

#endif

