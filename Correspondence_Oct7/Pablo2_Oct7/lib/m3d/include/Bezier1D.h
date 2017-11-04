#ifndef BEZIER1D_H
#define BEZIER1D_H

#include "Vector3D.h"

class Bezier1D
{
public:
	Bezier1D();
	Bezier1D(double * _ctrlPts, int _mesh = 0);
	Bezier1D(Vector3D _ctrlPts[4], int _mesh = 0);

	void newMesh(Vector3D _ctrlPts[4], int _mesh);

	Vector3D &getPoint(int i)
	{
//		assert(i <= mesh + 1 && i >= 0);
		return points[i];
	}

	Vector3D &getNormal(int i)
	{
//		assert(i <= mesh + 1 && i >= 0);
		return tangents[i];
	}

	Vector3D calcPoint(double t);
	Vector3D calcTangent(double t);

	int getMesh() { return mesh; }

	~Bezier1D();

private:
	void calcCoeffMat();
	void calcPoints();
	void initMesh();

	Vector3D ctrlPts[4];
	Vector3D coeffMat[4];
	Vector3D * points;
	Vector3D * tangents;

	int mesh;
};

#endif

