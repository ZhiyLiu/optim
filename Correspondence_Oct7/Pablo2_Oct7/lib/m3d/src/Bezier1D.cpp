#ifdef __APPLE__
#include <sys/malloc.h>
#else
#include <malloc.h>
#endif
#include "Bezier1D.h"

Bezier1D::Bezier1D()
{
	points = NULL;
	tangents = NULL;
}

Bezier1D::Bezier1D(double * _ctrlPts, int _mesh)
{
	points = NULL;
	tangents = NULL;

	mesh = _mesh;

	for(int i = 0; i <= 9; i += 3)
		ctrlPts[0].set(_ctrlPts + i);

	initMesh();
}

Bezier1D::Bezier1D(Vector3D _ctrlPts[4], int _mesh)
{
	points = NULL;
	tangents = NULL;

	mesh = _mesh;

	ctrlPts[0] = _ctrlPts[0];
	ctrlPts[1] = _ctrlPts[1];
	ctrlPts[2] = _ctrlPts[2];
	ctrlPts[3] = _ctrlPts[3];

	initMesh();
}

void Bezier1D::newMesh(Vector3D _ctrlPts[4], int _mesh)
{
	points = NULL;
	tangents = NULL;

	mesh = _mesh;

	ctrlPts[0] = _ctrlPts[0];
	ctrlPts[1] = _ctrlPts[1];
	ctrlPts[2] = _ctrlPts[2];
	ctrlPts[3] = _ctrlPts[3];

	initMesh();
}

void Bezier1D::initMesh()
{
	if(mesh > 0)
	{
		points = (Vector3D *) malloc(sizeof(Vector3D) * (mesh + 1));
		tangents = (Vector3D *) malloc(sizeof(Vector3D) * (mesh + 1));
	}

	calcCoeffMat();
	calcPoints();
}

void Bezier1D::calcCoeffMat()
{
	coeffMat[0] = ctrlPts[0];
	coeffMat[1] = 3 * (ctrlPts[1] - ctrlPts[0]);
	coeffMat[2] = 3 * (ctrlPts[0] + ctrlPts[2]) - 6 * ctrlPts[1];
	coeffMat[3] = ctrlPts[3] - ctrlPts[0] + 3 * (ctrlPts[1] - ctrlPts[2]);
}

void Bezier1D::calcPoints()
{
	double t;
	int i;

	double invMesh = 1.0 / (double) mesh;
	for(i = 0; i <= mesh; i++)
	{
		t = (double) i * invMesh;
		points[i] = coeffMat[0] + t * (coeffMat[1] + t * (coeffMat[2] + t * coeffMat[3]));
		tangents[i] = coeffMat[1] + t * (2 * coeffMat[2] + (3 * t) * coeffMat[3]);
	}
}

Vector3D Bezier1D::calcPoint(double t)
{
	return coeffMat[0] + t * (coeffMat[1] + t * (coeffMat[2] + t * coeffMat[3]));
}

Vector3D Bezier1D::calcTangent(double t)
{
	return coeffMat[1] + t * (2 * coeffMat[2] + (3 * t) * coeffMat[3]);
}

Bezier1D::~Bezier1D()
{
	if(points != NULL)
		free((void *) points);

	if(tangents != NULL)
		free((void *) tangents);
}



