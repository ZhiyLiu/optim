#ifdef __APPLE__
#include <sys/malloc.h>
#else
#include <malloc.h>
#endif
#include "Bezier2D.h"

Bezier2D::Bezier2D()
{
	points = NULL;
	normals = NULL;

	xMesh = 0;
	yMesh = 0;
}

Bezier2D::Bezier2D(double * _ctrlPts, int xmesh, int ymesh)
{
	double *ptr = _ctrlPts;

	int i, j;

	for(i = 0; i < 4; i++)
	{
		for(j = 0; j < 4; j++)
		{
			ctrlPts[i][j].set(ptr);
			ptr += 3;
		}
	}

	xMesh = xmesh;
	yMesh = ymesh;

	points = NULL;
	normals = NULL;

	initMesh();
}

Bezier2D::Bezier2D(Vector3D _ctrlPts[4][4], int xmesh, int ymesh)
{
	int i, j;

	for(i = 0; i < 4; i++)
	{
		for(j = 0; j < 4; j++)
			ctrlPts[i][j] = _ctrlPts[i][j];
	}

	xMesh = xmesh;
	yMesh = ymesh;

	points = NULL;
	normals = NULL;

	initMesh();
}

void Bezier2D::calcCoeffMat()
{
	coeffMat[0][0] = ctrlPts[0][0];
	coeffMat[0][1] = 3 * (ctrlPts[0][1] - ctrlPts[0][0]);
	coeffMat[0][2] = 3 * (ctrlPts[0][0] - 2 * ctrlPts[0][1] + ctrlPts[0][2]);
	coeffMat[0][3] = ctrlPts[0][3] - ctrlPts[0][0] + 3 * (ctrlPts[0][1] - ctrlPts[0][2]);

	coeffMat[1][0] = 3 * (ctrlPts[1][0] - ctrlPts[0][0]);
	coeffMat[1][1] = 9 * (ctrlPts[0][0] - ctrlPts[0][1] - ctrlPts[1][0] + ctrlPts[1][1]);
	coeffMat[1][2] = 9 * (ctrlPts[1][0] + ctrlPts[1][2] - ctrlPts[0][0] - ctrlPts[0][2]) +
		             18 * (ctrlPts[0][1] - ctrlPts[1][1]);
	coeffMat[1][3] = 3 * (ctrlPts[0][0] - ctrlPts[0][3] - ctrlPts[1][0] + ctrlPts[1][3]) +
		             9 * (ctrlPts[0][2] - ctrlPts[0][1] + ctrlPts[1][1] - ctrlPts[1][2]);

	coeffMat[2][0] = 3 * (ctrlPts[0][0] - 2 * ctrlPts[1][0] + ctrlPts[2][0]);
	coeffMat[2][1] = 9 * (ctrlPts[0][1] - ctrlPts[0][0] + ctrlPts[2][1] - ctrlPts[2][0]) +
		             18 * (ctrlPts[1][0] - ctrlPts[1][1]);
	coeffMat[2][2] = 9 * (ctrlPts[0][0] + ctrlPts[0][2] + ctrlPts[2][0] + ctrlPts[2][2]) +
		             -18 * (ctrlPts[0][1] + ctrlPts[1][0] + ctrlPts[1][2] + ctrlPts[2][1]) +
					 36 * ctrlPts[1][1];
	coeffMat[2][3] = 3 * (ctrlPts[0][3] - ctrlPts[0][0] + ctrlPts[2][3] - ctrlPts[2][0]) +
		             6 * (ctrlPts[1][0] - ctrlPts[1][3]) +
		             9 * (ctrlPts[0][1] - ctrlPts[0][2] + ctrlPts[2][1] - ctrlPts[2][2]) +
					 18 * (ctrlPts[1][2] - ctrlPts[1][1]);

	coeffMat[3][0] = ctrlPts[3][0] - ctrlPts[0][0] + 3 * (ctrlPts[1][0] - ctrlPts[2][0]);
	coeffMat[3][1] = 3 * (ctrlPts[0][0] - ctrlPts[0][1] + ctrlPts[3][1] - ctrlPts[3][0]) + 
		             9 * (ctrlPts[1][1] - ctrlPts[1][0] + ctrlPts[2][0] - ctrlPts[2][1]);
	coeffMat[3][2] = 3 * (ctrlPts[3][2] + ctrlPts[3][0] - ctrlPts[0][0] - ctrlPts[0][2]) +
		             6 * (ctrlPts[0][1] - ctrlPts[3][1]) +
					 9 * (ctrlPts[1][0] + ctrlPts[1][2] - ctrlPts[2][0] - ctrlPts[2][2]) +
					 18 * (ctrlPts[2][1] - ctrlPts[1][1]);
	coeffMat[3][3] = ctrlPts[0][0] - ctrlPts[0][3] + ctrlPts[3][3] - ctrlPts[3][0] +
		             3 * (ctrlPts[0][2] - ctrlPts[0][1] + ctrlPts[3][1] - ctrlPts[3][2] - 
					      ctrlPts[1][0] + ctrlPts[1][3] + ctrlPts[2][0] - ctrlPts[2][3]) +
					 9 * (ctrlPts[1][1] - ctrlPts[1][2] + ctrlPts[2][2] - ctrlPts[2][1]);
}

void Bezier2D::initMesh()
{
	calcCoeffMat();

	if(xMesh > 0 && yMesh > 0)
	{
		// Allocate the array of points and normals
		int size = (xMesh + 1) * (yMesh + 1);

		points = (Vector3D *) malloc(sizeof(Vector3D) * size);
		normals = (Vector3D *) malloc(sizeof(Vector3D) * size);

		calcPoints();
	}
}

void Bezier2D::calcPoints()
{
	Vector3D temp[4];

	// Indexing variables
	int i, j;
	int index = 0;
	int dIndex = 0;

	// Points to be evaluated
	double u, v;

	// Derivatives at a point
	Vector3D du, dv;

	Vector3D * dTemp = new Vector3D[(yMesh + 1) * 3];

	double invXMesh = 1.0 / (double) xMesh;
	double invYMesh = 1.0 / (double) yMesh;

	// Calculate v components of du
	for(i = 0; i <= yMesh; i++)
	{
		v = (double) i * invYMesh;

		dTemp[dIndex] = coeffMat[0][1] + v * (coeffMat[1][1] + v * (coeffMat[2][1] + v * coeffMat[3][1]));
		dTemp[dIndex+1] = coeffMat[0][2] + v * (coeffMat[1][2] + v * (coeffMat[2][2] + v * coeffMat[3][2]));
		dTemp[dIndex+2] = coeffMat[0][3] + v * (coeffMat[1][3] + v * (coeffMat[2][3] + v * coeffMat[3][3]));

		dIndex += 3;
	}

	// Calculate the points
	for(i = 0; i <= xMesh; i++)
	{
		u = (double) i * invXMesh;

		temp[0] = coeffMat[0][0] + u * (coeffMat[0][1] + u * (coeffMat[0][2] + u * coeffMat[0][3]));
		temp[1] = coeffMat[1][0] + u * (coeffMat[1][1] + u * (coeffMat[1][2] + u * coeffMat[1][3]));
		temp[2] = coeffMat[2][0] + u * (coeffMat[2][1] + u * (coeffMat[2][2] + u * coeffMat[2][3]));
		temp[3] = coeffMat[3][0] + u * (coeffMat[3][1] + u * (coeffMat[3][2] + u * coeffMat[3][3]));

		dIndex = 0;
		for(j = 0; j <= yMesh; j++)
		{
			v = (double) j * invYMesh;

			points[index] = temp[0] + v * (temp[1] + v * (temp[2] + v * temp[3]));

			dv = temp[1] + v * (2 * temp[2] + (3 * v) * temp[3]);

			du = dTemp[dIndex] + u * (2 * dTemp[dIndex+1] + (3 * u) * dTemp[dIndex+2]);
			dIndex += 3;

			normals[index] = du.cross(dv);
			normals[index].normalize();

			index++;
		}
	}

	delete [] dTemp;
}

Vector3D Bezier2D::calcPoint(double u, double v)
{
	Vector3D temp[4];

	temp[0] = coeffMat[0][0] + u * (coeffMat[0][1] + u * (coeffMat[0][2] + u * coeffMat[0][3]));
	temp[1] = coeffMat[1][0] + u * (coeffMat[1][1] + u * (coeffMat[1][2] + u * coeffMat[1][3]));
	temp[2] = coeffMat[2][0] + u * (coeffMat[2][1] + u * (coeffMat[2][2] + u * coeffMat[2][3]));
	temp[3] = coeffMat[3][0] + u * (coeffMat[3][1] + u * (coeffMat[3][2] + u * coeffMat[3][3]));

	return temp[0] + v * (temp[1] + v * (temp[2] + v * temp[3]));
}

Vector3D Bezier2D::calcDU(double u, double v)
{
	Vector3D temp[3];

	temp[0] = coeffMat[0][1] + v * (coeffMat[1][1] + v * (coeffMat[2][1] + v * coeffMat[3][1]));
	temp[1] = coeffMat[0][2] + v * (coeffMat[1][2] + v * (coeffMat[2][2] + v * coeffMat[3][2]));
	temp[2] = coeffMat[0][3] + v * (coeffMat[1][3] + v * (coeffMat[2][3] + v * coeffMat[3][3]));

	return temp[0] + u * (2 * temp[1] + (3 * u) * temp[2]);
}

Vector3D Bezier2D::calcDV(double u, double v)
{
	Vector3D temp[3];

	temp[0] = coeffMat[1][0] + u * (coeffMat[1][1] + u * (coeffMat[1][2] + u * coeffMat[1][3]));
	temp[1] = coeffMat[2][0] + u * (coeffMat[2][1] + u * (coeffMat[2][2] + u * coeffMat[2][3]));
	temp[2] = coeffMat[3][0] + u * (coeffMat[3][1] + u * (coeffMat[3][2] + u * coeffMat[3][3]));

	return temp[0] + v * (2 * temp[1] + (3 * v) * temp[2]);
}

Vector3D Bezier2D::calcNormal(double u, double v)
{
	Vector3D du = calcDU(u, v);
	Vector3D dv = calcDV(u, v);

	return du.cross(dv);
}

Bezier2D::~Bezier2D() 
{
	if(points != NULL)
		free((void *) points);

	if(normals != NULL)
		free((void *) normals);
}



