#include <math.h>
#include "Matrix4D.h"

using namespace std;

Matrix4D::Matrix4D()
{
    setToIdentity();
}

Matrix4D::Matrix4D(double * mat)
{
    int i;

    for(i = 0; i < 16; i++)
        matrix[i] = mat[i];
}

void Matrix4D::setData(double * mat)
{
    int i;

    for(i = 0; i < 16; i++)
        matrix[i] = mat[i];
}

Matrix4D::Matrix4D(float * mat)
{
    int i;

    for(i = 0; i < 16; i++)
        matrix[i] = mat[i];
}

void Matrix4D::setData(float * mat)
{
    int i;

    for(i = 0; i < 16; i++)
        matrix[i] = mat[i];
}

void Matrix4D::print(ostream & out) const
{
    out << matrix[0] << '\t' << matrix[1] << '\t' << matrix[2] << '\t' << matrix[3] << endl;
    out << matrix[4] << '\t' << matrix[5] << '\t' << matrix[6] << '\t' << matrix[7] << endl;
    out << matrix[8] << '\t' << matrix[9] << '\t' << matrix[10] << '\t' << matrix[11] << endl;
    out << matrix[12] << '\t' << matrix[13] << '\t' << matrix[14] << '\t' << matrix[15] << endl;
}

void Matrix4D::setToIdentity()
{
    matrix[0] = 1.0;
    matrix[1] = 0.0;
    matrix[2] = 0.0;
    matrix[3] = 0.0;

    matrix[4] = 0.0;
    matrix[5] = 1.0;
    matrix[6] = 0.0;
    matrix[7] = 0.0;

    matrix[8] = 0.0;
    matrix[9] = 0.0;
    matrix[10] = 1.0;
    matrix[11] = 0.0;

    matrix[12] = 0.0;
    matrix[13] = 0.0;
    matrix[14] = 0.0;
    matrix[15] = 1.0;
}

void Matrix4D::translate(const Vector3D &v)
{
    matrix[3] = matrix[3] + v.getX();
    matrix[7] = matrix[7] + v.getY();
    matrix[11] = matrix[11] + v.getZ();
}

void Matrix4D::rotate(const Quat &q)
{
    double mat[16];

    for(int i = 0; i < 16; i++)
        mat[i] = matrix[i];

    q.buildRotMatrix(&(matrix[0]));

    multiply(&(mat[0]));
}

void Matrix4D::scale(double s)
{
    matrix[0] = matrix[0] * s;
    matrix[1] = matrix[1] * s;
    matrix[2] = matrix[2] * s;
    matrix[3] = matrix[3] * s;

    matrix[4] = matrix[4] * s;
    matrix[5] = matrix[5] * s;
    matrix[6] = matrix[6] * s;
    matrix[7] = matrix[7] * s;

    matrix[8] = matrix[8] * s;
    matrix[9] = matrix[9] * s;
    matrix[10] = matrix[10] * s;
    matrix[11] = matrix[11] * s;
}

void Matrix4D::decompose(Quat & rotation, double & scale, Vector3D & translation) const
{
	double m[16];
	int i;

	for (i = 0; i < 16; i++)
		m[i] = matrix[i];

	translation = Vector3D(m[3], m[7], m[11]);
	m[3] = 0.0;
	m[7] = 0.0;
	m[11] = 0.0;

	double d = det();
	scale = pow(d, 1.0/3.0);

	for (i = 0; i < 11; i++)
		m[i] /= scale;

	rotation = Quat();
	rotation.matrixToQuat(m);
}

double Matrix4D::det() const
{
	return matrix[0]*(matrix[5]*matrix[10] - matrix[6]*matrix[9])
		- matrix[1]*(matrix[4]*matrix[10] - matrix[6]*matrix[8])
		+ matrix[2]*(matrix[4]*matrix[9] - matrix[5]*matrix[8]);
}

inline void swap(double & a, double & b)
{
	double d = a;
	a = b;
	b = d;
}

void Matrix4D::transpose()
{
	swap(matrix[1], matrix[4]);
	swap(matrix[2], matrix[8]);
	swap(matrix[3], matrix[12]);
	swap(matrix[6], matrix[9]);
	swap(matrix[7], matrix[13]);
	swap(matrix[11], matrix[14]);
}

void Matrix4D::multiply(double * mat)
{
    double matrix2[16];

    for(int i = 0; i < 16; i++)
        matrix2[i] = matrix[i];

    matrix[0] = matrix2[0] * mat[0] + matrix2[1] * mat[4] + matrix2[2] * mat[8] + matrix2[3] * mat[12];
    matrix[1] = matrix2[0] * mat[1] + matrix2[1] * mat[5] + matrix2[2] * mat[9] + matrix2[3] * mat[13];
    matrix[2] = matrix2[0] * mat[2] + matrix2[1] * mat[6] + matrix2[2] * mat[10] + matrix2[3] * mat[14];
    matrix[3] = matrix2[0] * mat[3] + matrix2[1] * mat[7] + matrix2[2] * mat[11] + matrix2[3] * mat[15];

    matrix[4] = matrix2[4] * mat[0] + matrix2[5] * mat[4] + matrix2[6] * mat[8] + matrix2[7] * mat[12];
    matrix[5] = matrix2[4] * mat[1] + matrix2[5] * mat[5] + matrix2[6] * mat[9] + matrix2[7] * mat[13];
    matrix[6] = matrix2[4] * mat[2] + matrix2[5] * mat[6] + matrix2[6] * mat[10] + matrix2[7] * mat[14];
    matrix[7] = matrix2[4] * mat[3] + matrix2[5] * mat[7] + matrix2[6] * mat[11] + matrix2[7] * mat[15];

    matrix[8] = matrix2[8] * mat[0] + matrix2[9] * mat[4] + matrix2[10] * mat[8] + matrix2[11] * mat[12];
    matrix[9] = matrix2[8] * mat[1] + matrix2[9] * mat[5] + matrix2[10] * mat[9] + matrix2[11] * mat[13];
    matrix[10] = matrix2[8] * mat[2] + matrix2[9] * mat[6] + matrix2[10] * mat[10] + matrix2[11] * mat[14];
    matrix[11] = matrix2[8] * mat[3] + matrix2[9] * mat[7] + matrix2[10] * mat[11] + matrix2[11] * mat[15];

    matrix[12] = matrix2[12] * mat[0] + matrix2[13] * mat[4] + matrix2[14] * mat[8] + matrix2[15] * mat[12];
    matrix[13] = matrix2[12] * mat[1] + matrix2[13] * mat[5] + matrix2[14] * mat[9] + matrix2[15] * mat[13];
    matrix[14] = matrix2[12] * mat[2] + matrix2[13] * mat[6] + matrix2[14] * mat[10] + matrix2[15] * mat[14];
    matrix[15] = matrix2[12] * mat[3] + matrix2[13] * mat[7] + matrix2[14] * mat[11] + matrix2[15] * mat[15];
}

void Matrix4D::multiply(const Matrix4D & mat)
{
    double matrix2[16];

    for(int i = 0; i < 16; i++)
        matrix2[i] = matrix[i];

    matrix[0] = matrix2[0] * mat.matrix[0] + matrix2[1] * mat.matrix[4] + matrix2[2] * mat.matrix[8] + matrix2[3] * mat.matrix[12];
    matrix[1] = matrix2[0] * mat.matrix[1] + matrix2[1] * mat.matrix[5] + matrix2[2] * mat.matrix[9] + matrix2[3] * mat.matrix[13];
    matrix[2] = matrix2[0] * mat.matrix[2] + matrix2[1] * mat.matrix[6] + matrix2[2] * mat.matrix[10] + matrix2[3] * mat.matrix[14];
    matrix[3] = matrix2[0] * mat.matrix[3] + matrix2[1] * mat.matrix[7] + matrix2[2] * mat.matrix[11] + matrix2[3] * mat.matrix[15];

    matrix[4] = matrix2[4] * mat.matrix[0] + matrix2[5] * mat.matrix[4] + matrix2[6] * mat.matrix[8] + matrix2[7] * mat.matrix[12];
    matrix[5] = matrix2[4] * mat.matrix[1] + matrix2[5] * mat.matrix[5] + matrix2[6] * mat.matrix[9] + matrix2[7] * mat.matrix[13];
    matrix[6] = matrix2[4] * mat.matrix[2] + matrix2[5] * mat.matrix[6] + matrix2[6] * mat.matrix[10] + matrix2[7] * mat.matrix[14];
    matrix[7] = matrix2[4] * mat.matrix[3] + matrix2[5] * mat.matrix[7] + matrix2[6] * mat.matrix[11] + matrix2[7] * mat.matrix[15];

    matrix[8] = matrix2[8] * mat.matrix[0] + matrix2[9] * mat.matrix[4] + matrix2[10] * mat.matrix[8] + matrix2[11] * mat.matrix[12];
    matrix[9] = matrix2[8] * mat.matrix[1] + matrix2[9] * mat.matrix[5] + matrix2[10] * mat.matrix[9] + matrix2[11] * mat.matrix[13];
    matrix[10] = matrix2[8] * mat.matrix[2] + matrix2[9] * mat.matrix[6] + matrix2[10] * mat.matrix[10] + matrix2[11] * mat.matrix[14];
    matrix[11] = matrix2[8] * mat.matrix[3] + matrix2[9] * mat.matrix[7] + matrix2[10] * mat.matrix[11] + matrix2[11] * mat.matrix[15];

    matrix[12] = matrix2[12] * mat.matrix[0] + matrix2[13] * mat.matrix[4] + matrix2[14] * mat.matrix[8] + matrix2[15] * mat.matrix[12];
    matrix[13] = matrix2[12] * mat.matrix[1] + matrix2[13] * mat.matrix[5] + matrix2[14] * mat.matrix[9] + matrix2[15] * mat.matrix[13];
    matrix[14] = matrix2[12] * mat.matrix[2] + matrix2[13] * mat.matrix[6] + matrix2[14] * mat.matrix[10] + matrix2[15] * mat.matrix[14];
    matrix[15] = matrix2[12] * mat.matrix[3] + matrix2[13] * mat.matrix[7] + matrix2[14] * mat.matrix[11] + matrix2[15] * mat.matrix[15];
}

Vector3D Matrix4D::multiplyVector(const Vector3D & v) const
{
    Vector3D retVector;

    double w = matrix[12] * v.getX() + matrix[13] * v.getY() + matrix[14] * v.getZ() + matrix[15];

    retVector.setX((matrix[0] * v.getX() + matrix[1] * v.getY() + matrix[2] * v.getZ() + matrix[3]) / w);
    retVector.setY((matrix[4] * v.getX() + matrix[5] * v.getY() + matrix[6] * v.getZ() + matrix[7]) / w);
    retVector.setZ((matrix[8] * v.getX() + matrix[9] * v.getY() + matrix[10] * v.getZ() + matrix[11]) / w);

    return retVector;
}

Matrix4D & Matrix4D::operator = (const Matrix4D & mat)
{
    int i;

    for(i = 0; i < 16; i++)
        matrix[i] = mat.matrix[i];

    return (*this);
}

Matrix4D::Matrix4D(const Matrix4D & mat)
{
    int i;

    for(i = 0; i < 16; i++)
        matrix[i] = mat.matrix[i];
}


