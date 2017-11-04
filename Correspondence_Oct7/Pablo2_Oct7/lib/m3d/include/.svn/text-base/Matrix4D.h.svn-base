#ifndef MATRIX_4D_H
#define MATRIX_4D_H


#include "Vector3D.h"
#include "Quat.h"



class Matrix4D
{
public:
    Matrix4D();
    Matrix4D(double * mat);
    Matrix4D(float * mat);

    void setData(double * mat);
    void setData(float * mat);
    const double * getData() { return &(matrix[0]); }

    double getElement(int i) const { return matrix[i]; }

    void print(std::ostream & out = std::cout) const;

	// Functions to change the matrix
    void translate(const Vector3D & v);
    void rotate(const Quat & q);
    void scale(double s);
    void setToIdentity();
	void transpose();

	// This assumes the matrix is a similarity transform.  If the returned
	// rotation, scale, and translation are applied to an identity matrix
	// in that order, they will produce the original composed matrix.
	void decompose(Quat & rotation, double & scale, Vector3D & translation) const;

    void multiply(double * mat);
    void multiply(const Matrix4D & mat);

	// Apply the matrix to a vector
    Vector3D multiplyVector(const Vector3D & v) const;

	// Determinant; only correct if the matrix is a similarity transform
	double det() const;

    Matrix4D & operator = (const Matrix4D & mat);
    Matrix4D(const Matrix4D & mat);

private:
    double matrix[16];
};


#endif

