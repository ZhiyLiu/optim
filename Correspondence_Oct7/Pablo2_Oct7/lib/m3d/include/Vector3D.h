/*---------------------------------------------------------------------------
Class definitions for Vector3D

Author:			Tom Fletcher
Version:		1.0
Last Modified:	6/06/99

Modified by Dibyendu on 06/03/2011 to include rotation operators
---------------------------------------------------------------------------*/

#ifndef VECTOR3D_H
#define VECTOR3D_H

#include <iostream>

#include "Mathdefs.h"

class Vector3D
{
public:
	// Constructor to create a zero vector
	Vector3D();

	// Constructor takes three coordinates
	Vector3D(double X, double Y, double Z);

	// Constructor takes a double vector
	Vector3D(double vect[3]);

	void print(std::ostream & out = std::cout) const;

	// Get and set functions
	void set(double X, double Y, double Z)
	{
		x = X;
		y = Y;
		z = Z;
	} 

	void set(double v[3])
	{
		x = v[0];
		y = v[1];
		z = v[2];
	}

	void set(const Vector3D & v)
	{
		x = v.x;
		y = v.y;
		z = v.z;
	}

	void setX(double X) {x = X;}
	void setY(double Y) {y = Y;}
	void setZ(double Z) {z = Z;}

	double getX() const {return x;}
	double getY() const {return y;}
	double getZ() const {return z;}

	// Arithmetic operations
	Vector3D& operator = (const Vector3D &v)
	{
		x = v.x;
		y = v.y;
		z = v.z;

		return (*this);
	}

	void operator += (const Vector3D &v)
	{
		x += v.x;
		y += v.y;
		z += v.z;
	}

	void operator += (const double v)
	{
		x += v;
		y += v;
		z += v;
	}

	void operator -= (const Vector3D &v)
	{
		x -= v.x;
		y -= v.y;
		z -= v.z;
	}

	void operator -= (const double v)
	{
		x -= v;
		y -= v;
		z -= v;
	}

	void operator *= (const double k)
	{
		x *= k;
		y *= k;
		z *= k;
	}

	void operator /= (const double k)
	{
		x /= k;
		y /= k;
		z /= k;
	}

	// Multiply by a 4x4 matrix
	void multByMatrix4(float * mat);

	// Compute the cross product with another vector
	Vector3D cross(const Vector3D &v) const;

	// replace this by element-wise minimum with 'a'
	void minPoint(const Vector3D &a);
	void maxPoint(const Vector3D &a);

	// Compute the element-wise product with another vector
	Vector3D vprod(const Vector3D &v) const
	{
		Vector3D ret;

		ret.x = x * v.x;
		ret.y = y * v.y;
		ret.z = z * v.z;

		return ret;
	}

	// Compute the element-wise division by another vector
        // Warning: no test for a zero divisior is made
	Vector3D vdiv(const Vector3D &v) const
	{
		Vector3D ret;

		ret.x = x / v.x;
		ret.y = y / v.y;
		ret.z = z / v.z;

		return ret;
	}

	// squared distance
	double distSquare(const Vector3D &v) const
	{
		return (x-v.x)*(x-v.x)+(y-v.y)*(y-v.y)+(z-v.z)*(z-v.z);
	}

    //  Convert the vector to its absolute value
    void abs();
	void round();

	// Compute the norm of the vector
	double norm() const;

	// Normalize the vector, returns old norm
	double normalize();

	// Dibyendu
	// Rotate the vector, using spherical coordinates, by (delTheta, delPhi)
	// theta = angle with z-axis, phi = angle (of x-y projection) with x-axis

	void rotateBy( double delTheta, double delPhi );

	// Declare our operations to be friends (so we can access private data)
	friend Vector3D operator + (const Vector3D &v1, const Vector3D &v2);
	friend Vector3D operator - (const Vector3D &v1, const Vector3D &v2);
	friend Vector3D operator - (const Vector3D &v);
	friend Vector3D operator * (const Vector3D &v, const double scalar);
	friend Vector3D operator * (const double scalar, const Vector3D &v);
    friend Vector3D operator / (const Vector3D &v, const double scalar);
	friend double operator * (const Vector3D &v1, const Vector3D &v2);
	friend bool operator == (const Vector3D &v1, const Vector3D &v2);
	friend bool operator != (const Vector3D &v1, const Vector3D &v2);

private:
	double x, y, z;
};

// Vector operations
Vector3D operator + (const Vector3D &v1, const Vector3D &v2);
Vector3D operator - (const Vector3D &v1, const Vector3D &v2);
Vector3D operator - (const Vector3D &v);
Vector3D operator * (const Vector3D &v, const double scalar);
Vector3D operator * (const double scalar, const Vector3D &v);
Vector3D operator / (const Vector3D &v, const double scalar);

// Dot product
double operator * (const Vector3D &v1, const Vector3D &v2);

bool operator == (const Vector3D &v1, const Vector3D &v2);
bool operator != (const Vector3D &v1, const Vector3D &v2);


inline std::ostream& operator<<(std::ostream& os, const Vector3D& v)
{ v.print(os); return os; }
#endif

