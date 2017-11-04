/*---------------------------------------------------------------------------
Class definitions for Quat

Author:			Tom Fletcher
Version:		1.0
Last Modified:	6/08/99
---------------------------------------------------------------------------*/

#ifndef QUAT_H
#define QUAT_H

#include "Vector3D.h"


class Quat
{
public:
	// Creates a quaternion with no rotation (<0, 0, 0>, 1)
	Quat();

	// Creates the quaternion (<X, Y, Z>, W)
	Quat(double X, double Y, double Z, double W);

	// Creates the quaternion (V, W)
	Quat(const Vector3D & V, double W);

	// Creates a quaternion from a unit frame
	// (represents rotation of frame from standard frame)
	Quat(const Vector3D & X, const Vector3D & Y, const Vector3D & Z);

	// Copy constructor
	Quat(const Quat & q);

	void print(std::ostream & out = std::cout) const;

	// get and set functions
	double getX() const {return x;}
	double getY() const {return y;}
	double getZ() const {return z;}
	double getW() const {return w;}

	Vector3D getVector() const {return Vector3D(x, y, z);}

	void set(double X, double Y, double Z, double W)
	{
		x = X;
		y = Y;
		z = Z;
		w = W;
	}

	void setX(double X) {x = X;}
	void setY(double Y) {y = Y;}
	void setZ(double Z) {z = Z;}
	void setW(double W) {w = W;}

	// Sets the quaternion as rotation of an angle about an axis
	void setAxisAngle(Vector3D axis, double angle);
	void getAxisAngle(Vector3D & axis, double & angle);

	// Assignment operator
	Quat & operator = (const Quat & q)
	{
		x = q.x;
		y = q.y;
		z = q.z;
		w = q.w;

		return (*this);
	}

	// Compute the quaternion's norm (NOT the same as length)
	double norm() const;

	// Compute the quaternion's length (sqrt of norm)
	double length();

	// Set to unit length
	void normalize();

	// Returns the conjugate quaternion
	Quat conj() const;

	// Rotates the vector using the quaternion
	void rotateVector(Vector3D & v) const;

	// Builds a rotation matrix from the current quaternion
	void buildRotMatrix(double * m) const;

	// Builds a quaternion from a rotation matrix
	void matrixToQuat(double m[4][4]);
	void matrixToQuat(double * m);

	/** 
	 * Build a quaternion from a set of frame vectors and saves the result
	 * in itself.
	 *
	 * @param	x	The rotated x axis of the new frame
	 * @param	y	The rotated y axis of the new frame
	 * @return	none
	 */
	void buildFromFrame( Vector3D x, Vector3D y);

	// Dot product
	double dot(const Quat & q) const;

	// Spherical linear interpolation
	Quat slerp(double t, const Quat & q0, const Quat & q1);

	// Exponential and log map for quaternion
	Quat expMap(Vector3D v);
	Vector3D logMap();
	void neg();

	// Declare operations to be friends (so we can access private data)
	friend Quat operator + (const Quat & q1, const Quat & q2);
	friend Quat operator - (const Quat & q1, const Quat & q2);
    friend Quat operator - (const Quat & q);
	friend Quat operator * (const Quat & q1, const Quat & q2);
	friend Quat operator * (const Quat & q, double scalar);
	friend Quat operator * (double scalar, const Quat & q);

private:
	double x, y, z, w;
};

// Quaternion operations
Quat operator + (const Quat & q1, const Quat & q2);
Quat operator - (const Quat & q1, const Quat & q2);
Quat operator - (const Quat & q);
Quat operator * (const Quat & q1, const Quat & q2);
Quat operator * (const Quat & q, double scalar);
Quat operator * (double scalar, const Quat & q);
inline std::ostream& operator<<(std::ostream& os, const Quat& q)
{ q.print(os); return os; }

#endif

