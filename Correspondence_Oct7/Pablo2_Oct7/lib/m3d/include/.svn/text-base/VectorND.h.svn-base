/*---------------------------------------------------------------------------
Class definitions for VectorND

Author:			Qiong Han
Last Modified:	10/06/03
---------------------------------------------------------------------------*/

#ifndef VectorND_H
#define VectorND_H

#include <iostream>

#include "Mathdefs.h"

class VectorND
{
public:

	// Constructors
	VectorND();
	VectorND(int sz);
	VectorND(const VectorND & v);

	~VectorND();

	void print(std::ostream & out = std::cout) const;

	// Get and set functions
	void setSize(int sz) {
		if (size > 0)
			delete [] vec;
		size = sz;
		vec = new double[size];
	}

	void set(double * v) {
		for (int i = 0; i < size; i++)
			vec[i] = v[i];
	}

	void set(const VectorND & v) {
		set(v.vec);
	}

	void set(int dex, double E) {
		if(size > 0 && dex >= 0 && dex < size)
			vec[dex] = E;
	}

	void set(int startDex, const Vector3D & v) {
		if (size > 0 && startDex >= 0 && startDex + 3 <= size) {
			vec[startDex] = v.getX();
			vec[startDex+1] = v.getY();
			vec[startDex+2] = v.getZ();
		}
	}

	void set(int startDex, VectorND & v) {
		if (size > 0 && startDex >= 0 && startDex + v.size <= size) {
			for (int i = 0; i < v.size; i++)
				vec[startDex+i] = v.get(i);
		}
	}

	void set(double c) {
		for (int i = 0; i < size; i++)
			vec[i] = c;
	}

	double * get() const { return vec; }
	double get(int dex) const { return vec[dex]; }

	int getSize() const { return size; }

	// Arithmetic operations
	// These assume that the two vectors are of equal length

	VectorND & operator = (const VectorND & v) {
		for (int i = 0; i < size; i++)
			vec[i] = v.get(i);
		return (*this);
	}

	void operator += (const VectorND &v) {
		for (int i = 0; i < size; i++)
			vec[i] += v.get(i);
	}

	void operator -= (const VectorND &v) {
		for (int i = 0; i < size; i++)
			vec[i] -= v.get(i);
	}

	void operator *= (const double k) {
		for (int i = 0; i < size; i++)
			vec[i] *= k;
	}

	void operator /= (const double k) {
		for (int i = 0; i < size; i++)
			vec[i] /= k;
	}

	// Multiply by a NxN matrix
	void multByMatrix(float * mat);

	// Compute the cross product with another vector
	// Currently not implemented
	//VectorND cross(const VectorND &v) const;

	// Compute the element-wise product with another vector
	VectorND vprod(VectorND & v) const {
		VectorND ret(size);

		for (int i = 0; i < size; i++)
			ret.set(i, vec[i]*v.get(i));

		return ret;
	}

	// Compute the element-wise division by another vector
	// Warning: no test for a zero divisior is made
	VectorND vdiv(VectorND & v) const {
		VectorND ret(size);

		for (int i = 0; i < size; i++)
			ret.set(i, vec[i]/v.get(i));

		return ret;
	}

    //  Convert the vector to its absolute value
    void abs();

	// Compute the norm of the vector
	double norm() const;
	double normSquare() const;

	// Normalize the vector, returns old norm
	double normalize();

	// Declare our operations to be friends (so we can access private data)
	friend VectorND operator + (const VectorND &v1, const VectorND &v2);
	friend VectorND operator - (const VectorND &v1, const VectorND &v2);
	friend VectorND operator - (const VectorND &v);
	friend VectorND operator * (const VectorND &v, const double scalar);
	friend VectorND operator * (const double scalar, const VectorND &v);
    friend VectorND operator / (const VectorND &v, const double scalar);
	friend double operator * (const VectorND &v1, const VectorND &v2);
	friend bool operator == (const VectorND &v1, const VectorND &v2);
	friend bool operator != (const VectorND &v1, const VectorND &v2);

private:
	int size;
	double * vec;
};

// Vector operations
VectorND operator + (const VectorND &v1, const VectorND &v2);
VectorND operator - (const VectorND &v1, const VectorND &v2);
VectorND operator - (const VectorND &v);
VectorND operator * (const VectorND &v, const double scalar);
VectorND operator * (const double scalar, const VectorND &v);
VectorND operator / (const VectorND &v, const double scalar);

// Dot product
double operator * (const VectorND &v1, const VectorND &v2);

bool operator == (const VectorND &v1, const VectorND &v2);
bool operator != (const VectorND &v1, const VectorND &v2);


#endif

