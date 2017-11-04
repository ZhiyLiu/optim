/*---------------------------------------------------------------------------
Member functions for class VectorND

Author:			Qiong Han
Last Modified:	10/06/03
--------------------------------------------------------*/
#include <math.h>

#include "Vector3D.h"
#include "VectorND.h"

#include <string.h>

using namespace std;

VectorND::VectorND()
{
	size=0;
	vec=NULL;
}

VectorND::VectorND(int sz)
{
	size = sz;
	vec = new double[size];
}

VectorND::~VectorND()
{
	if (size > 0)
		delete [] vec;
}

VectorND::VectorND(const VectorND & v)
{
	size = v.size;
	vec = new double[size];
	memcpy(vec, v.vec, sizeof(double)*size);
}

void VectorND::print(ostream & out) const
{
	int i;
	cout << "<";
	for(i=0; i<size-1; i++)
		out << vec[i] << ", ";
	cout << vec[size-1] << ">\n";
}

void VectorND::multByMatrix(float * mat)
{
	//VectorND newV(size);
	double *newX=new double[size];
	int i,
		j,
		dex;

	for(i=0; i<size; i++)
	{
		double tempX=0;
		dex=i;
		for(j=0; j<size; j++)
		{
			tempX+=*(mat+dex)*vec[j];
			dex+=size+1;
		}
		dex+=size+1;
		tempX+=*(mat+dex);
		newX[i]=tempX;
	}

	/*
	reference from Vector3D
	newX = *(mat    ) * x + *(mat + 4) * y + *(mat + 8) * z + *(mat + 12);
	newY = *(mat + 1) * x + *(mat + 5) * y + *(mat + 9) * z + *(mat + 13);
	newZ = *(mat + 2) * x + *(mat + 6) * y + *(mat + 10) * z + *(mat + 14);
	*/

	set(newX);
	delete []newX;
}

// Computes the cross product with another vector v
// Currently not implemented
/*
VectorND VectorND::cross(const VectorND &v) const
{
	VectorND vect;

	vect.x = y * v.z - z * v.y;
	vect.y = z * v.x - x * v.z;
	vect.z = x * v.y - y * v.x;

	return vect;
}
*/

//  Convert the vector to its absolute value
void VectorND::abs()
{
	for (int i = 0; i < size; i++)
		if (vec[i] < 0.0) vec[i] = -vec[i];
}

// Computes the norm (square of the dot product of the vector w/ itself)
double VectorND::norm() const
{
	return sqrt((*this) * (*this));
}

double VectorND::normSquare() const
{
	return (*this) * (*this);
}

// Sets vector to a unit length, and returns old norm
double VectorND::normalize()
{
	double length = norm();

	if(length == 0.0)
		return 0.0;

	int i;
	for(i=0; i<size; i++)
		vec[i]/=length;

	return length;
}

// Vector add
VectorND operator + (const VectorND &v1, const VectorND &v2)
{
	VectorND vect;

	if(v1.size==v2.size)
	{
		int i;
		vect.setSize(v1.size);
		for(i=0; i<vect.getSize(); i++)
			vect.set(i, v1.vec[i]+v2.vec[i]); //v1.get(i)+v2.get(i));
	}
	return vect;
}

// Vector subtract
VectorND operator - (const VectorND &v1, const VectorND &v2)
{
	VectorND vect;

	if(v1.size==v2.size)
	{
		int i;
		vect.setSize(v1.size);
		for(i=0; i<vect.getSize(); i++)
			vect.set(i, v1.vec[i]-v2.vec[i]);	//get(i)-v2.get(i));
	}

	return vect;
}

// Vector negate
VectorND operator - (const VectorND &v)
{
	VectorND vect;

	int i;
	vect.setSize(v.size);
	for(i=0; i<vect.size; i++)
		vect.set(i, -v.vec[i]);

	return vect;
}

// Vector scale
VectorND operator * (const VectorND &v, const double scalar)
{
	VectorND vect;

	int i;
	vect.setSize(v.size);
	for(i=0; i<vect.size; i++)
		vect.set(i, scalar*v.vec[i]);

	return vect;
}

// Another vector scale
VectorND operator * (const double scalar, const VectorND &v)
{
	VectorND vect;

	int i;
	vect.setSize(v.size);
	for(i=0; i<vect.size; i++)
		vect.set(i, scalar*v.vec[i]);

	return vect;
}

VectorND operator / (const VectorND &v, const double scalar)
{
    VectorND vect;

	int i;
	vect.setSize(v.size);
	for(i=0; i<vect.size; i++)
		vect.set(i, scalar/v.vec[i]);

    return vect;
}

// Dot product
double operator * (const VectorND &v1, const VectorND &v2)
{
	int i;
	if(v1.size!=v2.size)
		return 0;

	double sum=0;
	for(i=0; i<v1.size; i++)
		sum+=v1.vec[i]*v2.vec[i];
	return sum;
}

bool operator == (const VectorND &v1, const VectorND &v2)
{
	int i;
	VectorND vect;

	if(v1.size!=v2.size)
		return false;
	else
	{
		vect.setSize(v1.size);
		for(i=0; i<vect.size; i++)
			if(v1.vec[i]!=v2.vec[i])
				return false;
		return true;
	}
}

bool operator != (const VectorND &v1, const VectorND &v2)
{
	return !(v1 == v2);
}


