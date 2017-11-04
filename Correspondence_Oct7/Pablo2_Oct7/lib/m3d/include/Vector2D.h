#ifndef VECTOR2D_H
#define VECTOR2D_H

#include <iostream>

#include "Mathdefs.h"

class Vector2D
{
public:
    // Constructor to create a zero vector
    Vector2D();

    // Constructor takes three coordinates
    Vector2D(double X, double Y);

    // Constructor takes a double vector
    Vector2D(double vect[2]);

    void print(std::ostream & out = std::cout) const;

    // Get and set functions
    void set(double X, double Y)
    {
        x = X;
        y = Y;
    }

    void set(double v[2])
    {
        x = v[0];
        y = v[1];
    }

	void set(const Vector2D & v)
	{
		x = v.x;
		y = v.y;
	}

    void setX(double X) {x = X;}
    void setY(double Y) {y = Y;}

    double getX() const {return x;}
    double getY() const {return y;}

    // Arithmetic operations
    Vector2D& operator = (const Vector2D &v)
    {
        x = v.x;
        y = v.y;

        return (*this);
    }

    void operator += (const Vector2D &v)
    {
        x += v.x;
        y += v.y;
    }

    void operator -= (const Vector2D &v)
    {
        x -= v.x;
        y -= v.y;
    }

    void operator *= (const double k)
    {
        x *= k;
        y *= k;
    }

    void operator /= (const double k)
    {
        x /= k;
        y /= k;
    }

    // Compute the element-wise product with another vector
    Vector2D vprod(const Vector2D &v) const
    {
        Vector2D ret;

        ret.x = x * v.x;
        ret.y = y * v.y;

        return ret;
    }

    // Compute the element-wise division by another vector
    // Warning: no test for a zero divisior is made
    Vector2D vdiv(const Vector2D &v) const
    {
		Vector2D ret;

		ret.x = x / v.x;
		ret.y = y / v.y;

		return ret;
    }

    //  Convert the vector to its absolute value
    void abs();

    // Compute the norm of the vector
    double norm() const;

    // Normalize the vector, returns old norm
    double normalize();

    // Declare our operations to be friends (so we can access private data)
    friend Vector2D operator + (const Vector2D &v1, const Vector2D &v2);
    friend Vector2D operator - (const Vector2D &v1, const Vector2D &v2);
    friend Vector2D operator - (const Vector2D &v);
    friend Vector2D operator * (const Vector2D &v, const double scalar);
    friend Vector2D operator * (const double scalar, const Vector2D &v);
    friend Vector2D operator / (const Vector2D &v, const double scalar);
    friend double operator * (const Vector2D &v1, const Vector2D &v2);
    friend bool operator == (const Vector2D &v1, const Vector2D &v2);
    friend bool operator != (const Vector2D &v1, const Vector2D &v2);

private:
    double x, y, z;
};

// Vector operations
Vector2D operator + (const Vector2D &v1, const Vector2D &v2);
Vector2D operator - (const Vector2D &v1, const Vector2D &v2);
Vector2D operator - (const Vector2D &v);
Vector2D operator * (const Vector2D &v, const double scalar);
Vector2D operator * (const double scalar, const Vector2D &v);
Vector2D operator / (const Vector2D &v, const double scalar);

// Dot product
double operator * (const Vector2D &v1, const Vector2D &v2);

bool operator == (const Vector2D &v1, const Vector2D &v2);
bool operator != (const Vector2D &v1, const Vector2D &v2);


#endif

