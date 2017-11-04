#include <math.h>

#include "Vector2D.h"

using namespace std;

Vector2D::Vector2D()
{
    x = 0;
    y = 0;
}

Vector2D::Vector2D(double X, double Y)
{
    x = X;
    y = Y;
}

Vector2D::Vector2D(double vect[3])
{
    x = vect[0];
    y = vect[1];
}

void Vector2D::print(ostream & out) const
{
    out << '<' << x << ", " << y << ">\n";
}

//  Convert the vector to its absolute value
void Vector2D::abs()
{
	if (x < 0.0) x = -x;
	if (y < 0.0) y = -y;
}

// Computes the norm (square of the dot product of the vector w/ itself)
double Vector2D::norm() const
{
    return sqrt((*this) * (*this));
}

// Sets vector to a unit length, and returns old norm
double Vector2D::normalize()
{
    double length = norm();

    if(length == 0.0)
        return 0.0;

    x /= length;
    y /= length;

    return length;
}

// Vector add
Vector2D operator + (const Vector2D &v1, const Vector2D &v2)
{
    Vector2D vect;

    vect.x = v1.x + v2.x;
    vect.y = v1.y + v2.y;

    return vect;
}

// Vector subtract
Vector2D operator - (const Vector2D &v1, const Vector2D &v2)
{
    Vector2D vect;

    vect.x = v1.x - v2.x;
    vect.y = v1.y - v2.y;

    return vect;
}

// Vector negate
Vector2D operator - (const Vector2D &v)
{
    Vector2D vect;

    vect.x = -(v.x);
    vect.y = -(v.y);

    return vect;
}

// Vector scale
Vector2D operator * (const Vector2D &v, const double scalar)
{
    Vector2D vect;

    vect.x = scalar * v.x;
    vect.y = scalar * v.y;

    return vect;
}

// Another vector scale
Vector2D operator * (const double scalar, const Vector2D &v)
{
    Vector2D vect;

    vect.x = scalar * v.x;
    vect.y = scalar * v.y;

    return vect;
}

Vector2D operator / (const Vector2D &v, const double scalar)
{
    Vector2D vect;

    vect.x = v.x / scalar;
    vect.y = v.y / scalar;

    return vect;
}

// Dot product
double operator * (const Vector2D &v1, const Vector2D &v2)
{
    return ((v1.x * v2.x) + (v1.y * v2.y));
}

bool operator == (const Vector2D &v1, const Vector2D &v2)
{
    if(v1.x != v2.x)
        return false;
    if(v1.y != v2.y)
        return false;

    return true;
}

bool operator != (const Vector2D &v1, const Vector2D &v2)
{
    return !(v1 == v2);
}


