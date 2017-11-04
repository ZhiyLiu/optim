#ifndef M3D_BOUNDING_SPHERE_H
#define M3D_BOUNDING_SPHERE_H

#include "M3DObject.h"

// This class is used to compute an approximation to a bounding sphere
// for an m-rep object or figure.  The level of subdivision of the
// surface will determine how good the approximation is.

class M3DBoundingSphere
{
public:
    M3DBoundingSphere(const M3DObject & object, int subdivisions = 1);
    M3DBoundingSphere(const M3DFigure & figure, int subdivisions = 1);

    Vector3D & getCenter() { return center; }
    double getRadius() { return radius; }

private:
    Vector3D center;
    double radius;
};

#endif

