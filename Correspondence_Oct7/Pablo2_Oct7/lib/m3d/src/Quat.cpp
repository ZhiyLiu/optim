/*---------------------------------------------------------------------------
Member functions for class Quat

Author:         Tom Fletcher
Version:        1.0
Last Modified:  6/08/99
---------------------------------------------------------------------------*/
#include <math.h>

#include "Quat.h"

using namespace std;

const double epsilon = 1.0e-03;
#ifndef EPSILON
#define EPSILON 1e-8
#endif 

// Creates a quaternion with no rotation (<0, 0, 0>, 1)
Quat::Quat()
{
    x = 0;
    y = 0;
    z = 0;
    w = 1;
}

// Creates the quaternion (<X, Y, Z>, W)
Quat::Quat(double X, double Y, double Z, double W)
{
    x = X;
    y = Y;
    z = Z;
    w = W;
}

// Creates the quaternion (V, W)
Quat::Quat(const Vector3D & V, double W)
{
    x = V.getX();
    y = V.getY();
    z = V.getZ();
    w = W;
}

// Creates a quaternion from a unit frame
// (represents rotation of frame from standard frame)
Quat::Quat(const Vector3D & X, const Vector3D & Y, const Vector3D & Z)
{
    double rotMat[4][4];

#if 1
    rotMat[0][0] = X.getX();
    rotMat[0][1] = X.getY();
    rotMat[0][2] = X.getZ();
    rotMat[1][0] = Y.getX();
    rotMat[1][1] = Y.getY();
    rotMat[1][2] = Y.getZ();
    rotMat[2][0] = Z.getX();
    rotMat[2][1] = Z.getY();
    rotMat[2][2] = Z.getZ();
#else
    rotMat[0][0] = X.getX();
    rotMat[1][0] = X.getY();
    rotMat[2][0] = X.getZ();
    rotMat[0][1] = Y.getX();
    rotMat[1][1] = Y.getY();
    rotMat[2][1] = Y.getZ();
    rotMat[0][2] = Z.getX();
    rotMat[1][2] = Z.getY();
    rotMat[2][2] = Z.getZ();
#endif

    matrixToQuat(rotMat);
}

// Copy constructor
Quat::Quat(const Quat & q)
{
    x = q.x;
    y = q.y;
    z = q.z;
    w = q.w;
}

void Quat::print(ostream & out) const
{
    out << "<" << x << ", " << y << ", " << z << ", " << w << ">" << endl;
}

// Sets the quaternion as rotation of an angle about an axis
void Quat::setAxisAngle(Vector3D axis, double angle)
{
    double sinHalfAngle = sin(angle / 2.0);

    axis.normalize();

    x = axis.getX() * sinHalfAngle;
    y = axis.getY() * sinHalfAngle;
    z = axis.getZ() * sinHalfAngle;

    w = cos(angle / 2.0);
}

// Gets the quaternion as rotation of an angle about an axis
void Quat::getAxisAngle(Vector3D & axis, double & angle)
{
    axis.set(x, y, z);
    axis.normalize();
    angle = 2.0*acos(w);
}

// Compute the quaternion's norm (NOT the same as length; 
// this is different from usual mathematical terminology)
// Who picked this terminology?  Mathematically, the norm 
// is the same as the length.  This is the squared norm.
double Quat::norm() const
{
    return ((x * x) + (y * y) + (z * z) + (w * w));
}

// Compute the quaternion's length (sqrt of norm)
double Quat::length()
{
    return sqrt(norm());
}

// Set to unit length
void Quat::normalize()
{
    double len = length();

//  assert(len != 0);

    x /= len;
    y /= len;
    z /= len;
    w /= len;
}

// Returns the conjugate quaternion
Quat Quat::conj() const
{
    Quat qRet;

    double n = norm();

    qRet.x = -x / n;
    qRet.y = -y / n;
    qRet.z = -z / n;
    qRet.w = w / n;

    return qRet;
}

static int errorCount = 10;

// Rotates the vector using the quaternion
void Quat::rotateVector(Vector3D & v) const
{
    Quat qVect;

#ifdef _DEBUG
    // Make sure we have a unit quaternion
    double len = x*x + y*y + z*z + w*w;
    if (len >= 1.0 + epsilon || len <= 1.0 - epsilon)
    {
        if (--errorCount > 0)
        {
            cerr << "Attempting rotation of vector by non-unit "
                 << "quaternion (squared length=" << len << ") : ";
            this->print();
        }
        return;
    }
#endif

    qVect = (*this) * Quat(v, 0) * (*this).conj();

    v = qVect.getVector();
}

// Builds a rotation matrix
void Quat::buildRotMatrix(double * m) const
{
    // column 0
    *(  m) = 1.0 - 2.0 * (y * y + z * z);
    *(++m) = 2.0 * (x * y - z * w);
    *(++m) = 2.0 * (z * x + y * w);
    *(++m) = 0.0;

    // column 1
    *(++m) = 2.0 * (x * y + z * w);
    *(++m) = 1.0 - 2.0 * (z * z + x * x);
    *(++m) = 2.0 * (y * z - x * w);
    *(++m) = 0.0;

    // column 2
    *(++m) = 2.0 * (z * x - y * w);
    *(++m) = 2.0 * (y * z + x * w);
    *(++m) = 1.0 - 2.0 * (y * y + x * x);
    *(++m) = 0.0;

    *(++m) = 0.0;
    *(++m) = 0.0;
    *(++m) = 0.0;
    *(++m) = 1.0;
}

// Builds a quaternion from a rotation matrix
void Quat::matrixToQuat(double m[4][4])
{
    double tr, s, q[4];
    int i, j, k;

    int nxt[3] = { 1, 2, 0 };

    tr = m[0][0] + m[1][1] + m[2][2];

    if (1.0 + tr > 0.0) {
        s = sqrt(tr + 1.0);
        w = s / 2.0;
        s = 0.5 / s;
        x = (m[1][2] - m[2][1]) * s;
        y = (m[2][0] - m[0][2]) * s;
        z = (m[0][1] - m[1][0]) * s;
    }
    else {
        i = 0;
        if(m[1][1] > m[0][0])
            i = 1;

        if(m[2][2] > m[i][i])
            i = 2;

        j = nxt[i];
        k = nxt[j];

        s = sqrt((m[i][i] - (m[j][j] + m[k][k])) + 1.0);

        q[i] = s * 0.5;

        if (s != 0.0)
            s = 0.5 / s;

        q[3] = (m[j][k] - m[k][j]) * s;
        if( q[3] < -1e-5 || q[3] > 1e5 ) {
            cout << "ERROR: q[3] is non-zero, check the code!!! q[3] = " << q[3] << " "
                 << __FILE__ << ":" << __LINE__ << "\n";
        }
        q[j] = (m[i][j] + m[j][i]) * s;
        q[k] = (m[i][k] + m[k][i]) * s;

        x = q[0];
        y = q[1];
        z = q[2];
        w = q[3];
    }

}

// Builds a quaternion from a rotation matrix
void Quat::matrixToQuat(double * m)
{
    double tr, s, q[4];
    int i, j, k;

    int nxt[3] = { 1, 2, 0 };

    tr = m[0] + m[5] + m[10];

    if (1 + tr > 0.0) {
        s = sqrt(tr + 1.0);
        w = s / 2.0;
        s = 0.5 / s;
        x = (m[6] - m[9]) * s;
        y = (m[8] - m[2]) * s;
        z = (m[1] - m[4]) * s;
    }
    else {
        i = 0;
        if(m[5] > m[0])
            i = 1;

        if(m[10] > m[i + 4*i])  // m(i, j) is m[i + 4*j]
            i = 2;

        j = nxt[i];
        k = nxt[j];

        s = sqrt((m[i + 4*i] - (m[j + 4*j] + m[k + 4*k])) + 1.0);

        q[i] = s * 0.5;

        if (s != 0.0)
            s = 0.5 / s;

        q[3] = (m[j + 4*k] - m[k + 4*j]) * s;
        if( q[3] < -1e-5 || q[3] > 1e5 ) {
            cout << "ERROR: q[3] is non-zero, check the code!!! q[3] = " << q[3] << " "
                 << __FILE__ << ":" << __LINE__ << "\n";
        }
        q[j] = (m[i + 4*j] + m[j + 4*i]) * s;
        q[k] = (m[i + 4*k] + m[k + 4*i]) * s;

        x = q[0];
        y = q[1];
        z = q[2];
        w = q[3];
    }
}

/** 
 * Builds a quaternion from a set of frame vectors and saves the result
 * in itself. This code is a translated copy of the matlab function
 * quatFromFrame
 *
 * @param   x   The rotated x axis of the new frame
 * @param   y   The rotated y axis of the new frame
 * @return  none
 */
void Quat::buildFromFrame(Vector3D x, Vector3D y)
{
    x.normalize();
    y.normalize();
    const Vector3D z    = x.cross(y);

    const double T  = x.getX() + y.getY() + z.getZ();
    double cos_phi2;
    double qx,qy,qz;
    double S;
    // @Jeong
    // NOTE: T+1 == 0 when rotation by 180 degree == flip 
    // In converting rotation matrix to quaternion,
    // we need to handle this case specially.   
    if ( (1+T) > EPSILON ) {
        cos_phi2 = sqrt(T+1)/2;
        qx  = (y.getZ() - z.getY()) / (4*cos_phi2);
        qy  = (z.getX() - x.getZ()) / (4*cos_phi2);
        qz  = (x.getY() - y.getX()) / (4*cos_phi2);
    }
    else {  // Flip
        cos_phi2 = 0.0;
        if( 1 + x.getX() > 0.0 ) {
            qx  = sqrt( 0.5*(1+x.getX()));
            S   = 0.5/qx;
            qy  = y.getX() * S;
            qz  = z.getX() * S;
        }
        else if( 1 + y.getY() > 0.0 ) {
            qy  = sqrt( 0.5*(1+y.getY()));
            S   = 0.5/qy;
            qx  = y.getY() * S;
            qz  = z.getY() * S;
        }
        else if( 1 + z.getZ() > 0.0 ) {
            qz  = sqrt( 0.5*(1+z.getZ()));
            S   = 0.5/qz;
            qx  = z.getX() * S;
            qz  = z.getY() * S;
        }
        else {
            cerr << "Quat::buildFromFrame(...) : This error should not happen" << endl;
            cos_phi2    = 0.0;
            qx  = 0.0;
            qy  = 0.0;
            qz  = 0.0;
            x.print(cerr);
            y.print(cerr);
            cerr << flush;
        }
    }
    this->x = qx;
    this->y = qy;
    this->z = qz;
    this->w = cos_phi2;
    normalize();
}



// Dot product
double Quat::dot(const Quat & q) const
{
    return (x * q.x) + (y * q.y) + (z * q.z) + (w * q.w);
}

// Spherical Linear Interpolation
Quat Quat::slerp(double t, const Quat & q0, const Quat & q1)
{
    // cosine of angle between quats
    double cosAngle = fabs(q0.dot(q1));

    // Clamp between [-1, 1]
    if(cosAngle < -1.0)
        cosAngle = -1.0;
    else if(cosAngle > 1.0)
        cosAngle = 1.0;

    double angle = acos(cosAngle);
    double sinAngle = sin(angle);
    if(sinAngle < epsilon)
        return q0;

    double sinInv = 1.0 / sinAngle;

    return ((sin((1 - t) * angle) * q0) + (sin(t * angle) * q1)) * sinInv;
}

// Quaternion addition
Quat operator + (const Quat & q1, const Quat & q2)
{
    return Quat(q1.x + q2.x, q1.y + q2.y, q1.z + q2.z, q1.w + q2.w);
}

// Quaternion subtraction
Quat operator - (const Quat & q1, const Quat & q2)
{
    return Quat(q1.x - q2.x, q1.y - q2.y, q1.z - q2.z, q1.w - q2.w);
}

// Quaternion negation
Quat operator - (const Quat & q)
{
    return Quat(-q.x, -q.y, -q.z, -q.w);
}

// Quaternion multiplication
Quat operator * (const Quat & q1, const Quat & q2)
{
    Quat qRet;

    qRet.x = (q1.x * q2.w) + (q1.w * q2.x) + (q1.y * q2.z) - (q1.z * q2.y);
    qRet.y = (q1.y * q2.w) + (q1.w * q2.y) + (q1.z * q2.x) - (q1.x * q2.z);
    qRet.z = (q1.z * q2.w) + (q1.w * q2.z) + (q1.x * q2.y) - (q1.y * q2.x);
    qRet.w = (q1.w * q2.w) - (q1.x * q2.x) - (q1.y * q2.y) - (q1.z * q2.z);

    return qRet;
}

// Scales a quaternion
Quat operator * (const Quat & q, double scalar)
{
    Quat qRet;

    qRet.x = q.x * scalar;
    qRet.y = q.y * scalar;
    qRet.z = q.z * scalar;
    qRet.w = q.w * scalar;

    return qRet;
}

// Scales a quaternion
Quat operator * (double scalar, const Quat & q)
{
    Quat qRet;

    qRet.x = q.x * scalar;
    qRet.y = q.y * scalar;
    qRet.z = q.z * scalar;
    qRet.w = q.w * scalar;

    return qRet;
}

// Exponential and log map for quoternion
Quat Quat::expMap(Vector3D v)
{
    Quat qRet;
    double halfTheta = v.norm()/2;
    if (halfTheta == 0) {
        qRet.setW(1);
        qRet.setX(0);
        qRet.setY(0);
        qRet.setZ(0);
    }
    else {
        v.normalize();
        v = v*sin(halfTheta);
        qRet.setW(cos(halfTheta));
        qRet.setX(v.getX());
        qRet.setY(v.getY());
        qRet.setZ(v.getZ());
    }
    return qRet;
}

Vector3D Quat::logMap()
{
    Vector3D vRet;
    Quat qNew = *this;
    double halfTheta;

    if (norm() != 1)
        qNew.normalize();

    halfTheta=acos(qNew.getW());

    if (halfTheta == 0) {
        vRet.set(0, 0, 0);
    }
    else {
        vRet.set(qNew.getX(), qNew.getY(), qNew.getZ());
        vRet = vRet*(2.0*halfTheta/sin(halfTheta));
    }
    return vRet;
}

void Quat::neg()
{
    w =- w;
    x =- x;
    y =- y;
    z =- z;
}



