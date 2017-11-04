/*
 * This file has been modified. The Geodesic class is
 * deprecated and the GeodesicSym class is implemented
 * in the framework of the primitive classes.
 *
 * This file provides utility functions for shape space
 * operations.
 *
 * @author Rohit Saboo
 */
#ifndef GEODESIC_H
#define GEODESIC_H

#include <math.h>
#include "Vector2D.h"
#include "Vector3D.h"
#include "VectorND.h"
#include "Quat.h"

namespace ShapeSpace {

/**
 * namespace RP1.
 * real projective space 1
 * Provides the log map and exp for RP1.
 */
namespace RP1 {
	/**
	 * Takes the log map of x.
	 *
	 * @param	x	an angle lying in (0,PI).
	 * @return	log map of x which lies in (-inf,+inf).
	 */
	inline double Log(const double x) {
		return tan(x-R_HALF_PI);
	}
	/**
	 * Returns the exp map of x.
	 *
	 * @param	x	the log map of an angle lying between (-inf,+inf).
	 * @return	an angle lying between (0,PI).
	 */
	inline double Exp(const double x) {
		return atan(x) + R_HALF_PI;
	}
}

/**
 * namespace S2.
 * Spherical space embedded in three dimensions.
 * Provides the log and exp map for S2 and some
 * helper functions.
 */
namespace S2 {
	Quat rotationFromOrigin(const Vector3D& n);
	Quat rotationToOrigin(const Vector3D& n);
	/**
	 * Returns the log map of n.
	 *
	 * @param	n	A point on a unit sphere.
	 * @return	A point on the tangent space.
	 */
	Vector2D Log(const Vector3D& n);
	/**
	 * Returns the exponential map of v.
	 *
	 * @param	v	A point in tangent space
	 * @return	A point on the unit Sphere
	 */
	Vector3D Exp(const Vector2D & v);

	/**
	 * Computes the average of several
	 * different points on a unit sphere.
	 * The points should *be* on a unit sphere.
	 * They are not normalized.
	 *
	 * @param	points	An array of points on the unit sphere.
	 * @param	mNum	Number of points.
	 * @param	weights	The weights for the average.
	 *					This cannot be NULL.
	 * @return	A normalized Vector3D instance representing
	 *			the mean of the points passed in.
	 */
	 Vector3D mean(const Vector3D* points, const int mNum,
	 	const double* weights);
}


};


#ifndef BINARY
#include "M3DPrimitive.h"

/**
 * @deprecated
 * Geodesic is an ancient class which does statistics in lie space. It has
 * been superseded by M3DPGA suit of classes. @see M3DPGA.
 */
class Geodesic
{

public:
	// By default, ignorePrimType is true, which means all the atoms are
	// treated as standard ones.
	// To include end atoms, set ignorePrimType as false when instantializing
	// For certain operantions, all the atoms involved should be in the same
	// type, otherwise all will be considered as standard atoms
	Geodesic(bool _ignorePrimType = true)
	{
		ignorePrimType=_ignorePrimType;
	}

	// Exponential and Log map for an atom
	bool atomExpMap(VectorND &v, M3DPrimitive *mNew);
	bool atomLogMap(M3DPrimitive *m, VectorND *v);

	// Inverse of an atom
	bool atomInverse(M3DPrimitive* mOld, M3DPrimitive *mNew);

	// Compose two atoms by the group operantions
	bool atomCompose(M3DPrimitive* m1, M3DPrimitive* m2, M3DPrimitive *mNew);

	// Average of an array of atoms
	bool atomAverage(int mNum, M3DPrimitive **mList, M3DPrimitive *mNew);

	// Interpolation between atoms m1/m2 by t in [0, 1]
	bool atomInterp(double t, M3DPrimitive *m1, M3DPrimitive *m2, M3DPrimitive *mNew);

	// Geodesic distance between atoms m1/m2
	double atomDistance(M3DPrimitive* m1, M3DPrimitive* m2);
	double atomSquareDistance(M3DPrimitive* m1, M3DPrimitive* m2);

private:

	bool ignorePrimType;
};

#endif


#endif
