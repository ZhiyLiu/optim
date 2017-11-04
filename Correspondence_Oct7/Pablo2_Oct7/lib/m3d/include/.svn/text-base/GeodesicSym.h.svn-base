/*
 * This file is deprecated, do not use for any new development.
 *  - rrs
 */
#ifndef GEODESICSYM_H
#define GEODESICSYM_H

/*
    Class definitions for Geodesic

    Author:  Qiong Han
    Last Modified:  10/07/03
*/

#include <math.h>
#include "Vector2D.h"
#include "Vector3D.h"
#include "VectorND.h"
#include "Quat.h"


#ifdef BINARY
// #define OLD_BINARY_PABLO
#endif


class M3DSymPrimitive {

public:

	M3DSymPrimitive() {
		elongation = 1.0;
	}
	M3DSymPrimitive(Vector3D _x, double _r, Vector3D _n0, Vector3D _n1,
		double _elongation)
	{
		x = _x;
		r = _r;
		n0 = _n0;
		n1 = _n1;
		elongation = 1.0;
	}
	M3DSymPrimitive(Vector3D _x, double _r, Vector3D _n0, Vector3D _n1) {
		x = _x;
		r = _r;
		n0 = _n0;
		n1 = _n1;
		elongation = 1.0;
	}

    Vector3D x;         // Primitive medial site position
    double r;           // Primitive radius
	Vector3D n0, n1;
	double elongation;
};


/**
 * @deprecated
 * Geodesic is an ancient class which does statistics in lie space. It has
 * been superseded by M3DPGA suit of classes. @see M3DPGA.
 */
class GeodesicSym
{

public:

	// By default, ignore is true, which means all the atoms are treated
	// as standard ones.  To include end atoms, set ignore to false when
	// instantiating.  For certain operations, all the atoms involved
	// should be the same type, otherwise all will be considered to be
	// standard atoms.
	GeodesicSym(bool ignore = true) {
		ignorePrimType = ignore;
	}

	// Compose the difference mean to the mean M-rep
	bool atomCrossingSpaceCompose(M3DPrimitive * prim, M3DPrimitive * deltaPrim);

	// Align atom according to its neighboring atoms 
	bool atomsAlignment(int nPrims, M3DPrimitive ** mSet1, M3DPrimitive ** mSet2,
		Vector3D & tran, Quat & q, double & s, bool alignAll = false);

	// Exponential and Log map for an atom
	bool atomExpMap(VectorND & v, M3DPrimitive * mNew);
	bool atomLogMap(M3DPrimitive * m, VectorND * v);

	// Exponential and Log map for a sym atom
	bool symAtomExpMap(VectorND & v, M3DSymPrimitive * mNew);
	bool symAtomLogMap(M3DSymPrimitive * m, VectorND * v);
	Vector2D sphereLogMap(Vector3D n);
	Vector3D sphereExpMap(Vector2D v);

	Quat rotationFromOrigin(Vector3D n);
	Quat rotationToOrigin(Vector3D n);

	// Inverse of an atom
	bool atomInverse(M3DPrimitive * mOld, M3DPrimitive * mNew);

	// Compose two atoms by the group operations
	bool atomCompose(M3DPrimitive * m1, M3DPrimitive * m2, M3DPrimitive * mNew);

	// Average of an array of atoms.
	// If weights is NULL, 1/mNum will be the default weight.
	bool atomAverage(int mNum, M3DPrimitive ** mList, M3DPrimitive * mNew,
		double * weights = NULL);

	// Interpolation between atoms m1/m2 by t in [0, 1]
	bool atomInterp(double t, M3DPrimitive * m1, M3DPrimitive * m2,
		M3DPrimitive * mNew);

	// Geodesic distance between atoms m1/m2
	double atomDistance(M3DPrimitive * m1, M3DPrimitive * m2);
#ifdef OLD_BINARY_PABLO
	// This function was in Binary Pablo v.14.  It was deprecated on 28 June 2006
	// by Qiong Han and A. G. Gash.
	double atomSquareDistance(M3DPrimitive* m1, M3DPrimitive* m2, double *xD = NULL,
		double *rD = NULL, double *eD = NULL, double *n0D = NULL, double *n1D = NULL);
#else
	double atomSquareDistance(M3DPrimitive * m1, M3DPrimitive * m2,
		double * radius = NULL);
#endif
	double atomSetSquareDistance(int nPrims, M3DPrimitive ** mSet1,
		M3DPrimitive ** mSet2, double *radius = NULL);

	M3DSymPrimitive atomToSymAtom(M3DPrimitive * m);
	bool symAtomToAtom(M3DSymPrimitive * sM, M3DPrimitive * m);

	bool applyVector(VectorND & f, int nPrims, M3DPrimitive ** mSetBase,
		M3DPrimitive ** mSet);

private:

	bool ignorePrimType;

};


#endif

