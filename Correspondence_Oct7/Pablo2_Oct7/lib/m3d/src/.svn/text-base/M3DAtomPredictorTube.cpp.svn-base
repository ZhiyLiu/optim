#include <cmath>
#include <vector>
#include "Mathdefs.h"
#include <cstdio>
#include <ctime>
#include <sys/timeb.h>

#include "M3DAtomPredictorTube.h"

using std::vector;

const double dt		= 0.02;	// A step along the medial curve, should actually be called half_du
//const double step	= 0.1;	// Should always be >= 2*dt
const double step	= 0.1;	// Should always be >= 2*dt


double rSrad( const double phi, double dt, 
	const M3DTubePrimitive* prev, const M3DTubePrimitive* prim, const M3DTubePrimitive* next)
{
	assert( prev != NULL || next != NULL );
	if( prev == NULL ) {
		prev= prim;
		dt	/= 2;
	}
	else if( next == NULL ) {
		next= prim;
		dt	/= 2;
	}
	Vector3D U	= prim->getNormalizedYPhi(phi);
	Vector3D Sm, Sp;
	//
	// First find dS by dt.
	//
	Sm	= prev->getYPhi(phi);
	Sp	= next->getYPhi(phi);
	Vector3D dSbydu	= (Sp-Sm)/(2*dt);

	//
	// Now project along spoke S onto the tangent.
	// First remove any error components along the binormal direction.
	// The binormal here is T x U. T is dX/du
	//
	// FIXME: Note this way of getting differences for X is not right, we should
	// have used hermite in the first place to interpolate the atoms.
	// - update: this is an approximation as we are using Hermite now.
	const Vector3D T	= ((next->getX() - prev->getX()) / (2.0*dt));
	// bPerp will be zero iff the cone angle is zero or pi, but a 0/pi cone angle is not something valid.
	Vector3D bPerp	= T.cross(U);
	bPerp.normalize();
	const double errComponent	= dSbydu*bPerp;
	dSbydu	= dSbydu - errComponent * bPerp;

	//
	// Now dSbydu = r * dUbydu + drbydu * U
	// and dUbydu = a0 * U - a1 * dxbydu
	// a1 is Srad.
	// |U| = 1 => U*U = 1 and U * dUbydu = 0
	// So a0 = Srad * (dxbydu * U)
	// So dUbydu = Srad * ((dxbydu*U) * U - dxbydu )
	// and dSbydu = rSrad * ((dxbydu*U)*U - dxbydu ) + drbydu * U
	// rSrad is rKr here.
	//

	const double drbydu	= (next->getRPhi(phi) - prev->getRPhi(phi))/(2*dt);

	const Vector3D dr	= ((T*U)*U - T);
	const Vector3D nr	= (dSbydu - drbydu * U);
	//cout << nr.getX()/dr.getX() << " " << nr.getY()/dr.getY() << " " << nr.getZ()/dr.getZ() << " " << nr.norm()/dr.norm() << endl;
	const double rKr	= nr.norm() / dr.norm();

	return rKr;
}


double computeRSradPenalty( const M3DTubePrimitive* prev,
						   const M3DTubePrimitive* prim,
						   const M3DTubePrimitive* next,
						   const double pnorm,
						   const double threshold)
{
	//
	// First interpolate and find out the spokes in the neighboring region
	//
	double sigmarKrp	= 0.0;

	//
	// Find out the correct value for phi - one where the bending is the most.
	// This is given by finding out the plane in which the tangents lie
	// and then working over both the values of the phis.
	//
	double phi	= 0.0;
	vector<double> phis( prim->getNumberOfSpokes(), 0.0 );
	int iphi;

	int phiCount;
/*
	phiCount = 2;
	{
		//
		// Find the difference in the tangent vectors.
		//
		const Vector3D dT	= prim->getB() - prev->getB();
		//
		// Use this to compute a normal in the direction of maximum bending of the curve.
		// (basically second derivative of the curve)
		//
		Vector3D n	= dT - (dT*prim->getB()) * prim->getB();
		// Rotate into space of atom's frame.
		prim->getQ().conj().rotateVector(n);
		// Now find the angle made by n with the +ve y-axis in the yz plane.
		phis[0]	= atan2(n.getY(), n.getZ());
	}
	{
		//
		// Find the difference in the tangent vectors.
		//
		const Vector3D dT	= next->getB() - prim->getB();
		//
		// Use this to compute a normal in the direction of maximum bending of the curve.
		// (basically second derivative of the curve)
		//
		Vector3D n	= dT - (dT*prim->getB()) * prim->getB();
		// Rotate into space of atom's frame.
		prim->getQ().conj().rotateVector(n);
		// Now find the angle made by n with the +ve y-axis in the yz plane.
		phis[1]	= atan2(n.getY(), n.getZ());
	}
*/
	// The above heuristic only makes sense for tubes not quasi-tubes.
	phiCount	= prim->getNumberOfSpokes();
	for( iphi = 0; iphi != phiCount; iphi++ ) {
		phis[iphi]	= 2.0*R_PI*iphi/phiCount;
	}

	for( iphi = 0; iphi != phiCount; iphi++ ) {
		phi	= phis[iphi];
		const double rKr = rSrad( phi, dt, prev, prim, next );
		sigmarKrp	+= pow( (rKr < threshold) ? 0.0 : rKr - threshold, pnorm);
	}
	return sigmarKrp / phiCount;
}

double M3DAtomPredictorTube::getFigureRSradPenalty(
	M3DFigure *_fig, int atomId, RSRAD_PENALTY_TYPE pelType,
	const double pnorm, const double threshold)
{
	M3DTubeFigure* fig	= dynamic_cast<M3DTubeFigure*>(_fig);
	//
	// First make sure that the phi re-alignment is done.
	// DO NOT DO any realignment here, it should be done outside,
	// so that all penalties see a consistent effect.
	// (This is done in the deformation optimizer now)
	//
	fig->fixGlobalConsistency();

	double penalty	= 0.0;
	//
	// Step over all atoms finding srad penalty for them.
	//
	M3DTubePrimitive prim1, prim2, prim3;

	int totalPositions	= 0;
	double u;
	for( double i = 2.0*dt; i < (fig->getPrimitiveCount()-1)-2.0*dt; i += step, totalPositions++) {
		u	= i-dt;
		fig->atomAtCoordinates(&u, prim1);
		u	= i;
		fig->atomAtCoordinates(&u, prim2);
		u	= i+dt;
		fig->atomAtCoordinates(&u, prim3);
		penalty	+= computeRSradPenalty( &prim1, &prim2, &prim3, pnorm, threshold );
	}
	return pow(penalty / totalPositions, 1.0/pnorm);
}


/*
double computeRSradPenalty( const M3DTubePrimitive* prev,
						   const M3DTubePrimitive* prim,
						   const M3DTubePrimitive* next )
{
	//
	// First interpolate and find out the spokes in the neighboring region
	//
	const double dt		= 0.02;

	double sigmarKrp	= 0.0;

	//
	// Find out the correct value for phi - one where the bending is the most.
	// This is given by finding out the plane in which the tangents lie
	// and then working over both the values of the phis.
	//
	double phi	= 0.0, phis[] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
	int count	= 0;
	int iphi;
	if(prev != NULL) {
		//
		// Find the difference in the tangent vectors.
		//
		const Vector3D dT	= prim->getB() - prev->getB();
		//
		// Use this to compute a normal in the direction of maximum bending of the curve.
		// (basically second derivative of the curve)
		//
		Vector3D n	= dT - (dT*prim->getB()) * prim->getB();
		// Rotate into space of atom's frame.
		prim->getQ().conj().rotateVector(n);
		// Now find the angle made by n with the +ve y-axis in the yz plane.
		phis[count++]	= atan2(n.getY(), n.getZ());
	}
	if(next != NULL) {
		//
		// Find the difference in the tangent vectors.
		//
		const Vector3D dT	= next->getB() - prim->getB();
		//
		// Use this to compute a normal in the direction of maximum bending of the curve.
		// (basically second derivative of the curve)
		//
		Vector3D n	= dT - (dT*prim->getB()) * prim->getB();
		// Rotate into space of atom's frame.
		prim->getQ().conj().rotateVector(n);
		// Now find the angle made by n with the +ve y-axis in the yz plane.
		phis[count++]	= atan2(n.getY(), n.getZ());
	}
	int phiCount	= count;

	// Well, let's just discard what we did up here and use a sampling of phis along the
	// circumference.
	phiCount	= 8;
	for( iphi = 0; iphi != phiCount; iphi++ ) {
		phis[iphi]	= 2.0*R_PI*iphi/phiCount;
	}
	for( iphi = 0; iphi != phiCount; iphi++ ) {
		phi	= phis[iphi];
		Vector3D U	= prim->getNormalizedYPhi(phi);
		Vector3D spokes[2];
		spokes[0]	= U;
		const double weights[]	= { 1.0-dt, dt};
		Vector3D Sm, Sp;
		//
		// First find dS by dt.
		//
		if( prev == NULL ) {
			Sm	= prim->getR() * U;
		}
		else {
			spokes[1]	= prev->getNormalizedYPhi(phi);
			Sm	= exp( log(prim->getR()) * weights[0] + log(prev->getR()) * weights[1]) *
				ShapeSpace::S2::mean( spokes, 2, weights );
		}
		if( next == NULL ) {
			Sp	= prim->getR() * prim->getNormalizedYPhi(phi);
		}
		else {
			spokes[1]	= next->getNormalizedYPhi(phi);
			Sp	= exp( log(prim->getR()) * weights[0] + log(next->getR()) * weights[1]) *
				ShapeSpace::S2::mean( spokes, 2, weights );
		}
		Vector3D dSbydu	= (Sp-Sm)/(count*dt);

		//
		// Now project along spoke S onto the tangent.
		// First remove any error components along the binormal direction.
		// The binormal here is T x U.
		//
		const Vector3D T	= prim->getB() * (prim->getX() - ((next!=NULL) ? next->getX() : prev->getX()) ).norm();
		// bPerp will be zero iff the cone angle is zero or pi, but a 0/pi cone angle is not something valid.
		const Vector3D bPerp	= T.cross(U);
		dSbydu	= dSbydu - (dSbydu * bPerp) * bPerp;
		//
		// and then find a from dS/dt = a U + Kr T
		// T . N = 0 => a = dS/dt . N / U . N
		// Here, the N used is U - U . T
		//

		const Vector3D N	= U - (U * T) * T;
		const double a		= (dSbydu * N) / (U * N);

		dSbydu	= dSbydu - (a * U);
		const double rKr	= dSbydu.norm()/T.norm();
		sigmarKrp	+= pow(rKr, pnorm);
	}
	return sigmarKrp / count;
}

double M3DAtomPredictorTube::getFigureRSradPenalty(
	M3DFigure *_fig, int atomId, RSRAD_PENALTY_TYPE pelType)
{
	M3DTubeFigure* fig	= dynamic_cast<M3DTubeFigure*>(_fig);
	//
	// First do the phi re-alignment.
	// DO NOT DO any realignment here, it should be done outside,
	// so that all penalties see a consistent effect.
	// (This is done in the deformation optimizer now)
	//
	//fig->fixGlobalConsistency();

	double penalty	= 0.0;
	//
	// Step over all atoms finding srad penalty for them.
	//
	for( int i = 1; i < fig->getPrimitiveCount()-1; ++i) {
		penalty	+= computeRSradPenalty(
			dynamic_cast<M3DTubePrimitive*>(fig->getPrimitivePtr(i-1)),
			dynamic_cast<M3DTubePrimitive*>(fig->getPrimitivePtr(i)),
			dynamic_cast<M3DTubePrimitive*>(fig->getPrimitivePtr(i+1)) );
	}
	//
	// Do differently for end atoms.
	//
	penalty += computeRSradPenalty( NULL, 
		dynamic_cast<M3DTubePrimitive*>(fig->getPrimitivePtr(0)),
		dynamic_cast<M3DTubePrimitive*>(fig->getPrimitivePtr(1)) );
	penalty += computeRSradPenalty(
		dynamic_cast<M3DTubePrimitive*>(fig->getPrimitivePtr( fig->getPrimitiveCount()-2)),
		dynamic_cast<M3DTubePrimitive*>(fig->getPrimitivePtr( fig->getPrimitiveCount()-1)),
		NULL );
	return pow(penalty / fig->getPrimitiveCount(), 1.0/pnorm);
}
*/

