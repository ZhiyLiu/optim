/*-----------------------------------------------------------
To calculate the predicted atom by the neighboring atoms
of an atom in a m-rep figure
by Qiong Han, 12/30/2005
-----------------------------------------------------------*/

#include <cmath>
#include <cstdio>
#include <ctime>
#include <sys/timeb.h>

#include "M3DQuadFigure.h"
#include "M3DAtomPredictorQuad.h"


#define CUBIC_HERMITE_INTERP
//#define BLERP

#define USE_RSRAD
#ifdef USE_RSRAD
	#define CALCULATE_RSRAD_DIRECTLY
	#define LEGALIZE_ILLEGAL_SRAD
#endif

#define USE_NEW_R_INTERPOLATION


//#define ATOM_PREDICTOR_TESTING
//#define VERBOSE


// debug: speed testing
#define DEBUGGING_TOOL
#define RENDER_ENTIRE_INTERPOLATED_FIGURE
#define INTERPOLATION_IMPROVEMENT


void M3DAtomPredictorQuad::setInitialFigures(M3DFigure * _frechet, M3DFigure * _fig)
{
	frechetMeanFig = _frechet;
	refFig = _fig;

	if(refPredFig!=NULL) {
		delete refPredFig;
	}
	if(refDiffFig!=NULL) {
		delete refDiffFig;
	}
	refPredFig = dynamic_cast<M3DQuadFigure*>(_fig)->clone();
	refDiffFig = dynamic_cast<M3DQuadFigure*>(_fig)->clone();
	tmpFig = dynamic_cast<M3DQuadFigure*>(_fig)->clone();

	if(radii != NULL)
		delete []radii;
	radii = new double[_fig->getPrimitiveCount()];
}

// Estimate spoke derivatives for internal atoms
#define EPSILON 1E-4
#define ZERO_THRESHOLD 1E-2
bool M3DAtomPredictorQuad::getDervs(M3DQuadFigure * _candFig, int row, int col,
	Vector3D U1[], Vector3D Xv[], double Rv[], Vector3D dSdv[2][2])
{

	// step 1
	//	|			1
	//	V		0	c	2	-> v1
	//	v2			3
	int				i,
					rowC = _candFig->getRowCount(),
					colC = _candFig->getColumnCount();
	double			v[2];
	M3DQuadEndPrimitive	mTmp[4];
	M3DPrimitive	*mCur,
					*mNgh[4];

	mCur = _candFig->getPrimitivePtr(row, col);
	U1[0] = mCur->getNormalizedY0();
	U1[1] = mCur->getNormalizedY1();
	U1[2] = mCur->getN();

	//double divider;

	if(row == 0) {
		if(col == 0){
			// corner atom [0, 0]
			v[0] = 1.0-EPSILON;
			v[1] = EPSILON;

			mNgh[0] = _candFig->getPrimitivePtr(1, 0);
			mNgh[2] = _candFig->getPrimitivePtr(0, 1);
			mNgh[3] = _candFig->getPrimitivePtr(1, 1);
			mTmp[0].atomInterp(    v[0], mNgh[0], mCur);
			mTmp[2].atomInterp(1.0-v[0], mCur, mNgh[2]);
			mTmp[1].atomInterp(v[1], mCur, mNgh[3]);
			
			for(i=0; i<2; i++){
				if(i == 0) {
					Xv[i] = (mTmp[2].getX()-mTmp[i].getX())/EPSILON/1.414;
					dSdv[0][i] = (mTmp[2].getY0()-mTmp[i].getY0())/EPSILON/1.414;
					dSdv[1][i] = (mTmp[2].getY1()-mTmp[i].getY1())/EPSILON/1.414;
					
					Rv[i] = (mTmp[2].getR()-mTmp[i].getR())/EPSILON/1.414;
				}
				else {
					Xv[i] = (mTmp[i].getX()-mCur->getX())/EPSILON/1.414;
					dSdv[0][i] = (mTmp[i].getY0()-mCur->getY0())/EPSILON/1.414;
					dSdv[1][i] = (mTmp[i].getY1()-mCur->getY1())/EPSILON/1.414; 

					Rv[i] = (mTmp[i].getR()-mCur->getR())/EPSILON/1.414;
				}
			}
		}
		else if(col == colC-1) {
			// corner atom [0, colC-1]
			v[0] = 1.0-EPSILON;
			v[1] = EPSILON;
			
			mNgh[0] = _candFig->getPrimitivePtr(0, colC-2);
			mNgh[2] = _candFig->getPrimitivePtr(1, colC-1);
			mNgh[3] = _candFig->getPrimitivePtr(1, colC-2);
			mTmp[0].atomInterp(    v[0], mNgh[0], mCur);
			mTmp[2].atomInterp(1.0-v[0], mCur, mNgh[2]);
			mTmp[1].atomInterp(v[1], mCur, mNgh[3]);
			
			for(i=0; i<2; i++){
				if(i == 0) {
					Xv[i] = (mTmp[2].getX()-mTmp[i].getX())/EPSILON/1.414;
					dSdv[0][i] = (mTmp[2].getY0()-mTmp[i].getY0())/EPSILON/1.414;
					dSdv[1][i] = (mTmp[2].getY1()-mTmp[i].getY1())/EPSILON/1.414;

					Rv[i] = (mTmp[2].getR()-mTmp[i].getR())/EPSILON/1.414;
				}
				else {
					Xv[i] = (mTmp[i].getX()-mCur->getX())/EPSILON/1.414;
					dSdv[0][i] = (mTmp[i].getY0()-mCur->getY0())/EPSILON/1.414;
					dSdv[1][i] = (mTmp[i].getY1()-mCur->getY1())/EPSILON/1.414;

					Rv[i] = (mTmp[i].getR()-mCur->getR())/EPSILON/1.414;
				}
			}
		}
		else {
			// end atoms at the top edge
			v[0] = 1.0-EPSILON;
			v[1] = EPSILON;
			
			mNgh[0] = _candFig->getPrimitivePtr(0, col-1);
			mNgh[2] = _candFig->getPrimitivePtr(0, col+1);
			mNgh[3] = _candFig->getPrimitivePtr(1, col);
			mTmp[0].atomInterp(v[0],     mNgh[0], mCur);
			mTmp[2].atomInterp(1.0-v[0], mCur, mNgh[2]);
			mTmp[3].atomInterp(v[1],     mCur, mNgh[3]);
			
			for(i=0; i<2; i++){
				if(i == 0) {
					Xv[i] = (mTmp[2].getX()-mTmp[i].getX())/EPSILON/2.0;
					dSdv[0][i] = (mTmp[2].getY0()-mTmp[i].getY0())/EPSILON/2.0;
					dSdv[1][i] = (mTmp[2].getY1()-mTmp[i].getY1())/EPSILON/2.0;

					Rv[i] = (mTmp[2].getR()-mTmp[i].getR())/EPSILON/2.0;
				}
				else {
					Xv[i] = (mTmp[3].getX()-mCur->getX())/EPSILON;
					dSdv[0][i] = (mTmp[3].getY0()-mCur->getY0())/EPSILON;
					dSdv[1][i] = (mTmp[3].getY1()-mCur->getY1())/EPSILON;

					Rv[i] = (mTmp[3].getR()-mCur->getR())/EPSILON;
				}
			}
		}
	}
	else if(row == rowC-1){
		if(col == 0){
			// corner atom [rowC-1, 0]
			v[0] = 1.0-EPSILON;
			v[1] = 1.0-EPSILON;
			
			mNgh[0] = _candFig->getPrimitivePtr(rowC-2, 0);
			mNgh[2] = _candFig->getPrimitivePtr(rowC-1, 1);
			mNgh[1] = _candFig->getPrimitivePtr(rowC-2, 1);
			mTmp[0].atomInterp(    v[0], mNgh[0], mCur);
			mTmp[2].atomInterp(1.0-v[0], mCur, mNgh[2]);
			mTmp[1].atomInterp(v[1], mNgh[1], mCur);
			
			for(i=0; i<2; i++){
				if(i == 0) {
					Xv[i] = (mTmp[2].getX()-mTmp[i].getX())/EPSILON/1.414;
					dSdv[0][i] = (mTmp[2].getY0()-mTmp[i].getY0())/EPSILON/1.414;
					dSdv[1][i] = (mTmp[2].getY1()-mTmp[i].getY1())/EPSILON/1.414;

					Rv[i] = (mTmp[2].getR()-mTmp[i].getR())/EPSILON/1.414;
				}
				else {
					Xv[i] = (mCur->getX()-mTmp[i].getX())/EPSILON/1.414;
					dSdv[0][i] = (mCur->getY0()-mTmp[i].getY0())/EPSILON/1.414;
					dSdv[1][i] = (mCur->getY1()-mTmp[i].getY1())/EPSILON/1.414;

					Rv[i] = (mCur->getR()-mTmp[i].getR())/EPSILON/1.414;
				}
			}
		}
		else if(col == colC-1) {
			// corner atom [rowC-1, colC-1]
			v[0] = 1.0-EPSILON;
			v[1] = 1.0-EPSILON;
			
			mNgh[0] = _candFig->getPrimitivePtr(rowC-1, colC-2);
			mNgh[2] = _candFig->getPrimitivePtr(rowC-2, colC-1);
			mNgh[1] = _candFig->getPrimitivePtr(rowC-2, colC-2);
			mTmp[0].atomInterp(    v[0], mNgh[0], mCur);
			mTmp[2].atomInterp(1.0-v[0], mCur, mNgh[2]);
			mTmp[1].atomInterp(v[1], mNgh[1], mCur);
			
			for(i=0; i<2; i++){
				if(i == 0) {
					Xv[i] = (mTmp[2].getX()-mTmp[i].getX())/EPSILON/1.414;
					dSdv[0][i] = (mTmp[2].getY0()-mTmp[i].getY0())/EPSILON/1.414;
					dSdv[1][i] = (mTmp[2].getY1()-mTmp[i].getY1())/EPSILON/1.414;

					Rv[i] = (mTmp[2].getR()-mTmp[i].getR())/EPSILON/1.414;
				}
				else {
					Xv[i] = (mCur->getX()-mTmp[i].getX())/EPSILON/1.414;
					dSdv[0][i] = (mCur->getY0()-mTmp[i].getY0())/EPSILON/1.414;
					dSdv[1][i] = (mCur->getY1()-mTmp[i].getY1())/EPSILON/1.414;

					Rv[i] = (mCur->getR()-mTmp[i].getR())/EPSILON/1.414;
				}
			}
		}
		else {
			// end atoms at the bottom edge
			v[0] = 1.0-EPSILON;
			v[1] = 1.0-EPSILON;
			
			mNgh[0] = _candFig->getPrimitivePtr(row, col-1);
			mNgh[2] = _candFig->getPrimitivePtr(row, col+1);
			mNgh[1] = _candFig->getPrimitivePtr(row-1, col);
			mTmp[0].atomInterp(v[0],     mNgh[0], mCur);
			mTmp[2].atomInterp(1.0-v[0], mCur, mNgh[2]);
			mTmp[1].atomInterp(v[1],     mNgh[1], mCur);
			
			for(i=0; i<2; i++){
				if(i == 0) {
					Xv[i] = (mTmp[2].getX()-mTmp[i].getX())/EPSILON/2.0;
					dSdv[0][i] = (mTmp[2].getY0()-mTmp[i].getY0())/EPSILON/2.0;
					dSdv[1][i] = (mTmp[2].getY1()-mTmp[i].getY1())/EPSILON/2.0;

					Rv[i] = (mTmp[2].getR()-mTmp[i].getR())/EPSILON/2.0;
				}
				else {
					Xv[i] = (mCur->getX()-mTmp[1].getX())/EPSILON;
					dSdv[0][i] = (mCur->getY0()-mTmp[1].getY0())/EPSILON;
					dSdv[1][i] = (mCur->getY1()-mTmp[1].getY1())/EPSILON;

					Rv[i] = (mCur->getR()-mTmp[1].getR())/EPSILON;
				}
			}
		}
	}
	else if(col == 0) {
		// end atoms at the left edge
		v[0] = EPSILON;
		v[1] = 1.0-EPSILON;
		
		mNgh[2] = _candFig->getPrimitivePtr(row, col+1);
		mNgh[1] = _candFig->getPrimitivePtr(row-1, col);
		mNgh[3] = _candFig->getPrimitivePtr(row+1, col);
		mTmp[2].atomInterp(v[0],     mCur, mNgh[2]);
		mTmp[1].atomInterp(v[1],     mNgh[1], mCur);
		mTmp[3].atomInterp(1.0-v[1], mCur, mNgh[3]);
		
		for(i=0; i<2; i++){
			if(i == 0) {
				Xv[i] = (mTmp[2].getX()-mCur->getX())/EPSILON;
				dSdv[0][i] = (mTmp[2].getY0()-mCur->getY0())/EPSILON;
				dSdv[1][i] = (mTmp[2].getY1()-mCur->getY1())/EPSILON;

				Rv[i] = (mTmp[2].getR()-mCur->getR())/EPSILON;
			}
			else {
				Xv[i] = (mTmp[3].getX()-mTmp[1].getX())/EPSILON/2.0;
				dSdv[0][i] = (mTmp[3].getY0()-mTmp[1].getY0())/EPSILON/2.0;
				dSdv[1][i] = (mTmp[3].getY1()-mTmp[1].getY1())/EPSILON/2.0;

				Rv[i] = (mTmp[3].getR()-mTmp[1].getR())/EPSILON/2.0;
			}
		}
	}
	else if(col == colC-1) {
		// end atoms at the right edge
		v[0] = 1.0-EPSILON;
		v[1] = 1.0-EPSILON;
		
		mNgh[0] = _candFig->getPrimitivePtr(row, col-1);
		mNgh[1] = _candFig->getPrimitivePtr(row-1, col);
		mNgh[3] = _candFig->getPrimitivePtr(row+1, col);
		mTmp[0].atomInterp(v[0],     mNgh[0], mCur);
		mTmp[1].atomInterp(v[1],     mNgh[1], mCur);
		mTmp[3].atomInterp(1.0-v[1], mCur, mNgh[3]);
		
		for(i=0; i<2; i++){
			if(i == 0) {
				Xv[i] = (mCur->getX()-mTmp[i].getX())/EPSILON;
				dSdv[0][i] = (mCur->getY0()-mTmp[i].getY0())/EPSILON;
				dSdv[1][i] = (mCur->getY1()-mTmp[i].getY1())/EPSILON;

				Rv[i] = (mCur->getR()-mTmp[i].getR())/EPSILON;
			}
			else {
				Xv[i] = (mTmp[3].getX()-mTmp[1].getX())/EPSILON/2.0;
				dSdv[0][i] = (mTmp[3].getY0()-mTmp[1].getY0())/EPSILON/2.0;
				dSdv[1][i] = (mTmp[3].getY1()-mTmp[1].getY1())/EPSILON/2.0;

				Rv[i] = (mTmp[3].getR()-mTmp[1].getR())/EPSILON/2.0;
			}
		}
	}
	else {
		// internal atoms
		v[0] = 1.0-EPSILON;
		v[1] = 1.0-EPSILON;
		
		mNgh[0] = _candFig->getPrimitivePtr(row, col-1);
		mNgh[2] = _candFig->getPrimitivePtr(row, col+1);
		mNgh[1] = _candFig->getPrimitivePtr(row-1, col);
		mNgh[3] = _candFig->getPrimitivePtr(row+1, col);
		mTmp[0].atomInterp(v[0],     mNgh[0], mCur);
		mTmp[2].atomInterp(1.0-v[0], mCur, mNgh[2]);
		mTmp[1].atomInterp(v[1],     mNgh[1], mCur);
		mTmp[3].atomInterp(1.0-v[1], mCur, mNgh[3]);
		
		for(i=0; i<2; i++){
			if(i == 0) {
				Xv[i] = (mTmp[2].getX()-mTmp[0].getX())/EPSILON/2.0;
				dSdv[0][i] = (mTmp[2].getY0()-mTmp[0].getY0())/EPSILON/2.0;
				dSdv[1][i] = (mTmp[2].getY1()-mTmp[0].getY1())/EPSILON/2.0;

				Rv[i] = (mTmp[2].getR()-mTmp[0].getR())/EPSILON/2.0;
			}
			else {
				Xv[i] = (mTmp[3].getX()-mTmp[1].getX())/EPSILON/2.0;
				dSdv[0][i] = (mTmp[3].getY0()-mTmp[1].getY0())/EPSILON/2.0;
				dSdv[1][i] = (mTmp[3].getY1()-mTmp[1].getY1())/EPSILON/2.0;

				Rv[i] = (mTmp[2].getR()-mTmp[0].getR())/EPSILON/2.0;
			}
		}
	}

	for(i=0; i<2; i++)
		Xv[i] -= mCur->getN() * (mCur->getN() * Xv[i]);

	return true;
}	// getDervs


// Calculate the "rKrel" from spoke derivatives of internal atoms
bool M3DAtomPredictorQuad::getRKrel(Vector3D * U1, Vector3D * Xv,
	Vector3D dSdv[2][2], double rKrel[], RSRAD_PENALTY_TYPE pelType)
{	
	// step 2
	
	int			i;
	Vector3D	Rt,
				U2;

	if(pelType == PEL_END) {
		for(i=0; i<4; i++) {
			rKrel[i] = -1.0;
		}
		return true;
	}

	for(i=0; i<2; i++) {
		Rt = Xv[0];
		U2 = dSdv[i][1].cross(U1[i]);	//Xv[0].cross(U1[i]);	//dSdv[i][1].cross(U1[i]);	//Xv[1].cross(U1[i]);
		//U2 = U2.cross(U1[i]);
		rKrel[2*i]   = -(dSdv[i][0] * U2) / (Rt * U2);

		Rt = Xv[1];
		U2 = dSdv[i][0].cross(U1[i]);	//Xv[1].cross(U1[i]);	//dSdv[i][0].cross(U1[i]);	//Xv[0].cross(U1[i]);
		//U2 = U2.cross(U1[i]);
		rKrel[2*i+1] = -(dSdv[i][1] * U2) / (Rt * U2);
	}

	return true;
}	// getRKrel


// Calculate the rSrad matrix from spoke derivatives of internal atoms
bool M3DAtomPredictorQuad::getRSrad(Vector3D * U1,
	Vector3D * Xv, double Rv[], Vector3D dSdv[2][2], Matrix rSrad[2])
{
	// step 2
	
	int			i;
//	double		Rv[2];
	Vector3D	tmp[2],
				dU1du[2], 
				dU1dv[2];

	Matrix		leftSide(2, 3), rightSide(3, 2), 
				Il(2, 2), Ir(2, 2), IrInv(2, 2),
				UTU(3, 3), I(3, 3),
				U(1, 3), UT(3, 1);

	dU1du[0] = dSdv[0][0];
	dU1dv[0] = dSdv[0][1];
	dU1du[1] = dSdv[1][0];
	dU1dv[1] = dSdv[1][1];

	for(i=0; i<2; i++) {

		I.setAll(0);
		I(0, 0) = 1;
		I(1, 1) = 1;
		I(2, 2) = 1;

		U(0, 0) = U1[i].getX();
		U(0, 1) = U1[i].getY();
		U(0, 2) = U1[i].getZ();
		UT = U.t();
		UTU = UT * U;
		UTU = UTU - I;

	//	tmp[0] = Xv[0] * (UTU);
		U(0, 0) = Xv[0].getX();
		U(0, 1) = Xv[0].getY();
		U(0, 2) = Xv[0].getZ();
		U = U * UTU;
		tmp[0].set(U(0, 0), U(0, 1), U(0, 2));

	//	tmp[1] = Xv[1] * (UTU);
		U(0, 0) = Xv[1].getX();
		U(0, 1) = Xv[1].getY();
		U(0, 2) = Xv[1].getZ();
		U = U * UTU;
		tmp[1].set(U(0, 0), U(0, 1), U(0, 2));

		dU1du[i] -= Rv[0] * U1[i];
		dU1dv[i] -= Rv[1] * U1[i];

		leftSide(0, 0)  = dU1du[i].getX();
		leftSide(0, 1)  = dU1du[i].getY();
		leftSide(0, 2)  = dU1du[i].getZ();
		leftSide(1, 0)  = dU1dv[i].getX();
		leftSide(1, 1)  = dU1dv[i].getY();
		leftSide(1, 2)  = dU1dv[i].getZ();
		
		// Q is [ tmp[0]
		//        tmp[1] ]
		// rightSide is Q^T
		rightSide(0, 0) = tmp[0].getX();
		rightSide(1, 0) = tmp[0].getY();
		rightSide(2, 0) = tmp[0].getZ();
		rightSide(0, 1) = tmp[1].getX();
		rightSide(1, 1) = tmp[1].getY();
		rightSide(2, 1) = tmp[1].getZ();

		Il = leftSide * rightSide;

		Ir(0, 0)        = tmp[0]*tmp[0];
		Ir(0, 1)        = tmp[0]*tmp[1];
		Ir(1, 0)        = tmp[1]*tmp[0];
		Ir(1, 1)        = tmp[1]*tmp[1];
		Ir.inverse(IrInv);
		
		Il              = Il*IrInv;
		//rSrad.setSize(2, 2);
		rSrad[i] = Il.t();	
	} // for(i=0, 1), both sides

	return true;
}	// getRSrad

// Calculate the rSrad matrix from spoke derivatives of internal atoms by 
// the compatibility condition
bool M3DAtomPredictorQuad::getRSradCompatibility(Vector3D * U1,
	Vector3D * Xv, Vector3D dSdv[2][2], Matrix rSrad[2])
{
	
	int			i;
	double		RvOld[2];
	Vector3D	tmp[2],
				dU1du[2], 
				dU1dv[2];
	Matrix		leftSide(2, 3), rightSide(3, 2), 
				IlOld(2, 2), Ir(2, 2), IrInv(2, 2);

	dU1du[0] = dSdv[0][0];
	dU1dv[0] = dSdv[0][1];
	dU1du[1] = dSdv[1][0];
	dU1dv[1] = dSdv[1][1];

	for(i=0; i<2; i++) {
		RvOld[0]    = -Xv[0] * U1[i];
		RvOld[1]    = -Xv[1] * U1[i];

		tmp[0] = -RvOld[0]*U1[i] - Xv[0];
		tmp[1] = -RvOld[1]*U1[i] - Xv[1];

		dU1du[i] -= RvOld[0] * U1[i];
		dU1dv[i] -= RvOld[1] * U1[i];

		leftSide(0, 0)  = dU1du[i].getX();
		leftSide(0, 1)  = dU1du[i].getY();
		leftSide(0, 2)  = dU1du[i].getZ();
		leftSide(1, 0)  = dU1dv[i].getX();
		leftSide(1, 1)  = dU1dv[i].getY();
		leftSide(1, 2)  = dU1dv[i].getZ();
		
		// Q is [ tmp[0]
		//        tmp[1] ]
		// rightSide is Q^T
		rightSide(0, 0) = tmp[0].getX();
		rightSide(1, 0) = tmp[0].getY();
		rightSide(2, 0) = tmp[0].getZ();
		rightSide(0, 1) = tmp[1].getX();
		rightSide(1, 1) = tmp[1].getY();
		rightSide(2, 1) = tmp[1].getZ();
		
		IlOld = leftSide * rightSide;

		Ir(0, 0)        = tmp[0]*tmp[0];
		Ir(0, 1)        = tmp[0]*tmp[1];
		Ir(1, 0)        = tmp[1]*tmp[0];
		Ir(1, 1)        = tmp[1]*tmp[1];
		Ir.inverse(IrInv);
		
		IlOld              = IlOld*IrInv;
		//rSrad.setSize(2, 2);
		rSrad[i] = IlOld.t();
	} // for(i=0, 1), both sides

	return true;
}	// getRSradCompatibility


// estimate spoke derivatives for internal atoms
bool M3DAtomPredictorQuad::getEndDervs(M3DQuadFigure * _candFig, int atomId,
	Vector3D U1[], Vector3D &Rt, Vector3D dSdt[])
{
	// [0]---[2]---[1]

	// step 1

	int	i, numNeighbors;
	M3DQuadEndPrimitive	endM[2];
	M3DQuadEndPrimitive	* mCur;
	M3DPrimitive * mNgh[3];
	Vector3D Xt[3], tmp[3];

	mCur = dynamic_cast<M3DQuadEndPrimitive *>(_candFig->getPrimitivePtr(atomId));
	if (mCur->type() != M3D_END_PRIMITIVE)
		return false;

	// a local frame {U1[0], U1[1], U1[2]} determined by each end atom
	// U1[0] is in the same direction as B (the bisector)
	U1[2] = mCur->getBPerp();	//U1[0].cross(U1[1]);

	_candFig->getAtomsNeighbors(atomId, numNeighbors, mNgh, EDGES_HAVE_2_NEIGHBORS);

	Xt[2] = (mNgh[1]->getX() - mNgh[0]->getX()) / 2.0;
	mNgh[2] = mCur;
	for(i=2; i<3; i++) {
		// project each Xt[i] onto the tagent plane 
		// at each atom's medial sheet
		tmp[0] = mNgh[i]->getN();
		Xt[i] -= tmp[0] * (tmp[0] * Xt[i]);
	}
	Rt = Xt[2];

	endM[0].atomInterp(1.0-EPSILON, mNgh[0], mNgh[2]);
	endM[1].atomInterp(EPSILON,     mNgh[2], mNgh[1]);

	dSdt[0] = (endM[1].getY0() - endM[0].getY0()) / 2.0 / EPSILON;
	dSdt[1] = (endM[1].getY1() - endM[0].getY1()) / 2.0 / EPSILON;
	dSdt[2] = (endM[1].getExtendedB() - endM[0].getExtendedB()) / 2.0 / EPSILON;
	
	return true;
}	// getEndDervs

// calculate the rKrel (relative radial curvature) from spoke derivatives of end atoms
bool M3DAtomPredictorQuad::getEndRKrel(Vector3D U1[], Vector3D Rt, Vector3D dSdt[], 
	double rKrel[], RSRAD_PENALTY_TYPE pelType)
{

	//	step 2

	int		i;

	if(pelType == PEL_INTERNAL) {
		for(i=0; i<3; i++)
			rKrel[i] = -1.0;
	}

	for(i=0; i<3; i++)
		rKrel[i] = - (dSdt[i] * U1[2]) / (Rt * U1[2]);

	return true;
}	// getEndRKrel


#define RKSRAD_4_DEGENERATE_RSRAD .55555
//#define MAX_ALLOWED_RKSRAD 2.00001
//#define MAX_ALLOWED_RKREL 2.00001
double M3DAtomPredictorQuad::getRSradPenalty(Matrix * rSrad, double * rKrel, 
	double r, double n, double a, double r3, double ne, double ae)
{
	// step 3: 
	//		calcualte the eigenvalues and eigenvectors of the rSrad matrix
	//		the penalty is calculated as a function of the larger real 
	//		eigenvalue of rSrad
	double		penalty = 1e10;

	int			i = 0, 
				k = 0;
	bool		bilerpSrad = false;
	Vector2D	eigenVectr[2];
	Vector		L(4);
	Matrix		LI(2, 2), V(2, 2), VInv(2, 2);

	if(a <= 0)
		a = 3.0;
	if(n <= 0)
		n = 3.0;

	if(rSrad == NULL)
		L(0) = *rKrel;	//penalize the rKrel for end atoms;
	else {
		if(4*(*rSrad)(1, 0)*(*rSrad)(0, 1)+
			((*rSrad)(0, 0)-(*rSrad)(1, 1))*((*rSrad)(0, 0)-(*rSrad)(1, 1)) <= 0) {
			// rSrad is degenerate: either with duplicate or complex
			// eigenvalues, which means a maximal penalty is needed
			//bilerpSrad = true;
#ifdef VERBOSE
			printf("DEGENERATE rSrad: duplicate or complex eigenvalues @ (%d, %d)!\n", k, i);
			rSrad->print();
#endif
			// penalize the case when the rSrad is not decomposable
			rSrad->factorEV(L, V, NON_SYM);
			//L(0) = RKSRAD_4_DEGENERATE_RSRAD;
			if(L(2)<-ZERO_THRESHOLD || L(2)>ZERO_THRESHOLD)
				L(0) = -sqrt(L(0)*L(0) + L(2)*L(2));
		/*	
			L(0) = sqrt(L(0)*L(0) + L(2)*L(2));
			if(L(1) > 0)
				L(0) *= 1.0;
			else
				L(0) *= -1.0;
		*/
		} 
		else {
			rSrad->factorEV(L, V, NON_SYM);
			if(L(0)<L(1))
				L(0) = L(1);
		}
	}

	double t;

	double	catsX2[10] = { -1.0,  -.5,  0.0,  0.5,  1.0,   1.5,   2.0,    2.5,     3.0};
	double	catsP2[10] = {  0.0,  0.2,  1.5,  8.0, 32.0, 128.0, 512.0, 2048.0,  8192.0};
	double	catsN2[10] = {  0.0,  0.5,  2.0,  8.0, 32.0, 128.0, 512.0, 2048.0,  8192.0};
	int		nCat2 = 8,
			cat2  = 0;

	double	catsX[10] = { -.5,  0.0,  0.5,  1.0,   1.5,   2.0,   2.5,    3.0};
	double	catsP[10] = { 0.0,  0.2,  1.5,  8.0,  32.0, 128.0, 512.0, 2048.0};
	double	catsN[10] = { 0.0,  0.5,  2.0,  8.0,  32.0, 128.0, 512.0, 2048.0};
	int		nCat = 7,
			cat  = 0;

	{
		// keep the same response function for crest reagion
		if(L(0) <= catsX[cat] || L(0) >= catsX[nCat]) {
			if (L(0) <= catsX[cat])
				penalty = 0;
			else {
				t = L(0) - catsX[nCat];
				penalty = catsP[nCat] + catsN[nCat]*t;
			}
		}
		else {
			while (L(0) > catsX[cat+1]) {
				cat++;
			}
			t = (L(0)-catsX[cat]) / (catsX[cat+1] - catsX[cat]);
			penalty =MyMath::HH0(t)* catsP[cat] +
				MyMath::HH1(t)* catsP[cat+1] +
				MyMath::HH2(t)* catsN[cat] +
				MyMath::HH3(t)* catsN[cat+1];
		}
	}

#ifdef RSRAD_GATHERING
		printf("%10.5f\n", L(0));
#endif

	return penalty;
}


//#define WORST_SPOKE_PENALTY_FOR_EACH_ATOM
double M3DAtomPredictorQuad::getAtomPenalty(M3DQuadFigure * _candFig, int atomId, 
											RSRAD_PENALTY_TYPE pelType,
											const double pnorm, const double threshold)
{
	double			penalty = 1e10,
					rKrel[6];
	M3DQuadEndPrimitive	mTmp[4];
	M3DPrimitive	*mNew = NULL;

	double			Rv[2];

	Vector3D		U1[3],
					Xv[2],
					Rt,
					dSdv[2][2],
					dSdt[6];
	Matrix			rSrad[2];

	int				i,
					//spokeNum,
					rowC = _candFig->getRowCount(),
					colC = _candFig->getColumnCount(),
					col = atomId % colC,
					row = atomId / colC;

	bool			realEndAtom = (_candFig->getPrimitivePtr(atomId)->type() == M3D_END_PRIMITIVE) && 
								  (pelType != PEL_INTERNAL);

	// calculate for one atom only
	
	getDervs(_candFig, row, col, U1, Xv, Rv, dSdv);

	for(i=0; i<2; i++)
		rSrad[i].setSize(2, 2);
	getRSrad(U1, Xv, Rv, dSdv, rSrad);

	penalty = 0;
	if( (pelType != PEL_END) 
	  ) {
		// exclude corner atoms from rSrad penalty
		// only rKrel penalty is used due to the stability

#ifdef WORST_SPOKE_PENALTY_FOR_EACH_ATOM
		penTmp  = getRSradPenalty(&rSrad[0]);
		if(penalty < penTmp)
			penalty = penTmp;
		penTmp  = getRSradPenalty(&rSrad[1]);
		if(penalty < penTmp)
			penalty = penTmp;
#else
		double temp = getRSradPenalty(&rSrad[0]);
		temp = (temp > threshold) ? temp - threshold : 0.0;
		penalty  = pow(temp, pnorm);
		temp = getRSradPenalty(&rSrad[1]);
		temp = (temp > threshold) ? temp - threshold : 0.0;
		penalty  += pow(temp, pnorm);
#endif
	}

	if(realEndAtom) {
		// for doubleSided corner atom penalties
		// U1[0], U1[1] hold the two Rt's
		// dSdt[0-2] and dSdt[3-5] hold the corresponding
		// spoke derivs
		bool doubleSidedEndAtom = true;
		getEndDervs(_candFig, atomId, U1, Rt, dSdt);
		getEndRKrel(U1, Rt, dSdt, rKrel);
		for(i=0; i<3; i++) {
#ifdef WORST_SPOKE_PENALTY_FOR_EACH_ATOM
			penTmp  = getRSradPenalty(NULL, &rKrel[i]);
			if(penalty < penTmp)
				penalty = penTmp;
#else
			double temp = getRSradPenalty(NULL, &rKrel[i]);
			temp = (temp > threshold) ? temp - threshold : 0.0;
			penalty  += pow(temp, pnorm);
#endif
		}

#ifndef WORST_SPOKE_PENALTY_FOR_EACH_ATOM
		if( (pelType != PEL_END) 
		  ) {
			penalty /= 5.0;
		}
		else
			penalty /= 3.0;
#endif
	}
	else
#ifndef WORST_SPOKE_PENALTY_FOR_EACH_ATOM
		penalty /= 2.0;
#else
		;
#endif	

	return penalty;
}


//#define WORST_SPOKE_PENALTY_FOR_FIGURE
double M3DAtomPredictorQuad::getFigureRSradPenalty(M3DFigure * _qcandFig, 
	int atomId, RSRAD_PENALTY_TYPE pelType, const double pnorm, const double threshold)
{
	int i;
	double penalty = 1e10;

	M3DQuadFigure* _candFig	= dynamic_cast<M3DQuadFigure*>(_qcandFig);
	if (atomId == -1) {
		// calculate for the entire figure
		// as the averaged rSrad penalty for all atoms

#ifdef WORST_SPOKE_PENALTY_FOR_FIGURE
		penalty = -1e5;
#else
		penalty = 0;
#endif

		for(i=0; i<_candFig->getPrimitiveCount(); i++) {
#ifdef WORST_SPOKE_PENALTY_FOR_FIGURE
			double penTmp   = getAtomPenalty(_candFig, i, pelType, pnorm, threshold);
			if(penalty < penTmp)
				penalty = penTmp;
#else
			penalty += getAtomPenalty(_candFig, i, pelType, pnorm, threshold);
#endif
		}

#ifndef WORST_SPOKE_PENALTY_FOR_FIGURE
		penalty /= _candFig->getPrimitiveCount();
#endif
	}
	else {
		// calculate for one atom only
		penalty = getAtomPenalty(_candFig, atomId, pelType, pnorm, threshold);
	}	// atomId \in [0, primCount-1]

#ifdef WORST_SPOKE_PENALTY_FOR_FIGURE
	return penalty;
#else
	return pow(penalty, 1.0/pnorm);
#endif
}

bool M3DAtomPredictorQuad::getFigureRSradValues(M3DQuadFigure * _candFig,
 int atomId, double rS[], RSRAD_PENALTY_TYPE pelType)
{
	double rKrel[3];
	M3DQuadEndPrimitive	mTmp[4];
	double		Rv[2];
	Vector3D	U1[3], Xv[2], Rt, dSdv[2][2], dSdt[3];
	Matrix		rSrad[2];

	M3DPrimitive * mNew = NULL;
	double penalty = 1e10;

	if (atomId < 0 || atomId >= _candFig->getPrimitiveCount())
		return false;

	int				i,
					//spokeNum,
					rowC = _candFig->getRowCount(),
					colC = _candFig->getColumnCount(),
					col = atomId % colC,
					row = atomId / colC;


	bool			realEndAtom = (_candFig->getPrimitivePtr(atomId)->type() == M3D_END_PRIMITIVE) && 
								  (pelType != PEL_INTERNAL);

	Vector		L(4);
	Matrix		LI(2, 2), V(2, 2), VInv(2, 2);

	if(_candFig==NULL || atomId<0 || atomId>=_candFig->getPrimitiveCount() || rS==NULL)
		return false;

	for(i=0; i<2; i++)
		rSrad[i].setSize(2, 2);

	getDervs(_candFig, row, col, U1, Xv, Rv, dSdv);
	getRSrad(U1, Xv, Rv, dSdv, rSrad);

	for(i=0; i<2; i++) {
		if(4*(rSrad[i])(1, 0)*(rSrad[i])(0, 1)+
		  ((rSrad[i])(0, 0)-(rSrad[i])(1, 1))*((rSrad[i])(0, 0)-(rSrad[i])(1, 1)) <= 0) {
			// rSrad is degenerate: either with duplicate or complex
			// eigenvalues
#ifdef VERBOSE
			printf("DEGENERATE rSrad: duplicate or complex eigenvalues @ (%d, %d)!\n", k, i);
			rSrad[i].print();
#endif
			// penalize the case when the rSrad is not decomposable
			rSrad[i].factorEV(L, V, NON_SYM);
			//L(0) = RKSRAD_4_DEGENERATE_RSRAD;
			if(L(2)<-ZERO_THRESHOLD || L(2)>ZERO_THRESHOLD)
				L(0) = -sqrt(L(0)*L(0) + L(2)*L(2));
		}
		else {
			rSrad[i].factorEV(L, V, NON_SYM);
			if (L(0)<L(1))
				L(0) = L(1);
		}

		if( (pelType != PEL_END) 
		  ) {
			rS[i] = L(0);
		}
		else
			rS[i] = -1.0;
	}

	if(realEndAtom) {
		getEndDervs(_candFig, atomId, U1, Rt, dSdt);
		getEndRKrel(U1, Rt, dSdt, rKrel);
		for(i=0; i<3; i++) {
			rS[i+2] = rKrel[i];
		}
	}

	return true;
}

