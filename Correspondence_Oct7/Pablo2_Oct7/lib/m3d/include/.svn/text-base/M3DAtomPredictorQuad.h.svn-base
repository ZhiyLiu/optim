#ifndef M3D_ATOM_PREDICTOR_QUAD_H
#define M3D_ATOM_PREDICTOR_QUAD_H

#include "M3DAtomPredictor.h"

class M3DAtomPredictorQuad : public M3DAtomPredictor {

protected:

	double			*radii;				// the atom radii in the inital Frechet mean figure
	M3DFigure		*frechetMeanFig,	// a pointer to the Frechet mean
					*refFig,			// a pointer to the initial reference figure

					*refPredFig,		// the predicted (by interpolation) figure by the reference
					*refDiffFig,		// the difference figure
					*tmpFig;			// temporary storage

public:

	M3DAtomPredictorQuad() {
		radii = NULL;

		frechetMeanFig = NULL;
		refFig = NULL;

		refPredFig = NULL;
		refDiffFig = NULL;
		tmpFig = NULL;
	}
    ~M3DAtomPredictorQuad() {
		if(radii != NULL)
			delete []radii;

		if(refPredFig != NULL)
			delete refPredFig;
		if(refDiffFig != NULL)
			delete refDiffFig;
		if(tmpFig != NULL)
			delete tmpFig;
	}
	
	// setup the reference figure and calculate the initial residue atoms
	void setInitialFigures(M3DFigure *_frechet, M3DFigure *_fig);

	// given the new figure, return the penalty for the atom with index atomId
	// the penalty is the geodesic distance between the two residue atoms
	double getAtomPenaltyWithoutOpt(M3DFigure *_candFig, int atomId) {
		return 0;
	}

	// calculate the predicted atom at (.5, .5), without any optimization
	void getPredictedAtomWithoutOpt(int numOfPredictors, 
									M3DPrimitive **predictors, 
									M3DPrimitive *predicted) {
	}

	// calculate the predicted atom from the neighboring atoms
	// either store the predicted result in figRes or back into fig if figRes is NULL
	void getPredictedAtomWithinFigure(M3DFigure *_fig, int atomId, 
		M3DFigure *_figRes = NULL, int type = 0) {
	}

	// for internal atoms
	bool getDervs(M3DQuadFigure * _candFig, int row, int col,
					Vector3D U1[], Vector3D Xv[], double Rv[], Vector3D dSdv[2][2]);
	// calculate the "rKrel" from spoke derivatives of internal atoms
	bool getRKrel(Vector3D * U1, Vector3D *Xv, Vector3D dSdv[2][2], double rKrel[],
					RSRAD_PENALTY_TYPE pelType = PEL_INTERNAL);
	bool getRSrad(Vector3D * U1, Vector3D *Xv, double Rv[], Vector3D dSdv[2][2], Matrix rSrad[]);
	bool getRSradCompatibility(Vector3D * U1, Vector3D * Xv, Vector3D dSdv[2][2], Matrix rSrad[]);

	// for end atoms
	bool getEndDervs(M3DQuadFigure* _candFig, int atomId,
						Vector3D U1[], Vector3D &Rt, Vector3D dSdt[]);
	bool getEndRKrel(Vector3D U1[], Vector3D Rt, Vector3D dSdt[], double rSrel[],
						RSRAD_PENALTY_TYPE pelType = PEL_INTERNAL);

	// the rSrad/rKe penalty response function
	double getRSradPenalty(Matrix* rSrad = NULL, double *rSrel = NULL, 
								 double r = 1.0, double n = 2.2, double a = 2.0,
								 double r3 = 1.0, double ne = 3.0, double ae = 3.0);

	// type - of penalty
	// 0 -	ignore rKrel pel
	// 1 -	ignore rSrad pel
	// 2 -	include both

	// calculate the rSrad penalty one atome at a time
	double getAtomPenalty(M3DQuadFigure *_candFig, int atomId, 
							RSRAD_PENALTY_TYPE pelType,
							const double pnorm, const double threshold );

	// return the rSrad penalty for any given atom
	// or teh average rSrad penalty for all atoms in the figure (when atom Id == -1)
	double getFigureRSradPenalty(M3DFigure *_candFig, int atomId,
									RSRAD_PENALTY_TYPE pelType,
									const double pnorm, const double threshold);

	bool getFigureRSradValues(M3DQuadFigure *_candFig, int atomId, double rS[], 
								RSRAD_PENALTY_TYPE pelType = PEL_INTERNAL);

//	void getInterpolatedSpokesNew(M3DQuadFigure *figure, int markedPrimId, 
//									int resR, int resC, int type = 0);
//	M3DQuadFigure* getInterpolatedFigure(M3DQuadFigure *figure, 
//											int lvl, int type = 0);

};

#endif
