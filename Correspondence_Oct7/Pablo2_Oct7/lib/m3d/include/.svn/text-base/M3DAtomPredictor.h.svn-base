#ifndef M3D_ATOM_PREDICTOR_H
#define M3D_ATOM_PREDICTOR_H

#include "M3DFigure.h"
#include "matrix.h"
#include "Vector3D.h"
#include "GeodesicSym.h"

// handle weight tuple
class WeightTupleNew
{
private:
	int size;
	double *w;

	void allocateAndInit(int _size) {
		if(w!=NULL) {
			delete []w;
			w = NULL;
		}
		size = _size;
		w = new double[size];
		memset(w, 0, sizeof(double)*size);
	}

public:
	WeightTupleNew(int _size = 4) {
		w = NULL;
		if(_size<=0) {
			size = 0;
		}
		else {
			allocateAndInit(_size);
		}
	}

	~WeightTupleNew() {
		if(w!=NULL)
			delete []w;
	}


	int getSize() const {
		return size;
	}
	double *getWeights() {
		return w;
	}

	bool validDex(int dex) { 
		return (dex>=0 && dex<=size-1); 
	}

	void setSize(int _size) {
		if(_size<=0) {
			if(w!=NULL)
				delete []w;
			size = 0;
		}
		else {
			allocateAndInit(_size);
		}
	}
	void set(double *_w) {
		memcpy(w, _w, sizeof(double)*size);
	}
	void set(double _w, int dex) {
		if(validDex(dex))
			w[dex] = _w;
	}

	//copy constructor
	WeightTupleNew(const WeightTupleNew &_w) {
		w = NULL;
		allocateAndInit(_w.getSize());
		memcpy(w, _w.w, sizeof(double)*size);
	}

	// Arithmetic operations
	WeightTupleNew& operator = (const WeightTupleNew &_w) {
		memcpy(w, _w.w, sizeof(double)*size);
		return (*this);
	}

	// Arithmetic operations
	void operator += (const WeightTupleNew &_w) {
		for(int i=0; i<size; i++)
			w[i] += _w.w[i];
	}

	void operator *= (const double s) {
		for(int i=0; i<size; i++)
			w[i] *= s;
	}
	WeightTupleNew operator * (const double s) {
		WeightTupleNew wNew(size);
		for(int i=0; i<size; i++)
			wNew.w[i] = w[i]*s;
		return wNew;
	}
	void operator /= (const double s) {
		if(s!=0)
			for(int i=0; i<size; i++)
				w[i] /= s;
	}
	WeightTupleNew operator / (const double s) {
		WeightTupleNew wNew(size);
		if(s!=0)
			for(int i=0; i<size; i++)
				wNew.w[i] = w[i]/s;
		return wNew;
	}

	// This function is not defined.
	//friend WeightTupleNew operator + (WeightTupleNew &w1, WeightTupleNew &w2);
};



// generate some of the weights
class MyMath 
{
public:
static void blerpWeights(WeightTupleNew * w, WeightTupleNew * weights,
	int paramNum, double * params)
{
	WeightTupleNew w1(w->getSize()), w2(w->getSize());
	
	if(paramNum == 2) {
	/*
	w1 =  (weights[3]*params[0]);
	w1 += (weights[0]*(1.0-params[0]));
	w2 =  (weights[2]*params[0]);
	w2 += (weights[1]*(1.0-params[0]));
	
	  *w =  (w2*params[1]);
	  *w += (w1*(1.0-params[1]));
		*/
		w->set((1.0-params[0])*(1.0-params[1]), 0);
		w->set((1.0-params[0])*(    params[1]), 1);
		w->set((    params[0])*(    params[1]), 2);
		w->set((    params[0])*(1.0-params[1]), 3);
	}
	else {
		*w =  (weights[0]*(1.0-params[0]));
		*w += (weights[1]*params[0]);
	}
}

static double HH0(double t) {
	return (t*t * ( 2.0*t-3) + 1.0);	// p0
}

static double HH1(double t) {
	return (t*t * (-2.0*t+3)      );	// p1
}

static double HH2(double t) {
	return (t * (t * (t-2.0) + 1.0));	// deriv 0
}

static double HH3(double t) {
	return (t * (t * (t - 1.0) )  );	// deriv 1
}

static void hermiteWeights(int type, WeightTupleNew * w, double t,
	WeightTupleNew * weights, double * params = NULL) 
{
//	WeightTupleNew w1(w->getSize()), w2(w->getSize());
	
	if(type == 1) {
		// for the interpolation of hub position
		w->set(HH0(t), 0);
		w->set(HH1(t), 1);
		w->set(HH2(t), 2);
		w->set(HH3(t), 3);
	}
	else if(type==2 && params!=NULL) {
		// for the interpolation of spoke length and orientation
		*w =  weights[0] * (HH0(params[0]) * HH0(params[1]));
		*w += weights[1] * (HH0(params[0]) * HH1(params[1]));
		*w += weights[2] * (HH1(params[0]) * HH1(params[1]));
		*w += weights[3] * (HH1(params[0]) * HH0(params[1]));

		blerpWeights(w, weights, 2, params);
	}
	else if(type==3) {
		w->set(HH0(t), 0);
		w->set(HH1(t), 1);
	}
	else {
		w->set(0, 0);
		w->set(0, 1);
		w->set(0, 2);
		w->set(0, 3);
	}
}

static double HHDeriv0(double t) {
	return 6*t*( t - 1);
}

static double HHDeriv1(double t) {
	return 6*t*(-t + 1);
}

static double HHDeriv2(double t) {
	return t*(3*t - 4) +1;
}

static double HHDeriv3(double t) {
	return t*(3*t - 2);
}

static void hermiteWeightsDeriv(int type, WeightTupleNew * w, double t,
	WeightTupleNew * weights, double * params = NULL)
{
//	WeightTupleNew w1(w->getSize()), w2(w->getSize());
	
	if(type == 1) {
		// for the interpolation of hub position
		w->set(HHDeriv0(t), 0);
		w->set(HHDeriv1(t), 1);
		w->set(HHDeriv2(t), 2);
		w->set(HHDeriv3(t), 3);
	}
	else if(type==2 && params!=NULL) {
		// for the interpolation of spoke length and orientation
		*w =  weights[0] * (HHDeriv0(params[0]) * HHDeriv0(params[1]));
		*w += weights[1] * (HHDeriv0(params[0]) * HHDeriv1(params[1]));
		*w += weights[2] * (HHDeriv1(params[0]) * HHDeriv1(params[1]));
		*w += weights[3] * (HHDeriv1(params[0]) * HHDeriv0(params[1]));

//		blerpWeights2(w, weights, 2, params);
	}
	else if(type==3) {
		w->set(HHDeriv0(t), 0);
		w->set(HHDeriv1(t), 1);
	}
	else {
		w->set(0, 0);
		w->set(0, 1);
		w->set(0, 2);
		w->set(0, 3);
	}
}

static void evaluateHub(int type, WeightTupleNew * w, Vector3D * P,
	double * params, Vector3D * hub)
{
	int i, j;
	WeightTupleNew weights;	
	Vector3D Q[4];

	// type = 0/1/2 to evaluate cubic Hermite patch of X/Xu/Xv
	if(type == 0 || type == 1 || type == 2)
	{
		// type == 0 to evaluate X
		if(type==1)
			// evaluate Xu
			hermiteWeightsDeriv(1, &weights, params[0], w);
		else
			hermiteWeights(1, &weights, params[0], w);
		for(i=0; i<4; i++) {
			Q[i].set(0, 0, 0);
			for(j=0; j<4; j++)
				Q[i] += (weights.getWeights())[j] * P[i*4+j];
		}
		if(type==2)
			// evaluate Xv
			hermiteWeightsDeriv(1, &weights, params[1], w);
		else
			hermiteWeights(1, &weights, params[1], w);
		hub->set(0, 0, 0);
		for(i=0; i<4; i++)
			*hub += Q[i] * (weights.getWeights())[i];
		//hermiteWeights2(2, &w, 0.0, weights, params);
		//			geo.atomAverageNew(numOfPredictors, predictors, &mNew, w.getWeights(), 1);
		//			mNew.setX(hubPos.getX(), hubPos.getY(), hubPos.getZ());
	}
	else if(type == 3) {
		hermiteWeights(1, &weights, params[0], w);
		hub->set(0, 0, 0);
		for(i=0; i<4; i++)
			*hub += P[i] * (weights.getWeights())[i];
	}
}

static void evaluateRad(int type, WeightTupleNew * w, double * rP,
	double * params, double * rad)
{
	int i, j;
	WeightTupleNew weights;	
	double rQ[4];

	{
		// type == 0 to evaluate R
		if(type==1)
			// evaluate Ru
			hermiteWeightsDeriv(1, &weights, params[0], w);
		else
			hermiteWeights(1, &weights, params[0], w);
		for(i=0; i<4; i++) {
			rQ[i]=0;
			for(j=0; j<4; j++)
				rQ[i] += (weights.getWeights())[j] * rP[i*4+j];
		}
		
		if(type==2)
			// evaluate Rv
			hermiteWeightsDeriv(1, &weights, params[1], w);
		else
			hermiteWeights(1, &weights, params[1], w);
		*rad=0;
		for(i=0; i<4; i++)
			*rad += rQ[i] * (weights.getWeights())[i];
		//hermiteWeights2(2, &w, 0.0, weights, params);
		//			geo.atomAverageNew(numOfPredictors, predictors, &mNew, w.getWeights(), 1);
		//			mNew.setX(hubPos.getX(), hubPos.getY(), hubPos.getZ());
	}	
}

static void evaluateSpoke(int type, WeightTupleNew * w, Vector3D * P,
	double * rP, double * params, M3DSymPrimitive * atomSym)
{
	int i;			//, j;

//	Vector3D Q[4];
//	double	rQ[4];

	Vector3D X, Xu, Xv, gR, N;
	double   R, Ru, Rv, cosA;

	WeightTupleNew WEIGHTS[4];
	for(i=0; i<4; i++) {
		WEIGHTS[i].setSize(4);
		WEIGHTS[i].set(1, i);
	}

	evaluateHub(0, WEIGHTS, P, params, &X);
	evaluateHub(1, WEIGHTS, P, params, &Xu);
	evaluateHub(2, WEIGHTS, P, params, &Xv);
	evaluateRad(0, WEIGHTS, rP, params, &R);
	evaluateRad(1, WEIGHTS, rP, params, &Ru);
	evaluateRad(2, WEIGHTS, rP, params, &Rv);

	Matrix Xuv(3, 2), Ix(2, 2), IxInv(2, 2);
	Vector gradR(3), Ruv(2);
	Xuv(0, 0) = Xu.getX();
	Xuv(0, 1) = Xv.getX();
	Xuv(1, 0) = Xu.getY();
	Xuv(1, 1) = Xv.getY();
	Xuv(2, 0) = Xu.getZ();
	Xuv(2, 1) = Xv.getZ();
	Ix(0, 0)  = Xu*Xu;
	Ix(0, 1)  = Xu*Xv;
	Ix(1, 0)  = Xv*Xu;
	Ix(1, 1)  = Xv*Xv;
	Ix.inverse(IxInv);
	Ruv(0)			 = Ru;
	Ruv(1)			 = Rv;
	gradR = Xuv * IxInv * Ruv;

	gR.set(gradR(0), gradR(1), gradR(2));
	cosA = gR.norm();
	N    = Xu.cross(Xv);
	N.normalize();
	atomSym->x  = X;
	atomSym->r  = R;
	if(cosA<=1.0){
	}
	else {
		cosA = 1;
	}

	atomSym->n0 = -gR + sqrt(1-cosA*cosA)*N;
	atomSym->n1 = -gR - sqrt(1-cosA*cosA)*N;
}

static void evaluateQuadInterpBSplineCurve(double t, Vector3D * p,
	Vector3D * hub)
{
	int i;
	double q[3];

	hub->set(0, 0, 0);
	q[0] = 1 + (-4)*t + ( 4)*t*t;
	q[1] =     ( 4)*t + (-4)*t*t;
	q[2] =     (-1)*t + ( 2)*t*t;
	for(i=0; i<3; i++)
		*hub += q[i] * p[i];
}
static void evaluateLagrangeCurve(int orderN, double t, Vector3D * p,
	Vector3D * hub)
{
	int	    i, j;
	double  grid[100],
			q[100],
			dex;

	hub->set(0, 0, 0);

	dex = t/(1.0/orderN);
	if(dex-(int)(dex) == 0) {
		// at a grid point
		*hub = p[(int)(dex)];
	}
	else {
		double  divider  = 1.0,
				dividentCommon = 1.0;

		for(i=0; i<orderN+1; i++)
			grid[i] = i*1.0/orderN;

		for(j=0; j<orderN+1; j++)
			dividentCommon *= t-grid[j];

		for(i=0; i<orderN+1; i++) {
			divider = 1;
			for(j=0; j<orderN+1; j++) {
				if(j==i)
					continue;
				divider *= (grid[i]-grid[j]);
			}

			q[i] = dividentCommon/divider/(t-grid[i]);
		}

		for(i=0; i<orderN+1; i++)
			*hub += q[i] * p[i];
	}
}

};


enum RSRAD_PENALTY_TYPE {
	PEL_INTERNAL,
	PEL_END,
	PEL_BOTH
};


class M3DAtomPredictor {
  public:
	virtual double getFigureRSradPenalty(M3DFigure *_candFig, int atomId,
								RSRAD_PENALTY_TYPE pelType,
								const double pnorm, const double threshold) = 0;
};


#include "M3DAtomPredictorQuad.h"
#include "M3DAtomPredictorTube.h"
#endif
