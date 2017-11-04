/*
    Member functions for class GeodesicSym

    Author:         Qiong Han
    Last Modified:  06/04
*/

#include <math.h>

#include "M3DPrimitive.h"
#include "M3DEndPrimitive.h"
#include "M3DFigure.h"
#include "GeodesicSym.h"
#include "utility.h"

using namespace std;

const double epsilon = 1.0e-03;
#define MEAN_ITERATION_THRESHOLD			1.0e-11
#define STANDARD_ATOM_LOG_MAP_VECTOR_SIZE	8
#define END_ATOM_LOG_MAP_VECTOR_SIZE		9
#define NUMERICAL_ZERO 1.0e-10
#define MAX_ITERATION_NUMBER 200

//#define debug

#define BUG_FIX

Vector2D GeodesicSym::sphereLogMap(Vector3D n)
{
	double d;

// This is only a way to avoid getting math errors.  There is an underlying
// problem that has not been resolved yet.  AGG.
#ifdef BUG_FIX
	d = n.getX();
	if (d > 1.0)
		d = 1.0;
	d = acos(d);
#else
	d = acos(n.getX());
#endif
//	if(d == 1.2096279737202726) {
//		cout << "Got it\n";
//	}
	if(d == 0)
		return Vector2D(0.0, 0.0);

	return (d/sin(d))*Vector2D(n.getY(), n.getZ());
}

Vector3D GeodesicSym::sphereExpMap(Vector2D v)
{
	double d, alpha;

	d = v.norm();
	if(d == 0)
		return Vector3D(1.0, 0.0, 0.0);

	alpha = sin(d)/d;
	return Vector3D(cos(d), alpha * v.getX(), alpha * v.getY());
}

bool GeodesicSym::atomCrossingSpaceCompose(M3DPrimitive * prim,  M3DPrimitive* deltaPrim)
{
	// Convert the  deltaprim in to symatom.. from lie group in to sysmetric space 

	M3DSymPrimitive symPrim ,symDeltaPrim;

    Quat q0,q1,q0_delta,q1_delta,Q0,Q1;

	symPrim= atomToSymAtom(prim);  

	symDeltaPrim=atomToSymAtom(deltaPrim);
	q0_delta=rotationFromOrigin(symDeltaPrim.n0);
	q1_delta=rotationFromOrigin(symDeltaPrim.n1);  

	Vector3D N0,N1;

	//rotate the orignial normal with the delta angle ,get the new symPrim
	q0_delta.rotateVector(symPrim.n0);
	q1_delta.rotateVector(symPrim.n1);

	symAtomToAtom(&symPrim,prim);

	return true;
}

bool GeodesicSym::applyVector(VectorND & f, int nPrims, M3DPrimitive ** mSetBase,
	M3DPrimitive ** mSet) 
{
	int			i;

	Vector3D	tran;
	Quat		q;
	double		s;

	tran.set(0, 0, 0);
	q.set(1, 0, 0, 0);
	s = 1.0;

	if(f.getSize() <= 0)
		return false;
	else if(f.getSize() == 1) {
		s = exp(f.get(0));
	}
	else if(f.getSize() == 3) {
		q.expMap(Vector3D(f.get(0), f.get(1), f.get(2)));
	}
	else if(f.getSize() == 4) {
		q = q.expMap(Vector3D(f.get(0), f.get(1), f.get(2)));
		s = exp(f.get(3));
	}
	else if(f.getSize() == 6) {
		tran.set(f.get(0), f.get(1), f.get(2));
		q = q.expMap(Vector3D(f.get(3), f.get(4), f.get(5)));
	}
	else if(f.getSize() == 7) {
		tran.set(f.get(0), f.get(1), f.get(2));
		q = q.expMap(Vector3D(f.get(3), f.get(4), f.get(5)));
		s = exp(f.get(6));
	}
	else 
		return false;

	Vector3D	cog,
				deltaX;

	cog.set(0, 0, 0);
	for(i=0; i<nPrims; i++) {
		cog += mSetBase[i]->getX();
	}
	cog /= nPrims;

	for(i=0; i<nPrims; i++) {
		*mSet[i] = *mSetBase[i];
		deltaX = mSet[i]->getX() - cog;

		if(f.getSize() > 1) {
			q.rotateVector(deltaX);
			mSet[i]->rotateBy(q);
		}
		if(f.getSize() == 1 ||
		   f.getSize() == 4 ||
		   f.getSize() == 7   ) {
			deltaX *= s;
			mSet[i]->scaleBy(s);
		}

		mSet[i]->setX( cog + deltaX );

		if(f.getSize() == 6 || f.getSize() == 7)
			mSet[i]->translateBy(tran);
	}

	return true;
}

double GeodesicSym::atomSetSquareDistance(int nPrims, M3DPrimitive **mSet1, 
										M3DPrimitive **mSet2, double *radius) 
{
	int i;
	double dis = 0;

	for (i = 0; i < nPrims; i++) {
		if (radius != NULL)
			dis += atomSquareDistance(mSet1[i], mSet2[i], &radius[i]);
		else
			dis += atomSquareDistance(mSet1[i], mSet2[i], NULL);
	}
	return dis;
}

bool GeodesicSym::atomsAlignment(int nPrims, M3DPrimitive** mSet1, M3DPrimitive** mSet2,
	Vector3D &tran, Quat &q, double &s, bool alignAll)
{
	int				i;
	M3DPrimitive	**mSetTmp = new M3DPrimitive* [nPrims];

	Vector3D	cogs[2];

	for(i=0; i<2; i++)
		cogs[i].set(0, 0, 0);


	for(i=0; i<nPrims; i++) {
		if(mSet1[i]->type() == M3D_END_PRIMITIVE)
			mSetTmp[i] = new M3DQuadEndPrimitive();
		else
			mSetTmp[i] = new M3DQuadPrimitive();

		*mSetTmp[i] = *mSet1[i];

		cogs[0] += mSet1[i]->getX();
		cogs[1] += mSet2[i]->getX();
	}

	for(i=0; i<2; i++)
		cogs[i] /= nPrims;

	tran = cogs[1] - cogs[0];

	int				iterNum = 0;
	double			prevDis,
					dis,
					disPlus, disMinus;

	// no need to optimize translation
	// opt on rotation (+scale) only
	VectorND		f,
					deriv;

	if(alignAll) {
		f.setSize(7);
		deriv.setSize(7);
	}
	else {
		f.setSize(6);
		deriv.setSize(6);
	}

	f.set(0.0);
	deriv.set(0.0);

	// translation is already determined by the difference 
	// between the two COGs

	// the rest (rotation/scaling) is calculated by minimizing
	// the sumed square geodesic distance between 2 sets of 
	// atoms, using gradient descent
	f.set(0, tran);

	dis = atomSetSquareDistance(nPrims, mSet1, mSet2);

	do {
		prevDis = dis;
		for(i=3; i<f.getSize(); i++) {
			f.set(i, f.get(i) + epsilon);
			applyVector(f, nPrims, mSet1, mSetTmp);
			disPlus = atomSetSquareDistance(nPrims, mSetTmp, mSet2);

			f.set(i, f.get(i) - epsilon*2.0);
			applyVector(f, nPrims, mSet1, mSetTmp);
			disMinus = atomSetSquareDistance(nPrims, mSetTmp, mSet2);

			f.set(i, f.get(i) + epsilon);

			deriv.set(i, (disPlus-disMinus) / 2.0 / epsilon );
		}

		double	stepSize = 20.0;
		int		trialNum = 0;
		f -= deriv * stepSize;
		applyVector(f, nPrims, mSet1, mSetTmp);
		dis = atomSetSquareDistance(nPrims, mSetTmp, mSet2);

		while(dis >= prevDis && trialNum < 40) {
			// recover the old "f" first
			f += deriv * stepSize;

			// cut the step size to half
			stepSize /= 2.0;
			f -= deriv * stepSize;
			applyVector(f, nPrims, mSet1, mSetTmp);
			dis = atomSetSquareDistance(nPrims, mSetTmp, mSet2);
			trialNum ++;
		}

		if(dis >= prevDis) {
			dis = prevDis;
			f  += deriv * stepSize;
		}

		iterNum++;
	}
	while(prevDis-dis > 1.0e-18 && 
		  iterNum < MAX_ITERATION_NUMBER           );

	if(f.getSize() <= 0)
		return false;
	else if(f.getSize() == 1) {
		s = exp(f.get(0));
	}
	else if(f.getSize() == 3) {
		q = q.expMap( Vector3D(f.get(0), f.get(1), f.get(2)) );
	}
	else if(f.getSize() == 4) {
		q = q.expMap( Vector3D(f.get(0), f.get(1), f.get(2)) );
		s = exp(f.get(3));
	}
	else if(f.getSize() == 6) {
		tran.set(f.get(0), f.get(1), f.get(2));
		q = q.expMap( Vector3D(f.get(3), f.get(4), f.get(5)) );
		s = 1.0;
	}
	else if(f.getSize() == 7) {
		tran.set(f.get(0), f.get(1), f.get(2));
		q = q.expMap( Vector3D(f.get(3), f.get(4), f.get(5)) );
		s = exp(f.get(6));
	}
	else 
		return false;

	return true;
}

// another possible bug in M3DPGA.cpp?
Quat GeodesicSym::rotationFromOrigin(Vector3D n)
{
	Vector2D logN;
	Vector3D axis;
	double angle;
	Quat q;

	logN = sphereLogMap(n);
	axis.set(0.0, n.getZ(), -n.getY());
	axis.normalize();
	angle = logN.norm();

	q.setAxisAngle(axis, -angle);
	return q;
}

Quat GeodesicSym::rotationToOrigin(Vector3D n)
{
	Vector2D logN;
	Vector3D axis;
	double angle;
	Quat q;

	logN = sphereLogMap(n);
	axis.set(0.0, n.getZ(), -n.getY());
	axis.normalize();
	angle = logN.norm();

	q.setAxisAngle(axis, angle);
	return q;
}

/* Exponential and Log map for atoms
   1.	Only the atom (X, R, Q, Theta) are updated here.
	Other fields including selected/type/hinge should be updated somewhere
	else.
   2.	All atoms are considered as standard atoms by default.
	Special treatment is needed for end atoms by instantialize GeodesicSym
	with 'FALSE'.
*/
bool GeodesicSym::atomLogMap(M3DPrimitive *m, VectorND *v)
{
	// v has different sizes for different atom type
	if(!ignorePrimType && m->type()==M3D_END_PRIMITIVE)
	{
		if(v==NULL)
			v=new VectorND(END_ATOM_LOG_MAP_VECTOR_SIZE);
		v->set(8, log((dynamic_cast<M3DQuadEndPrimitive *>(m))->getElongation()));
	}
	else
	{
		if(v==NULL)
			v=new VectorND(STANDARD_ATOM_LOG_MAP_VECTOR_SIZE);
	}

	v->set(0, m->getX()); 
	v->set(3, log(m->getR()));
	v->set(4, m->getQ().logMap());
	v->set(7, m->getTheta());
	return true;
}/*
bool GeodesicSym::symAtomLogMap(M3DSymPrimitive *sM, VectorND *v)
{
	// v has different sizes for different atom type
	if(!ignorePrimType)
	{
		if(v==NULL)
			v=new VectorND(END_ATOM_LOG_MAP_VECTOR_SIZE);
		v->setX(8, log(sM->elongation));
	}
	else
	{
		if(v==NULL)
			v=new VectorND(STANDARD_ATOM_LOG_MAP_VECTOR_SIZE);
	}

	v->setX(0, sM->x); 
	v->setX(3, log(sM->r));
	v->setX(4, sphereLogMap(sM->n0));
	v->setX(6, sphereLogMap(sM->n1));
	return true;
}
*/
bool GeodesicSym::atomExpMap(VectorND &v, M3DPrimitive *mNew)
{
	Vector3D logQ;
	Quat q;

	// to treat standard/end atoms differently
	if(!ignorePrimType && v.getSize()==END_ATOM_LOG_MAP_VECTOR_SIZE)
	{
		if(mNew==NULL)
			mNew=new M3DQuadEndPrimitive();
		(dynamic_cast<M3DQuadEndPrimitive *>(mNew))->setElongation(exp(v.get(8)));
	}
	else
	{
		if(mNew==NULL)
			mNew=new M3DQuadPrimitive();
	}

	logQ.set(v.get(4), v.get(5), v.get(6));
	mNew->setX(v.get(0), v.get(1), v.get(2));
	mNew->setR(exp(v.get(3)));
	mNew->setQ(q.expMap(logQ));
	mNew->setTheta(v.get(7));

	return true;
}

bool GeodesicSym::symAtomExpMap(VectorND &v, M3DSymPrimitive *mNew)
{
	//Vector2D logN;
	//Quat q;

	// to treat standard/end atoms differently
	if(!ignorePrimType && v.getSize()==END_ATOM_LOG_MAP_VECTOR_SIZE)
	{
		if(mNew==NULL)
			mNew=new M3DSymPrimitive();
		mNew->elongation=exp(v.get(8));
	}
	else
	{
		if(mNew==NULL)
			mNew=new M3DSymPrimitive();
	}

	//logQ.set(v.get(4), v.get(5), v.get(6));
	mNew->x=Vector3D(v.get(0), v.get(1), v.get(2));
	mNew->r=exp(v.get(3));
	//logN.set(v.get(4), v.get(5);
	mNew->n0=sphereExpMap(Vector2D(v.get(4), v.get(5)));
	mNew->n1=sphereExpMap(Vector2D(v.get(6), v.get(7)));

	return true;
}

bool GeodesicSym::atomInverse(M3DPrimitive *m, M3DPrimitive *mNew)
{
	if(!ignorePrimType && m->type()==M3D_END_PRIMITIVE)
	{
		if(mNew==NULL)
			mNew=new M3DQuadEndPrimitive();
		(dynamic_cast<M3DQuadEndPrimitive *>(mNew))->setElongation(1.0/
			 (dynamic_cast<M3DQuadEndPrimitive *>(m))->getElongation());
	}
	else
	{
		if(mNew==NULL)
			mNew=new M3DQuadPrimitive();
	}
	mNew->setX(-m->getX());
	mNew->setR(1.0/m->getR());
	mNew->setQ(m->getQ().conj());
	mNew->setTheta(-m->getTheta());

	return true;
}

/*
void regularizeQuaternion(M3DPrimitive *mNew)
{
	Quat qTmp;
	qTmp=mNew->getQ();
	if(qTmp.getW()<0 && fabs(qTmp.getW())<NUMERICAL_ZERO)
		qTmp.setW(0);
	mNew->setQ(qTmp);
}
*/
bool GeodesicSym::atomCompose(M3DPrimitive *m1, M3DPrimitive *m2, M3DPrimitive *mNew)
{
	double thetaTmp;

	if (m2==NULL)
		cout << "Error: m2 is NULL in GeodesicSym::atomCompose()" << endl;

	if(!ignorePrimType && m1->type()==m2->type() && m1->type()==M3D_END_PRIMITIVE)
	{
		if(mNew==NULL)
			mNew=new M3DQuadEndPrimitive();
		(dynamic_cast<M3DQuadEndPrimitive *>(mNew))->setElongation((dynamic_cast<M3DQuadEndPrimitive *>(m1))->getElongation()*
												 (dynamic_cast<M3DQuadEndPrimitive *>(m2))->getElongation());
	}
	else
	{
		if(mNew==NULL)
			mNew=new M3DQuadPrimitive();
	}
	mNew->setX(m1->getX()+m2->getX());
	mNew->setR(m1->getR()*m2->getR());
//	mNew->setQ(m1->getQ()*m2->getQ());
	mNew->setQ(m2->getQ()*m1->getQ());
	//cout << "m2: norm" << m2->getQ().norm() << endl;
	//cout << "m1: norm" << m1->getQ().norm() << endl;

//	regularizeQuaternion(mNew);

	if(mNew->getQ().getW()<0)
	{
#ifdef debug
		cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
#endif
		Quat qTmp;
		qTmp=mNew->getQ();
		qTmp.neg();
		mNew->setQ(qTmp);
	}
	thetaTmp=m1->getTheta()+m2->getTheta();

	mNew->setTheta(thetaTmp);

	return true;
}

bool GeodesicSym::atomAverage(int mNum, M3DPrimitive **mList, M3DPrimitive *mNew, double *weights)
{
	int i;
	double delta=100.0;
	M3DPrimitive	*mInv=NULL,
					*mTmp=NULL;
	VectorND		vectSum,
					*vectTmp;
	bool			allEndPrims=true,
					_ignorePrimType=ignorePrimType;

	double *realWeights=new double[mNum];
	if(weights==NULL) {
		for(i=0; i<mNum; i++)
			realWeights[i]=1.0/mNum;
	}
	else {
		double totalWeights=0;
		for(i=0; i<mNum; i++) {
			// invalid weight here!
//			if(weights[i]<0)
//				return false;
			totalWeights+=weights[i];
		}
		for(i=0; i<mNum; i++)
			realWeights[i]=weights[i]/totalWeights;
	}

	for(i=0; i<mNum; i++)
	{
		if(mList[i]->type()!=M3D_END_PRIMITIVE)
		{
			allEndPrims=false;
			break;
		}
	}

	if(!ignorePrimType && allEndPrims)
	{
		if(mNew==NULL)
			mNew=new M3DQuadEndPrimitive();
		*mNew=*mList[0];

		if(mInv==NULL)
			mInv=new M3DQuadEndPrimitive();
		if(mTmp==NULL)
			mTmp=new M3DQuadEndPrimitive();

		vectSum.setSize(END_ATOM_LOG_MAP_VECTOR_SIZE);
		vectTmp=new VectorND(END_ATOM_LOG_MAP_VECTOR_SIZE);
	}
	else
	{
		if(mNew==NULL)
			mNew=new M3DQuadPrimitive();
		*mNew=*mList[0];

		if(mInv==NULL)
			mInv=new M3DQuadPrimitive();
		if(mTmp==NULL)
			mTmp=new M3DQuadPrimitive();

		vectSum.setSize(STANDARD_ATOM_LOG_MAP_VECTOR_SIZE);
		vectTmp=new VectorND(STANDARD_ATOM_LOG_MAP_VECTOR_SIZE);

		// in case of endPrimAvg==false, i.e. mixed prim types in mList
		// enforce the program to treat all atoms as standard ones
		ignorePrimType=true;
	}

	Vector2D *n0Log=NULL, *n1Log=NULL;
	//Vector2D uTempN0, uTempN1;
	Vector3D uTempN0, uTempN1;
	M3DSymPrimitive *sM=NULL;
	M3DSymPrimitive sMAvg;

	sM=new M3DSymPrimitive[mNum];
	sMAvg.x.set(0, 0, 0);
	sMAvg.r=0;
	sMAvg.elongation=0;

	n0Log=new Vector2D[mNum];
	n1Log=new Vector2D[mNum];
	uTempN0.set(0, 0, 0);
	uTempN1.set(0, 0, 0);
	for(i=0; i<mNum; i++) {
		sM[i] = atomToSymAtom(mList[i]);
		sMAvg.x+=sM[i].x*realWeights[i];
		sMAvg.r+=log(sM[i].r)*realWeights[i];
		sMAvg.elongation+=log(sM[i].elongation)*realWeights[i];

//		n0Log[i]=sphereLogMap(sM[i].n0);
//		n1Log[i]=sphereLogMap(sM[i].n1);
//		uTempN0+=n0Log[i]*realWeights[i];
//		uTempN1+=n1Log[i]*realWeights[i];
		uTempN0+=sM[i].n0;
		uTempN1+=sM[i].n1;
	}
	uTempN0.normalize();
	uTempN1.normalize();
	sMAvg.r = exp(sMAvg.r);
	sMAvg.elongation = exp(sMAvg.elongation);
	Quat qN0, qN1;
	qN0=rotationToOrigin(uTempN0);
	qN1=rotationToOrigin(uTempN1);
	qN0=qN0.conj();
	qN1=qN1.conj();


	Vector3D *rotatedData=NULL;
	rotatedData=new Vector3D[mNum];
	Vector2D logMean;
	Vector3D axis;
	double angle;

	delta = 1.0;
	// average up n0s
	while(delta>MEAN_ITERATION_THRESHOLD) {
		logMean.set(0, 0);
		for(i=0; i<mNum; i++) {
			rotatedData[i]=sM[i].n0;
			qN0.conj().rotateVector(rotatedData[i]);
			n0Log[i]=sphereLogMap(rotatedData[i]);
			logMean+=n0Log[i]*realWeights[i];
		}
		uTempN0=sphereExpMap(logMean);
		qN0.rotateVector(uTempN0);

		axis.set(0.0, -uTempN0.getZ(), uTempN0.getY());
		axis.normalize();
		angle=acos(uTempN0.getX());
		qN0.setAxisAngle(axis, angle);
		delta=logMean.norm();
	}
	sMAvg.n0.set(1, 0, 0);
	qN0.rotateVector(sMAvg.n0);

	// average up n1s
	delta = 1.0;
	while(delta>MEAN_ITERATION_THRESHOLD) {
		logMean.set(0, 0);
		for(i=0; i<mNum; i++) {
			rotatedData[i]=sM[i].n1;
			qN1.conj().rotateVector(rotatedData[i]);
			n1Log[i]=sphereLogMap(rotatedData[i]);
			logMean+=n1Log[i]*realWeights[i];
		}
		uTempN1=sphereExpMap(logMean);
		qN1.rotateVector(uTempN1);

		axis.set(0.0, -uTempN1.getZ(), uTempN1.getY());
		axis.normalize();
		angle=acos(uTempN1.getX());
		qN1.setAxisAngle(axis, angle);
		delta=logMean.norm();
	}
	sMAvg.n1.set(1, 0, 0);
	qN1.rotateVector(sMAvg.n1);

	symAtomToAtom(&sMAvg, mNew);

	ignorePrimType=_ignorePrimType;
	if(mInv!=NULL)
		delete mInv;
	if(mTmp!=NULL)
		delete mTmp;
	if(vectTmp!=NULL)
		delete vectTmp;
	if(sM!=NULL)
		delete[] sM;
	if(n0Log!=NULL)
		delete[] n0Log;
	if(n1Log!=NULL)
		delete[] n1Log;
	return true;
}


/*
 * If you modify this function, also modify
 * M3D[Quad|Tube]Primitive::atomInterp(...) accordingly.
 */
bool GeodesicSym::atomInterp(double t, M3DPrimitive *m1, M3DPrimitive *m2, M3DPrimitive *mNew)
{
	M3DPrimitive *mTmp=NULL;
	VectorND vTmp;

	if(!ignorePrimType && m1->type()==m2->type() && m1->type()==M3D_END_PRIMITIVE)
	{
		if(mNew==NULL)
			mNew=new M3DQuadEndPrimitive();
		mTmp=new M3DQuadEndPrimitive();
		vTmp.setSize(END_ATOM_LOG_MAP_VECTOR_SIZE);
	}
	else
	{
		if(mNew==NULL)
			mNew=new M3DQuadPrimitive();
		mTmp=new M3DQuadPrimitive();
		vTmp.setSize(STANDARD_ATOM_LOG_MAP_VECTOR_SIZE);
	}

	M3DSymPrimitive sM1, sM2, deltaSM, interpolatedSM;
	Vector3D rotatedSpoke;

	sM1=atomToSymAtom(m1);
	sM2=atomToSymAtom(m2);
	deltaSM.x=(sM2.x-sM1.x)*t;
	deltaSM.r=(log(sM2.r)-log(sM1.r))*t;
	deltaSM.elongation=(log(sM2.elongation)-log(sM1.elongation))*t;

	interpolatedSM.x=sM1.x+deltaSM.x;
	interpolatedSM.r=exp(log(sM1.r)+deltaSM.r);
	interpolatedSM.elongation=exp(log(sM1.elongation)+deltaSM.elongation);

	Quat q=rotationToOrigin(sM1.n0);
	q=q.conj();
	rotatedSpoke = sM2.n0;
	q.conj().rotateVector(rotatedSpoke);
	Vector2D logDataN0=sphereLogMap(rotatedSpoke);
	logDataN0*=t;
	interpolatedSM.n0=sphereExpMap(logDataN0);
	q.rotateVector(interpolatedSM.n0);

	q=rotationToOrigin(sM1.n1);
	q=q.conj();
	rotatedSpoke = sM2.n1;
	q.conj().rotateVector(rotatedSpoke);
	Vector2D logDataN1=sphereLogMap(rotatedSpoke);
	logDataN1*=t;
	interpolatedSM.n1=sphereExpMap(logDataN1);
	q.rotateVector(interpolatedSM.n1);

	symAtomToAtom(&interpolatedSM, mNew);

/*
	atomInverse(m1, mTmp);
	atomCompose(mTmp, m2, mTmp);
	atomLogMap(mTmp, &vTmp);
	vTmp*=t;
	atomExpMap(vTmp, mTmp);
	atomCompose(m1, mTmp, mNew);
*/
	if(m1->isSelected())
		mNew->select();
	else
		mNew->deselect();
	mNew->toggleHinge(m1->isHinge());

	if(mTmp!=NULL)
		delete mTmp;
	return true;
}

double GeodesicSym::atomDistance(M3DPrimitive *m1, M3DPrimitive *m2)
{
	return sqrt(atomSquareDistance(m1, m2));
}

#ifdef OLD_BINARY_PABLO

// Special version for Binary Pablo
double GeodesicSym::atomSquareDistance(M3DPrimitive *m1, M3DPrimitive *m2,
	double *xD, double *rD, double *eD, double *n0D, double *n1D)
{
	M3DPrimitive *mTmp=NULL;
	VectorND vTmp;
	double distance;

	if(!ignorePrimType && m1->type()==m2->type() && m1->type()==M3D_END_PRIMITIVE)
	{
		mTmp=new M3DQuadEndPrimitive();
		vTmp.setSize(END_ATOM_LOG_MAP_VECTOR_SIZE);
	}
	else
	{
		mTmp=new M3DQuadPrimitive();
		vTmp.setSize(STANDARD_ATOM_LOG_MAP_VECTOR_SIZE);
	}

	M3DSymPrimitive sM1, sM2, deltaSM;
	sM1=atomToSymAtom(m1);
	sM2=atomToSymAtom(m2);
	deltaSM.x=sM2.x-sM1.x;
	deltaSM.r = (log(sM2.r)-log(sM1.r)) * 0.03162;
	deltaSM.elongation=log(sM2.elongation)-log(sM1.elongation);

	Quat q=rotationToOrigin(sM1.n0);
	q.normalize();
	//Quat qTmp=rotationFromOrigin(sM1.n0);
	q=q.conj();
	(q.conj()).rotateVector(sM2.n0);
	Vector2D logDataN0=sphereLogMap(sM2.n0)*m1->getR();
	q=rotationToOrigin(sM1.n1);
	q.normalize();
	q=q.conj();
	(q.conj()).rotateVector(sM2.n1);
	Vector2D logDataN1=sphereLogMap(sM2.n1)*m1->getR();

	double sumSquare=(deltaSM.x.norm())*(deltaSM.x.norm());
	sumSquare+=deltaSM.r*deltaSM.r+deltaSM.elongation*deltaSM.elongation;
	sumSquare+=logDataN0.norm()*logDataN0.norm();
	sumSquare+=logDataN1.norm()*logDataN1.norm();

	distance=sumSquare;

	if (xD != NULL) {
		(*xD)  += (deltaSM.x.norm())*(deltaSM.x.norm());
		(*rD)  += deltaSM.r*deltaSM.r;
		(*eD)  += deltaSM.elongation*deltaSM.elongation;
		(*n0D) += logDataN0.norm()*logDataN0.norm();
		(*n1D) += logDataN1.norm()*logDataN1.norm();
	}

//	atomInverse(m1, mTmp);
//	atomCompose(mTmp, m2, mTmp);
//	atomLogMap(mTmp, &vTmp);
//	distance=vTmp.normSquare();

	if(mTmp!=NULL)
		delete mTmp;
	return distance;
}

#else
/*
 * If you modify this function, also modify
 * M3D[Quad|Tube]Primitive::atomDistance(...) accordingly.
 */
double GeodesicSym::atomSquareDistance(M3DPrimitive *m1, M3DPrimitive *m2, double *radius)
{
	double effectiveR = 1.0;
	if(radius == NULL)
		effectiveR = sqrt(m1->getR()*m2->getR());
	else
		effectiveR = *radius;

	M3DPrimitive *mTmp=NULL;
	VectorND vTmp;
	double distance;

	if(!ignorePrimType && m1->type()==m2->type() && m1->type()==M3D_END_PRIMITIVE) {
		mTmp=new M3DQuadEndPrimitive();
		vTmp.setSize(END_ATOM_LOG_MAP_VECTOR_SIZE);
	}
	else {
		mTmp=new M3DQuadPrimitive();
		vTmp.setSize(STANDARD_ATOM_LOG_MAP_VECTOR_SIZE);
	}

	M3DSymPrimitive sM1, sM2, deltaSM;
	sM1=atomToSymAtom(m1);
	sM2=atomToSymAtom(m2);
	if(!ignorePrimType && m1->type()==m2->type() && m1->type()==M3D_END_PRIMITIVE) {
		sM1.elongation = (dynamic_cast<M3DQuadEndPrimitive *>(m1))->getElongation();
		sM2.elongation = (dynamic_cast<M3DQuadEndPrimitive *>(m2))->getElongation();
	}

	deltaSM.x = sM2.x-sM1.x;
	deltaSM.r = (log(sM2.r)-log(sM1.r))*effectiveR;
	deltaSM.elongation = (log(sM2.elongation)-log(sM1.elongation))*effectiveR;

	Quat q=rotationToOrigin(sM1.n0);
	q.normalize();
//	q=q.conj();
//	(q.conj()).rotateVector(sM2.n0);
	q.rotateVector(sM2.n0);
	Vector2D logDataN0 = sphereLogMap(sM2.n0)*effectiveR;
//	logDataN0 = sphereLogMap(sM2.n0)*sqrt(m1->getR()*m2->getR());

	q=rotationToOrigin(sM1.n1);
	q.normalize();
//	q=q.conj();0
//	(q.conj()).rotateVector(sM2.n1);
	q.rotateVector(sM2.n1);
	Vector2D logDataN1 = sphereLogMap(sM2.n1)*effectiveR;
//	logDataN1 = sphereLogMap(sM2.n1)*sqrt(m1->getR()*m2->getR());

	double sumSquare = 0;
	sumSquare  = deltaSM.x.norm()*deltaSM.x.norm();
	sumSquare += deltaSM.r*deltaSM.r;
	sumSquare += deltaSM.elongation*deltaSM.elongation;
	sumSquare += logDataN0.norm()*logDataN0.norm();
	sumSquare += logDataN1.norm()*logDataN1.norm();

	distance=sumSquare;

//	atomInverse(m1, mTmp);
//	atomCompose(mTmp, m2, mTmp);
//	atomLogMap(mTmp, &vTmp);
//	distance=vTmp.normSquare();

	if(mTmp!=NULL)
		delete mTmp;
	return distance;
}

#endif

M3DSymPrimitive GeodesicSym::atomToSymAtom(M3DPrimitive *m) {
	M3DSymPrimitive sM;
	sM.x = m->getX();
	sM.r = m->getR();
	sM.n0 = m->getNormalizedY0();
	sM.n1 = m->getNormalizedY1();

	if(!ignorePrimType && m->type()==M3D_END_PRIMITIVE) {
		sM.elongation = (dynamic_cast<M3DQuadEndPrimitive *>(m))->getElongation();
	}

	return sM;
}

bool GeodesicSym::symAtomToAtom(M3DSymPrimitive *sM, M3DPrimitive *m)
{
	Quat q;
	double theta;

	m->setX(sM->x);
	m->setR(sM->r);

    Vector3D b, n, bPerp;

    // Choose n to be the difference of the two spokes
	// (n is in the same direction as Y1, not Y0?)

    if (sM->n1 == sM->n0) {
		printf("Two spokes identical in GeodesicSym::symAtomToAtom!\n");

        // Handle the case where n0 = n1
        if ( fabs(sM->n0.getX()) < fabs(sM->n0.getY()) )
            if ( fabs(sM->n0.getX()) < fabs(sM->n0.getZ()) )
                n.set(1.0, 0.0, 0.0);
            else
                n.set(0.0, 0.0, 1.0);
        else
            if ( fabs(sM->n0.getY()) > fabs(sM->n0.getZ()) )
                n.set(0.0, 0.0, 1.0);
            else
                n.set(0.0, 1.0, 0.0);
		n = n - (n * sM->n1) * sM->n1;
    }
    else
        n = sM->n1 - sM->n0;
    n.normalize();

    b = sM->n0 + sM->n1;
    // Handle the case where n0 = -n1
    if (b.norm() == 0) {
		printf("Two spokes back to back in GeodesicSym::symAtomToAtom!\n");

		// If n0 and n1 are back to back,
		// pick a set of local frames
        if ( fabs(sM->n0.getX()) < fabs(sM->n0.getY()) )
            if ( fabs(sM->n0.getX()) < fabs(sM->n0.getZ()) )
                b.set(1.0, 0.0, 0.0);
            else
                b.set(0.0, 0.0, 1.0);
        else
            if ( fabs(sM->n0.getY()) > fabs(sM->n0.getZ()) )
                b.set(0.0, 0.0, 1.0);
            else
                b.set(0.0, 1.0, 0.0);
    }
    b = b - (b * n) * n;
    b.normalize();

    bPerp = b.cross(n);

    theta = acos(b * sM->n0);
    q = Quat(b, n, bPerp);

	m->setTheta(theta);
	m->setQ(q);

	if (! ignorePrimType && m->type() == M3D_END_PRIMITIVE)
		(dynamic_cast<M3DQuadEndPrimitive *>(m))->setElongation(sM->elongation);

	return true;
}



