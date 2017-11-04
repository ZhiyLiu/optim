/*
 *
 */

#include "Geodesic.h"
#include "M3DFigure.h"

using namespace std;

#define MEAN_ITERATION_THRESHOLD			1.0e-04

#define BUG_FIX

/**
 * Returns the log map of the vector n which is a point on a unit sphere.
 */
Vector2D ShapeSpace::S2::Log(const Vector3D& n)
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
	if(d == 0)
		return Vector2D(0.0, 0.0);

	return (d/sin(d))*Vector2D(n.getY(), n.getZ());
}


Vector3D ShapeSpace::S2::Exp(const Vector2D& v)
{
	double d, alpha;

	d = v.norm();
	if(d == 0)
		return Vector3D(1.0, 0.0, 0.0);

	alpha = sin(d)/d;
	return Vector3D(cos(d), alpha * v.getX(), alpha * v.getY());
}

// another possible bug in M3DPGA.cpp?
Quat ShapeSpace::S2::rotationFromOrigin(const Vector3D& n)
{
	Vector2D logN;
	Vector3D axis;
	double angle;
	Quat q;

	logN = Log(n);
	axis.set(0.0, n.getZ(), -n.getY());
	axis.normalize();
	angle = logN.norm();

	q.setAxisAngle(axis, -angle);
	return q;
}

Quat ShapeSpace::S2::rotationToOrigin(const Vector3D& n)
{
	Vector2D logN;
	Vector3D axis;
	double angle;
	Quat q;

	logN = Log(n);
	axis.set(0.0, n.getZ(), -n.getY());
	axis.normalize();
	angle = logN.norm();

	q.setAxisAngle(axis, angle);
	return q;
}

Vector3D ShapeSpace::S2::mean(const Vector3D* points, const int mNum, const double* weights)
{
	int i;
	Vector2D pointLog;
	Vector3D rotatedPoint;
	Vector3D meanPt;
	Vector2D logMean;
	Vector3D axis;
	double angle;

	for( i = 0; i != mNum; ++i ) {
		meanPt	+= points[i];
	}
	meanPt.normalize();
	Quat qPt	= rotationToOrigin(meanPt);
	qPt	= qPt.conj();

	double delta	= 1.0;
	while( delta > MEAN_ITERATION_THRESHOLD ) {
		logMean.set( 0.0, 0.0 );

		for( i = 0; i != mNum; ++i ) {
			rotatedPoint	= points[i];
			qPt.conj().rotateVector( rotatedPoint );
			pointLog		= Log( rotatedPoint );
			if( weights == NULL ) {
				pointLog	*= 1.0/mNum;
			}
			else {
				pointLog	*= weights[i];
			}
			logMean	+= pointLog;
		}
		meanPt	= Exp(logMean);
		qPt.rotateVector( meanPt );

		axis.set( 0.0, -meanPt.getZ(), meanPt.getY() );
		axis.normalize();
		angle	= acos( meanPt.getX() );
		qPt.setAxisAngle( axis, angle );
		delta	= logMean.norm();
	};

	meanPt.set( 1.0, 0.0, 0.0 );
	qPt.rotateVector( meanPt );
	return meanPt;
}


/*
 * The following code is deprecated and no longer used.
 */
#ifndef BINARY

const double epsilon = 1.0e-03;
#define STANDARD_ATOM_LOG_MAP_VECTOR_SIZE	8
#define END_ATOM_LOG_MAP_VECTOR_SIZE		9
#define NUMERICAL_ZERO 1.0e-10

//#define debug

/* Exponential and Log map for atoms
   1.	Only the atom (X, R, Q, Theta) are updated here.
	Other fields including selected/type/hinge should be updated somewhere
	else.
   2.	All atoms are considered as standard atoms by default.
	Special treatment is needed for end atoms by instantialize Geodesic
	with 'FALSE'.
*/
bool Geodesic::atomLogMap(M3DPrimitive *m, VectorND *v)
{
	// v has different sizes for different atom type
	if(!ignorePrimType && typeid(*m) == typeid(M3DQuadEndPrimitive))
	{
		if(v==NULL)
			v=new VectorND(END_ATOM_LOG_MAP_VECTOR_SIZE);
		v->set(8, log((dynamic_cast<M3DEndPrimitive *>(m))->getElongation()));
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
}

bool Geodesic::atomExpMap(VectorND &v, M3DPrimitive *mNew)
{
	Vector3D logQ;
	Quat q;

	// to treat standard/end atoms differently
	if(!ignorePrimType && v.getSize()==END_ATOM_LOG_MAP_VECTOR_SIZE)
	{
		if(mNew==NULL)
			mNew=new M3DQuadEndPrimitive();
		(dynamic_cast<M3DEndPrimitive *>(mNew))->setElongation(exp(v.get(8)));
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

bool Geodesic::atomInverse(M3DPrimitive *m, M3DPrimitive *mNew)
{
	if(!ignorePrimType && typeid(*m) == typeid(M3DQuadEndPrimitive))
	{
		if(mNew==NULL)
			mNew=new M3DQuadEndPrimitive();
		(dynamic_cast<M3DEndPrimitive *>(mNew))->setElongation(1.0/
			 (dynamic_cast<M3DEndPrimitive *>(m))->getElongation());
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

void regularizeQuaternion(M3DPrimitive *mNew)
{
	Quat qTmp;
	qTmp=mNew->getQ();
	if(qTmp.getW()<0 && fabs(qTmp.getW())<NUMERICAL_ZERO)
		qTmp.setW(0);
	mNew->setQ(qTmp);
}

bool Geodesic::atomCompose(M3DPrimitive *m1, M3DPrimitive *m2, M3DPrimitive *mNew)
{
	double thetaTmp;
	if(!ignorePrimType && typeid(*m1) == typeid(*m2) && typeid(*m1) == typeid(M3DQuadEndPrimitive))
	{
		if(mNew==NULL)
			mNew=new M3DQuadEndPrimitive();
		(dynamic_cast<M3DEndPrimitive *>(mNew))->setElongation((dynamic_cast<M3DEndPrimitive *>(m1))->getElongation()*
												 (dynamic_cast<M3DEndPrimitive *>(m2))->getElongation());
	}
	else
	{
		if(mNew==NULL)
			mNew=new M3DQuadPrimitive();
	}
	mNew->setX(m1->getX()+m2->getX());
	mNew->setR(m1->getR()*m2->getR());
	mNew->setQ(m1->getQ()*m2->getQ());

	regularizeQuaternion(mNew);

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

bool Geodesic::atomAverage(int mNum, M3DPrimitive **mList, M3DPrimitive *mNew)
{
	int i;
	double delta=1.0;
	M3DPrimitive	*mInv=NULL,
					*mTmp=NULL;
	VectorND		vectSum,
					*vectTmp;
	bool			allEndPrims=true,
					_ignorePrimType=ignorePrimType;

	for(i=0; i<mNum; i++)
	{
		if(typeid(*mList[i]) != typeid(M3DQuadEndPrimitive))
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

	while(delta>MEAN_ITERATION_THRESHOLD)
	{
		vectSum.set(0.0);
#ifdef debug
		vectSum.print();
		mNew->print();
#endif

		atomInverse(mNew, mInv);

#ifdef debug
		mInv->print();
#endif

		for(i=0; i<mNum; i++)
		{
			atomCompose(mInv, mList[i], mTmp);
#ifdef debug
			mTmp->print();
#endif

			atomLogMap(mTmp, vectTmp);

#ifdef debug
			vectTmp->print();
#endif

			vectSum+=*vectTmp;

#ifdef debug
			vectSum.print();
#endif
		}
		vectSum/=mNum;
		atomExpMap(vectSum, mTmp);
		atomCompose(mNew, mTmp, mNew);
		atomLogMap(mTmp, vectTmp);
		delta=vectTmp->norm();
		// here vectSum.norm() should be equal to delta
	}

	ignorePrimType=_ignorePrimType;
	if(mInv!=NULL)
		delete mInv;
	if(mTmp!=NULL)
		delete mTmp;
	if(vectTmp!=NULL)
		delete vectTmp;
	return true;
}

bool Geodesic::atomInterp(double t, M3DPrimitive *m1, M3DPrimitive *m2, M3DPrimitive *mNew)
{
	M3DPrimitive *mTmp=NULL;
	VectorND vTmp;

	if(!ignorePrimType && typeid(*m1) == typeid(*m2) && typeid(*m1) == typeid(M3DQuadEndPrimitive))
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

	atomInverse(m1, mTmp);
	atomCompose(mTmp, m2, mTmp);
	atomLogMap(mTmp, &vTmp);
	vTmp*=t;
	atomExpMap(vTmp, mTmp);
	atomCompose(m1, mTmp, mNew);

	if(mTmp!=NULL)
		delete mTmp;
	return true;
}

double Geodesic::atomDistance(M3DPrimitive *m1, M3DPrimitive *m2)
{
	return sqrt(atomSquareDistance(m1, m2));
}

double Geodesic::atomSquareDistance(M3DPrimitive *m1, M3DPrimitive *m2)
{
	M3DPrimitive *mTmp=NULL;
	VectorND vTmp;
	double distance;

	if(!ignorePrimType && typeid(*m1)==typeid(*m2) && typeid(*m1) == typeid(M3DQuadEndPrimitive))
	{
		mTmp=new M3DQuadEndPrimitive();
		vTmp.setSize(END_ATOM_LOG_MAP_VECTOR_SIZE);
	}
	else
	{
		mTmp=new M3DQuadPrimitive();
		vTmp.setSize(STANDARD_ATOM_LOG_MAP_VECTOR_SIZE);
	}

	atomInverse(m1, mTmp);
	atomCompose(mTmp, m2, mTmp);
	atomLogMap(mTmp, &vTmp);
	distance=vTmp.normSquare();

	if(mTmp!=NULL)
		delete mTmp;
	return distance;
}

#endif

