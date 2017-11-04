#include "GeodesicDistanceFunction.h"
#include "ConjugateGradientMethod.h"
#include <iostream>

using std::cout;
using std::endl;
using std::vector;


GeodesicDistanceFunction::GeodesicDistanceFunction(
	primList ref, primList target, SimTransShapeSpace* space)
{
	refPrims = ref;
	targetPrims = target;
	shapeSpace = space; 
	transWeight = 1;
	scaleWeight = 1;	
	quatWeight = 1;
	thetaWeight = 1;
}

SimilarityTransform3D GeodesicDistanceFunction::createSimTrans(const Vector &v)
{
	return shapeSpace->createSimTrans(v);
}

vector<M3DPrimitive*> GeodesicDistanceFunction::createCandidateAtoms(const Vector& v)
{
	SimilarityTransform3D xform = createSimTrans(v);
	primList candPrims;
	primList::iterator refIter = refPrims.begin();
	for ( ; (refIter != refPrims.end()); refIter++)
	{
		M3DPrimitive *candPrim = 0;

		if (*refIter) {

			candPrim = (*refIter)->copyPtr();
			// rotate
			candPrim->rotateBy(xform.getRotation());
			Vector3D temp = candPrim->getX();
			xform.getRotation().rotateVector(temp);
			candPrim->setX(temp);

			//scale
			candPrim->scaleBy(xform.getScale());
			candPrim->setX(xform.getScale() * candPrim->getX());

			//translate
			candPrim->setX(xform.getTranslation() + candPrim->getX());
		} // if(*refIter)
		candPrims.push_back(candPrim);
	}
	return candPrims;
}
	
double GeodesicDistanceFunction::evaluate(const Vector& v)
{
	double dist = 0;
	primList candPrims = createCandidateAtoms(v);
	primList::iterator targetIter = targetPrims.begin(); 
	primList::iterator candIter = candPrims.begin(); 

	for ( ;
		(candIter != candPrims.end()) && (targetIter != targetPrims.end());
		candIter++, targetIter++)
	{

		if ((*candIter == 0) || (*targetIter == 0)) {
			continue;
		}

		M3DPrimitive* candPrim = *candIter;
		M3DPrimitive* targetPrim = *targetIter;

		double thisDist = 0;
		thisDist += transWeight * pow((candPrim->getX() - targetPrim->getX()).norm(), 2);

		// dibyendu

		// if this is a quadPrimitive, then take an average of the difference of the two radii separately

		if( ( dynamic_cast<M3DQuadPrimitive *>(candPrim) ) != NULL && 
			( dynamic_cast<M3DQuadPrimitive *>(targetPrim) ) != NULL ) {

			double dLogR0 = log( (dynamic_cast<M3DQuadPrimitive *>(candPrim))->getR0() ) - 
							log( (dynamic_cast<M3DQuadPrimitive *>(targetPrim))->getR0() ) ;

			double dLogR1 = log( (dynamic_cast<M3DQuadPrimitive *>(candPrim))->getR1() ) - 
							log( (dynamic_cast<M3DQuadPrimitive *>(targetPrim))->getR1() ) ;

			thisDist += scaleWeight * pow( 0.5 * (dLogR0 + dLogR1), 2);
		}
		else
			thisDist += scaleWeight * pow(log(candPrim->getR()) - log(targetPrim->getR()), 2);

		thisDist += quatWeight * ((candPrim->getQ().logMap()) - (targetPrim->getQ().logMap())).norm();

		thisDist += thetaWeight * pow(cos(candPrim->getTheta()) - cos(targetPrim->getTheta()), 2);

		dist += thisDist;
		//  cout << "Evaluate: prim: " << i << " val is " << dist << endl << endl;
	}
	return dist;
}

// This really invokes conjugate gradient
Vector GeodesicDistanceFunction::doOptimization() 
{
	int numParams = shapeSpace->parameterCount();
	Vector start(numParams);
	for (int idx = 0; idx < numParams; idx++) {
		start(idx) = 0;
	}

	NumericalFunction wrapper(*this, 1e-6);

	int iter = 0;
	ConjugateGradientMethod opt(wrapper, start);
	while (! opt.isFinished()) {
		Vector best = opt.getBestEverX();
		double v = opt.getBestEverValue();
		cout << "Iter: " << iter << " bestX: ";
		best.print();
		cout << " bestVal: " << v << endl;
		opt.performIteration();
		iter++;
	}
	Vector best = opt.getBestEverX();
	return best;
}

