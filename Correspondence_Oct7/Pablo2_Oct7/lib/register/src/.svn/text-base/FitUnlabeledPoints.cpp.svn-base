#include <cstdlib>
#include <vector>
#include <fstream>
#include <iostream>
#include <sstream>
#include "matrix.h"
#include "DistanceToPointSetFunction.h"
#include "M3DObject.h"
#include "ConjugateGradientMethod.h"
#include "SimilarityTransform3D.h"
#include "SimTransShapeSpace.h"
#include "M3DPGAStats.h"
#include "FitUnlabeledPoints.h"

using std::vector;
using std::ios;


FitUnlabeledPoints::FitUnlabeledPoints(int level, int figureId,
	int order, double closeness, double outerLimit, double innerLimit)
{
	closenessLimit = closeness;
	figure_id = figureId;
	pg_order = order;
	convergence_limit = outerLimit;
	fit_limit = innerLimit;
	surf_level = level;
	bestObj = NULL;
	pga_coefs = NULL;
	xform = NULL;
}

FitUnlabeledPoints::~FitUnlabeledPoints()
{
	if (xform)
		delete xform;
	if (bestObj)
		delete bestObj;
	if (pga_coefs)
		delete pga_coefs;
}

bool FitUnlabeledPoints::fit(M3DObject * model, const char * filename,
	Image3D * image, PointFitType type)
{
	DistanceToPointSetFunction::SimTransType fitType;

    if (model == NULL)
		return false;

    if (filename == NULL || *filename == '\0')
		return false;

    // Configure the distance function
	if (type == ScaleTrans)
		fitType = DistanceToPointSetFunction::Both;
	else if (type == Trans)
		fitType = DistanceToPointSetFunction::Translation;
	else if (type == Scale || type == ScalePGA)
		fitType = DistanceToPointSetFunction::Scale;
	else	// if (type == Similarity || type = SimPGA)
		fitType = DistanceToPointSetFunction::Full;
	// Note that fitType == Translation is not currently being used
    DistanceToPointSetFunction f(filename, fitType, surf_level,
		figure_id, image);
	return runPointsFit(f, model, image, type);
}

bool FitUnlabeledPoints::fit(M3DObject * model, double * x, double * y,
	double * z, int numPoints, Image3D * image, PointFitType type)
{
	int i;
	DistanceToPointSetFunction::SimTransType fitType;

    if (model == NULL)
		return false;

    // Configure the distance function
	Vector3D * points = new Vector3D[numPoints];
	for (i = 0; i < numPoints; i++)
		points[i] = Vector3D(x[i], y[i], z[i]);
	if (type == ScaleTrans)
		fitType = DistanceToPointSetFunction::Both;
	else if (type == Trans)
		fitType = DistanceToPointSetFunction::Translation;
	else if (type == Scale || type == ScalePGA)
		fitType = DistanceToPointSetFunction::Scale;
	else	// if (type == Similarity || type = SimPGA)
		fitType = DistanceToPointSetFunction::Full;
	// Note that fitType == Translation is not currently being used
    DistanceToPointSetFunction f(points, numPoints, fitType, surf_level,
		figure_id, image);
	delete [] points;
	return runPointsFit(f, model, image, type);
}

bool FitUnlabeledPoints::runPointsFit(DistanceToPointSetFunction & f,
	M3DObject * model, Image3D * image, PointFitType type)
{
	int i, numPGs;
    int iter;
    double finalScore, lastScore;

	f.setModel(model);
    f.setClosenessLimit(closenessLimit);
    f.setPGAOrder(pg_order);
	if (pg_order >= 0) {
		numPGs =
			(pg_order >= 0) ? model->getPGAStats()->getPGDataPtr(pg_order)->numPGs : 0;
		if (pga_coefs)
			delete pga_coefs;
		pga_coefs = new Vector(numPGs);   // Create a zeroed vector
		f.setPGACoefs(*pga_coefs);
    }
    f.setLevel(surf_level);

    int numParams = f.parameterCount();

    Vector postSimTransCoefs(numParams);
    for (i = 0; i < numParams; i++)
		postSimTransCoefs(i) = 0.0;
    f.setPostSimTrans(postSimTransCoefs);   // Sets the length of the vector in f

    if (globalVerbosity >= 1) {
		cout << "Point fitting parameters:\n";
		cout << "  Level = " << surf_level << '\n';
		cout << "  Closeness = " << closenessLimit << '\n';
		cout << "  PG order = " << pg_order << '\n';
		cout << "  Number of PG coefficients = " << numPGs << '\n';
		cout << "  Number of parameters = " << numParams << endl;
    }

    Vector best;
    NumericalFunction wrapper(f, fit_limit);

    double simTransScore = 1.0e9; // Value of objective function after simtrans
    // This is used as stopping criteria for alternating the optimization between
	// fitting by similarity transformation and fitting by PGA deformation.  In
	// particular, when the PGA fit doesn't improve on result of the similarity
	// transformation fit, stop.

    if (type >= PGA)	// Both PGA and sim trans being done
        finalScore = 1.0e9;	// Force the while condition below to fail
    else
        finalScore = 0.0;

    int grandIter = 0;

    // Alternate between optimizing over similarity transform and PGA coefficients
    do {
        lastScore = finalScore;
        if (type != PGA) {
		    // Similarity transform iteration
		    f.setStage(DistanceToPointSetFunction::SimT);
		    iter = 0;
		    ConjugateGradientMethod opt(wrapper, postSimTransCoefs);
		    while (! opt.isFinished()) {
			    best = opt.getBestEverX();
			    double v = opt.getBestEverValue();
			    if (globalVerbosity >= 1) {
				    cout << "[POST-ST] Iter: " << grandIter << ", " << iter << " bestX: ";
				    best.print();
				    cout << " bestVal: " << v << endl;
			    }
			    opt.performIteration();
			    iter++;
		    }
		    postSimTransCoefs = opt.getBestEverX();
		    finalScore = opt.getBestEverValue();
		    simTransScore = finalScore;
		    f.setPostSimTrans(postSimTransCoefs);

		    if (globalVerbosity >= 1)
			    cout << "Iter: " << iter << " bestX: " << postSimTransCoefs(0)
				    << " bestVal: " << finalScore << endl;
		    else if (globalVerbosity == 0)
			    cout << "Result after " << iter << " sim trans iterations = " << finalScore << endl;
        }

		if (type >= PGA && pg_order >= 0) {  // PGA iteration
			f.setStage(DistanceToPointSetFunction::PGA);
			iter = 0;
			ConjugateGradientMethod opt(wrapper, *pga_coefs);
			while (! opt.isFinished()) {
				best = opt.getBestEverX();
				double v = opt.getBestEverValue();
				if (globalVerbosity >= 1) {
					cout << "[PG] Iter: " << grandIter << ", " << iter << " bestX: ";
					best.print();
					cout << " bestVal: " << v << endl;
				}
				opt.performIteration();
				iter++;
			}
			*pga_coefs = opt.getBestEverX();
			finalScore = opt.getBestEverValue();
			if (globalVerbosity >= 1)
				cout << "Iter: " << iter << " bestX: " << (*pga_coefs)(0) << " bestVal: "
					<< finalScore  << endl;
			else if (globalVerbosity == 0)
				cout << "Result after " << iter << " PGA iterations = " << finalScore << endl;
			f.setPGACoefs(*pga_coefs);
		}

		grandIter++;

		if (type >= PGA && globalVerbosity >= 1)
			// Both PGA and sim trans being done
			cout << "[ALTERNATION] difference is: " << (simTransScore - finalScore)
				<< endl;

    } while ((lastScore - finalScore) > convergence_limit);

	// Save the computed optimal transform for function result()
	if (! xform)
	    xform = new SimilarityTransform3D;
	*xform = f.getSpace()->createSimTrans(postSimTransCoefs);
	// Save the PGA-deformed model for function result()
	if (type != PGA || pg_order >= 0)									// AGG: Why not always do it?
		bestObj = f.createTargetObject(best);

    return true;
}

void FitUnlabeledPoints::result(SimilarityTransform3D & outputXform)
{
	outputXform = *xform;
}

void FitUnlabeledPoints::result(std::vector<double> & pgaCoefs)
{
	pgaCoefs.clear();
	if (pga_coefs != NULL) {
		for (int i = 0; i < pga_coefs->size(); i++)
			pgaCoefs.push_back((*pga_coefs)(i));
	}
}

