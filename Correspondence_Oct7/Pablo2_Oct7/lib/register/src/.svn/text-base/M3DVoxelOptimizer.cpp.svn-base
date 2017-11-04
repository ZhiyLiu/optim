#include <fstream>
#include <math.h>
#include <time.h>
#include "OptimizerBase.h"
#include "M3DVoxelOptimizer.h"

//#define DEBUG


const double MIN_PENALTY_WEIGHT = 0;
const double MIN_PENALTY_WEIGHT_EXP = 0.0;
const double MAX_PENALTY_WEIGHT_EXP = 10.0;
const double PENALTY_WEIGHT_POWER = 0.4;

const double TOLERENCE = 1e-4; //1e-5;

using namespace std;
using namespace ThallCode;

M3DVoxelOptimizer::M3DVoxelOptimizer()
{
	type = BoundaryDispl;

	targetObject = NULL;
	surfaceLevel = 0;
	figureId = 0;
	finished = false;

	numPts = 0;
	Bpnt = NULL;
	neighbors = NULL;
	match = NULL;

	cumChange = 0.0;
	lastPenalty = 0.0;
	lastMatchValue = 0.0;

    // Hardcoded to match initial slider value
    setPenaltyWeight(0.5);
    setConstraintsPenaltyWeight(0.5);

	// parameters
	sigma1 = 0.3;
	sigma2 = 1.0;
}

M3DVoxelOptimizer::~M3DVoxelOptimizer()
{
	clearWorkSpace();
	deleteTargetObject();
}

void M3DVoxelOptimizer::clearWorkSpace()
{
	if (Bpnt != NULL) 
	{
		delete [] Bpnt;
		Bpnt = NULL;
	}

	if (neighbors != NULL)
	{
		delete [] neighbors;
		neighbors = NULL;
	}
}

void M3DVoxelOptimizer::deleteTargetObject()
{
	if (targetObject != NULL) {
		delete targetObject;
	targetObject = NULL;
	}
}

void M3DVoxelOptimizer::setPenaltyWeight(double w)
{
	penaltyWeight = w;
}

void M3DVoxelOptimizer::setConstraintsPenaltyWeight(double w)
{
    double penaltyExp = MIN_PENALTY_WEIGHT_EXP +
        (MAX_PENALTY_WEIGHT_EXP - MIN_PENALTY_WEIGHT_EXP) * (w - 1);
    constraintsPenaltyWeight = MIN_PENALTY_WEIGHT + pow(10.0, penaltyExp);
	//cout << "Constraints penalty weight is " << constraintsPenaltyWeight << endl;
}

M3DObject * M3DVoxelOptimizer::getReferenceObject()
{
	if (match != NULL)
		return match->getReferenceObject();
	else
		return NULL;
}

bool M3DVoxelOptimizer::initialize(Match * matchPtr, M3DObject * object, 
												  int figId, int level, bool useWindowing)
{
#ifdef DEBUG
	cout << "M3DVoxelOptimizer::initialize()\n";
#endif
	figureId = figId;
	surfaceLevel = level;
	windowIntensities = useWindowing;

	match = matchPtr;
    if (getReferenceObject() == NULL) {
		cout << "Match not initialized.  Quitting BoundaryOpt::init()." << endl;
        return false;
	}

	dmask = match -> getDMask(figureId);
	dmask->setDmaskType(match->getType());
	// The dmask actually contains all the types  of filters (currently
	// Gaussian and NegGaussian).  A member array of dmask contains
	// the preferred filter type for every point that will be called for in
	// Optimization.  As of 3/03, I'm thinking level 4 all the time.

	if (dmask == NULL || dmask->getMaskSize() == 0)
	{
		cout << "Error: selected template type not available." << endl;
		return false;
	}

	targetImage = match -> getTargetImage();
	if (targetImage == NULL)
	{
		cout << "Error: no image information." << endl;
		return false;
	}

	if (object == NULL) {
		cout << "Error: no object information in BDOpt::init()." << endl;
		return false;
	}

	deleteTargetObject();
	targetObject = object->assign();

	if (targetObject == NULL) {
		cout << "Error: null object(s) in BOpt:init." << endl;
		return false;
	}

	initializeBoundary(targetObject);

	M3DQuadFigure * figure = dynamic_cast<M3DQuadFigure*>(targetObject -> getFigurePtr(figureId));
	if (figure == NULL) {
		cout << "Error: null figure in BOpt::init()." << endl;
		return false;
	}

	dmask->initializeSurfaceEls(figure, level);

	boundary = figure -> getBoundaryPtr();
	if (boundary == NULL)
	{
		cout << "Error: null boundary in BdOpt::init." << endl;
		return false;
	}

	boundary -> initDisplacements(surfaceLevel, numPts);
	Displacements * dsp = boundary -> getDisplacements(surfaceLevel);
	if (dsp == NULL)
	{
		cout << "Error: null displacement info in BdOpt::init." << endl;
		return false;
	}

	dVals = dsp -> getVals();
	if (dVals == NULL)
	{
		cout << "Error: null displacement array in BdOpt::init." << endl;
		return false;
	}

	finished = false;

	return true;
}

void M3DVoxelOptimizer::initializeBoundary(M3DObject * object)
{
	if (object == NULL) {
		std::cout << "Null object in SubdivBoundary::init." << endl;
		return;
	}

	M3DFigure * figure = object -> getFigurePtr(figureId);
	if (figure == NULL) {
		std::cout << "Null figure in BOpt::initBd." << endl;
		return;
	}
	Xferlist * xferList = convertM3DtoXfer(figure);

	clearWorkSpace();

	Pointlist_server2 * pList = new Pointlist_server2;
	Tileneighbors * vn = NULL;

	pList -> init(xferList);
	pList -> ComputeSubdivPointCloud(surfaceLevel);
	pList -> copySubdivBoundaryInfo(&numPts, &Bpnt);

	int num_n = 0;
	pList -> subdivvertexneighbors(&num_n, &vn);

	if (num_n == numPts && vn != NULL) 
	{
		// Get neighbor list for each vertex:
		neighbors = new NeighborList[numPts];
		int i, j, tmp;

		for (i=0; i<numPts; i++)
		{
			for (j=0; j<MAX_NEIGHBOR_COUNT; j++)
			{
				tmp = vn[i][j];
				if (tmp >= 0 && tmp<numPts) //if (tmp >= 0)
					neighbors[i][j] = tmp;
				else
				{
					if(tmp>=numPts)
						cout << "Invalid vertex neighbor index " << tmp \
						<< " in M3DVoxelOptimizer::initializeBoundary(...)" << endl;
					break;
				}	// 2003/07/07
			}
			neighbors[i].set_degree(j);
			neighbors[i][j] = neighbors[i][0];
		}
	}
	else
		cout << "Error in vertex neighbor computation in M3DBdOpt::initBd." << endl;

	delete[] (xferList -> atomlist);
	delete xferList;
	delete pList;
}

bool M3DVoxelOptimizer::performIterations(int nIterations)
{
	int i;
	cumChange = 0.0;

#ifdef DEBUG
	cout << "M3DVoxelOptimizer::performIterations()\n";
#endif
	if (Bpnt == NULL || neighbors == NULL)
	{
		cout << "Error: no boundary information in BdOpt::performIts." << endl;
		return false;
	}

	// Weights
	double alpha = penaltyWeight;

	// Alpha_one is a scalar to make the new values fit in the optimization that was
	// already being done.  The beta is so that the steps taken are a decent size.
	double alpha_one = 10000.0;
	double beta = .001;  // hard coded, 0.01

	double cutoff = dmask->getCutoff();
	int dmaskSize = dmask->getMaskSize();
	double dsize = 2 * cutoff / dmaskSize;

	double *profile = new double[dmaskSize + 2];
	double sum, sum_sq, rms, mu;


	for (int currentIteration = 0; currentIteration < nIterations; currentIteration++)
	{
		double delx, dely, delz, dv, val1;
		int j;

		// Create a random permutation
		int *perm = new int[numPts];
		int index, temp;
		for(i = 0; i < numPts; i++)
			perm[i] = i;
		for(i = 0; i < numPts; i++)
		{
			index = i + (rand() % (numPts - i));
			temp = perm[i];
			perm[i] = perm[index];
			perm[index] = temp;
		}

		for (index = 0; index < numPts; index++)
		{
			i = perm[index];

			Vector3D Pi(Bpnt[i].pnt);
			Vector3D Ni(Bpnt[i].norm);
			double r = Bpnt[i].rad;

			DerivMaskElement * delm;
			double area = dmask -> getSurfaceElement(surfaceLevel, i) -> area;
			Vector3D coord;
			double mVal = 0.0, 
				   pVal = 0.0;

			// We need to use mean zero, rms 1 intensity values along the
			// profile, so that the responses calculated are all in [-1, 1].
			// We don't actually need the mean zero, since the filter (dmask)
			// is already that.

			sum = 0.0; sum_sq = 0.0;
			for (j = 0; j < dmaskSize; j++)
			{
				delm = dmask -> getMaskElement(i,j);
				coord = Pi + Ni * (dVals[i] + delm->dist) * r;
				targetImage->modelToImageCoordinates(coord);

				if (windowIntensities)
					profile[j] = targetImage->getWindowedInterpolatedVoxelValue(
						coord.getX(), coord.getY(), coord.getZ());
				else
					profile[j] = targetImage->getInterpolatedVoxelValue(
						coord.getX(), coord.getY(), coord.getZ());

				sum += profile[j];
				sum_sq += profile[j]*profile[j];

			}

			coord = Pi + Ni * (dVals[i] - cutoff) * r;
			targetImage->modelToImageCoordinates(coord);
			if (windowIntensities)
				profile[dmaskSize] = targetImage->getWindowedInterpolatedVoxelValue(
					coord.getX(), coord.getY(), coord.getZ());
			else
				profile[dmaskSize] = targetImage->getInterpolatedVoxelValue(
					coord.getX(), coord.getY(), coord.getZ());

			sum += profile[dmaskSize];
			sum_sq += profile[dmaskSize]*profile[dmaskSize];

			coord = Pi + Ni * (dVals[i] + cutoff) * r;
			targetImage->modelToImageCoordinates(coord);
			if (windowIntensities)
				profile[dmaskSize+1] = targetImage->getWindowedInterpolatedVoxelValue(
					coord.getX(), coord.getY(), coord.getZ());
			else
				profile[dmaskSize+1] = targetImage->getInterpolatedVoxelValue(
					coord.getX(), coord.getY(), coord.getZ());

			sum += profile[dmaskSize+1];
			sum_sq += profile[dmaskSize+1]*profile[dmaskSize+1];

			rms = sum_sq - (1.0/(dmaskSize+2))*sum*sum;
			rms = rms/(dmaskSize+1);
			// We divide the computed response with this before
			// determining the motion of optimization.
			if(rms <= R_SMALL_TOLERANCE)
				rms = 0.0;
			else
				rms = sqrt(rms);

			//if (rms < 10) 
			//	double qq = 0;

			//if (dmask->getTemplateType(i) == 2)
			//	int qq = 0;

			for (j = 0; j < dmaskSize; j++)
			{
				delm = dmask -> getMaskElement(i,j);

				val1 = delm->dTempVal;

				mVal += -val1*profile[j];
			}

			mVal *= dsize*area;

			val1 = dmask->getCutoffValue(0,i);
			mVal += -val1 * profile[dmaskSize] * area;

			val1 = dmask->getCutoffValue(1,i);
			mVal += val1 * profile[dmaskSize+1] * area;
			mVal /= rms;

			// If the rms/mean(profile) is small, then we're in a very flat
			// area and we shouldn't guess on movement.

			mu = sum/(dmaskSize+2.0);
			mu = rms/mu;

			if (mu < .01)
			{
				if (mu > .004)
					mVal *= (1.0 - exp(-200000*(mu-.004)*(mu-.004)));
				else
					mVal = 0;
			}

			pVal += -dVals[i] / sigma1 / sigma1;
			int c = neighbors[i].get_degree();
			for (j = 0; j < c; j += 2)
			{
				int ij = neighbors[i][j];
				delx = Bpnt[i].pnt[0] - Bpnt[ij].pnt[0];
				dely = Bpnt[i].pnt[1] - Bpnt[ij].pnt[1];
				delz = Bpnt[i].pnt[2] - Bpnt[ij].pnt[2];
				dv   = delx*delx + dely*dely + delz*delz;

				pVal += -(dVals[i]*r - dVals[ij]*Bpnt[ij].rad) / (sigma2*sigma2*dv);
			}
			pVal /= numPts;

			double change = beta*((1-alpha)*alpha_one*mVal + alpha*pVal) ;

			// If the rms for mVal normalization was too small, then
			// don't move at all.
			if (rms <= R_SMALL_TOLERANCE)
				change = beta*(alpha*pVal);

			dVals[i] += change;

			cumChange = max(cumChange, fabs(change));
		}

		delete [] perm;

	}

	lastObjectiveFunctionValue = (1 - penaltyWeight) * computeLastMatchValue() + 
									penaltyWeight * computeLastPenaltyValue();

	if (cumChange < TOLERENCE)
		finished = true;
	return true;
}

double M3DVoxelOptimizer::computeLastMatchValue()
{
	int i, j;
	double imgVal;

	if (Bpnt == NULL)
	{
		cout << "Error: no boundary information in BdOpt::compMatch." << endl;
		return 0.0;
	}

	int dmaskSize = dmask->getMaskSize();
	double cutoff = dmask->getCutoff();
	lastMatchValue = 0.0;

	double mu = 0.0, totval2Sqr = 0.0, rms;
//	cout << numPts;
	for (i = 0; i < numPts; i++)
	{
		DerivMaskElement * delm;
		double area = dmask->getSurfaceElement(surfaceLevel, i)->area;
		double r = Bpnt[i].rad;
		double vol = 2 * area * r * cutoff / dmaskSize;
		Vector3D Pi(Bpnt[i].pnt);
		Vector3D Ni(Bpnt[i].norm);

		for (j = 0; j < dmaskSize; j++)
		{
			delm = dmask -> getMaskElement(i, j);
			Vector3D coord = Pi + Ni * (dVals[i] + delm->dist) * r;
			targetImage->modelToImageCoordinates(coord);

			if (windowIntensities)
				imgVal = targetImage->mapDisplayToRelative(
					targetImage->getWindowedInterpolatedVoxelValue(
					coord.getX(), coord.getY(), coord.getZ()));
			else
				imgVal = targetImage->mapDisplayToRelative(
					targetImage->getInterpolatedVoxelValue(
					coord.getX(), coord.getY(), coord.getZ()));

			mu += imgVal * vol;
			totval2Sqr += imgVal * imgVal * vol;
		}
	}

	double maskVolume = dmask -> getTotalVolume();
//	cout << "\nvol: " << maskVolume << ' '; 
    mu /= maskVolume;
    rms = (totval2Sqr / maskVolume) - (mu * mu);

//	cout << "\nmu: " << mu << ", val2sqr: " << totval2Sqr << ", rms: " << rms << '\n';
    if(rms <= R_SMALL_TOLERANCE)
        rms = 0.0;
    else
        rms = sqrt(rms);

    if(rms <= R_SMALL_TOLERANCE)
        return 0.0;

	for (i = 0; i < numPts; i++)
	{
		DerivMaskElement * delm;
		double area = dmask->getSurfaceElement(surfaceLevel, i)->area;
		double r = Bpnt[i].rad;
		double vol = 2 * area * r * cutoff / dmaskSize;
		Vector3D Pi(Bpnt[i].pnt);
		Vector3D Ni(Bpnt[i].norm);

		for (j = 0; j < dmaskSize; j++)
		{
			delm = dmask -> getMaskElement(i,j);
			Vector3D coord = Pi + Ni * (dVals[i] + delm->dist) * r;
			targetImage->modelToImageCoordinates(coord);

			if (windowIntensities)
				imgVal = targetImage->mapDisplayToRelative(
					targetImage->getWindowedInterpolatedVoxelValue(
					coord.getX(), coord.getY(), coord.getZ()));
			else
				imgVal = targetImage->mapDisplayToRelative(
					targetImage->getInterpolatedVoxelValue(
					coord.getX(), coord.getY(), coord.getZ()));

			lastMatchValue += -delm->tempVal * (imgVal-mu)/rms * vol;
		}
	}

	lastMatchValue /= maskVolume;
//	cout << "\nmatch: " << lastMatchValue << '\n';
	return lastMatchValue;
}

double M3DVoxelOptimizer::computeLastPenaltyValue()
{
	int i, j;
	double delx, dely, delz, dv2;

	if (Bpnt == NULL || neighbors == NULL)
	{
		cout << "Error: no boundary information in BdOpt::compPenalty." << endl;
		return 0.0;
	}

	lastPenalty = 0.0;

	for (i = 0; i < numPts; i++)
	{
		lastPenalty += dVals[i]*dVals[i] / (2*sigma1*sigma1);

		for (j = 0; j < neighbors[i].get_degree(); j += 2)
		{
			int ij = neighbors[i][j];
			delx = Bpnt[i].pnt[0] - Bpnt[ij].pnt[0];
			dely = Bpnt[i].pnt[1] - Bpnt[ij].pnt[1];
			delz = Bpnt[i].pnt[2] - Bpnt[ij].pnt[2];
			dv2  = delx*delx + dely*dely + delz*delz;

			lastPenalty += (dVals[i]*Bpnt[i].rad - dVals[ij]*Bpnt[ij].rad) * 
				(dVals[i]*Bpnt[i].rad - dVals[ij]*Bpnt[ij].rad) / (2*sigma2*sigma2*dv2);
		}

		double area = dmask->getSurfaceElement(surfaceLevel, i)->area;
		lastPenalty *= area;
	}

	return lastPenalty;
}



