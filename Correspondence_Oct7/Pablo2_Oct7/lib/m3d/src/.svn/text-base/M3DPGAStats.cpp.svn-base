#include <math.h>
#include <fstream>

#include "matrix.h"
#include "M3DFigurePredictor.h"
#include "M3DPGA.h"
#include "M3DPGAStats.h"

#include <typeinfo>

/*
	Ja-Yeon 3/10/2005

	Currently, optimization of multi-objects using PG's based on shape space prediction 
	has not been implemented. 

*/

//#define DEBUG

using namespace std;

M3DPGAStats::M3DPGAStats()
{
	meanObj = NULL;
	intermedObj = NULL;
	latestObj = NULL;
	meanDiffObj = NULL;
	type = Scaled;	// Default is "Scaled"

	// CL: initialize best sim trans
	bestSimilarityTransform.setToIdentity();
}

M3DPGAStats::~M3DPGAStats()
{
	clear();
}

M3DPGAStats::M3DPGAStats(const M3DPGAStats & ps)
{
	meanObj = NULL;
	intermedObj = NULL;
	latestObj = NULL;
	meanDiffObj = NULL;
	type = ps.type;
	bestSimilarityTransform = ps.bestSimilarityTransform;
	PGs = ps.PGs;
	PGSets = ps.PGSets;
	Predictions = ps.Predictions;
}

void M3DPGAStats::clear()
{
	if (intermedObj != NULL)
		delete intermedObj;
	intermedObj = NULL;
	if (latestObj != NULL)
		delete latestObj;
	latestObj = NULL;
	if (meanDiffObj != NULL)
		delete meanDiffObj;
	meanDiffObj = NULL;

	// meanObj does not belong to this class
	meanObj = NULL;

	resetPGs();
}

void M3DPGAStats::resetPGs()
{
	int numPGSets = PGSets.size();
	for (int i = 0; i < numPGSets; i++) {
		if (! PGSets.empty()) {
			if (PGSets[i] != NULL)
				delete PGSets[i];
		}
		if (! PGs.empty()) {
			if (PGs[i] != NULL)
				delete PGs[i];
		}
		if (! Predictions.empty()) {
			if (i > 0 && i < numPGSets - 1 && Predictions[i - 1] != NULL)
				delete Predictions[i - 1];
		}
	}

	PGSets.clear();
	PGs.clear();
	Predictions.clear();
}

void M3DPGAStats::resetObjsToMean()
{
	if (intermedObj != NULL) {
		delete intermedObj;
		intermedObj = NULL;
	}
	if (latestObj != NULL) {
		delete latestObj;
		latestObj = NULL;
	}

	if (meanObj == NULL)
		return;

    intermedObj = meanObj->assign();
	latestObj = meanObj->assign();
}

PGData * M3DPGAStats::getPredictionPtr(int order)
{
	if (Predictions.empty())
		return NULL;
	else
		return Predictions[order];
}

void M3DPGAStats::setIntermedObj(M3DObject * obj)
{
	if (intermedObj != NULL)
		delete intermedObj;

	intermedObj = obj->assign();
}

void M3DPGAStats::setLatestObj(M3DObject * obj)
{
	if (latestObj != NULL)
		delete latestObj;

	latestObj = obj->assign();
}

int M3DPGAStats::findPGSetId(int figureId)
{
	// For single figure object, there will be always only one pgSet.
	if (meanObj->getFigureTreeCount() == 1) 
		return 0;

	// For ensemble object
	for (int i = 1; i < PGSets.size(); i++) {
		PGSet * pgSet = getPGSetPtr(i);

		for (int j = 0; j < pgSet->getNumFigs() ; j++) {
			if (pgSet->figIndexes[j] == figureId)
				return i;
		}
	}

	cout << "Warning: the pgSet of figure " << figureId << " is absent.\n";
	return -1;
}

void M3DPGAStats::removeLastPGSet()
{
	int numPGSets = getNumOfPGSets();

	// PGData to remove
	int last = --numPGSets;
	PGData *pgData = getPGDataPtr(last);
	int lenPG = pgData->lenPG;
	int i,j;

	
	if (numPGSets == 1) //single figure - remove PG Stats
		clear();
	else {
		//first set for ensemble object
		PGSet  * pgSet  = getPGSetPtr(0);
		PGData * pgData = getPGDataPtr(0);
		
		pgSet->figIndexes.pop_back();

		int keep = pgData->lenPG - lenPG;
		pgData->lenMean = keep;
		pgData->lenPG = keep;
		for(i=0; i<pgData->numPGs; i++){
			double *v = new double[keep];
			for (j=0; j<keep; j++){
				v[j] = (pgData->pgVec[i])[j];
			}
			delete pgData->pgVec[i];
			pgData->pgVec[i] = v;
		}
		
		if (pgData->meanRes != NULL){
			double *v = new double[keep];
			for (j=0; j<keep; j++){
				v[j] = pgData->meanRes[j];
			}
			delete pgData->meanRes;
			pgData->meanRes = v;
		}
		
		//Remove the last pg set and data
		pgSet = getPGSetPtr(last);
		if (pgSet != NULL)
			delete pgSet;
		pgData = getPGDataPtr(last);
		if (pgData != NULL)
			delete pgData;
		PGSets.pop_back();
		PGs.pop_back();
	}
}

int M3DPGAStats::totalNumParams(int depth)
{
	int totalParams = 0;
	for (int i = 0; i <= depth; i++) {
		PGData * pg = getPGDataPtr(i);
		if (pg == NULL)
			return -1;
		totalParams += pg->numPGs;
	}
	return totalParams;
}


/*	This function is to apply mean of residue to 'latestObj' (which reflects
	all the  deformations made up to 'order-1') and stored back to 'obj'.
	Normally, 'obj' is passed in from P3DControl class.
*/
bool M3DPGAStats::applyMeanResidue(M3DObject * obj, int order)
{	
	return applyMeanResidue(obj, latestObj, order);
}


/*	This function is to deform 'intermedObj (which reflects all the deformations
	made up to 'order-1' and mean of residue) and stored back in 'obj'.
	Normally, 'obj' is passed in from P3DControl class.
*/
bool M3DPGAStats::doPGADeform(M3DObject * obj, std::vector<double> vals,
							  int order, bool predict)
{
	int numParams,
		i;

	PGData *pgData = getPGDataPtr(order);
	if (pgData == NULL) {
		cout << "Null PG in M3DPGAStats::doPGADeform(..)\n";
		return false;
	}

	// Compare the size of vals and number of pgs 
	numParams = pgData->numPGs;
	if (vals.size() < numParams)
		for (i = vals.size(); i < numParams; i++)
			vals.push_back(0.0);

	if (obj == NULL)
		return false;

	return pgaDeform(obj, intermedObj, vals, order, predict);
}


/*------------------------- PRIVATE MEMBER FUNCTIONS ------------------------------------*/

// Called ONLY ONCE when the PGA statistics file is read.
// Convert meanDiff (PG[0].meanRes) into M3DObject. 
void M3DPGAStats::convertMeanDiff()
{
	int numSets = getNumOfPGSets();
	if (numSets < 0) {
		cout<<"Error: No geometric statistics "<< endl;
		return;
	}

	if (meanObj != NULL) {
		meanDiffObj = meanObj->assign();
	}
	else {
		cout<< "NULL mean object when M3DPGAStats::convertMeanDiff() is called" << endl;
		return;
	}

	int nFigs = meanDiffObj->getFigureCount();
	int nPrims = 0;
	M3DFigure *figPtr = NULL;

	if (nFigs != PGSets[0]->getNumFigs()) {
		cout << "Error: Not all figures are included in global statistics" << endl;
		return;
	}

	PGData *globalPG = getPGDataPtr(0);	// global statistics
	// FIXME: will break tubes
	if ((globalPG->lenMean)/NUM_PGA_PARAMS != meanDiffObj->getPrimitiveCount()) {
		cout << "Error: Number of primitives in statistics and model are not same" << endl;
		return;
	}

	double *meanPtr = globalPG->meanRes;

	int count = 0;

	for (int i=0; i<nFigs; i++) {

		figPtr = meanDiffObj->getFigurePtr(i);
		nPrims = figPtr->getPrimitiveCount();
		// FIXME: will break tubes
		count = i*NUM_PGA_PARAMS;

		M3DPrimitive *prim =NULL;
		Vector3D n0, n1;
		Quat q; 
		double theta;

		for (int j=0; j<nPrims; j++) {
			// FIXME: will break tubes
			int offset = count + j*NUM_PGA_PARAMS;

			prim = figPtr->getPrimitivePtr(j);

			Vector3D pos(meanPtr[offset], meanPtr[offset+1], meanPtr[offset+2]);
			prim->setX(pos);

			prim->setR(exp(meanPtr[offset+3]));


			n0 = ShapeSpace::S2::Exp(Vector2D(meanPtr[offset+4], meanPtr[offset+5]));
			n1 = ShapeSpace::S2::Exp(Vector2D(meanPtr[offset+6], meanPtr[offset+7]));

			symToLieAtom(q, theta, n1, n0);

			prim->setQ(q);
			prim->setTheta(theta);
		}		
	}
}

// For adaptive case, convertMeanDiff() should be called 
// before rescale() is called.

void M3DPGAStats::rescale()
{
	M3DObject * mean;


	int numSets = getNumOfPGSets();
	if (numSets < 0) {
		cout << "Error: No geometric statistics" << endl;
		return;
	}

	if (type) {	// scaled
		cout << "The PG statistics are already mean(R) scaled" << endl;
		return; 
	}

	if (meanObj->isAdaptive())
		mean = meanDiffObj;	// At this point, meanDiffObj != NULL  
	else
		mean = meanObj;

	if (mean == NULL) {
		cout << "No mean model for PGA statistics rescaling" << endl;
		return; 
	}

	PGSet *pgSet = NULL;
	PGData *pgData = NULL;

	int nFigs, nPrims, nAugs, figIdx, primIdx;
	M3DFigure *figPtr = NULL;
	M3DPrimitive *primPtr = NULL;

	int i, j, k, l;

	double meanR = 1.0;

	for (i=0; i<numSets; i++)
	{
		pgSet = getPGSetPtr(i);
		pgData = getPGDataPtr(i);

		int offset = 0;

		nFigs = pgSet->getNumFigs();
		for (j=0; j<nFigs ; j++) {

			figIdx = pgSet->figIndexes[j];
			figPtr = mean->getFigurePtr(figIdx);
			if (figPtr == NULL)
				continue;

			nPrims = figPtr->getPrimitiveCount();			

			for (k=0; k<nPrims; k++) {

				primPtr = figPtr->getPrimitivePtr(k);
				meanR = primPtr->getR();

				int nPGs = pgData->numPGs;
				// FIXME: will break tubes
				for (l=0; l<nPGs; l++) {
					pgData->pgVec[l][offset+k*NUM_PGA_PARAMS+3] /= meanR;  // radius
					pgData->pgVec[l][offset+k*NUM_PGA_PARAMS+4] /= meanR;  // n0
					pgData->pgVec[l][offset+k*NUM_PGA_PARAMS+5] /= meanR;
					pgData->pgVec[l][offset+k*NUM_PGA_PARAMS+6] /= meanR;  // n1
					pgData->pgVec[l][offset+k*NUM_PGA_PARAMS+7] /= meanR;
					pgData->pgVec[l][offset+k*NUM_PGA_PARAMS+8] /= meanR;  // elongation
				}	
			}
			// FIXME: will break tubes
			offset += nPrims*NUM_PGA_PARAMS;
		}

		nAugs = pgSet->getNumAugs();

		for (j=0; j<nAugs; j++)
		{
			figIdx = (pgSet->augmentations[j])->figIndex;
			nPrims = (pgSet->augmentations[j])->getNumPrims();

			for (k=0; k<nPrims; k++) {

				primIdx = (pgSet->augmentations[j])->primIndexes[k];
				primPtr = mean->getPrimitivePtr(figIdx, primIdx);
				if (primPtr == NULL)
					continue;

				meanR = primPtr->getR();

				int nPGs = pgData->numPGs;
				for (l=0; l<nPGs; l++) {
					pgData->pgVec[l][offset+3] /= meanR;  // radius
					pgData->pgVec[l][offset+4] /= meanR;  // n0
					pgData->pgVec[l][offset+5] /= meanR;
					pgData->pgVec[l][offset+6] /= meanR;  // n1
					pgData->pgVec[l][offset+7] /= meanR;
					pgData->pgVec[l][offset+8] /= meanR;  // elongation
				}
				// FIXME: will break tubes
				offset += NUM_PGA_PARAMS;
			}
		}
	}
}


//PRIVATE MEMBER FUNCTIONS -----------------------------------------------------------------------------------

/*	This function is to deform 'startObj' by principal geodesics at 'order' level.
	Normally, 'startObj' will be 'latestObj' in class when this function is called and 
	will be called ONLY ONCE before doing any pga deformation.
*/
bool M3DPGAStats::applyMeanResidue(M3DObject * obj, M3DObject *startObj, int order)
{
	int size,
		lenMean, i;

	PGData *pgData = getPGDataPtr(order);

    if (obj == NULL || startObj == NULL || pgData == NULL)
        return false;

	size = pgData->lenPG;
	lenMean = pgData->lenMean;

#ifdef DEBUG
	cout << "In M3DPGAStats::applyMeanResidue(),  " << "lenPG = " << size << endl;
#endif

	// Note :length of mean residue = length of principal geodesic component
	double *meanRes = new double[size];
	for (i = 0; i < size; i++) {
		if (lenMean == 0)
			meanRes[i] = 0.0;
		else if (lenMean == size)
			meanRes[i] = pgData->meanRes[i];
		else {	// This should not happen!!
			cout << "Error: check PGs\n";
			delete [] meanRes;
			return false;
		}
	}

	if (! pgaDeformToSet(obj, startObj, meanRes, order, NULL, false))
	{
		delete [] meanRes;
		cout << "Error: check pgaDeformToSet()\n";
		return false;
	}

	// Set intermediate object : mean + all changes up to
	// 'order-1' level + mean residue at 'order' level
	setIntermedObj(obj);
	setLatestObj(obj);

    delete [] meanRes;
    return true;
}



/*	This function is to deform 'startObj'.
	'vals' contains coefficients for principal geodesics in pgSet[order].
	Normally, startObj will be intermedObj.
*/
bool M3DPGAStats::pgaDeform(M3DObject * obj, M3DObject *startObj,
	std::vector<double> vals, int order, bool predict)
{
	int numParams,
		size,
		i, j;

	if (order < 0 || order >= getNumOfPGSets())
		return false;

#ifdef DEBUG
	cout << "order = " << order << '\n';
#endif

	PGData *pgData = getPGDataPtr(order);
	if (obj == NULL || startObj == NULL  || pgData == NULL)
        return false;

	// Compare the number of principal geodesic components and the size of vals
    if (pgData->numPGs > vals.size())
        numParams = vals.size();
    else
        numParams = pgData->numPGs;

	size = pgData->lenPG;


	double * v1 = new double[size];
    for (i = 0; i < size; i++)
        v1[i] = 0.0;

    for (i = 0; i < numParams; i++)
        for (j = 0; j < size; j++)
            v1[j] += vals[i] * (pgData->pgVec[i])[j];

	if (! pgaDeformToSet(obj, startObj, v1, order, &vals, predict))
	{	
		delete [] v1;
		cout << "Error: check pgaDeformToSet()\n";
		return false;
	}

	// Set intermediate object : mean + all changes up to
	// order-1 + mean residue of order 
	setLatestObj(obj);

    delete [] v1;

	return true;
}

/*
	This function is to apply 'v1' to every atom indicated in 'PGSets[order]' by
	calling pgaDeformToPrim() for every atom.
	'size' : size of 'v1'
	'order' : which set
*/
bool M3DPGAStats::pgaDeformToSet(M3DObject * obj, M3DObject * startObj, const double * v1, 
								 int order, const std::vector<double> * vals,  bool predict)
{
//	int countPrims; 
	int	numFigs, figId, i, j, lenVec;

    M3DPrimitive * prim,
                 * meanPrim; 

	PGSet *pgSet = getPGSetPtr(order);
	PGData *pgData = getPGDataPtr(order);

    if (obj == NULL || startObj == NULL || v1 == NULL || pgData == NULL)
        return false;

	numFigs = pgSet->getNumFigs();

	// 1) Check if the number of atoms matches the length of principal geodesics

	// Length of principal geodesic component = sum(over all primitives) {prim.NUM_PGA_PARAMS}
	// + sum(over all augmented primitives) {prim.NUM_PGA_PARAMS}
	static bool checked = false;
	if (!checked) {
		lenVec	= 0;
		for(i = 0; i != numFigs; ++i ) {
			figId = pgSet->figIndexes[i];
			lenVec	+= (obj->getFigurePtr(figId)->getPrimitiveCount())*(obj->getFigurePtr(figId)->NUM_PGA_PARAMS());
		}
		
		for(i = 0; i != pgSet->getNumAugs(); i++) {
			AugmentedAtoms * aug = pgSet->augmentations[i];
			if (aug == NULL) {
				return false;
			}
			
			figId = aug->figIndex;
			M3DFigure *fig = obj->getFigurePtr(figId);
			
			// FIXME: will break tubes 
			// JJ: This should fix breaking tube
			if (typeid(*fig) == typeid(M3DQuadFigure)) {
				lenVec	+= aug->getNumPrims() * fig->NUM_PGA_PARAMS();
			}
			else if (typeid(*fig) == typeid(M3DTubeFigure)){
				lenVec	+= aug->getNumPrims() * fig->NUM_PGA_PARAMS();
			}		
		}
		
		if (lenVec != pgData->lenPG) {
			std::cerr << "pgaDeformToSet: length mismatch " << __FILE__ << ":" << __LINE__ << std::endl;
			std::cerr << "numFigs = " << numFigs << ",order = " << order << ",lenVec = " << lenVec << ", pgdata-length = " << pgData->lenPG << std::endl;
			return false;
		}
		checked = true;
	}

	// 2) Deform each atom in whole figures in pgSet[order]

	//countPrims	= 0;
	
	int offset	= 0;
	for (i = 0; i < numFigs ; i++)
	{	
		figId = pgSet->figIndexes[i];
		M3DFigure *fig = obj->getFigurePtr(figId);

		for (j = 0; j < fig->getPrimitiveCount() ; j++)
		{
			prim = obj->getPrimitivePtr(figId, j);
			meanPrim = startObj->getPrimitivePtr(figId, j);
			//int offset = (countPrims+j)*NUM_PGA_PARAMS;

			pgaDeformToPrim(prim, meanPrim, v1+offset);
			offset	+= fig->NUM_PGA_PARAMS();
		}

		//
		// Correct any consistency issues.
		//
		fig->fixGlobalConsistency();

		//countPrims += fig->getPrimitiveCount();
	}

	// 3) Deform atoms in augmented figures in pgSet[order]
	for (i = 0; i < pgSet->getNumAugs() ; i++)
	{
		AugmentedAtoms *aug = pgSet->augmentations[i];		
		figId = aug->figIndex;
		
		for (j = 0; j < aug->getNumPrims() ; j++)
		{
			int primId = aug->primIndexes[j];
			prim = obj->getPrimitivePtr(figId, primId);
			meanPrim = startObj->getPrimitivePtr(figId, primId);
			pgaDeformToPrim(prim, meanPrim, v1+offset);
			offset	+= obj->getFigurePtr(figId)->NUM_PGA_PARAMS();			
		}	
	}  
	
	return true;
}

bool M3DPGAStats::pgaDeformToPrim(M3DPrimitive *prim, M3DPrimitive *meanPrim, const double *v )
{
    Vector3D axis;
    double theta;

	if (prim == NULL || meanPrim == NULL || v == NULL)
		return false;
	// TODO: It would be nice if this code would go into some common place.
	if( typeid(*prim) == typeid(M3DQuadPrimitive) ||typeid(*prim) ==typeid(M3DQuadEndPrimitive))  {
		const Vector3D pos(v[0], v[1], v[2]);

		prim->setX(meanPrim->getX() + pos);
		prim->setR(meanPrim->getR() * exp(v[3]));

		Vector3D n0 = meanPrim->getNormalizedY0();
		Vector3D n1 = meanPrim->getNormalizedY1();
		Quat q0 = ShapeSpace::S2::rotationFromOrigin(n0);
		Quat q1 = ShapeSpace::S2::rotationFromOrigin(n1);

		n0 = ShapeSpace::S2::Exp(Vector2D(v[6], v[7]));
		n1 = ShapeSpace::S2::Exp(Vector2D(v[4], v[5]));
		q0.rotateVector(n0);
		q1.rotateVector(n1);

		Quat q;
		symToLieAtom(q, theta, n1, n0);

		prim->setQ(q);
		prim->setTheta(theta);

		if (prim->type() == M3D_END_PRIMITIVE)
			(dynamic_cast<M3DQuadEndPrimitive *>(prim))->setElongation((dynamic_cast<M3DQuadEndPrimitive *>(meanPrim))
			->getElongation() * exp(v[8]));
	}
	else if( typeid(*prim) == typeid(M3DTubePrimitive) ||typeid(*prim) ==typeid(M3DTubeEndPrimitive))  {
		const Vector3D pos(v[0], v[1], v[2]);

		// Now, rotate the spokes
		Vector3D t	= ShapeSpace::S2::Exp(Vector2D(v[4], v[5]));
		// This should by definition be the same as t.
		//Vector3D tnew(1.0, 0.0, 0.0);
		//ShapeSpace::S2::rotationFromOrigin(t).rotateVector(tnew);
		Vector3D n(0.0, 1.0, 0.0);
		ShapeSpace::S2::rotationFromOrigin(t).rotateVector(n);
		Vector3D b	= t;
		meanPrim->getQ().rotateVector(t);
		meanPrim->getQ().rotateVector(n);

		theta	= meanPrim->getTheta();
		// FIXME: verify RP1 exp and log maps.
		theta	= ShapeSpace::RP1::Exp( ShapeSpace::RP1::Log(theta) + v[6] );

		// First add the translation
		prim->setX(meanPrim->getX() + pos);
		// Next multiply the radius
		prim->setR(meanPrim->getR() * exp(v[3]));
		Quat q;
		q.buildFromFrame(t,n);
		// Successive accumulation of error leads to q,t,n,b becoming non-unit.
		q.normalize();
		prim->setQ( q );
		prim->setTheta(theta);
		// 
		if (prim->type() == M3D_END_PRIMITIVE)
			(dynamic_cast<M3DTubeEndPrimitive *>(prim))->setElongation((dynamic_cast<M3DTubeEndPrimitive *>(meanPrim))
				->getElongation() * exp(v[7]));
	}
	return true;
}

int M3DPGAStats::primitiveCount(int order)
{
	if (order > PGSets.size()) {
		cout << "Error: PGSet ID is not within range\n";
		return -1;
	}

	PGSet * pgSet = PGSets[order];

	int i;
	int numPrims;

	numPrims = 0;
	if (meanObj != NULL) {
		for (i = 0; i < pgSet->getNumFigs(); i++) {
			int figId = pgSet->figIndexes[i];
			M3DFigure * fig = meanObj->getFigurePtr(figId);

			if (fig != NULL)
				numPrims += fig->getPrimitiveCount(); 
		}
	}

	for (i = 0; i < pgSet->getNumAugs(); i++) {
		AugmentedAtoms * augs = pgSet->augmentations[i];
		numPrims += augs->getNumPrims();
	}

	return numPrims;
}

// FIXME: will break tubes

// Find coeffient of pgs in the remainder set and store it in 'remVals'
void M3DPGAStats::prediction(std::vector<double> &remVals, M3DObject * obj,
	std::vector<double> vals, int order)
{
	if (order<1 || order>getNumOfPGSets() - 1)
		return;

	PGSet *pgSet = getPGSetPtr(order);
	PGData *pgData = getPGDataPtr(order);
	PGData *predData = getPredictionPtr(order-1);

	int nPGs, i;

	// Compare the number of principal geodesic components and the size of vals
    if (pgData->numPGs > vals.size())
        nPGs = vals.size();
    else
        nPGs = pgData->numPGs;

	int nAugs = pgSet->getNumAugs();
	if (nAugs < 1)	// Do nothing
		return;	

	// XXX For now, assuming that, in object level, each set has only one figure
	// extraction of part of PG's and mean's in target and remainder sets
	int	figId = pgSet->figIndexes[0];  

	M3DFigure *fig = obj->getFigurePtr(figId);
	int offset = fig->getPrimitiveCount()*NUM_PGA_PARAMS;

	int nAugAtoms = 0;
	for (i = 0; i < nAugs ; i++) {
		AugmentedAtoms *aug = pgSet->augmentations[i];
		if (aug == NULL)
			return;
		nAugAtoms += aug->getNumPrims();
	}

	Matrix bMean(nAugAtoms*NUM_PGA_PARAMS, 1);
	Matrix bPG(nAugAtoms*NUM_PGA_PARAMS, nPGs);	// PG's in target set corresponding 
												// to augmented atoms

	for (int c=0; c<nPGs ; c++) {
		for (int r=0; r<nAugAtoms*NUM_PGA_PARAMS ; r++) {
			bPG(r, c) = (pgData->pgVec[c])[offset+r];	// at the end of pg vectors
			if (! c)
				bMean(r, c) = (pgData->meanRes)[offset+r];
		}
	}

	Matrix rMean(nAugAtoms*NUM_PGA_PARAMS, 1);
	Matrix rPG(nAugAtoms*NUM_PGA_PARAMS, nPGs);	// PG's in remainder set corresponding 
												// to augmented atoms
	int nAtoms = 0;
	for (i = 0; i < nAugs ; i++) {
		AugmentedAtoms *aug = pgSet->augmentations[i];
		if (aug == NULL)
			return;

		for (int j=0; j<aug->getNumPrims(); j++) {
			int primId = aug->primIndexes[j];

			for (int c = 0; c < nPGs; c++)
				for (int r=0 ; r<NUM_PGA_PARAMS ; r++) {
					rPG(j*NUM_PGA_PARAMS+r, c) = (predData->pgVec[c])[(nAtoms+primId)*NUM_PGA_PARAMS + r];
					if (! c)
						rMean(j*NUM_PGA_PARAMS+r, c) = (predData->meanRes)[(nAtoms+primId)*NUM_PGA_PARAMS + r];
				}			
		}

		M3DFigure *fig = obj->getFigurePtr(aug->figIndex);
		nAtoms += fig->getPrimitiveCount();	
	}

	Matrix beta(nPGs, 1);
	for (i=0 ; i<nPGs ; i++)   
		beta(i, 0) = vals[i];

	Matrix b, A, invA, alpha;
	b = rPG.t()*(bMean - rMean + bPG*beta);
	A = rPG.t()*rPG;
	A.inverse(invA);
	alpha = invA*b;

	for (i=0 ; i<nPGs ; i++)
		remVals.push_back(alpha(i, 0));
}

	// Deform atoms in augmented figures in pgSet[order]
	/* 
		Ja-Yeon 3/10/05
		Currently, for loop below works only when numAugs == 1 
		due to prediction.
	*/
	/*numAugs = pgSet->getNumAugs();
	offset	= 0;
	for (i = 0; i < numAugs ; i++)
	{
		AugmentedAtoms *aug = pgSet->augmentations[i];
		if (aug == NULL) {
			return false;
		}

		figId = aug->figIndex;
		for (j = 0; j < aug->getNumPrims() ; j++)
		{
			int primId = aug->primIndexes[j];
			prim = obj->getPrimitivePtr(figId, primId);
			meanPrim = startObj->getPrimitivePtr(figId, primId);
			// FIXME: augmentation?
			//int offset = (countPrims+j)*NUM_PGA_PARAMS;

			pgaDeformToPrim(prim, meanPrim, v1+offset, obj->getFigurePtr(figId));
			offset	+= obj->getFigurePtr(figId)->NUM_PGA_PARAMS();			
		}	
		//countPrims += aug->getNumPrims();

	/*
		Ja-Yeon 3/10/05
		Currently, prediction based on approximate similarity transform 
		works only when the augmented atoms are from one figure.
		Take that part out of Qiong's code???
	*/
	/*	if (predict)
		{
			//===========================================================
			// Applying prediction based on shape space
			//===========================================================

			//Note: 0< order< getNumOfPGSets()-1 here
			//because numAugs = 0 at ensemble stage of order 0.

			bool doSimPrediction = false;

			if (Predictions.empty())
				doSimPrediction = true;

			if (! doSimPrediction) 
			{
				PGData *predData = getPredictionPtr(order-1);
				if (predData == NULL) {
					cout << "Error: M3DPGAStats::doPGADeformToSet()" << endl;
					return false;
				}

				int size = predData->lenPG;
				int numParams = predData->numPGs;

				std::vector<double> remVals;
				prediction(remVals, obj, (*vals), order);

				double * v2 = new double[size];  // Components for PGs in remainder 
				for (i = 0; i < size; i++)
					v2[i] = 0.0;

				for (i = 0; i < numParams; i++)
					for (j = 0; j < size; j++)
						v2[j] += (*vals)[i] * (predData->pgVec[i])[j];

					
					// Ja-Yeon 3/10/05
					//	Need to rewrite to cover general cases.
					//	Check prediction() also!!!
					
					int objId = pgSet->figIndexes[0];  
					int countPrimsInRemainder = 0;
					int nFigs = obj->getFigureCount();

					// Bladder 0 => then prostate 1, rectum 2		AGG:	THIS SHOULD BE ORGAN-INDEPENDENT
					for (int k = objId + 1; k < nFigs; k++) {
						M3DFigure *fig = obj->getFigurePtr(k);	
						for (int l = 0; l < fig->getPrimitiveCount() ; l++) {
							prim = obj->getPrimitivePtr(k, l);
							meanPrim = startObj->getPrimitivePtr(k, l);
							// FIXME: I've no clue about how to handle this - rrs
							// FIXME: will break tubes
							int offset = (countPrimsInRemainder+l)*NUM_PGA_PARAMS;
							pgaDeformToPrim(prim, meanPrim, v2+offset, fig);
						}
						countPrimsInRemainder += fig->getPrimitiveCount();
					}

					delete [] v2;					
			}

			//==============================================================
			// Applying prediction based on approximate similarity transform
			//==============================================================

			else if (doSimPrediction)
			{
				//Predict the remaining atoms of the figure of figId
				M3DQuadFigure *newFigPtr = dynamic_cast<M3DQuadFigure*>(obj->getFigurePtr(figId));
				M3DQuadFigure *figPtr = dynamic_cast<M3DQuadFigure*>((startObj->getFigurePtr(figId))->assign());
				int linkCount = aug->getNumPrims();
				int *linkIds = new int[linkCount];
				for (i = 0; i < linkCount; i++)
					linkIds[i] = aug->primIndexes[i];

				M3DFigurePredictor figPredictor(linkCount, linkIds, figPtr, newFigPtr);
				obj->replaceFigure(figId, figPtr);

				delete [] linkIds;

				//Ja-Yeon 11/2/04 : propagate similarity transformation to other figures
				const SimilarityTransform * simTrans = figPredictor.bestTransform();
				int numPGSets = getNumOfPGSets();
				for (int k = order; k < numPGSets; k++)
				{
					PGSet *chkSet = PGSets[k];
					for (int l = 0; l < chkSet->getNumFigs(); l++)
					{
						if (chkSet->figIndexes[l] == figId)
						{
							numAugs = chkSet->getNumAugs();
							for (int m = 0; m < numAugs ; m++)
							{	
								AugmentedAtoms *aug = chkSet->augmentations[m];
								if (aug == NULL)  return false;
								int otherFigId = aug->figIndex;

								int linkCount = aug->getNumPrims();
								int *linkIds = new int[linkCount];
								for (int n = 0; n < linkCount; n++)
								{
									linkIds[n] = aug->primIndexes[n];
									prim =  obj->getPrimitivePtr(otherFigId, linkIds[n]);

									//apply similarity transform
									Vector3D diff = prim->getX() - simTrans->COR;
									Vector temp(3, diff.getX(), diff.getY(), diff.getZ());
									temp = simTrans->scale*((Matrix)(simTrans->rotM)*temp);
									Vector3D temp2(temp(0), temp(1), temp(2));
									prim->setX(temp2 + simTrans->COR + simTrans->trans);
									prim->setR(prim->getR()*simTrans->scale);

									prim->rotateBy(simTrans->rotQ);

									if (prim->type()==M3D_END_PRIMITIVE)
										(dynamic_cast<M3DEndPrimitive *>(prim))->setElongation(
											(dynamic_cast<M3DEndPrimitive *>(prim))
											->getElongation()*simTrans->scale);
								}	

								//Predict the remaining atoms of the figure of figId
								M3DQuadFigure *newFigPtr = dynamic_cast<M3DQuadFigure*>(obj->getFigurePtr(otherFigId));
								M3DQuadFigure *figPtr = dynamic_cast<M3DQuadFigure *>(startObj->getFigurePtr(otherFigId)->assign());

								M3DFigurePredictor figPredictor(linkCount, linkIds, figPtr, newFigPtr);
								obj->replaceFigure(otherFigId, figPtr);

								delete [] linkIds;
							}							
						}
					}
				}
			} // END if (doSimPrediction)

		}  // END if (predict)

	}  //END for (i = 0; i < numAugs ; i++)
	*/

