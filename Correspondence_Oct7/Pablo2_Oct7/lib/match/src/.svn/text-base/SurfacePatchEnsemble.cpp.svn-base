#include <iostream>
#include <fstream>
#include "SurfacePatchEnsemble.h"
#include "P3DControl.h"
#include "utility.h"

#include <algorithm>

#include <typeinfo>

using namespace std;

static bool printTime = false;

bool lessProfileComponent(ProfileComponent c1, ProfileComponent c2)
{
	return c1.weight < c2.weight;
}

bool lessQuantileComponent(QuantileComponent c1, QuantileComponent c2)
{
	return c1.value < c2.value;
}

/************************************************************************************
  SurfacePatchEnsemble
************************************************************************************/

SurfacePatchEnsemble::SurfacePatchEnsemble(HistogramMatchType hmtIn,
										   const char *filename,
										   const char * QFDataFilename)
{
	hmt = hmtIn;

			
	// Load the registry
	try {
		reg.readFromFile(filename);
	}
	catch (RException excp) {
		excp.print(std::cout);
		cout << "ERROR: SurfacePatchEnsemble file, " << filename 
             << ", not found !!OR!! ERROR loading.\n";
	}

	int version = reg.getIntValue("config.version", -1);
	if (version != 7)
		printf("WARNING: Incorrect surface patch ensemble version loaded. Current version: 7\n");

	//Do everything other than QFData during initialize -- when have object

	numFigures = 0;
	numPatches = 0;
	pList = NULL;

	steps = NULL;
	checkImage = NULL;
	mean = NULL;
	scaleVector = NULL;

	numFeatures = 1;
	//featureVectors = NULL;
	//featureVectorLength = 0;
	
	numValues = 0;
	numVectors = 0;
	vectors = NULL;
	vectorLengths = NULL;
	vectorOffsets = NULL;
	vectorVariances = NULL;
	vectorResidualID = NULL;
	vectorAtomID = NULL;
	vectorFigureID = NULL;

	numResiduals = 0;
	residualVariances = NULL;
	residualMasks = NULL;
	residualMaskLengths = NULL;

	initialized = false;

	//Load QFData if exists
	QFData = NULL;
	QFDataDim[0] = 0; QFDataDim[1] = 0; QFDataDim[2] = 0;
	QFDataMargin[0] = 0; QFDataMargin[1] = 0; QFDataMargin[2] = 0;

	if(QFDataFilename != NULL){

		FILE *fid = fopen(QFDataFilename, "rb");
		if(fid == NULL){
			printf("ERROR: QFDataFilename could not be opened for reading: %s\n", QFDataFilename);
			return;
		}

		fread(&numFeatures, sizeof(int), 1, fid);
		fread(QFDataDim,    sizeof(int), 3, fid);
		fread(QFDataMargin, sizeof(int), 3, fid);

		printf("QFData: numFeatures = %d, dimensions(%d, %d, %d), margin(%d, %d, %d)\n", numFeatures, 
               QFDataDim[0], QFDataDim[1], QFDataDim[2], QFDataMargin[0], QFDataMargin[1], QFDataMargin[2]);

		int numData = QFDataDim[0] * QFDataDim[1] * QFDataDim[2] * numFeatures;

		float *QFDataTemp = new float[numData];		
		QFData           = new double[numData];		

		//fread(QFData, sizeof(double), numData, fid);
		int numRead = safeFread(QFDataTemp, sizeof(float) * numData, fid);

		if(numRead != sizeof(float) * numData){
			printf("ERROR reading all of QFData\n");
			fclose(fid);
			exit(1);
		}
		fclose(fid);

		for(int i = 0; i < numData; i++){
			QFData[i] = (double) QFDataTemp[i];
		}
		delete [] QFDataTemp;
		QFDataTemp = NULL;
	}
}


void SurfacePatchEnsemble::setMean(double *newMean)
{
	for (int n=0; n < numValues; n++)
		mean[n] = newMean[n];
}

void SurfacePatchEnsemble::initialize(const M3DObject * object, Image3D * refImage){


	int n, i, j, k;

	if(initialized){
		printf("ERROR: SPE already initialized!");
		return;
	}

	initialized = true;

	// Sampling Variables
	histogramPatchType = intToHistogramPatchType( 
        reg.getIntValue("config.histogramPatchType", 0));
	inDistance = reg.getDoubleValue("config.inDistance", 1.0);
	outDistance = reg.getDoubleValue("config.outDistance", 1.0);
	inGaussianSTD = reg.getDoubleValue("config.inGaussianSTD", 0.5);
	outGaussianSTD = reg.getDoubleValue("config.outGaussianSTD", 0.5);

	// Variables that determine the regions
	numInProfileSteps        = reg.getIntValue("config.numInProfileSteps", 4);
	numOutProfileSteps       = reg.getIntValue("config.numOutProfileSteps", 6);
	voxelProfileSurfaceLevel = reg.getIntValue("config.voxelProfileSurfaceLevel", 4);
	usePixelNeighbors        = reg.getIntValue("config.usePixelNeighbors", 1);

	// Local region definition
	localPatchSurfaceLevel	= reg.getIntValue("config.localPatchSurfaceLevel", 0);
	localPatchRadius				= reg.getDoubleValue("config.localPatchRadius", 1.5);
	localPatchSTD						= reg.getDoubleValue("config.localPatchSTD", 0.75);

	// Gas and Bone defaults
	useLowHigh							= reg.getIntValue("config.useLowHigh", 1);
	lowThreshold						= reg.getIntValue("config.lowThreshold",   800);
	highThreshold						= reg.getIntValue("config.highThreshold", 1200);

	// Feature collection
	numLowBins							= reg.getIntValue("config.numLowBins",   1);
	numMidBins							= reg.getIntValue("config.numMidBins", 100);
	numHighBins							= reg.getIntValue("config.numHighBins",  1);

	//neighborhoodSize        = reg.getIntValue("config.neighborhoodSize",  1);
	int configNumFeatures             = reg.getIntValue("config.numFeatures",  0);
	if(configNumFeatures == 0)
		reg.setIntValue("config.numFeatures", numFeatures);
	else if(configNumFeatures != numFeatures){
		printf("ERROR: mismatch between QFData and config NUMFEATURES\n");
	}

	/*
	if(numFeatures == 0 || (numFeatures == 1 && neighborhoodSize == 1)){
		numFeatures = 1;
		reg.setIntValue("config.numFeatures", 1);
		if(neighborhoodSize != 1){ printf("ERROR: neighborhoodSize != 1 but no features are supplied. Why?\n"); }
		neighborhoodSize = 1;
		featureVectors = new double *[1];
		featureVectorLength = 1;
		featureVectors[0] = new double[1];
		featureVectors[0][0] = 1;
	}
	else{

		featureVectors = new double *[numFeatures];
		featureVectorLength = neighborhoodSize * neighborhoodSize * 1;

		for(i = 0; i < numFeatures; i++){
			//get partial vector and save length
			featureVectors[i] = reg.getDoubleArray("config.featureVectors[%d]", &n, i);
			if(featureVectorLength != n){ printf("ERROR: featureVectors and neighborhoodSize mismatch\n"); }
		}
	}
	*/

	// Don't need PCA settings

	///////////////////////////////////////////
	// Initialize
	///////////////////////////////////////////
	numFigures = object->getFigureCount();

	int f;

	if(histogramPatchType == GLOBAL_HPT) {
		numPatches = numFigures;
		for(f = 0; f < numFigures; f++) {

			//add ALL profiles to the patch with a weight of 1
			int numProfiles;
			Bpoint *Bpoints;
			ProfileComponent pc;
			pc.weight = 1;

			getBPoints(object, f, voxelProfileSurfaceLevel, numProfiles, Bpoints);
			int offset = f * 2 * (3 + numLowBins + numMidBins + numHighBins) * numFeatures;
			patches.push_back(new SurfacePatch(f, -1, offset, numFeatures));
			for(int p = 0; p < numProfiles; p++){
				pc.index = p;
				patches[f]->profiles.push_back(pc);
			}
		}
	} else {
		initializeLocalRegions(object, refImage);
	}

	numValues = numPatches * 2 * (3 + numLowBins + numMidBins + numHighBins) * numFeatures;
	values = new double [numValues];
	for(i = 0; i < numValues; i++){
		values[i] = 0;
	}
	
	///////////////////////////////////////////
	// Load Statistics if exist
	///////////////////////////////////////////

	if(reg.hasKey("stats")){

		mean = reg.getDoubleArray("stats.mean", &n);
		if(n != numValues){ printf("ERROR: mean incorrect size\n"); }

		scaleVector = reg.getDoubleArray("stats.scaleVector", &n);
		if(n != numValues){ printf("ERROR: scaleVector incorrect size\n"); }

		numVectors = reg.getIntValue("stats.numVectors", 0);

		vectorAtomID = reg.getIntArray("stats.vectorAtomID", &n);
		if(n != numVectors){ printf("ERROR: vectorAtomID incorrect size\n");}

		vectorFigureID = reg.getIntArray("stats.vectorFigureID", &n);
		if(n != numVectors){ printf("ERROR: vectorFigureID incorrect size\n");}

		//REMEMBER they are 1 based!
		vectorResidualID = reg.getIntArray("stats.vectorResidualID", &n);
		if(n != numVectors){ printf("ERROR: vectorResidualID incorrect size\n");}

		//get vector offsets -- REMEMBER they are 1 based!
		vectorOffsets = reg.getIntArray("stats.vectorOffsets", &n);
		if(n != numVectors){ printf("ERROR: vectorOffsets incorrect size\n");}

		vectors = new double *[numVectors];
		vectorLengths = new int [numVectors];

		for(i = 0; i < numVectors; i++){
			//get partial vector and save length
			vectors[i] = reg.getDoubleArray("stats.vectors[%d]", &vectorLengths[i], i);

			//Adjust 1 based numbers to zero based
			vectorResidualID[i]--;
			vectorOffsets[i]--;
		}

		vectorVariances = reg.getDoubleArray("stats.vectorVariances", &n);
		if(n != numVectors){ printf("ERROR: vectorVariances incorrect size\n");}

		numResiduals = reg.getIntValue("stats.numResiduals", 0);

		if(numResiduals > 1){
			residualVariances = reg.getDoubleArray("stats.residualVariances", &n);
 		  if(n != numResiduals){ printf("ERROR: residualVariances incorrect size\n");}
		}
		else{
			residualVariances = new double[1];
			residualVariances[0] = reg.getDoubleValue("stats.residualVariances", 0);
 		  if(residualVariances[0] == 0){ printf("ERROR: residualVariances[0] is 0\n");}
		}

		if(numResiduals > 0){

			// Need to compute the valid mask for each residual
			// Old code assumes that all entries of the vector are used with the residual
			// This was correct for jointMode=0 and "all joint" 
			// BUT INCORRECT for jointMode=1 when using gasBone or multiple features
			residualMasks = new bool * [numResiduals];
			residualMaskLengths = new int[numResiduals];
			for(i = 0; i < numResiduals; i++){
				residualMasks[i] = NULL;
				residualMaskLengths[i] = 0;

				//all vectors with the same residualID should be the same length. If not, error.
				for(j = 0; j < numVectors; j++){
					if(vectorResidualID[j] != i)
						continue;

					if(residualMasks[i] == NULL){
						residualMaskLengths[i] = vectorLengths[j];
						residualMasks[i] = new bool[residualMaskLengths[i]];

						for(k = 0; k < residualMaskLengths[i]; k++)
							residualMasks[i][k] = false;

					}

					if(vectorLengths[j] != residualMaskLengths[i]){
						printf("ERROR: SPE Read stats: vectors with same residualID have different lengths.\n");
						exit(1);
					}
					for(k = 0; k < residualMaskLengths[i]; k++){
						if(vectors[j][k] != 0.0)
							residualMasks[i][k] = true;
					}
				}

				//debug print residual mask
				/*
				printf("Residual %d:", i);
				for(k = 0; k < residualMaskLengths[i]; k++)
					printf(" %d", residualMasks[i][k]);
				*/
			}
		}
	}
	else{
		printf("WARNING: No SPE statistics supplied\n");
	}

	print();

	localLogFilename = 0;
}

void SurfacePatchEnsemble::print()
{
	printf("SURFACE PATCH ENSEMBLE SETTINGS:\n");
	printf("    matchType %d, numPatches %d\n", histogramMatchTypeToInt(hmt), numPatches);
	printf("  Region Definition\n");
	printf("    patchType %d, inDistance %.2lf, outDistance %.2lf, inGaussianSTD %.2lf, outGaussianSTD %.2lf\n",
							histogramPatchTypeToInt(histogramPatchType), inDistance, outDistance, inGaussianSTD, outGaussianSTD);
	printf("  Sampling variables\n");
	printf("    inSteps %d, outSteps %d, voxelSurfaceLevel %d, usePixelNeighbors %d\n",
							numInProfileSteps, numOutProfileSteps, voxelProfileSurfaceLevel, usePixelNeighbors);
	printf("  Local Region Definition\n");
	printf("    localSurfaceLevel %d, localRadius %.4lf, localSTD %.4lf\n", 
							localPatchSurfaceLevel, localPatchRadius, localPatchSTD);
	printf("  Gas and Bone\n");
	printf("    useLowHigh %d, lowThreshold %d, highThreshold %d\n", useLowHigh, lowThreshold, highThreshold);
	printf("  Bins\n");
	printf("    numLowBins %d, numMidBins %d, numHighBins %d\n", numLowBins, numMidBins, numHighBins);
	printf("  Statistics\n");
	printf("    numVectors %d, numResiduals %d\n", numVectors, numResiduals);
	printf("  QFData\n");
	printf("    numFeatures %d, dimensions (%d, %d, %d), margin (%d, %d, %d)\n", numFeatures, 
							QFDataDim[0], QFDataDim[1], QFDataDim[2], QFDataMargin[0], QFDataMargin[1], QFDataMargin[2]);
}

bool SurfacePatchEnsemble::save(const char * filename)
{

	// Adjust registry from creation to have current DATA folder
	// patchToFigure, patchToAtom, values 
	int i;
	int *patchToFigure = new int[numPatches];
	int *patchToAtom   = new int[numPatches];
	double *values2    = new double[numValues];
	for(i = 0; i < numPatches; i++){
		patchToFigure[i] = patches[i]->figureID;
		patchToAtom[i]   = patches[i]->atomID;
	}
	for(i = 0; i < numValues; i++){
		values2[i] = values[i];
	}
	
	// reg needs its own copy -- one isn't deleted later
	reg.setIntArray(   "matlabInfo.patchFigureID", numPatches, patchToFigure);
	reg.setIntArray(   "matlabInfo.patchAtomID",   numPatches, patchToAtom);
	reg.setDoubleArray("values",									 numValues,  values2);


	/*
	//Need to add back in ARRAYS in original reg (BOOO) - and reg needs its own copy
	for(f = 0; f < numFeatures; f++){
		double * v = new double [featureVectorLength];
		for(i = 0; i < featureVectorLength; i++){
			v[i] = featureVectors[f][i];
		}
		reg.setDoubleArray("config.featureVectors[%d]", featureVectorLength, v, f);
	}
	*/

	//Statistics if supplied will automatically be resaved! -- I don't think that is true

	bool retr = true;
	try {
		if (! reg.writeToFile(filename))
			cout << "Warning: Output SPE file may be corrupted" << endl;

	}
	catch (RException excp) {
		cout << excp.message() << endl;
		retr = false;
	}

	return retr;
}

SurfacePatchEnsemble::~SurfacePatchEnsemble()
{
	int f, p, i;
	//delete patches and pList and values
	if(values != NULL)
		delete [] values;

	if(pList != NULL){
		for(f = 0; f < numFigures; f++){
			if(pList[f] != NULL)
				delete pList[f];
		}
		delete [] pList;
	}

	for(p = 0; p < numPatches; p++){
		delete patches[p];
	}
	patches.clear();

	if (steps != NULL)
		delete [] steps;

	if (checkImage != NULL)
		delete checkImage;

	/*
	if(featureVectors != NULL){
		for(i = 0; i < numFeatures; i++){
			if(featureVectors[i] != NULL)
				delete [] featureVectors[i];
		}
		delete [] vectors;
	}
	*/
	if(QFData != NULL){
		delete [] QFData;
	}
	numFeatures = 0;

	if(mean != NULL)
		delete [] mean;

	if(scaleVector != NULL)
		delete [] scaleVector;

	if(vectors != NULL){
		for(i = 0; i < numVectors; i++){
			if(vectors[i] != NULL)
				delete [] vectors[i];
		}
		delete [] vectors;
	}

	numVectors = 0;

	if(vectorLengths != NULL)
		delete [] vectorLengths;

	if(vectorOffsets != NULL)
		delete [] vectorOffsets;

	if(vectorResidualID != NULL)
		delete [] vectorResidualID;

	if(vectorAtomID != NULL)
		delete [] vectorAtomID;

	if(vectorFigureID != NULL)
		delete [] vectorFigureID;

	if(vectorVariances!= NULL)
		delete [] vectorVariances;

	if(residualVariances!= NULL)
		delete [] residualVariances;

	if(residualMasks != NULL){
		for(i = 0; i < numResiduals; i++){
			if(residualMasks[i] != NULL)
				delete [] residualMasks[i];
		}
		delete [] residualMasks;
	}

	if(residualMaskLengths != NULL)
		delete [] residualMaskLengths;

	numResiduals = 0;
}

void SurfacePatchEnsemble::initializeLocalRegions(const M3DObject *object, Image3D *refImage)
{

	Bpoint *Bpoints;
	Bpoint *Bpoints2;
	int f, p, j;
	Vector3D p1, p2;
	double tempd;
	SurfacePatch *patch;
	int offset = 0;
	ProfileComponent pc;


	double localPatchRadiusSquared = localPatchRadius / refImage->getModelToWorldScale();
	localPatchRadiusSquared *= localPatchRadiusSquared;

	for(f = 0; f < numFigures; f++){

		int numFigPatches = 0;
		int numProfiles;

		//-> sets numFigPatches,
		// so change 
		getBPoints(object, f, localPatchSurfaceLevel, numFigPatches, Bpoints, true);

		////////////////////////////////////////////////////////////////////////////////
		//atomIDs - This should be the only non-tube compatible code
		//Need list of atomIDs given a list of level 0 boundary points
		const M3DFigure * figure = object->getFigurePtr(f);
		int *atomIDs = new int[numFigPatches];
		if( typeid(*figure) == typeid(M3DTubeFigure) ) {
			// tube figure?, atomid is nothing but u.
			// FIXME, TODO: for tubes atomId should be a combo of
			// atomId and spokeNumber.
			// (Maybe this variable should be renamed to spokeId)
			for(p = 0; p < numFigPatches; p++){
				atomIDs[p]	= int(Bpoints[p].u);
			}
		}
		else {
			const M3DQuadFigure* qFigure	= dynamic_cast<const M3DQuadFigure*>(figure);
			pList[f]->computeLvl0VertAtomList(qFigure->getRowCount(), qFigure->getColumnCount(), atomIDs);
		}
		////////////////////////////////////////////////////////////////////////////////

		getBPoints(object, f, voxelProfileSurfaceLevel, numProfiles, Bpoints2, false);

		numPatches += numFigPatches;
		for(p = 0; p < numFigPatches; p++){

#ifdef BINARY
			const bool discardEnds	= tuningWt(WeightMatchAtEnds) <= 1e-5;
#else
			const bool discardEnds	= int(tuningWt(IgnoreMatchAtEnds)) == 1;
#endif

			if( typeid(*figure) == typeid(M3DTubeFigure) &&
				discardEnds &&
				(Bpoints[p].u == 0 || Bpoints[p].u == figure->getPrimitiveCount() ) )
			{
				// Tube figure and patch corresponds to spokes from an end atom.
				numPatches--;
				continue;
			}

			patch = new SurfacePatch(f, atomIDs[p], offset, numFeatures);
			patches.push_back(patch);
			offset += 2 * (3 + numLowBins + numMidBins + numHighBins) * numFeatures;

			//Now find near profiles to each one
			p1.set(Bpoints[p].pnt);
			for (j = 0; j < numProfiles; j++) {
				p2.set(Bpoints2[j].pnt);
				tempd = p1.distSquare(p2);
				if (tempd <= localPatchRadiusSquared) {
					if (fabs(Bpoints[p].t) < 1 || Bpoints[p].t * Bpoints2[j].t > 0){
						pc.index = j;
						pc.weight = 1;
						if(localPatchSTD > 0.0)
							pc.weight = exp( -0.5 * tempd / (localPatchSTD * localPatchSTD));
						patch->profiles.push_back(pc);						
					}
				}
			}

			// Put the profiles in order, from center out so max. weight per voxel is used
			sort(patch->profiles.begin(), patch->profiles.end(), lessProfileComponent);

		}
		delete [] atomIDs;
		delete [] Bpoints;
	}
}

bool SurfacePatchEnsemble::getBPoints(const M3DObject *object, int figureID, int inSurfaceLevel,
                                      int &numPoints, Bpoint * &Bpoints, bool copy)
{

	int f;

	M3DFigure * figure = object->getFigurePtr(figureID);
	Xferlist * xferList = convertM3DtoXfer(figure);

	if(pList == NULL){
		pList = new ThallCode::Pointlist_server2 *[numFigures];
		for(f = 0; f < numFigures; f++){
			pList[f] = NULL;
		}
	}

	if (pList[figureID] == NULL) {
		//Requisite Incantations
		pList[figureID] = new ThallCode::Pointlist_server2;
		pList[figureID]->init(xferList);
		pList[figureID]->ComputeSubdivPointCloud(inSurfaceLevel);
	}
	else
		pList[figureID]->UpdateSubdivPointCloud(inSurfaceLevel, xferList);

	if(copy)
		pList[figureID]->copySubdivBoundaryInfo(&numPoints, &Bpoints);
	else
		pList[figureID]->subdivboundaryinfo(&numPoints, &Bpoints);

	delete [] (xferList->atomlist);
	delete xferList;

	//cout << numPoints << " boundary points used.\n";
	//boundary displacements not supported

	return true;
}

bool SurfacePatchEnsemble::update(Image3D * refImage,
								  const M3DObject * object,
								  int figureID, int atomID, 
								  const char * binaryImageFilename,
								  const char * patchImageFilename)
{
	printTime = printTime && ((double) rand() / RAND_MAX < 0.005);
	clock_t t1, t2, t3, t4;
	t1 = clock();

	returnWorstMatch = true;

	if(!initialized)
		initialize(object, refImage);

	int p, i, j, f;
	Image3D * patchImage = NULL;
	GreyValue * patchImageVoxels = NULL;
	GreyValue * patchIDValues = NULL;
	double worldToModelScale = 1 / refImage->getModelToWorldScale();

	// Clear Patches
	for(p = 0; p < numPatches; p++)
		patches[p]->clear();


	if (! object->verifyInBounds()) {
		cout << "Warning: This model is larger than or outside the unit cube --\n";
        //  << "    Not Erroring" << endl;
		return false;
	}


	// Init patch image
	if (patchImageFilename != NULL) {
		patchImage = new Image3D(refImage->getXDim(), refImage->getYDim(), refImage->getZDim());
        Vector3D temp = refImage->getWorldOrigin();
        patchImage->setSpacingAndOrigin(refImage->getXSpacing(), refImage->getYSpacing(), refImage->getZSpacing(), &temp);
		patchImage->clear();
		patchImageVoxels = patchImage->getVoxels();

		patchIDValues = new GreyValue[2 * numPatches];
		for (p = 0; p < 2 * numPatches; p++)
            patchIDValues[p] = 20000 + 30000 * rand();
	}

	// Initialize voxelCheck image
	if (checkImage == NULL) {
		checkImage = new Image3D(refImage->getXDim(), refImage->getYDim(), refImage->getZDim());
		checkImage->setSpacingAndOrigin(refImage->getXSpacing(), refImage->getYSpacing(), refImage->getZSpacing());
	}

	int startFigureID = figureID;
	int endFigureID = figureID;
	if(figureID == -1){
		startFigureID = 0;
		endFigureID = numFigures - 1;
	}

	t3 = clock();
	checkImage->clear(); 
	t4 = clock();
	if (printTime)
		printf("%.6lf seconds for check image clear\n", (t4-t3)/(double)CLOCKS_PER_SEC);

	bool useBinary = (binaryImageFilename != NULL);
	Image3D * binaryImage = NULL;

	// Load binary image (during training) to limit valid voxels
	if (useBinary) {
		RAWImageFile imageFile;

		binaryImage = imageFile.read(binaryImageFilename);

		if (binaryImage == NULL) {
			useBinary = false;
			cout << "Binary image " << binaryImageFilename << " FAILED to load.\n";
			return false;
		}
	}

	for(f = startFigureID; f <= endFigureID; f++){

		int numPoints, b;
		Bpoint * Bpoints;
		t3 = clock();
		getBPoints(object, f, voxelProfileSurfaceLevel, numPoints, Bpoints);
		t4 = clock();
		if (printTime)
			printf("%.6lf seconds for get Boundary Points\n", (t4-t3)/(double)CLOCKS_PER_SEC);


		int numNeighbors = 1;
		int xNeighbors[75] = {0, 0, -1,  1,  0, -1, -1,  1,  1,-2,-1, 0, 1, 2, 2, 2, 2, 2, 1, 0,-1,-2,-2,-2,-2, 0, 0, -1,  1,  0, -1, -1,  1,  1,-2,-1, 0, 1, 2, 2, 2, 2, 2, 1, 0,-1,-2,-2,-2,-2, 0, 0, -1,  1,  0, -1, -1,  1,  1,-2,-1, 0, 1, 2, 2, 2, 2, 2, 1, 0,-1,-2,-2,-2,-2};
		int yNeighbors[75] = {0,-1,  0,  0,  1, -1,  1, -1,  1, 2, 2, 2, 2, 2, 1, 0,-1,-2,-2,-2,-2,-2,-1, 0, 1, 0,-1,  0,  0,  1, -1,  1, -1,  1, 2, 2, 2, 2, 2, 1, 0,-1,-2,-2,-2,-2,-2,-1, 0, 1, 0,-1,  0,  0,  1, -1,  1, -1,  1, 2, 2, 2, 2, 2, 1, 0,-1,-2,-2,-2,-2,-2,-1, 0, 1};
		int zNeighbors[75] = {0, 0,  0,  0,  0,  0,  0,  0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1,  1,  1,  1,  1,  1,  1,  1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,-1,-1, -1, -1, -1, -1, -1, -1, -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1};

		if (usePixelNeighbors == 1)
			numNeighbors = 9;
		if (usePixelNeighbors == 2)
			numNeighbors = 75;


		// Compute normal steps
		int numProfileSteps = numInProfileSteps + numOutProfileSteps + 2;

		if (steps == NULL) {
			double modelInDistance  = inDistance * worldToModelScale;
			double modelOutDistance = outDistance * worldToModelScale;

			steps = new double[numProfileSteps];

			//inside steps - duplicate so always there when alternate
			steps[0] = 0.0;
			steps[1] = 0.0;

			if (numOutProfileSteps < numInProfileSteps) {
				cout << "Warning: numOutProfileSteps must be >= numInProfileSteps, changing OutSteps from "
                     << numOutProfileSteps << " to " << numInProfileSteps << ".\n";
				numOutProfileSteps = numInProfileSteps;
			}

			int offset;
			//printf("Offsets: ");
			for (i = 1; i <= numInProfileSteps; i++) {
				offset = 2 + ((i-1) / 2) * 4 + (i-1)%2;
				//printf(" %d", offset);
				steps[offset] = -1.0 * i / numInProfileSteps * modelInDistance;

				offset = 4 + ((i-1) / 2) * 4 + (i-1)%2;
				//printf(" %d", offset);
				steps[offset] = 1.0 * i / numOutProfileSteps * modelOutDistance;
			}

			//outside steps
			for (i = numInProfileSteps +1; i <= numOutProfileSteps; i++) {
				offset = 4 + ((numInProfileSteps-1) / 2) * 4 + (numInProfileSteps-1)%2 + i - numInProfileSteps;
				//printf(" %d", offset);
				steps[offset] = 1.0 * i / numOutProfileSteps * modelOutDistance;
			}
			//printf("\n");
			/*
              printf("%d Steps: ", numProfileSteps);
              for(i = 0; i < numProfileSteps; i++)
              printf("%.4lf ", steps[i]);
              printf("\n");
			*/

			cout << numPoints << " profiles x " << numProfileSteps/2 << " steps = "
                 << (numPoints * numProfileSteps/2) << endl;
		}

		// For each profile
		GreyValue * checkVoxels = checkImage->getVoxels();
		GreyValue * binaryVoxels;
		int index, vb;
		Vector3D coord, coord2, coordNew, coord3;
		int slices[3], slices2[3];
		int numPatchProfiles;
		bool inVoxel, inBounds;
		double *features;

		double profileWeight, pixelWeight, mahal;
		double modelInGaussianSTD  = inGaussianSTD  * worldToModelScale;
		double modelOutGaussianSTD = outGaussianSTD * worldToModelScale;
		double modelInGaussianSTDMultiplier = -0.5 / modelInGaussianSTD / modelInGaussianSTD;
		double modelOutGaussianSTDMultiplier = -0.5 / modelOutGaussianSTD / modelOutGaussianSTD;

		// Used to prune out pixels outside collar - for consistency but not needed
		double inMaxDistanceSquared = inDistance * inDistance * worldToModelScale * worldToModelScale;
		double outMaxDistanceSquared = outDistance * outDistance * worldToModelScale * worldToModelScale;

		if (useBinary)
			binaryVoxels = binaryImage->getVoxels();

		// Time main loop
		t3 = clock();

		//Local and global Regions
		for (p = 0; p < numPatches; p++) {
			
			if (f != patches[p]->figureID || (atomID != -1 && atomID != patches[p]->atomID))
				continue;

			numPatchProfiles = patches[p]->profiles.size();
			
			//printf("Patch %d with %d profiles\n", p, numPatchProfiles);
			//printf("%d Steps: ", numProfileSteps);
			//for (int tt = 0; tt < numProfileSteps; tt++)
			//	printf("%.3lf ", steps[tt]);
			//printf("\n");
			
			int numInBounds = 0, numInBounds2 = 0, numCheck = 0, numIn = 0, numOut = 0, numIn2 = 0, numOut2 = 0;

			for (j = 0; j < numProfileSteps/2; j++) {
				for(vb = 0; vb < numPatchProfiles; vb++) {

					profileWeight = patches[p]->profiles[vb].weight;

					b = patches[p]->profiles[vb].index;
					i = j*2 + b%2;
					Vector3D Pi(Bpoints[b].pnt);
					Vector3D Ni(Bpoints[b].norm);

					coord = Pi + Ni * steps[i];

					refImage->modelToImageCoordinates(coord, slices2);
					inBounds = refImage->clipToImage(slices2);
					if (!inBounds)
						continue;

					for (int ii = 0; ii < numNeighbors; ii++) {          
						slices[0] = slices2[0] + xNeighbors[ii];
						slices[1] = slices2[1] + yNeighbors[ii];
						slices[2] = slices2[2] + zNeighbors[ii];

						inBounds = refImage->imageCoordinatesInBounds(slices);

						if(!inBounds)
							continue;
						numInBounds++;

						index = refImage->imageIndex(slices);

						//check that the QFData is valid for this pixel (i.e. that its neighborhood is valid)
						inBounds = checkQFData(slices);
						if(!inBounds)
							continue;
						numInBounds2++;

						if (checkVoxels[index] == (p + 1))
							continue;
						checkVoxels[index] = (p + 1);
						numCheck++;

						refImage->imageToModelCoordinates(slices, coord2);
                        mahal = (Ni * (coord2 - Pi));
						inVoxel = (mahal <= 0);

						if (useBinary)
							inVoxel = (binaryVoxels[index] > 0);

						mahal *= mahal;
						pixelWeight = 1.0;

						if (inVoxel) {
							numIn++;
							if (mahal >= inMaxDistanceSquared)
								continue;
							numIn2++;
							if (inGaussianSTD > 0.0)
								pixelWeight = exp(mahal * modelInGaussianSTDMultiplier);

							features = getVoxelFeatures(refImage, index);
							addElement(p, features, profileWeight, pixelWeight, true);

							if (patchImageVoxels != NULL)
								patchImageVoxels[index] = patchIDValues[2*p];

						}
						else {
							numOut++;
							if (mahal >= outMaxDistanceSquared)
								continue;
							numOut2++;
							if (outGaussianSTD > 0.0)
								pixelWeight = exp(mahal * modelOutGaussianSTDMultiplier);							

							features = getVoxelFeatures(refImage, index);
							addElement(p, features, profileWeight, pixelWeight, false);

							if (patchImageVoxels != NULL)
								patchImageVoxels[index] = patchIDValues[2*p+1];
						}
					}
				}
			}
			if (patchImageFilename != NULL) {
				printf("numInBounds = %d, numInBounds2 = %d, numCheck = %d, numIn = %d, numIn2 = %d, numOut = %d, numOut2 = %d\n",
                       numInBounds, numInBounds2, numCheck, numIn, numIn2, numOut, numOut2);
			}
			if (!patches[p]->simplify(values, mean, numLowBins, numMidBins, numHighBins)) {
				cout << "ERROR: patches->simplify()" << endl;
				return false;
			}
		}

		t4 = clock();
		if (printTime)
			printf("%.6lf seconds for main loop\n", (t4-t3)/(double)CLOCKS_PER_SEC);

	}
	//save patchImage
	if (patchImage != NULL) {
		RAWImageFile imageFile;

		if (imageFile.write(patchImageFilename, *patchImage))
			cout << "Saved Patch Image: " << patchImageFilename << endl;
		else
			cout << "Patch Image failed to save: " << patchImageFilename << endl;

		delete patchImage;
	    delete [] patchIDValues;
		patchImage = NULL;
		patchIDValues = NULL;
	}

	if (binaryImage != NULL) {
		delete binaryImage;
		binaryImage = NULL;
	}

	t2 = clock();
	if (printTime)
		printf("%.6lf seconds for update()\n", (t2-t1)/(double)CLOCKS_PER_SEC);

	returnWorstMatch = false;
	return true;
}

bool SurfacePatchEnsemble::checkQFData(int * slices){
	//						inBounds = checkQFData(slices);

	if(QFData == NULL)
		return true;

	if(slices[0] - QFDataMargin[0] < 0 || slices[0] + QFDataMargin[0] >= QFDataDim[0])
		return false;
	if(slices[1] - QFDataMargin[1] < 0 || slices[1] + QFDataMargin[1] >= QFDataDim[1])
		return false;
	if(slices[2] - QFDataMargin[2] < 0 || slices[2] + QFDataMargin[2] >= QFDataDim[2])
		return false;

	return true;

}

double * SurfacePatchEnsemble::getVoxelFeatures(Image3D * refImage, int index){

	double * features;

	static double dummyFeature;

	if (QFData == NULL) {
		dummyFeature = (double) refImage->mapDisplayToActual(refImage->getVoxelValue(index));
		features = &dummyFeature;
		return features;
	}

	features = &(QFData[index * numFeatures]);

	//temp debug
	/*
	int i;
	static int numPrints = 0;
	if(numPrints == 0)
		printf("Extent: %d %d %d\n", refImage->getXDim(), refImage->getYDim(), refImage->getZDim());

	if(numPrints % 10 == 0 && numPrints < 1000){
		printf("(%d) -> %lf: ", index, (double) refImage->mapDisplayToActual(refVoxels[index]));
		for(i = 0; i < numFeatures; i++)
			printf("%lf ", features[i]);
		printf("\n");
	}
	numPrints++;
  */

	return features;
}

void SurfacePatchEnsemble::addElement(int patchID, double *features, double profileWeight, double pixelWeight, bool isInside){

	bool isLow = false;
	bool isHigh = false;


	//The definition of thresholds with multiple features needs work
	//This code assumes gas/bone is only turned on with single features


	//inside gas/bone should not use pixel weight? -- TESTING
	double weight = profileWeight * pixelWeight;

	if (useLowHigh && features[0] < lowThreshold) {
		isLow = true;
		//if(isInside)
		//	weight = profileWeight;
	}
	else if (useLowHigh && features[0] > highThreshold) {
		isHigh = true;
		//if(isInside)
		//	weight = profileWeight;
	}

	patches[patchID]->addElement(features, weight, isInside, isLow, isHigh);	

	//delete [] features;
	//features = NULL;
}

double SurfacePatchEnsemble::getMatchValue(int figureID, int atomID)
{
	double match = 0.0;
	int v, i, r;
	int numUsedVectors = 0;
	double m;
	static bool firstTime = false; //true

	if(!initialized){
		printf("ERROR: Getting SPE Match without being initialized\n");
		return 999;
	}

	if(returnWorstMatch)
		return 999;

	std::ofstream logfile;
	if (localLogFilename) {
	  cout << "Saving local image match to " << localLogFilename << endl;
	  logfile.open(localLogFilename);
  }

	//initialize Residuals
	bool   *residualInit     = new bool[numResiduals];
	double *residuals        = new double[numResiduals];
	int    *residualFigureID = new int[numResiduals];
	int    *residualAtomID   = new int[numResiduals];

	for(i = 0; i < numResiduals; i++){
		residualInit[i]     = false;
		residuals[i]        = 0.0;
		residualFigureID[i] = -1;
		residualAtomID[i]   = -1;
	}

	//Subtract the mean and scale
	for(v = 0; v < numValues; v++)
		values[v] = (values[v] - mean[v]) * scaleVector[v];

	// Apply PCA Vectors
	for(v = 0; v < numVectors; v++){
		if(figureID != -1 && vectorFigureID[v] != -1 && figureID != vectorFigureID[v])
			continue;
		if(atomID != -1 && vectorAtomID[v] != -1 && atomID != vectorAtomID[v])
			continue;

		numUsedVectors++;
		if(firstTime)
			printf("Num used vectors: %d\n", numUsedVectors);

		m = 0;
		for(i = 0; i < vectorLengths[v]; i++)
			m += values[vectorOffsets[v] + i] * vectors[v][i];

		m *= m;
	
		//setup Residual and subtract out amount traveled
		r = vectorResidualID[v];
		if(r >= 0){
			if(!residualInit[r]){
				residualInit[r] = true;
				residualFigureID[r] = vectorFigureID[v];
				residualAtomID[r]   = vectorAtomID[v];
				for(i = 0; i < vectorLengths[v]; i++)
					residuals[r] += values[vectorOffsets[v] + i] * values[vectorOffsets[v] + i] * residualMasks[r][i];
			}
			residuals[r] -= m;
		}

		//add mahal distance
		m /= vectorVariances[v];
		match += m;

		if(firstTime)
			printf("Vector %d: %lf\n", v, m);

    if (localLogFilename) {
			logfile << vectorFigureID[v] << ' ' << vectorAtomID[v] << ' ' << m << endl;
    }
	}

	// Apply Residuals
	for(v = 0; v < numResiduals; v++){
		if(!residualInit[v])
			continue;

		numUsedVectors++;

		m = residuals[v] / residualVariances[v];

		match += m;

		if(firstTime)
			printf("Residual %d: %lf\n", v, m);
    if (localLogFilename) {
			logfile << residualFigureID[v] << ' ' << residualAtomID[v] << ' ' << m << endl;
    }
	}

	delete [] residuals;
	delete [] residualInit;
	delete [] residualFigureID;
	delete [] residualAtomID;
	
	if(hmt == WEIGHTED_GAUSSIAN_HMT){
		
		double matchVariance = numUsedVectors;

		//Pablo Shouldn't need this info from the registry, but for now...
		// Weight needs to be adjusted for outside gas/bone that is downweighted
		int jointMode        = reg.getIntValue("config.jointMode", -1);
		double inLowWeight   = reg.getDoubleValue("config.inLowWeight",  1.0);
		double inHighWeight  = reg.getDoubleValue("config.inHighWeight",  1.0);
		double outLowWeight	 = reg.getDoubleValue("config.outLowWeight",1.0);
		double outHighWeight = reg.getDoubleValue("config.outHighWeight",1.0);

		if(jointMode <= 2){

			//COMPUTE number of patches for figureID
			int numFigPatches = 0;
			for(i = 0; i < numPatches; i++){
				if(patches[i]->figureID == figureID)
					numFigPatches++;
			}
			if(firstTime)
				printf("numUsedVectors: %d\n", numUsedVectors);

			//Need to change #vectors for when Weight is NON zero -- otherwise vector wasn't included
			if(useLowHigh == 1){
				if(inLowWeight > 0)
					matchVariance += numFigPatches * (inLowWeight - 1);
				if(inHighWeight > 0)
					matchVariance += numFigPatches * (inHighWeight - 1);
				if(outLowWeight > 0)
					matchVariance += numFigPatches * (outLowWeight - 1);
				if(outHighWeight > 0)
					matchVariance += numFigPatches * (outHighWeight - 1);
			}
		}
		
		if(firstTime)
			printf("Normalizing by: %lf\n", matchVariance);
		cout.flush();
		
		match = (match - matchVariance) / sqrt(2.0 * matchVariance);
	}
  
	if(firstTime){
		printf("Match: %lf\n", match);
		cout.flush();
	}
	firstTime = false;

  if (localLogFilename)
		logfile.close();

	return match;
}

int SurfacePatchEnsemble::histogramMatchTypeToInt(HistogramMatchType hmt)
{
	if (hmt == GAUSSIAN_HMT)
		return 0;
	if (hmt == WEIGHTED_GAUSSIAN_HMT)
		return 1;
	return -1;
}

HistogramMatchType SurfacePatchEnsemble::intToHistogramMatchType(int hmtInt)
{
	if (hmtInt == 0)
		return GAUSSIAN_HMT;
	if (hmtInt == 1)
		return WEIGHTED_GAUSSIAN_HMT;

	cout << "ERROR: SPE::intToHistogramMatchType: Unknown match type\n";
	return GAUSSIAN_HMT;
}

M3DObject * SurfacePatchEnsemble::getSingleFigureObject(const M3DObject * oldObject, int figureID)
{
	// Break the cast: this is OK, since the returned object is temporary
	M3DObject * object = ((M3DObject *) oldObject)->assign();

	// If fig != -1 generate new object of just that figure
	if (figureID >= 0) {
		int nFigures = object->getFigureCount();
		for (int i = 0; i < nFigures; i++) {
			if (i < figureID) {
				object->deleteFigure(i);
				figureID--;
				i--;
				nFigures--;
			}
			if (i > figureID) {
				object->deleteFigure(i);
				i--;
				nFigures--;
			}
		}	
	}

	return object;
}

int SurfacePatchEnsemble::histogramPatchTypeToInt(HistogramPatchType hpt)
{
	if (hpt == GLOBAL_HPT)
		return 0;
	if (hpt == LOCAL_HPT)
		return 1;
	return -1;
}

HistogramPatchType SurfacePatchEnsemble::intToHistogramPatchType(int hptInt) 
{
	if (hptInt == 0)
		return GLOBAL_HPT;
	if (hptInt == 1)
		return LOCAL_HPT;
	return GLOBAL_HPT;
}

/************************************************************************************
  SurfacePatch
************************************************************************************/
SurfacePatch::SurfacePatch(int inFigureID, int inAtomID, int inOffset, int inNumFeatures){
	figureID = inFigureID;
	atomID = inAtomID;
	offset = inOffset;
	numFeatures = inNumFeatures;

	profiles.clear();
	profiles.reserve(1000);

	mixtures = new QuantileMixture[2 * numFeatures];

	clear();
}

SurfacePatch::~SurfacePatch(){
	profiles.clear();
}

void SurfacePatch::clear()
{
	int i;
	for(i = 0; i < 2*numFeatures; i++){
		mixtures[i].lowFrequency = 0.0;
		mixtures[i].midFrequency = 0.0;
		mixtures[i].highFrequency = 0.0;

		mixtures[i].lowComponents.clear();
		mixtures[i].lowComponents.reserve(1000);
		mixtures[i].midComponents.clear();
		mixtures[i].midComponents.reserve(10000);
		mixtures[i].highComponents.clear();
		mixtures[i].highComponents.reserve(1000);
	}
}

bool SurfacePatch::addElement(double *features, double weight, bool isInside, bool isLow, bool isHigh){

	QuantileComponent qc;

	int mixtureOffset = 0;
	int i;

	if(!isInside){
		mixtureOffset = numFeatures;
	}


	//Here low and high checks are presumed -- and are based only on features[0]

	qc.weight = weight;
	for(i = 0; i < numFeatures; i++){

		qc.value = features[i];

		if(isLow){
			mixtures[i+mixtureOffset].lowFrequency += weight;
			mixtures[i+mixtureOffset].lowComponents.push_back(qc);
		}
		else if(isHigh){
			mixtures[i+mixtureOffset].highFrequency += weight;
			mixtures[i+mixtureOffset].highComponents.push_back(qc);
		}
		else{
			mixtures[i+mixtureOffset].midFrequency += weight;
			mixtures[i+mixtureOffset].midComponents.push_back(qc);
		}
	}

	return true;
}

// Convert rawValues into values
bool SurfacePatch::simplify(double *values, double * mean, int numLowBins, int numMidBins, int numHighBins)
{
	bool ret = true;
	int i;
	int currentOffset = offset;

	//sorting is slower but required if non-integer data. assume 1 feature = intensity = int. Otherwise sort.
	bool useSort = false;
	//if(numFeatures == 1) useSort = false;

	for(i = 0; i < 2*numFeatures; i++){

		//normalize the Frequencies -- cannot apply until after simply though!
		double total  =  mixtures[i].lowFrequency +  mixtures[i].midFrequency +  mixtures[i].highFrequency;
		
		ret = ret && simplifyHelper(values, mean, numLowBins, currentOffset, mixtures[i].lowFrequency, total, mixtures[i].lowComponents, useSort);
		currentOffset += 1 + numLowBins;

		ret = ret && simplifyHelper(values, mean, numMidBins, currentOffset, mixtures[i].midFrequency, total, mixtures[i].midComponents, useSort);
		currentOffset += 1 + numMidBins;

		ret = ret && simplifyHelper(values, mean, numHighBins,currentOffset, mixtures[i].highFrequency, total, mixtures[i].highComponents, useSort);
		currentOffset += 1 + numHighBins;
	}

	clear();
	return ret;
}

bool SurfacePatch::simplifyHelper(double *values, double *mean, int numBins, int offset, double count, double total, 
																	std::vector<QuantileComponent> &components, bool useSort){
																												
	// Construct quantiles - while *Frequencies are still weights
	int i, b;
	int numE = components.size();

	double wpb = count / numBins;	// Weight per Bin
	double currentWeight, partialWeight, w;
	int currentBin;
	
	values[offset] = count / total;

	offset++; //Nothing else cares about the frequency/count!

	//Zero out / initialize the current values
	for (b = 0; b < numBins; b++)
		values[offset + b] = 0;

	if (count == 0.0) {
		//If there is no weight and a MEAN is given -- set to the mean so no penalty
		if(mean != NULL){
			for (b = 0; b < numBins; b++)
				values[offset + b] = mean[b];
		}
		return true;
	}

	//do the easy case fast
	if(numBins == 1){
		double v = 0;
		for (i = 0; i < numE; i++) {
			v += components[i].value * components[i].weight;
		}
		values[offset] = v / wpb;
		return true;
	}

	if(useSort){
		sort(components.begin(), components.end(), lessQuantileComponent);
	
		// Fill test values 
		currentBin = 0;
		currentWeight = 0;
		for (i = 0; i < numE; i++) {
			w = components[i].weight;
			while (w > 0.0) {
				//if point falls entirely in bin
				if ((currentBin == numBins-1) || (currentWeight + w <= wpb)) {
					currentWeight += w;
					values[offset + currentBin] += w * components[i].value;
					w = 0.0;
				}
				else {
					// else point goes into next bins too
					partialWeight = wpb - currentWeight;
					values[offset + currentBin] += partialWeight * components[i].value;

					currentWeight = 0.0;
					w -= partialWeight;
					currentBin++;
				}
			}
		}
	}
	else{

		// NEW METHOD - only works for integer valued image intensities

		// Find the min and max values
		int v = (int) (components[0].value + 0.5);
		int minBin = v;
		int maxBin = v;
		for (i = 1; i < numE; i++) {
			v = (int) (components[i].value + 0.5);
			if (v < minBin)
				minBin = v; 
			else if (v > maxBin)
				maxBin = v;
		}

		int numFreqBins = maxBin - minBin + 1;

		// Create Frequency bins
		double *freqBins = new double[numFreqBins];
		for (i = 0; i < numFreqBins; i++)
			freqBins[i] = 0.0;

		// Fill frequency bins
		for (i = 0; i < numE; i++) {
			b = ((int)(components[i].value + 0.5)) - minBin;
			freqBins[b] += components[i].weight;
		}

		// Fill test values - a lot like old method
		currentBin = 0;
		currentWeight = 0;
		for (i = 0; i < numFreqBins; i++) {
			w = freqBins[i];
			while (w > 0.0) {
				//if point falls entirely in bin
				if ((currentBin == numBins-1) || (currentWeight + w <= wpb)) {
					currentWeight += w;
					values[offset + currentBin] += w * (minBin + i);
					w = 0.0;
				}
				else {
					// else point goes into next bins too
					partialWeight = wpb - currentWeight;
					values[offset + currentBin] += partialWeight * (minBin + i);

					currentWeight = 0.0;
					w -= partialWeight;
					currentBin++;
				}
			}
		}

		// Clear bins
		delete [] freqBins;
		freqBins = NULL;
	}

	// Normalize the bins back to values
	double wpbInv = 1.0 / wpb;
	for (b = 0; b < numBins; b++)
		values[offset + b] *= wpbInv;

	return true;
}
