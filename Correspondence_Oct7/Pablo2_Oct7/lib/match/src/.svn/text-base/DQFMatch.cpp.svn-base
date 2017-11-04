//Distance Quantile Function image match stuff is here and in the .h.

#include <iostream>
#include <fstream>
#include "DQFMatch.h"
#include "P3DControl.h"

#include <algorithm>


/********************************************************************************************************************************************************/



/********************************************************************************************************************************************************
Code for the DQuantileFunction class, for easy use with the DQFTraining and DQFMatch classes.
*********************************************************************************************************************************************************/




/********************************************************************************************************************************************************/




DQuantileFunction::DQuantileFunction()
{
	rawData.clear();
	setQFSize(10);

	useMyWeighting = false;
	dirty = true;
	qf = NULL;
}

DQuantileFunction::DQuantileFunction(int size)
{
	rawData.clear();
	qf = NULL;

	setQFSize(size);

	useMyWeighting = false;
	dirty = true;
}

DQuantileFunction::~DQuantileFunction()
{
	rawData.clear();
	delete []qf;
}

void DQuantileFunction::useDefaultWeighting(double umu, double uvar)
{
	useMyWeighting = true;
	mu = umu;
	var = uvar;
}

void DQuantileFunction::setQFSize(int size)
{
	qfSize = size; 
	dirty = true;

	if (qf != NULL) 
		delete []qf;

	qf = new double[qfSize];
	memset(qf, 0, sizeof(double)*qfSize);

}

void DQuantileFunction::addData(DDataElement element)
{
	rawData.push_back(element);
	dirty = true;
}

void DQuantileFunction::addData(double element)
{
	DDataElement de;

	de.element = element;

	//Now for the weight.  If the user sent a double and not a DDataElement, 
	//then they get constant weighting unless useMyWeighting is true, in
	//which case, the weight should be computed here.

	if (useMyWeighting)
	{
		//Here we use the weighting 
		
		if (mu <= 0)
		{
			if (element >= mu && element <= -mu)
				de.weight = 1.0;
			else
				de.weight = exp((element - mu)*(element - mu)/(-2.0*var));
		}
		else if (mu >= 0)
		{
			if (element <= mu && element >= -mu)
				de.weight = 1.0;
			else
				de.weight = exp((element - mu)*(element - mu)/(-2.0*var));
		}
			

	}
	else
	{
		de.weight = 1.0;
	}

	rawData.push_back(de);
	dirty = true;
}

//This must be defined for the sort of the vector.
bool operator<(const DDataElement& a, const DDataElement& b) 
{
    return a.element < b.element;
}

void DQuantileFunction::makeQF(int size)
{
	int i, numE = rawData.size();
	if (numE < qfSize)
		cerr << "Not enough data to makeQF." << endl;

	if (size == qfSize)
		memset(qf, 0, sizeof(double)*qfSize);
	else
		setQFSize(size);


	//The relatively slow way: This way loops through the rawData more
	//than once. A potentially much better way (assuming there are lots
	//of elements in rawData) would be to make a histogram at, say, 4*qfSize
	//resolution, and then loop through that.  But I can always come back.
	std::sort(rawData.begin(), rawData.end());

	//Now rawData is sorted.  Take the average of each proportion of weight.
	//First find the sum of the weights.
	double sum = 0.0, wpb = 1.0, wpbInv = 1.0;

	for (i = 0; i < numE; i ++)
	{
		sum += rawData[i].weight;
	}
//cout << sum << ".\n";
	if (sum < 1.0e-41)
	{
		cout << "Total weight quite small: " << sum << ".\n";

		//This basically means that the thousands of bone voxels are so far away
		//as to have a miniscule effect on the result.

		
	}
	//else
	//{
		wpb = sum/qfSize;
		wpbInv = 1.0/wpb;
	//}
	
	//loop through the data again.
	int currentBin = 0;
	double partialWeight, currentWeight = 0.0;
	i = 0;
	
	while (currentBin < qfSize && i < numE)
	{
		while (currentWeight + rawData[i].weight < wpb && i < numE)
		{
			qf[currentBin] += rawData[i].weight * rawData[i].element;
			currentWeight += rawData[i].weight;
			i++;
		}
		//get the partial weight and move on, unless i >= numE.
		if (i < numE)
		{
			partialWeight = wpb - currentWeight;
			qf[currentBin] += partialWeight * rawData[i].element;
			
			partialWeight = rawData[i].weight - partialWeight;
			//what if the element's weight was greater than that of a bin, for those
			//really far off voxels...
			while (partialWeight > wpb)
			{
				partialWeight -= wpb;
				qf[++currentBin] += wpb * rawData[i].element;
			}
			qf[++currentBin] += partialWeight * rawData[i].element;
			
			currentWeight = partialWeight;
			i ++;
		}
	}

	for (currentBin = 0; currentBin < qfSize; currentBin ++)
	{
		qf[currentBin] *= wpbInv;
	}
	
}

double * DQuantileFunction::getQF()
{
	if (!dirty)
		return qf;

	makeQF(qfSize);

	dirty = false;
	return qf;
}

double DQuantileFunction::getMean()
{
	if (dirty)
		makeQF(qfSize);

	double sum = 0.0;

	for (int i = 0; i < qfSize; i ++)
	{
		sum += qf[i];
	}

	return sum/qfSize;
}














/********************************************************************************************************************************************************/



/********************************************************************************************************************************************************
Code for the DQFTraining class, for creating the dqf images and model dqfs for training
*********************************************************************************************************************************************************/




/********************************************************************************************************************************************************/





DQFTraining::DQFTraining()
{
	//boundingBoxBuffer = 0.0;
	Xmin.set(0.0,0.0,0.0);
	Xmax.set(0.0,0.0,0.0);
	
	Registry reg;
	loadConfig(reg);
}

DQFTraining::DQFTraining(Registry & reg)
{
	Xmin.set(0.0,0.0,0.0);
	Xmax.set(0.0,0.0,0.0);

	loadConfig(reg);
}

void DQFTraining::loadConfig(Registry & reg)
{
	sWidth             = reg.getIntValue("dqfSWidth", 10);

	boundingBoxBuffer  = reg.getDoubleValue("dqfBoundingBoxBuffer", 5.0);

	weightingMu		   = reg.getDoubleValue("dqfWeightingMu", 2.0);
	weightingVar       = reg.getDoubleValue("dqfWeightingVar", 1.0);
	modelSubDivLevel   = reg.getIntValue("dqfModelSurfaceLevel", 0);

	imageDelta         = reg.getIntValue("dqfImageDelta", 0.0);
	imageMultiplier    = reg.getIntValue("dqfImageMultiplier", 10.0);
}

void DQFTraining::setDQFImageBoundingBox(Vector3D uXmin, Vector3D uXmax)
{
	Xmin = uXmin;
	Xmax = uXmax;
}

void DQFTraining::setBoundingBoxBuffer(double buffer)
{
	boundingBoxBuffer = buffer;
}

//This function will determine the size of a bounding box to be used
//in construction of the dqf image later.
void DQFTraining::setDQFImageBoundingBox(Image3D *image, const M3DObject * object, int figureID)
{

	double maxx,maxy,maxz,minx,miny,minz;

	M3DFigure * figure = object->getFigurePtr(figureID);
    
	Xferlist * xferList = convertM3DtoXfer(figure);

	ThallCode::Pointlist_server2 * pList = new ThallCode::Pointlist_server2;
	pList->init(xferList);
	pList -> ComputeSubdivPointCloud(0);


	int numPoints;
	Bpoint * Bpoints;
    pList->subdivboundaryinfo(&numPoints, &Bpoints);

	Vector3D Pi(Bpoints[0].pnt);
	image->modelToImageCoordinates(Pi);

	minx = maxx = Pi.getX();
	miny = maxy = Pi.getY();
	minz = maxz = Pi.getZ();

	for (int i = 1; i < numPoints; i ++)
	{
		Pi.set(Bpoints[i].pnt); //the figure boundary point.
		image->modelToImageCoordinates(Pi);

		if (Pi.getX() < minx)
			minx = Pi.getX();
		if (Pi.getX() > maxx)
			maxx = Pi.getX();

		if (Pi.getY() < miny)
			miny = Pi.getY();
		if (Pi.getY() > maxy)
			maxy = Pi.getY();

		if (Pi.getZ() < minz)
			minz = Pi.getZ();
		if (Pi.getZ() > maxz)
			maxz = Pi.getZ();

	}

	//now set the bounding box to be a bit larger than those points.
	//how much larger is from boundingBoxBuffer.  Also, boundingBoxBuffer is in cm, whereas
	//the Xmin and Xmax are in voxels.  So I have to convert, then add, then convert back...

	Xmin.set(minx, miny, minz);
	Xmax.set(maxx, maxy, maxz);


	//image->imageToWorldCoordinates(Xmin);
	//image->imageToWorldCoordinates(Xmax);

	//Xmin -= boundingBoxBuffer;
	//Xmax += boundingBoxBuffer;

	//image->worldToImageCoordinates(Xmin);
	//image->worldToImageCoordinates(Xmax);
	
	//Now fix in case the y's got switched...  
	//if (Xmax.getY() < Xmin.getY())
	//{
	//	miny = Xmax.getY();
	//	Xmax.setY(Xmin.getY());
	//	Xmin.setY(miny);
	//}


	Xmin.setX(Xmin.getX() - image->worldXToImageDistance(boundingBoxBuffer));
	Xmin.setY(Xmin.getY() - image->worldYToImageDistance(boundingBoxBuffer));
	Xmin.setZ(Xmin.getZ() - image->worldZToImageDistance(boundingBoxBuffer));

	Xmax.setX(Xmax.getX() + image->worldXToImageDistance(boundingBoxBuffer));
	Xmax.setY(Xmax.getY() + image->worldYToImageDistance(boundingBoxBuffer));
	Xmax.setZ(Xmax.getZ() + image->worldZToImageDistance(boundingBoxBuffer));


	//round to actual voxel coords.
	Xmin.round();
	Xmax.round();

	//Cut these coordinates to the image.
	if (!image->imageCoordinatesInBounds(Xmin.getX(), Xmin.getY(), Xmin.getZ()) ||
		!image->imageCoordinatesInBounds(Xmax.getX(), Xmax.getY(), Xmax.getZ()))
	{
		if (Xmin.getX() < 0)
			Xmin.setX(0);
		
		if (Xmax.getX() >= image->getXDim())
			Xmax.setX(image->getXDim()-1);

		if (Xmin.getY() < 0)
			Xmin.setY(0);
		
		if (Xmax.getY() >= image->getYDim())
			Xmax.setY(image->getYDim()-1);

		if (Xmin.getZ() < 0)
			Xmin.setZ(0);
		
		if (Xmax.getZ() >= image->getZDim())
			Xmax.setZ(image->getZDim()-1);
	}

	delete [] (xferList->atomlist);
	delete xferList;
	delete pList; 

	//Now Xmin and Xmax are image coords.

	cout << "DQF image bounding box, from "; 
	Xmin.print();
	cout << " to "; 
	Xmax.print();
	cout << ".\n";
}


void DQFTraining::floodFill(Image3D *image, int x, int y, int z, GreyValue &targetGrey, GreyValue &replacementGrey)
//This probably takes way too much memory.  The program just quits. I have no idea why. 
//So I'm messing with the implementation details until I get it to work. 
//got it to work: with while loops that don't belong in a recursive scheme. but anyway.
{
	if (x == 221 && y == 278 && z == 61)
		int qq = 0;

	if (!image->imageCoordinatesInBounds(x,y,z)) 
		return;

	if (image->getVoxelValue(x,y,z) != targetGrey) //&& image->getVoxelValue(x,y,z) != replacementGrey 
	{
		return;
	}
	
	
	int startx = x;
	while (image->imageCoordinatesInBounds(x,y,z) && image->getVoxelValue(x,y,z) == targetGrey)
	{
		image->setVoxel(x,y,z,replacementGrey);
		floodFill(image, x, y+1, z, targetGrey, replacementGrey);
		
		x ++;
	}
	int posx = x - 1;
	
	x = startx - 1;
	
	while (image->imageCoordinatesInBounds(x,y,z) && image->getVoxelValue(x,y,z) == targetGrey)
	{
		image->setVoxel(x,y,z,replacementGrey);
		floodFill(image, x, y+1, z, targetGrey, replacementGrey);
		
		x --;
	}
	int negx = x + 1;
	
	for (x = negx; x <= posx; x ++)
		floodFill(image, x, y-1, z, targetGrey, replacementGrey);
	
	
	x = startx;
	//floodFill(image, x+1, y, z, targetGrey, replacementGrey);
	//cout << "1";
	floodFill(image, x-1, y, z, targetGrey, replacementGrey);
	//cout << "2";
	floodFill(image, x, y+1, z, targetGrey, replacementGrey);
	//cout << "3";
	floodFill(image, x, y-1, z, targetGrey, replacementGrey);
	

	return;
}

int DQFTraining::conCompSize(Image3D *target, Image3D *seen, int x, int y, int z, int & maxComponentSize)
//Returns the size of the MAX_GREY_VALUE connected component including voxel x,y,z in the 
//z slice.
{
	if (seen->getVoxelValue(x,y,z))
		return 0;

	if (!target->getVoxelValue(x,y,z))
	{
		return 0;
	}
	else
	{
		seen->setVoxel(x,y,z,1);

		//Can't quit early, cause then seen isn't appropriately updated.
		int north = conCompSize(target, seen, x, y+1, z, maxComponentSize);
		//if (north > maxComponentSize)
		//	return maxComponentSize;

		int south = conCompSize(target, seen, x, y-1, z, maxComponentSize);
		//if (south > maxComponentSize)
		//	return maxComponentSize;

		int east = conCompSize(target, seen, x+1, y, z, maxComponentSize);
		//if (east > maxComponentSize)
		//	return maxComponentSize;

		int west = conCompSize(target, seen, x-1, y, z, maxComponentSize);
		//if (west > maxComponentSize)
		//	return maxComponentSize;

		return 1 + north + south + east + west;

	}
}

void DQFTraining::deleteBoneBits(Image3D *target, Image3D *seen, int slice, int maxComponentSize)
//this function calls conCompSize to get the size of a connected component, and floodFill to 
//delete small components. These functions work on a single slice, meaning the z just comes
//along. sorry.
{
	
	int i,j, compSize;
	GreyValue targetGrey = MAX_GREY_VALUE, replacementGrey = 0;

	for (j = 0; j < target->getYDim(); j ++)
	{
		for (i = 0; i < target->getXDim(); i ++)
		{
			//if (i == 373 && j == 260 && slice == 51)
			//{
			//	int qq = 0;
			//	RAWImageFile imageFile;
			//	imageFile.write("seendebug.raw3", *seen);
			//}

			
			if (target->getVoxelValue(i,j,slice) && !seen->getVoxelValue(i,j,slice))
			{
				compSize = conCompSize(target, seen, i,j,slice, maxComponentSize);
				if (compSize < maxComponentSize)
					floodFill(target, i,j,slice, targetGrey, replacementGrey);
			}
		}
	}
}

//A prettier version of the code to create a distance quantile image.
//the old p3dcontrol code is below and commented.
bool DQFTraining::createDQFImage(Image3D *image, const char * outputFilename, const M3DObject * object, int figureID)
{
	//if the object is given, then use it to determine the bounding box.
	if (object)
	{
		setDQFImageBoundingBox(image, object, figureID);
	}
	
	//Some Error checking.
	if (!image->imageCoordinatesInBounds(Xmin.getX(), Xmin.getY(), Xmin.getZ()) ||
		!image->imageCoordinatesInBounds(Xmax.getX(), Xmax.getY(), Xmax.getZ()))
	{
		cerr << "createDQFImage Error: the bounding box extends outside of the image." << endl;
		return false;
	}

	char * imageFilename = new char[13];

	if (outputFilename == NULL)  //make it from the image.
	{
		cerr << "createDQFImage Warning: no filename given, using default DQFIMAGE.raw3" << endl;

		strcpy(imageFilename, "DQFIMAGE.raw3");
	}


	RAWImageFile imageFile;

	//First, we'll declare and initialize the output images.
	Image3D * boneBoundaryImage = NULL;
	Image3D * gasBoundaryImage = NULL;
	Image3D * boneDistanceImage = NULL;
	Image3D * gasDistanceImage = NULL;

	GreyValue * imageVoxels = NULL;
	GreyValue * gbVoxels = NULL;
	
	boneBoundaryImage = new Image3D(image->getXDim(), image->getYDim(),
		image->getZDim());
	boneBoundaryImage->setSpacingAndOrigin(image->getXSpacing(),
		image->getYSpacing(), image->getZSpacing());
	boneBoundaryImage->clear();

	gasBoundaryImage = new Image3D(image->getXDim(), image->getYDim(),
		image->getZDim());
	gasBoundaryImage->setSpacingAndOrigin(image->getXSpacing(),
		image->getYSpacing(), image->getZSpacing());
	gasBoundaryImage->clear();

	boneDistanceImage = new Image3D(image->getXDim(), image->getYDim(),
		image->getZDim());
	boneDistanceImage->setSpacingAndOrigin(image->getXSpacing(),
		image->getYSpacing(), image->getZSpacing());
	boneDistanceImage->clear();

	gasDistanceImage = new Image3D(image->getXDim(), image->getYDim(),
		image->getZDim());
	gasDistanceImage->setSpacingAndOrigin(image->getXSpacing(),
		image->getYSpacing(), image->getZSpacing());
	gasDistanceImage->clear();

	//initalization is done.  


	int i, j, k, q, numVoxels;
	int xres = image->getXDim(), yres = image->getYDim(), zres = image->getZDim();
	int count1, count2;
	int numBoneVoxels;

	//Now, the boundary images are easiest.

	imageVoxels = image->getVoxels();
	numVoxels = xres*yres*zres;

	//First, bone boundary image.

	Image3D *blankImage = gasBoundaryImage;

	gbVoxels = blankImage->getVoxels();
	count1 = count2 = 0;

	for (i = 0; i < numVoxels; i ++)
	{
		if (image->mapDisplayToActual(imageVoxels[i]) > 1250)
		{
			gbVoxels[i] = MAX_GREY_VALUE;
			count1 ++;
		}
	}

	//Now we'll do a morphological open I think, with a 3x3 piece.  The idea is to make sure there aren't holes,
	//so that when we do a fill later (to determine the boundary voxels), the fill doesn't push into the bone, where 
	//it would run into inner boundaries we don't want.

	boneBoundaryImage->clear(MAX_GREY_VALUE);

	for (k = 0; k < zres; k ++)
	{
		for (j = 0; j < yres; j ++)
		{
			for (i = 0; i < xres; i ++)
			{
				//This serves as the kernel we're using for the morph.
				if (!blankImage->getVoxelValue(i,j,k) & !blankImage->getVoxelValue(i,j+1,k) & 
					!blankImage->getVoxelValue(i,j-1,k) & !blankImage->getVoxelValue(i+1,j,k) & 
					!blankImage->getVoxelValue(i-1,j,k))
				{
					boneBoundaryImage->setVoxel(i,j,k,0); boneBoundaryImage->setVoxel(i,j+1,k,0);
					boneBoundaryImage->setVoxel(i,j-1,k,0); boneBoundaryImage->setVoxel(i+1,j,k,0);
					boneBoundaryImage->setVoxel(i-1,j,k,0);
				}
			}
		}
	}

	//I'm commenting out the image writes, which were there for debug anyway.

	//if (imageFile.write("boneBoundaryImageBeforeFill.raw3", *boneBoundaryImage))
	//	cout << "Saved Image: boneBoundaryImageBeforeFill.raw3" << endl;
	//else
	//	cout << "Image failed to save: boneBoundaryImageBeforeFill.raw3" << endl;



	//A problem with the threshold image is that some bone bits are not bone, like prostate interior calcifications and
	//just random bits of noise.  I want to get rid of these.
	//Another step in bone processing. Here, I run some code that determines the size of a connected component of bone.  
	//If its size is too small (another hack parameter, like the open operator above), I'll delete it.
	//I might also use the training model. 
	//Again, its a per-slice operation

	blankImage->clear(); //used for deleteBoneBits.

	int componentSize = 50;  //get rid of bone bits smaller than 50 voxels.
	for (k = 0; k < zres; k ++)
		deleteBoneBits(boneBoundaryImage, blankImage, k, componentSize);

	//if (imageFile.write("boneBoundaryImageBeforeFillAfterBoneBits.raw3", *boneBoundaryImage))
	//	cout << "Saved Image: boneBoundaryImageBeforeFillAfterBoneBits.raw3" << endl;
	//else
	//	cout << "Image failed to save: boneBoundaryImageBeforeFill.raw3" << endl;



	//Now, boneBoundary ought to have fewer holes from the open and fewer non-bone bone bits.  
	//Now we can do a fill from the outside.  
	//I've never done a fill before, but it's straightforward I'm told...from Wiki, I'll
	//do the recursive flood fill...But it didn't run, probably because of limited stack size (maybe?).
	//so i made some modifications which makes it do more work seems to me, but finishes and does its job. 

	GreyValue targetGrey = 0, replacementGrey = 1;
	for (k = 0; k < zres; k ++)
		floodFill(boneBoundaryImage, 0,0,k, targetGrey, replacementGrey);
	//floodFill(boneBoundaryImage, 0,0,0, targetGrey, replacementGrey);

	//if (imageFile.write("boneBoundaryImageAfterFill.raw3", *boneBoundaryImage))
	//	cout << "Saved Image: boneBoundaryImageAfterFill.raw3" << endl;
	//else
	//	cout << "Image failed to save: boneBoundaryImageAfterFill.raw3" << endl;


	//Now I want to say if you're white (max value) and none of your neighbors is 1, then you're an interior
	//white voxel and should make yourself black.

	for (k = 0; k < zres; k ++)
	{
		for (j = 0; j < yres; j ++)
		{
			for (i = 0; i < xres; i ++)
			{
				//if all of the neighbors of (i,j,k) are on then set to zero.
				//This is pretty fast, since most voxels will quit early in the
				//if. In plane neighbors since the z spacing is so large.
				if (boneBoundaryImage->getVoxelValue(i,j,k) == MAX_GREY_VALUE && boneBoundaryImage->getVoxelValue(i,j+1,k) != 1 && 
					boneBoundaryImage->getVoxelValue(i,j-1,k) != 1 && boneBoundaryImage->getVoxelValue(i+1,j,k) != 1 && 
					boneBoundaryImage->getVoxelValue(i-1,j,k) != 1)//& boneBoundaryImage->getVoxelValue(i,j,k+1) & boneBoundaryImage->getVoxelValue(i,j,k-1))
				{
					boneBoundaryImage->setVoxel(i,j,k,0);
					//I used to set it to 2, but it's irrelevant, just extra info.
					count2 ++;
				}
			}
		}
	}

	//if (imageFile.write("boneBoundaryImageWith01MAX.raw3", *boneBoundaryImage))
	//	cout << "Saved Image: boneBoundaryImageWith01MAX.raw3" << endl;
	//else
	//	cout << "Image failed to save: boneBoundaryImage.raw3" << endl;

	//The boneBoundaryImage at this point is very interesting, giving us inside bone boundary, bone boundary, and outside 
	//bone label for every voxel.  Steve wants this used for negative versus positive distance at a voxel... 
	//So I need to keep it for the boneDistanceImage calculation... So copy it to blankImage, which I use as a
	//temp variable here.
	blankImage->clear();

	GreyValue *blnkVoxels = blankImage->getVoxels();

	gbVoxels = boneBoundaryImage->getVoxels();

	for (i = 0; i < numVoxels; i ++)
	{
		blnkVoxels[i] = gbVoxels[i];
	}

	//Okay, now turn boneBoundaryImage to just bone boundary and not bone boundary.
	for (i = 0; i < numVoxels; i ++)
	{
		if (gbVoxels[i] == 1 || gbVoxels[i] == 2)
		{
			gbVoxels[i] = 0;
		}
	}

	cout << count1 << " bone voxels, " << count1 - count2 << " bone boundary voxels. " << endl; 

	numBoneVoxels = count1 - count2;
	
	//if (imageFile.write("boneBoundaryImage.raw3", *boneBoundaryImage))
	//	cout << "Saved Image: boneBoundaryImage.raw3" << endl;
	//else
	//	cout << "Image failed to save: boneBoundaryImage.raw3" << endl;


	//Turns out it's easier to deal with all the bone stuff, then go after the gas.

	//I'm not completely satisfied with the boneBoundaryImage, since the labelled voxels are not 
	//only the bone exterior.  But we'll work with it anyway.
	//Not true, they are now with the open and fill having been implemented.  These ops are not applied to gas
	//though, so remember that.

	//We have to make a list of the positions of bone voxels, so as to speed the process. The
	//positions should be in world coordinates for speed, since we'll have to have the world
	//coords every time we want to compute a distance.

	std::vector<Vector3D *> boneCoords;
	Vector3D * vox;


	//Here I need to use the Xmin and Xmax for the bounding box.
	//z will go from Xmin.getZ() to Xmax.getZ().
	//y goes from Xmin.getY() to Xmax.getY().
	//x goes from Xmin.getX() to Xmax.getX().

	for (k = Xmin.getZ(); k <= Xmax.getZ(); k ++)
	{
		for (j = Xmin.getY(); j <= Xmax.getY(); j ++)
		{
			for (i = Xmin.getX(); i < Xmax.getX(); i ++)
			{
				if (boneBoundaryImage->getVoxelValue(i,j,k))
				{
					//Vector3D vox(i,j,k);
					vox = new Vector3D((double) i, (double) j, (double) k);
					image->imageToWorldCoordinates(*vox);
					boneCoords.push_back(vox);
				}
			}
		}
	}

	cout << "Bone voxels to look at : " << boneCoords.size() << endl;

	//Now for every voxel, we compute the distance-weighted average distance to the bone voxels.

	double sumd, tempd, sumweight;

	double maxdiffsofar = 0.0;


	//Time to finally make the distance image.
	vox = new Vector3D(0,0,0);
	Vector3D tempv(1,0,0);

	image->imageToWorldCoordinates(*vox);
	image->imageToWorldCoordinates(tempv);

	cout << "1 voxel distance is " << sqrt(vox->distSquare(tempv)) << ", should be .98mm?" << endl;


	//Now let's declare and populate a DQF image.
	//The constant 10 is sent as the size of the quantile functions...
	//Change that with a parameter of course.
	dqfImage = new DQFImage;
	
	dqfImage->setDims(Xmax.getX() - Xmin.getX() + 1, Xmax.getY() - Xmin.getY() + 1,
		Xmax.getZ() - Xmin.getZ() + 1, 10, image->getXDim(), image->getYDim(), image->getZDim());
	
	dqfImage->roi(Xmin, Xmax);

	//The 10.0 is hack.  
	dqfImage->setMapping(MAX_GREY_VALUE/imageMultiplier, imageDelta);


/*if (imageFile.write("BlankDQFImageBeforeSpacing.raw3", *dqfImage))
		cout << "Saved DQF image: " << "BlankDQFImageBeforeSpacing.raw3" << ".\n";
	else
	{
		cout << "DQF image failed to save: " << "BlankDQFImageBeforeSpacing.raw3" << ".\n";
		return false;
	}

if (imageFile.write("BlankDQFImageBeforeSpacing2.raw3", *dqfImage))
		cout << "Saved DQF image: " << "BlankDQFImageBeforeSpacing2.raw3" << ".\n";
	else
	{
		cout << "DQF image failed to save: " << "BlankDQFImageBeforeSpacing2.raw3" << ".\n";
		return false;
	}
	*/

	dqfImage->clear(0, false); 

	//if (imageFile.write("BlankDQFImageAfterClear.raw3", *dqfImage))
	//	cout << "Saved DQF image: " << "BlankDQFImageAfterClear.raw3" << ".\n";
	//else
	//{
	//	cout << "DQF image failed to save: " << "BlankDQFImageAfterClear.raw3" << ".\n";
	//	return false;
	//}

	Vector3D origin(image->getWorldOrigin());

	dqfImage->setSpacingAndOrigin(image->getXSpacing(), image->getYSpacing(), image->getZSpacing(), &origin);
	//This call doesn't do what you might expect (let you convert coords in the orignal image).  That is 
	//because some of the variables in the function are overwritten by DQFImage.  One must write this function
	//specifically for DQFImage, using the variables set aside for the embedding image.

/*	if (imageFile.write("BlankDQFImageAfterSpacing.raw3", *dqfImage))
		cout << "Saved DQF image: " << "BlankDQFImageAfterSpacing.raw3" << ".\n";
	else
	{
		cout << "DQF image failed to save: " << "BlankDQFImageAfterSpacing.raw3" << ".\n";
		return false;
	}
*/

	//Now we have specified the dqfImage, as in what the size is and the region of interest with 
	//respect to the orginal greyscale.  Now we can start setting voxels as quantile functions...

	DQuantileFunction *qf = new DQuantileFunction(dqfImage->getWidth());
	qf->useDefaultWeighting(weightingMu, weightingVar);
	qf->clearData();


	maxdiffsofar = 0.0;

	double recmu = 0.0;
	double madmu = 0.0;
	
	double * madqf;//this was causing memory leak = new double[dqfImage->getWidth()];
	double * recqf = new double[dqfImage->getWidth()];

	for (q = 0; q < dqfImage->getWidth(); q ++)
	{
		recqf[q] = 0.0;
	}
				
	for (k = Xmin.getZ(); k <= Xmax.getZ(); k ++)
	//for (k = 43; k <= 72; k ++)
	{
		cerr << ".";
	
		for (j = Xmin.getY(); j <= Xmax.getY(); j ++)
		//for (j = 222; j <= 314; j ++)
		{
			for (i = Xmin.getX(); i <= Xmax.getX(); i ++)
			//for (i = 213; i <= 316; i ++)
			{
				vox->set((double) i, (double) j, (double) k);
				//vox is the current voxel to consider.
				image->imageToWorldCoordinates(*vox);

				qf->clearData();
		
				sumd = 0.0;
				sumweight = 0.0;
				for (q = 0; q < boneCoords.size(); q ++) //for every bone voxel.
				{
					tempd = sqrt(vox->distSquare(*boneCoords[q]));

					qf->addData(tempd);

					/*if (tempd <= 1)
					{						
						sumd += tempd;
						sumweight += 1.0;
					}
					else
					{
						tempi = exp(-.5*(tempd - 1)*(tempd - 1));
						sumd += tempi*tempd;
						sumweight += tempi;
					}*/
					//sumd += vox->distSquare(*boneCoords[q]);
					//debug:
					//Vector3D v = *boneCoords[q];
					//tempi = vox->distSquare(*boneCoords[q]);
					//if (tempi < 1) 
					//	sumd = boneCoords.size();
				}
				
				//sumd /= sumweight;
				
				//Now use the inside outside boundary image stored in blankImage to negate some voxels.
				//if (blankImage->getVoxelValue(i,j,k) == 0) //then on bone inside, 
				//	sumd *= -1.0;
				
				
				//tempi = (double) MAX_GREY_VALUE * sumd / image->maxExtent();
				//tempi = (double) MAX_GREY_VALUE * sumd / (image->maxExtent() * image->maxExtent());
				//With neg and pos values now, let's plus minus the middle greyvalue.
				
				//if (sumd < 0)
				//	tempi = 0.5 * MAX_GREY_VALUE + MAX_GREY_VALUE * 5*sumd/image->maxExtent();
				//else
				//	tempi = 0.5 * MAX_GREY_VALUE + MAX_GREY_VALUE * sumd/image->maxExtent();
				
				//boneDistanceImage->setVoxel(i,j,k,(GreyValue) tempi);
				//boneDistanceImage->setVoxel(i,j,k,(GreyValue) 500);

				dqfImage->setVoxel(i-Xmin.getX(),j-Xmin.getY(),k-Xmin.getZ(), qf->getQF());
				//dqfImage->setVoxel(i-Xmin.getX(),j-Xmin.getY(),k-Xmin.getZ(), recqf);

				madqf = qf->getQF();
				dqfImage->getImageVoxelValue(i,j,k, recqf);

				recmu = 0.0;
				madmu = 0.0;

				//compute the means
				for (q = 0; q < dqfImage->getWidth(); q ++)
				{
					recmu += recqf[q];
					madmu += madqf[q];
				}

				recmu /= dqfImage->getWidth();
				madmu /= dqfImage->getWidth();

				if (fabs(madmu - recmu) > maxdiffsofar)
				{
					maxdiffsofar = fabs(madmu - recmu);
					//if (maxdiffsofar > .1)
					//	int qq = 0;
					cout << "Max diff ideal/image changed: now " << maxdiffsofar << " " << qf->getMean() << ".\n";
				}

				if (maxdiffsofar > qf->getMean())
				{
					cout << "\tdiff too big. voxel (" << i << ", " << j << ", " << k << ")\n";
				}
	
			}
		}
	}

	//if (imageFile.write("boneDistanceImage.raw3", *boneDistanceImage))
	//	cout << "Saved Image: boneDistanceImage.raw3" << endl;
	//else
	//	cout << "Image failed to save: boneDistanceImage.raw3" << endl;

	

	if (imageFile.write(outputFilename, *dqfImage))
		cout << "Saved DQF image: " << outputFilename << ".\n";
	else
	{
		cout << "DQF image failed to save: " << outputFilename << ".\n";
		return false;
	}

	//I need to debug the dqf image.  So, I'm going to save the first element into one of the other images so 
	//I can see if they line up.

/*	boneDistanceImage->clear();

	for (k = Xmin.getZ(); k <= Xmax.getZ(); k ++)
	{
		for (j = Xmin.getY(); j <= Xmax.getY(); j ++)
		{
			for (i = Xmin.getX(); i <= Xmax.getX(); i ++)
			{

				dqfImage->getImageVoxelValue(i,j,k, recqf);
				//Now write recqf[0] to the image, but need to turn it into a greyvalue

				//tempi = MAX_GREY_VALUE * (1.0 - recqf[0]/20.0);
				tempi = recqf[0]*MAX_GREY_VALUE/20.0;
				boneDistanceImage->setVoxel(i,j,k, (GreyValue) tempi);
			}
		}
	}	
	if (imageFile.write("debugDQFelement0.raw3", *boneDistanceImage))
		cout << "Saved Image: debugDQFelement1.raw3" << endl;
	else
		cout << "Image failed to save: debugDQFelement1.raw3" << endl;


	boneDistanceImage->clear();

	for (k = Xmin.getZ(); k <= Xmax.getZ(); k ++)
	{
		for (j = Xmin.getY(); j <= Xmax.getY(); j ++)
		{
			for (i = Xmin.getX(); i <= Xmax.getX(); i ++)
			{

				dqfImage->getImageVoxelValue(i,j,k, recqf);
				//Now write recqf[0] to the image, but need to turn it into a greyvalue

				//tempi = MAX_GREY_VALUE * (1.0 - recqf[4]/20.0);
				tempi = recqf[4]*MAX_GREY_VALUE/20.0;
				boneDistanceImage->setVoxel(i,j,k, (GreyValue) tempi);
			}
		}
	}	
	if (imageFile.write("debugDQFelement4.raw3", *boneDistanceImage))
		cout << "Saved Image: debugDQFelement4.raw3" << endl;
	else
		cout << "Image failed to save: debugDQFelement4.raw3" << endl;


	boneDistanceImage->clear();

	for (k = Xmin.getZ(); k <= Xmax.getZ(); k ++)
	{
		for (j = Xmin.getY(); j <= Xmax.getY(); j ++)
		{
			for (i = Xmin.getX(); i <= Xmax.getX(); i ++)
			{

				dqfImage->getImageVoxelValue(i,j,k, recqf);
				//Now write recqf[0] to the image, but need to turn it into a greyvalue

				//tempi = MAX_GREY_VALUE * (1.0 - recqf[9]/20.0);
				tempi = recqf[9]*MAX_GREY_VALUE/20.0;
				boneDistanceImage->setVoxel(i,j,k, (GreyValue) tempi);
			}
		}
	}	
	if (imageFile.write("debugDQFelement9.raw3", *boneDistanceImage))
		cout << "Saved Image: debugDQFelement9.raw3" << endl;
	else
		cout << "Image failed to save: debugDQFelement9.raw3" << endl;

*/
	return true;
}


//This function, which records a dqf file for a model and dqf image pair, has lots of comments related to the 
//fact that it used to be a part of hacky code to test the dqf stuff.  Also, the file it used to write contained
//all distances rather than the quantile function, and matlab was responsible for the processing.
bool DQFTraining::doModelDQFTraining(Image3D *image, const M3DObject * object, int figureID, const char *saveFilename, const char *dqfFilename)
{

	RAWImageFile imageFile;

	if (dqfFilename != NULL)
	{
		dqfImage = (DQFImage *) imageFile.read(dqfFilename);
	
		if (dqfImage == NULL)
		{
			cout << "Error doModelDQFTraining: DQF image failed to load.\n";
			return false;
		}
	}

	if (saveFilename == NULL)
	{
		saveFilename = new char[13];
		strcpy((char *) saveFilename, "MODEL_DQF.dqf");
	}

	//DQuantileFunction *qf = new DQuantileFunction(dqfImage->getWidth());
	//qf->useDefaultWeighting();
	//qf->clearData();



	//More stuff?  I guess.  I am going to have to make this stuff prettier soon, before it
	//gets out of hand.

	//But first, Let's make an image of actual DQFs, rather than the single bin distance-weighted
	//average that is the above boneDistanceImage.  I can't use the image3d class for this, since
	//every voxel will have, say an array of doubles (the quantile function on distance at that voxel.

	//But just to start, Eli thought an interesting view would be to record the dqfs for all the boundary
	//points on the object of interest.  Then we could visualize some interesting things.  So, in order
	//not to generate and image of dqfs for this, I'll just say for every boundary point...

	//And in order not to write out a file that's just too huge, I'll record only those within 5 cm maybe...
	
	//okay, get the points.


	M3DFigure * figure = object->getFigurePtr(figureID);
    
	Xferlist * xferList = convertM3DtoXfer(figure);

	ThallCode::Pointlist_server2 * pList = new ThallCode::Pointlist_server2;
	pList->init(xferList);
	pList -> ComputeSubdivPointCloud(modelSubDivLevel);

	//const char * filename = script.getStringValue("Model", NULL);
	//char * myfilename = new char[1000];

	//sprintf(myfilename, "%s", filename);

	//if (outputFilename)
	//	sprintf(myfilename, "%s", outputFilename);
	//else
	//	sprintf(myfilename, "%s", imageFilename);

	//myfilename[strlen(myfilename) - 4] = '\0';
	//sprintf(myfilename, "%s.dqf", myfilename);


	//FILE *fout = fopen(myfilename, "w");
	FILE *fout = fopen(saveFilename, "w");

	int numPoints;//, numDists;
	Bpoint * Bpoints;
    pList->subdivboundaryinfo(&numPoints, &Bpoints);

	//double *distances = new double[boneCoords.size()];
	
	double *recqf = new double[dqfImage->getWidth()];


	fprintf(fout, "%i %i ", numPoints, dqfImage->getWidth());

	

 	for (int i = 0; i < numPoints; i ++)
	{
		Vector3D Pi(Bpoints[i].pnt); //the figure boundary point.
		image->modelToImageCoordinates(Pi);

		//cout << "Pi printing:\n";
		//Pi.print();

		//Vector3D Pii(Bpoints[i].pnt);
		//image->modelToImageCoordinates(Pii);

		//cout << "Now Pii printing:\n";
		//Pii.print();

		//Vector3D difff(Pi - Pii);
		//difff.print();

		//now loop over all bone voxels and record the distances.

		//memset(distances, 0, boneCoords.size()*sizeof(double));
		//numDists = 0;

		//fill distances out, then write it in the file...

		//for (q = 0; q < boneCoords.size(); q ++) //for every bone voxel.
		//{
		//	tempd = sqrt(Pi.distSquare(*boneCoords[q]));
			//only record those within 4 cm. or 8 or something
		//	if (tempd <= 8)
		//	{
		//		distances[numDists] = tempd;
		//		numDists ++;
		//	}
		//}

		//qf->clearData();

		//now can write out the distances
		//fprintf(fout, "%i ", numDists);

		dqfImage->getInterpolatedImageVoxelValue(Pi.getX(), Pi.getY(), Pi.getZ(), recqf);

		for (int j = 0; j < dqfImage->getWidth(); j ++)
		{
			fprintf(fout, "%lf ", recqf[j]);
			//qf->addData(distances[j]);
		}

		//cout << "Quantile Function says mean is " << qf->getMean() << ".\n";

		//image->worldToImageCoordinates(Pi);

		//Using the interpolated voxel value reduces the error from .06 out of 4 to
		//.002 out of 4.  But it does involve a lot more accesses.  We want to be sure
		//to train and segment the same way. Let's stick with interpolated.
		//dqfImage->getInterpolatedImageVoxelValue(Pi.getX(), Pi.getY(), Pi.getZ(), recqf);
		//Pi.round();
		//dqfImage->getImageVoxelValue(Pi.getX(), Pi.getY(), Pi.getZ(), recqf);
		//sumd = 0.0;
		//for (j = 0; j < dqfImage->getWidth(); j ++)
		//	sumd += recqf[j];

		//sumd /= dqfImage->getWidth();
		//if (fabs(qf->getMean() - sumd) > maxdiffsofar)
		//{
		//	maxdiffsofar = fabs(qf->getMean() - sumd);
		//	cout << "Max diff model/image increased: now " << maxdiffsofar << " " << qf->getMean() << ".\n";
		//}

		//cout << "DQF image says mean is " << sumd << ".\n";

	}

	fclose(fout);

	delete [] (xferList->atomlist);
	delete xferList;
	delete pList; 
	
	//There, now use matlab to process all that raw info.

	return true;
}

















/********************************************************************************************************************************************************/



/********************************************************************************************************************************************************
Code for the DQFMatch class, which will compute the dqf match at target time.
*********************************************************************************************************************************************************/




/********************************************************************************************************************************************************/




DQFpca::DQFpca()
{
	neigs = 0;
	dim = 0;
	mu = NULL;
	var = NULL;
	eigs = NULL;
}

DQFpca::~DQFpca()
{
	delete []mu;
	delete []var;
	for (int i = 0; i < neigs; i ++)
		delete []eigs[i];
}

void DQFpca::readDQFpca(FILE * fin, int uwidth)
{
	//Here's how I imagine the format for a pca.  it's dim x (neigs + 2), just
	//like in matlab.  First dim elements is mu, next dim is var.  The number
	//of non-zero entries of var is neigs + 1, since the last non-zero is the 
	//residual.  That follows with dim elements, neigs times. And that's it.
	
	dim = uwidth;
	
	//dim is the same for all the points for a figure, so that shouldn't have to 
	//be read.  Knowing width I can just start.
	int i, j;

	mu = new double[dim];
	var = new double[dim];
	//don't know how big eigs is yet, have to read how many non-zeros are in vars.

	for (i = 0; i < dim; i ++)
		fscanf(fin, "%lf ", &mu[i]);

	neigs = 0;
	for (i = 0; i < dim; i ++)
	{
		fscanf(fin, "%lf ", &var[i]);

		if (var[i] > 0)
			neigs ++;
	}

	//Now, neigs - 1 is how many eigenmodes are given; the last non-zero is the 
	//residual error.

	neigs --;
	
	eigs = new double *[neigs];

	for (j = 0; j < neigs; j ++)
	{
		eigs[j] = new double[dim];

		for (i = 0; i < dim; i ++)
		{
			fscanf(fin, "%lf ", &eigs[j][i]);
		}
	}

	//There.  Now this particular DQFPca is populated.
		
}

double DQFpca::computeDQFMahal(double *qf)
{
	//this function uses the pca info to compute a mahalanobis distance
	//of the arg qf.  The normalization must happen elsewhere.
	//That's not true.  I'll do it here, since I know neigs.

	//First compute difference wrt mu;
	double *diff, *coeffs;
	int i, j;

	diff = new double[dim]; //the difference vector;
	coeffs = new double[neigs];
	//the coefficients of the projection.

	for (i = 0; i < dim; i++)
	{
		diff[i] = qf[i] - mu[i];
	}

	//Now compute the projections onto the eigenmodes.
	double distInPCASubspace = 0.0;

	for (j = 0; j < neigs; j ++)
	{
		coeffs[j] = 0.0;

		//Simple dot product.
		for (i = 0; i < dim; i ++)
		{
			coeffs[j] += diff[i]*eigs[j][i];
		}
		coeffs[j] *= coeffs[j];

		distInPCASubspace += coeffs[j];
		coeffs[j] /= var[j];
	}

	//Now compute the actual distance to the mean, and compare with the distInPCASubspace,
	//for the residual.  The match returned will be the sum of coeffs and this additional
	//term.

	double distInSpace = 0.0, mahal = 0.0;
	
	for (i = 0; i < dim; i ++)
	{
		distInSpace += diff[i]*diff[i];
	}

	for (j = 0; j < neigs; j ++)
		mahal += coeffs[j];

	mahal += (distInSpace - distInPCASubspace)/var[neigs];
	//Because distInPCASubspace is always less than distInSpace, don't have to check 
	//the sign.

	//Now, finally, make the result MU zero and unit STD/VAR, so to speak...if everything's 
	//Gaussian...It has to do with the fact that the Mahal distance of a Gaussian variable 
	//follows a Chi-square distribution with degree of freedom equal to the dimension of the
	//Gaussian. so expected value neigs + 1, with var = 2*(neigs + 1).

	mahal -= neigs + 1;
	mahal /= sqrt(2.0*(neigs + 1));

	return mahal;
}

int DQFpca::getNEigs()
{
	return neigs;
}

DQFMatch::DQFMatch()
{
	surfaceLevel = 0;
	numFigs = 0;
	dqfPca = NULL;
	dqfImage = NULL;
	numPoints = NULL;
	mu = 2.0;
	var = 1.0;
}

DQFMatch::~DQFMatch()
{
}

DQFMatch::DQFMatch(Registry & reg)
{
	surfaceLevel = 0;
	numFigs = 0;
	dqfPca = NULL;
	dqfImage = NULL;
	numPoints = NULL;
	mu = 2.0;
	var = 1.0;
	loadConfig(reg);
}

void DQFMatch::loadConfig(Registry & reg)
{
	surfaceLevel   = reg.getIntValue("dqfModelSurfaceLevel", 0);
}

void DQFMatch::initializeDQFMatch(Image3D *uimage, const M3DObject * object, const char ** dqfImageFilenames, const char ** dqfStatFilenames)
{
	image = uimage;

	RAWImageFile imageFile;



	numFigs = object->getFigureCount();


	numPoints = new int[numFigs];
	width = new int[numFigs];

	dqfPca = new DQFpca *[numFigs];
	dqfImage = new DQFImage *[numFigs];

	for (int f = 0; f < numFigs; f ++)
	{
		dqfImage[f] = NULL;
		dqfPca[f] = NULL;
		//const char * filename = P3DControl::calcDQFFilename((char *) imageFilename, f);

		dqfImage[f] = (DQFImage *) imageFile.read(dqfImageFilenames[f]);

		loadDQFpcasForFigure(dqfStatFilenames[f], f);
	}
}

//I can't write this until I know what the matlab ouput format is.  Hmmm....
//Now I'll go write that out.  See readDQFpca for more comments...
void DQFMatch::loadDQFpcasForFigure(const char * dqfStatFilename, int figureID)
{
	if (dqfPca[figureID] != NULL)
		delete [] dqfPca[figureID];


	FILE * fin = fopen(dqfStatFilename, "r");

	if (fin == NULL) //The stat file was not defined or did not exist.
 		return;

	fscanf(fin, "%i %i ", &numPoints[figureID], &width[figureID]);

	dqfPca[figureID] = new DQFpca[numPoints[figureID]];

	for (int i = 0; i < numPoints[figureID]; i ++)
	{
		dqfPca[figureID][i].readDQFpca(fin, width[figureID]);
	}

	//Now we've read the pcas for the points of this figure.  
}

double DQFMatch::computeDQFMatch(Image3D * image, const M3DObject * targetObject, int figureID)
{
	if (figureID < 0 || dqfImage[figureID] == NULL || dqfPca[figureID] == NULL) 
	{
		return 0.0;  
	}
	

	//The steps are to
/*	For each point on the surface of the figure.
		get the qf for that point.
		compute its match and add it into the running total.

    return the total.

	I'm going to do the same normalization that Eli does, that depends on the total variances


  */

	M3DFigure * figure = targetObject->getFigurePtr(figureID);
    
	Xferlist * xferList = convertM3DtoXfer(figure);

	ThallCode::Pointlist_server2 * pList = new ThallCode::Pointlist_server2;
	pList->init(xferList);
	pList -> ComputeSubdivPointCloud(surfaceLevel);
	//The suddivlevel should come from the dqfPca


	int numPoints;//, numDists;
	Bpoint * Bpoints;
    pList->subdivboundaryinfo(&numPoints, &Bpoints);

	double *recqf = new double[dqfImage[figureID]->getWidth()];


	double match = 0.0;

 	for (int i = 0; i < numPoints; i ++)
	{
		Vector3D Pi(Bpoints[i].pnt); //the figure boundary point.
		image->modelToImageCoordinates(Pi);

		dqfImage[figureID]->getInterpolatedImageVoxelValue(Pi.getX(), Pi.getY(), Pi.getZ(), recqf);

		//Now compute the match of recqf with respect to dqfPca[figureID][i].

		match += dqfPca[figureID][i].computeDQFMahal(recqf);
	}

	match /= sqrt((double) numPoints);

	delete [] (xferList->atomlist);
	delete xferList;
	delete pList; 


	return match;
}


