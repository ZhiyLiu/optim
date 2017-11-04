//Image Quantile Function image match stuff is here and in the .h.

// MSVC6: disable compiler warning about debugging symbol names too long
#pragma warning (disable:4786)

#include <iostream>
#include <fstream>
#include <algorithm>
#include "ObjectRelativeSampling.h"
#include "SurfacePatchEnsemble.h"
#include "P3DControl.h"
#include "assert.h"
#include <float.h>

#include <typeinfo>

/********************************************************************************************************************************************************/



/********************************************************************************************************************************************************
Code for the QuantileFunction class.
*********************************************************************************************************************************************************/
/*
-po -cs config-mPatch-training.txt -model 3109.fr04.bp.ma.m3d -image 3109.fr04.raw3 -mpConfig config-mPatch-training-GlobalInLocalOut2Regions.txt -mpOut 3109.fr04.bp.try-GlobalInt2Regions.mpf
-po -cs config-mPatch-training.txt -model 001-rectum.normalized.m3d -image 3106.fr01.raw3 -mpConfig config-mPatch-training-RectumGlobalInLocalOut2Regions.txt -mpOut 001.try-GlobalInt.mpf 

-po -cs config-mPatch-training.txt -model 001-rectum.normalized.m3d -image 3106.fr01.raw3 -mpConfig config-mPatch-training-RectumGlobalInLocalOut2Regions.txt -mpOut 001.try-GlobalInt.mpf 


-lv -po -ci -cs r:\proj\hn1\scripts\config-grey-hist-figure-local-MPF.txt -image r:\home\saboo\public\rectums\3106\data\images\gray\017.raw3 -outTile r:\home\saboo\public\rectums\3106\data\segmentations\mpfDebug-s2-align3\rectum\figure\loo-course-mpf\017-rectum.byu -tileUnit 1 -model r:\home\saboo\public\rectums\3106\data\segmentations\mpfDebug-s2-align3\rectum\contourInit\loo-course\final\noXform\017-rectum.m3d -PGAModel r:\home\saboo\public\rectums\3106\data\shapeStat\align3\rectum\loo-courseResidue\rectum.017SymMean.m3d -mpConfig r:\home\saboo\public\rectums\3106\data\imageStat\mpfDebug-s2\rectum\loo-course-mpf\mpfStatsDebug2Regions\017.try-GlobalInt-trainedMPF.mpf -outModel 017-rectum.m3d

-po -cs r:\proj\hn1\scripts\config-grey-hist-figure-local-MPF.txt -image  r:\home\saboo\public\rectums\3108\data\images\gray\001.raw3 -model  r:\home\saboo\public\rectums\3108\data\shapeStat\align3\rectum\loo-course\M3dNormalized\001-rectum.normalized.m3d -PGAModel  r:\home\saboo\public\rectums\3108\data\shapeStat\align3\rectum\loo-courseResidue\rectum.001SymMean.m3d -mpConfig  r:\home\saboo\public\rectums\3108\data\imageStat\mpfDebug-s2\rectum\loo-course-mpf\mpfStatsDebugLocalLocal\001.try-LocalLocal-trainedMPF.mpf -outModel tmp3108-001fromCorrect.m3d -figureIterations 1

-lv -po -ci -cs r:\proj\hn1\scripts\config-grey-hist-figure-local-MPF-ignoreEnds.txt -image r:\home\saboo\public\rectums\3108\data\images\gray\001.raw3 -outTile myRectum001.byu -tileUnit 1 -model r:\home\saboo\public\rectums\3108\data\segmentations\mpfDebug-s2-align3\rectum\contourInit\loo-course\final\noXform\001-rectum.m3d -PGAModel r:\home\saboo\public\rectums\3108\data\shapeStat\align3\rectum\loo-courseResidue\rectum.001SymMean.m3d -mpConfig r:\home\saboo\public\rectums\3108\data\imageStat\mpfDebug-s2\rectum\loo-course-mpf\mpfStatsDebug2Regions\001.try-2Regions-trainedMPF.mpf -outModel myRectum001.m3d

  -lv -po -ci -cs r:\proj\hn1\scripts\config-grey-hist-figure-local-MPF-ignoreEnds.txt -image r:\home\saboo\public\rectums\3108\data\images\gray\001.raw3 -outTile myRectum001.byu -tileUnit 1 -model r:\home\saboo\public\rectums\3108\data\segmentations\mpfDebug-s2-align3\rectum\contourInit\loo-course\final\noXform\001-rectum.m3d -PGAModel r:\home\saboo\public\rectums\3108\data\shapeStat\align3\rectum\loo-courseResidue\rectum.001SymMean.m3d -mpConfig 001.try-2Regions-trainedMPF.mpf -outModel myRectum001.m3d

  -lv -po -ci -cs r:\proj\prostate9\malePelvis1\tests\run21\pipeScripts\config\config-grey-hist-figure-local-MPF-augAtoms.txt  -transform r:\proj\prostate9\malePelvis1\pats\3106\adapt\xform\run400-bone\3106.fr17.simTrans.m3d -image r:\proj\prostate9\malePelvis1\pats\3106\adapt\image\run01-gray-origRes\3106.fr17.raw3 -model r:\proj\prostate9\malePelvis1\pats\3106\adapt\shapeStat\run313-bone-aligned\3106.fr17.bp.aa-SymMean.m3d -mpConfig r:\proj\prostate9\malePelvis1\pats\3106\adapt\imageStat\run575-BPR\mpfTrained-local0-gasBone\3106.fr17-trainedMPF.mpf -outModel 3106.fr17.bpr.OldProstate-noAugment-OldFunc.m3d

-lv -ci -cs h:\levy\mpf\config\config-mpf-training.txt -mpConfig h:\levy\mpf\config\config-mpf-training-local0.txt -po -image h:\levy\mpf\BWH_PNL_01\intRescaleEMS\BWH_PNL_01_MRI_RAI_regAtlas-adj.gipl -mpBinaryImage0 h:\levy\mpf\BWH_PNL_01\intRescaleEMS\BWH_PNL_01_seg_Lcaudate_RAI_regAtlas.gipl -figureId 0 -model  h:\levy\mpf\BWH_PNL_01\intRescaleEMS\BWH_PNL_01_roundTwoLeft.m3d -mpOut BWH_PNL_01_riqf.mpf



*/
/********************************************************************************************************************************************************/




/************************************************************************************************************


Defining some utilities.  I have a mixQF and linspace and interp1.  I don't know how I'll want them, so 
I put them here at the top.  The idea is that consider two regions and their associated qfs.  The qf of 
the pooled region is something.  If the regions have no voxels in common, then my mixQF should give the
right result.  But since these regions may partially overlap, the real pooled region qf and the mixed 
qf might be different. This is important since in my region growing, at training (matlab) I can't pool regions
and at target, I may want to combine and uncombine regions a lot.

  *****************************************************************************************/

void mixQF(const double *inqf0, const double *inqf1, int qf0Size, int qf1Size, double *qfInt, int qfIntSize, double alpha);
double * linspace(double xMin, double xMax, int numSteps);
double * interp1(const double *x, const double *Y, int size, double *xi, int outSize, double extrapval);

void mixQF(const double *inqf0, const double *inqf1, int qf0Size, int qf1Size, 
		   double *qfInt, int qfIntSize, double alpha)
{
	//invert, interpolate, invert;  assume qfs here.

	
	double *xOut = linspace((1.0/(2*qfIntSize)), 1.0-(1.0/(2*qfIntSize)), qfIntSize);
//pha = .001;
	//This was me testing out for really small alphas.
	//We have the bin values, which are for .005, .015, etc.  but we should also want
	//0 and 1.  So we extrapolate a little.
	//It didn't work like my matlab changes did.  
	
	double *x010 = linspace(-(1.0/(2*qf0Size)), 1.0+(1.0/(2*qf0Size)), qf0Size+2);
	double *x011 = linspace(-(1.0/(2*qf1Size)), 1.0+(1.0/(2*qf1Size)), qf1Size+2);
	qf0Size += 2; 
	qf1Size += 2;
	x010[0] = 0.0; x010[qf0Size-1] = 1.0;
	x011[0] = 0.0; x011[qf1Size-1] = 1.0;
	double sm = 1.0e-10;

	double *qf0 = new double[qf0Size];
	memcpy(qf0+1, inqf0, (qf0Size-2)*sizeof(double));
	qf0[0] = qf0[1]-.5*(qf0[2]-qf0[1])-sm;
	qf0[qf0Size-1] = qf0[qf0Size-2]+.5*(qf0[qf0Size-2]-qf0[qf0Size-3])+sm;

	double *qf1 = new double[qf1Size];
	memcpy(qf1+1, inqf1, (qf1Size-2)*sizeof(double));
	qf1[0] = qf1[1]-.5*(qf1[2]-qf1[1])-sm;
	qf1[qf1Size-1] = qf1[qf1Size-2]+.5*(qf1[qf1Size-2]-qf1[qf1Size-3])+sm;

	
	//We want the made up values for 0 and 1 included in the xi, to allow some
	//interpolation past the defined values.  It's helpful for neagtive alpha.
	double xMin = MIN(qf0[0], qf1[0]);
	double xMax = MAX(qf0[qf0Size-1], qf1[qf1Size-1]);

	int superSize = (int) ceil(5*(xMax - xMin)); //5*qfIntSize; //for more accuracy.
	double *xi = linspace(xMin, xMax, superSize);


	//Let's quit in the case of one or zero, for that first combination.
	if (alpha < sm && alpha > 0)
	{
		double * myqfInt = interp1(x010, qf0, qf0Size, xOut, qfIntSize, 0);
		memcpy(qfInt, myqfInt, qfIntSize*sizeof(double));
		delete [] myqfInt;
		delete []xi;
		delete []qf0;
		delete []qf1;
		delete []x010;
		delete []x011;
		delete []xOut;
		return;
	}
	else if (alpha > 1.0 - sm)
	{
		double * myqfInt = interp1(x011, qf1, qf1Size, xOut, qfIntSize, 0);
		memcpy(qfInt, myqfInt, qfIntSize*sizeof(double));
		delete [] myqfInt;
		delete []xi;
		delete []qf0;
		delete []qf1;
		delete []x010;
		delete []x011;
		delete []xOut;
		return;
	}

	//Remember, we're allowing for negative alphas.

	double *cdf0 = interp1(qf0, x010, qf0Size, xi, superSize, 0);
	double *cdf1 = interp1(qf1, x011, qf1Size, xi, superSize, 0);


	int i;
	//Check for zeros after ones.
	for (i = 1; i < superSize; i ++)
	{
		if (cdf0[i] == 0 && cdf0[i] < cdf0[i-1])
			cdf0[i] = 1.0;

		if (cdf1[i] == 0 && cdf1[i] < cdf1[i-1])
			cdf1[i] = 1.0;
	}


	//Now we can interpolate the cdf's.

	for (i = 0; i < superSize; i ++)
	{
		cdf0[i] = MIN(cdf0[i], 1.0);
		cdf1[i] = MIN(cdf1[i], 1.0);

		cdf0[i] = (1-alpha)*cdf0[i] + alpha*cdf1[i];

		cdf0[i] = MAX(cdf0[i], 0.0);
		cdf0[i] = MIN(cdf0[i], 1.0);
	}

	//Unlike the zeros, where we want the last zero, we want the first 1.
	//but interp will find the last, so we actually just make the cdf smaller.

	
	i = 0;
	while(i < superSize && cdf0[i] < 1.0-sm)
		i ++;
	
	
	int newSuperSize = i + 1;
	if (i == superSize)
		newSuperSize = superSize;


	//now invert.

	double *temp = interp1(cdf0, xi, newSuperSize, xOut, qfIntSize, 0);

	memcpy(qfInt, temp, qfIntSize*sizeof(double));

	
	//Now one more thing. The above method is simple interp, when the qf is
	//supposed to be the average of some amount of data. So here, I'll go 
	//through the data one more time, take the average of values within 
	//some domain.
/*
	double *xOut2 = linspace(0.0, 1.0, qfIntSize+1);

	i = 0;
	double summ;
	int cnt, cur = 0;
	while(i < newSuperSize)
	{
		summ = 0.0;
		cnt = 0;
		while(cdf0[i] >= xOut2[cur] && cdf0[i] < xOut2[cur+1] && i < newSuperSize)
		{
			summ += xi[i];
			cnt ++;
			i ++;
		}

		if (cnt)
			qfInt[cur] = summ/cnt;

		cur ++;
	}


	delete []xOut2;

*/


	//deletes.
	delete []xi;
	delete []qf0;
	delete []qf1;
	delete []x010;
	delete []x011;
	delete []xOut;
	delete []cdf0;
	delete []cdf1;
	delete []temp;

}

double * interp1(const double *x, const double *Y, int size, double *xi, int outSize, double extrapval)
{
	//Fill in yi at the xi given the current x and Y.
	//Basically, for every xi, find the before and after in x, and interpolate the associated 
	//Y appropriately.

	int i, before;

	double *yi = new double[outSize];

	for( i = 0; i < outSize; i ++ )
	{
		//Find xi[i] in x.
		if (xi[i] < x[0] || xi[i] > x[size-1])
		{
			yi[i] = extrapval;
			continue;
		}
		
		before = 0;

		//I changed the below < x[before+1] from <= because I wanted to go to the last
		//instance.  in matlab, I unique the vectors before sending to interp, but here
		//no need with this change.
		while (!(xi[i] >= x[before] && xi[i] < x[before+1]) && (before < size-1))
		{
			before ++;
		}

		//This is a singleton case.  if xi was greater than the greatest x, we would 
		//have already quit.  The above while stops and xi must be equal or less than
		//the greatest x.  This case checks equal basically.
		if (before == size-1)
		{
			yi[i] = Y[size-1];
			continue;
		}

		if (x[before] == x[before + 1])
		{
			while (x[before] == x[before + 1])
				before ++;

			yi[i] = Y[before];
			continue;
		}

		//now before is the left side.
		//simple linear interp.
		double alpha = (xi[i] - x[before]) / (fabs(x[before] - x[before + 1]));

		//so, simply
		yi[i] = (1.0-alpha)*Y[before] + alpha*Y[before+1];

	}

	return yi;

}

double * linspace(double xMin, double xMax, int numSteps)
{
	double * x = new double[numSteps];
	double stepSize = (xMax - xMin) / (numSteps - 1);


	for (int i = 0; i < numSteps; i ++)
	{
		x[i] = xMin + i*stepSize;
	}

	return x;
}





/********************************************************************************************************************************************************
Code for the QuantileFunction class.
*********************************************************************************************************************************************************/



QuantileFunction::QuantileFunction()
{
	rawData.clear();
	qf = NULL;
	qfSize = 0;
	setQFSize(10);

	useMyWeighting = false;
	dirty = true;

	lowWeight = .0;
	highWeight = .0;
	totWeight = .0;
}

QuantileFunction::QuantileFunction(int size)
{
	rawData.clear();
	qf = NULL;

	setQFSize(size);

	useMyWeighting = false;
	dirty = true;

	lowWeight = .0;
	highWeight = .0;
	totWeight = .0;
}

QuantileFunction::QuantileFunction(QuantileFunction q, int size = 10)
{
	rawData.clear();
	qf = NULL;
	setQFSize(size);

	int numE = q.rawData.size();
	for (int i = 0; i < numE; i ++)
	{
		rawData.push_back(q.rawData[i]);
	}

	useMyWeighting = false;
	dirty = true;

	lowWeight = .0;
	highWeight = .0;
	totWeight = .0;
	
	setCutOffs();
}

QuantileFunction::~QuantileFunction()
{
	rawData.clear();
	delete []qf;
}

void QuantileFunction::useDefaultWeighting(double umu, double uvar)
{
	useMyWeighting = true;
	mu = umu;
	var = uvar;
}

void QuantileFunction::setQFSize(int size)
{ 
	dirty = true;

	if (qfSize != size)
	{
		qfSize = size;

		if (qf != NULL)
		{
			delete []qf;
			qf = new double[qfSize];
		}
	}
	
	if (qf == NULL)
		qf = new double[qfSize];

	memset(qf, 0, sizeof(double)*qfSize);

}

void QuantileFunction::addData(ImageElement element)
{
	rawData.push_back(element);
	dirty = true;
}

void QuantileFunction::addData(double element, Vector3D coord)
{
	ImageElement de;

	de.element = element;
	de.coord = coord;

	//Now for the weight.  If the user sent a double and not a ImageElement, 
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

//This is weird.  I want to add an image element.  the weight of the image element
//comes from a product of the profile's weight and the distance along tau.  If I'm 
//trying to combine patches, and thus their qf's,  I want the weight of this element
//to be different. How? The easy way is to force a resample (higher up) of the patch
//of this quantile function.  The maybe faster way is to be sure to change the 
//element's weight correctly, relative to how the patches are added. Just know that
//the weight of the element should be changed to be simply the tau weight if the 
//profile that sampled it is now one of the interior profiles (and so woulld likely 
//have weight 1).  

//Final thoughts:  a simple drawing says that as long as the two profiles that sampled
//this coord are coming in from about the same direction (the surface is relatively
//low curvature), then the sum of the weights should be about the right thing to do.
//w1 = p1*t1, w2 = p2*t2.  if for the patch combine we sum p1 and p2, then great.

//You know, I'll probably remember this and just resample after the patch combine.
//To avoid the linear time single insert (quadratic total), I will sort the combined
//vector by coordinate, then combine equal coord elements in a single linear run
//through.

//Define the sort method.

bool cmpImageElements (const ImageElement &a, const ImageElement &b)
{
	if (a.coord.getZ() < b.coord.getZ())
		return true;
	else
	{
		if (a.coord.getZ() == b.coord.getZ())
		{
			if (a.coord.getY() < b.coord.getY())
				return true;
			else
			{
				if (a.coord.getY() == b.coord.getY())
				{
					if (a.coord.getX() < b.coord.getX())
						return true;
				}
			}
		}
		return false;
	}
}

void QuantileFunction::addQF(const QuantileFunction &q)
{

	int i, j, numE;// numNE = q.rawData.size(), numE = rawData.size();


	//The easy way to add one quantile function to another is to take a weighted 
	//average.  But how different is that from the right way, which is to mix the 
	//underlying raw data?  Not very if we're talking about a unimodal mix.  But, an 
	//additional confusion is that many of the voxels may be in common here, when 
	//they're assumed indpendent in the matlab training (when we take a simple
	//weighted average. 

	//The right way! Interpolate the CDF's, which are the inverses of the qfs.  Pretty
	//simple (see matlab/.../mixQF).  Just need an interp1 for C++. 
	//Did it!. Now I'll port it over.
	//I'm going to comment this stuff out if it works.
	//It works, enough anyway.  The mixed qf is not quite the real pooled qf, but 
	//since I only have the mixed qf in matlab, I should used mixedqfs here (well, in 
	//RegionShift anyway) as well.

	//Going to compare my mixed qf to the actual qf after pooling the data.

	//This qf mixing debuggery has been replaced by the larger scale addition and 
	//subtraction mixing commented out in RegionShift::getMatch.
/*
	double *qf0 = getQF();
	int qf0Size = qfSize;
	const double *qf1 = q.getQF();
	int qf1Size = q.getQFSize();
	double w0 = getMiddleWeight();
	double w1 = q.getMiddleWeight();
	double alpha = w1 / (w0 + w1);

	double *qfInt = new double[qfSize];
	mixQF(qf0, qf1, qf0Size, qf1Size, qfInt, qfSize, alpha);
	//Now I'm going to compare qfInt to getQF after this combining of data.

	FILE *fout = fopen("QFIntVSReal.txt", "w");
	//I'm going to write out the two qfs as they are now.  then the interp then the after.

	fprintf(fout, "%lf\n", alpha);
	for (i = 0; i < qf0Size; i++)
		fprintf(fout, "%lf ", qf0[i]);

	fprintf(fout, "\n");

	for (i = 0; i < qf1Size; i ++)
		fprintf(fout, "%lf ", qf1[i]);

	fprintf(fout, "\n");
*/	



	//insert the (perhaps) duplicate elements of the other qf.
	rawData.insert(rawData.end(), q.rawData.begin(), q.rawData.end());

	//Now sort on coordinate.
	sort(rawData.begin(), rawData.end(), cmpImageElements);

	//Now the elements are sorted on coord.  I can go through the vector once,
	//combining duplicate elements.

	numE = rawData.size();

	//std::vector<ImageElement>::iterator tmp; //the one to erase below.
	
	i = j = 0;

	while (i < numE)
	{
		j = i;
		i ++;
		while (i < numE && rawData[i].coord == rawData[j].coord)
		{
			//add the weight to this sample
			rawData[j].weight += rawData[i].weight;
			//delete the extraneous sample.
			rawData.erase(rawData.begin() + i);
			//reset size.
			numE--;
		}

		//set weight to one if it must be.
		if (rawData[j].weight > 1.0)
			rawData[j].weight = 1.0;
	
	}


	dirty = true;

	
	//More for the comparing of pooled to mixed qfs.
/*	
	//Now that we've combined the data, let's look at the quantile function.
	double *realQf = getQF();

	//I want to compare real with average qf's.  Write them out

	for (i = 0; i < qfSize; i ++)
		fprintf(fout, "%lf ", qfInt[i]);

	fprintf(fout, "\n");

	for (i = 0; i < qfSize; i ++)
		fprintf(fout, "%lf ", realQf[i]);

	fclose(fout);
	delete [] qfInt;
*/
	
	


}

void QuantileFunction::setCutOffs(double inLowCutOff, double inHighCutOff)
{
	lowCutOff = inLowCutOff;
	highCutOff = inHighCutOff;
}

//This must be defined for the sort of the vector.
bool operator<(const ImageElement& a, const ImageElement& b) 
{
    return a.element < b.element;
}

void QuantileFunction::makeQF(int size)
{
	int i, numE = rawData.size(), istart = 0, iend = 0;
	//if (numE < qfSize)
	//	cerr << "Not enough data to makeQF. numE " << numE << " and qfSize " << qfSize << ".\n";

	
	setQFSize(size);


	//The relatively slow way: This way loops through the rawData more
	//than once. A potentially much better way (assuming there are lots
	//of elements in rawData) would be to make a histogram at, say, 4*qfSize
	//resolution, and then loop through that.  But I can always come back.
	std::sort(rawData.begin(), rawData.end());

	//Now rawData is sorted.  Take the average of each proportion of weight.
	//First find the sum of the weights.
	double sum = 0.0, wpb = 1.0, wpbInv = 1.0;

	//While I'm computing the total weight, I can also update low and high weight.
	lowWeight = .0;
	highWeight = .0;
	totWeight = .0;

	if (numE == 0)
	{
		cout << "ERROR QuantileFunction: no data.\n";
	}

	for (i = 0; i < numE && rawData[i].element <= lowCutOff; i ++)
	{
		totWeight += rawData[i].weight;
	}
	lowWeight = totWeight;
	istart = i;

	for (; i < numE && rawData[i].element < highCutOff; i ++)
	{
		totWeight += rawData[i].weight;
	}
	iend = i;

	for (; i < numE; i ++)
	{
		highWeight += rawData[i].weight;
	}
	totWeight += highWeight;
	sum = totWeight - lowWeight - highWeight;

	//Remember, sum should be without the low and high.



//cout << sum << ".\n";
	if (sum < 1.0e-10 || istart == iend)
	{
		cout << "Total weight quite small: " << sum << ". low is " << lowWeight << " and high is " << highWeight 
			 << ". No quantile function possible.\n";

		//This basically means that the thousands of bone voxels are so far away
		//as to have a miniscule effect on the result.

		return;
	}
	else
	{
		wpb = sum/qfSize;
		wpbInv = 1.0/wpb;
	}
	
	//loop through the data again.
	int currentBin = 0;
	double partialWeight, currentWeight = 0.0;
	i = istart;


	//now make qf for the data in range
	while (currentBin < qfSize && i < iend)
	{
		while (i < iend && (currentWeight + rawData[i].weight < wpb))
		{
			qf[currentBin] += rawData[i].weight * rawData[i].element;
			currentWeight += rawData[i].weight;
			i++;
		}
		//get the partial weight and move on, unless i >= numE.
		if (i < iend)
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
			if (++currentBin < qfSize)
				qf[currentBin] += partialWeight * rawData[i].element;
			else
				assert(partialWeight < 1.0e-5);
			
			currentWeight = partialWeight;
			i ++;
		}
	}

	for (currentBin = 0; currentBin < qfSize; currentBin ++)
	{
		qf[currentBin] *= wpbInv;
	}
	
}

double * QuantileFunction::getQF(int size)
{
	if (!dirty)
		return qf;

	if (size > 0)
		makeQF(size);
	else
		makeQF(qfSize);

	dirty = false;
	return qf;
}

const double * QuantileFunction::getQF(int size) const
{
	if (!dirty)
		return qf;

	return NULL;
}


double QuantileFunction::getMean()
{
	if (dirty)
		makeQF(qfSize);

	double sum = 0.0;

	for (int i = 0; i < qfSize; i ++)
	{
		sum += qf[i];
	}

	dirty = false;
	return sum/qfSize;
}

double * QuantileFunction::getWeights()
{
	if (dirty)
		makeQF(qfSize);

	double *w = new double[3];

	w[0] = lowWeight;
	w[1] = highWeight;
	w[2] = totWeight;
	dirty = false;

	return w;
}






/********************************************************************************************************************************************************/



/********************************************************************************************************************************************************
Code for the PatchPCA class,.
*********************************************************************************************************************************************************/




/********************************************************************************************************************************************************/


PatchPCA::PatchPCA()
{
	mu = NULL;
	var = NULL;
	eigs = NULL;
	weights = NULL;

	neigs = 0;
	dim = 0;

	gasMatters = true;
	boneMatters = true;
}

PatchPCA::~PatchPCA()
{
	delete []weights;
	delete []mu;
	delete []var;
	for (int i = 0; i < neigs; i ++)
		delete []eigs[i];
}

void PatchPCA::setPca(double *inPca, int inDim, int length)
{
	int i,j;
	dim = inDim;

	mu = new double[dim];
	var = new double[dim];

	int c = 0;
	for (i = 0; i < dim; i ++)
	{
		mu[i] = inPca[c];
		c++;
	}

	neigs = 0;
	for (i = 0; i < dim; i ++)
	{
		var[i] = inPca[c];
		c++;

		if (var[i] > 0)
			neigs ++;
	}

	//What if there is no residue?
	if (length < (neigs + 2)* dim)
		neigs --; //for the residue.

	if (length < (neigs + 2)* dim)
	{
		cout << "ERROR PatchPCA::setPca:length " << length << " not appropriate, given dim " << dim 
			 << " and neigs " << neigs << ".Will probably seg fault.\n";
	}

	eigs = new double *[neigs];

	for (j = 0; j < neigs; j ++)
	{
		eigs[j] = new double[dim];

		for (i = 0; i < dim; i ++)
		{
			eigs[j][i] = inPca[c];
			c++;
		}
	}

	//Now, this pca object's pca component should be populated.
	delete []inPca; 
	//this is because the registry, where we got this, gave us control, and
	//the memory is about to stop being referenced.
}

double PatchPCA::computeMahal(double *qf, bool normalize)
{
	//this function uses the pca info to compute a mahalanobis distance
	//of the arg qf.  The normalization must happen elsewhere.
	//That's not true.  I'll do it here, since I know neigs.

	//If the qf is all zeros, just quit.  this is also a hack.just hope you don't 
	//have pca's with mean the zero qf.

	if (qf[0] == 0 && qf[dim-1] == 0)
		return 0.0;

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

	if (neigs < dim) //do residue.
	{
		double temp = sqrt(distInSpace) - sqrt(distInPCASubspace);
		mahal += temp*temp/var[neigs];
	}
	//Because distInPCASubspace is always less than distInSpace, don't have to check 
	//the sign.

	//Now, finally, make the result MU zero and unit STD/VAR, so to speak...if everything's 
	//Gaussian...It has to do with the fact that the Mahal distance of a Gaussian variable 
	//follows a Chi-square distribution with degree of freedom equal to the dimension of the
	//Gaussian. so expected value neigs + 1, with var = 2*(neigs + 1).

	if (normalize)
	{
		mahal -= neigs;
		if (neigs < dim)
		{
			mahal -= 1;
			mahal /= sqrt(2.0*(neigs + 1));
		}
		else
			mahal /= sqrt(2.0*neigs);
	}

	//deletions
	delete [] diff;
	delete [] coeffs;

	return mahal;
}

double PatchPCA::computeMatch(double *qf, double *w, bool normalize)
{
	//This is where we define the low and high matches.  These ought to be worth about nothing,
	//since they're so variable. But that's what the gasMatters and boneMatters tells us.

	//These matches are really penalties...
	
	//Definitely want the mahal.
	double mahalMatch = computeMahal(qf, normalize);


	//That was easy.
	double gasMatch = 0.0, boneMatch = 0.0;
	int count = 1;

	//w should be [lowWeight highWeight totWeight], and I may assume more stuff later. 
	//if gasMatters, then we want to penalize the presence of gas, or at least relative 
	//to the training.  But on the inside of the rectum for example, we don't want to 
	//penalize gas. We don't want to say more relative gas is better either, cause it'll
	//try to shrink itself.  So just give the best match.

	//You also see that weights is used as [lowMean highMean lowVar highVar];

	gasMatch = computeGasMatch(w, normalize);
	if (gasMatters)
		count ++;

	boneMatch = computeBoneMatch(w, normalize);
	if (boneMatters)
		count ++;


	//just return something...

	if (normalize)
		return (mahalMatch + gasMatch + boneMatch) / sqrt(1.0*count);
	else
		return (mahalMatch + gasMatch + boneMatch);
}

//The code that uses this, in MultiPatch will assume the dimension is 1 for 
//normalization purposes.
double PatchPCA::computeGasMatch(double *w, bool normalize)
{
	double gasMatch = 0.0;
	if (w == NULL) return gasMatch;

	if (gasMatters)
	{
		gasMatch = ((w[0]/w[2]) - weights[0]);
		gasMatch = gasMatch*gasMatch / weights[2];
		
		if (normalize)
		{
			gasMatch --;
			gasMatch /= sqrt(2.0);
		}
	}

	return gasMatters*gasMatch;
}

//The code that uses this, in MultiPatch will assume the dimension is 1 for 
//normalization purposes.
double PatchPCA::computeBoneMatch(double *w, bool normalize)
{
	double boneMatch = 0.0;
	if (w == NULL) return boneMatch;

	if (boneMatters)
	{
		boneMatch = ((w[1]/w[2]) - weights[1]);
		boneMatch = boneMatch*boneMatch / weights[3];
		
		if (normalize)
		{
			boneMatch --;
			boneMatch /= sqrt(2.0);
		}
	}

	return boneMatters*boneMatch;
}

void PatchPCA::setWeights(double *inWeights, int length)
{//I do the copy just in case this memory addr is better suited than wherever that came from.
	//Right now I'll work like there's [lowMean highMean lowVar highVar totMean totVar] in 
	//here.  maybe more later.

	if (weights != NULL)
		delete []weights;

	weights = new double[length];
	for (int i = 0; i < length; i ++)
	{
		weights[i] = inWeights[i];
	}

	//Also, have a minimum variance to avoid the zeros.
	weights[2] = (weights[2] < 1.0e-9) ? 1.0e-9 : weights[2];
	weights[3] = (weights[3] < 1.0e-9) ? 1.0e-9 : weights[3];

	//The above should have been taken care of by the matlab.  Anyway.


	delete []inWeights; 
}

double ** PatchPCA::projectData(double ** qfs, int num)
{

	double ** res = new double *[num];
	int i, j, p;

	double *diff, *coeffs;

	diff = new double[dim]; //the difference vector;
	coeffs = new double[neigs];
	//the coefficients of the projection.
	
	for (p = 0; p < num; p ++)
	{
		res[p] = new double[neigs];

		for (i = 0; i < dim; i++)
		{
			//diff[i] = qfs[i][p] - mu[i];
			//Going to swap i and p, so that it's easier
			//to build the qfs that I send here.
			diff[i] = qfs[p][i] - mu[i];
		}
		
		for (j = 0; j < neigs; j ++)
		{
			res[p][j] = 0.0;
			
			//Simple dot product.
			for (i = 0; i < dim; i ++)
			{
				res[p][j] += diff[i]*eigs[j][i];
			}
			
		}
	}

	return res;
}





/********************************************************************************************************************************************************/



/********************************************************************************************************************************************************
Code for the Patch class,.
*********************************************************************************************************************************************************/




/********************************************************************************************************************************************************/



Patch::Patch(int inFigureID, int inSpokeID, Bpoint bp, SamplingConfiguration *inSc)
{
	figureID = inFigureID;
	sc = inSc;

	spokeIDs.clear();
	anchorPoints.clear();

	spokeIDs.push_back(inSpokeID);
	anchorPoints.push_back(bp);

	data.setQFSize(sc->emSize);
	data.setCutOffs(sc->lowCutOff, sc->highCutOff);

	pca = NULL;
}

Patch::Patch()
{
	//Nothing to set really, I think anyway.
	spokeIDs.clear();
	anchorPoints.clear();

	sc = NULL;
	pca = NULL;
}

Patch::~Patch()
{
	clear();
	profiles.clear();
	spokeIDs.clear();
	anchorPoints.clear();

	data.clearData();
	//delete data;

	delete pca;
}

void Patch::clear()
{
	spokeIDs.clear();
	anchorPoints.clear();

	profiles.clear();
	data.clearData();
}


void Patch::sampleYourself(Image3D *image, Bpoint *Bpoints, Image3D *patchImage, GreyValue g, Image3D *binImage)
{
	int i;
	Vector3D coord, coord2;
	double mahal;
	double maxDistanceSquared = sc->samplingDistance * sc->samplingDistance;
	if (!sc->isRRelative)
		maxDistanceSquared /= (image->getModelToWorldScale() * image->getModelToWorldScale());

	double mahalSTDMult = 1.0 / (-2.0 * sc->distanceSTD * sc->distanceSTD);

	if (sc->distanceSTD <= 0.0)
		mahalSTDMult = 0.0;

	//First, clear the qf data. This is what we're going to repopulate now anyway.
	data.clearData();

	ImageElement ie;

	std::set<ProfileComponent, compareProfileComponents>::iterator it;

	std::vector<double>::iterator tauIt;
	std::vector<double>::iterator tauWeightIt;

	for (i = 0, it = profiles.begin(); i < profiles.size(); i ++, it++)
	{
		//Bpoint mybp = Bpoints[it->index];
		//profiles[i] is the index into Bpoints and the along surface weight of the profile.
		Vector3D Pi(Bpoints[it->index].pnt);
		Vector3D Ni(Bpoints[it->index].norm);

		for (tauIt = sc->taus.begin(), tauWeightIt = sc->tauWeights.begin(); tauIt != sc->taus.end(); tauIt++, tauWeightIt++)
		{
			//The last bit is whether this is r-relative.
			//That business has to be done now, at sample time. if r-relative, then
			//mult by r, if not, then mult by worldToModelScale (which has been done).  Again, see 
			//SurfacePatchEnsemble and Eli's code.  weights will be the same regardless.

			if (sc->isRRelative)
				coord = Pi + Ni * Bpoints[it->index].rad * (*tauIt);
			else
				coord = Pi + Ni * (*tauIt);


			image->modelToImageCoordinates(coord);

			coord.round();


			if (!image->clipToImage(coord))
				continue;

			//if ((!sc->isInside && patchImage->getVoxelValue(coord.getX(), coord.getY(), coord.getZ()))
			//	|| (sc->isInside && patchImage->getVoxelValue(coord.getX(), coord.getY(), coord.getZ()) ) )
			if (patchImage->getVoxelValue(coord.getX(), coord.getY(), coord.getZ()) == g)
				continue;

			//there's a rought bit here.  My inside/outside regions are not paired, as in Eli.  So when he's 
			//checking inside versus outside, if the binImage is defined he simply uses that to decide, where
			//my inside patch for example will never sample outside of the model.  
			//In order to improve that (do it the wrong way Eli does, see SurfacePatchEnsemble::update), I 
			//will change my MultiPatch::updateTauAndWeights to sample a little on the wrong side as well 
			//as the correct side given inside/outside.  Then here, if there's no binary, use that fact.

			coord2 = coord;
			image->imageToModelCoordinates(coord2);
			mahal = (Ni * (coord2 - Pi));
			//Now coord2 is in model space and mahal is the normal projected distance from center of the voxel.  

			if (binImage != NULL)
			{	
				if ((sc->isInside && !binImage->getVoxelValue(coord.getX(), coord.getY(), coord.getZ()))
					|| (!sc->isInside && binImage->getVoxelValue(coord.getX(), coord.getY(), coord.getZ())))
					continue;
			}
			else
			{ 
				//This ensures that if we don't have the binary, then the tau tells us whether to 
				//sample or not.  if inside patch don't sample if tau > 0, for example.
				//if ((sc->isInside && *tauIt > 1.0e-5) || (!sc->isInside && *tauIt < -1.0e-5))
				//	continue;

				//We're going to do the normal dot product check that Eli does 
				//for inside/outside, instead of the one above.  These are usually highly weighted voxels on the
				//boundary, where the commented out code above and this provide different results.

				if ((sc->isInside && mahal > 0) || (!sc->isInside && mahal <= 0))
					continue;
			}

			//patchImage->setVoxel(coord.getX(), coord.getY(), coord.getZ(), g); moved down


			ie.element = (double) image->mapDisplayToActual(image->getVoxelValue(coord.getX(), coord.getY(), coord.getZ()));
			ie.coord = coord;

			//ie.weight = it->weight * (*tauWeightIt);
			//mahal is a distance in model coordinates and the weight of the voxel should be based on this, and
			//not on the predefined tau step that caught this voxel.

			//The r-relative stuff is confusing kind of.  If I take out the r in mahal, then I can run the same check 
			//against maxDist, but then the mahalSTDMult (from distanceSTD) doesn't make any sense...So I'm not 
			//going to take out the radius, but I'll check against a different maxDist
			//if (sc->isRRelative)
			//	mahal /= Bpoints[it->index].rad;

			mahal *= mahal;

			double radd = Bpoints[it->index].rad;

			if (sc->isRRelative && mahal >= maxDistanceSquared * radd * radd)
					continue;
			else if (mahal >= maxDistanceSquared)
					continue;

			patchImage->setVoxel(coord.getX(), coord.getY(), coord.getZ(), g);
			//I put this here, cause it should happen after all possibility of quiting.

			//We have to account for r-relative again and the distanceSTD.  WAIT! No.  the mahal is an r multiplied
			//dist in the model frame, and the distanceSTD is r-relative too if r-relative is on.  So nothing needs to 
			//be done here.

			ie.weight = it->weight * exp(mahal * mahalSTDMult);
			//It could be useful to compare exp(mahal * mahalSTDMult) with (*tauWeightIt), as they should
			//have a cos relationship with the two being equal if the voxel center is on the line Ni
			//(which in turn should be degenerately rare). 

#ifdef WIN32
			if (_isnan(ie.weight))
#else
			if (isnan(ie.weight))
#endif
			{
				ie.weight = 0.0;
			}



			data.addData(ie);
			//if (ie.weight < 1.0)
			//	double qq = it->weight;
		}
	}//do i agree with the totalWeights that come out?

	//cout << "Patch " << spokeIDs[0] << " sampled. Size " << data.getSize() << "\n";

	//int ss = data.getSize();

	return;
}


bool operator<(const ProfileComponent& a, const ProfileComponent& b) 
{
    return a.index < b.index;
}


//To add a patch to another, you want.
void Patch::addPatch(const Patch &p)
{
	int i, j;

	//Patches should only be added if their SamplingConfigurations are the same.
	//I'm not puting that check here, but it should be here, something like:
	//if (sc->distanceSTD != p.sc->distanceSTD || sc->distance != p.sc->distance)
	//	error


	//add anchorPoints and spokeIDs.

	if (sc == NULL)
		setSC(p.sc);

	if (sc != p.sc)
		cout << "ERROR: Patches with different sampling configs cannot combine.\n";


	int numNS = p.anchorPoints.size(), numS = anchorPoints.size();
	
	for (i = 0; i < numNS; i ++)
	{
		for (j = 0; j < numS && p.anchorPoints[i] != anchorPoints[j]; j ++);

		if (j == numS)
		{
			spokeIDs.push_back(p.spokeIDs[i]);
			anchorPoints.push_back(p.anchorPoints[i]);
		}
	}


	//add the profiles that need to be added.
	//Also, I set the weight of an in common profile to be the sum
	//of the two weights with a max of 1.0.
	int numNE = p.profiles.size(), numE = profiles.size();

	std::set<ProfileComponent, compareProfileComponents>::const_iterator it;
	std::set<ProfileComponent, compareProfileComponents>::iterator found;

	for (i = 0, it = p.profiles.begin(); i < numNE; i ++, it ++)
	{
		found = profiles.find(*it);
		if (found != profiles.end())
		{
			ProfileComponent tempP = *found;
			profiles.erase(found);
			
			if (tempP.weight + it->weight <= 1.0)
				tempP.weight += it->weight;
			else
				tempP.weight = 1.0;

			profiles.insert(tempP);
		}
		else
			profiles.insert(*it);
	}

/*
	for (i = 0; i < numNE; i ++)
	{
		for (j = 0; j < numE && p.profiles[i].index != profiles[j].index; j ++);

		if (j == numE)
		{
			profiles.push_back(p.profiles[i]);
			numE++;
		}else
		{
			if (profiles[j].weight + p.profiles[i].weight < 1.0)
				profiles[j].weight += p.profiles[j].weight;
			else
				profiles[j].weight = 1.0;
		}
	}
*/
	//Lastly, don't forget to add the quantile functions of the two patches.  
	//This is really what it's about.  The above is adding the ability to 
	//sample the combined patch.  In optimization, for a given prospective
	//model, the idea will be to add a neighboring patch to this one, when 
	//the patches have already been sampled.  
	data.addQF(p.data);

}

double Patch::getProfileWeight(int index)
{
/*	std::set<ProfileComponent, compareProfileComponents>::const_iterator it;
	for (it = profiles.begin(); it != profiles.end(); it ++)
	{
		if (it->index == index)
			return it->weight;
	}

	return 0.0;
*/
	//second way, much faster.
	std::set<ProfileComponent, compareProfileComponents>::const_iterator found;
	ProfileComponent pc;
	pc.index = index;

	found = profiles.find(pc);
	if (found != profiles.end())
	{
		return found->weight;
	}
	
	return 0.0;

}





/********************************************************************************************************************************************************/



/********************************************************************************************************************************************************
Code for the RegionShift class,.
*********************************************************************************************************************************************************/




/********************************************************************************************************************************************************/

RegionShift::RegionShift(PatchPCA *inProjPCA, std::vector< PatchPCA * > inTypePCAs, 
						 double *inP, int priorSize, SamplingConfiguration *inSc)
{
	int i, j, ind;

	projPCA = inProjPCA;
	sc = inSc;

	priors = new double[priorSize];
	memcpy(priors, inP, sizeof(double)*priorSize);

	typePCAs.clear();

	typePCAs.assign(inTypePCAs.begin(), inTypePCAs.end());

	for (i = 0; i < typePCAs.size(); i ++)
	{		
		//Now for error catching purposes, let's see if the sizes and all that match up.
		//With the combined regions, I don't use projPCA so it's trash.  just ignore the 
		//error.
		if (typePCAs[i]->getDim() != projPCA->getNEigs())
		{
			//printf("MultiPatch. RegionShift error: proj pca neigs (%d) not equal to dim (%d) of type pca %d", 
			//	projPCA->getNEigs(), typePCAs[i]->getDim(), i);
		}
	}
	
	numTypes = typePCAs.size();

	numPts = priorSize / numTypes;  //Should be a whole number, so I'll do integer division here.

	initModel = new int[numPts];
	finalModel = new int[numPts];
	permutation = new int[numPts];
	weights = new double[numPts];
	memset(weights, 0, numPts*sizeof(double));

	numModes = 0;

	//int howmuch[] = {0,0,0,0};
	double max;
	for (i = 0; i < numPts; i ++)
	{
		//independent bit, just initialize permutation.
		permutation[i] = i;

		max = 0.0;
		for (j = 0; j < numTypes; j ++)
		{
			ind = numTypes*i + j;
			if (priors[ind] > max)
			{
				max = priors[ind];
				initModel[i] = j;
			}
		}
		//howmuch[initModel[i]]++;
	}
}

void RegionShift::buildAdjacency(const M3DObject * object, int f, int inSurfL)
{
	surfL = inSurfL;

	ThallCode::Pointlist_server2 pList;
    M3DQuadFigure * figure;
	int npoints;
    Bpoint * pointset;
    int i, j;
	int num_n = 0;
	int max_neighbors = 8;
	
    
	figure = (M3DQuadFigure *) object->getFigurePtr(f);
	if (figure == NULL)
		return;
	
	// Generate the surfaces, for subdivlevels 0 through 3.
	Xferlist * xferList = convertM3DtoXfer(figure);
	SubdivBoundary * boundary = figure->getBoundaryPtr();
	
	pList.init(xferList);
	pList.ComputeSubdivPointCloud(surfL);
	pList.ComputeSubdivBoundaryTiles(surfL);
	
	pList.subdivboundaryinfo(&npoints, &pointset);

	if (npoints != numPts)
	{
		printf("RegionShift error: number of points given by prior (%d) not consistent with \n"
			"that given by adjacency (%d) at surface level %d.", numPts, npoints, surfL);
	}
	
	adj = new int*[numPts];
	
	ThallCode::Tileneighbors * vn = NULL;
	pList.subdivvertexneighbors(&num_n, &vn);
	
	//Now output the neighbor indices for every point.
	for (i = 0; i < npoints; i ++)
	{
		adj[i] = new int[max_neighbors];
		
		for (j = 0; j < max_neighbors; j ++)
		{
			adj[i][j] = vn[i][j];
		}
	}
	
	delete [] (xferList->atomlist);
	delete xferList;
	
	
}

RegionShift::~RegionShift()
{
	delete [] initModel;
	delete [] finalModel;
	delete [] priors;
	delete [] permutation;

	delete projPCA;

	int i;
	for (i = 0; i < typePCAs.size(); i ++)
		delete typePCAs[i];

	typePCAs.clear();

	for (i = 0; i < numPts; i ++)
		delete [] adj[i];
}

double RegionShift::getMatch(double ** qfs, double &totModes)
{
	//Let's start by projecting the data.
	double ** cX = projPCA->projectData(qfs, numPts);
	
	srand( time(NULL) );

	double beforeMatch = computeShiftedMatch(initModel, cX);

	doShifting(cX);
	//memcpy(finalModel, initModel, sizeof(int)*numPts);

	double theMatch = computeShiftedMatch(finalModel, cX);

	printf("RegionShift match: from %lf to %lf.\n", beforeMatch, theMatch);

	totModes = 1.0*numModes;

	
	/*
	//Debug.  Just going to combine patches a little, to see whether my mixQF is
	//working and how accrate it is if regions overlap.

	Patch * patch = new Patch();
	double **componentqfs = new double*[6];
	double **interpolants = new double*[6];
	double **intermediates = new double*[6];	

	double mytotalweights[6];
	double myweights[6];
	int n = 0, pt = 15, cur = 0;

	//Let's take patch zero and its neighbors and combine.
	patch->addPatch(*patches[pt]);
	patch->setSC(patches[pt]->getSC());
	int mysize = patch->getEMSize();

	double *qfInt = new double[mysize];
	memset(qfInt, 0, mysize*sizeof(double));

	componentqfs[cur] = new double[mysize];
	double *temp = patches[pt]->getQF();
	myweights[cur] = patches[pt]->getMiddleWeight();
	memcpy(componentqfs[cur], temp, mysize*sizeof(double));

	mytotalweights[cur] = patch->getMiddleWeight();
	interpolants[cur] = new double[mysize];
	double alpha = 1.0;
	mixQF(qfInt, componentqfs[cur], mysize, mysize, interpolants[cur], mysize, alpha);

	intermediates[cur] = new double[mysize];
	temp = patch->getQF();
	memcpy(intermediates[cur], temp, mysize*sizeof(double));
	cur ++;
	int i; 
	double weightsofar = 0.0;


	while (adj[pt][n] != -1)
	{
		if (finalModel[adj[pt][n]] == finalModel[pt])
		{
			patch->addPatch(*patches[adj[pt][n]]);
			
			
			componentqfs[cur] = new double[mysize];
			temp = patches[adj[pt][n]]->getQF();
			myweights[cur] = patches[adj[pt][n]]->getMiddleWeight();
			memcpy(componentqfs[cur], temp, mysize*sizeof(double));
			
			mytotalweights[cur] = patch->getMiddleWeight();
			interpolants[cur] = new double[mysize];
			alpha = myweights[cur];
			weightsofar = 0.0;
			for (i = cur; i >= 0; i --) weightsofar += myweights[i];
			alpha /= weightsofar;
			mixQF(interpolants[cur-1], componentqfs[cur], mysize, mysize, interpolants[cur], mysize, alpha);
			
			intermediates[cur] = new double[mysize];
			temp = patch->getQF();
			mytotalweights[cur] = patch->getMiddleWeight();
			memcpy(intermediates[cur], temp, mysize*sizeof(double));
			cur ++;
		
		}
		n++;
	}

	//Now debug if negative mixtures works. Go through the neighbors again,
	//and subtract them from the total.  And positive for that matter, just the whole
	//thing.
	n = 0;


	FILE *fout = fopen("QFMixturesAndNegativeMixtures.txt", "w");
		
	for (n = 0; n < cur; n++)
		fprintf(fout, "%lf ", myweights[n]);
	fprintf(fout, "\n");
	for (n = 0; n < cur; n++)
		fprintf(fout, "%lf ", mytotalweights[n]);
	fprintf(fout, "\n\nThe local regional qfs\n\n");

	//The local regional quantile functions.
	for (n = 0; n < cur; n++)
	{
		for (i = 0; i < mysize; i ++)
		{
			fprintf(fout, "%lf ", componentqfs[n][i]);
		}
		fprintf(fout, "\n");
	}
	

	fprintf(fout, "\n\nThe interpolants I compute\n\n");
	//Now the interpolants. These are NOT the actual ones, those come next.
	for (n = 0; n < cur; n++)
	{
		for (i = 0; i < mysize; i ++)
		{
			fprintf(fout, "%lf ", interpolants[n][i]);
		}
		fprintf(fout, "\n");
	}


	fprintf(fout, "\n\nThe intermediates that the interpolants should be.\n\n");
	//Now the intermediates. These are the actual ones, not the interpolants.
	for (n = 0; n < cur; n++)
	{
		for (i = 0; i < mysize; i ++)
		{
			fprintf(fout, "%lf ", intermediates[n][i]);
		}
		fprintf(fout, "\n");
	}

	fprintf(fout, "\n\nThe out mixed, compare to the intermediates.\n");
	//Now out mix the quantile functions and write them.

	//I wouldn't have mytotalweights, the accurate ones, without 
	//building the larger region.  So compute the ones to use here.  
	//Just a cumsum of myweights. Also, I use interpolants, not intermediates.
	double cumsummyweights[6];
	cumsummyweights[0] = myweights[0];
	for (n = 1; n < cur; n++)
		cumsummyweights[n] = myweights[n] + cumsummyweights[n-1];


	for (n = cur-1; n >= 0; n --)
	{
		//compute the alpha to rid the intermediate of the component.
		alpha = myweights[n] / cumsummyweights[n];
		alpha = -1.0*alpha / (1.0 - alpha);

		//Now mix.
		mixQF(interpolants[n], componentqfs[n], mysize, mysize, qfInt, mysize, alpha);

		//Now write.
		for (i = 0; i < mysize; i ++)
		{
			fprintf(fout, "%lf ", qfInt[i]);
		}
		fprintf(fout, "\n");
	}



	fclose(fout);
	//This is the end of the qf mixing debuggery.
	*/


	return theMatch;
}

void RegionShift::doShifting(double ** cX)
{
	//X here is assumed to be the dim of the typePCAs by numPts.

	int pt, i, n;
	
	memcpy(finalModel, initModel, sizeof(int)*numPts);
	
	int numChanged = 1, cur = 0, bestType, totalChanged = 0, temp;
	bool onBorder;

	double bestMatch;
	
	int *possibilities = new int[numTypes];
	double *scores = new double[numTypes];
	
	while (numChanged > 0)
	{
		numChanged = 0;
		//Get a permutation.
		Permute();
		
		for (i = 0; i < numPts; i ++)
		{
			//if none of the adjacent verts to this point are different from it, then 
			//don't do anything.
			
			memset(possibilities, 0, numTypes*sizeof(int));
			memset(scores, 0, numTypes*sizeof(double));
			
			pt = permutation[i];
			
			possibilities[finalModel[pt]] = 1;
			

			//Here we loop through the neighbors, determining which types amongst our
			//neighbors are at once different from us (finalModel[pt]) and possible for 
			//us to be.
			n = 0; 
			onBorder = false;

			while (adj[pt][n] != -1)
			{
				temp = finalModel[adj[pt][n]]; //the neighboring type.
				if (temp != finalModel[pt]	//is not us.
					&& priors[pt*numTypes + temp] > 0 ) //but we can be it.
				{
					onBorder = true;
					possibilities[temp] = 1;
				}
				n ++;
			}
			
			if (!onBorder)
				continue;
			
			
			//Now we're going to score the point with respect to the possibilities.
			
			for (n = 0; n < numTypes; n ++)
			{
				if (possibilities[n])
					scores[n] = typePCAs[n]->computeMahal(cX[pt], true);
			}
			//We normalize since the decision on whether to change is based on an
			//energy including other terms. These terms should be commensurate, so
			//that the terms can offset one another.
			
			
			//Now we do our energy function thing.  
			//First stab: find the best score for non-zero prior probability.  The non-zero
			//prior was in the checking agains the neighbors.
			
			bestMatch = 1.0e10; bestType = -1;
			for (n = 0; n < numTypes; n ++)
			{
				if (possibilities[n] && scores[n] < bestMatch)
				{
					bestType = n;
					bestMatch = scores[n];
				}
			}
			
			if (bestType == -1)
			{
				printf("RegionShift: bestType is bad, bestMatch found is %lf.\n", bestMatch); 
			}

			if (bestType != finalModel[pt]) //we're going to change the type of a point.
			{
				finalModel[pt] = bestType;
				numChanged ++;
				totalChanged ++;
				
			}
			
		}
		
		cur ++;
	}

	//Okay, done.  So now finalModel contains the best type per point, modulo the 
	//region shifting constraints.

}

//This is like the getMatch, but with larger combined region image information. 
//Here is what I did that is not good.  I was doing a shifting region match like
//in getMatch, where every local patch belonged to one of the typePCAs, and we
//re-assigned local patches according to whether those patches could be reassigned 
//and whether that would improve the overall image match.  That is shifting regions.
//This has the typePCAs not be the pcas of the clustered local patches, but of larger
//scale regions consisting of those patches assigned to a single type in the training.
//The typePCAs then specify the statistics of the pooled regions.  
//This match uses this information.  The idea is to start with the expected patch 
//assignments. Then change the patch assignments according to whether the patch can
//be re-assigned (that's the priors), whether subtracting that patch qf from the
//larger region to which is belongs would improve that region, and whether adding the
//patch qf to newly assigned larger region would improve that region as well.
double RegionShift::getCombinedMatch(double **qfs, double &totModes)
{
	//Let's start by projecting the data.
	//No projections here, as opposed to the original
	//double ** cX = projPCA->projectData(qfs, numPts);

	srand( time(NULL) );

	computeWeights();
	//now weights is filled out.

	double beforeMatch = computeShiftedMatchCombined(initModel, qfs);

	doShiftingCombined(qfs);
	//memcpy(finalModel, initModel, sizeof(int)*numPts);

	double theMatch = computeShiftedMatchCombined(finalModel, qfs);

	printf("Combined RegionShift match: from %lf to %lf.\n", beforeMatch, theMatch);

	totModes += 1.0*numModes;


	return theMatch;
}


//This function fills out a weights array, that will be used for computing mixing
//alphas...
void RegionShift::computeWeights()
{
	//Go through the patches and get the weights.
	for (int i = 0; i < patches.size(); i ++)
	{
		weights[i] = patches[i]->getMiddleWeight();
	}
}

//This is different from the above doShifting. This combines the qfs
//belonging to the same cluster and uses the typePCAs, which in this 
//case are stats for the combined regions.
/*
I'm also adding the combineQFs function, to build the larger region qfs.
This is tough: the regions change when we add or subtract a patch, and that
changes the alphas, which depend on the order of the additions...
I see. When/if a patch is reassigned, the changed larger regions must be
recomputed.
No that's not so. I can subtract and add in any order. So build the larger
qfs once at the beginning. finalModel stores the data needed to build the 
larger qfs
  */

void RegionShift::doShiftingCombined(double ** qfs)
{
	//X here is assumed to be the dim of the typePCAs by numPts.

	int pt, i, n, fmpt;
	double match0, match1, match2, match3;
	
	memcpy(finalModel, initModel, sizeof(int)*numPts);
	
	int numChanged = 1, cur = 0, bestType, totalChanged = 0, temp;
	int size = sc->emSize;
	bool onBorder;

	double bestMatch, alpha, beta, changeScore, priorScore;

	double *wts = new double[numTypes];
	double **lqfs = new double*[numTypes];
	double *tempqf = new double[size];
	double *tempqf1 = new double[size];
	for (i = 0; i < numTypes; i ++)
	{
		lqfs[i] = new double[size];
		memset(lqfs[i], 0, size*sizeof(double));
	}

	double *lqfScores = new double[numTypes];
	memset(lqfScores, 0, numTypes*sizeof(double));
	
	int *possibilities = new int[numTypes];
	double *scores = new double[numTypes];


	//I want to constrain the motion of the smaller regions, but only in the 
	//local shifting case, not the level 0 clustered stuff I'm doing.  I'll just look
	//and see if the number of points is high, then let's constrain the shifting.
	//see RegionShiftCombined in matlab. I need the adjprior 5 x numTypes.
//	cerr << "-001\n";
	double **adjPrior =  new double *[6];
	double lac = 0.0, lap = 0.0;
	double *adjpt = new double[numTypes];
	for (i = 0; i < 6; i ++)
	{
		adjPrior[i] = new double[numTypes];
		memset(adjPrior[i], 0, numTypes*sizeof(double));
	}
//cerr << "-001b\n";
	for (pt = 0; pt < numPts; pt ++)
	{
		adjPrior[1][finalModel[pt]] += 1;
		n = 0;
		while (adj[pt][n] != -1)
		{
			temp = finalModel[adj[pt][n]]; //the neighboring type.
			if (temp == finalModel[pt])	//is us.
			{
				adjPrior[2][finalModel[pt]] += 1;
			}
			n ++;
		}
	}
//cerr << "-001c\n";
	for (i = 0; i < numTypes; i ++)
	{
		adjPrior[3][i] = adjPrior[2][i] / adjPrior[1][i];
		adjPrior[4][i] = (adjPrior[3][i] / 14.0);
		adjPrior[4][i] *= adjPrior[4][i];
		adjPrior[5][i] = 0;
	}

	//1 is number of pts of this type.
	//2 is total adjacency of the type
	//3 is the average adjacency. 4 is the hardcoded variance allowed. 
	//5 is the penalty associated with the current configuration. 
	//1, 2, and 5 will change every time a pt switches allegiance.
	//So now I have the adjPrior as in matlab. It controls the deformation of the 
	//larger region using a gaussian on the average adjacency of the pt's making
	//up the region.


	int *changed = new int[numTypes];
	memset(changed, 1, numTypes*sizeof(int));
	
//	cerr << "-002\n";


	//debug: why is numChanged almost always 6?
	//It was because I was using i in two nested for loops. oops.  
	//But the below is nice as prints the model after each iteration through
	//the surface.
	//	cout << "\n\n";
	//	for (i = 0; i < numPts; i ++)
	//	{
	//		cout << finalModel[i] << " ";
	//		if ((i % 34) == 33)
	//			cout << endl;
	//	}
		//cout << "\n\n";


	while (numChanged > 0 && cur < 10) //force stopping.
	{

		//Redo the weights and make the larger region qfs. 
		//This is like matlab's combineQFs.
		//But only the ones we have to. Looked at the changed.
		for (n = 0; n < numTypes; n ++)
		{
			if (changed[n]) //non-zero means.
			{
				wts[n] = 0.0;
				memset(lqfs[n], 0, size*sizeof(double));
			}
		}

		for (i = 0; i < numPts; i ++)
		{
			fmpt = finalModel[i];
			
			if (changed[fmpt])
			{
				wts[fmpt] += weights[i];
				alpha = weights[i] / wts[fmpt];

				mixQF(lqfs[fmpt], qfs[i], size, size, tempqf, size, alpha);
				memcpy(lqfs[fmpt], tempqf, size*sizeof(double));
			}
			
		}

		for (n = 0; n < numTypes; n++)
		{
			if (changed[n]) //non-zero means.
			{
				lqfScores[n] = typePCAs[n]->computeMahal(lqfs[n], true);
			}
		}
		memset(changed, 0, numTypes*sizeof(int));


		/*
		memset(wts, 0, numTypes*sizeof(double));
		for (n = 0; n < numTypes; n ++)
			memset(lqfs[n], 0, size*sizeof(double));

		for (i = 0; i < numPts; i ++)
		{
			wts[finalModel[i]] += weights[i];
			alpha = weights[i] / wts[finalModel[i]];

			mixQF(lqfs[finalModel[i]], qfs[i], size, size, tempqf, size, alpha);
			memcpy(lqfs[finalModel[i]], tempqf, size*sizeof(double));
		}

		for (n = 0; n < numTypes; n ++)
		{
			lqfScores[n] = typePCAs[n]->computeMahal(lqfs[n], true);
		}
		*/
		//Okay, now we have the combined qfs (lqfs) and wts stores the total weights
		//contributing to the larger qfs.



		numChanged = 0;
		//Get a permutation.
		//Permute();
		
		for (i = 0; i < numPts; i ++)
		{
			//if none of the adjacent verts to this point are different from it, then 
			//don't do anything.
			
			memset(possibilities, 0, numTypes*sizeof(int));
			memset(scores, 0, numTypes*sizeof(double));
			
			pt = permutation[i];

			fmpt = finalModel[pt];
			
			//possibilities[fmpt] = 1;
			//Want to ensure that the point is not by itself, so don't set that here.
			

			//Here we loop through the neighbors, determining which types amongst our
			//neighbors are at once different from us (finalModel[pt]) and possible for 
			//us to be.
			n = 0; 
			onBorder = false;
			memset(adjpt, 0, numTypes*sizeof(double));

			while (adj[pt][n] != -1)
			{
				temp = finalModel[adj[pt][n]]; //the neighboring type.
				if (temp != fmpt)	//is not us.
				{	
					if (priors[pt*numTypes + temp] > 0 ) //but we can be it.
					{
						onBorder = true;
						possibilities[temp] = 1;
					}
				}
				else //neighbor is like us, so we're not alone.
				{
					possibilities[fmpt] = 1;
				}
				
				if (priors[pt*numTypes + temp] > 0)
				{
					//Additional bit to control the motion through the adjPrior, the average adjacency
					//penalty.
					adjpt[temp] += 1; 
				}
				n ++;
			}
			
			if (!onBorder)
				continue;
			
			//if (!possibilities[fmpt])
			//	int qq = 1;
			
			//Now we're going to score the point with respect to the possibilities.

			//match0 = typePCAs[fmpt]->computeMahal(lqfs[fmpt], true);
			match0 = lqfScores[fmpt];
			//the match of the currently assigned region for this point.
			
			//Subtract pt's qf from the current region.
			alpha = weights[pt] / wts[fmpt];
			if (alpha > .99)
			{
				memset(tempqf1,0, size*sizeof(double));
			}
			else
			{
				beta = -alpha / (1.0-alpha);
				mixQF(lqfs[fmpt], qfs[pt], size, size, tempqf1, size, beta);
			}
			
			match3 = typePCAs[fmpt]->computeMahal(tempqf1, true);
			//the match of the currently assigned region if we subtract this pt's qf.

			changeScore = (match3 - match0);// + log(priors[pt*numTypes + fmpt]);

			if (numTypes > 20)
			{
				//Additional bit to control the motion through the adjPrior, the average adjacency
				//penalty. Use the computed adjacency of this point. adjpt.
				lac = (adjPrior[2][fmpt]-2*adjpt[fmpt]) / (adjPrior[1][fmpt] - 1); //new average adj.
				lac = (lac - adjPrior[3][fmpt])*(lac - adjPrior[3][fmpt]) / adjPrior[4][fmpt]; //new penalty.
				lac = lac - adjPrior[5][fmpt]; //diff with current penalty.
			}
			else
			{
				priorScore = (wts[fmpt]-weights[pt])*log(priors[pt*numTypes + fmpt]);
			}
			

			
			for (n = 0; n < numTypes; n ++)
			{
				//Unlike in getMatch above, here we have to ask about the larger regions.
				//We need to know how the match changes if pt is assigned to a new region,
				//which involves knowing the current larger qfs, versus two new larger qfs.
				//one the region currently assigned, but without this qf, the other the 
				//potential region assignment, but with this qf.

				if (possibilities[n])
				{


					//If the possible type is the type I am already, I shouldn't add this qf
					//to it (it's already been added).  
					if (n != finalModel[pt])
					{
						//match1 = typePCAs[n]->computeMahal(lqfs[n], true);
						match1 = lqfScores[n];
						//the current match of the neighboring region.
						
						//Add pt's qf to the neighboring region. 
						alpha = weights[pt] / (weights[pt] + wts[n]);
						mixQF(lqfs[n], qfs[pt], size, size, tempqf, size, alpha);
						match2 = typePCAs[n]->computeMahal(tempqf, true);


						//The two matches of these scenarios are what to consider.  How is match2
						//an improvement over match1 and how is match3 an improvement over match0.

						//scores[n] = (wts[n]*(match2 - match1) + wts[fmpt]*changeScore) / (wts[fmpt] + wts[n]);  
						
						scores[n] = (wts[n]+weights[pt])*(match2 - match1) + (wts[fmpt]-weights[pt])*changeScore;  

						//scores[n] = priors[pt*numTypes + n]*(wts[n]+weights[pt])*(match2 - match1) + priors[pt*numTypes + fmpt]*(wts[fmpt]-weights[pt])*changeScore; 

						//- log(priors[pt*numTypes + n]);
						//because a smaller number is better, if match2 < match1, then pt's qf improves
						//the neighbor region.  If match3 < match0, then excluding pt's qf from the current
						//region improves the current region. Both terms may come out negative (then 
						//a definite improvement) but they may also compete. That's why they're normalized.

						if (numTypes > 20) //Arbitrary check, don't want to do this for clustered shifting, but yes for local shifting.
						{	
							//additional control of the deformation according to average adjacency.
							lap = (adjPrior[2][n]+2*adjpt[n]) / (adjPrior[1][n] + 1); //new average adj.
							lap = (lap - adjPrior[3][n])*(lap - adjPrior[3][n]) / adjPrior[4][n]; //new penalty.
							lap = lap - adjPrior[5][n]; //diff with current penalty.
							scores[n] =  (match2 - match1) + changeScore + lac + lap;
						}else
						{
							//scores[n] += priorScore - (wts[n]+weights[pt])*log(priors[pt*numTypes + n]);
							//if the region-size-weighted prior for the patch is higher as the current type
							//than as the prospective type, this addtional term is positive, otherwise negative.
							scores[n] *= priors[pt*numTypes + n];
						}
					}
					else
					{
						scores[n] = 0.0;//- 2.0 * log(priors[pt*numTypes + n]);
					}
				}
			}
			//We normalize since the decision on whether to change is based on an
			//energy including other terms. These terms should be commensurate, so
			//that the terms can offset one another.
			

			//Now we do our energy function thing.  
			//First stab: find the best score for non-zero prior probability.  The non-zero
			//prior was in the checking agains the neighbors.
			
			bestMatch = 1.0e10; bestType = -1;
			for (n = 0; n < numTypes; n ++)
			{
				if (possibilities[n] && scores[n] < bestMatch)
				{
					bestType = n;
					bestMatch = scores[n];
				}
			}
			
			if (bestType == -1)
			{
				printf("RegionShift: bestType is bad, bestMatch found is %lf. Keeping current assignment %d\n", bestMatch, fmpt); 
				bestType = fmpt;
			}

			if (fmpt == 45 || bestType == 45)
				int qq = 1;

			if (bestType != fmpt) //we're going to change the type of a point.
			{
				
				//need to change the current larger qfs to reflect this change.  maybe a better way, cause here we'll have 
				//	to recompute the larger qfs, but we already did in the match2 and match3, though they would have been
				//	overwritten by other possibilities.  so maybe integrate this part with the above... but then numChanged 
				//	may change more than once at a single point during the optimization...extra qfs for storing the results?
				//	no point.  just redo the computation.  also changes wts.
				alpha = weights[pt] / (weights[pt] + wts[bestType]);
				mixQF(lqfs[bestType], qfs[pt], size, size, tempqf, size, alpha);
				memcpy(lqfs[bestType], tempqf, size*sizeof(double));
				lqfScores[bestType] = typePCAs[bestType]->computeMahal(lqfs[bestType], true);
				
				alpha = weights[pt] / wts[fmpt];
				if (alpha > .99)
				{
					memset(tempqf1, 0, size*sizeof(double));
				}
				else
				{
					beta = -alpha / (1.0-alpha);
					mixQF(lqfs[fmpt], qfs[pt], size, size, tempqf1, size, beta);
				}

				memcpy(lqfs[fmpt], tempqf1, size*sizeof(double));
				lqfScores[fmpt] = typePCAs[fmpt]->computeMahal(lqfs[fmpt], true);

				wts[fmpt] -= weights[pt];
				wts[bestType] += weights[pt];

			//	cout << "winning score " << scores[bestType] << " against " << scores[finalModel[pt]] << endl;


				finalModel[pt] = bestType;
				numChanged ++;
				changed[fmpt] += 1;
				changed[bestType] += 1;
				totalChanged ++;

				//cout << "\n" << pt << "\t" << fmpt << "\t" << bestType;


				//also change the adjPrior stuff to be right, changing 1,2,5 as described above.
				lac = (adjPrior[2][fmpt]-2*adjpt[fmpt])/(adjPrior[1][fmpt]-1);
				adjPrior[5][fmpt] = (lac - adjPrior[3][fmpt])*(lac - adjPrior[3][fmpt]) / adjPrior[4][fmpt];
				adjPrior[1][fmpt] = adjPrior[1][fmpt]-1;
				adjPrior[2][fmpt] = adjPrior[2][fmpt] - 2*adjpt[fmpt];
				
				lap = (adjPrior[2][bestType]+2*adjpt[bestType])/(adjPrior[1][bestType]+1);
				adjPrior[5][bestType] = (lap - adjPrior[3][bestType])*(lac - adjPrior[3][bestType]) / adjPrior[4][bestType];
				adjPrior[1][bestType] = adjPrior[1][bestType] + 1;
				adjPrior[2][bestType] = adjPrior[2][bestType] + 2*adjpt[bestType];




				//I'm going to redo the lqfs more often:
				for (n = 0; n < numTypes; n ++)
				{
					if (changed[n] > 10) 
					{
						wts[n] = 0.0;
						memset(lqfs[n], 0, size*sizeof(double));


						for (int iter = 0; iter < numPts; iter ++)
						{
							if (finalModel[iter] == n)
							{
								wts[n] += weights[iter];
								alpha = weights[iter] / wts[n];
								
								mixQF(lqfs[n], qfs[iter], size, size, tempqf, size, alpha);
								memcpy(lqfs[n], tempqf, size*sizeof(double));
							}
							
						}

						lqfScores[n] = typePCAs[n]->computeMahal(lqfs[n], true);
						changed[n] = 1;
					}
				}



			}//if bestType equals finalModel[pt] then we don't need to change the larger qfs or the wts.
		}
		
		cur ++;
		cout << " " << numChanged << " ";


		//debug: why is numChanged almost always 6?
		//cout << "\n\n";
		//for (i = 0; i < numPts; i ++)
		//{
		//	cout << finalModel[i] << " ";
		//	if ((i % 34) == 33)
		//		cout << endl;
		//}
		//cout << "\n\n";


	}

	//Okay, done.  So now finalModel contains the best type per point, modulo the 
	//region shifting constraints.

	//cerr << "001\n";
	delete [] possibilities;
	delete [] scores;
	delete [] wts;
	delete [] tempqf;
	delete [] tempqf1;
	for (i = 0; i < numTypes; i ++)
		delete [] lqfs[i];
	delete [] lqfs;
	//cerr << "002\n";
	for (i = 0; i < 6; i ++)
		delete [] adjPrior[i];
	delete [] adjPrior;
	//cerr << "003\n";
	delete [] lqfScores;
	delete [] changed;
	//cerr << "004\n";

}

void RegionShift::Permute(int n)
{
	int temp, r1, r2, i;

	//for (i = 0; i < numPts; i ++)
	//{
	//	if (permutation[i] < 0 || permutation[i] >= numPts)
	//	{
	//		printf("ShiftingRegion error: permutation contains bad value %d at %d", permutation[i], i);
	//	}
	//}

	//Swap random entries n times.
	for (i = 0; i < n; i ++)
	{
		r1 = (int) numPts*(1.0*rand())/RAND_MAX;
		r2 = (int) numPts*(1.0*rand())/RAND_MAX;

		if (r1 >= numPts)
			r1 = numPts-1;
		if (r2 >= numPts)
			r2 = numPts-1;

		temp = permutation[r1];
		permutation[r1] = permutation[r2];
		permutation[r2] = temp;
	}
}

double RegionShift::computeShiftedMatch(const int *typePerPoint, double ** cX, bool normalize)
{
	double match = 0.0;

	//First compute the mahal scores over all points, then normalize with the modes.
	int i;
	numModes = 0;

	for (i = 0; i < numPts; i ++)
	{
		match += typePCAs[typePerPoint[i]]->computeMahal(cX[i], false);
		numModes += typePCAs[typePerPoint[i]]->getNEigs();
	}

	if (normalize)
		match = (match - numModes) / sqrt(2.0 * numModes);


	return match;
}

double RegionShift::computeShiftedMatchCombined(const int *typePerPoint, double ** qfs, bool normalize)
{
	double match = 0.0;

	//We build up all the larger region qfs, then compute the mahals over the larger regions.
	int i, n, size = sc->emSize;
	numModes = 0;

	double *wts = new double[numTypes];
	double **lqfs = new double*[numTypes];
	double *tempqf = new double[size];
	double alpha;

	double * sumpriors = new double[numTypes];
	int *totalpriors = new int[numTypes];

	for (n = 0; n < numTypes; n ++)
	{
		lqfs[n] = new double[size];
		memset(lqfs[n], 0, size*sizeof(double));
		wts[n] = 0.0;

		sumpriors[n] = 0.0;
		totalpriors[n] = 0;
	}

	//Build the larger scale quantile functions.
	for (i = 0; i < numPts; i ++)
	{
		wts[typePerPoint[i]] += weights[i];
		alpha = weights[i] / wts[typePerPoint[i]];

		//if (typePerPoint[i] == 16)
		//	int qqq = 0;
		
		mixQF(lqfs[typePerPoint[i]], qfs[i], size, size, tempqf, size, alpha);
		memcpy(lqfs[typePerPoint[i]], tempqf, size*sizeof(double));

		//if (priors[i*numTypes + typePerPoint[i]] != 1.0)
		//{
			sumpriors[typePerPoint[i]] += priors[i*numTypes + typePerPoint[i]];
			totalpriors[typePerPoint[i]] += 1;
		//}
	}

	//Compute the larger regions' matches.
	for (n = 0; n < numTypes; n ++)
	{
		double tempmatch = typePCAs[n]->computeMahal(lqfs[n], false);

		//I want to penalize models that are far from what the prior prescribes.  I had this problem
		//where a bladder had been init'd as too small, so bladder intensity, which looks like intestine
		//and therefore shows up on the outside occasionally, was on the outside, and the model, with a
		//non-zero prior for this type, said okay. But if we looked to see that it was dark type .92 of 
		//the time, and intestine .08 of the time, maybe we wouldn't be stuck with a delineation that was
		//so unlikely.  Now the two choices are to add the prior into the region shifting itself, or to 
		//add it to a place like here, to multiply the overall match penalty by a factor >= 1, depending 
		//on the prior at the points. this is above as well, when I compute the scores in the region
		//shift.

		//if (totalpriors[n] != 0) //if it is zero, the mahal is really bad.
		//	match += tempmatch/(sumpriors[n]/totalpriors[n]);
		//else
			match += tempmatch;


		numModes += typePCAs[n]->getNEigs() + 1; //1 for the residue.
	}

	if (normalize)
		match = (match - numModes) / sqrt(2.0 * numModes);

	delete [] wts;
	for (i = 0; i < numTypes; i ++)
		delete [] lqfs[i];
	delete [] lqfs;
	delete [] tempqf;

	delete [] sumpriors;
	delete [] totalpriors;


	return match;
}

//Gets the gas bone match for the larger scale regions. One shouldn't believe the level 2 (very local)
//gas bone stuffs.  So this takes the weighted average of the l2 patch gas bones and uses the stats
//given by training, which were stored in typePCA, see readConfig in MultiPatch.
double RegionShift::getCombinedGasBoneMatch(double &totModes)
{
	double theMatch = 0;
	int n, p;

	double ** wts = new double*[numTypes];
	for (n = 0; n < numTypes; n ++)
	{
		wts[n] = new double[3];
		memset(wts[n], 0, 3*sizeof(double));
	}

	//We need to loop over all points, adding their contributions to the relevant larger scale regions.
	double *w;
	for (p = 0; p < numPts; p ++)
	{
		w = patches[p]->getWeights();
		wts[finalModel[p]][0] += w[0];
		wts[finalModel[p]][1] += w[1];
		wts[finalModel[p]][2] += w[2];
	}

	//So now wts is filled out with the small patch low, high and total weights. now compute the matches against
	//the typePCAs.


	for (n = 0; n < numTypes; n++)
	{
		if (wts[n][2] > 0)
		{
			theMatch += typePCAs[n]->computeGasMatch(wts[n], false);
			theMatch += typePCAs[n]->computeBoneMatch(wts[n], false);
			
			totModes += typePCAs[n]->getGasBoneMatters();
		}
	}

	for (n = 0; n < numTypes; n ++)
		delete [] wts[n];
	delete [] wts;


	return theMatch;
}

/********************************************************************************************************************************************************/



/********************************************************************************************************************************************************
Code for the MultiPatch class,.
*********************************************************************************************************************************************************/




/********************************************************************************************************************************************************/


MultiPatch::MultiPatch(const char * filename)
{
	readSamplingConfiguration(filename);

	pList = NULL;

	patchImage = NULL;

	randColors = NULL;

	loggingFileName = 0;
}


MultiPatch::MultiPatch(const M3DObject * object, Image3D * refImage, const char * filename)
{
	
	readSamplingConfiguration(filename);

	pList = NULL;

	patchImage = NULL;

	randColors = NULL;

	initializePatches(object, refImage);

	loggingFileName = 0;
}

void MultiPatch::readSamplingConfiguration(const char * filename)
{
	try {
		reg.readFromFile(filename);
	}
	catch (RException excp) {
		excp.print(std::cout);
		cout << "ERROR: MultiPatch file, " << filename << ", not found !!OR!! ERROR loading:\n";
		excp.print(cout);
	}
	

	//save the name of the sampling config file.
	strcpy(scFilename, filename);


	//create sampling config for inside and outside of each figure, and initialize (and instantiate) patches too.
	//samplingConfigs should be numFigures x 2, where each figure has an inside and outside patch sampling scheme.

	numFigures = reg.getIntValue("config.numFigures", -1);
	//This needs to be defined in the script, if we're going to call it before we have an object.

	samplingConfigs = new SamplingConfiguration[2*numFigures];


	char startStr[1024];
	char tmpStr[1024];
	char scndStr[1024];
	//std::vector< SamplingConfiguration> scs;
	int i, j, p, delinLength, gmLength, bmLength, side;

	//samplingConfigs = new SamplingConfiguration*[numFigures];


	startStr[0] = tmpStr[0] = scndStr[0] = '\0';
	strcat(startStr, "config.figure[%d]");

	for (i = 0; i < numFigures; i ++)
	{
		//first, inside, then outside.
		//samplingConfigs[i] = new SamplingConfiguration[2];

		for (side = 0; side <= 1; side ++)
		{
			sprintf(tmpStr, startStr, i);
			if (side == 0)
				sprintf(tmpStr, "%s.inside", tmpStr);
			else
				sprintf(tmpStr, "%s.outside", tmpStr);
			
			//Now sample without having to think about the side again.
			samplingConfigs[2*i+side].delins = NULL;
			
			samplingConfigs[2*i+side].patchType					= intToPatchType(reg.getIntValue("%s.patchType", 0, tmpStr));
			samplingConfigs[2*i+side].numPatches					= reg.getIntValue("%s.numPatches", -1, tmpStr);
			
			if (samplingConfigs[2*i+side].patchType == DEFINED_PATCHES)
			{
				//read in the gasMatters and boneMatters patches.
				if (samplingConfigs[2*i+side].numPatches > 1)
				{
					samplingConfigs[2*i+side].gasMatters				= reg.getDoubleArray("%s.gasMatters", &gmLength, tmpStr);
					samplingConfigs[2*i+side].boneMatters			= reg.getDoubleArray("%s.boneMatters", &bmLength, tmpStr);
				}else
				{
					samplingConfigs[2*i+side].gasMatters			=  new double(reg.getDoubleValue("%s.gasMatters", 1.0, tmpStr));
					samplingConfigs[2*i+side].boneMatters			=  new double(reg.getDoubleValue("%s.boneMatters", 1.0, tmpStr));
				}
				
				
				//read in an array that represents the delineation on the surface.  
				samplingConfigs[2*i+side].delins = new double *[samplingConfigs[2*i+side].numPatches];
				
				for(j = 0; j < samplingConfigs[2*i+side].numPatches; j ++)
				{
					//read in another patch lower left, upper right, should be six element array.
					sprintf(scndStr, "%s.patch[%d]", tmpStr, j);
					samplingConfigs[2*i+side].delins[j]			= reg.getDoubleArray(scndStr, &delinLength);
					if (delinLength != 6)
						cout << "ERROR: MultiPatch : fig " << i << ", inside, delin_patch " << j << " not read properly.\n";
					
					cout << "Defined Patch: ";
					for (p = 0; p < 6; p ++)
					{
						cout << samplingConfigs[2*i+side].delins[j][p] << " ";
					}
					cout << endl;
				}
			}
			else //not defined patches, so just one gasMatters and boneMatters value.
			{
				samplingConfigs[2*i+side].gasMatters			=  new double(reg.getDoubleValue("%s.gasMatters", 1.0, tmpStr));
				samplingConfigs[2*i+side].boneMatters			=  new double(reg.getDoubleValue("%s.boneMatters", 1.0, tmpStr));
			}
			
			
			samplingConfigs[2*i+side].samplingDistance				= reg.getDoubleValue("%s.samplingDistance", 1.0, tmpStr);
			samplingConfigs[2*i+side].distanceSTD					= reg.getDoubleValue("%s.distanceSTD", 0.5, tmpStr);
			samplingConfigs[2*i+side].isRRelative					= (bool) reg.getBooleanValue("%s.isRRelative", false, tmpStr);
			
			samplingConfigs[2*i+side].numSteps						= reg.getIntValue("%s.numSteps", 6, tmpStr);
			samplingConfigs[2*i+side].localPatchSurfaceLevel		= reg.getIntValue("%s.localPatchSurfaceLevel", 0, tmpStr);
			samplingConfigs[2*i+side].samplingProfileSurfaceLevel	= reg.getIntValue("%s.samplingProfileSurfaceLevel", 4, tmpStr);
			
			samplingConfigs[2*i+side].localPatchRadius				= reg.getDoubleValue("%s.localPatchRadius", 1.5, tmpStr);
			samplingConfigs[2*i+side].localPatchSTD					= reg.getDoubleValue("%s.localPatchSTD", 0.0, tmpStr);
			
			samplingConfigs[2*i+side].lowCutOff						= reg.getDoubleValue("%s.lowCutOff", 800.0, tmpStr);
			samplingConfigs[2*i+side].highCutOff					= reg.getDoubleValue("%s.highCutOff", 1200.0, tmpStr);
			
			samplingConfigs[2*i+side].emSize						= reg.getIntValue("%s.emSize", 100, tmpStr);

			if (side == 0)
				samplingConfigs[2*i+side].isInside						= true;
			else
				samplingConfigs[2*i+side].isInside						= false;
			

		}//side

		
	}//figure

	//if the file has stats, then populate the patchPCAs vector.
	//in saveMP you see that numPatches is set correctly, so when then is run at target time, 
	//the keys will be correct in the config.


	if(reg.hasKey("stats"))
	{
		std::vector< PatchPCA *> tempPatchPCAs; //to be pushed onto the patchPCAs data member.
		int pcaSize, weightSize, otherSize;

		startStr[0] = tmpStr[0] = scndStr[0] = '\0';
		strcat(startStr, "stats.figure[%d]");

		//patchPCAs.assign(2*numFigures, std::vector< PatchPCA *> );
		//regionShifts.assign(2*numFigures, NULL);
		patchPCAs.resize(2*numFigures);
		regionShifts.resize(2*numFigures);

		for (i = 0; i < numFigures; i ++)
		{
			//In the default stats, there is a pca for each patch.  In match_mle, the pcas are
			//listed up front (as matchType's) and numPatches means the number of them.  Then,
			//there is not a list of patch[i], since there is no particular appearance model at
			//at patch.


			//first inside, then out.

			for (side = 0; side <= 1; side ++)
			{
				sprintf(tmpStr, startStr, i);

				if (side == 0)
					sprintf(tmpStr, "%s.inside", tmpStr);
				else
					sprintf(tmpStr, "%s.outside", tmpStr);

				//Now tmpStr gets to the side of the figure: stats.figure[1].outside

				
				samplingConfigs[2*i+side].statType = intToStatType(reg.getIntValue("%s.statType", 0, tmpStr));
				samplingConfigs[2*i+side].numPatches = reg.getIntValue("%s.numPatches", 0, tmpStr);
				
				if (samplingConfigs[2*i+side].statType == MATCH_DEFAULT || samplingConfigs[2*i+side].statType == MATCH_MLE)
				{
					
					
					if (samplingConfigs[2*i+side].statType == MATCH_DEFAULT)
						sprintf(scndStr, "%s.patch", tmpStr);
					else if (samplingConfigs[2*i+side].statType == MATCH_MLE)
						sprintf(scndStr, "%s.matchType", tmpStr);
					
					
					for (p = 0; p < samplingConfigs[2*i+side].numPatches; p ++)
					{
						PatchPCA *pca = new PatchPCA();
						tempPatchPCAs.push_back(pca); //put it in the list.
						
						double * d = reg.getDoubleArray("%s[%d].pca", &pcaSize, scndStr, p);
						pca->setPca(d, samplingConfigs[2*i+side].emSize, pcaSize);
						
						d = reg.getDoubleArray("%s[%d].weights", &weightSize, scndStr, p);
						pca->setWeights(d, weightSize);
						
						if (samplingConfigs[2*i+side].patchType == DEFINED_PATCHES && samplingConfigs[2*i+side].numPatches > 1)
							pca->setGasBoneMatters(samplingConfigs[2*i+side].gasMatters[p], samplingConfigs[2*i+side].boneMatters[p]);
						else
							pca->setGasBoneMatters(*samplingConfigs[2*i+side].gasMatters, *samplingConfigs[2*i+side].boneMatters);
						
						//that's all for this patch pca .
					}//patch
					
					
					
					//patchPCAs.push_back(tempPatchPCAs);
					//push back doesn't work if the statType is shifting for any side of any figure,
					//in that you'll push the interior on, but then the next push is another interior.
					patchPCAs.insert(patchPCAs.begin() + 2 * i + side, tempPatchPCAs);
					
					//now clear the tempPatches vector.
					tempPatchPCAs.clear();
				}
				else if (samplingConfigs[2*i+side].statType == MATCH_SHIFTING)
				{
					
					PatchPCA *projpca = new PatchPCA();

					double * d = reg.getDoubleArray("%s.projPCA", &pcaSize, tmpStr);
					if (d != NULL)
						projpca->setPca(d, samplingConfigs[2*i+side].emSize, pcaSize);	


					int nclust = reg.getIntValue("%s.numClust", 0, tmpStr);

					for (p = 0; p < nclust; p ++)
					{
						PatchPCA *pca = new PatchPCA();
						tempPatchPCAs.push_back(pca);

						d = reg.getDoubleArray("%s.typePCA[%d]", &pcaSize, tmpStr, p);
						pca->setPca(d, samplingConfigs[2*i+side].emSize, pcaSize);

						if (nclust > 0) //pick up the weights as well.
						{
							d = reg.getDoubleArray("%s.typeWeights[%d]", &weightSize, tmpStr, p);
							pca->setWeights(d, weightSize);
							
							if (samplingConfigs[2*i+side].patchType == DEFINED_PATCHES && samplingConfigs[2*i+side].numPatches > 1)
								pca->setGasBoneMatters(samplingConfigs[2*i+side].gasMatters[p], samplingConfigs[2*i+side].boneMatters[p]);
							else
								pca->setGasBoneMatters(*samplingConfigs[2*i+side].gasMatters, *samplingConfigs[2*i+side].boneMatters);
						}
					}


					//Read in the prior, and construct the region shifting object.
					double *pr = reg.getDoubleArray("%s.priors", &otherSize, tmpStr);

					RegionShift * rS = new RegionShift(projpca, tempPatchPCAs, pr, otherSize, 
						&samplingConfigs[2*i+side]);

					regionShifts.insert(regionShifts.begin() + 2 * i + side, rS);

					//now clear the tempPatches vector.
					tempPatchPCAs.clear();

					//Now get gas and bone from the individual patches.
					sprintf(scndStr, "%s.patch", tmpStr);

					for (p = 0; p < samplingConfigs[2*i+side].numPatches; p ++)
					{
						PatchPCA *pca = new PatchPCA();
						tempPatchPCAs.push_back(pca); //put it in the list.
					
						d = reg.getDoubleArray("%s[%d].weights", &weightSize, scndStr, p);
						pca->setWeights(d, weightSize);
						
						if (samplingConfigs[2*i+side].patchType == DEFINED_PATCHES && samplingConfigs[2*i+side].numPatches > 1)
							pca->setGasBoneMatters(samplingConfigs[2*i+side].gasMatters[p], samplingConfigs[2*i+side].boneMatters[p]);
						else
							pca->setGasBoneMatters(*samplingConfigs[2*i+side].gasMatters, *samplingConfigs[2*i+side].boneMatters);
						
						//that's all for this patch pca .
					}//patch

					patchPCAs.insert(patchPCAs.begin() + 2 * i + side, tempPatchPCAs);

					tempPatchPCAs.clear();

				}

			}//side

			
		}//figure

		/*int my0 = patchPCAs[0].size();
		int my1 = patchPCAs[1].size();
		int my2 = patchPCAs[2].size();
		int my3 = patchPCAs[3].size();
		int my4 = patchPCAs[4].size();
		int my5 = patchPCAs[5].size();
		*/

	}//if stats
}

MultiPatch::~MultiPatch()
{
	int i, j, f;

	for (i = 0; i < patches.size(); i ++)
	{
		for (j = 0; j < patches[i].size(); j ++)
		{
			delete patches[i].at(j);
		}
	}
	patches.clear();


	for (i = 0; i < patchPCAs.size(); i ++)
	{
		for (j = 0; j < patchPCAs[i].size(); j ++)
		{
			delete patchPCAs[i].at(j);
		}
	}
	patchPCAs.clear();


	delete [] samplingConfigs;


	if(pList != NULL){
		for(f = 0; f < numFigures; f++){
			if(pList[f] != NULL)
				delete pList[f];
		}
		delete [] pList;
	}

	delete patchImage;

	if (randColors != NULL)
	{
		for (f = 0; f < numFigures; f++)
		{
			if (randColors[f] != NULL)
				delete randColors[f];
		}
		delete [] randColors;
	}
}


void MultiPatch::initializePatches(const M3DObject * object, Image3D * refImage)
{
	int i; 

	numFigures = object->getFigureCount();

	for (i = 0; i < numFigures; i ++)
	{
		initializeFigurePatches(object, refImage, i, 0);  //inside.
		initializeFigurePatches(object, refImage, i, 1);  //outside.
	}
}


void MultiPatch::updateTauAndWeights(int figureID, int whichSide, Image3D *refImage)
{
	//access samplingConfig to get the sampling info.
	
	SamplingConfiguration *sc = &samplingConfigs[2 * figureID + whichSide];
	sc->taus.clear();
	sc->tauWeights.clear();

	//See Patch::sampleYourself and the note I made there.  Basically, when using the binary
	//to sample the regions, the binary is all eli uses, whereas my patch sampling from the
	//boundary out or in never gets the chance to ask for out voxels that happen to be in the
	//binary (for an inside patch) or the other way.  So here I am, forcing the wrong sampling.

	double wrongSamplingMult = 0.5;


	int n = (int) floor((1.0 + wrongSamplingMult) * sc->numSteps);
	double first, second;
	
	if (whichSide == 0)
	{
		first = -1 * sc->samplingDistance;
		second = wrongSamplingMult * sc->samplingDistance;
	}else
	{
		first = -1.0 * wrongSamplingMult * sc->samplingDistance;
		second = sc->samplingDistance;
	}

	//if not r-relative, then mult by worldToModelScale. Again, see 
	//SurfacePatchEnsemble and Eli's code. weights will be the same regardless.

	if (!sc->isRRelative)
	{
		double temp = refImage->getModelToWorldScale();
		first /= temp;
		second /= temp;
		sc->distanceSTD /= temp;
	}


	//So a simple line from first to second in n steps. The weights
	//can also be computed here.
	double thisWeight;
	double step = (second - first) / (n - 1);
	double sqStd = sc->distanceSTD;
	sqStd *= sqStd;

	sc->taus.push_back(first);
	if (sqStd > 0)
		sc->tauWeights.push_back(exp( -0.5 * (sc->taus[0] * sc->taus[0]) / sqStd));
	else
		sc->tauWeights.push_back(1.0);


	int i;
	for (i = 1; i < n; i ++)
	{
		//sc->taus.push_back(sc->taus[i-1] + step);
		sc->taus.push_back(first + i*step);
		thisWeight = 1.0;
		if (sqStd > 0)
			thisWeight = exp( -0.5 * (sc->taus[i] * sc->taus[i]) / sqStd);
		sc->tauWeights.push_back(thisWeight);
	}

	for (i = 0; i < n; i ++)
	{
		cout << sc->taus[i] << " ";
	}
	cout << endl << endl;
	for (i = 0; i < n; i ++)
	{
		cout << sc->tauWeights[i] << " ";
	}
	cout << endl << endl;

}

void MultiPatch::initializeFigurePatches(const M3DObject *object, Image3D *refImage, 
										 int figureID, int whichSide)
{
	Bpoint *Bpoints;
	Bpoint *Bpoints2;
	int p, j;
	Vector3D p1, p2;
	double tempd;
	int offset = 0;
	ProfileComponent pc;

	std::vector< Patch *> tempPatches;

	updateTauAndWeights(figureID, whichSide, refImage);

	SamplingConfiguration *sc = &samplingConfigs[2 * figureID + whichSide];

	//for match_shifting we have to build the adjacency as well.
	if (sc->statType == MATCH_SHIFTING)
		regionShifts[2 * figureID + whichSide]->buildAdjacency(object, figureID, sc->localPatchSurfaceLevel);

	
	double localPatchRadiusSquared = sc->localPatchRadius / 
										refImage->getModelToWorldScale();
	localPatchRadiusSquared *= localPatchRadiusSquared;

	double mylocalPatchSTD = sc->localPatchSTD / refImage->getModelToWorldScale();
	
	int numFigPatches = 0;
	int numProfiles;
	

	//-> sets numFigPatches: if delin patches then this number is wrong, is really the number of 
	// spokes returned for that surface level.
	// so change 
	getBPoints(object, figureID, whichSide, sc->localPatchSurfaceLevel,
		       numFigPatches, Bpoints, true);
	



	//if sc.delins == DEFINED_PATCHES, then these points above are spokeID's and anchorPoints, while
	//the profiles will be all the higher level points. 
	//OOH!.  I should just make the local patches with the voxel level profiles.  Then, if I'm doing
	//defined patches, combine them using the addPatch functions.

	//So, first make the local patches, with one anchorpoint and spokeID per patch.

	getBPoints(object, figureID, whichSide, sc->samplingProfileSurfaceLevel,
				numProfiles, Bpoints2, false);


	const M3DFigure * figure = object->getFigurePtr(figureID);

	/*double *profileWeights = new double[numProfiles];
	int *profilePatch = new int[numProfiles];
	//For debugging, I just want to see if the pattern of profile weights is wavy...
	for (j = 0; j < numProfiles; j ++)
	{
		profileWeights[j] = .0;
		profilePatch[j] = 0;
	}

	double *whichPatch = new double[numFigPatches];
	for (p = 0; p < numFigPatches; p ++)
		whichPatch[p] = (.25 * rand())/RAND_MAX;
	*/
	//The above code is for the different colored patches.



	//See the most inside loop below for more info.  I just wanted to not do the division in every 
	//loop.
	double littleSTD = mylocalPatchSTD / 5;  //just enough to get the surounding areas of a spoke.
	double littleSTDsq = littleSTD*littleSTD;

	//littleSTD = 0;



	for(p = 0; p < numFigPatches; p++)
	{
		//int(tuningWt(IgnoreMatchAtEnds)) == 1 in place of true below.

		if( typeid(*figure) == typeid(M3DTubeFigure) &&
#ifndef BINARY
			int(tuningWt(IgnoreMatchAtEnds)) == 1 &&
#endif
			(Bpoints[p].u == 0 || Bpoints[p].u == figure->getPrimitiveCount()-1 ) )
			// uncomment to inject bug to match what was in the spe code
			//(Bpoints[p].u == 0 || Bpoints[p].u == figure->getPrimitiveCount()-0 ) )
		{
			// Tube figure and patch corresponds to spokes from an end atom. Don't use them.
			continue;
		}
		
		//Patch * patch = new Patch(figureID, p, Bpoints[p], &sc);
		//patch->setSC(&sc);

		tempPatches.push_back(new Patch(figureID, p, Bpoints[p], sc));
		//patches.push_back(patch);

		//Now find near profiles to each one
		p1.set(Bpoints[p].pnt);
		for (j = 0; j < numProfiles; j++) {
			p2.set(Bpoints2[j].pnt);
			tempd = p1.distSquare(p2);
			if (tempd <= localPatchRadiusSquared) {
				if (fabs(Bpoints[p].t) < 1 || Bpoints[p].t * Bpoints2[j].t > 0){
					pc.index = j;
					pc.weight = 1;
					//this bit is rough here, where we're trying to assign a weight to this profile
					//given the patch we want to add it to.  I was getting little rings of low weight 
					//near patch center in defined regions, since other spokes weren't close enough by
					//to bring the weight of this profile to 1 (see addPatch above where we add profile weights.)
					//So, here I want a little more flat function than one that falls off like a Gaussian.
					//A problem with the below.  this 
					if(mylocalPatchSTD > 0.0)
					{
						if (tempd > littleSTDsq)
						{
							double tempsqrt = sqrt(tempd);
							pc.weight = exp( -0.5 * (tempsqrt - littleSTD)*(tempsqrt - littleSTD) / (mylocalPatchSTD * mylocalPatchSTD));
						}
					}
					//patch->addProfile(pc);	
					//Debug on profile weights;
					//if(pc.weight > profileWeights[j])
					//{
					//	profileWeights[j] = pc.weight;
					//	profilePatch[j] = p;
					//}


					tempPatches.back()->addProfile(pc);
				}
			}
		}
	}
	



/*All this code below is for different colored patches.
	//int bladder[] = {30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 60, 61, 62, 63, 64, 65, 77};
		int bladder[] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76};
		int rectum[] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 25, 26, 27, 28, 29, 30, 31, 33, 34, 35, 36, 37, 38, 39, 41, 42, 43, 44, 45, 46, 47, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138};


	//int rectum[] = {39, 47, 55, 63, 71, 79, 87};

	if (figureID == 0)
	{
		for (j = 0; j < numProfiles; j ++)
		{
			for (p = 0; p < 61; p ++)
			{
				if (profilePatch[j] == bladder[p])
					profilePatch[j] = -1;
			}
		}
	}
	else if (figureID == 2)
	{
		for (j = 0; j < numProfiles; j ++)
		{
			for (p = 0; p < 135; p ++)
			{
				if (profilePatch[j] == rectum[p])
					profilePatch[j] = -1;
			}
		}
	}
	
	char mypatchname[1000];
	sprintf(mypatchname, "patchPerProfile_%d.txt", figureID);
	FILE *fout = fopen(mypatchname, "w");

	for (j = 0; j < numProfiles; j ++)
	{
		if (profilePatch[j] >= 0)
			fprintf(fout, "%lf ", whichPatch[profilePatch[j]]);
		else
			fprintf(fout, "%lf ", 0.5);
	}

	fclose(fout);

*/




	//Debug with profile weights.  Just output the file.  Actually I do this below in the 
	//delin patches case.  oops.



	//Delete the vector of patches (patch pointers actually) that is 
	//at the position where this list needs to be.
	//This code is pretty much untested, as patches is usually being constructed 
	//for the first time here.
	if (patches.size() > 2 * figureID + whichSide)
	{
		for (j = 0; j < patches[2 * figureID + whichSide].size(); j ++)
		{
			delete patches[2 * figureID + whichSide].at(j);
		}
		patches.erase(patches.begin() + 2 * figureID + whichSide);
	}
	
	//Now, we either insert the tempPatches vector we just made, or combine those with the 
	//the delineations defined in the sc.  Also, in here, once we have the patches, attach
	//their PatchPCA elements.

	if (sc->patchType == LOCAL_PATCHES)
	{
		//First, attach the PatchPCA elements.  Since all of this is done by value of a pointer,
		//I could wait til after the insertion, but that's less intuitive.
		if(reg.hasKey("stats"))
		{
			if (sc->statType == MATCH_DEFAULT || sc->statType == MATCH_SHIFTING)
			{
				for (j = 0; j < tempPatches.size(); j ++)
				{
					tempPatches[j]->setPCA(patchPCAs[2 * figureID + whichSide][j]);
				}
			}
			else if (sc->statType == MATCH_MLE)
			{
				//Do nothing, since we don't need a pca per patch.  They were all listed up front.
				tempPatches[j]->setPCA(NULL);
			}
		}

		patches.insert(patches.begin() + 2 * figureID + whichSide, tempPatches);

		if (sc->statType == MATCH_SHIFTING)
			regionShifts[2 * figureID + whichSide]->setPatches(tempPatches);


		/*
		if (whichSide)
		{
			FILE * fout;
			char fileName[1024];
			for (p = 0; p < numFigPatches; p ++)
			{
				sprintf(fileName, "profileWeightForPatch%d_Fig%d.txt", p, figureID);
				fout = fopen(fileName, "w");
				
				for (j = 0; j < numProfiles; j++)
				{
					fprintf(fout, "%lf ", tempPatches[p]->getProfileWeight(j));
				}
				
				
				fclose(fout);
			}
		}

		int ww;
		ww = 0;
		*/
		

	}
	else if (sc->patchType == DEFINED_PATCHES)
	{


		//Now I should reconstitute the tempPatches if I need to (if the sampling config requires it,
		//as in a rectum (just an out of the blue example).
		int curDelinCount;
		int *matching = new int[sc->numPatches];
		
		//the u,v,t of lower left and the uvt of upper right are the six items of any sc.delims[i] vector.
		
		//build an additional vector of patch pointers
		std::vector< Patch * > delinPatches;
		for (p = 0; p < sc->numPatches; p ++)
		{
			delinPatches.push_back(new Patch());
			delinPatches[p]->setFigureID(figureID);
			delinPatches[p]->setSC(sc);
			if(reg.hasKey("stats"))
			{
				if (sc->statType == MATCH_DEFAULT)
					delinPatches[p]->setPCA(patchPCAs[2 * figureID + whichSide][p]);
			}
		}
	
		bool isTube = (bool) (0 != (typeid(*figure) == typeid(M3DTubeFigure)));

		for (p = 0; p < tempPatches.size(); p ++)
		{
			//Determine the delineations of the current Patch, still referenced by Bpoints
			//if (p == 59)
			//	int qq = 0;

			//if we're dealing with tubes the below bp may not be the same as Bpoints[p], if the
			//ends are not included.
			Bpoint bp = tempPatches[p]->getFirstAnchorPoint();
			curDelinCount = whichDelin(bp, sc, matching, isTube);

			if (curDelinCount == 0)
			{
				cout << "WARNING: MultiPatch Delineations, patch " << tempPatches[p]->getFirstSpokeID() << " not assigned for figure " << figureID 
					 << ". Bpoint (" << bp.u << ", " << bp.v << " " << bp.t << ")\n";

			}

			while (curDelinCount > 0)
			{
				delinPatches[matching[curDelinCount-1]]->addPatch(*tempPatches[p]);
				curDelinCount --;
			}
			//Now we've added this local patch to all the defined patches it belongs to.

			//Should destroy the patches in tempPatches.
			delete tempPatches[p];
		}


		patches.insert(patches.begin() + 2 * figureID + whichSide, delinPatches);


		//for debug, just write a file with the profile weights of the delin patches.
		/*FILE * fout;
		char fileName[1024];
		for (p = 0; p < sc->numPatches; p ++)
		{
			sprintf(fileName, "profileWeightForPatch%d.txt", 2*figureID + p);
			fout = fopen(fileName, "w");

			for (j = 0; j < numProfiles; j++)
			{
				fprintf(fout, "%lf ", delinPatches[p]->getProfileWeight(j));
			}


			fclose(fout);
		}*/



		tempPatches.clear();
		//Delineated patches finished;		

	}
	else
	{
		//Flexible regions not yet defined.

	}
	
	

	delete [] Bpoints;








/*

	////////////////////////////////////////////////////////////////////////////////
	//atomIDs - This should be the only non-tube compatible code
	//Need list of atomIDs given a list of level 0 boundary points
	const M3DFigure * figure = object->getFigurePtr(figureID);
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
		
		if( typeid(*figure) == typeid(M3DTubeFigure) &&
			int(tuningWt(IgnoreMatchAtEnds)) == 1 &&
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
	*/
}

bool MultiPatch::getBPoints(const M3DObject *object, int figureID, int whichSide,
							int inSurfaceLevel, int &numPoints, Bpoint * &Bpoints, bool copy)
{

	int f;
	//int inSurfaceLevel = samplingConfigs[2*figureID + whichSide].samplingProfileSurfaceLevel;

	M3DFigure * figure = object->getFigurePtr(figureID);
	Xferlist * xferList = convertM3DtoXfer(figure);

	if(pList == NULL){
		pList = new ThallCode::Pointlist_server2 *[2 * numFigures];
		for(f = 0; f < 2 * numFigures; f++){
			pList[f] = NULL;
		}
	}

	if (pList[2 * figureID + whichSide] == NULL) {
		//Requisite Incantations
		pList[2 * figureID + whichSide] = new ThallCode::Pointlist_server2;
		pList[2 * figureID + whichSide]->init(xferList);
		pList[2 * figureID + whichSide]->ComputeSubdivPointCloud(inSurfaceLevel);
	}
	else
		pList[2 * figureID + whichSide]->UpdateSubdivPointCloud(inSurfaceLevel, xferList);

	if(copy)
		pList[2 * figureID + whichSide]->copySubdivBoundaryInfo(&numPoints, &Bpoints);
	else
		pList[2 * figureID + whichSide]->subdivboundaryinfo(&numPoints, &Bpoints);

	delete [] (xferList->atomlist);
	delete xferList;

	//cout << numPoints << " boundary points used.\n";
	//boundary displacements not supported

	return true;
}


PatchType MultiPatch::intToPatchType(int n)
{
	if (n == 0)
		return LOCAL_PATCHES;

	if (n == 1)
		return DEFINED_PATCHES;


	return FLEXIBLE_PATCHES;
}
	

int MultiPatch::patchTypeToInt(PatchType p)
{
	if (p == LOCAL_PATCHES)
		return 0;

	if (p == DEFINED_PATCHES)
		return 1;
	
	return 2;
}

StatType MultiPatch::intToStatType(int n)
{
	if (n == 0)
		return MATCH_DEFAULT;

	if (n == 1)
		return MATCH_MLE;

	if (n == 2)
		return MATCH_SHIFTING;

	return MATCH_UNKNOWN;
}

int MultiPatch::statTypeToInt(StatType s)
{
	if (s == MATCH_DEFAULT)
		return 0;

	if (s == MATCH_MLE)
		return 1;

	if (s == MATCH_SHIFTING)
		return 2;

	return 3;
}

int MultiPatch::whichDelin(Bpoint bp, SamplingConfiguration *sc, int *matching, bool isTube)
{
	//This code will return the row index of the delineation defined in
	//sc for the Bpoint bp, based on u,v,t and that delins are specified by
	//lower left, upper right in u,v,t. 
	//OLD: Really, returns the first index that
	//matches (since it's up to the user to define the patches well).
	//NEW: It returns all the indices that match.
	//Also, u and v are switched somewhere else in pablo.  if you select the 
	//last atom through the gui, it'll say u,v is 4,5 for example, but down
	//here, that's 5,4.  God damnit. So I switch them here.
	int i = 0, count = 0;

	double u = bp.u, v = bp.v, t = bp.t;

	if (!isTube)
	{
		v = u;
		u = bp.v;
	}
	//The below check is not valid in v for tubes if you want more than all v.  

	while (i < sc->numPatches)
	{
		if (	(u >= sc->delins[i][0] && u <= sc->delins[i][3])
			&&  (v >= sc->delins[i][1] && v <= sc->delins[i][4])
			&&  (t >= sc->delins[i][2] && t <= sc->delins[i][5]))
		{
			matching[count] = i;
			count ++;
		}

		i++;
	}

	if (count == 0)
		return -1;  //no matching defined patch found for this bpoint.
	
	return count;
}


void MultiPatch::sampleYourself(const M3DObject * object, Image3D * image, int f, Image3D *binImage, bool debug)
{
	//for each figure, for each of inside and outside, update Bpoints and ask
	//ask each patch to sample itself.

	Bpoint *Bpoints;
	int numProfiles, j, p, inOut, nump0 = patches[2 * f].size(), nump1 = patches[2 * f + 1].size(); //f

	if(patchImage == NULL)
	{
		patchImage = new Image3D(image->getXDim(), image->getYDim(), image->getZDim());
		patchImage->setSpacingAndOrigin(image->getXSpacing(), image->getYSpacing(), image->getZSpacing());
	}
	else if (  patchImage->getXDim() != image->getXDim()
			|| patchImage->getYDim() != image->getYDim()
			|| patchImage->getZDim() != image->getZDim() )
	{
		delete patchImage;

		patchImage = new Image3D(image->getXDim(), image->getYDim(), image->getZDim());
		patchImage->setSpacingAndOrigin(image->getXSpacing(), image->getYSpacing(), image->getZSpacing());
	}
	
	patchImage->clear();

	//I'm going to create a unique array of grey values.  a problem popped up because if two patches on the same
	//side had the same spokeID[0], then the patchImage would force the second patch to be undersampled, in that
	//that number (spokeID[0] + 1) is used to note voxels that have been sampled.  So I create this set so that no
	//two patches have the same identifier in the patchImage.
	if (randColors == NULL)
	{
		randColors = new GreyValue *[numFigures];
		for (j = 0; j < numFigures; j ++)
			randColors[j] = NULL;
	}
	
	if (randColors[f] == NULL)
	{
		//Then I need to instantiate the array.
		randColors[f] = new GreyValue[nump0 + nump1];

		bool good; 
		GreyValue g;
		
		for (p = 0; p < nump0 + nump1; p ++)
		{
			good = false;
			while (!good)
			{
				good = true;
				g = 3000 * (1.0* rand())/RAND_MAX;
				j = 0;
				while (j < p && randColors[f][j] != g)
					j++;
				
				if (j == p)
					randColors[f][j] = g;
				else
					good = false;
			}
		}

	}	


	//Now, ready to sample.

	//for (f = 0; f < numFigures; f ++)
	//{
	for (inOut = 0; inOut <= 1; inOut++)
	{
		int sL = samplingConfigs[2*f + inOut].samplingProfileSurfaceLevel;
		getBPoints(object, f, inOut, sL, numProfiles, Bpoints, false);

		//Vector3D Pi(Bpoints[0].pnt);
		//image->modelToImageCoordinates(Pi);

		
		for (j = 0; j < patches[2 * f + inOut].size(); j ++)
		{
			patches[2 * f + inOut][j]->sampleYourself(image, Bpoints, patchImage, randColors[f][inOut*nump0 + j], binImage);
			//cout << "\tLater, size is " << patches[2 * f + inOut][j]->getDataSize() << ".\n";
		}
		
	}

	//I've sampled in and out.  now write out the image.

	if (debug)
	{
		RAWImageFile imageFile;
		char fileName[1024];
		
		sprintf(fileName, "patchImage%d.raw3", f);
		
		if (imageFile.write(fileName, *patchImage))
			cout << "Saved Image: " << fileName << "." << endl;
		else
			cout << "Image failed to save: " << fileName << "." << endl;
	}


	//}

}

//STOUGH::This function should be modified to do the FLEXIBLE patchType. That way, to
//the rest of pablo, it's just getting the match, but in here there will be another
//optimization, that of the regional delineation. 

//Scratch that.  I did modify it, now flexibleMatch needs to be filled out.
double MultiPatch::getYourMatch(const M3DObject * object, Image3D * image, int f, Image3D *binImage)
{
	double theMatch = 0.0;

	//First, sample ourselves.
	sampleYourself(object, image, f, binImage);

	//Now get the matches. I want inside and outside to have equal weight, but if there
	//are only two inside and ~70 outside, then that's not fair...
	double *theMatches = new double[2];
	double *total = new double[2];

	total[0] = total[1] = 0;
	theMatches[0] = theMatches[1] = 0.0;

	int inOut, p, totalCount = 0;
    // levy hackery to create a log file of the format
	// FIGURE IN/OUT  PATCH   THE_MATCH   TOTAL_MATCH
	std::ofstream out;
    if (loggingFileName) {
	  out.open(loggingFileName);
	  out << "figure#  IN(0)/OUT(1)  patch#  MPF  #eigModes" << endl;
    }
	
	for (inOut = 0; inOut <= 1; inOut++)
	{
		if (samplingConfigs[2*f + inOut].patchType == FLEXIBLE_PATCHES)
		{
			flexibleMatch(object, image, f, inOut, theMatches[inOut], total[inOut], binImage);
		}
		else
		{
			//In the MATCH_DEFAULT case, each patch has a pca to inquire about, but in
			//the MATCH_MLE case, a patch can be used against any of the pca's for that 
			//side for that figure. So call a new match function to take care of that 
			//case.
			if (samplingConfigs[2*f + inOut].statType == MATCH_MLE)
			{
				MLEMatch(object, image, f, inOut, theMatches[inOut], total[inOut], binImage);
			}
			else if (samplingConfigs[2*f + inOut].statType == MATCH_SHIFTING)
			{
				//build the qfs.
				int nPatches = patches[2 * f + inOut].size(); // should be also samplingConfigs[2*f + inOut].numPatches;
				double **qfs = new double *[nPatches];
				for (p = 0; p < nPatches; p ++)
				{
					qfs[p] = patches[2 * f + inOut][p]->getQF();
				}

				//theMatches[inOut] = regionShifts[2 * f + inOut]->getMatch(qfs, total[inOut]);
				theMatches[inOut] = regionShifts[2 * f + inOut]->getCombinedMatch(qfs, total[inOut]);

				delete [] qfs;

				//Now get gas and bone, which are ignored in regionshifting for now.
				if (regionShifts[2*f + inOut]->getNumTypes() > 0) //must be the local shifting.
				{
					//Ask the region Shift object for the gas bone match, since we want it computed
					//at the larger regional scale.  
					double thisMatch = regionShifts[2 * f + inOut]->getCombinedGasBoneMatch(total[inOut]);
					theMatches[inOut] += thisMatch;
				}
				else
				{
					for (p = 0; p < patches[2 * f + inOut].size(); p ++)
					{
						double thisMatch = patches[2 * f + inOut][p]->getYourGas(false);
						theMatches[inOut] += thisMatch;
						//false means do not normalize.
						thisMatch = patches[2 * f + inOut][p]->getYourBone(false);
						theMatches[inOut] += thisMatch;
						
						
						total[inOut]  += patches[2 * f + inOut][p]->getGasBoneMatters();
					}
				}



			}
			else //just do the default.
			{
				for (p = 0; p < patches[2 * f + inOut].size(); p ++)
				{
					double thisMatch = patches[2 * f + inOut][p]->getYourMatch(false);
					theMatches[inOut] += thisMatch;
					//false means do not normalize.
					
					
					double thisTotal = patches[2 * f + inOut][p]->getNEigs() + 1; //1 for the residue.
					thisTotal += patches[2 * f + inOut][p]->getGasBoneMatters(); 
					total[inOut]  += thisTotal;

					if (f == 1)
						int qqq = 0;


					if (out) {
						out << f << " " <<  inOut << " " << p << " " << thisMatch << " " << thisTotal << endl;
					}	
				}
			}//default match
		}
		
		//We'll normalize after everything, like Eli.
		//theMatches[inOut] /= sqrt((double) total);
		//for the normalization thingy.
	}

	//the sum of two unit normals is normal with variance 2.
	
	//theMatch = (theMatches[0] + theMatches[1]) / sqrt(2.0);

	//Now I'm going to normalize the whole thing.  Most of the commented out code above is for 
	//normalizing each side separately or for normalizing each match (qf, gas, bone) separately.

	theMatch = (theMatches[0] + theMatches[1]) - (total[0] + total[1]);
	theMatch /= sqrt(2.0*(total[0] + total[1]));

    if (out) { out.close(); }
	  

	delete [] theMatches;
	delete [] total;

	return theMatch;
}

//This match is meant to be for an augmented set of atoms.  So, for example, you want a prostate
//with some bladder and rectum atom image match, not just prostate.  If you have a model for each
//spoke of each figure, for inside and for outside, then we can do this. 
//Apparently, we're only concerned with the exteriors of the augmenting spokes, but that should
//not matter, I need to keep it for when Steve says wait we do want both sides. So, two new
//parameters to replace the figureID f...
//WAIT:  I'll do it better, if f is >= 0, then use it to compute the whole match of that one figure, 
//and use the insides and outsides on the other figures.
//insides and outsides are each numfigs x something, where something is the indices of the spokes to 
//include in the match.
double MultiPatch::getYourMatch(const M3DObject * object, Image3D * image, int f, 
		std::vector<std::vector< int> > insides, std::vector<std::vector< int> > outsides, Image3D * binImage)
{
	double theMatch = 0.0;
	int fiter, p;

	//First, sample ourselves. But we have to do all the figures,
	for (fiter = 0; fiter < numFigures; fiter ++)
		sampleYourself(object, image, fiter, binImage);


	//Now we'll get the matches. if f is >= 0, it's assumed a figure number, and we'll compute
	//the whole match of that figure, so that vector in insides and outsides doesn't have to 
	//be just everything. Actually, I'll change the appropriate insides and outsides, so that 
	//I can run the match computation the same.
	if (f >= 0)
	{
		insides.erase(insides.begin() + f);
		outsides.erase(outsides.begin() + f);

		std::vector< int > tempSpokes;
		
		//use all the patches available for that figure.
		
		for(p = 0; p < patches[2 * f + 0].size(); p ++)
		{
			tempSpokes.push_back(p);
		}
		insides.insert(insides.begin() + f, tempSpokes);

		tempSpokes.clear();

		for(p = 0; p < patches[2 * f + 1].size(); p ++)
		{
			tempSpokes.push_back(p);
		}
		outsides.insert(outsides.begin() + f, tempSpokes);

		tempSpokes.clear();
	}


	//So now, we loop over each figure, over all the elements of that index of 
	//insides and outsides, gathering the matches.

	double * theMatches = new double[2*numFigures];
	double * total = new double[2*numFigures];
	double thisMatch, thisTotal;

	for (fiter = 0; fiter < numFigures; fiter++)
	{
		
		theMatches[2*fiter + 0] = 0.0;
		total[2*fiter + 0] = 0.0;
		
		//First, insides
		if (insides[fiter].size() == 0)
			break;
		
		//if I should count some inside patches, then loop through them.
		for(p = 0; p < insides[fiter].size(); p ++)
		{
			thisMatch = patches[2 * fiter + 0][insides[fiter][p]]->getYourMatch(false);
			theMatches[2 * fiter + 0] += thisMatch;
			//false means do not normalize.
			
			
			thisTotal = patches[2 * fiter + 0][insides[fiter][p]]->getNEigs() + 1; //1 for the residue.
			thisTotal += patches[2 * fiter + 0][insides[fiter][p]]->getGasBoneMatters(); 
			total[2 * fiter + 0]  += thisTotal;
		}
		
		//Then outsides
		theMatches[2*fiter + 1] = 0.0;
		total[2*fiter + 1] = 0.0;
		
		if (outsides[fiter].size() == 0)
			break;
		
		//if I should count some inside patches, then loop through them.
		int myp;
		for(p = 0; p < outsides[fiter].size(); p ++)
		{
			myp = outsides[fiter][p];
			thisMatch = patches[2 * fiter + 1][myp]->getYourMatch(false);
			theMatches[2 * fiter + 1] += thisMatch;
			//false means do not normalize.
			
			
			thisTotal = patches[2 * fiter + 1][myp]->getNEigs() + 1; //1 for the residue.
			thisTotal += patches[2 * fiter + 1][myp]->getGasBoneMatters(); 
			total[2 * fiter + 1]  += thisTotal;
		}
		
	}
	
	//From the other match function.  I'll do this with all figures.
	//theMatch = (theMatches[0] + theMatches[1]) - (total[0] + total[1]);
	//theMatch /= sqrt(2.0*(total[0] + total[1]));
	double totalSum = 0.0, totalMatch = 0.0;
	for (int i = 0; i < 2*numFigures; i ++)
	{
		totalMatch += theMatches[i];
		totalSum += total[i];
	}

	theMatch = totalMatch - totalSum;
	theMatch /= sqrt(2.0 * totalSum);

	delete [] theMatches;
	delete [] total;

	return theMatch;
}


//STOUGH: fill this out with an optimization of the regional delineation and subsequent 
//matches.
void MultiPatch::flexibleMatch(const M3DObject * object, Image3D * image, int f, int whichSide,
								 double &theMatch, double &total, Image3D *binImage)
{
	theMatch = 0.0;
	total = 1.0;
}

//In this match, for each patch we look through the patchPCA's for the best match, and 
//keep that one.
void MultiPatch::MLEMatch(const M3DObject * object, Image3D * image, int f, int whichSide,
								 double &theMatch, double &total, Image3D *binImage)
{
	int p, q;
	double thisMatch, thisTotal, tempMatch;
	
	theMatch = total = 0.0;

	for (p = 0; p < patches[2 * f + whichSide].size(); p ++)
	{
		thisMatch = 1.0e100; thisTotal = 1.0;

		//Now loop through the patchPCAs, searching for the best match.
		for (q = 0; q < patchPCAs[2 * f + whichSide].size(); q ++)
		{
			tempMatch = patchPCAs[2 * f + whichSide][q]->computeMatch(patches[2 * f + whichSide][p]->getQF(),
																patches[2 * f + whichSide][p]->getWeights(),
																false);

			if (tempMatch < thisMatch)
			{
				thisMatch = tempMatch;
				thisTotal = patchPCAs[2 * f + whichSide][q]->getNEigs() + 1; //1 for the residue.
				thisTotal += patchPCAs[2 * f + whichSide][q]->getGasBoneMatters(); 
			}
		}

		theMatch += thisMatch;
		total += thisTotal;
	}

}


//Again, as per Eli's save mechanism, such that the data can be read by matlab. 
//It's integrated now so that the output of training, the input to matlab, output
//of matlab and input to seg time is all the same format, that of Registry.
bool MultiPatch::saveMP(const char * saveFilename)
{
	// Adjust registry from creation, adding training case info, namely, the qf's and high
	// low weights of the patches.  Basically, record numFigures, then for each figure, 
	// record for inside and outside, numpatches and for each patch record spokeID and 
	// quantile function. whew.  Unlike Eli's concatenated vectors, this will have a lot 
	// of whitespace overhead, but much easier to read.

	//everything recorded here is inside a training subfolder.

	reg.setIntValue("training.numFigures", numFigures);


	int j, f;

	//double weights[3];

	char startStr[1024];
	char tmpStr[1024];
	char scndStr[1024];



	//if the patchType is DEFINED, we have to reset the delins. In all cases.  we 
	//need to reset the numPatches, since in testing it determines how many pca's 
	//we read etc.
	startStr[0] = tmpStr[0] = scndStr[0] = '\0';
	strcat(startStr, "config.figure");
	for (f = 0; f < numFigures; f ++)
	{
		sprintf(tmpStr, "%s[%d].inside", startStr, f);
		//now tmpStr is "config.figure[0].inside";

		reg.setIntValue("%s.numPatches", patches[2 * f].size(), tmpStr);

		if (samplingConfigs[2 * f].patchType == DEFINED_PATCHES)
		{
			//set the gasMatters and boneMatters arrays.
			reg.setDoubleArray("%s.gasMatters", samplingConfigs[2*f].numPatches, samplingConfigs[2*f].gasMatters, tmpStr);
			reg.setDoubleArray("%s.boneMatters", samplingConfigs[2*f].numPatches, samplingConfigs[2*f].boneMatters, tmpStr);

			for(j = 0; j < samplingConfigs[2*f].numPatches; j ++)
			{
				sprintf(scndStr, "%s.patch[%d]", tmpStr, j);
				reg.setDoubleArray(scndStr, 6, samplingConfigs[2*f].delins[j]);
			}
		}
		else
		{
			reg.setDoubleValue("%s.gasMatters", *samplingConfigs[2*f].gasMatters, tmpStr);
			reg.setDoubleValue("%s.boneMatters", *samplingConfigs[2*f].boneMatters, tmpStr);
		}

		sprintf(tmpStr, "%s[%d].outside", startStr, f);
		//now tmpStr is "config.figure[0].outside";

		reg.setIntValue("%s.numPatches", patches[2 * f + 1].size(), tmpStr);

		if (samplingConfigs[2*f + 1].patchType == DEFINED_PATCHES)
		{
			//set the gasMatters and boneMatters arrays.
			reg.setDoubleArray("%s.gasMatters", samplingConfigs[2*f+1].numPatches, samplingConfigs[2*f+1].gasMatters, tmpStr);
			reg.setDoubleArray("%s.boneMatters", samplingConfigs[2*f+1].numPatches, samplingConfigs[2*f+1].boneMatters, tmpStr);

			for(j = 0; j < samplingConfigs[2*f + 1].numPatches; j ++)
			{
				sprintf(scndStr, "%s.patch[%d]", tmpStr, j);
				reg.setDoubleArray(scndStr, 6, samplingConfigs[2*f + 1].delins[j]);
			}
		}
		else
		{
			reg.setDoubleValue("%s.gasMatters", *samplingConfigs[2*f+1].gasMatters, tmpStr);
			reg.setDoubleValue("%s.boneMatters", *samplingConfigs[2*f+1].boneMatters, tmpStr);
		}
	}


	//Now we're going to set the training info in the same registry.

	startStr[0] = tmpStr[0] = scndStr[0] = '\0';
	strcat(startStr, "training.figure");

	//sprintf(scndStr, "%s.inside.patch[%d]", tmpStr, j);

	for (f = 0; f < numFigures; f ++)
	{

		sprintf(tmpStr, "%s[%d].inside", startStr, f);
		//now tmpStr is "training.figure[0].inside";

		reg.setIntValue("%s.numPatches", patches[2*f].size(), tmpStr);

		
		for (j = 0; j < patches[2*f].size(); j ++)
		{
			reg.setIntArray("%s.patch[%d].spokeIDs", patches[2*f][j]->numSpokeIDs(), 
				              patches[2*f][j]->getSpokeIDs(), tmpStr, j);
			// training.figure[0].inside.patch[0].spokeIDs       #spokeIDs           ref to array of spokeIDs
			
			// now qf.

			reg.setDoubleArray("%s.patch[%d].qf", patches[2*f][j]->getEMSize(), patches[2*f][j]->getQF(), tmpStr, j);

			//now the weights.

			reg.setDoubleArray("%s.patch[%d].weights", 3, patches[2*f][j]->getWeights(), tmpStr, j); 

		}




		sprintf(tmpStr, "%s[%d].outside", startStr, f);
		//now tmpStr is "training.figure[0].outside";

		reg.setIntValue("%s.numPatches", patches[2*f + 1].size(), tmpStr);
		
		for (j = 0; j < patches[2*f + 1].size(); j ++)
		{
			reg.setIntArray("%s.patch[%d].spokeIDs", patches[2*f + 1][j]->numSpokeIDs(), 
				              patches[2*f + 1][j]->getSpokeIDs(), tmpStr, j);
			// training.figure[0].inside.patch[0].spokeIDs       #spokeIDs           ref to array of spokeIDs
			
			// now qf.

			reg.setDoubleArray("%s.patch[%d].qf", patches[2*f + 1][j]->getEMSize(), patches[2*f + 1][j]->getQF(), tmpStr, j);

			//now the weights.

			reg.setDoubleArray("%s.patch[%d].weights", 3, patches[2*f + 1][j]->getWeights(), tmpStr, j); 

		}
	}


	//matlab does the stat saving, there may be more info to save here later.


	bool retr = true;
	try {
		if (! reg.writeToFile(saveFilename))
			cout << "Warning: Output MultiPatch file may be corrupted" << endl;

	}
	catch (RException excp) {
		cout << excp.message() << endl;
		retr = false;
	}

	return retr;
}

