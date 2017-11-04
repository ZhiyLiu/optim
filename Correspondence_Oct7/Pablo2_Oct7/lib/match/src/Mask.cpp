#include <math.h>
#include <iostream>
#include "Mask.h"
#include "DistanceMap3D.h"

//#define DEBUG

const int MASK_MEDIAL_SUBDIVISIONS = 7;


using namespace std;



/*  This function will now fill the mask elements array where
    each tile's array is separately rms 1 and mean 0.  This
    allows for the locally varying match templates of PROFILE_MATCH,
	now obsolete.
*/
Mask::Mask(M3DObject *referenceObject, int _figureId, float _cutoff, MatchType type,
           Image3D *trainingImage, bool useWindowing, int surfaceLevel,
		   const char * templateFilename, const char * profileFilename,
		   SimpleMaskFile * simpleMask)
{

    Vector3D coord,
			 vec0,
             vec1,
             vec2,
             vec3,
             tmpvec;
    int i, j;
    MaskElement * newElement;
    M3DFigure * figure;
    DistanceMap3D distanceMap;
    Bpoint * bptlistptr = NULL;

    int numPts,nb1, nb2, nb3;
    double val1,
           totval1,
           totval1sq,
           area,
           radius,
           volume,
           dist,
           c1,
           c2,
           mu,
           rms;
    bool inBlendRegion;
    bool useDistanceMap;


	MatchType curtype = type;
	ok = true;
    figureId = _figureId;

    if (figureId < 0 || figureId >= referenceObject->getFigureCount())
    {
        cout << "Invalid figure id in Mask." << endl;
		ok = false;
        return;
    }

    cutoff = _cutoff;
    numstep = MASK_NUM_STEPS; // HARDCODED
    sigma = cutoff / 2.0;
	notch_sigma = .6*sigma;

#ifndef BINARY
    SamplesPerPoint = 2*numstep + 1;
#endif

    double sigmaSq = sigma * sigma;
    double sigmaCube = sigmaSq * sigma;
	double sigmaFifth = sigmaSq*sigmaCube;

	double nsigmaSq = notch_sigma * notch_sigma;
    double nsigmaCube = nsigmaSq * notch_sigma;
	double nsigmaFifth = nsigmaSq*nsigmaCube;

	// qiong: the DistanceMap3D has been updated for multi-figure object
    // HACK: WE ARE ONLY DOING SINGLE-OBJECT MULTIPLE-FIGURE
    // OR MULTIPLE-OBJECT SINGLE-FIGURE
    if (referenceObject->getFigureTreeCount() == 1 &&
       referenceObject->getFigureCount() > 1)
    {
        useDistanceMap = false;
    }
    else
        useDistanceMap = false;

	// the boundary tiles used to
	int num_n = 0;
	ThallCode::Tileneighbors * vn = NULL;
	figure = referenceObject->getFigurePtr(figureId);
    if (figure == NULL) {
		cout << "No figure in Mask" << endl;
		ok = false;
        return;
	}

	referenceObjectPtr = referenceObject;
	stackedImageMask = figure->getStackedImageMask();

	//Set the figureStats pointer if we have figure stats, rather than the 
	//profile file from earlier.  In that the profile files can be merged 
	//with the figures they correspond to in the m3d and read at that time,
	//this is the same functionality with less user input.  
	figureStats = figure->getFigureStatsPtr();

	Xferlist * xferList = convertM3DtoXfer(figure);
	ThallCode::Pointlist_server2 * pList = new ThallCode::Pointlist_server2;
	pList->init(xferList);
	pList->ComputeSubdivPointCloud(surfaceLevel);
	pList->subdivboundaryinfo(&numPts, &bptlistptr);
	pList->subdivvertexneighbors(&num_n, &vn);

#ifndef BINARY
	// the new neighboring vertices list
	ThallCode::NeighborVertices *nVerts = NULL;
	pList->subdivNeighborVertices(surfaceLevel, &num_n, &nVerts);
#endif

	int vertexToRemoveNum;
	short *vertexMapToRemove=NULL;
	if(referenceObject==NULL) {
		cout << "No reference object" << endl;
		ok = false;
        return;
	}

	// alright alright, still a hack to deal w/ 2-figure object, multi-level is yet to come
	ThallCode::Pointlist_serverB *pListB=new ThallCode::Pointlist_serverB;
	if (referenceObject->getFigureTreeCount() == 1 &&
		referenceObject->getFigureCount() > 1)
	{
		// future work to do: change referenceObject to the current object tree!
		pListB->init(referenceObject, surfaceLevel);	//pListB->init(referenceObject, surfaceLevel);
		pListB->ComputeFigureMapToRemove(&vertexToRemoveNum, &vertexMapToRemove, surfaceLevel, _figureId);
	}

	if (num_n != numPts || vn == NULL) {
		printf("Error in vertex neighbor computation in Mask::Mask.\n");
		ok = false;
		return;
	}

	// Here we read the template information, so that what is a gaussian derivative
	// is whatever is decided outside of the program, as well for neg gauss and notch.
	double long_light_dark[11], long_dark_light[11], long_notch[11];

	if (templateFilename != NULL) {
		FILE *filterfile;
		double tempd;
		filterfile = fopen(templateFilename, "r");

		if (filterfile == NULL) {
			printf("Mask.cpp: template file not found.\n");
			templateFilename = NULL;
		}
		else {
			// Read in the info.
			fscanf(filterfile, "long_light_dark\n");
			for (i = 0; i < 11; i ++) {
				fscanf(filterfile, "%lf\n", &tempd);
				long_light_dark[i] = tempd;
			}
			fscanf(filterfile, "long_dark_light\n");
			for (i = 0; i < 11; i ++) {
				fscanf(filterfile, "%lf\n", &tempd);
				long_dark_light[i] = tempd;
			}
			fscanf(filterfile, "long_notch\n");
			for (i = 0; i < 11; i ++) {
				fscanf(filterfile, "%lf\n", &tempd);
				long_notch[i] = tempd;
			}
			fclose(filterfile);
		}
	}

	totalVolume = 0.0;
	totval1 = 0;
	totval1sq = 0;

#ifndef BINARY

	//If SIMPLE_MASK_MATCH, then we can set the dist incrementer here.
	double tau_start, sinterval;

	if (type == SIMPLE_MASK_MATCH && simpleMask != NULL)
	{
		tau_start = simpleMask->getTau(figureId, 0);
		sinterval = (simpleMask->getTau(figureId, 1) - tau_start)/(simpleMask->getNumSamples(figureId)-1);

		if (numPts > simpleMask->getNumPoints(figureId)) {
			cout << "Incorrect number of points in simple mask for current surface level" << endl;
			ok = false;
			return;
		}
	}

#endif

	for (i = 0; i < numPts; i++)
	{
		vec0.set(bptlistptr[i].pnt[0], bptlistptr[i].pnt[1], bptlistptr[i].pnt[2]);
		area = 0;

#ifndef BINARY
		for (j = 0; j < nVerts[i].nNumber*2; j += 2)
#else
		for (j = 0; j < MAX_NEIGHBOR_COUNT; j += 2)
#endif
		{
#ifndef BINARY
			nb1 = j;
			nb2 = (j+1) % (nVerts[i].nNumber*2);
			nb3 = (j+2) % (nVerts[i].nNumber*2);

			vec1.set(nVerts[i][nb1].pnt[0], nVerts[i][nb1].pnt[1], nVerts[i][nb1].pnt[2]);
			vec2.set(nVerts[i][nb2].pnt[0], nVerts[i][nb2].pnt[1], nVerts[i][nb2].pnt[2]);
			vec3.set(nVerts[i][nb3].pnt[0], nVerts[i][nb3].pnt[1], nVerts[i][nb3].pnt[2]);

#else
			nb1 = vn[i][j];
			nb2 = vn[i][j+1];
			nb3 = vn[i][(j+2)%MAX_NEIGHBOR_COUNT];
			if (nb1 < 0 || nb1 >= numPts || nb2 <0 || nb2 >= numPts || nb3 < 0|| nb3 >= numPts)
			{
				if(nb1>=numPts || nb2>=numPts || nb3>=numPts)
					cout << "Invalid vertex neighbor index " << nb1 << "," << nb2 << "," << nb3 \
					<< " in Mask::Mask(...)" << endl;
				break;
			}	// 2003/07/07

			vec1.set(bptlistptr[nb1].pnt[0], bptlistptr[nb1].pnt[1], bptlistptr[nb1].pnt[2]);
			vec2.set(bptlistptr[nb2].pnt[0], bptlistptr[nb2].pnt[1], bptlistptr[nb2].pnt[2]);
			vec3.set(bptlistptr[nb3].pnt[0], bptlistptr[nb3].pnt[1], bptlistptr[nb3].pnt[2]);

			vec2 = .25 * (vec0 + vec1 + vec2 + vec3);
			vec1 = 0.5 * (vec0 + vec1);
			vec3 = 0.5 * (vec0 + vec3);
#endif

			tmpvec = vec1 - vec2;
			area += 0.5 * (tmpvec.cross(vec3 - vec2)).norm();

			tmpvec = vec3 - vec0;
			area += 0.5 * (tmpvec.cross(vec1 - vec0)).norm();
		}

        radius = bptlistptr[i].rad;

#ifndef BINARY
        if (type == SIMPLE_MASK_MATCH && simpleMask != NULL)
		{
			double extent = simpleMask->getTau(figureId, 1) - simpleMask->getTau(figureId, 0);
			volume = area * extent * radius / (simpleMask->getNumSamples(figureId));
		}
		else
#endif
			volume = area * cutoff * radius / (2 * numstep + 1);

		if(vertexMapToRemove!=NULL)
			inBlendRegion=(vertexMapToRemove[i]==1);
		else
			inBlendRegion=false;

        // March along +/-numstep subdiv of the cutoff in rnormalized terms
        // along the normal to the 1st point
        // Do the first triangle
		curtype = type;
		if (inBlendRegion)
		{
	        for (j = -numstep; j <= numstep; j++)
			{
				// in fact these initialization of the mask contents could be removed
				dist = ((double) j) * cutoff / numstep;
				if (dist < -1.0) {
					printf("Dist < -1\n");
					ok = false;
					return;
				}
				coord.set(bptlistptr[i].norm);
				coord *= (dist * radius);
				coord += vec0;
				newElement = new MaskElement;
				newElement->x = coord.getX();
				newElement->y = coord.getY();
				newElement->z = coord.getZ();
				newElement->weight = 1.0;
				newElement->figureId = figureId;
				newElement->u = bptlistptr[i].u;
				newElement->v = bptlistptr[i].v;
				newElement->t = bptlistptr[i].t;
				newElement->dist = dist;
				newElement->volelem = volume;
				maskElement.insert(maskElement.end(), newElement);

				// the only thing matters might be just this line
				newElement->inBlendRegion = inBlendRegion;
			}	// if (in the blended region)
		}
#ifdef BINARY
		else
#else
		else if (type != SIMPLE_MASK_MATCH)
#endif
		{
			for (j = -numstep; j <= numstep; j++)
			{	
				dist = ((double) j) * cutoff / numstep;

				if (dist < -1.0) {
					printf("Dist < -1\n");
					ok = false;
					return;
				}

				coord.set(bptlistptr[i].norm);
				coord *= (dist * radius);
				coord += vec0;

				newElement = new MaskElement;
				newElement->x = coord.getX();
				newElement->y = coord.getY();
				newElement->z = coord.getZ();
				newElement->weight = 1.0;
				newElement->figureId = figureId;
				newElement->u = bptlistptr[i].u;
				newElement->v = bptlistptr[i].v;
				newElement->t = bptlistptr[i].t;
				newElement->dist = dist;
				newElement->volelem = volume;
				newElement->inBlendRegion = inBlendRegion;

				// Now do the different template types
				if (curtype == GAUSSIAN_DERIVATIVE_MATCH)
				{
					c1 = -1.0 / sigmaCube;
					c2 = -0.5 / sigmaSq;
					if (templateFilename != NULL)
						val1 = long_light_dark[j+numstep];
					else val1 = c1 * dist * exp(c2 * dist * dist);

					if (!(figure->hasPositivePolarity()))
						val1 = -val1;
				}

				// Dark to light gaussian.  Almost identical to Gaussian.
				if (curtype == NEG_GAUSSIAN_DERIVATIVE_MATCH)		
				{
					c1 = -1.0 / sigmaCube;
					c2 = -0.5 / sigmaSq;
					if (templateFilename != NULL)
						val1 = long_dark_light[j+numstep];
					else val1 = -c1 * dist * exp(c2 * dist * dist);
				}

#ifdef BINARY
				if (curtype == SIMPLE_MASK_MATCH && simpleMask != NULL)
				{
					if (! simpleMask->getValue(figureId, i, j + numstep, val1))
						return;
				}
#endif

				newElement->tempVal = val1;

				maskElement.insert(maskElement.end(), newElement);
				// Need this for subtracting out the mean
				totval1 += val1*volume;
				totalVolume += newElement->volelem;
				totval1sq += val1 * val1 * volume;
			}
		}
#ifndef BINARY
		else {
			//type is SIMPLE_MASK_MATCH.  In this version, every figure 
			//can have different tau's and samples, so we must load the 
			//specific information. dist is the thing that needs to be
			//different from above.

			for (j = 0; j < simpleMask->getNumSamples(figureId); j ++)
			{
				dist = tau_start + ((double) j)*sinterval;

				if (dist < -1.0) {
					printf("Dist < -1\n");
					ok = false;
					return;
				}
				coord.set(bptlistptr[i].norm);
				coord *= (dist * radius);
				coord += vec0;

				newElement = new MaskElement;
				newElement->x = coord.getX();
				newElement->y = coord.getY();
				newElement->z = coord.getZ();
				newElement->weight = 1.0;
				newElement->figureId = figureId;
				newElement->u = bptlistptr[i].u;
				newElement->v = bptlistptr[i].v;
				newElement->t = bptlistptr[i].t;
				newElement->dist = dist;
				newElement->volelem = volume;
				newElement->inBlendRegion = inBlendRegion;

				if (! simpleMask->getValue(figureId, i, j, val1)) {
					ok = false;
					return;
				}
				newElement->tempVal = val1;

				maskElement.insert(maskElement.end(), newElement);
				// Need this for subtracting out the mean
				totval1 += val1*volume;
				totalVolume += newElement->volelem;
				totval1sq += val1 * val1 * volume;
			}
		}
#endif
    }	// for

	/*  Here we normalize the entire mask (which is the concatenation of the template
		profiles for each point, which may also be called the template stencil).  This 
		is done in a way such that when image match is computed later on, the maximal
		possible image match of a consistently normalized target stencil (normalized
		in the same way) is 1.0.
	*/

	//if profile match, then don't do the normalization, because they're
	//already done.
	mu = totval1 / totalVolume;
	rms = (totval1sq / totalVolume) - (mu * mu);

	if(rms <= R_SMALL_TOLERANCE)
		rms = 0.0;
	else
		rms = sqrt(rms);

	// Now need to normalize the template value for each figure
//	double tt;
	for(i = 0; i < maskElement.size(); i++)
	{
//		tt = (maskElement[i]->tempVal - mu) / rms;
		maskElement[i]->tempVal = (maskElement[i]->tempVal - mu) / rms;
	}

#ifndef BINARY
	if (type == SIMPLE_MASK_MATCH)
	{
		SamplesPerPoint = simpleMask->getNumSamples(figureId);
	}
#endif

    delete pList;
	delete pListB;
    delete [] xferList->atomlist;
    delete xferList;
	delete []vertexMapToRemove;
}

// This function expects an index into the maskElements array, and 
// returns a corresponding rms value, in the profile match case.
double Mask::getRms(int i)
{
	int ind = (int) i/SamplesPerPoint;

	if (figureStats != NULL)
		return figureStats->getRms(ind);
	else
		return mt.getRms(ind);
}

double Mask::getMeanOffset(int i) 
{ 
	if (figureStats != NULL)
		return figureStats->getMeanOffset(i);
	else
		return mt.getMeanOffset(i); 
}

double Mask::getMeanStd(int i) 
{ 
	if (figureStats != NULL)
		return figureStats->getMeanStd(i);
	else
		return mt.getMeanStd(i); 
}

double Mask::getPMatchVar() 
{ 
	if (figureStats != NULL)
		return figureStats->getProfMatchVar();
	else
		return mt.getPMatchVar(); 
}

double Mask::getMMatchVar() 
{ 
	if (figureStats != NULL)
		return figureStats->getMeanMatchVar();
	else
		return mt.getMMatchVar(); 
}

#ifndef BINARY

double Mask::getCenterIntensity()
{
	if (figureStats != NULL)
	{
		printf("FigureStats center intensity not defined.\n");
		return 1.0;  //CHANGE THIS STOUGH
	}
	else
		return mt.getCenterIntensity();
}

#endif

Mask::~Mask()
{
    int size;
    int i;

    size = maskElement.size();

    for(i = 0; i < size; i++)
        delete maskElement[i];
}

void Mask::setVals(int _figureId, double _cutoff, double _sigma,
        bool _includeInterior, const std::vector<MaskElement *> & elements)
{
    int size;
    int i;

    figureId = _figureId;
    cutoff = _cutoff;
    sigma = _sigma;
    includeInterior = _includeInterior;

    size = maskElement.size();
    for(i = 0; i < size; i++)
        delete maskElement[i];
    maskElement.erase(maskElement.begin(), maskElement.end());

    size = elements.size();
    totalVolume = 0.0;
    for(i = 0; i < size; i++)
    {
        maskElement.push_back(elements[i]);

        if(elements[i] != NULL)
            totalVolume += elements[i]->volelem;
    }
}

void Mask::calcPositions(M3DObject * object)
{
    int i, j, size;
    M3DFigure * figure;
    ThallCode::Pointlist_server2 * pList;
    Bpoint * bptlistptr,
           * tmpPtr;
    Xferlist * xferList;
    int ntiles;
    Vector3D coord;
    double radius,
           dist;

	printf("Mask.cpp: calcPositions not modified to use point surfaces.\n");

    if(figureId < 0 || figureId >= object->getFigureCount())
    {
        cout << "Invalid figure id in Mask." << endl;
        return;
    }

    numstep = MASK_NUM_STEPS; // HARDCODED

    figure = object->getFigurePtr(figureId);

    if(figure == NULL)
        return;

    xferList = convertM3DtoXfer(figure);

    pList = new ThallCode::Pointlist_server2;
    pList->init(xferList);
    pList->ComputeSubdivBoundaryTiles(MATCH_POINT_CLOUD_SUBDIVISIONS);
    pList->subdivtileinfo(&ntiles, &bptlistptr);

    tmpPtr = bptlistptr;
    int index = 0;
    size = maskElement.size();
    for(i = 0; i < ntiles; i++, tmpPtr += 4)
    {
        radius = tmpPtr[0].rad;

        // March along +/-numstep subdiv of the cutoff in rnormalized terms
        // along the normal to the 1st point
        // Do the first triangle
        for(j = -numstep; j <= numstep && index < size; j++)
        {
            dist = ((double) j) * cutoff / numstep;

            coord.set(tmpPtr[0].norm);
            coord *= (dist * radius);
            coord += Vector3D(tmpPtr[0].pnt);

            if(maskElement[index] == NULL)
            {
                cout << "Mask::calcPositions(): invalid mask." << endl;
                break;
            }

            maskElement[index]->x = coord.getX();
            maskElement[index]->y = coord.getY();
            maskElement[index]->z = coord.getZ();

            index++;
        }
    }

    delete pList;
    delete [] xferList->atomlist;
    delete xferList;
}

DMask::DMask(M3DObject * object, int _figureId, double _cutoff, MatchType _type,
           int subdivLevel, const char * templateFilename, const char * profileFilename)
{
	int i, k;

	M3DFigure * figure = object->getFigurePtr(_figureId);
	if (figure == NULL)
		return;

	figureId = _figureId;

	type = _type;

	// Here pick up the dmask to use from templateFilename.
    double light_dark[10] = {0};
    double dark_light[10] = {0};
    double notch[10] = {0};
	double long_light_dark[11], long_dark_light[11], long_notch[11];

	if (templateFilename != NULL) {
		FILE *filterfile;
		double tempd;

		filterfile = fopen(templateFilename, "r");

		if (filterfile == NULL) {
			printf("Mask.cpp: template file not found.\n");
			templateFilename = NULL;
		}
		else {
			// Get the mask creation numbers.
			fscanf(filterfile, "long_light_dark\n");
			for (i = 0; i < 11; i ++) {
				fscanf(filterfile, "%lf\n", &tempd);
				long_light_dark[i] = tempd;
			}
			fscanf(filterfile, "long_dark_light\n");
			for (i = 0; i < 11; i ++) {
				fscanf(filterfile, "%lf\n", &tempd);
				long_dark_light[i] = tempd;
			}
			fscanf(filterfile, "long_notch\n");
			for (i = 0; i < 11; i ++) {
				fscanf(filterfile, "%lf\n", &tempd);
				long_notch[i] = tempd;
			}

			// Now read the info.
			fscanf(filterfile, "light_dark\n");
			for (i = 0; i < 10; i ++) {
				fscanf(filterfile, "%lf\n", &tempd);
				light_dark[i] = tempd;
			}
			fscanf(filterfile, "dark_light\n");
			for (i = 0; i < 10; i ++) {
				fscanf(filterfile, "%lf\n", &tempd);
				dark_light[i] = tempd;
			}
			fscanf(filterfile, "notch\n");
			for (i = 0; i < 10; i ++) {
				fscanf(filterfile, "%lf\n", &tempd);
				notch[i] = tempd;
			}
		}
	}

	// Here we also need to pick up the profile from the figure.
	numstep = MASK_NUM_STEPS;	// Hard coded
	double dist, val1, val2, c1, c2;

	cutoff = _cutoff;
    sigma = cutoff / 2; 
	notch_sigma = .6*sigma;

    double sigmaSq = sigma * sigma;
    double sigmaCube = sigmaSq * sigma;
	double sigmaFifth = sigmaSq*sigmaCube;

	double nsigmaSq = notch_sigma * notch_sigma;
    double nsigmaCube = nsigmaSq * notch_sigma;
	double nsigmaFifth = nsigmaSq*nsigmaCube;

	int index;

	double dld[10], ddl[10], dnotch[10];

	dld[0] = (light_dark[1] - light_dark[0])/0.06;
	ddl[0] = (dark_light[1] - dark_light[0])/0.06;
	dnotch[0] = (notch[1] - notch[0])/0.06;

	for (k = 1; k < 9; k ++) {
		dld[k] = (light_dark[k+1] - light_dark[k-1])/0.12;
		ddl[k] = (dark_light[k+1] - dark_light[k-1])/0.12;
		dnotch[k] = (notch[k+1] - notch[k-1])/0.12;
	}

	dld[9] = (light_dark[9] - light_dark[8])/0.06;
	ddl[9] = (dark_light[9] - dark_light[8])/0.06;
	dnotch[9] = (notch[9] - notch[8])/0.06;


	// With statistically determined ideal templates,
	// these dervative computations are no longer analytic
	// (except for the never used abs. gauss.).  Take the 
	// forward and backward difference at the ends, and the
	// symmetric difference in the middle, of light_dark, 
	// dark_light, and notch.  The distances are then .06 first
	// and last, with .12 in the middle.

	for (k = -numstep-1; k < numstep+1; k++)
	{
		index = k + numstep;

		dist = ((double) k + 0.5) * cutoff / numstep;

		if (k == -numstep-1) dist = -cutoff;
		if (k == numstep) dist = cutoff;

		// Depending on the type of match desired.
		switch (type)
		{
			case BINARY_IMAGE_MATCH:
			case GAUSSIAN_DERIVATIVE_MATCH:
				c1 = -1.0 / sigmaCube;
				c2 = -0.5 / sigmaSq;
				if (templateFilename != NULL) {
					if (k == -numstep - 1)
						val1 = long_light_dark[0];
					else if (k == numstep)
						val1 = long_light_dark[10];
					else {
						val1 = light_dark[index];
						val2 = dld[index];
					}
				} else {
					val1 = c1 * dist * exp(c2 * dist * dist);
					val2 = c1 * exp(c2* dist * dist) - dist / sigmaSq * val1;
				}

				if (! (figure->hasPositivePolarity())) 
				{
					val1 = -val1;
					val2 = -val2;
				}
				break;
			case NEG_GAUSSIAN_DERIVATIVE_MATCH:
				c1 = -1.0 / sigmaCube;
				c2 = -0.5 / sigmaSq;
				if (templateFilename != NULL) {
					if (k == -numstep - 1)
						val1 = long_dark_light[0];
					else if (k == numstep)
						val1 = long_dark_light[10];
					else {
						val1 = dark_light[index];
						val2 = ddl[index];
					}
				}
				else {
					val1 = -(c1 * dist * exp(c2 * dist * dist));
					val2 = -1.0*c1 * exp(c2* dist * dist) - dist / sigmaSq * val1;
				}
				break;
		}

		int not_used;
		not_used = 0;

		if (k== -numstep-1)
			cutoffValues[0] = val1;
		else if (k == numstep)
			cutoffValues[1] = val1;
		else
		{
			DerivMaskElement * newDElement = new DerivMaskElement;
			newDElement -> dist = dist;
			newDElement -> tempVal = val1;
			newDElement -> dTempVal = val2;
			dmaskElements.insert(dmaskElements.end(), newDElement);
		}
	}

	initializeSurfaceEls(figure, subdivLevel);
}

DMask::~DMask()
{
	if (!dmaskElements.empty()) {
		for( int i = 0; i != dmaskElements.size(); ++i ) {
			delete dmaskElements[i];
		}
		dmaskElements.clear();
	}

#ifndef BINARY
	if (! dmaskElementsAll) {
		for (int i = 0; i < mt.getNumTemplateTypes(); i ++)
			delete [] dmaskElementsAll[i];
		delete [] dmaskElementsAll;
	}
#endif

	for (int l = 0; l < MAX_SUBDIV_LEVEL; l++)
		if (!surfaceElements[l].empty()) {
			for( int i = 0; i != surfaceElements[l].size(); ++i ) {
				delete surfaceElements[l][i];
			}
			surfaceElements[l].clear();
		}
}

/*  RMS and mean zero must be done per boundary point anymore (if there
    there are different templates locally).
*/
void DMask::initializeSurfaceEls(M3DFigure * figure, int level)
{
    Bpoint * bptlistptr;
    int numPts,
		i, j, 
		nb1, nb2, nb3;
    Vector3D vec0,
             vec1,
             vec2,
             vec3,
             tmpvec;
    double area,
           mu,
           rms;

	int nMaskEl = getMaskSize();

	if (!surfaceElements[level].empty()) {
		for( i = 0; i != surfaceElements[level].size(); ++i ) {
			delete surfaceElements[level][i];
		}
		surfaceElements[level].clear();
	}

    Xferlist * xferList = convertM3DtoXfer(figure);

    ThallCode::Pointlist_server2 * pList = new ThallCode::Pointlist_server2;
    pList->init(xferList);
	pList -> ComputeSubdivPointCloud(level);
	pList -> subdivboundaryinfo(&numPts, &bptlistptr);

	int num_n = 0;
	ThallCode::Tileneighbors * vn = NULL;
	pList -> subdivvertexneighbors(&num_n, &vn);

#ifndef BINARY
	// the new neighboring vertices list for the area calculation
	ThallCode::NeighborVertices *nVerts = NULL;
	pList->subdivNeighborVertices(level, &num_n, &nVerts);
#endif

	if (num_n != numPts || vn == NULL) 
	{
		printf("Error in vertex neighbor computation in DMask::initSurfEl.\n");
		return;
	}

	totalVolume = 0;

	for (i = 0; i < numPts; i++)
	{
		vec0.set(bptlistptr[i].pnt[0], bptlistptr[i].pnt[1], bptlistptr[i].pnt[2]);
		area = 0;

#ifdef BINARY
		for (j = 0; j < MAX_NEIGHBOR_COUNT; j += 2)
#else
		for (j = 0; j < nVerts[i].nNumber*2; j += 2)
#endif
		{
#ifndef BINARY
			nb1 = j;
			nb2 = (j+1) % (nVerts[i].nNumber*2);
			nb3 = (j+2) % (nVerts[i].nNumber*2);

			vec1.set(nVerts[i][nb1].pnt[0], nVerts[i][nb1].pnt[1], nVerts[i][nb1].pnt[2]);
			vec2.set(nVerts[i][nb2].pnt[0], nVerts[i][nb2].pnt[1], nVerts[i][nb2].pnt[2]);
			vec3.set(nVerts[i][nb3].pnt[0], nVerts[i][nb3].pnt[1], nVerts[i][nb3].pnt[2]);

#else
			nb1 = vn[i][j];
			nb2 = vn[i][j+1];
			nb3 = vn[i][(j+2)%MAX_NEIGHBOR_COUNT];
			if (nb1 < 0 || nb1 >= numPts || nb2 <0 || nb2 >= numPts || nb3 < 0|| nb3 >= numPts)
			{
				if(nb1>=numPts || nb2>=numPts || nb3>=numPts)
					cout << "Invalid vertex neighbor index " << nb1 << "," << nb2 << "," << nb3 \
						<< " in DMask::initializeSurfaceEls(...)" << endl;
				break;
			}	// 2003/07/07

			vec1.set(bptlistptr[nb1].pnt[0], bptlistptr[nb1].pnt[1], bptlistptr[nb1].pnt[2]);
			vec2.set(bptlistptr[nb2].pnt[0], bptlistptr[nb2].pnt[1], bptlistptr[nb2].pnt[2]);
			vec3.set(bptlistptr[nb3].pnt[0], bptlistptr[nb3].pnt[1], bptlistptr[nb3].pnt[2]);

			vec2 = .25 * (vec0 + vec1 + vec2 + vec3);
			vec1 = 0.5 * (vec0 + vec1);
			vec3 = 0.5 * (vec0 + vec3);
#endif

			tmpvec = vec1 - vec2;
			area += 0.5 * (tmpvec.cross(vec3 - vec2)).norm();

			tmpvec = vec3 - vec0;
			area += 0.5 * (tmpvec.cross(vec1 - vec0)).norm();
		}

		SurfaceElement * newSurfEl = new SurfaceElement;
		newSurfEl->area = area;

        surfaceElements[level].insert(surfaceElements[level].end(), newSurfEl);

		double radius = bptlistptr[i].rad;
		double volelem = area * 2 * cutoff * radius;
		totalVolume += volelem;
	}

    delete pList;

	double sum = 0.0, sum_sq = 0.0;

	for (j = 0; j < nMaskEl; j++)
	{
		// The first parameter to getMaskElement doesn't
		// matter if type isn't PROFILE_MATCH (obsolete).
		DerivMaskElement * dmaskEl = getMaskElement(0,j);
		sum += dmaskEl->tempVal;
		sum_sq += dmaskEl->tempVal*dmaskEl->tempVal;
	}

	mu = sum/(1.0*nMaskEl);

	rms = sum_sq - (1.0/nMaskEl)*sum*sum;
	rms = (1.0/(nMaskEl-1.0))*rms;

	if(rms <= R_SMALL_TOLERANCE)
		rms = 0.0;
	else
		rms = sqrt(rms);

	double debug_d;

	for (j = 0; j < nMaskEl; j++)
	{
		DerivMaskElement * dmaskEl = getMaskElement(0,j);
		dmaskEl->tempVal = (dmaskEl->tempVal - mu) / rms;
		dmaskEl->dTempVal = dmaskEl->dTempVal / rms;
		debug_d = 0;
	}

	cutoffValues[0] = (cutoffValues[0] - mu) / rms;
	cutoffValues[1] = (cutoffValues[1] - mu) / rms;

    delete [] xferList->atomlist;
    delete xferList;
}

DerivMaskElement * DMask::getMaskElement(int i, int j) 
{ 
		return dmaskElements[j]; 
}

double DMask::getCutoffValue(int j, int i)
{
		return cutoffValues[j];
}

int DMask::getTemplateType(int i)
{
	return mt.getProfileType(i);
}


