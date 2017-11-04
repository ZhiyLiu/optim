#include "DistanceToPointSetFunction.h"

#define D_XFERLIST
#define D_POINTLIST_SERVER2
#define D_LINALG
#include "Shapedepend.h" 
#include "M3DPGAStats.h"
#include <iostream>
#include <fstream>
#include "P3DControl.h"
#include "SimTransShapeSpace.h"
#include "ScaleTransShapeSpace.h"
#include "OctTreeTauBand.h"

#include <algorithm>

using std::ifstream;
using std::cout;
using std::endl;
using std::vector;
using std::sort;


DistanceToPointSetFunction::DistanceToPointSetFunction(
    const char * pointsFileName, SimTransType simType, int level,
	int figureId, Image3D * img, int padding)
{
	pgaOrder = -1;
    image = img;
	slicePadding = padding;
	if (pointsFileName)
        if (! loadPointSet(pointsFileName)) {
            simType = Fail;
            return;
        }

	subdivLevel = level;
	figure_id = figureId;
	sizeWeight = 0.0;
	if (simType == Full)
		space = new SimTransShapeSpace;
	else {
		space = new ScaleTransShapeSpace;
		if (simType == Scale)
			((ScaleTransShapeSpace *) space)->setType(ScaleTransShapeSpace::Scale);
		else if (simType == Translation)
			((ScaleTransShapeSpace *) space)->setType(ScaleTransShapeSpace::Translation);
	}

	setClosenessLimit(0.04);
	
	optimizationStep = None;
	referenceObject = NULL;
}

DistanceToPointSetFunction::DistanceToPointSetFunction(
    Vector3D * pointSet, int numPoints, SimTransType simType, int level,
	int figureId, Image3D * img, int padding)
{
	pgaOrder = -1;
    image = img;
	slicePadding = padding;
	if (pointSet)
        if (! addPointSet(pointSet, numPoints)) {
            simType = Fail;
            return;
        }

	subdivLevel = level;
	figure_id = figureId;
	sizeWeight = 0.0;
	if (simType == Full)
		space = new SimTransShapeSpace;
	else {
		space = new ScaleTransShapeSpace;
		if (simType == Scale)
			((ScaleTransShapeSpace *) space)->setType(ScaleTransShapeSpace::Scale);
		else if (simType == Translation)
			((ScaleTransShapeSpace *) space)->setType(ScaleTransShapeSpace::Translation);
	}
	setClosenessLimit(0.04);
	optimizationStep = None;
	referenceObject = NULL;
}

DistanceToPointSetFunction::~DistanceToPointSetFunction()
{
	delete space;
}

bool DistanceToPointSetFunction::loadPointSet(const char * pointsFileName)
{
	double x, y, z;

	points.clear();
	ifstream in(pointsFileName);
	in >> x >> y >> z;
	if (! in) {	
		cout << "Error: Cannot read points from file " << pointsFileName << endl;
		return false;
	}
	
	while (in) {
		in >> x >> y >> z;
		if (in) {
			Vector3D pnt(x, y, z);
			points.push_back(pnt);
		}
	}

	adjustPointSet();
	normalizePointSet();
    return true;
}

bool DistanceToPointSetFunction::addPointSet(Vector3D * pointSet,
	int numPoints)
{
    if (pointSet == NULL) return false;
	points.clear();
	for (int i = 0; i < numPoints; i++)
		points.push_back(pointSet[i]);

	adjustPointSet();
	normalizePointSet();
    return true;
}

// Convert the specified points from world to model coordinates
void DistanceToPointSetFunction::adjustPointSet() {  
	if (image) {
		for ( vector<Vector3D>::iterator it = points.begin(); 
			it != points.end();
			it++) {
				image->worldToModelCoordinates(*it);
		}		
	}

	// Levy: Because we are no longer assuming 3 designated
	// axial slices,  the adjustment to add an extra slice 
	// thickness no longer makes sense and has been removed.

	// Keeping track of the extent of the bounding box so that we can 
	// normalize coordinates for tau-band calculations.
	minX = points[0].getX();
    minY = points[0].getY();
    minZ = points[0].getZ();;
	maxX = minX;
    maxY = minY;
    maxZ = minZ;

    int numPoints = points.size();
    // AGG: This loop can be speeded up, but it is probably not worth the effort
	for (int i = 1; i < numPoints; i++) {
		if (points[i].getX() < minX) { minX = points[i].getX(); }
		else if (points[i].getX() > maxX) { maxX = points[i].getX(); }
		
		if (points[i].getY() < minY) { minY = points[i].getY(); }
		else if (points[i].getY() > maxY) { maxY = points[i].getY(); }
		
		if (points[i].getZ() < minZ) { minZ = points[i].getZ(); }
		else if (points[i].getZ() > maxZ) { maxZ = points[i].getZ(); }
	}
}

void DistanceToPointSetFunction::normalizePointSet()
{
	normalizedPoints.clear();

	double lenX = maxX - minX;
	double lenY = maxY - minY;
	double lenZ = maxZ - minZ;

	for (vector<Vector3D>::iterator it = points.begin();
		it != points.end(); it++) {

			Vector3D temp((it->getX() - minX) / lenX,
				(it->getY() - minY) / lenY,
				(it->getZ() - minZ) / lenZ);
			normalizedPoints.push_back(temp);
			
	}
}

M3DObject * DistanceToPointSetFunction::createTargetObject(const Vector & v)
{
#ifdef UNFLIPPED
	const Vector3D center(0.0, 1.0, 0.0);
#else
	const Vector3D center(0.0, 0.0, 0.0);
#endif

	M3DObject * xformed = referenceObject->assign();
	// Note: assuming model atoms are already selected

    // Perform a PGA deformation
	if (pgaOrder >= 0) {
		vector<double> coefs;

		M3DPGAStats * stats = referenceObject->getPGAStats();
		stats->setIntermedObj(xformed);

		if (optimizationStep == PGA) {
			for (int i = 0; i < v.size(); i++)
				coefs.push_back(v(i));
			// Apply the deformation
			stats->doPGADeform(xformed, coefs, pgaOrder, false /*no prediction*/);
		}
		else if (pgaCoefs.size() > 0)
			stats->doPGADeform(xformed, pgaCoefs, pgaOrder, false /*no prediction*/);
	}

    // Apply the optimized post simtrans
	SimilarityTransform3D varXform;
    if (optimizationStep == SimT)
		varXform = space->createSimTrans(v);
	else
		varXform = space->createSimTrans(postSimTransCoefs);
	xformed->applySimilarity(varXform, center);

	return xformed;
}

double DistanceToPointSetFunction::evaluate(const Vector & v)
{
	M3DObject * obj = createTargetObject(v);
	double dist = computeDistanceSqr(obj);
	delete obj;
	return dist;
}

double DistanceToPointSetFunction::computeDistanceSqr(M3DObject * obj)
{
	M3DFigure * fig = obj->getFigurePtr(figure_id);
	Xferlist * xferList = convertM3DtoXfer(fig);
	ThallCode::Pointlist_server2 pList;
	pList.init(xferList);
	pList.ComputeSubdivPointCloud(subdivLevel);

	int numPoints;
	Bpoint * bPoints = 0;
	pList.subdivboundaryinfo(&numPoints, &bPoints);

	if(numPoints > 0) {
		OctTreeTauBand mrepOctTree(bPoints[0].pnt);
		
		for (int i = 1; i<numPoints; i++) {
			mrepOctTree.insert(bPoints[i].pnt);
		}

		mrepOctTree.normalize();
		
		mrepOctTree.xLimit = true;
		mrepOctTree.yLimit = true;
		mrepOctTree.zLimit = true;

		double totalDistSqr = 0;
		vector<Vector3D>::iterator pIt = points.begin();
		vector<Vector3D>::iterator nIt = normalizedPoints.begin();
		while ((pIt != points.end()) && 
			(nIt != normalizedPoints.end())) {
				totalDistSqr += mrepOctTree.getDistanceSqr(*pIt,
					*nIt, sqrtClosenessLimit);
				pIt++;
				nIt++;
		}

		return (totalDistSqr); // returning sum, not distance since match will tweak it / (points.size());

	}
	return 0.0;
}

int DistanceToPointSetFunction::parameterCount()
{
    return space->parameterCount();
}

void DistanceToPointSetFunction::setPGACoefs(const Vector & pg)
{
    pgaCoefs.clear();
    for (int i = 0; i < pg.size(); i++) {
        pgaCoefs.push_back(pg(i));
    }
}

