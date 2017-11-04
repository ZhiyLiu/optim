#include "M3DFigure.h"
#include "ControlParms.h"
#include "utility.h"

#include <typeinfo>

//#define DEBUG


const float DEFAULT_FIGURE_COLOR[3] = {0.2f, 0.2f, 0.8f};

int M3DFigure::defaultFigureTolerance = 0;


using namespace std;


M3DFigure::M3DFigure()
{
#ifdef DEBUG
    cout << "M3DFigure::M3DFigure()" << endl;
#endif

    properties = new struct FigProps_t;
    properties->counter = 1;

    properties->name = NULL;

    properties->modified = false;

    properties->color[0] = DEFAULT_FIGURE_COLOR[0];
    properties->color[1] = DEFAULT_FIGURE_COLOR[1];
    properties->color[2] = DEFAULT_FIGURE_COLOR[2];

    properties->positiveSpace = true;
    properties->positivePolarity = true;
    properties->visibility = true;
    properties->tolerance = getDefaultSurfaceTolerance();

    markedLandmark = -1;
#ifdef LM_METHOD_OBJ
    for (int i = 0; i < 99; i++)
    {
        landmarkAtomIndices[i] = -1;
        landmarkAtomTs[i] = 0;		// default is 'b spoke'
        landmarkAtomUs[i] = -1;		// default is invalid
        landmarkAtomVs[i] = -1;		// default is invalid
        landmarkNames[i] = NULL;
    }
    landmarkCount = 0;
#endif
    figureStats = NULL;
    stackedImageBit = 0;
}

M3DFigure::M3DFigure(const M3DFigure & fig)
{
#ifdef DEBUG
	cout << "M3DFigure::M3DFigure(const M3DFigure & fig)" << endl;
#endif

	properties = new struct FigProps_t;
	properties->counter = 1;

	// Copy the name
    if (fig.properties->name != NULL)
    {
        properties->name = new char[strlen(fig.properties->name) + 1];
        strcpy(properties->name, fig.properties->name);
    }
    else
        properties->name = NULL;

    properties->modified = fig.properties->modified;

    properties->color[0] = fig.properties->color[0];
    properties->color[1] = fig.properties->color[1];
    properties->color[2] = fig.properties->color[2];

    properties->positiveSpace = fig.properties->positiveSpace;
    properties->positivePolarity = fig.properties->positivePolarity;
	properties->visibility = fig.properties->visibility;
	properties->tolerance = fig.properties->tolerance;

	ifconstraints = fig.ifconstraints;
	inverse_constraints = fig.inverse_constraints;

	landmarks = fig.landmarks;
#ifdef LM_METHOD_OBJ
	memcpy(landmarkNames, fig.landmarkNames, sizeof(landmarkNames));
	memcpy(landmarkAtomIndices, fig.landmarkAtomIndices, sizeof(landmarkAtomIndices));
	memcpy(landmarkAtomTs, fig.landmarkAtomTs, sizeof(landmarkAtomTs));
	memcpy(landmarkAtomUs, fig.landmarkAtomUs, sizeof(landmarkAtomUs));
	memcpy(landmarkAtomVs, fig.landmarkAtomVs, sizeof(landmarkAtomVs));
	landmarkCount = fig.landmarkCount;
#else
	landmarkNames = fig.landmarkNames;
#endif
	markedLandmark = fig.markedLandmark;

	figureStats = fig.figureStats;	// Just copy the pointer
    stackedImageBit = fig.stackedImageBit;
}

M3DFigure & M3DFigure::copy(const M3DFigure & fig)
{
#ifdef DEBUG
	cout << "M3DFigure::copy()" << endl;
#endif
	if (&fig == this)
		return *this;

	if (--properties->counter == 0) {
		if (properties->name != NULL)
			delete [] properties->name;
		delete properties;
	}
    properties = fig.properties;
    properties->counter++;

	ifconstraints = fig.ifconstraints;
	inverse_constraints = fig.inverse_constraints;

	landmarks = fig.landmarks;
#ifdef LM_METHOD_OBJ
	memcpy(landmarkNames, fig.landmarkNames, sizeof(landmarkNames));
	memcpy(landmarkAtomIndices, fig.landmarkAtomIndices, sizeof(landmarkAtomIndices));
	memcpy(landmarkAtomTs, fig.landmarkAtomTs, sizeof(landmarkAtomTs));
	memcpy(landmarkAtomUs, fig.landmarkAtomUs, sizeof(landmarkAtomUs));
	memcpy(landmarkAtomVs, fig.landmarkAtomVs, sizeof(landmarkAtomVs));
	landmarkCount = fig.landmarkCount;
#else
	landmarkNames = fig.landmarkNames;
#endif
	markedLandmark = fig.markedLandmark;

	// FIXME: SHALLOW COPY
	// will cause problem in destructor.
	figureStats = fig.figureStats;	// Just copy the pointer
	stackedImageBit = fig.stackedImageBit;

	return *this;
}

M3DFigure::~M3DFigure()
{
#ifdef DEBUG
	cout << "M3DFigure::~M3DFigure(): this = " << this << endl;
#endif
	if (--properties->counter == 0) {
		if (properties->name != NULL)
			delete [] properties->name;
		delete properties;
	}

	landmarks.clear();
	if (figureStats != NULL)
		delete figureStats;
}

void M3DFigure::print(bool dump, ostream & out, int markedPrimIndex) const
{
    if (properties->name != NULL)
        out << "Figure name: \"" << properties->name << "\"\n";
    else
        out << "Unnamed Figure\n";

    out << "Type: " << typeid(*this).name();
	out << "Visibility: " << getVisibility() << '\t';
	out << "Number of constraints: " << ifconstraints.size() << '\n';
	out << "Number of landmarks: " << getLandmarkCount() << endl;
}

char * M3DFigure::copyName() const
{
    char * retStr;

    if (properties->name == NULL)
        return NULL;

    retStr = new char[strlen(properties->name) + 1];
    strcpy(retStr, properties->name);

    return retStr;
}

void M3DFigure::setName(const char * newName)
{

    if (properties->name != NULL)
        delete [] properties->name;

    if(newName != NULL)
    {
        properties->name = new char[strlen(newName) + 1];
        strcpy(properties->name, newName);
    }
    else
        properties->name = NULL;
}

#ifndef LM_METHOD_OBJ

const char * M3DFigure::getLandmarkName(int index) const
{
	if (3*index >= landmarks.size()) {
		cout << "Invalid landmark index\n";
		return NULL;
	}

    return landmarkNames[index];
}

bool M3DFigure::setLandmarkName(int index, const char * newName)
{
	cout << "ERR: setLandmarkName" << endl;

	char * lname;

	if (3*index >= landmarks.size()) {
		cout << "Invalid landmark index\n";
		return false;
	}

    if(landmarkNames[index] != NULL)
        delete [] landmarkNames[index];

    if(newName != NULL)
    {
        lname = new char[strlen(newName) + 1];
        strcpy(lname, newName);
    }
    else
        lname = NULL;
	landmarkNames[index] = lname;
	return true;
}

#endif

void M3DFigure::setFigureStatsPtr(M3DFigureStats * newFigStatPtr)
{
	figureStats = newFigStatPtr;
}

M3DFigureStats * M3DFigure::getFigureStatsPtr()
{
	return figureStats;
}

void M3DFigure::addLandmark(Vector3D & landmark)
{
#ifndef LM_METHOD_OBJ
	cout << "ERR: addLandmark" << endl;
	// Landmarks are always added to the end of the list
	landmarks.push_back(landmark.getX());
	landmarks.push_back(landmark.getY());
	landmarks.push_back(landmark.getZ());
	landmarkNames.push_back(NULL);
#endif
}

void M3DFigure::translateLandmarks(const Vector3D & vTrans)
{
#ifndef LM_METHOD_OBJ
	for (int i = 0; i < landmarks.size(); i += 3) {
		landmarks[i] += vTrans.getX();
		landmarks[i + 1] += vTrans.getY();
		landmarks[i + 2] += vTrans.getZ();
	}
#endif
}

void M3DFigure::rotateLandmarks(const Quat &q, const Vector3D &vCenter)
{
#ifndef LM_METHOD_OBJ
	for (int i = 0; i < landmarks.size(); i += 3) {
		Vector3D l(landmarks[i], landmarks[i + 1], landmarks[i + 2]);
		l -= vCenter;
		q.rotateVector(l);
		l += vCenter;
		landmarks[i] = l.getX();
		landmarks[i + 1] = l.getY();
		landmarks[i + 2] = l.getZ();
	}
#endif
}

void M3DFigure::scaleLandmarks(double scalefact, const Vector3D &vCenter)
{
#ifndef LM_METHOD_OBJ
	for (int i = 0; i < landmarks.size(); i += 3) {
		Vector3D l(landmarks[i], landmarks[i + 1], landmarks[i + 2]);
		l -= vCenter;
		l *= scalefact;
		l += vCenter;
		landmarks[i] = l.getX();
		landmarks[i + 1] = l.getY();
		landmarks[i + 2] = l.getZ();
	}
#endif
}

#ifdef LM_METHOD_OBJ

Vector3D M3DFigure::getLandmark(int index) {
	Vector3D l;
	return l;
}

#else

Vector3D M3DFigure::getLandmark(int index) {
	if (3*index >= landmarks.size()) {
		Vector3D l(0.0, 0.0, 0.0);
		cout << "Invalid landmark index\n";
		return l;
	}
	Vector3D l(landmarks[3*index], landmarks[3*index + 1], landmarks[3*index + 2]);
	return l;
}

void M3DFigure::clearLandmarks()
{
	if (landmarks.size() == 0)
		return;

	landmarks.clear();
	landmarkNames.clear();
	markedLandmark = -1;
}

bool M3DFigure::clearLandmark(int index)
{
	int i, s;

	s = landmarks.size();
	if (index >= s/3 || index < 0)
		return false;
	if (s == 3) {
		markedLandmark = -1;
		landmarks.clear();
		landmarkNames.clear();
		return true;
	}

	if (index == markedLandmark)
		markedLandmark = -1;
	else if (markedLandmark == s/3 - 1)
		markedLandmark = index;

	// Must retain the landmark ordering, because of pairing with image landmarks
	for (i = 3*index; i < landmarks.size() - 3; i++)
		landmarks[i] = landmarks[i + 3];
	landmarks.pop_back();
	landmarks.pop_back();
	landmarks.pop_back();

	for (i = index; i < landmarkNames.size() - 1; i++)
		landmarkNames[i] = landmarkNames[i + 1];
	landmarkNames.pop_back();

	return true;
}

bool M3DFigure::markLandmark(int index)
{
	int s = landmarks.size();
	if (index >= s/3 || index < 0)
		return false;
	markedLandmark = index;
	return true;
}

#endif

#ifdef LM_METHOD_OBJ

// Assumes that landmark's position is at the tip of the specifed spoke.
// defaults to bisector spoke with elongation
void M3DFigure::addLandmark(int atomIndex, const char *name, double atomT)
{
    if (typeid(*this) == typeid(M3DQuadFigure)) {
	  M3DQuadFigure* quad = dynamic_cast<M3DQuadFigure*>(this);
	  int u, v;
	  // u is col, v is row
	  quad->figuralCoordinates(atomIndex, v, u);
	  addLandmark(u, v, name, atomT);
	}   
    else if (typeid(*this) == typeid(M3DTubeFigure)) {
      cout << "ERROR: M3DFigure::addLandmark [" << name << "] not enough information to specify tube landmark" << endl;
	}
}

void M3DFigure::addLandmark(double atomU, double atomV, const char *name, double atomT)
{
    int atomIndex;
    if (typeid(*this) == typeid(M3DQuadFigure)) {
	  M3DQuadFigure* quad = dynamic_cast<M3DQuadFigure*>(this);
	  // u is col, v is row
	  atomIndex =  quad->indexOfAtom((int) atomV, (int) atomU);
    }
    else if (typeid(*this) == typeid(M3DTubeFigure)) {
	  atomIndex= (int) atomU;
    }

	// Landmarks are always added to the end of the list
	landmarkAtomIndices[landmarkCount] = atomIndex;
	landmarkAtomTs[landmarkCount] = atomT;
	landmarkAtomUs[landmarkCount] = atomU;
	landmarkAtomVs[landmarkCount] = atomV;
	landmarkNames[landmarkCount] = (char*)name;
	landmarkCount++;
}

	
// Find the atomIndex in landmarkAtomIndices[] corresponding to an atom name
// in landmarkNames[]
int	M3DFigure::findLandmarkAtomIndexByName(const char *name)
{
	if (name && *name)
	{
		for (int lm = 0; lm < landmarkCount; lm++)
		{
			char *myName = landmarkNames[lm];
			if (myName && *myName && !strcmp(myName, name))
				return landmarkAtomIndices[lm];
		}
	}
	return -1;
}

// Find the landmarkIndex corresponding to an atom name
// in landmarkNames[]
int	M3DFigure::findLandmarkIndexByName(const char *name)
{
	if (name && *name)
	{
		for (int lm = 0; lm < landmarkCount; lm++)
		{
			char *myName = landmarkNames[lm];
			if (myName && *myName && !strcmp(myName, name))
				return lm;
		}
	}
	return -1;
}

// Find the atom's name in landmarkNames[] corresponding to an atom index
// in landmarkAtomIndices[]
const char *M3DFigure::findLandmarkAtomNameByIndex(const int atomIndex)
{
	if (atomIndex >= 0)
	{
		for (int lm = 0; lm < landmarkCount; lm++)
		{
			int myAtomIndex = landmarkAtomIndices[lm];
			if (myAtomIndex >= 0 && myAtomIndex == atomIndex)
				return landmarkNames[lm];
		}
	}
	return NULL;
}

void M3DFigure::clearLandmarks()
{
	landmarks.clear();
	for (int i=0; i < landmarkCount; i++)
	{
		landmarkNames[i] = NULL;
		landmarkAtomIndices[i] = -1;
		landmarkAtomTs[i] = 0;
	}
	markedLandmark = -1;
	landmarkCount = 0;
}

#endif

void M3DFigure::setStackedImageNumber(int num)
{
	if (num < 0 || num > 15)
		stackedImageBit = 0;
	else
		stackedImageBit = (unsigned short) 0x1 << num;
}


void M3DFigure::applySimilarity(const SimilarityTransform3D & transform,
	const Vector3D & center)
{
    Quat q = transform.getRotation();
    rotateBy(q, center);
    scaleBy(transform.getScale(), center);
    translateBy(transform.getTranslation());
    M3DTubeFigure* tubeFigure	= dynamic_cast<M3DTubeFigure*>(this);
    if( tubeFigure != NULL ) {
        // Special tube phi business for tubes.
        // Get base atom from the tube, and do the phi alignment.
        M3DPrimitive* prim = tubeFigure->getPrimitivePtr( 
            tubeFigure->getBaseAtomIndex() );
        Quat q;
        const double phi = transform.getTubePhi();
        const Vector3D newN = prim->getN() * cos(phi) + 
                              prim->getBPerp() * sin(phi);
        q.buildFromFrame( prim->getB(), newN );
        prim->setQ(q);
        tubeFigure->fixGlobalConsistency();
    }
}


// REGULARITY CALCULATION

// PND is a primitive neighborhood definition applicable only to [Quad|Tube]Meshes right now, 
// eventually we'll extend it to tubes and tri-meshes and whatever else...

double M3DFigure::dist2FromFigure(M3DFigure * figure, DistanceType DT, bool verbose)
{
	// eg, Average Geodesic Distance is calculated as:
	// R = 1/N * Sum over all atoms (geodescic distance (reference atom, candidate atom));

    M3DFigure * figure1 = this;
    M3DFigure * figure2 = figure;

	int nPrims = figure1->getPrimitiveCount();  // Had better be the same in the candidate!
	double dist=0, sum_dist_squared = 0;

	for (int i=0; i<nPrims; i++)
	{
		M3DPrimitive * m1 = figure1->getPrimitivePtr(i);
		M3DPrimitive * m2 = figure2->getPrimitivePtr(i);
		if (! m1->isSelectedForRegularity()) continue;
		//if (! m2->isSelectedForRegularity()) continue;

		if (DT == GEODESIC_DIST || DT == AVE_GEODESIC_DIST)
		{
			dist=m1->atomDistance(m2);
		} else
		{
			// EUCLIDEAN (XYZ) DISTANCE CALCULATION
			Vector3D m1_to_m2 = m2->getX() - m1->getX();
			dist = m1_to_m2.norm();
		}
		sum_dist_squared += dist * dist;
	}

	// Average, if requested
	double result;
	if (DT == AVE_GEODESIC_DIST || DT == AVE_EUCLIDEAN_DIST) 
	{
		double ave_dist_squared = sum_dist_squared / (float) nPrims;
		result = ave_dist_squared;
	}
	else result = sum_dist_squared;


	return result;
}

// Returns the figure such that each atom is the geodesic average of its
// neighbors
M3DFigure* M3DFigure::fromAveOfNeighbors(PrimNeighborhoodDefn PND)
{
	// We'll write into this figure and return it
	M3DFigure * aveFigure = clone();
	M3DPrimitive * atom;
	M3DPrimitive * neighbors[4];  // No atom has more than 4 neighbors
	int numNeighbors;
	M3DPrimitive * new_atom;

	int nPrims = getPrimitiveCount();

	for (int i = 0; i < nPrims; i++) {
		atom = getPrimitivePtr(i);
		new_atom = aveFigure->getPrimitivePtr(i);

		getAtomsNeighbors(i, numNeighbors, neighbors, PND);
		// Cast this figure to its real type and populate these values.
		// This puts the number of neighbors into numNeighbors and puts
		// pointers to them into neighbors

		if (numNeighbors >= 2) {
			// Put new_atom to the predicted position
			new_atom->atomAverage(numNeighbors, neighbors);
		}
		//else {
			// A corner or something with zero or one neighbor, ie, pinned
			//*new_atom	= *atom;	// We have already cloned, redundant.
		//}
	}

	return aveFigure;
}

M3DFigure * M3DFigure::fractionalStepToFigure(M3DFigure * targets[], int weights[],
	int numTargets, DistanceType DT)
{
	GeodesicSym geo; // for internal atoms 
	GeodesicSym geoE(false); // for end atoms 

	// We'll write into this figure and return it
	M3DFigure * aveFigure = clone();
	//M3DQuadFigure * aveFigure = dynamic_cast<M3DQuadFigure*>( assign());

	int numPrims = getPrimitiveCount();

	// Determine total weighting
	int sumWeights = 0;
	for (int i = 0; i < numTargets; i++) {
		sumWeights += weights[i];
	}

	// For each atom
	for (int j = 0; j < numPrims; j++) {
		Vector3D sum_over_targets;
		sum_over_targets.set(0,0,0);

		// GEODESIC AVERAGING
		if (DT == GEODESIC_DIST || DT == AVE_GEODESIC_DIST) {

			M3DPrimitive * atom = aveFigure->getPrimitivePtr(j);
			M3DPrimitive ** targetAtoms = new M3DPrimitive *[sumWeights];

			int offset = 0; // extra atoms from weighting
			for (int i = 0; i < numTargets; i++) {
				for (int k = 0; k < weights[i]; k++) {
					targetAtoms[i+k+offset] = targets[i]->getPrimitivePtr(j);
				}
				offset += weights[i] - 1;  // If it's zero, we need to step back one i
			}

			// Average into the ptr to this atom
			if(atom->type() != M3D_END_PRIMITIVE) {
				geo.atomAverage(sumWeights, targetAtoms, atom);
			}
			else {
				geoE.atomAverage(sumWeights, targetAtoms, atom);
			}
		}
		else {		// EUCLIDEAN AVERAGING
			// Average each target
			for (int i = 0; i < numTargets; i++) {
				// Euclidean (XYZ) distance calculation
				Vector3D this_target_xyz = targets[i]->getPrimitivePtr(j)->getX();

				// Add this weighted target value into the sum
				sum_over_targets += this_target_xyz * weights[i];
			}

			// Set this atom's X vector to the average
			Vector3D ave_over_targets = sum_over_targets / sumWeights;  // Average the coordinates
			aveFigure->getPrimitivePtr(j)->setX(ave_over_targets.getX(),
				ave_over_targets.getY(), ave_over_targets.getZ());
		}
	}
	return aveFigure;
}


// Persistent object support :
M3DFigure* M3DFigure::readFigure(int figureIndex, Registry& registry )
{
    const char * figureTypeStr;

    figureTypeStr = registry.getStringValue("model.figure[%d].type", NULL, figureIndex);
    if(figureTypeStr == NULL) {
        return NULL;
	}
	else if( strcmp(figureTypeStr, M3DTubeFigure::friendlyFigureName) == 0 ) {
		return M3DTubeFigure::readFigure(figureIndex, registry);
    }
    else if( strcmp(figureTypeStr, M3DQuadFigure::friendlyFigureName) == 0 ) {
        return M3DQuadFigure::readFigure(figureIndex, registry);
    }
	else {
		cerr << "Unsupported figure type <" << figureTypeStr << "> being read\n";
    	return NULL;
	}
}

void M3DFigure::resampleForRegularSpacing()
{
	std::cout << typeid(*this).name() << "::resample() is not implemented" << std::endl;
}


