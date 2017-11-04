#include <stdio.h>
#include <stdarg.h>
#include <sys/stat.h>
#include <time.h>
#include <string>

// dibyendu
#include <vector>
#include <fstream>

#include "Registry.h"
#include "Vector3D.h"
#include "Image3D.h"
#include "SimilarityTransform3D.h"
#include "M3DPGAStats.h"
#include "M3DPGAPrimitiveStats.h"

// dibyendu - cpns
#include "M3DCPNSStats.h"

#include "pablo_version.h"
#include "M3DObject.h"
#include "WorldSystem.h"
#include "M3DObjectFile.h"
#include "utility.h"


//#define DEBUG

#include <typeinfo>

#ifdef BINARY
extern int globalVerbosity;
#endif


using namespace std;


static bool skipSettingPGAMean = false;    // Reentry control flag


M3DObject * M3DObjectFile::read(const char * filename, int & markedPrimitiveId,
	SimilarityTransform3D * xform, WorldSystem * world, const char * pgaFilename)
{
    //cout << "M3D Read..." << endl;
    M3DObject * object;
    int numFigures,
        numFigureTrees,
        i;
    M3DFigureTreeNode * treeNode;
	M3DObjectFile pgaObjectFile;
	M3DObject * pgaObject;
	M3DPGAStats * pgaStats;
	M3DPGAPrimitiveStats * pgaAtomStats;

	// dibyendu cpns
	M3DCPNSStats * cpnsStats ;
//    int numFigureLinks;

    if (filename == NULL)
        return NULL;

	if (pgaFilename != NULL) {
		// This leaves the PGA object pointers in pgaObject set to
		// NULL, after the following read is completed, so when it
		// is deleted, they will persist in object.
		skipSettingPGAMean = true;
		pgaObject = pgaObjectFile.read(pgaFilename);
		skipSettingPGAMean = false;
		if (pgaObject == NULL)
			return NULL;
	}
	else
		pgaObject = NULL;

    if (! fileExists(filename))
    {
        cout << "Not a valid file: " << filename << endl;
        return NULL;
    }

	//cout << "Going to read from file" << filename << endl ;

	try {
		registry.readFromFile(filename);
	}
	catch (RException excp) {
		excp.print(cout);
		return NULL;
	}	

	//registry.writeToFile( "mean_model_registry.txt" ) ;	

#undef DEBUG_DIBYENDU_REGISTRY
#ifdef DEBUG_DIBYENDU_REGISTRY

	registry.writeToFile( "../../bin/mean_model_registry.txt" ) ;

#endif 

    object = new M3DObject;
#ifdef BINARY
    object->setName(registry.getStringValue("model.name", "noName"));
#else
    object->setName(registry.getStringValue("model.name", NULL));
#endif

    numFigures = registry.getIntValue("model.figureCount", 0);
    numFigureTrees = registry.getIntValue("model.figureTrees.count", 0);
    markedPrimitiveId = registry.getIntValue("model.markedPrimitiveId", -1);

	bool needToSaveBack	= false;

    // Add figures with default figure trees
    if(numFigureTrees == 0)
    {
        for(i = 0; i < numFigures; i++) {
			M3DFigure* figure	= M3DFigure::readFigure(i,registry);
			if( typeid(*figure) == typeid(M3DTubeFigure) ) {
//				needToSaveBack = needToSaveBack | (dynamic_cast<M3DTubeFigure*>(figure))->orientTubeTangents();
			}
			// FIXME @see read in figure
			figure->setFigureStatsPtr(readFigureStats("model.figure[%d]", i));
            object->addFigure(figure);
		}
    }
    // Add figures without default trees
    else
    {
        for(i = 0; i < numFigures; i++) {
			M3DFigure* figure	= M3DFigure::readFigure(i,registry);
			if( typeid(*figure) == typeid(M3DTubeFigure) ) {
//				needToSaveBack = needToSaveBack | (dynamic_cast<M3DTubeFigure*>(figure))->orientTubeTangents();
			}
			figure->setFigureStatsPtr(readFigureStats("model.figure[%d]",i));
            object->addFigureWithoutTree(figure);
		}
    }

    for(i = 0; i < numFigureTrees; i++)
    {
        treeNode = readFigureTree("model.figureTrees.tree[%d]", i);
        object->addFigureTree(treeNode);
    }



    //-----------------------------------
   /* for(i = 0; i < numFigures; i++) {
        cout<<"Liyun, deltaU is ----:  "<<registry.getFloatValue("model.primitive[%d][%d].deltaU",i)<<endl;
        cout<<"Liyun, deltaV is ----:  "<<registry.getFloatValue("model.primitive[%d][%d].deltaV",i)<<endl;
    }*/






	if (xform != NULL) {
		xform->tie(&registry);
		(void) xform->readSimilarity(NULL);
	}

	// Read PGA data or use that provided in pgaObject
	if (pgaObject == NULL) {
		pgaStats = readPGAStats();
		pgaAtomStats = readPrimitivePGAStats(object);
	}
	else {
		pgaStats = pgaObject->getPGAStats();
		pgaAtomStats = pgaObject->getAtomPGAStats();
	}

	// Dibyendu - cpns
	// Read CPNS data from the registry if it is available in the model file

    //cout << "Preparing to read CPNS stats from registry" << endl ;
	
	cpnsStats = readCPNSStats() ;

    //cout << "CPNS stats read successfully from registry" << endl ;
    //cout << cpnsStats << endl;

#undef DEBUG_DIBYENDU_CPNSSTATS

#ifdef DEBUG_DIBYENDU_CPNSSTATS

	const char* debugModelFilename = "../../bin/test_cpns_eigenDeform.m3d" ;		
	if( cpnsStats != NULL && object != NULL ) {

		double scores[10] ;

		for( int i = 0 ; i < 10 ; i++ )
			// scores[i] = 3 * sqrt( cpnsStats->getEigenValue(i) ) ;	
			scores[i] = 0.0 ;
		
		//M3DObject * newObj11 = cpnsStats->eigenmodeDeformMean( scores ) ;

		M3DObject * newObj11 = cpnsStats->eigenmodeDeformMean( scores, object ) ;

		if( newObj11 == NULL ) 
			cout << "Error in eigenmodeDeformMean( ) ! NULL pointer returned" << endl ;
		else {
			M3DObjectFile mFile1 ;
			mFile1.write( debugModelFilename, *newObj11, 0, 0, -1, 0, 0, 0 ) ;
		}		
	}

#endif 


	// Set the pga pointers in object
	if (pgaStats != NULL)
		object->setPGAStats(pgaStats);
	if (pgaAtomStats != NULL)
		object->setAtomPGAStats(pgaAtomStats);

	// dibyendu - cpns
	// Set the CPNS pointers in object	

	if( cpnsStats != NULL )
		object->setCPNSStats( cpnsStats ) ;

    //cout << "CPNS Stats has been set in the s-rep object" << endl ;

	// These lines must be before setting the mean object
	object->that = object->assign();
	object->that->that = object->that;

	if (! skipSettingPGAMean) {   // Skip, if only reading to get PGA data
		if (pgaStats != NULL) {
			pgaStats->setMeanObj(object->loadedObject());	// object->that
			pgaStats->resetObjsToMean();
		}

		if (pgaAtomStats != NULL)
 			pgaAtomStats->setMeanObj(object->loadedObject());	// object->that

		if (pgaStats != NULL) {

			// JJ : 11/4/05
			// Order matters: first decide what type of model it is,
			// then rescale (b/c the radius)
			if (pgaStats->getPGDataPtr(0)->lenMean > 0) {
				pgaStats->convertMeanDiff();
				object->setModelType(M3DObject::Adaptive);
				object->that->setModelType(M3DObject::Adaptive);
				cout << "Model's PGA statistics are for adaptive optimization\n";
			}
			else {
				object->setModelType(M3DObject::NotAdaptive);
				object->that->setModelType(M3DObject::NotAdaptive);
				cout << "Model's PGA statistics are for non-adaptive optimization\n";				
			}
			pgaStats->rescale();
			pgaStats->type = M3DPGAStats::Scaled;
		}
	}

	if (pgaObject != NULL) {
		// Disconnect the pga variables from the old object, before deleting it
		if (pgaObject->getPGAStats() != NULL)
			pgaObject->setPGAStats(NULL, true);
		if (pgaObject->getAtomPGAStats() != NULL)
			pgaObject->setAtomPGAStats(NULL, true);
		delete pgaObject;
	}

	// This does not set the world in the object, but loads the world variable from
	// the m3d file.
	if (world != NULL) {
		world->clear();
		if (readWorld(&world->origin, &world->bound, &world->spacing, &world->imagePath,
			&world->imageModTime))
		{
			world->genConversions();
			world->ok = true;
		}
	}

	/* A value of -1 for the flipped variable is used in RadOnc to
	   indicate that the world coordinate system is right handed.
	   In Comp. Sci, the value of 1, or no value, defaulting to 0,
	   is used to indicate that the world system has been flipped
	   along the Y-axis, and is thus left handed.
	*/
    int orientation = registry.getIntValue("coordSystem.yDirection", 0);
	if (orientation >= 0)
		object->flipped = false;
	else
		object->flipped = true;
	object->that->flipped = object->flipped;

	if (needToSaveBack)
		std::cout << "Warning: Updates were applied; the model should be saved now";

    //cout << "Object successfully read" << endl ;

    return object;
}

bool M3DObjectFile::write(const char * filename, M3DObject & object,
	SimilarityTransform3D * xform, bool asMatrix, int markedPrimitiveId,
	Image3D * image, const char * imageFileName, bool savePGA)
{
    int numFigures,
        i;
    const char * name;
    char figureStr[1024];
	const char * order[4] = {"pabloVersion", "coordSystem", "model", NULL};

//    M3DFiguralLinksCollection figuralLinks;
//    int numFigureLinks;

    int numFigureTrees;
	M3DPGAStats * pgaStats;
	M3DPGAPrimitiveStats * pgaPrimitiveStats;

    if (filename == NULL)
        return false;

    registry.setStringValue("pabloVersion", revision);

    if (object.orientation())
		registry.setIntValue("coordSystem.yDirection", -1);
	else
		registry.setIntValue("coordSystem.yDirection", 1);

    name = object.getName();
    if (name != NULL)
        registry.setStringValue("model.name", name);

    numFigures = object.getFigureCount();
    registry.setIntValue("model.figureCount", numFigures);

	if (image != NULL)
		writeWorld(image, imageFileName);

    for(i = 0; i < numFigures; i++) {
        sprintf(figureStr, "model.figure[%d]", i);
		M3DFigure* figure	= object.getFigurePtr(i);
        figure->writeFigure(figureStr, registry);
		// FIXME @see M3D*Figure::writeFigure
		// Write boundary information
		//writeBoundary(figureStr, figure->getBoundaryPtr());
	}

//    figuralLinks = object.getFiguralLinks();
//    numFigureLinks = figuralLinks.getLinkCount();

//    registry.setIntValue("model.figuralLinks.count", numFigureLinks);
//    for(i = 0; i < numFigureLinks; i++)
//        writeFiguralLink("model.figuralLinks.figuralLink[%d]", figuralLinks.getLink(i), i);

    numFigureTrees = object.getFigureTreeCount();
    registry.setIntValue("model.figureTrees.count", numFigureTrees);
	if (markedPrimitiveId >= 0)
		registry.setIntValue("model.markedPrimitiveId", markedPrimitiveId);
    for(i = 0; i < numFigureTrees; i++)
        writeFigureTree("model.figureTrees.tree[%d]", object.getFigureTreeRoot(i), i);

	if (xform != NULL) {
		xform->tie(&registry);
		xform->writeSimilarity(NULL, asMatrix);
	}

	if (savePGA) {
		pgaStats = object.getPGAStats();
		if (pgaStats != NULL) 
			writePGAStats(pgaStats);
		
		pgaPrimitiveStats = object.getAtomPGAStats();
	   if (pgaPrimitiveStats != NULL)
			writePrimitivePGAStats(pgaPrimitiveStats);	
	}
	else
		cout << "PGA statistics not saved" << endl;

	registry.ordering(order);
	bool retr = true;
	try {
		if (! registry.writeToFile(filename))
			cout << "Warning: model file may be corrupted" << endl;
	}
	catch (RException excp) {
		cout << excp.message() << endl;
		retr = false;
	}
	registry.clear();

	return retr;
}

// dibyendu
// Function to write the distances at all the spoke ends of an M3D object (quad figures only)

//bool writeDistancesAtSpokeEnds( M3DObject& object, Match * const match, const char* filename ) {
//
//	// Collect the distance data into a vector first
//
//	/* 	If no. of rows = n, no. of cols = m
//
//	Data per object:	
//
//		No. of end atoms = 2 * n  + 2 * (m-2)
//		No. of std atoms = m * n - ( 2 * n  + 2 * (m-2) )
//
//		No. of entities per std atom = 2 
//		No. of entities per end atom = 3 
//
//	*/ 
//
//	int nEndPrimitives = 2 * nRows + 2 * ( nCols - 2 ) ;
//
//	int nStdPrimitives = nRows * nCols - nEndPrimitives ;
//
//	int nEntries = 2 * nStdPrimitives + 3 * nEndPrimitives ;
//
//	// declaring the array of values
//
//	vector <double> distVector ;
//
//	// for every object, take the first figure, get the distances and put them in the vector
//
//	M3DFigure * figThis = object.getFigurePtr(0) ;
//
//	if( ( dynamic_cast <M3DTubeFigure*> (figThis) ) != NULL ) {
//
//		cout << "This would not work with Tube figures" << endl ;
//		return (0) ;
//	}
//
//	int nPrims = figThis->getPrimitiveCount() ;	
//
//	Vector3D thisSpokeEnd ;
//	double   thisDist ;
//
//	for( int ii = 0 ; ii < nPrims ; ii ++ ) {
//
//		M3DPrimitive * primThis = figThis->getPrimitivePtr( ii ) ;
//
//		// spoke 0
//
//		thisSpokeEnd = primThis->getX() + primThis->getY0() ;
//
//		distVector.push_back( match->binaryDistanceMap->getDistance(thisSpokeEnd) ) ;
//
//		// spoke 1
//
//		thisSpokeEnd = primThis->getX() + primThis->getY1() ;
//
//		distVector.push_back( match->binaryDistanceMap->getDistance(thisSpokeEnd) ) ;
//
//		// spoke End
//
//		if( primThis->type() == M3D_END_PRIMITIVE ) {		
//
//			thisSpokeEnd = primThis->getX() + (dynamic_cast <M3DQuadEndPrimitive*> (primThis))->getYEnd() ;
//
//			distVector.push_back( match->binaryDistanceMap->getDistance(thisSpokeEnd) ) ;
//
//		}		
//
//	}
//
//	// open the text file for writing 
//
//	std::ofstream f;
//
//	f.open( filename,  ios::trunc);
//
//	if (!f) {
//		string s("Unable to open " + string(filename) + " for writing spoke end distances!" ) ;
//		cout << s << endl ;
//		return(0) ;
//	}	
//
//	for( int ii = 0 ; ii < distVector.size() ; ii ++ ) {		
//		f << setw(20) << distVector[ii] << endl ;		
//
//	f.close() ;
//
//	cout << filename << " has been written !" << endl ;
//
//	return( 1 ) ;
//
//}



/*	Write the world coordinates of model space into the model file.  This
	can be used to convert atom coordinates to world space as follows:
		world_coords = model_coords*(bound - origin) + origin,
	where each variable is a 3-vector and * is element-by-element
	multiplication.  If available, the voxel spacing, image file name,
	and last modification time are also provided in the output.

	The last argument is assumed to be a full pathname.  On Windows PC's,
	the pathname may contain embedded spaces, which will not be escaped.
*/
void M3DObjectFile::writeWorld(const Image3D * image, const char * imageFileName)
{
    char format[] = "model.world.%s";
	struct stat buf;
	char tbuf[27];

	Vector3D origin = image->getModelOrigin();
	registry.setDoubleValue(format, origin.getX(), "origin.x");
	registry.setDoubleValue(format, origin.getY(), "origin.y");
	registry.setDoubleValue(format, origin.getZ(), "origin.z");

	Vector3D bound = image->getModelBound();
	registry.setDoubleValue(format, bound.getX(), "bound.x");
	registry.setDoubleValue(format, bound.getY(), "bound.y");
	registry.setDoubleValue(format, bound.getZ(), "bound.z");

	registry.setDoubleValue(format, image->getXSpacing(), "spacing.x");
	registry.setDoubleValue(format, image->getYSpacing(), "spacing.y");
	registry.setDoubleValue(format, image->getZSpacing(), "spacing.z");

	if (imageFileName == NULL)
		return;
	if (0 == strlen(imageFileName))
		return;

	// Write the (full) image file path into the model file
	registry.setStringValue(format, imageFileName, "imagePath");

    if (stat(imageFileName, &buf) != 0)
		return;
	strcpy(tbuf, ctime(&buf.st_mtime));
	tbuf[24] = '\0';
	registry.setStringValue(format, tbuf, "imageModTime");
}

/*	Read the world coordinates of model space from the model file.  This
	can be used to place the model in an image having a different world
	coordinate system.
*/
bool M3DObjectFile::readWorld(Vector3D * origin, Vector3D * bound, Vector3D * spacing,
	string * imagePath, string * imageModTime)
{
    char format[] = "model.world.%s";
	double x, y, z;
	const char * str;

	if (! registry.hasKey(format, "origin"))
		return false;

	x = registry.getDoubleValue(format, 0.0, "origin.x");
	y = registry.getDoubleValue(format, 0.0, "origin.y");
	z = registry.getDoubleValue(format, 0.0, "origin.z");
	*origin = Vector3D(x, y, z);

	x = registry.getDoubleValue(format, 0.0, "bound.x");
	y = registry.getDoubleValue(format, 0.0, "bound.y");
	z = registry.getDoubleValue(format, 0.0, "bound.z");
	*bound = Vector3D(x, y, z);

	x = registry.getDoubleValue(format, 0.0, "spacing.x");
	y = registry.getDoubleValue(format, 0.0, "spacing.y");
	z = registry.getDoubleValue(format, 0.0, "spacing.z");
	*spacing = Vector3D(x, y, z);

	str = registry.getStringValue(format, "", "imagePath");
	*imagePath = string(str);

	str = registry.getStringValue(format, "", "imageModTime");
	*imageModTime = string(str);

	return true;
}

/*
M3DFiguralLink * M3DObjectFile::readFiguralLink(const char * regStr, ...)
{
    char newStr[1024];
    char childStr[1024];

    M3DFiguralLink * newLink;
    int i;

    va_list val;
    va_start(val, regStr);
    vsprintf(newStr, regStr, val);
    va_end(val);

    strcpy(childStr, newStr);

    strcat(newStr, ".%s");
    strcat(childStr, ".childLink[%d]");

    newLink = new M3DFiguralLink;

    newLink->setParentId(registry.getIntValue(newStr, -1, "parentId"));

    int childCount = registry.getIntValue(newStr, 0, "childCount");
    for(i = 0; i < childCount; i++)
        newLink->addChildLink(readChildFigureLink(childStr, i));

    return newLink;
}

M3DChildFigureLink * M3DObjectFile::readChildFigureLink(const char * regStr, ...)
{
    char newStr[1024];

    M3DChildFigureLink * childLink;

    va_list val;
    va_start(val, regStr);
    vsprintf(newStr, regStr, val);
    va_end(val);

    strcat(newStr, ".%s");

    childLink = new M3DChildFigureLink;

    childLink->id = registry.getIntValue(newStr, -1, "id");
    childLink->blendExtent = registry.getDoubleValue(newStr, 1.0, "blendExtent");
    childLink->blendAmount = registry.getDoubleValue(newStr, 0.0, "blendAmount");

    return childLink;
}

void M3DObjectFile::writeFiguralLink(const char * regStr, M3DFiguralLink * link, ...)
{
    char newStr[1024];
    char childStr[1024];

    int childCount;
    int i;

    va_list val;
    va_start(val, link);
    vsprintf(newStr, regStr, val);
    va_end(val);

    strcpy(childStr, newStr);

    strcat(newStr, ".%s");
    strcat(childStr, ".childLink[%d]");

    childCount = link->getChildCount();
    registry.setIntValue(newStr, link->getParentId(), "parentId");
    registry.setIntValue(newStr, childCount, "childCount");

    for(i = 0; i < childCount; i++)
        writeChildLink(childStr, link->getChildLink(i), i);
}

void M3DObjectFile::writeChildLink(const char * regStr, M3DChildFigureLink * link, ...)
{
    char newStr[1024];

    va_list val;
    va_start(val, link);
    vsprintf(newStr, regStr, val);
    va_end(val);

    strcat(newStr, ".%s");

    registry.setIntValue(newStr, link->id, "id");
    registry.setIntValue(newStr, link->blendAmount, "blendAmount");
    registry.setIntValue(newStr, link->blendExtent, "blendExtent");
}
*/
M3DFigureTreeNode * M3DObjectFile::readFigureTree(const char * regStr, ...)
{
    char newStr[1024];
    char linkStr[1024];
    char childStr[1024];

    M3DFigureTreeNode * figureTreeNode;

    int linkCount,
        childCount;
    int i;
    M3DPrimitiveLinkInfo * link;


    va_list val;
    va_start(val, regStr);
    vsprintf(newStr, regStr, val);
    va_end(val);

    strcat(newStr, ".%s");

    figureTreeNode = new M3DFigureTreeNode;

    linkCount = registry.getIntValue(newStr, 0, "linkCount");
    childCount = registry.getIntValue(newStr, 0, "childCount");
    figureTreeNode->setFigureId(registry.getIntValue(newStr, -1, "figureId"));
    figureTreeNode->setBlendAmount(registry.getDoubleValue(newStr, 0.0, "blendAmount"));
    figureTreeNode->setBlendExtent(registry.getDoubleValue(newStr, 0.0, "blendExtent"));

    figureTreeNode->setAttachmentMode(
		(M3DFigureTreeNode::SubfigureAttachment_t) registry.getIntValue(newStr,
		(int) M3DFigureTreeNode::UNATTACHED, "attachmentMode"));

    char * str = strrchr(newStr, '.');
    str[0] = '\0';

    for(i = 0; i < linkCount; i++)
    {
        sprintf(linkStr, "%s.link[%d]", newStr, i);
        strcat(linkStr, ".%s");

        link = new M3DPrimitiveLinkInfo;
        link->primitiveId = registry.getIntValue(linkStr, -1, "primitiveId");
        link->u = registry.getDoubleValue(linkStr, 0.0, "u");
        link->v = registry.getDoubleValue(linkStr, 0.0, "v");
        link->t = registry.getDoubleValue(linkStr, 0.0, "t");

        figureTreeNode->addLink(link);
    }

    for(i = 0; i < childCount; i++)
    {
		M3DFigureTreeNode * child;

        sprintf(childStr, "%s.child[%d]", newStr, i);
        child = readFigureTree(childStr);
		if (child->getAttachmentMode() == M3DFigureTreeNode::UNATTACHED)
			child->setAttachmentMode(M3DFigureTreeNode::PROTRUDE);	// Default for subfigures
        figureTreeNode->addChild(child);
    }

    return figureTreeNode;
}

void M3DObjectFile::writeFigureTree(const char * regStr, M3DFigureTreeNode * root, ...)
{
    char newStr[1024];
    char linkStr[1024];
    char childStr[1024];

    int linkCount,
        childCount;
    int i;

    M3DPrimitiveLinkInfo * link;

    va_list val;
    va_start(val, root);
    vsprintf(newStr, regStr, val);
    va_end(val);

    if(root == NULL)
        return;

    strcat(newStr, ".%s");

    linkCount = root->getLinkCount();
    childCount = root->getChildCount();

    registry.setIntValue(newStr, childCount, "childCount");
    registry.setIntValue(newStr, linkCount, "linkCount");
    registry.setIntValue(newStr, root->getFigureId(), "figureId");
    registry.setDoubleValue(newStr, root->getBlendAmount(), "blendAmount");
    registry.setDoubleValue(newStr, root->getBlendExtent(), "blendExtent");

	registry.setIntValue(newStr, (int) root->getAttachmentMode(), "attachmentMode");

    char * str = strrchr(newStr, '.');
    str[0] = '\0';

    for(i = 0; i < linkCount; i++)
    {
        sprintf(linkStr, "%s.link[%d]", newStr, i);
        strcat(linkStr, ".%s");

        link = root->getLink(i);
        if(link == NULL)
            continue;

        registry.setIntValue(linkStr, link->primitiveId, "primitiveId");
        registry.setDoubleValue(linkStr, link->u, "u");
        registry.setDoubleValue(linkStr, link->v, "v");
        registry.setDoubleValue(linkStr, link->t, "t");
    }

    for(i = 0; i < childCount; i++)
    {
        sprintf(childStr, "%s.child[%d]", newStr, i);
        writeFigureTree(childStr, root->getChild(i));
    }
}

// Read boundary displacement info from registry
SubdivBoundary * M3DObjectFile::readBoundary(const char * regStr, ...)
{
	char newStr[1024], baseStr[1024], tmpStr[1024];
	//int level = 2; //  Used to work at level 2 only.  Now works at all levels.

    va_list val;
    va_start(val, regStr);
    vsprintf(baseStr, regStr, val);
    va_end(val);

	SubdivBoundary * bd = new SubdivBoundary();

	for (int level = 0; level < MAX_SUBDIV_LEVEL; level++)
	{
		sprintf(newStr, "%s.subdivLevel[%d]", baseStr, level);
		strcpy(tmpStr, newStr);

		strcat(newStr, ".%s");
		int numPts = registry.getIntValue(newStr, 0, "numPoints");

		if (numPts == 0)	// no boundary information at current level
			continue;

		Displacements * disp = new Displacements(level, numPts);

		strcpy(newStr, tmpStr);
		strcat(newStr, ".vals.dsp[%d]");

		disp->setLevel(level);
		disp->setNumPts(numPts);

		for (int i = 0; i < numPts; i++)
			disp->setVal(i, registry.getDoubleValue(newStr, 0.0, i));

		bd -> setDisplacements(disp);
	}

	return bd;
}

// Write boundary displacement info to registry
void M3DObjectFile::writeBoundary(const char * regStr, SubdivBoundary * boundary)
{
	char newStr[1024], baseStr[1024];
	//int level = 2; // Used to work at level 2 only.  Now works at all levels.

	if (boundary == NULL)
		return;

	for (int level = 0; level < MAX_SUBDIV_LEVEL; level++)
	{
		Displacements * disp = boundary -> getDisplacements(level);

		if (disp == NULL)	// no boundary information available at current level
			continue;

		int numPts = disp -> getNumPts();
		if (numPts == 0)
			continue;

		sprintf(baseStr, "%s.displacements.subdivLevel[%d]", regStr, level);
		strcpy(newStr, baseStr);
		registry.setIntValue(strcat(newStr, ".numPoints"), numPts); 

		double * dVals = disp -> getVals();
		if (dVals == NULL)
			continue;

		strcpy(newStr, baseStr);
		strcat(newStr, ".vals.dsp[%d]");

		for (int i = 0; i < numPts; i++)
			registry.setDoubleValue(newStr, dVals[i], i);
	}
}

// Write the current figure statistics to the m3d file.
M3DFigureStats * M3DObjectFile::readFigureStats(const char * regStr, ...)
{
	char newStr[1024], baseStr[1024], tmpStr[1024];

    va_list val;
    va_start(val, regStr);
    vsprintf(baseStr, regStr, val);
    va_end(val);

	sprintf(newStr, "%s.trainingStats", baseStr);
	strcpy(tmpStr, newStr);

	// Load some singleton information.
	int numPoints = registry.getIntValue("%s.%s", 0, newStr, "numPoints");
	int numTemplates = registry.getIntValue("%s.%s", 0, newStr, "numTemplates");
	int dimension = registry.getIntValue("%s.%s", 0, newStr, "dimension");
	double pvar = registry.getDoubleValue("%s.%s", 0, newStr, "profileMatchVariance");
	double mvar = registry.getDoubleValue("%s.%s", 0, newStr, "meanMatchVariance");

	if (numPoints == 0)
	    return NULL;	// There was no information about templates.

	M3DFigureStats * figstats = new M3DFigureStats;

	// Save the information so far.
	figstats->setnumPoints(numPoints);
	figstats->setnumTemplates(numTemplates);
	figstats->setdimension(dimension);
	figstats->setProfileMatchVar(pvar);
	figstats->setMeanMatchVar(mvar);

	// Now the further nested data.
	int i, j;

	// The templates.
	sprintf(newStr, "%s.templates", tmpStr);
	strcat(newStr, ".type[%d].val[%d]");

	for (i = 0; i < numTemplates; i ++)
	{
		for (j = 0; j < dimension; j ++)
		{
			figstats->setTemplate(i,j, registry.getDoubleValue(newStr, 0.0, i, j));
		}
	}

	// The template types.
	sprintf(newStr, "%s.templateTypes", tmpStr);
	strcat(newStr, ".type[%d]");

	for (i = 0; i < numPoints; i ++)
	{
		figstats->setTemplateType(i, registry.getIntValue(newStr, 0, i));
	}

	// The intensity normalization rmses
	sprintf(newStr, "%s.rmses", tmpStr);
	strcat(newStr, ".val[%d]");

	for (i = 0; i < numPoints; i ++)
	{
		figstats->setRmses(i, registry.getDoubleValue(newStr, 1.0, i));
	}

	// The mean offset (typical profile intensity at a point).
	sprintf(newStr, "%s.meanOffsets", tmpStr);
	strcat(newStr, ".val[%d]");

	for (i = 0; i < numPoints; i ++)
	{
		figstats->setMeanOffset(i, registry.getDoubleValue(newStr, 1.0, i));
	}

	// The mean intensity standard deviations per point.
	sprintf(newStr, "%s.meanStds", tmpStr);
	strcat(newStr, ".val[%d]");

	for (i = 0; i < numPoints; i ++)
	{
		figstats->setMeanStd(i, registry.getDoubleValue(newStr, 1.0, i));
	}

	return figstats;
}

// Write the statistics information into the m3d file, if the data exist.
void M3DObjectFile::writeFigureStats(const char * regStr, M3DFigureStats * figstats)
{
	char newStr[1024], baseStr[1024];

	if (figstats->getnumPoints() == 0) return;

	sprintf(baseStr, "%s.trainingStats", regStr);
	strcpy(newStr, baseStr);

	// Set the singleton information.
	strcat(newStr, ".%s");
	registry.setIntValue(newStr, figstats->getnumPoints(), "numPoints");
	registry.setIntValue(newStr, figstats->getnumTemplates(), "numTemplates");
	registry.setIntValue(newStr, figstats->getdimension(), "dimension");
	registry.setDoubleValue(newStr, figstats->getProfMatchVar(), "profileMatchVariance");
	registry.setDoubleValue(newStr, figstats->getMeanMatchVar(), "meanMatchVariance");

	// Now the nested data.
	int i, j;

	// The templates.
	strcpy(newStr, baseStr);
	strcat(newStr, ".templates.type[%d].val[%d]");

	for (i = 0; i < figstats->getnumTemplates(); i ++)
	{
		for (j = 0; j < figstats->getdimension(); j ++)
			registry.setDoubleValue(newStr, figstats->getTemplateValbyType(i, j), i, j);
	}

	// The template types.
	strcpy(newStr, baseStr);
	strcat(newStr, ".templateTypes.type[%d]");

	for (i = 0; i < figstats->getnumPoints(); i ++)
	{
		registry.setIntValue(newStr, figstats->getTemplateType(i), i);
	}

	// The intensity normalization rmses.
	strcpy(newStr, baseStr);
	strcat(newStr, ".rmses.val[%d]");

	for (i = 0; i < figstats->getnumPoints(); i ++)
	{
		registry.setDoubleValue(newStr, figstats->getRms(i), i);
	}

	// The mean offset (typical profile intensity at a point).
	strcpy(newStr, baseStr);
	strcat(newStr, ".meanOffsets.val[%d]");

	for (i = 0; i < figstats->getnumPoints(); i ++)
	{
		registry.setDoubleValue(newStr, figstats->getMeanOffset(i), i);
	}

	// The mean intensity standard deviations per point.
	strcpy(newStr, baseStr);
	strcat(newStr, ".meanStds.val[%d]");

	for (i = 0; i < figstats->getnumPoints(); i ++)
	{
		registry.setDoubleValue(newStr, figstats->getMeanStd(i), i);
	}
}

// =========================   Start of PGA Functions  =========================

// Read info about augmented atoms in each principal geodesic statistic
AugmentedAtoms * M3DObjectFile::readAugmentedAtoms(const char * regStr, ...)
{
	char newStr[1024];
	int numPrims, i;

	va_list val;
	va_start(val, regStr);
	vsprintf(newStr, regStr, val);
	va_end(val);

	strcat(newStr,".%s");

	AugmentedAtoms * augPtr = new AugmentedAtoms;

	numPrims = registry.getIntValue(newStr, 0, "numPrims");
	augPtr->figIndex = registry.getIntValue(newStr, 0, "figIndex");

	strcat(newStr,"[%d]");

	for(i=0; i < numPrims; i++)
		augPtr->primIndexes.push_back(registry.getIntValue(newStr, 0, "primIndex", i));

	return augPtr;
}

// Read a set(header info) of principal geodesic statistics
PGSet * M3DObjectFile::readPGSet(const char * regStr, ...)
{
	char newStr[1024];
	int numFigs, numAugs, i;

	PGSet * setPtr = new PGSet;

	va_list val;
	va_start(val, regStr);
	vsprintf(newStr, regStr, val);
	va_end(val);

	strcat(newStr,".%s");

	numFigs = registry.getIntValue(newStr,-1, "numFigs");
	numAugs = registry.getIntValue(newStr,-1, "numAugs");

	//NOT Found
	if (numFigs == -1) {
		delete setPtr;
		return NULL;
	}

	setPtr->setName(registry.getStringValue(newStr, NULL, "name"));
	strcat(newStr,"[%d]");

	for (i=0; i < numFigs; i++)
	{
		int temp = registry.getIntValue(newStr, 0, "figIndex", i);
		setPtr->figIndexes.push_back(registry.getIntValue(newStr, 0, "figIndex", i));
	}

	for (i=0; i< numAugs; i++)
		setPtr->augmentations.push_back(readAugmentedAtoms(newStr, "AugmentedAtoms", i));

	return setPtr;
}

// Read principal geodesic data corresponding to setPtr
PGData * M3DObjectFile::readPGData(const char * regStr, ...)
{	
	char newStr[1024];
	int lenPG, numPGs, lenMean, i;

	va_list val;
	va_start(val, regStr);
	vsprintf(newStr, regStr, val);
	va_end(val);

	strcat(newStr,".%s");

	PGData *pg = new PGData();

	pg->meanRes = registry.getDoubleArray(newStr, &lenMean, "mean");
	pg->lenMean = lenMean;

#ifdef DEBUG
	cerr << "length of mean is " << lenMean << " (in readPGData())" << endl;
#endif

    pg->lenPG = lenPG = registry.getIntValue(newStr, -1, "PGLength");
	pg->numPGs = numPGs = registry.getIntValue(newStr, -1, "numPGs");

	//length of mean and PG should be the same!!
	if(lenMean != lenPG && lenMean != 0)
		cerr << "Error : Length of PG is not equal to the length of mean." << endl;

	//Not FOUND
	if (lenPG == -1) {
		delete pg;
		return NULL;
	}

	strcat(newStr, "[%d]");
	for(i = 0; i < numPGs; i++)
    {
		int PGLength;
		double *pgVec= registry.getDoubleArray(newStr, &PGLength, "PG", i);
		if (PGLength != lenPG) {
			cerr << "Error : Length of PG does not match.\n" << endl;
			return pg;
		}
		pg->pgVec.push_back(pgVec);
    }

	return pg;
}

// Read a file that contains principal geodesic statistics
M3DPGAStats * M3DObjectFile::readPGAStats()
{
	int numSets, i;
	M3DPGAStats * pgaStats;

	numSets = registry.getIntValue("PGAStats.PGSets.numStats", -1);
	if (numSets < 0) {
		//No statistics
		return NULL;
	}

	pgaStats = new M3DPGAStats;

	for (i = 0; i < numSets; i++) {
		// Should be done in this order
		// 1. get the header info
		// 2. get the pg data corresponding to the header
		PGSet * setPtr = readPGSet("PGAStats.PGSets.Set[%d]", i);
		if (setPtr != NULL)
			pgaStats->PGSets.push_back(setPtr);

		PGData * pgPtr = readPGData("PGAStats.Set[%d]", i);

		if (pgPtr != NULL)
			pgaStats->PGs.push_back(pgPtr);

		if (i > 0 && i < numSets - 1) {
			pgPtr = readPGData("PGAStats.Set[%d].Prediction", i);
			if (pgPtr != NULL)
				pgaStats->Predictions.push_back(pgPtr);
		}
	}

	int type = registry.getIntValue("PGAStats.type", 0);
	if (type == 1)	// mean(r)-scaled
		pgaStats->type = M3DPGAStats::Scaled;
	else	// Not mean(r)-scaled
		pgaStats->type = M3DPGAStats::NotScaled;

	int scaled = registry.getIntValue("PGAStats.scaled", 0);
	if (scaled == 1) // mean(r) scaled
		pgaStats->type = M3DPGAStats::Scaled;
	else    // Not mean(r) scaled
		pgaStats->type = M3DPGAStats::NotScaled;

#ifdef DEBUG
	cerr << "name of object = " << pgaStats->meanObj->getName() << endl;
	cerr << "num of figures = " << pgaStats->meanObj->getFigureTreeCount() << endl;
#endif

	return pgaStats;
}

// dibyendu
// Read CPNS infomation from the registry
M3DCPNSStats * M3DObjectFile::readCPNSStats() {

	//cout << "Reading CPNS..." << endl;

	M3DCPNSStats * cpnsStats ;	

	int numEigenmodes	= registry.getIntValue("CPNSStats.nEigenmodes", -1) ;
	int numSpokes		= registry.getIntValue("CPNSStats.nSpokes", -1) ;
	int eigenVectorLen  = registry.getIntValue("CPNSStats.eigenVectorLength", -1) ;

	int numRows			= registry.getIntValue("CPNSStats.nAtomRows", -1) ;
	int numCols			= registry.getIntValue("CPNSStats.nAtomCols", -1) ;

	// CPNS statistics are not in this model file
	if( numEigenmodes < 0 || numSpokes < 0 || eigenVectorLen < 0 || numRows < 0 || numCols < 0 )
		return NULL ;	

	// allocate space for a new M3DCPNSStats object

	cpnsStats = new M3DCPNSStats( numSpokes, numEigenmodes, eigenVectorLen, numRows, numCols ) ;

	// -----------------  Read scaleShape --------------------------------	

	double _scaleShape	= registry.getDoubleValue("CPNSStats.scaleShape", -1) ;

	if( ! cpnsStats->setScaleShape( _scaleShape ) ) {
		cout << "Error in M3DObjectFile::readCPNSStats() ! scaleShape incorrect" << endl ;
		delete cpnsStats ; 
		cpnsStats = NULL ;
		return NULL ;
	}

	cout << "scaleShape set" << endl ;

	// -----------------  Read meanShape --------------------------------	

	int mLen = 3 ;

	double * _meanShape = registry.getDoubleArray("CPNSStats.meanShape", &mLen ) ;

	if( ! cpnsStats->setMeanShape( _meanShape ) ) {
		cout << "Error in M3DObjectFile::readCPNSStats() ! meanShape NULL pointer" << endl ;
		delete cpnsStats ; 
		cpnsStats = NULL ;
		return NULL ;
	}	
	
	cout << "meanShape set" << endl ;

	// -----------------  Read PNSShape --------------------------------

#define DEBUG_DIBYENDU_CPNS_READ

#ifdef DEBUG_DIBYENDU_CPNS_READ

	
	M3DPNSTransform * _PNSShape = readPNSShape() ;	// donot delete this - its owned by cpnsStats

	// cout << "PNSShape read successfully" << endl ;

	if( _PNSShape == NULL ) {
		cout << "Error in M3DObjectFile::readCPNSStats() ! PNSShape not read correctly." << endl ;
		delete cpnsStats ; 
		cpnsStats = NULL ;
		return NULL ;
	}
	else {
		//cout << "Attempting to set PNSShape" << endl ;
		cpnsStats->setPNSShape( _PNSShape ) ;
		// delete _PNSShape ;
		_PNSShape = NULL ;
	}

	cout << "PNSShape set" << endl ;

#endif // DEBUG_DIBYENDU_CPNS_READ


	
	// -----------------  Read PNSSpokes --------------------------------

	vector <M3DPNSTransform *> _PNSSpoke ;
	
	

	if( ! readPNSSpokes( _PNSSpoke, numSpokes ) ) {
		cout << "Error in M3DObjectFile::readCPNSStats() ! PNSSpoke[] not read correctly." << endl ;
		delete cpnsStats ; 
		cpnsStats = NULL ;
		return NULL ;
	}
	else {
		if( ! cpnsStats->setPNSSpoke( _PNSSpoke ) ) {
			cout << "Error in M3DObjectFile::readCPNSStats() ! PNSSpoke[] not set correctly." << endl ;
			delete cpnsStats ; 
			cpnsStats = NULL ;
			return NULL ;
		}
	}

	cout << "PNSSpoke set" << endl ;

	// -----------------  Read scaleSpokes --------------------------------

	double * _scaleSpoke = registry.getDoubleArray("CPNSStats.scaleSpoke", &numSpokes ) ;

	if( _scaleSpoke == NULL ) {
		cout << "Error in M3DObjectFile::readCPNSStats() ! scaleSpoke not read correctly" << endl ;
		delete cpnsStats ; 
		cpnsStats = NULL ;
		return NULL ;
	}
	else if( ! cpnsStats->setScaleSpoke( _scaleSpoke ) ) {
		cout << "Error in M3DObjectFile::readCPNSStats() ! scaleSpoke NULL pointer" << endl ;
		delete cpnsStats ; 
		cpnsStats = NULL ;
		return NULL ;
	}	

	cout << "scaleSpoke set" << endl ;

	// -----------------  Read eigenValues --------------------------------

	double * _eigenValue = registry.getDoubleArray("CPNSStats.eigenValue", &numEigenmodes ) ;

	if( _eigenValue == NULL ){ 
		cout << "Error in M3DObjectFile::readCPNSStats() ! eigenValues not read correctly from registry" << endl ;
		delete cpnsStats ; 
		cpnsStats = NULL ;
		return NULL ;
	}

	cout << "eigenValues set" << endl ;

	// -----------------  Read eigenVectors --------------------------------

	vector <double *> _eigenVector ;

	if( ! readCPNSEigenVectors( _eigenVector, numEigenmodes, eigenVectorLen ) ) {
		cout << "Error in M3DObjectFile::readCPNSStats() ! eigenVector not read correctly from registry" << endl ;
		delete cpnsStats ; 
		cpnsStats = NULL ;
		return NULL ;
	}
	else if( ! cpnsStats->setEigenModes( _eigenValue,_eigenVector ) ) {		
		cout << "Error in M3DObjectFile::readCPNSStats() ! eigenModes not set correctly" << endl ;
		delete cpnsStats ; 
		cpnsStats = NULL ;
		return NULL ;
	}	

	cout << "eigenVectors set" << endl ;

	// -----------------  write CPNS Stats to an ASCII file --------------------------------

#undef DEBUG_DIBYENDU_CPNS_ASCII

#ifdef DEBUG_DIBYENDU_CPNS_ASCII

	if( cpnsStats->writeToFile( "../../bin/cpns_output.txt" ) )
		cout << "CPNS stats has been written to ../../bin/cpns_output.txt" << endl ;

#endif 

	//cout << "Preparing to return cpnsStats from M3DObjectFile::readCPNSStats" << endl ;

	cout << "CPNS Read" << endl;

#ifdef DEBUG_DIBYENDU_CPNS_READ

	return cpnsStats ;
	

#else // DEBUG_DIBYENDU_CPNS_READ

	return NULL ;

#endif // DEBUG_DIBYENDU_CPNS_READ

	return cpnsStats;

}

M3DPNSTransform * M3DObjectFile::readPNSShape() {

	//cout << "M3DObjectFile::readPNSShape started" << endl ;

	int nDimsThis = registry.getIntValue( "CPNSStats.PNSShape.nDims", -1 ) ;

	if( nDimsThis < 0 )
		return NULL ;

	// ---------------  set nDims  ------------------------------

	M3DPNSTransform * pnsShape = new M3DPNSTransform( nDimsThis ) ;	

	//cout << "new PNSShape allocated" << endl ;

	// ---------------  set sphereAxis (v[1], v[2] ... v[d]) ------------------------------

	char fieldName[100] ; 

	int lenThis ; 

	double * sphereAxisThis = NULL ;


	for( int nd = 0 ; nd < nDimsThis ; nd++ ) {	

		if( nd == (nDimsThis-1) )	// length is 1 for the last axis (instead of being 2)
			lenThis = 1 ;
		else
			lenThis = ( nDimsThis + 1 ) - nd ;

		sprintf( fieldName, "CPNSStats.PNSShape.sphereAxis[%d]", nd) ;		

		// cout << fieldName << nd << endl ;
	
		sphereAxisThis = registry.getDoubleArray( (const char*) fieldName, &lenThis ) ;		

		//cout << "sphereAxisThis[" << nd << "] has been read from registry successfully" << endl ;

		if( sphereAxisThis != NULL ) {

			if( ! pnsShape->setSphereAxis( nd, sphereAxisThis, lenThis ) ) {
				cout << "Error setting sphereAxis[" << nd << "] in M3DObjectFile::readPNSShape()" << endl ;
				delete [] sphereAxisThis ;
				sphereAxisThis = NULL ;
				return NULL ;
			}

			delete [] sphereAxisThis ;		

			sphereAxisThis = NULL ;
		}
		else {	
			cout << "Error reading sphereAxis[" << nd << "] from registry in M3DObjectFile::readPNSShape()" << endl ;
			return NULL ;
		}

	}

	//cout << "sphereAxis has been read" << endl ;

	// ---------------  set sphereDist (r[1], r[2] ... r[d-1]) ------------------------------

	int aLength = nDimsThis-1 ;

	double * _sphereDist = registry.getDoubleArray( "CPNSStats.PNSShape.sphereDist", &aLength ) ;

	if( _sphereDist != NULL ) {

		if( ! pnsShape->setSphereDist( _sphereDist, aLength ) ) {
			delete [] _sphereDist ;
			return NULL ;
		}

		delete [] _sphereDist ;
		_sphereDist = NULL ;
	}
	else {
		cout << "Error reading sphereDist from registry in M3DObjectFile::readPNSShape()" << endl ;
		return NULL ;
	}

	//cout << "sphereDist has been read" << endl ;

	//cout << "M3DObjectFile::readPNSShape completed successfully" << endl ;

	return pnsShape ;

	//return NULL ;

}

bool M3DObjectFile::readPNSSpokes( std::vector <M3DPNSTransform *> & PNSSpoke, int nSpokes ) {

	
	PNSSpoke.reserve( nSpokes ) ;

	// all 2-D spoke PNS transformations

	M3DPNSTransform * pnsSpokeThis = NULL ;

	double * _sphereDist = NULL ;

	double * sphereAxisThis = NULL ;

	char fieldName[100] ;
	
	int lenSphereDist = 1 ;

	for( int ns = 0 ; ns < nSpokes ; ns ++) {

		pnsSpokeThis = new M3DPNSTransform( 2 ) ;

	

		// ---------------  set sphereDist (r[1]) ------------------------------

		sprintf( fieldName, "CPNSStats.PNSSpoke[%d].sphereDist", ns) ;

		_sphereDist = registry.getDoubleArray( (const char*)fieldName, &lenSphereDist ) ;

		if( _sphereDist != NULL ) {

			if( ! pnsSpokeThis->setSphereDist( _sphereDist, lenSphereDist ) ) {
				cout << "Error in M3DObjectFile::readPNSSpokes !" << endl ;
				cout << "Error in M3DPNSTransform::setSphereDist() for spoke no." << ns << endl ;
				delete [] _sphereDist ;
				PNSSpoke.clear() ;
				return 0 ;
			}

			delete [] _sphereDist ;
			_sphereDist = NULL ;
		}
		else {
			cout << "Error reading sphereDist for spoke no." << ns << " from registry in M3DObjectFile::readPNSSpokes()" << endl ;
			PNSSpoke.clear() ;
			return 0 ;
		}		
		
		// ---------------  set sphereAxis (v[1], v[2]) ------------------------------		
		
		

		for( int i = 0 ; i < 2 ; i++ ) {			

			int lenThis = 3 ;

			if( i == 1 )		// length is 1 for the last axis (instead of being 2)
				lenThis = 1 ;			
			
			sprintf( fieldName, "CPNSStats.PNSSpoke[%d].sphereAxis[%d]", ns, i) ;			

			sphereAxisThis = registry.getDoubleArray( (const char*) fieldName, &lenThis ) ;

			if( sphereAxisThis != NULL ) {

				if( ! pnsSpokeThis->setSphereAxis( i, sphereAxisThis, lenThis ) ) {
					cout << "Error setting sphereAxis[" << i << "] in spoke no." << ns << " in M3DObjectFile::readPNSSpokes()" << endl ;
					delete [] sphereAxisThis ;
					sphereAxisThis = NULL ;
					PNSSpoke.clear() ;
					return(0) ;
				}
				delete [] sphereAxisThis ;
				sphereAxisThis = NULL ;
			}
			else {	
				cout << "Error reading sphereAxis[" << i << "] in spoke no. " << ns << " from registry in M3DObjectFile::readPNSSpokes()" << endl ;
				PNSSpoke.clear() ;
				return(0) ;
			}
		}

		PNSSpoke.push_back( pnsSpokeThis ) ;

		pnsSpokeThis = NULL ;

	} // end of ns (spoke) loop
	return(1) ;
}


bool M3DObjectFile::readCPNSEigenVectors( vector <double *> & _eigenVectors, int numEigenmodes, int eigenVectorLen ) {

	char fieldName[100] ;

	double * thisEigenVector ;

	for( int i = 0 ; i < numEigenmodes ; i ++ ) {

		sprintf( fieldName, "CPNSStats.eigenVector[%d]", i) ;			

		thisEigenVector = registry.getDoubleArray( (const char*) fieldName, & eigenVectorLen ) ;

		if( thisEigenVector == NULL ) {
			cout << "Error in M3DObjectFile::readCPNSEigenVectors() ! eigenVector[" << i << "] could not be read from registry" << endl ;
			return(0) ;
		}

		_eigenVectors.push_back( thisEigenVector ) ;

		thisEigenVector = NULL ;

	}

	return(1) ;

}
void M3DObjectFile::writeAugmentations(AugmentedAtoms * aug, const char * regStr, ...)
{
	char newStr[1024];
	int numPrims, i;

	va_list val;
	va_start(val, regStr);
	vsprintf(newStr, regStr, val);
	va_end(val);

	strcat(newStr, ".%s");

	registry.setIntValue(newStr, aug->figIndex, "figIndex");
	numPrims = aug->getNumPrims();
	registry.setIntValue(newStr, numPrims, "numPrims");

	strcat(newStr, "[%d]");
	for (i = 0; i < numPrims; i++)
		registry.setIntValue(newStr, aug->primIndexes[i], "primIndex", i);
}

/*
	Write the residue statistics of multi-object that pgaStat has into a file 
*/
void M3DObjectFile::writePGSet(PGSet * set, const char * regStr)
{
	char newStr[1024];
	int numFigs, numAugs, i;

	strcpy(newStr, regStr);
	strcat(newStr, ".%s");

	numFigs = set->getNumFigs();
	registry.setIntValue(newStr, numFigs, "numFigs");
	numAugs = set->getNumAugs();
	registry.setIntValue(newStr, numAugs, "numAugs");
	if(set->name != NULL)
		registry.setStringValue(newStr, set->name, "name");

	strcat(newStr, "[%d]");
	for (i = 0; i < numFigs; i++)
		registry.setIntValue(newStr, set->figIndexes[i], "figIndex", i);

	for (i = 0; i < numAugs; i++)
		writeAugmentations(set->augmentations[i], newStr, "AugmentedAtoms", i);
}

void M3DObjectFile::writePGData(PGData * dataPtr, const char * regStr)
{
	char newStr[1024];
	int numPGs, PGlength, i;

	if (dataPtr == NULL)
		return;

	strcpy(newStr, regStr);
	strcat(newStr, ".%s");

	numPGs = dataPtr->numPGs;
	PGlength = dataPtr->lenPG;
	registry.setIntValue(newStr, numPGs, "numPGs");
	registry.setIntValue(newStr, PGlength, "PGLength");

	if (dataPtr->meanRes != NULL) {
		int lenMean = dataPtr->lenMean;
		double * mean  = new double[lenMean];
		for (i = 0; i < lenMean; i++)
			mean[i] = dataPtr->meanRes[i];

		registry.setDoubleArray(newStr, lenMean, mean, "mean");
	}

	strcat(newStr, "[%d]");

	for (i = 0; i < numPGs; i++) {
		if (dataPtr->pgVec[i] != NULL) {
			double * pgVec  = new double[PGlength];
			for (int j = 0; j < PGlength; j++)
				pgVec[j] = (dataPtr->pgVec[i])[j];

			registry.setDoubleArray(newStr, PGlength, pgVec, "PG", i);
		}
	}
}

void M3DObjectFile::writeAtomPGData(PGData * dataPtr, const char * regStr)
{
	char newStr[1024];

	int numPGs, PGlength, i;

	if (dataPtr == NULL)
		return;

	strcpy(newStr, regStr);
	strcat(newStr, ".%s");



	numPGs = dataPtr->numPGs;
	PGlength = dataPtr->lenPG;
	registry.setIntValue(newStr, numPGs, "numPGs");
	registry.setIntValue(newStr, PGlength, "PGLength");

	if (dataPtr->meanRes != NULL) {
		int lenMean = dataPtr->lenMean;
		double * mean  = new double[lenMean];
		for (i = 0; i < lenMean; i++)
			mean[i] = dataPtr->meanRes[i];

		registry.setDoubleArray(newStr, lenMean, mean, "mean");
	}

	strcat(newStr, "[%d]");

	for (i = 0; i < numPGs; i++) {
		if (dataPtr->pgVec[i] != NULL) {
			double * pgVec  = new double[PGlength];
			for (int j = 0; j < PGlength; j++)
				pgVec[j] = (dataPtr->pgVec[i])[j];

			registry.setDoubleArray(newStr, PGlength, pgVec, "PG", i);
		}
	}
}

void M3DObjectFile::writePGAStats(M3DPGAStats * pgaStats)
{
	int i, numSets;
    char setStr[1024];

    if (pgaStats == NULL)
		return;

	registry.setIntValue("PGAStats.scaled", (int) pgaStats->isScaled());

	numSets = pgaStats->getNumOfPGSets();
    registry.setIntValue("PGAStats.PGSets.numStats", numSets);
    for(i = 0; i < numSets; i++) {
        sprintf(setStr, "PGAStats.PGSets.Set[%d]", i);
        writePGSet(pgaStats->getPGSetPtr(i), setStr);

		sprintf(setStr, "PGAStats.Set[%d]", i);
		writePGData(pgaStats->getPGDataPtr(i), setStr); 

		if (i > 0 && i < numSets - 1) {
			sprintf(setStr, "PGAStats.Set[%d].Prediction", i);
			writePGData(pgaStats->getPredictionPtr(i - 1), setStr);
		}
	}
}

// Read both mean difference and the PGAs for atoms
M3DPGAPrimitiveStats * M3DObjectFile::readPrimitivePGAStats(M3DObject * object)
{
	M3DPGAPrimitiveStats * pgaPrimitiveStats; 
	int numFigs, figNum, figureId;
    int numAtoms, atomNum, atomId;

	numFigs = registry.getIntValue("atomPGAStats.numFigs", -1);
	if (numFigs == -1)
		return NULL;

    pgaPrimitiveStats = new M3DPGAPrimitiveStats;
	pgaPrimitiveStats->setNumFigs(numFigs);

	for (figNum = 0; figNum< numFigs; figNum++) {
		figureId = registry.getIntValue("atomPGAStats.figure[%d].figIndex", -1, figNum);
	
		
		M3DFigure * figure = object->getFigurePtr(figureId);

        // The figure with the "figureID" has "numAtoms[figNum]" atoms
        numAtoms = registry.getIntValue("atomPGAStats.figure[%d].numAtoms", -1, figNum);
		pgaPrimitiveStats->setNumAtoms(numAtoms);

		M3DPrimitive * primitivePtr;
		for (atomNum = 0; atomNum < numAtoms; atomNum++) {
		    atomId = registry.getIntValue("atomPGAStats.figure[%d].deltaAtom[%d]",
				figNum, atomNum);

			if( typeid(*figure) == typeid(M3DQuadFigure) ){
				primitivePtr = M3DQuadPrimitive::readPrimitive(registry, "atomPGAStats.figure[%d].deltaAtom[%d].deltaMeanAtom", figNum, atomNum);
			}

			if( typeid(*figure) == typeid(M3DTubeFigure) ){
				int spokeNum = dynamic_cast<const M3DTubeFigure*>(figure)->getNumberOfSpokes();
				primitivePtr = M3DTubePrimitive::readPrimitive(registry, spokeNum, "atomPGAStats.figure[%d].deltaAtom[%d].deltaMeanAtom", figNum, atomNum);
			}
			
			pgaPrimitiveStats->setDeltaMeanPrimitives(primitivePtr);
			
			PGData * pgPtr = readPGData(
				"atomPGAStats.figure[%d].deltaAtom[%d].deltaAtomPGAStats", figNum, atomNum);

			if (atomNum == 0) {	// Ensemble PG stats 
				// Non adaptive case
				if (pgPtr->lenMean != 0) {
					cout << "Warning: Model's Atom-PGA statistics are for adaptive optimization:\n";
					cout << "    Discarding Atom-PGA statistics\n";
					delete pgaPrimitiveStats;
					delete pgPtr;
					return NULL;
				}
			}
			if (pgPtr != NULL){
				pgPtr->atomIndex = atomId;
				pgPtr->figureIndex = figureId;
				pgaPrimitiveStats->setPGs(pgPtr);
			}
			else {
				delete pgaPrimitiveStats;
				return NULL;
			}		  
		}
	}
	return pgaPrimitiveStats;
}

void M3DObjectFile::writePrimitivePGAStats(M3DPGAPrimitiveStats * atomPgaStats)
{
	int numFigs, figIndex;
	int numAtoms, atomIndex;
    char setStr[1024];
	char tempStr[1024];
	char temp1Str[1024];
	PGData * dataPtr;

    if (atomPgaStats == NULL)
		return;

	numFigs = atomPgaStats->getNumFigs();
    registry.setIntValue("atomPGAStats.numFigs", numFigs);

	for (int i = 0; i < numFigs; i++) {        
		numAtoms = (atomPgaStats->getNumAtoms())[i];

		figIndex = (atomPgaStats->getPGDataPtr(i,0))->figureIndex;				
		registry.setIntValue("atomPGAStats.figure[%d].figIndex", figIndex, i);
        registry.setIntValue("atomPGAStats.figure[%d].numAtoms", numAtoms, i);

		for (int j = 0;j < numAtoms; j++) {			
			sprintf(setStr, "atomPGAStats.figure[%d].deltaAtom[%d]", i, j);

			dataPtr = atomPgaStats->getPGDataPtr(i, j);			
			atomIndex = dataPtr->atomIndex;

			strcpy(tempStr, setStr);
			strcat(tempStr, ".atomIndex");
			registry.setIntValue(tempStr,atomIndex);

			if (	(atomPgaStats->getDeltaMeanPrimitives())[i * numAtoms + j] !=NULL){
				strcpy(temp1Str, setStr);
				strcat(temp1Str, ".deltaMeanAtom");
				(atomPgaStats->getDeltaMeanPrimitives())[i * numAtoms + j]->writePrimitive(registry,temp1Str);
			}
			strcat(setStr, ".deltaAtomPGAStats");
			writeAtomPGData(dataPtr, setStr); 
		}		
	}
}

PGSet * M3DObjectFile::readPrimitivePGSet(const char *regStr, ...)
{
	char newStr[1024];	
	PGSet * setPtr = new PGSet;
	va_list val;

	va_start(val, regStr);
	vsprintf(newStr, regStr, val);
	va_end(val);

	strcat(newStr,".%s");		

	setPtr->setName(registry.getStringValue(newStr, NULL, "name"));

	strcat(newStr,"[%d]");		

	return setPtr;	
}

// =========================   End of PGA Functions  =========================





