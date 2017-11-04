#ifndef M3DOBJECTFILE_H
#define M3DOBJECTFILE_H

#include "Registry.h"
#include "M3DObject.h"
//#include "Match.h"
#include "M3DCPNSStats.h"
#include "M3DPNSTransform.h"


class WorldSystem;
class Image3D;
class PGData;
class PGSet;
class AugmentedAtoms;


class M3DObjectFile
{
public:

    M3DObjectFile() {}
    ~M3DObjectFile() {}

    M3DObject * read(const char * filename, int & markedPrimitiveId,
		SimilarityTransform3D * xform = NULL, WorldSystem * world = NULL,
		const char * pgaFilename = NULL);
    M3DObject * read(const char * filename, SimilarityTransform3D *
		xform = NULL, WorldSystem * world = NULL,
		const char * pgaFilename = NULL)
	{
		int mark;
		return read(filename, mark, xform, world, pgaFilename);
	}
    bool write(const char * filename, M3DObject & object,
		SimilarityTransform3D * xform = NULL, bool asMatrix = false,
		int markedPrimitiveId = -1, Image3D * image = NULL,
		const char * imageFileName = NULL, bool savePGA = false);

	// dibyendu
	// bool writeDistancesAtSpokeEnds( M3DObject& object, Match * const match, const char* filename ) ;

	Registry * getRegistry() { return & registry; }

protected:
    void writeWorld(const Image3D * image, const char * imageFileName);
	bool readWorld(Vector3D * origin, Vector3D * bound, Vector3D * spacing,
		std::string * imagePath, std::string * imageModTime);

    M3DFigureTreeNode * readFigureTree(const char * regStr, ...);
    void writeFigureTree(const char * regStr, M3DFigureTreeNode * root, ...);

	SubdivBoundary * readBoundary(const char * regStr, ...);
	void writeBoundary(const char * regStr, SubdivBoundary * boundary);

	M3DFigureStats * readFigureStats(const char * regStr, ...);
	void writeFigureStats(const char * regStr, M3DFigureStats * FigureStats);

	M3DPGAStats * readPGAStats();

	// dibyendu - cpns
	M3DCPNSStats *		readCPNSStats() ;	// read CPNSstatistics from the registry
	M3DPNSTransform *	readPNSShape() ;	// read the PNSShape part of CPNSStats
	bool readPNSSpokes( std::vector <M3DPNSTransform *> & _PNSSpoke, int nSpokes ) ;	// read the PNSSpoke[1:n] from CPNSStats
	bool readCPNSEigenVectors( vector <double *> & _eigenVectors, int numEigenmodes, int eigenVectorLen ) ;

	PGData * readPGData(const char * regStr, ...);
	PGSet * readPGSet(const char * regStr, ...);
	AugmentedAtoms * readAugmentedAtoms(const char * regStr, ...);

	PGSet * readPrimitivePGSet(const char * regStr,...);
	M3DPGAPrimitiveStats * readPrimitivePGAStats(M3DObject * object);

	void writeAugmentations(AugmentedAtoms * aug, const char * regStr, ...);
	void writePGSet(PGSet * set, const char * regStr);
	void writePGData(PGData * dataPtr, const char * regStr);
	void writeAtomPGData(PGData * dataPtr, const char * regStr);
	void writePGAStats(M3DPGAStats * pgaStats);

	// dibyendu - cpns
	// void writeCPNSStats( M3DCPNSStats * cpnsStats ) ;

	void writePrimitivePGAStats(M3DPGAPrimitiveStats * pgaStats);

    Registry registry;

private:

	M3DObjectFile(const M3DObjectFile & of);	// Not implemented
	M3DObjectFile & operator=(const M3DObjectFile & of) const;	// Not implemented
};

#endif

