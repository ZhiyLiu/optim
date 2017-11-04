#ifndef M3D_PGA_H
#define M3D_PGA_H

#include "Vector2D.h"
#include "SimilarityTransform3D.h"
#include "M3DObject.h"


extern const int NUM_PGA_PARAMS;

struct PGAData
{
    PGAData()
    {
        mean = NULL;
        v = NULL;
        numVectors = 0;
        vectorSize = 0;
    }

    ~PGAData()
    {
        if(mean != NULL)
            delete mean;
        if(v != NULL)
            delete v;
    }

    int numVectors;
    int vectorSize;
    M3DObject * mean;
    double * v;
};

class M3DPGA
{
public:

	M3DPGA(M3DObject * object = NULL);

    void initializeMean(M3DObject * object);
#ifndef BINARY
    bool readPGAFile(const char * pgaFilename);
#endif
	virtual bool pgaDeform(M3DObject * object, std::vector<double> & vals);

    int componentCount() const { return pgaData.numVectors; }

	SimilarityTransform3D & getHandPlacementTransform() {
		return handPlacementTransform;
	}
	void setHandPlacementTransform(const SimilarityTransform3D & transform)
	{
		handPlacementTransform = transform;
	}

protected:

	/**
	 * @deprecated
	 * This function is deprecated. DO NOT USE anymore.
	 */
	void symToLieAtom(Quat & q, double & theta, Vector3D n0, Vector3D n1);

    PGAData pgaData;
	SimilarityTransform3D handPlacementTransform;

};

class AugmentedAtoms 
{
	public:
		AugmentedAtoms() { figIndex = 0; }
		~AugmentedAtoms() {}

		int getNumPrims() { return primIndexes.size(); }	// Atoms that are augmented in the figure of figIndex
		int figIndex;

		std::vector<int> primIndexes;
};

// Header information for each set of pg's
class PGSet
{
	public:
		PGSet() {
			order = 0; 
			name = NULL;
		}

		~PGSet() {
			int numAugs = augmentations.size();
			for (int i = 0; i < numAugs; i++) {
				if (augmentations[i] != NULL)
					delete augmentations[i];
				augmentations[i] = NULL;
			}
			augmentations.clear();

			if (name != NULL)
				delete [] name;
		}

		int getNumFigs() { return figIndexes.size(); };  // Number of figures
		int getNumAugs() { return augmentations.size(); }; // Number of augmentations
		void setName(const char * newName) {	
			if (name != NULL)
				delete [] name;
			if (newName != NULL) {
				name = new char[strlen(newName) + 1];
				strcpy(name, newName);
			}
			else
				name = NULL;
		};

		int order;	// Order of PG applied in optimization stage : ensemble(0) object(1..N)
		char * name;

		//std::vector<int> primitiveIndexes;
		std::vector<int> figIndexes;
		std::vector<AugmentedAtoms *> augmentations;	// Augmented atoms 
};

// Actual principal geodesics of corresponding PG header
class PGData
{
	public:
		PGData() {
			numPGs = 0;
			lenPG = 0;
			lenMean = 0;
			meanRes = NULL;
		}

		~PGData() {
			for (int i = 0; i < numPGs; i++) {
				if (pgVec[i] != NULL)
					delete [] pgVec[i];
				pgVec[i] = NULL;
			}
			pgVec.clear();

			if (meanRes != NULL)
				delete [] meanRes;
		}

		int lenMean; // Length of mean
		double * meanRes; 
		int numPGs;  // Number of modes [Note: pgVec.size() == numPGs]
		int lenPG;   // Length of principal geodesics [Note: lenMean == lenPG]

		int atomIndex;	// Used only for atom PGA data
		int figureIndex;	// Used only for atom PGA data

		std::vector<double *> pgVec;
};


#endif

