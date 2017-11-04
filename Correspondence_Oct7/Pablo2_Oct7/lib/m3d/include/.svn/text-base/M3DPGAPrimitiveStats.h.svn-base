#ifndef M3D_PRIMITIVE_PGA_STATS_H
#define M3D_PRIMITIVE_PGA_STATS_H


#include "M3DPGAStats.h"
#include "matrix.h"


// Atom PGA 


class PGData;

class M3DPGAPrimitiveStats : public M3DPGAStats 
{
//friend class M3DObjectFile;

public:
	void resetPGs();
	M3DPGAPrimitiveStats();
	
	M3DPGAPrimitiveStats * copyPtr(){	
		return new M3DPGAPrimitiveStats(*this);
	}

	virtual ~M3DPGAPrimitiveStats();

	PGData * getPGDataPtr(int figId, int atomId) ;
	M3DPrimitive *getDeltaMeanPrimitivePtr(int figId, int atomId);
	M3DObject * getAlignReferenceObj(){return alignReferenceObj;};

    virtual bool pgaDeform(M3DPrimitive * beforePrim, const Vector & vals, int figId, int atomId,
		M3DPrimitive * afterPrim);

	virtual bool doPGADeform(M3DPrimitive * beforePrim, const Vector & val, int figId, int atomId,
		M3DPrimitive * afterprim);

	// PGA Deformation applied to mean primitive
//	virtual bool pgaDeformToPrim(M3DPrimitive * prim, M3DPrimitive * meanPrim, const double * v);
	
	int getNumAtoms (int figId) {return numAtoms[figId];}

	void setMeanObj(M3DObject * obj) { alignReferenceObj = obj; }
	
	std::vector <int> getNumAtoms(){ return numAtoms; }
	void setNumAtoms(int numAtom){ numAtoms.push_back(numAtom); }
	
	int getNumFigs(){ return numFigs; }
    void setNumFigs( int figNum){ numFigs = figNum; }

	void setDeltaMeanPrimitives( M3DPrimitive * primitivePtr){ deltaMeanPrimitives.push_back(primitivePtr); }
    std::vector <M3DPrimitive *>  getDeltaMeanPrimitives(){ return deltaMeanPrimitives; }

    std::vector<PGData *> getPGs(){ return PGs;}
	void setPGs(PGData * pg){ PGs.push_back(pg); }

protected:
	//std::vector<PGSet *> PGSets;        // PGSets[i] has header info of PGs[i] 
	std::vector<PGData *> PGs;	// The i-th primitive PGAstats
	std::vector <M3DPrimitive *> deltaMeanPrimitives;	// The i-th delta mean atom
	std:: vector <int> numAtoms;	// Number of atoms for each figure
	int numFigs;
	M3DObject * alignReferenceObj;

};



#endif

