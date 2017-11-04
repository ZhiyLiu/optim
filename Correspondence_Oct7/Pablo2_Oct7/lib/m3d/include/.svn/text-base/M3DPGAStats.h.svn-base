#ifndef M3D_PGA_STATS_H
#define M3D_PGA_STATS_H

/*
 Class for hierarchical geometric statistics.
 Last modified : 3/9/05
*/

#include "M3DPGA.h"

class M3DObject;


// This class stores sets of global/residue principle geodesic statistics of multi-objects
class M3DPGAStats : public M3DPGA
{
	friend class M3DObjectFile;
//	friend class M3DPGAStatsFile;

public:

	enum geomStat_t { NotScaled = 0, Scaled };


	M3DPGAStats();
	~M3DPGAStats();
	M3DPGAStats(const M3DPGAStats &);

	void clear();
	void resetPGs();
	void resetObjsToMean();

	void setMeanObj(M3DObject * obj) { meanObj = obj; }
	void setLatestObj(M3DObject * obj);
	void setIntermedObj(M3DObject * obj);
	M3DObject * getMeanObjPtr()		{  return meanObj;  }
	M3DObject * getIntermedObjPtr()	{  return intermedObj;  }
	M3DObject * getLatestObjPtr()	{  return latestObj;  }

	int totalNumParams(int depth);	// Total no. of pg's used in all levels
	int getOrder(PGSet set)  {  return set.order;  }
	int getNumOfPGSets()	 {  return PGSets.size();  }

	int findPGSetId(int figureId);  // Order of PGData to use for figure stage optimization

	PGSet * getPGSetPtr(int order)   {  return PGSets[order]; }
	PGData * getPGDataPtr(int order) {  return PGs[order];  }
	PGData * getPredictionPtr(int order);
	void removeLastPGSet();

	SimilarityTransform3D & getBestTransform() {  return bestSimilarityTransform;  }
	void setBestTransform(const SimilarityTransform3D & transform) {
		bestSimilarityTransform = transform;
	}

	//wrapper 
	//object, vals, order are passed in from P3DControl
	bool applyMeanResidue(M3DObject * obj, int order);
	bool doPGADeform(M3DObject * obj, std::vector<double> vals, int order, bool predict = false);

	geomStat_t isScaled() { return type; }

	void rescale();
	void convertMeanDiff();

protected:

	virtual bool applyMeanResidue(M3DObject * obj, M3DObject * startObj, int order);

	// 1/20/05 -----------clean up up to here
	virtual bool pgaDeform(M3DObject * obj, M3DObject * startObj, std::vector<double> vals, int order,
		bool predict = false);
	virtual	bool pgaDeformToSet(M3DObject * obj, M3DObject * startObj, const double * v1, int order,
		const std::vector<double> * vals = NULL, bool predict = false);

	// PGA Deformation applied to mean primitive
	virtual	bool pgaDeformToPrim(M3DPrimitive * prim, M3DPrimitive * meanPrim, const double *v);

	int primitiveCount(int order);

	//1/19/05 prediction based on shape space
	void prediction(std::vector<double> &remVals, M3DObject * obj, std::vector<double> vals, int order);

	std::vector<PGSet *> PGSets;        // PGSets[i] has header info of PGs[i] 
	std::vector<PGData *> PGs; 
	std::vector<PGData *> Predictions;  // corresponding to PGs[1 ... number of PGs-1]

	//Initially, intermedObj == latestObj == meanObj(but separate copies)
	M3DObject * meanObj;	//mean ensemble object, For adaptive case, it is day 1 object currently.
	M3DObject * latestObj;	//intermedObj + most recently deformation made onto intermedObj in each stage
	M3DObject * intermedObj; //latestObj + meanResidue at each level

	// For the adaptive case
	M3DObject * meanDiffObj;

	// Counterintuitively, scaled puts PGAs in (xbar, r, theta1, theta2, eta) space, while
	//     Unscaled puts PGAs in (xbar, rbar*r, rbar*theta1, rbar*theta2, rbar*eta) space.
	//     The scale factor is 1/rbar, where rbar is the mean of r.
	geomStat_t type;	

	// CL: this is optimized in the first step, and kept for use in later optimization steps.
	SimilarityTransform3D bestSimilarityTransform;

	M3DPGAStats & operator=(const M3DPGAStats &);	// Not implemented
};

#endif
