
#include <math.h>
#include <fstream>
#include <iostream>
#include "matrix.h"
#include "M3DPGA.h"
#include "M3DPGAPrimitiveStats.h"


using namespace std;


M3DPGAPrimitiveStats::M3DPGAPrimitiveStats()
{
	numFigs = 0;
}

/*
// not tested yet
M3DPGAPrimitiveStats::M3DPGAPrimitiveStats( const M3DPGAPrimitiveStats * primitiveStats){

     PGs = primitiveStats->PGs;	// The i-th primitive PGAstats
	deltaMeanPrimitives = primitiveStats->deltaMeanPrimitives;	// The i-th delta mean atom
	numAtoms = primitiveStats->numAtoms;	// Number of atoms of the i-th figure
	numFigs = primitiveStats->numFigs;
 
}
*/


M3DPGAPrimitiveStats::~M3DPGAPrimitiveStats()
{
  resetPGs();
}

void M3DPGAPrimitiveStats::resetPGs()
{
	int i;
	int numPGs = PGs.size();
	
	for (i = 0; i < numPGs; i++) {	
		if (! PGs.empty()) {
			if (PGs[i] != NULL)
				delete PGs[i];
		}
	}
	
	PGs.clear();	
	
}

PGData * M3DPGAPrimitiveStats::getPGDataPtr(int figId, int atomId) { 
		int i, order;
		order = atomId;
		for (i = 0; i < figId;i ++)
		    order += numAtoms[i] ;		
		 return PGs[order]; 
	}

M3DPrimitive * M3DPGAPrimitiveStats ::getDeltaMeanPrimitivePtr(int figId, int atomId)
{
	int i, order;
		order = atomId;
		for (i = 0; i < figId ; i ++)
		    order += numAtoms[i] ;
	return deltaMeanPrimitives[order];
}

bool M3DPGAPrimitiveStats::doPGADeform(M3DPrimitive * beforePrim, const Vector &val,
	int figId, int atomId, M3DPrimitive * afterprim)
{
	PGData * pgData = getPGDataPtr(figId, atomId);
	if (pgData == NULL) {
		cout << "Null PG in M3DPGAStats::doPGADeform(..)\n";
		return false;
	}

	// Compare the size of vals and number of pgs 
	/*	numParams = pgData->numPGs;
	if (val.size() < numParams)
		for (i = val.size(); i < numParams ; i++)
			val.push_back(0.0);
	*/
	if (beforePrim == NULL)
		return false;

	return pgaDeform(beforePrim, val, figId, atomId, afterprim);
}

bool M3DPGAPrimitiveStats::pgaDeform(M3DPrimitive * beforePrim,
	const Vector & vals, int figId, int atomId, M3DPrimitive * afterPrim)
{
	int numParams,
		size,
		i, j;

	if (figId< 0)
		return false;
#ifdef DEBUG

	cout << "x0=" << vals(0) << "  x1=" << vals(1) << 
		"  x2=" << vals(2) << "  x3=" << vals(3) << "  x4=" << vals(4) << 
        "  x5=" << vals(5) << "  x6=" << vals(6) << "  x7=" << vals(7) <<
		 '\n';
#endif


	PGData * pgData = getPGDataPtr(figId, atomId);
	if (beforePrim == NULL || afterPrim == NULL  || pgData == NULL )
        return false;

	// Compare the number of principal geodesic components and the size of vals
    if (pgData->numPGs > vals.size())
        numParams = vals.size();
    else
        numParams = pgData->numPGs;

	size = pgData->lenPG;

	double * v1 = new double[size];
    for (i = 0; i < size; i++)
        v1[i] = 0.0;

	for (i = 0; i < numParams; i++)
        for(j = 0; j < size; j++) {
			//cout << "i: " << i << " j: "<<j << "v1[j]:" <<v1[j]<<
			//"vals(i)"<<vals(i)<<" pg: " << pgData->pgVec[i][j] <<endl;

            v1[j] += (vals(i) * ((pgData->pgVec[i])[j]));
			//cout << "after: v1[j]:" << v1[j] << endl;
		}

	if (! pgaDeformToPrim(afterPrim, beforePrim, v1)) {	
		delete [] v1;
		cout << "Error: check pgaDeformToSet()\n";
		return false;
	}

    delete [] v1;

	return true;
}

