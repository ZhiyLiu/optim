// SubdivBoundary.cpp: implementation of the SubdivBoundary class.
//
//////////////////////////////////////////////////////////////////////

#include "SubdivBoundary.h"
#include <float.h>

#include <iostream>
#include <fstream>

#include <string.h>

using namespace std;

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////
SubdivBoundary::SubdivBoundary()
{
	for (int i=0; i<MAX_SUBDIV_LEVEL; i++)
		displacements[i] = NULL; // intialize to NULL;
}

///////////////////////////////////////////////////////////////
// Copy constructor.  Only copies the displacement information
///////////////////////////////////////////////////////////////
SubdivBoundary::SubdivBoundary(SubdivBoundary * const refBoundary)
{
	for (int l=0; l<MAX_SUBDIV_LEVEL; l++)
	{
		if (refBoundary->displacements[l] != NULL)
			displacements[l] = new Displacements(refBoundary->displacements[l]);
		else
			displacements[l] = NULL;
	}
}

SubdivBoundary::~SubdivBoundary()
{
	for (int l=0; l < MAX_SUBDIV_LEVEL; l++)
	{
		if (displacements[l] != NULL)
			delete displacements[l];
		displacements[l] = NULL;
	}
}

//////////////////////////////////////////////////////////////////////
// Member functions
//////////////////////////////////////////////////////////////////////
Displacements * SubdivBoundary::getDisplacements(int level)
{
	if (level >= 0 && level < MAX_SUBDIV_LEVEL)
		return displacements[level];
	else 
		return NULL;
}

void SubdivBoundary::initDisplacements(int level, int nPts)
{
	if (level < 0 && level >= MAX_SUBDIV_LEVEL)
		return;

	if (displacements[level] != NULL)
	{
		double * vals = displacements[level] -> getVals();
		if (vals != NULL)
			delete [] vals;
	}

	displacements[level] = new Displacements(level, nPts);
}

void SubdivBoundary::setDisplacements(Displacements * ds)
{
	if (ds == NULL)
		return;

	int level = ds->getLevel();
	displacements[level] = ds;
}

/****************************************************************/
Displacements::Displacements(const Displacements * dp)
{
	figureId = dp -> figureId;
	level = dp -> level;
	numPts = dp -> numPts;

	if (dp -> vals != NULL)
	{
		vals = new double[numPts];
		memcpy(vals, dp -> vals, numPts*sizeof(double));
	}
	else
		vals = NULL;
}

Displacements::Displacements(int l, int nPts)
{
	figureId = 0;
	level = l;
	numPts = nPts;

	if (numPts > 0)
	{
		vals = new double[nPts];
		for (int i=0; i < numPts; i++)
			vals[i] = 0.0;
	}
	else 
		vals = NULL;
}

Displacements::~Displacements()
{
	if (vals != NULL)
		delete [] vals;
}


