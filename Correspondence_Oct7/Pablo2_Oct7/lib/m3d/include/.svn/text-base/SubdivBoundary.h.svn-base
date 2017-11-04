// SubdivBoundary.h: interface for the SubdivBoundary class.
//
//////////////////////////////////////////////////////////////////////

#ifndef SUBDIV_BOUNDARY_H
#define SUBDIV_BOUNDARY_H

#define MAX_SUBDIV_LEVEL 6

class Displacements
{
private:
	int figureId;
	int level;			// subdivision level
	int numPts;			// number of points
	double * vals;		// array of displacements

public:
	Displacements(const Displacements * dp);
	Displacements(int level, int nPts);
	~Displacements();

	void setLevel(int _level) { level = _level; }
	void setNumPts(int n) { numPts = n; }
	void setVal(int index, double val) { vals[index] = val; }

	int getLevel() const { return level; }
	int getNumPts() const { return numPts; }
	double * getVals() const { return vals; }
};

class SubdivBoundary  
{
public:
	SubdivBoundary();
	SubdivBoundary(SubdivBoundary * const refBoundary);
	~SubdivBoundary();

	Displacements * getDisplacements(int level);
	void initDisplacements(int level, int nPts);
	void setDisplacements(Displacements * ds);

private:
	Displacements * displacements[MAX_SUBDIV_LEVEL];
};

#endif

