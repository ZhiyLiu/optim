#ifndef INTERFIGURAL_CONSTRAINTS_H
#define INTERFIGURAL_CONSTRAINTS_H

#include <vector>

/*  Class InterfiguralConstraints.

    This is a container class containing two equal length lists, a list
	of figures and a list of corresponding distances.  The listed figures,
	considered to be governed figures, are constrained, in the registration
	and deformation stages of the program, to be the stored distance from
	the figure containing the instance of this class, which is termed the
	governor.

    The distance associated with each constrained (governed) figure implies
	a set of boundary points on the surface of the governed figure, all of
	which fall within the specified distance.

	Rendering the constrained objects (the "partial" rendering code in
	M3DObjectSurfaceRenderer in seurat) consists of looping over the lists
	and checking the distances against a cutoff distance to decide whether
	or not to render the boundary points.

    Every M3DFigure contains one instance of this class for use when
	it is a governor.

*/

class InterfiguralConstraints
{

public:

    InterfiguralConstraints() { }
    InterfiguralConstraints(const InterfiguralConstraints & ifc);
    InterfiguralConstraints & operator= (const InterfiguralConstraints & ifc);
    virtual ~InterfiguralConstraints() { }

	void clear();
	int addFigure(int figureID, float dist);
	bool deleteFigure(int figureID);
	bool updateFigure(int figureID, float dist);
	bool changeFigureId(int oldID, int newID);
	void remapFigureIds(int * map, int len);

	int size() const { return figures.size(); }
	float distance(int index) const { return dists[index]; }
	int figure(int index) const { return figures[index]; }

	void print(const char * msg = NULL);

protected:

    std::vector<int> figures;
    std::vector<float> dists;
};

#endif

