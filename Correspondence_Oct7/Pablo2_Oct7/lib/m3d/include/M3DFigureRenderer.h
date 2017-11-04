#ifndef M3D_FIGURE_RENDERER_H
#define M3D_FIGURE_RENDERER_H

#include "M3DPrimitiveRenderer.h"
#include "M3DFigure.h"
#include <vector>

// return a vector of landmark indices that involve this atom.
void findAtomsLandmarkIndices(M3DFigure *fig, int atomIndex, std::vector<int> &lmIndices);


class M3DFigureRenderer
{
public:
    M3DFigureRenderer(M3DFigure* _fig = NULL);
    virtual ~M3DFigureRenderer();

    void setStartingSelectionId(int id) { startingSelectionId = id; }
    void setPrimitiveRenderer(M3DPrimitiveRenderer * pr)
    {
        primitiveRenderer = pr;
    }

    void setPointScaleFactor(double factor)
    {
        if(primitiveRenderer != NULL) 
            primitiveRenderer->setPointScaleFactor(factor);
    }

    void setColor(const float * color);

    virtual void render() = 0;
    virtual void renderLandmarks() = 0;
    virtual void renderFigureName();

protected:
    M3DFigure * fig;
	float lineWidth;

    M3DPrimitiveRenderer * primitiveRenderer;

    int startingSelectionId;
	bool meshLineType;
	float rgb[3];
	// indexed by [atomIndex][instance]
	//  allows up to 3 named tips per atom
	char *tipNames[MAX_ATOMS_IN_FIGURE][3];	// up to 200 atoms
	double tipTs[MAX_ATOMS_IN_FIGURE][3];	// ...and their spoke-end selectors;
							// LATER: this may grow to a full (u,v,t) coord,
							//  or perhaps some other coord system
};

#endif

