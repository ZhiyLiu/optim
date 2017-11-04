#ifndef M3D_TUBE_FIGURE_RENDERER_H
#define M3D_TUBE_FIGURE_RENDERER_H

#include "M3DFigureRenderer.h"

class M3DTubeFigureRenderer : public M3DFigureRenderer
{
protected:
	double tipVs[200][3];	// ... and the position on the circumference selector.
public:
    M3DTubeFigureRenderer(M3DTubeFigure * fig = NULL);
    ~M3DTubeFigureRenderer();

//	void setFigurePtr(M3DTubeFigure * fig);

    virtual void render();
    virtual void renderLandmarks();
};

#endif

